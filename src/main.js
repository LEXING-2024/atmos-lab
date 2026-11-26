/**
 * 浏览器演示主程序。
 *
 * 只做三件事：驱动时间推进、渲染、把控件接到模型上。
 * 所有物理与数值都在 src/model、src/physics 里，换界面不影响它们。
 */

import { SCALE_PRESETS, presetById, buildGrid, buildModel, measureThroughput, autoPlan } from './model/plan.js';
import { createRenderer, FIELD_LIST } from './viz/render.js';
import { createTurbulence } from './physics/turb.js';
import { remap2D } from './core/remap.js';

const canvas = document.getElementById('view');
const renderer = createRenderer(canvas);

const state = {
  presetId: 'cloud',
  level: 1,
  speed: 20,
  field: 'cloud',
  scale: 1,
  paused: false,
  rh: 0.78,
  tsfc: 299.5,
  sensible: 0,
  latent: 0,
  amp: 3,
  radius: 4000,
  zc: 2000,
  micro: true,
  filter: 0.02,
  showAxes: true,
  budgetMs: 14,
};

let model = null;
let fps = 0;
let frameMs = 0;
let lastT = 0;
let hudT = 0;
let stepsLastFrame = 0;
const perf = { cellStepsPerMs: 0, measured: false };

const presetInit = {
  macro: { amp: 4, radius: 80000, zc: 4000 },
  meso: { amp: 3, radius: 20000, zc: 2500 },
  cloud: { amp: 3, radius: 4000, zc: 2000 },
  les: { amp: 1.5, radius: 400, zc: 250 },
};

/* URL 参数：?preset=cloud&level=1&speed=20&field=cloud&warm=1200
 * warm 是"预积分秒数"——重置后先离线跑这么久，省得等它慢慢长起来；
 * 一个 advance(warm, 0) 调用用的是模型自己的 dt（最大 8 s），比按帧推进快几十倍。
 *
 * 注意：这里只解析、不落地。真正的赋值放在启动时（applyUrlState），
 * 因为 bindSlider 初始化时会把滑块里的 DOM 值写回 state，那会覆盖这里的设定。 */
const urlq = new URLSearchParams(typeof location !== 'undefined' ? location.search : '');
let warmUp = Math.max(0, Number(urlq.get('warm')) || 0);

function applyUrlState() {
  if (urlq.has('preset') && SCALE_PRESETS.some((p) => p.id === urlq.get('preset'))) {
    state.presetId = urlq.get('preset');
  }
  if (urlq.has('level')) state.level = Math.max(0, Math.min(2, Number(urlq.get('level'))));
  if (urlq.has('speed')) state.speed = Number(urlq.get('speed')) || state.speed;
  if (urlq.has('field') && FIELD_LIST.some((f) => f.id === urlq.get('field'))) state.field = urlq.get('field');

  const preset = presetById(state.presetId);
  const ini = presetInit[preset.id] || presetInit.cloud;
  state.amp = ini.amp; state.radius = ini.radius; state.zc = ini.zc;
  state.micro = preset.micro !== false;

  const sf = preset.surface || {};
  state.sensible = sf.sensible || 0;
  state.latent = sf.latent || 0;
  if (urlq.has('sensible')) state.sensible = Number(urlq.get('sensible'));
  if (urlq.has('latent')) state.latent = Number(urlq.get('latent'));
  if (urlq.has('amp')) state.amp = Number(urlq.get('amp'));
  if (urlq.has('radius')) state.radius = Number(urlq.get('radius'));
}

/* ------------------------------ 建模型 ------------------------------ */

function currentOverrides() {
  const preset = presetById(state.presetId);
  const sf = preset.surface || {};
  return {
    base: { TSfc: state.tsfc, RHsfc: state.rh },
    micro: { enable: state.micro },
    turb: {
      surface: {
        sensible: state.sensible,
        latent: state.latent,
        // 必须给条带中心：xCenter 为 null 时 facCol 全域为 1，等于均匀加热，
        // 只有湍流没有环流，看不出任何中尺度结构。
        xCenter: preset.Lx / 2,
        xWidth: (sf.xWidthKm || 4) * 1000,
      },
    },
    top: { filter: state.filter },
  };
}

function resetModel() {
  const preset = presetById(state.presetId);
  const grid = buildGrid(preset, state.level);
  model = buildModel(grid, preset, currentOverrides());
  model.initAtmosphere({
    amplitude: state.amp, radius: state.radius, zCenter: state.zc, noise: 0.02,
  });
  if (warmUp > 0) {
    model.advance(warmUp, 0);   // 一次调用用模型自己的 dt，比按帧推进快几十倍
    warmUp = 0;
  }
  stepsLastFrame = 0;
}

/** 只改分辨率：保留场做守恒重映射 */
function changeLevel(level) {
  const preset = presetById(state.presetId);
  const grid = buildGrid(preset, level);
  if (!model) { resetModel(); return; }
  const m = buildModel(grid, preset, currentOverrides());
  model = remapInto(model, m);
}

/** 用质量加权重映射把旧场搬到新网格 */
function remapInto(oldModel, newModel) {
  const g0 = oldModel.grid, g1 = newModel.grid;
  const b0 = oldModel.base, b1 = newModel.base;
  const tmp = new Float64Array(g1.nCell);
  const wsrc = new Float64Array(g0.nCell);
  for (const f of ['u', 'th', 'qv', 'qc', 'qr']) {
    for (let k = 0; k < g0.nz; k++) {
      const off = g0.nx * k;
      for (let i = 0; i < g0.nx; i++) wsrc[off + i] = b0.rho0[k] * oldModel.S[f][off + i];
    }
    remap2D(wsrc, g0.nx, g0.nz, g1.nx, g1.nz, tmp);
    for (let k = 0; k < g1.nz; k++) {
      const off = g1.nx * k;
      for (let i = 0; i < g1.nx; i++) newModel.S[f][off + i] = tmp[off + i] / b1.rho0[k];
    }
  }
  const wSrc = new Float64Array(g0.nW), wOut = new Float64Array(g1.nW);
  for (let k = 0; k <= g0.nz; k++) {
    const off = g0.nx * k;
    for (let i = 0; i < g0.nx; i++) wSrc[off + i] = b0.rhoW[k] * oldModel.S.w[off + i];
  }
  remap2D(wSrc, g0.nx, g0.nz + 1, g1.nx, g1.nz + 1, wOut);
  for (let k = 0; k <= g1.nz; k++) {
    const off = g1.nx * k;
    for (let i = 0; i < g1.nx; i++) newModel.S.w[off + i] = wOut[off + i] / b1.rhoW[k];
  }
  newModel.time = oldModel.time;
  newModel.precipAccum = oldModel.precipAccum;
  newModel.waterInput = oldModel.waterInput;
  newModel.initMass = oldModel.initMass;
  return newModel;
}

/* ------------------------------ 控件 ------------------------------ */

function markActive(container, idx) {
  const kids = container.querySelectorAll('button');
  kids.forEach((b, i) => b.classList.toggle('on', i === idx));
}

const scalesEl = document.getElementById('scales');
SCALE_PRESETS.forEach((p, idx) => {
  const b = document.createElement('button');
  b.textContent = p.name;
  b.title = p.note;
  b.onclick = () => {
    state.presetId = p.id;
    markActive(scalesEl, idx);
    buildLevelButtons();
    const ini = presetInit[p.id];
    state.amp = ini.amp; state.radius = ini.radius; state.zc = ini.zc;
    state.micro = p.micro !== false;
    document.getElementById('micro').checked = state.micro;
    document.getElementById('presetNote').textContent = p.note;
    const sf = p.surface || {};
    state.sensible = sf.sensible || 0;
    state.latent = sf.latent || 0;
    syncSurface();
    syncSliders();
    resetModel();
    updateInfo();
  };
  scalesEl.appendChild(b);
});
markActive(scalesEl, SCALE_PRESETS.findIndex((p) => p.id === state.presetId));

const levelsEl = document.getElementById('levels');
function buildLevelButtons() {
  levelsEl.innerHTML = '';
  const preset = presetById(state.presetId);
  preset.levels.forEach((lv, idx) => {
    const b = document.createElement('button');
    b.textContent = lv.label;
    if (idx === state.level) b.classList.add('on');
    b.onclick = () => {
      state.level = idx;
      changeLevel(idx);
      markActive(levelsEl, idx);
      updateInfo();
    };
    levelsEl.appendChild(b);
  });
}
buildLevelButtons();

const speeds = [1, 5, 20, 60, 240];
const speedEl = document.getElementById('speeds');
speeds.forEach((s, idx) => {
  const b = document.createElement('button');
  b.textContent = `${s}×`;
  if (s === state.speed) b.classList.add('on');
  b.onclick = () => { state.speed = s; markActive(speedEl, idx); };
  speedEl.appendChild(b);
});

const fieldsEl = document.getElementById('fields');
FIELD_LIST.forEach((f, idx) => {
  const b = document.createElement('button');
  b.textContent = f.name;
  if (f.id === state.field) b.classList.add('on');
  b.onclick = () => { state.field = f.id; markActive(fieldsEl, idx); };
  fieldsEl.appendChild(b);
});

const bindSlider = (id, key, fmt, onChange) => {
  const el = document.getElementById(id);
  const out = document.getElementById(`${id}Out`);
  const upd = (initial) => {
    state[key] = Number(el.value);
    if (out) out.textContent = fmt(state[key]);
    // 初始化时不能触发 onChange：那会在启动阶段多重建一次模型，
    // 把 URL 里的预积分（warm）用在一次随后被丢弃的模型上。
    if (!initial && onChange) onChange();
  };
  el.addEventListener('input', () => upd(false));
  upd(true);
};

bindSlider('rh', 'rh', (v) => `${(v * 100).toFixed(0)}%`, () => {
  // 湿度改的是基态，必须重建
  resetModel();
});
bindSlider('sensible', 'sensible', (v) => `${v.toFixed(0)} W/m²`, applyLive);
bindSlider('latent', 'latent', (v) => `${v.toFixed(0)} W/m²`, applyLive);
bindSlider('amp', 'amp', (v) => `${v.toFixed(1)} K`);
bindSlider('radius', 'radius', (v) => `${(v / 1000).toFixed(1)} km`);

document.getElementById('micro').onchange = (e) => { state.micro = e.target.checked; applyLive(); };
document.getElementById('pause').onclick = () => {
  state.paused = !state.paused;
  document.getElementById('pause').textContent = state.paused ? '继续' : '暂停';
};
document.getElementById('reset').onclick = () => resetModel();
document.getElementById('probe').onclick = () => runProbe();
document.getElementById('auto').onclick = () => autoPick();

function applyLive() {
  if (!model) return;
  model.cfg.micro.enable = state.micro;
  model.cfg.turb.surface.sensible = state.sensible;
  model.cfg.turb.surface.latent = state.latent;
  model.turb = createTurbulence(model.grid, model.cfg.turb);
}

function syncSliders() {
  document.getElementById('amp').value = state.amp;
  document.getElementById('ampOut').textContent = `${state.amp.toFixed(1)} K`;
  const r = Math.max(20, state.radius);
  document.getElementById('radius').value = r;
  document.getElementById('radiusOut').textContent = `${(r / 1000).toFixed(1)} km`;
}

/** 把地表通量的实际值写回滑块（换档时会跟着变，所见即所用） */
function syncSurface() {
  document.getElementById('sensible').value = state.sensible;
  document.getElementById('sensibleOut').textContent = `${state.sensible.toFixed(0)} W/m²`;
  document.getElementById('latent').value = state.latent;
  document.getElementById('latentOut').textContent = `${state.latent.toFixed(0)} W/m²`;
}

/* --------------------------- 性能探测/定档 --------------------------- */

function runProbe() {
  const preset = presetById(state.presetId);
  const t = measureThroughput(preset, state.level, 120, currentOverrides());
  perf.cellStepsPerMs = t.cellStepsPerMs;
  perf.measured = true;
  const dt = model ? Math.max(model.lastDt, 1) : 4;
  const plan = autoPlan({ cellStepsPerMs: t.cellStepsPerMs }, preset, {
    frameMs: 16, simSpeed: state.speed, dt,
  });
  document.getElementById('probeOut').textContent =
    `实测 ${(t.cellStepsPerMs / 1000).toFixed(2)}k 格点·步/ms ｜ ${t.msPerStep.toFixed(2)} ms/步 ｜ 当前倍速建议 ${plan.label}（每帧 ${plan.stepsPerFrame} 步，约 ${plan.estMsPerFrame.toFixed(1)} ms）`;
  return plan;
}

function autoPick() {
  const plan = runProbe();
  if (plan.levelIndex !== state.level) {
    state.level = plan.levelIndex;
    changeLevel(state.level);
    markActive(levelsEl, state.level);
  }
  updateInfo();
}

/* ------------------------------ 主循环 ------------------------------ */

function updateInfo() {
  if (!model) return;
  const g = model.grid;
  const preset = presetById(state.presetId);
  document.getElementById('info').textContent =
    `${preset.name}档 ｜ 域 ${(g.Lx / 1000).toFixed(g.Lx > 1e5 ? 0 : 1)} km × ${(g.H / 1000).toFixed(1)} km ｜ `
    + `Δx ${g.dx.toFixed(g.dx < 100 ? 1 : 0)} m Δz ${g.dz.toFixed(g.dz < 100 ? 1 : 0)} m ｜ `
    + `格点 ${g.nx}×${g.nz} = ${g.nCell}`;
}

function draw() {
  renderer.draw(model, { field: state.field, viz: presetById(state.presetId).viz, showAxes: state.showAxes });
}

function frame(t) {
  const dtReal = lastT ? Math.min(0.25, (t - lastT) / 1000) : 0;
  lastT = t;
  fps = fps ? fps * 0.9 + (1 / Math.max(dtReal, 1e-3)) * 0.1 : 1 / Math.max(dtReal, 1e-3);
  frameMs = frameMs ? frameMs * 0.9 + dtReal * 1000 * 0.1 : dtReal * 1000;

  if (!state.paused && model) {
    const want = dtReal * state.speed;
    const r = model.advance(want, state.budgetMs);
    stepsLastFrame = r.steps;
  }

  draw();

  if (t - hudT > 250) {
    hudT = t;
    updateHud();
  }
  requestAnimationFrame(frame);
}

function updateHud() {
  if (!model) return;
  const g = model.grid, d = model.diagnostics();
  const drift = (d.totalWater + model.precipAccum * g.Lx - model.waterInput - model.initMass)
    / (model.initMass || 1) * 100;
  const simMin = model.time / 60;
  const el = document.getElementById('hud');
  el.innerHTML = [
    ['FPS', fps.toFixed(0)],
    ['帧耗时', `${frameMs.toFixed(1)} ms`],
    ['Δt', `${model.lastDt.toFixed(2)} s`],
    ['步/帧', stepsLastFrame],
    ['模拟时间', `${simMin.toFixed(1)} min`],
    ['吞吐', `${(g.nCell * stepsLastFrame / Math.max(frameMs, 1) / 1000).toFixed(2)}k 格点·步/ms`],
    ['max w', `${d.maxW.toFixed(1)} m/s`],
    ['云顶', `${(d.cloudTop / 1000).toFixed(1)} km`],
    ['降水率', `${d.precipRate.toFixed(1)} mm/h`],
    ['累计降水', `${model.precipAccum.toFixed(2)} mm`],
    ['散度残差', d.divResidual.toExponential(1)],
    ['水量漂移', `${drift.toFixed(2)} %`],
  ].map(([k, v]) => `<div class="row"><span>${k}</span><b>${v}</b></div>`).join('');
}

/* 点击画面注入热泡 */
canvas.addEventListener('pointerdown', (e) => {
  if (!model) return;
  const p = renderer.pick(model, e.clientX, e.clientY, state.showAxes);
  model.addBubble({ x: p.x, z: p.z, amplitude: state.amp, radius: state.radius });
});

/* 启动 */
applyUrlState();            // 必须在 bindSlider 初始化之后，否则会被 DOM 默认值覆盖
markActive(scalesEl, SCALE_PRESETS.findIndex((p) => p.id === state.presetId));
buildLevelButtons();        // 档位标签跟着预设走，URL 换档后要重建
markActive(speedEl, speeds.indexOf(state.speed));
markActive(fieldsEl, FIELD_LIST.findIndex((f) => f.id === state.field));
syncSurface();
syncSliders();
resetModel();
document.getElementById('micro').checked = state.micro;
document.getElementById('presetNote').textContent = presetById(state.presetId).note;
updateInfo();
requestAnimationFrame((t) => { lastT = t; frame(t); });
