/**
 * 三维演示主程序：三个正交切片 + 侧栏控制。
 *
 * 与二维版共用同一套"档 → 网格 → 模型"的构造方式，
 * 只是模型换成 AnelasticModel3，视野换成三个坐标面。
 */

import { SCALE_PRESETS3, presetById3, buildGrid3, buildModel3, measureThroughput3, autoPlan3 } from './model/plan3.js';
import { FIELD_LIST3 } from './viz/slices.js';
import { createSliceView } from './viz/slices.js';

const $ = (id) => document.getElementById(id);

const state = {
  presetId: 'conv3d',
  level: 0,
  speed: 20,
  field: 'cloud',
  paused: false,
  micro: true,
  rh: 0.78,
  sensible: 0,
  latent: 0,
  amp: 3,
  radius: 2000,
  wind: 0,
  shear: 0,
  follow: true,
  j: 0,
  i: 0,
  k: 0,
  pending: 0,
};

let model = null;
let views = null;
let lastWall = performance.now();
let stepMs = 0;
let lastStepMs = 0;
let stepsLastFrame = 0;

/* ----------------------------- URL 参数 ----------------------------- */

const urlq = new URLSearchParams(location.search);

/* ----------------------------- 模型构造 ----------------------------- */

function preset() { return presetById3(state.presetId); }

function makeModel() {
  const p = preset();
  const grid = buildGrid3(p, state.level);
  const wx = (p.surface?.xWidthKm ?? 2) * 1000;
  const wy = (p.surface?.yWidthKm ?? p.surface?.xWidthKm ?? 2) * 1000;
  const overrides = {
    base: { TSfc: 299.5, RHsfc: state.rh },
    micro: { enable: state.micro },
    turb: {
      surface: {
        sensible: state.sensible, latent: state.latent,
        xCenter: grid.Lx * 0.5, yCenter: grid.Ly * 0.5, xWidth: wx, yWidth: wy,
      },
    },
  };
  const m = buildModel3(grid, p, overrides);
  const b = p.bubble ?? { amplitude: 3, radius: 2000, zCenter: 1200 };
  m.initAtmosphere({
    amplitude: state.amp, radius: state.radius, zCenter: b.zCenter, noise: 0.02,
    wind: { u0: state.wind, shearU: state.shear / 1000 },
  });
  return m;
}

function resetModel() {
  model = makeModel();
  state.j = Math.floor(model.grid.ny / 2);
  state.i = Math.floor(model.grid.nx / 2);
  state.k = Math.floor(model.grid.nz / 3);
  state.pending = 0;
  syncSliceRanges();
  markActive();
}

/* ----------------------------- 控件 ----------------------------- */

function markActive() {
  for (const btn of document.querySelectorAll('#scales3 button')) btn.classList.toggle('on', btn.dataset.id === state.presetId);
  const p = preset();
  for (const btn of document.querySelectorAll('#levels3 button')) {
    btn.classList.toggle('on', Number(btn.dataset.lv) === state.level);
  }
  for (const btn of document.querySelectorAll('#speeds3 button')) btn.classList.toggle('on', Number(btn.dataset.sp) === state.speed);
  for (const btn of document.querySelectorAll('#fields3 button')) btn.classList.toggle('on', btn.dataset.id === state.field);
  $('note3').textContent = p.note;
  $('info3').textContent = `${p.levels[state.level].label} · ${(model ? model.grid.nCell : 0).toLocaleString()} 格点 · Δx ${(p.Lx / p.levels[state.level].nx).toFixed(0)} m`;
}

function buildButtons() {
  const box = $('scales3');
  for (const p of SCALE_PRESETS3) {
    const b = document.createElement('button');
    b.textContent = p.name;
    b.dataset.id = p.id;
    b.addEventListener('click', () => {
      state.presetId = p.id;
      state.level = 0;
      const bb = p.bubble ?? {};
      state.amp = bb.amplitude ?? 3;
      state.radius = bb.radius ?? 2000;
      state.sensible = p.surface?.sensible ?? 0;
      state.latent = p.surface?.latent ?? 0;
      syncSliders();
      buildLevelButtons();
      resetModel();
    });
    box.appendChild(b);
  }
}

function buildLevelButtons() {
  const box = $('levels3');
  box.innerHTML = '';
  preset().levels.forEach((lv, i) => {
    const b = document.createElement('button');
    b.textContent = lv.label;
    b.dataset.lv = String(i);
    b.addEventListener('click', () => { state.level = i; resetModel(); });
    box.appendChild(b);
  });
}

function buildSpeedButtons() {
  const box = $('speeds3');
  const opts = [5, 20, 60, 180];
  for (const sp of opts) {
    const b = document.createElement('button');
    b.textContent = `${sp}×`;
    b.dataset.sp = String(sp);
    b.addEventListener('click', () => { state.speed = sp; markActive(); });
    box.appendChild(b);
  }
}

function buildFieldButtons() {
  const box = $('fields3');
  for (const f of FIELD_LIST3) {
    const b = document.createElement('button');
    b.textContent = f.name;
    b.dataset.id = f.id;
    b.addEventListener('click', () => { state.field = f.id; markActive(); });
    box.appendChild(b);
  }
}

function syncSliders() {
  $('rh').value = String(state.rh);
  $('sensible').value = String(state.sensible);
  $('latent').value = String(state.latent);
  $('amp').value = String(state.amp);
  $('radius').value = String(state.radius);
  $('wind').value = String(state.wind);
  $('shear').value = String(state.shear);
  $('rhOut').textContent = state.rh.toFixed(2);
  $('sensibleOut').textContent = `${state.sensible} W/m²`;
  $('latentOut').textContent = `${state.latent} W/m²`;
  $('ampOut').textContent = `${state.amp.toFixed(1)} K`;
  $('radiusOut').textContent = `${(state.radius / 1000).toFixed(1)} km`;
  $('windOut').textContent = `${state.wind} m/s`;
  $('shearOut').textContent = `${state.shear} m/s·km⁻¹`;
  $('follow3').checked = state.follow;
}

function syncSliceRanges() {
  const g = model.grid;
  const mk = (id, max, val) => { const el = $(id); el.max = String(max); el.value = String(val); };
  mk('sliceY', g.ny - 1, state.j);
  mk('sliceX', g.nx - 1, state.i);
  mk('sliceZ', g.nz - 1, state.k);
  updateSliceLabels();
}

function updateSliceLabels() {
  const g = model.grid;
  $('sliceYOut').textContent = `${(g.yc(state.j) / 1000).toFixed(1)} km`;
  $('sliceXOut').textContent = `${(g.xc(state.i) / 1000).toFixed(1)} km`;
  $('sliceZOut').textContent = `${(g.zc(state.k) / 1000).toFixed(1)} km`;
}

/**
 * 三个切面跟随云体：取凝结物柱积分最大的 y / x，以及 qc 最多的层。
 * 环境风把云吹走之后，固定切面会切到空白处，看起来"什么都没长"。
 */
let followTick = 0;
function trackCloud(force = false) {
  if (!state.follow || !model) return;
  followTick++;
  if (!force && followTick % 8) return;
  const g = model.grid, S = model.S;
  const rowY = new Float64Array(g.ny);
  const rowX = new Float64Array(g.nx);
  let bestK = 0, bestLayer = -1, total = 0;
  for (let k = 0; k < g.nz; k++) {
    let layer = 0;
    for (let j = 0; j < g.ny; j++) {
      const off = g.nx * (j + g.ny * k);
      let s = 0;
      for (let i = 0; i < g.nx; i++) s += S.qc[off + i];
      rowY[j] += s;
      layer += s;
    }
    for (let i = 0; i < g.nx; i++) {
      let s = 0;
      for (let j = 0; j < g.ny; j++) s += S.qc[i + g.nx * (j + g.ny * k)];
      rowX[i] += s;
    }
    total += layer;
    if (layer > bestLayer) { bestLayer = layer; bestK = k; }
  }
  if (total < 1e-4) return;                    // 没有云就不动切面
  let bj = 0, bi = 0;
  for (let j = 1; j < g.ny; j++) if (rowY[j] > rowY[bj]) bj = j;
  for (let i = 1; i < g.nx; i++) if (rowX[i] > rowX[bi]) bi = i;
  state.j = bj; state.i = bi; state.k = bestK;
  $('sliceY').value = String(bj);
  $('sliceX').value = String(bi);
  $('sliceZ').value = String(bestK);
  updateSliceLabels();
}

/* ----------------------------- 主循环 ----------------------------- */

function stepSim() {
  if (!model) return;
  const dtNow = model.suggestDt();
  state.pending += state.speed * Math.min(80, performance.now() - lastWall) / 1000;
  if (state.pending > 60) state.pending = 60;     // 卡顿保护
  const t0 = performance.now();
  stepsLastFrame = 0;
  while (state.pending >= dtNow && stepsLastFrame < 8 && performance.now() - t0 < 45) {
    model.step(dtNow);
    state.pending -= dtNow;
    stepsLastFrame++;
  }
  stepMs = performance.now() - t0;
  if (stepsLastFrame > 0) lastStepMs = stepMs / stepsLastFrame;
}

function drawAll() {
  if (!model) return;
  trackCloud();
  const p = preset();
  const viz = p.viz;
  views.xz.draw(model, { field: state.field, plane: 'xz', index: state.j, viz });
  views.yz.draw(model, { field: state.field, plane: 'yz', index: state.i, viz });
  views.xy.draw(model, { field: state.field, plane: 'xy', index: state.k, viz });
}

let hudTick = 0;
function hud() {
  if (!model) return;
  if (hudTick++ % 6) return;
  const d = model.diagnostics();
  const g = model.grid;
  const mm = Math.floor(d.time / 60), ss = Math.floor(d.time % 60);
  $('hud3').innerHTML = [
    ['模拟时间', `${mm}:${String(ss).padStart(2, '0')}`],
    ['垂直速度', `${d.maxW.toFixed(1)} m/s`],
    ['水平风速', `${d.maxU.toFixed(1)} m/s`],
    ['云顶', d.cloudTop > 0 ? `${(d.cloudTop / 1000).toFixed(1)} km` : '—'],
    ['降水', `${d.precipRate.toFixed(1)} mm/h`],
    ['累计降水', `${d.precipAccum.toFixed(1)} mm`],
    ['CFL', d.cfl.toFixed(2)],
    ['每步耗时', `${lastStepMs.toFixed(1)} ms`],
    ['本帧步数', `${stepsLastFrame}`],
    ['网格', `${g.nx}×${g.ny}×${g.nz}`],
  ].map(([k, v]) => `<div class="row"><span>${k}</span><b>${v}</b></div>`).join('');
}

function frame() {
  const now = performance.now();
  if (!state.paused) stepSim();
  drawAll();
  hud();
  lastWall = now;
  requestAnimationFrame(frame);
}

/* ----------------------------- 交互 ----------------------------- */

function onCanvasClick(ev) {
  if (!model) return;
  const cv = ev.currentTarget;
  const rect = cv.getBoundingClientRect();
  const dpr = cv.width / rect.width;
  const xPix = (ev.clientX - rect.left) * dpr;
  const yPix = (ev.clientY - rect.top) * dpr;
  const padL = Math.round(38 * dpr), padB = Math.round(20 * dpr), padT = Math.round(16 * dpr);
  const plotW = Math.max(1, cv.width - padL);
  const plotH = Math.max(1, cv.height - padB - padT);
  const simX = Math.max(0, Math.min(model.grid.Lx, (xPix - padL) / plotW * model.grid.Lx));
  const simZ = Math.max(0, Math.min(model.grid.H, (1 - (yPix - padT) / plotH) * model.grid.H));
  model.addBubble({ x: simX, y: model.grid.yc(state.j), z: simZ, amplitude: state.amp, radius: state.radius });
}

function applyUrlState() {
  if (urlq.has('preset')) state.presetId = presetById3(urlq.get('preset')).id;
  if (urlq.has('level')) state.level = Number(urlq.get('level'));
  if (urlq.has('speed')) state.speed = Number(urlq.get('speed'));
  if (urlq.has('field')) state.field = urlq.get('field');
  if (urlq.has('rh')) state.rh = Number(urlq.get('rh'));
  if (urlq.has('sensible')) state.sensible = Number(urlq.get('sensible'));
  if (urlq.has('latent')) state.latent = Number(urlq.get('latent'));
  if (urlq.has('amp')) state.amp = Number(urlq.get('amp'));
  if (urlq.has('radius')) state.radius = Number(urlq.get('radius'));
  if (urlq.has('wind')) state.wind = Number(urlq.get('wind'));
  if (urlq.has('shear')) state.shear = Number(urlq.get('shear'));
  if (urlq.has('micro')) state.micro = urlq.get('micro') !== '0';
}

/* ----------------------------- 启动 ----------------------------- */

function boot() {
  applyUrlState();
  buildButtons();
  buildSpeedButtons();
  buildFieldButtons();
  buildLevelButtons();
  syncSliders();

  views = {
    xz: createSliceView($('viewXZ')),
    yz: createSliceView($('viewYZ')),
    xy: createSliceView($('viewXY')),
  };

  $('rh').addEventListener('input', () => { state.rh = Number($('rh').value); $('rhOut').textContent = state.rh.toFixed(2); resetModel(); });
  $('sensible').addEventListener('input', () => { state.sensible = Number($('sensible').value); $('sensibleOut').textContent = `${state.sensible} W/m²`; resetModel(); });
  $('latent').addEventListener('input', () => { state.latent = Number($('latent').value); $('latentOut').textContent = `${state.latent} W/m²`; resetModel(); });
  $('amp').addEventListener('input', () => { state.amp = Number($('amp').value); $('ampOut').textContent = `${state.amp.toFixed(1)} K`; });
  $('radius').addEventListener('input', () => { state.radius = Number($('radius').value); $('radiusOut').textContent = `${(state.radius / 1000).toFixed(1)} km`; });
  $('wind').addEventListener('input', () => { state.wind = Number($('wind').value); $('windOut').textContent = `${state.wind} m/s`; resetModel(); });
  $('shear').addEventListener('input', () => { state.shear = Number($('shear').value); $('shearOut').textContent = `${state.shear} m/s·km⁻¹`; resetModel(); });

  $('sliceY').addEventListener('input', () => { state.follow = false; $('follow3').checked = false; state.j = Number($('sliceY').value); updateSliceLabels(); });
  $('sliceX').addEventListener('input', () => { state.follow = false; $('follow3').checked = false; state.i = Number($('sliceX').value); updateSliceLabels(); });
  $('sliceZ').addEventListener('input', () => { state.follow = false; $('follow3').checked = false; state.k = Number($('sliceZ').value); updateSliceLabels(); });
  $('follow3').addEventListener('change', () => { state.follow = $('follow3').checked; });

  $('micro3').addEventListener('change', () => { state.micro = $('micro3').checked; resetModel(); });
  $('pause3').addEventListener('click', () => { state.paused = !state.paused; $('pause3').textContent = state.paused ? '继续' : '暂停'; });
  $('reset3').addEventListener('click', () => resetModel());

  $('probe3').addEventListener('click', () => {
    const p = preset();
    const r = measureThroughput3(p, state.level, 120);
    $('probeOut3').textContent = `实测 ${r.msPerStep.toFixed(1)} ms/步（${Math.round(r.cellStepsPerMs)} 格点·步/ms）`;
  });
  $('auto3').addEventListener('click', () => {
    const p = preset();
    const probe = measureThroughput3(p, state.level, 100);
    const plan = autoPlan3(probe, p, { frameMs: 33, simSpeed: state.speed, dt: model ? model.suggestDt() : 6 });
    state.level = plan.levelIndex;
    resetModel();
    buildLevelButtons();
    markActive();
    $('probeOut3').textContent = `定档 ${plan.label}（上限 ${plan.maxCells} 格点，预计 ${plan.estMsPerFrame.toFixed(0)} ms/帧）`;
  });

  $('viewXZ').addEventListener('click', onCanvasClick);

  resetModel();
  if (urlq.has('warm')) {
    const w = Number(urlq.get('warm'));
    if (w > 0) {
      const r = model.advance(w, 0);
      if (r.steps > 0) lastStepMs = r.ms / r.steps;
    }
  }
  trackCloud(true);
  $('micro3').checked = state.micro;
  lastWall = performance.now();
  requestAnimationFrame(frame);
}

if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', boot);
else boot();
