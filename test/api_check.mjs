/** 临时：按 API 文档逐条核对签名与返回值 */
import { Grid } from '../src/core/grid.js';
import { AnelasticModel, DEFAULT_CONFIG, resizeModel } from '../src/model/anelastic.js';
import { makeBase } from '../src/model/base.js';
import { SCALE_PRESETS, presetById, buildGrid, buildModel, measureThroughput, autoPlan, frameCost } from '../src/model/plan.js';
import { colorize, FIELD_LIST, DEFAULT_VIZ } from '../src/viz/render.js';
import { createMicrophysics } from '../src/physics/micro.js';
import { repairMoisture, clampMoisture, maxAbs } from '../src/model/state.js';
import { qsat, esat, exner, satIncrement, CONST } from '../src/physics/thermo.js';
import { remapUniform, remap2D } from '../src/core/remap.js';
import { createPoissonSolver } from '../src/core/poisson.js';

const ok = [];
const bad = [];
const chk = (name, cond, extra = '') => (cond ? ok : bad).push(`${name}${extra ? ' — ' + extra : ''}`);

// 1. Grid
const g = new Grid({ nx: 128, nz: 64, Lx: 32000, H: 16000 });
chk('Grid 字段', g.nCell === 8192 && g.nW === 128 * 65 && Math.abs(g.dz - 250) < 1e-9);
chk('Grid.sIdx/wIdx', g.sIdx(3, 2) === 3 + 256 && g.wIdx(0, 64) === 8192);
chk('Grid.wrap 负数', g.wrap(-1) === 127);
let threw = false;
try { new Grid({ nx: 100, nz: 32, Lx: 1, H: 1 }); } catch { threw = true; }
chk('Grid 非 2 的幂抛错', threw);

// 2. base / thermo
const base = makeBase(g, {});
chk('makeBase 长度', base.p0.length === 64 && base.rhoW.length === 65);
chk('makeBase 默认 lapse=K/m', Math.abs((base.T0[0]) - (299.5 - 125 * 0.007)) < 0.05, `T0[0]=${base.T0[0].toFixed(2)}`);
chk('esat/qsat/exner', esat(300) > 3000 && qsat(300, 1e5) > 0.01 && exner(1e5) === 1);
chk('satIncrement 凝结为正', satIncrement(300, 1e5, 0.03, 0) > 0, `qs=${qsat(300, 1e5).toFixed(4)}`);
chk('satIncrement 蒸发为负', satIncrement(300, 1e5, 0.01, 1e-3) < 0);
chk('CONST', CONST.g === 9.81 && CONST.p00 === 100000);

// 3. 泊松
const solver = createPoissonSolver(g);
const out = new Float64Array(g.nCell);
chk('solver 暴露三件', !!solver.fft && typeof solver.project === 'function' && typeof solver.computeDiv === 'function');
chk('computeDiv 返回标量', typeof solver.computeDiv(new Float64Array(g.nCell), new Float64Array(g.nW), base, out) === 'number');

// 4. remap
const src = new Float64Array([1, 1, 1, 1]);
const dst = new Float64Array(2);
remapUniform(src, 4, 2, dst);
chk('remapUniform 送 1→1', dst[0] === 1 && dst[1] === 1);
chk('remap2D 可调用', (() => { const o = new Float64Array(4); remap2D(new Float64Array(4), 2, 2, 1, 4, o); return o.length === 4; })());

// 5. 四档 × 三档全部可构造
for (const p of SCALE_PRESETS) {
  for (let i = 0; i < 3; i++) {
    const grid = buildGrid(p, i);
    const m = buildModel(grid, p, {
      turb: { surface: { sensible: p.surface?.sensible || 0, latent: 0, xCenter: grid.Lx / 2, xWidth: 1e4 } },
    });
    chk(`${p.id} L${i} 构造`, m.grid.nx === p.levels[i].nx);
  }
}

// 6. 模型 API
const m = new AnelasticModel(g, { base: { TSfc: 299.5, RHsfc: 0.78 }, micro: { enable: true } });
m.initAtmosphere({ amplitude: 3, radius: 4000, zCenter: 2000, noise: 0.02 });
const r = m.advance(60, 0);
chk('advance 返回字段', typeof r.steps === 'number' && typeof r.dt === 'number' && typeof r.simAdvanced === 'number' && typeof r.ms === 'number');
chk('advance 推进量', Math.abs(r.simAdvanced - 60) < 1e-6, `simAdvanced=${r.simAdvanced}`);
const d = m.diagnostics();
const keys = ['time', 'maxW', 'maxU', 'cloudTop', 'precipRate', 'precipAccum', 'cfl', 'divResidual', 'totalWater', 'waterInput', 'cells'];
chk('diagnostics 字段齐', keys.every((k) => k in d));
const p1 = m.probe(16000, 2500);
chk('probe 字段', ['u', 'w', 'th', 'qv', 'qc', 'qr', 'T'].every((k) => k in p1));
m.addBubble({ x: 8000, z: 1500, amplitude: 1, radius: 1000 });
const snap = m.snapshot();
chk('snapshot 字段', ['u', 'w', 'th', 'qv', 'qc', 'qr', 'time', 'precipAccum'].every((k) => k in snap));
chk('snapshot 不泄内部引用', snap.th !== m.S.th);
m.restore(snap);
chk('suggestDt 合理', m.suggestDt() > 0.05 && m.suggestDt() <= 8);

// 7. resizeModel 默认模式 = reset
const m2 = resizeModel(m, new Grid({ nx: 256, nz: 96, Lx: 32000, H: 16000 }));
chk('resizeModel 默认 reset（time 归零）', m2.time === 0, `time=${m2.time}`);
const m3 = resizeModel(m, new Grid({ nx: 256, nz: 96, Lx: 32000, H: 16000 }), { mode: 'conservative' });
chk('resizeModel conservative 继承 time', Math.abs(m3.time - m.time) < 1e-9);

// 8. plan
const probeRes = measureThroughput(presetById('cloud'), 0, 60, {});
chk('measureThroughput 字段', typeof probeRes.cellStepsPerMs === 'number' && probeRes.steps >= 3);
const plan = autoPlan(probeRes, presetById('cloud'), { frameMs: 16, simSpeed: 20, dt: 4 });
chk('autoPlan 字段', ['levelIndex', 'label', 'cells', 'maxCells', 'stepsPerFrame', 'estMsPerFrame', 'estMsPerStep', 'grid'].every((k) => k in plan));
chk('frameCost', frameCost({ cells: 8192, cellStepsPerMs: 2000, dt: 4, simSpeed: 20, frameMs: 16 }).msPerFrame > 0);
chk('presetById 回退', presetById('nope').id === 'cloud');
chk('预设 surface 默认', SCALE_PRESETS.find((x) => x.id === 'macro').surface.sensible === 250);

// 9. micro / state
const micro = createMicrophysics(g, { enable: false });
chk('micro disable 返回 0', micro.apply(m.S, m.base, 1) === 0);
chk('repairMoisture 返回结构', (() => { const res = repairMoisture(m.S, g.nx, g.nz); return 'repaired' in res && 'lost' in res; })());
chk('clampMoisture 返回数值', typeof clampMoisture(m.S) === 'number');
chk('maxAbs', maxAbs(new Float64Array([-2, 1])) === 2);

// 10. viz
const img = colorize(m, 'w', presetById('cloud').viz);
chk('colorize 结构', img.nx === g.nx && img.nz === g.nz && img.data.length === g.nx * g.nz * 4);
chk('colorize 未知 field 不抛错', colorize(m, 'nope').data.length === g.nx * g.nz * 4);
chk('FIELD_LIST 6 项 / DEFAULT_VIZ 5 键', FIELD_LIST.length === 6 && Object.keys(DEFAULT_VIZ).length === 5);
chk('DEFAULT_CONFIG evapTau', DEFAULT_CONFIG.micro.evapTau === 1200);

console.log(`通过 ${ok.length} 项`);
if (bad.length) {
  console.log(`不符 ${bad.length} 项：`);
  for (const b of bad) console.log('  ✗ ' + b);
  process.exit(1);
}
