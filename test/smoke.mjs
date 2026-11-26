/**
 * 无头验证。在 Node 中运行，检查：
 *   散度残差、水量守恒、CFL、稳定性、云顶与最大垂直速度
 * 用法：node test/smoke.mjs
 */

import { Grid } from '../src/core/grid.js';
import { AnelasticModel, resizeModel } from '../src/model/anelastic.js';

const fmt = (v, d = 2) => (Number.isFinite(v) ? v.toFixed(d) : String(v));

function runCase(name, o) {
  const grid = new Grid({ nx: o.nx, nz: o.nz, Lx: o.Lx, H: o.H });
  const m = new AnelasticModel(grid, o.opts || {});
  m.initAtmosphere(o.init || {});

  const rows = [];
  const total = o.simSeconds;
  const chunk = o.chunk ?? 600;
  const t0 = Date.now();
  let worst = { cfl: 0, div: 0, nan: false, maxW: 0, maxTop: 0 };

  for (let t = 0; t < total; t += chunk) {
    const dt = o.dt ?? m.suggestDt();
    m.advance(Math.min(chunk, total - t), 0);
    const d = m.diagnostics();
    worst.cfl = Math.max(worst.cfl, d.cfl);
    worst.div = Math.max(worst.div, d.divResidual);
    worst.maxW = Math.max(worst.maxW, d.maxW);
    worst.maxTop = Math.max(worst.maxTop, d.cloudTop);
    if (!Number.isFinite(d.maxW) || !Number.isFinite(d.totalWater)) worst.nan = true;
    rows.push(d);
    void dt;
  }
  const ms = Date.now() - t0;

  const last = rows[rows.length - 1];
  const waterNow = last.totalWater + m.precipAccum * grid.Lx - m.waterInput;
  const drift = m.initMass > 0 ? (waterNow - m.initMass) / m.initMass : 0;

  console.log(`\n=== ${name} ===`);
  console.log(`网格 ${o.nx}×${o.nz}  Δx=${fmt(grid.dx, 1)}m Δz=${fmt(grid.dz, 1)}m  dt=${fmt(m.lastDt, 2)}s`);
  console.log(`耗时 ${ms} ms  |  ${fmt(ms / (total / m.lastDt), 2)} ms/步  |  ${fmt(total / m.lastDt, 0)} 步`);
  console.log(`maxW=${fmt(last.maxW)}  maxU=${fmt(last.maxU)}  云顶=${fmt(last.cloudTop / 1000, 1)}km`);
  console.log(`全程峰值：maxW=${fmt(worst.maxW)} m/s  最高云顶=${fmt(worst.maxTop / 1000, 1)} km`);
  console.log(`降水率=${fmt(last.precipRate, 2)} mm/h  累计=${fmt(last.precipAccum, 2)} mm`);
  console.log(`CFL(峰值)=${fmt(worst.cfl, 3)}  散度残差(峰值)=${worst.div.toExponential(2)}`);
  console.log(`水量守恒漂移=${(drift * 100).toFixed(3)}%   NaN=${worst.nan}`);

  const ok = Number.isFinite(last.maxW)
    && worst.div < 1e-8
    && worst.cfl < 0.8
    && Math.abs(drift) < 0.05
    && last.maxW < 80;
  console.log(ok ? '判定: 通过' : '判定: 失败');
  return { m, ok, last, drift, worst };
}

/* 1) 静力平衡：无扰动，1 小时应当保持静止 */
runCase('静力平衡（无扰动 1h）', {
  nx: 64, nz: 48, Lx: 32000, H: 12000,
  opts: { micro: { enable: false }, turb: { surface: { sensible: 0, latent: 0 } } },
  init: { amplitude: 0, noise: 0 },
  simSeconds: 3600, chunk: 900,
});

/* 2) 干动力：热泡上升，不含水汽相变 */
runCase('干热泡（无微物理 1h）', {
  nx: 128, nz: 64, Lx: 32000, H: 16000,
  opts: {
    micro: { enable: false },
    base: { TSfc: 299.5, RHsfc: 0.05 },
    turb: { surface: { sensible: 0, latent: 0 } },
  },
  init: { amplitude: 2.5, radius: 2000, zCenter: 1500, noise: 0.02 },
  simSeconds: 3600, chunk: 900,
});

/* 3) 湿对流：完整物理，出现深对流与降水 */
const wet = runCase('湿对流（完整物理 2h）', {
  nx: 128, nz: 64, Lx: 32000, H: 16000,
  opts: {
    base: { TSfc: 299.5, RHsfc: 0.78 },
    turb: { surface: { sensible: 0, latent: 0 } },
  },
  init: { amplitude: 3, radius: 4000, zCenter: 2000, noise: 0.02 },
  simSeconds: 7200, chunk: 900,
});

/* 4) 分辨率切换：守恒重映射应保持水量与稳定性 */
{
  const old = wet.m;
  const newGrid = new Grid({ nx: 256, nz: 96, Lx: 32000, H: 16000 });
  const m2 = resizeModel(old, newGrid, { mode: 'conservative' });
  const metric = (mm) => mm.totalWater() + mm.precipAccum * newGrid.Lx - mm.waterInput;
  const before = metric(m2);
  m2.advance(600, 0);
  const d = m2.diagnostics();
  const water = metric(m2);
  console.log('\n=== 重映射（128×64 → 256×96，再跑 10 分钟）===');
  console.log(`重映射后水量 ${fmt(before / 1e6, 4)} ×10⁶ kg/m   10 分钟后 ${fmt(water / 1e6, 4)} ×10⁶ kg/m   变化 ${((water / before - 1) * 100).toFixed(3)}%`);
  console.log(`maxW=${fmt(d.maxW)}  云顶=${fmt(d.cloudTop / 1000, 1)}km  散度残差=${d.divResidual.toExponential(2)}`);
  console.log(Number.isFinite(d.maxW) && d.divResidual < 1e-8 ? '判定: 通过' : '判定: 失败');
}

/* 5) 快照 / 恢复一致性 */
{
  const g = new Grid({ nx: 128, nz: 64, Lx: 32000, H: 16000 });
  const m = new AnelasticModel(g, {});
  m.initAtmosphere({ amplitude: 2.5, radius: 2000 });
  m.advance(600, 0);
  const snap = m.snapshot();
  m.advance(600, 0);
  const a = m.diagnostics().maxW;
  m.restore(snap);
  m.advance(600, 0);
  const b = m.diagnostics().maxW;
  console.log('\n=== 快照/恢复 ===');
  console.log(`maxW(连续) = ${fmt(a)}   maxW(恢复后重算) = ${fmt(b)}   差 ${fmt(Math.abs(a - b), 4)}`);
  console.log(Math.abs(a - b) < 1e-6 ? '判定: 通过（逐位复现）' : '判定: 失败');
}
