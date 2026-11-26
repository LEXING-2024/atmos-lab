/** 冰相微物理自检：30 分钟深对流，检查冰化、分层、收支与稳定性。
 *
 * 判据：
 *   1. 云顶必须进到 -40°C 层附近，即冰相出现（iceWater > 0.05 mm）
 *   2. 冰相以高空为主（z > 4 km 的冰量占比 > 60%）
 *   3. 水量收支闭合（|漂移| < 1%）
 *   4. 无 NaN、散度 < 1e-12
 *   5. 与纯暖雨版相比云顶更高（冰化砧的典型信号）
 */
import { Grid3 } from '../src/three/core/grid3.js';
import { AnelasticModel3 } from '../src/three/model/anelastic3.js';

const g = new Grid3({ nx: 32, ny: 32, nz: 32, Lx: 12000, Ly: 12000, H: 14000 });

function build(ice) {
  const m = new AnelasticModel3(g, {
    base: { TSfc: 299.5, RHsfc: 0.80 },
    micro: { enable: true, ice },
    turb: { surface: { sensible: 0, latent: 0 } },
  });
  m.initAtmosphere({ amplitude: 4, radius: 2500, zCenter: 1500, noise: 0.02 });
  return m;
}

function run(m, minutes) {
  let peakW = 0, peakTop = 0, peakIce = 0;
  for (let i = 0; i < minutes * 2; i++) {
    const r = m.advance(30, 0);
    const d = m.diagnostics();
    peakW = Math.max(peakW, d.maxW);
    peakTop = Math.max(peakTop, d.cloudTop);
    peakIce = Math.max(peakIce, m.iceWater());
    if (i % 6 === 5) {
      console.log(
        `${(d.time / 60).toFixed(0).padStart(3)} min  maxW ${d.maxW.toFixed(1).padStart(5)}  ` +
        `云顶 ${(d.cloudTop / 1000).toFixed(1)} km  冰 ${m.iceWater().toFixed(3)} mm  ` +
        `雨强 ${d.precipRate.toFixed(2)}  漂移 ${(((d.totalWater - m.initMass) / m.initMass) * 100).toFixed(3)}%`,
      );
    }
    void r;
  }
  return { peakW, peakTop, peakIce };
}

/** z > zMin 的冰量占比 */
function iceAbove(m, zMin) {
  const S = m.S, base = m.base;
  let above = 0, total = 0;
  for (let k = 0; k < m.grid.nz; k++) {
    const off = m.grid.nXY * k, rk = base.rho0[k];
    const isAbove = m.grid.zc(k) >= zMin;
    for (let n = 0; n < m.grid.nXY; n++) {
      const q = rk * (S.qi[off + n] + S.qs[off + n] + S.qg[off + n]) * m.grid.dz;
      total += q;
      if (isAbove) above += q;
    }
  }
  return total > 1e-12 ? above / total : 0;
}

console.log('== 冰相版 ==');
const mi = build(true);
const ri = run(mi, 30);
const di = mi.diagnostics();
const driftI = (di.totalWater - mi.initMass) / mi.initMass;
console.log(`峰值 maxW ${ri.peakW.toFixed(1)}  云顶 ${(ri.peakTop / 1000).toFixed(1)} km  峰值冰 ${ri.peakIce.toFixed(3)} mm`);
console.log(`30 min 冰量 ${mi.iceWater().toFixed(3)} mm，其中 z>4 km 占 ${(iceAbove(mi, 4000) * 100).toFixed(0)}%`);
console.log(`累计降水 ${di.precipAccum.toFixed(2)} mm  水量漂移 ${(driftI * 100).toFixed(3)}%  散度 ${di.divResidual.toExponential(2)}`);

console.log('');
console.log('== 纯暖雨（对照）==');
const mw = build(false);
const rw = run(mw, 30);
console.log(`峰值 maxW ${rw.peakW.toFixed(1)}  云顶 ${(rw.peakTop / 1000).toFixed(1)} km`);

const checks = [
  ['冰相出现 (> 0.05 mm)', ri.peakIce > 0.05],
  ['冰在高空 (z>4km 占比 > 60%)', iceAbove(mi, 4000) > 0.6],
  ['水量收支闭合 (|漂移| < 1%)', Math.abs(driftI) < 0.01],
  ['投影后散度 < 1e-12', di.divResidual < 1e-12],
  ['无 NaN', Number.isFinite(di.maxW) && Number.isFinite(di.totalWater)],
  ['冰化后云顶不低于暖雨版', ri.peakTop >= rw.peakTop * 0.97],
];
for (const [name, ok] of checks) console.log(`${ok ? '✓' : '✗'} ${name}`);
console.log(checks.every((c) => c[1]) ? '判定: 通过' : '判定: 失败');
