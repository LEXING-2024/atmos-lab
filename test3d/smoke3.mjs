/** 三维湿对流冒烟：30 分钟演化、守恒、散度、峰值强度
 *
 * 判定口径：单泡闭域无持续强迫，30 分钟时已进入衰减期（maxW ≈ 2.7），
 * 用瞬时值判断强弱会误报。改为记录全程峰值，并以峰值、云顶、
 * 水量收支、散度四项一起判。
 */
import { Grid3 } from '../src/three/core/grid3.js';
import { AnelasticModel3 } from '../src/three/model/anelastic3.js';

const g = new Grid3({ nx: 32, ny: 32, nz: 24, Lx: 12000, Ly: 12000, H: 9000 });
const m = new AnelasticModel3(g, {
  base: { TSfc: 299.5, RHsfc: 0.78 },
  micro: { enable: true },
});
m.initAtmosphere({ amplitude: 3, radius: 2000, zCenter: 1200, noise: 0.02 });

const t0 = Date.now();
let peakW = 0;
let peakU = 0;
let peakTop = 0;
let steps = 0;
for (let i = 0; i < 10; i++) {
  const r = m.advance(180, 0);
  steps += r.steps;
  const d = m.diagnostics();
  peakW = Math.max(peakW, d.maxW);
  peakU = Math.max(peakU, d.maxU);
  peakTop = Math.max(peakTop, d.cloudTop);
  console.log(
    `${(d.time / 60).toFixed(0).padStart(3)} min  maxW ${d.maxW.toFixed(2).padStart(6)}  ` +
      `云顶 ${(d.cloudTop / 1000).toFixed(1)} km  雨强 ${d.precipRate.toFixed(2)} mm/h  ` +
      `漂移 ${(((d.totalWater - m.initMass) / m.initMass) * 100).toFixed(3)}%`,
  );
}
const d = m.diagnostics();
const used = Date.now() - t0;
const drift = (d.totalWater - m.initMass) / m.initMass;

console.log('');
console.log(`总计 ${(d.time / 60).toFixed(1)} min，${steps} 步，${(used / steps).toFixed(1)} ms/步`);
console.log(`峰值 maxW ${peakW.toFixed(2)} m/s  峰值 maxU ${peakU.toFixed(2)}  云顶 ${(peakTop / 1000).toFixed(1)} km`);
console.log(`累计降水 ${d.precipAccum.toFixed(3)} mm  水量漂移 ${(drift * 100).toFixed(3)}%  散度残差 ${d.divResidual.toExponential(2)}`);

const checks = [
  ['峰值 w > 10 m/s', peakW > 10],
  ['云顶 > 6 km', peakTop > 6000],
  ['水量收支闭合 (|漂移| < 1%)', Math.abs(drift) < 0.01],
  ['投影后散度 < 1e-12', d.divResidual < 1e-12],
  ['未出 NaN', Number.isFinite(d.maxW) && Number.isFinite(d.totalWater)],
];
for (const [name, ok] of checks) console.log(`${ok ? '✓' : '✗'} ${name}`);
console.log(checks.every((c) => c[1]) ? '判定: 通过' : '判定: 失败');
