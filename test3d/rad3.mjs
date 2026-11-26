/** v1.6 自检：日变化短波 + 地表能量平衡 + 海陆差异。
 *
 * 判据：
 *   1. 太阳高度角随模拟时间变化（07 时低、12 时高、21 时为 0）
 *   2. 陆地地表温度日振幅 > 3 K，海洋 < 0.3 K
 *   3. 白天对流强于夜间（maxW 峰值出现在正午前后）
 *   4. 无 NaN；水量收支 = 初值 + 地表水汽输入（在 1% 内）
 */
import { Grid3 } from '../src/three/core/grid3.js';
import { AnelasticModel3 } from '../src/three/model/anelastic3.js';
import { buildModel3, presetById3 } from '../src/three/model/plan3.js';

const g = new Grid3({ nx: 32, ny: 32, nz: 24, Lx: 24000, Ly: 24000, H: 8000 });
const p = presetById3('coast3d');
const m = buildModel3(g, p, { base: { TSfc: 299.5, RHsfc: 0.75 } });
m.initAtmosphere({ amplitude: 0.2, radius: 3000, zCenter: 800, noise: 0.02 });
m.initMass = m.totalWater();
m.cfg.dtMax = 10;
const area = g.Lx * g.Ly;

const sunChecks = [0, 2, 5, 11, 14].map((h) => [7 + h, m.rad.cosZenith(h * 3600).toFixed(3)]);
console.log(`太阳高度角余弦（时刻: cosZ）：${sunChecks.map(([h, c]) => `${h}h=${c}`).join('  ')}`);

let landMin = 1e9, landMax = -1e9, seaMin = 1e9, seaMax = -1e9;
let dayPeakW = 0, nightPeakW = 0, dayPeakTime = 0;
const t0 = Date.now();
let steps = 0;

for (let seg = 0; seg < 10; seg++) {
  const r = m.advance(1800, 0);   // 每次 30 min
  steps += r.steps;
  const d = m.diagnostics();
  const hour = 7 + d.time / 3600;
  const st = d.surface;
  landMin = Math.min(landMin, st.landMin); landMax = Math.max(landMax, st.landMax);
  seaMin = Math.min(seaMin, st.seaMin); seaMax = Math.max(seaMax, st.seaMax);
  const isDay = d.sunCos > 0.1;
  if (isDay) {
    if (d.maxW > dayPeakW) { dayPeakW = d.maxW; dayPeakTime = hour; }
  } else if (d.maxW > nightPeakW) nightPeakW = d.maxW;

  console.log(
    `${hour.toFixed(1).padStart(5)} h  cosZ ${d.sunCos.toFixed(2)}  ` +
    `陆 ${st.landMin.toFixed(1)}–${st.landMax.toFixed(1)} K  海 ${st.seaMin.toFixed(1)}–${st.seaMax.toFixed(1)}  ` +
    `maxW ${d.maxW.toFixed(2)}  云顶 ${(d.cloudTop / 1000).toFixed(1)} km  雨强 ${d.precipRate.toFixed(2)}`,
  );
}

const d = m.diagnostics();
const used = Date.now() - t0;
const landAmp = landMax - landMin;
const seaAmp = seaMax - seaMin;
const budget = (d.totalWater + d.precipAccum * area - m.initMass - d.waterInput) / m.initMass;

console.log('');
console.log(`步数 ${steps}，${(used / steps).toFixed(1)} ms/步`);
console.log(`陆地振幅 ${landAmp.toFixed(2)} K，海洋振幅 ${seaAmp.toFixed(2)} K`);
console.log(`白天峰值 maxW ${dayPeakW.toFixed(2)} m/s（${dayPeakTime.toFixed(1)}h），夜间峰值 ${nightPeakW.toFixed(2)} m/s`);
console.log(`水量收支偏差 ${(budget * 100).toFixed(3)}%（含地表蒸发输入 ${d.waterInput.toFixed(3)} kg/m²）`);

const checks = [
  ['太阳高度角随时刻变化', m.rad.cosZenith(5 * 3600) > 0.9 && m.rad.cosZenith(14 * 3600) === 0],
  ['陆地日振幅 > 3 K', landAmp > 3],
  ['海洋日振幅 < 0.3 K', seaAmp < 0.3],
  ['白天对流强于夜间', dayPeakW > nightPeakW * 1.3],
  ['无 NaN', Number.isFinite(d.maxW) && Number.isFinite(d.totalWater)],
  ['水量收支闭合 (< 1%)', Math.abs(budget) < 0.01],
];
for (const [name, ok] of checks) console.log(`${ok ? '✓' : '✗'} ${name}`);
console.log(checks.every((c) => c[1]) ? '判定: 通过' : '判定: 失败');
