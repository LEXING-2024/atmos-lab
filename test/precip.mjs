import { Grid } from '../src/core/grid.js';
import { AnelasticModel } from '../src/model/anelastic.js';

const g = new Grid({ nx: 128, nz: 64, Lx: 32000, H: 16000 });
const m = new AnelasticModel(g, { turb: { surface: { sensible: 100, latent: 150, xWidth: 6000 } } });
m.initAtmosphere({ amplitude: 2.5, radius: 4000, zCenter: 1500, noise: 0.02 });

const w0 = m.totalWater();
console.log(' t(min) maxW  云顶km  地表雨强mm/h  累计mm  总水(x1e6)  守恒偏差');
for (let s = 1; s <= 24; s++) {
  m.advance(300, 0);
  const d = m.diagnostics();
  const loss = m.precipAccum * g.Lx;
  const drift = (d.totalWater + loss - m.waterInput - w0) / w0 * 100;
  console.log(`  ${String(s * 5).padStart(3)}  ${d.maxW.toFixed(1).padStart(5)}  ${(d.cloudTop / 1000).toFixed(1).padStart(5)}  ${d.precipRate.toFixed(1).padStart(11)}  ${m.precipAccum.toFixed(2).padStart(6)}  ${(d.totalWater / 1e6).toFixed(3)}  ${drift.toFixed(3)}%`);
}
