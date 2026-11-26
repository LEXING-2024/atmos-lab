import { Grid } from '../src/core/grid.js';
import { AnelasticModel } from '../src/model/anelastic.js';

const g = new Grid({ nx: 128, nz: 64, Lx: 32000, H: 16000 });

const cases = [
  ['A amp3 r4000 zc2000', { amplitude: 3, radius: 4000, zCenter: 2000 }, { filter: 0.02 }],
  ['B amp3 r5000 RH0.85', { amplitude: 3, radius: 5000, zCenter: 2000 }, { filter: 0.02, base: { RHsfc: 0.85 } }],
  ['C amp2.5 r4000 +SH100/LH150', { amplitude: 2.5, radius: 4000, zCenter: 2000 }, { filter: 0.02, turb: { surface: { sensible: 100, latent: 150, xWidth: 6000 } } }],
  ['D amp4 r3000 zc1500', { amplitude: 4, radius: 3000, zCenter: 1500 }, { filter: 0.02 }],
];

for (const [name, init, cfg] of cases) {
  const m = new AnelasticModel(g, cfg);
  m.initAtmosphere(init);
  let peakW = 0, peakQc = 0, peakTop = 0, precipMax = 0;
  const t0 = Date.now();
  for (let s = 0; s < 24; s++) {
    m.advance(300, 0);
    const d = m.diagnostics();
    if (d.maxW > peakW) peakW = d.maxW;
    if (d.cloudTop > peakTop) peakTop = d.cloudTop;
    if (d.precipRate > precipMax) precipMax = d.precipRate;
    for (let n = 0; n < m.S.qc.length; n++) if (m.S.qc[n] > peakQc) peakQc = m.S.qc[n];
  }
  const mm = m.precipAccum;
  console.log(`${name.padEnd(30)} 峰值 maxW=${peakW.toFixed(1).padStart(5)} m/s  峰值qc=${(peakQc * 1000).toFixed(2)} g/kg  最高云顶=${(peakTop / 1000).toFixed(1)} km  峰值降水=${precipMax.toFixed(1)} mm/h  累计=${mm.toFixed(2)} mm  用时${((Date.now() - t0) / 1000).toFixed(1)}s`);
}
