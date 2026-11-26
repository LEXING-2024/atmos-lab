import { Grid } from '../src/core/grid.js';
import { AnelasticModel } from '../src/model/anelastic.js';

const g = new Grid({ nx: 128, nz: 64, Lx: 32000, H: 16000 });
const m = new AnelasticModel(g, {
  base: { TSfc: 299.5, RHsfc: 0.78 },
  turb: { surface: { sensible: 0, latent: 0 } },
});
m.initAtmosphere({ amplitude: 2.5, radius: 2000, zCenter: 1500, noise: 0.02 });

console.log('t(min)  maxW   cloudTop  maxQc(g/kg)  precip(mm/h)  water drift');
const init = m.initMass;

for (let t = 0; t < 120; t += 5) {
  m.advance(300, 0);
  let maxQc = 0;
  for (let n = 0; n < m.S.qc.length; n++) if (m.S.qc[n] > maxQc) maxQc = m.S.qc[n];
  const d = m.diagnostics();
  const water = d.totalWater + m.precipAccum * 1000 * g.Lx;
  const drift = (water - init) / init;
  console.log(
    `${(t + 5).toString().padStart(5)}  ${d.maxW.toFixed(2).padStart(6)}  ` +
    `${(d.cloudTop / 1000).toFixed(1).padStart(7)}km  ` +
    `${(maxQc * 1000).toFixed(2).padStart(9)}  ` +
    `${(d.precipRate * 3600 * 1000).toFixed(1).padStart(11)}  ` +
    `${(drift * 100).toFixed(2).padStart(8)}%  dt=${m.lastDt.toFixed(1)}`,
  );
}
