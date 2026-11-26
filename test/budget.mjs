import { Grid } from '../src/core/grid.js';
import { AnelasticModel } from '../src/model/anelastic.js';

const g = new Grid({ nx: 128, nz: 64, Lx: 32000, H: 16000 });

const cases = [
  ['无地表通量 amp3 r4000', { amplitude: 3, radius: 4000, zCenter: 2000 }, {}],
  ['有感热无潜热 amp3 r4000', { amplitude: 3, radius: 4000, zCenter: 2000 }, { turb: { surface: { sensible: 100, latent: 0 } } }],
  ['仅潜热 amp3 r4000', { amplitude: 3, radius: 4000, zCenter: 2000 }, { turb: { surface: { sensible: 0, latent: 150, xWidth: 6000 } } }],
];

for (const [name, init, cfg] of cases) {
  const m = new AnelasticModel(g, cfg);
  m.initAtmosphere(init);
  const w0 = m.totalWater();
  m.advance(3600, 0);
  const d = m.diagnostics();
  const loss = m.precipAccum * g.Lx;
  console.log(`${name.padEnd(24)} 初始 ${(w0 / 1e6).toFixed(4)}  末 ${(d.totalWater / 1e6).toFixed(4)}  降水 ${(loss / 1e6).toFixed(4)}  地表输入 ${(m.waterInput / 1e6).toFixed(4)}  残差 ${(((d.totalWater + loss - m.waterInput) / w0 - 1) * 100).toFixed(3)}%  云顶 ${(d.cloudTop / 1000).toFixed(1)}km`);
}
