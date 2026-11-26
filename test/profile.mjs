import { Grid } from '../src/core/grid.js';
import { AnelasticModel } from '../src/model/anelastic.js';

const g = new Grid({ nx: 128, nz: 64, Lx: 32000, H: 16000 });
const m = new AnelasticModel(g, {});
m.initAtmosphere({ amplitude: 3, radius: 4000, zCenter: 2000 });
for (let i = 0; i < 20; i++) m.step(4);

const n = 60;
const time = (label, fn) => {
  const t0 = performance.now();
  for (let i = 0; i < n; i++) fn();
  const ms = (performance.now() - t0) / n;
  console.log(`${label.padEnd(24)} ${ms.toFixed(3)} ms  (${(g.nCell / ms / 1000).toFixed(0)}k 格点/ms)`);
};

time('整步 step()', () => m.step(4));
time('仅 computeTendencies', () => m.computeTendencies(m.S, m.T));
time('仅 project', () => m.solver.project(m.S.u, m.S.w, 4, m.pi, m.base));
time('仅 filter+project', () => m.applyFilter(0.02));
time('仅 micro', () => m.micro.apply(m.S, m.base, 4));
