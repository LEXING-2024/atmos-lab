import { Grid } from '../src/core/grid.js';
import { createPoissonSolver } from '../src/core/poisson.js';
import { makeBase } from '../src/model/base.js';

for (const [nx, nz] of [[64, 48], [128, 64], [64, 32]]) {
  const g = new Grid({ nx, nz, Lx: 32000, H: 12000 });
  const base = makeBase(g, {});
  const solver = createPoissonSolver(g);

  let seed = 7;
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff - 0.5; };

  const u = new Float64Array(nx * nz);
  const w = new Float64Array(nx * (nz + 1));
  for (let n = 0; n < u.length; n++) u[n] = rnd() * 10;
  for (let k = 1; k < nz; k++) for (let i = 0; i < nx; i++) w[i + nx * k] = rnd() * 5;
  for (let i = 0; i < nx; i++) { w[i] = 0; w[i + nx * nz] = 0; }

  const pi = new Float64Array(nx * nz);
  const res = solver.project(u, w, 1.0, pi, base);
  console.log(`nx=${nx} nz=${nz}  divBefore=${res.divBefore.toExponential(2)}  divAfter=${res.divAfter.toExponential(2)}`);
}
