import { Grid } from '../src/core/grid.js';
import { createPoissonSolver } from '../src/core/poisson.js';
import { makeBase } from '../src/model/base.js';

const nx = 64, nz = 48;
const g = new Grid({ nx, nz, Lx: 32000, H: 12000 });
const base = makeBase(g, {});
const solver = createPoissonSolver(g);
const { dx, dz } = g;

let seed = 7;
const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff - 0.5; };

const u = new Float64Array(nx * nz);
const w = new Float64Array(nx * (nz + 1));
for (let n = 0; n < u.length; n++) u[n] = rnd() * 10;
for (let k = 1; k < nz; k++) for (let i = 0; i < nx; i++) w[i + nx * k] = rnd() * 5;

const div0 = new Float64Array(nx * nz);
solver.computeDiv(u, w, base, div0);

const pi = new Float64Array(nx * nz);
solver.project(u, w, 1.0, pi, base);

// 用严格的离散算子复算 L(π)
let maxErr = 0, maxL = 0;
for (let k = 0; k < nz; k++) {
  for (let i = 0; i < nx; i++) {
    const n = i + nx * k;
    const ip = i + 1 === nx ? 0 : i + 1;
    const im = i === 0 ? nx - 1 : i - 1;
    let L = (pi[ip + nx * k] - 2 * pi[n] + pi[im + nx * k]) / (dx * dx);

    if (k === 0) L += (pi[i + nx * 1] - pi[i]) / (dz * dz);
    else if (k === nz - 1) L += (pi[i + nx * (nz - 2)] - pi[i + nx * (nz - 1)]) / (dz * dz);
    else L += (pi[i + nx * (k + 1)] - 2 * pi[n] + pi[i + nx * (k - 1)]) / (dz * dz);

    const rhs = div0[n];
    const err = Math.abs(L - rhs);
    if (err > maxErr) maxErr = err;
    if (Math.abs(rhs) > maxL) maxL = Math.abs(rhs);
  }
}
console.log(`max|L(pi) - div| = ${maxErr.toExponential(3)}   max|div| = ${maxL.toExponential(3)}`);

// 逐列检查：误差最大的 k
const errByK = new Float64Array(nz);
for (let k = 0; k < nz; k++) {
  for (let i = 0; i < nx; i++) {
    const n = i + nx * k;
    const ip = i + 1 === nx ? 0 : i + 1;
    const im = i === 0 ? nx - 1 : i - 1;
    let L = (pi[ip + nx * k] - 2 * pi[n] + pi[im + nx * k]) / (dx * dx);
    if (k === 0) L += (pi[i + nx] - pi[i]) / (dz * dz);
    else if (k === nz - 1) L += (pi[i + nx * (nz - 2)] - pi[i + nx * (nz - 1)]) / (dz * dz);
    else L += (pi[i + nx * (k + 1)] - 2 * pi[n] + pi[i + nx * (k - 1)]) / (dz * dz);
    const e = Math.abs(L - div0[n]);
    if (e > errByK[k]) errByK[k] = e;
  }
}
console.log('误差按层（前 5 层与后 5 层）：',
  Array.from(errByK.slice(0, 5)).map((v) => v.toExponential(1)),
  '...',
  Array.from(errByK.slice(-5)).map((v) => v.toExponential(1)));
