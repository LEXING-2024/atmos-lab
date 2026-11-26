/** 三维投影精度：随机散度场 → 投影 → 剩余散度 + 严格离散 Laplacian 复算 */
import { Grid3 } from '../src/three/core/grid3.js';
import { createPoissonSolver3 } from '../src/three/core/poisson3.js';
import { makeBase } from '../src/model/base.js';

const g = new Grid3({ nx: 32, ny: 32, nz: 24, Lx: 12000, Ly: 12000, H: 9000 });
const base = makeBase(g, {});
const solver = createPoissonSolver3(g);
const { nx, ny, nz, dx, dy, dz } = g;
const nXY = nx * ny;

let seed = 11;
const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff - 0.5; };

const u = new Float64Array(g.nCell);
const v = new Float64Array(g.nCell);
const w = new Float64Array(g.nW);
for (let n = 0; n < u.length; n++) { u[n] = rnd() * 8; v[n] = rnd() * 8; }
for (let k = 1; k < nz; k++) for (let n = 0; n < nXY; n++) w[k * nXY + n] = rnd() * 4;

const div0 = new Float64Array(g.nCell);
solver.computeDiv(u, v, w, base, div0);

const pi = new Float64Array(g.nCell);
const res = solver.project(u, v, w, 1.0, pi, base);
console.log(`divBefore ${res.divBefore.toExponential(3)}  divAfter ${res.divAfter.toExponential(3)}`);

// 严格离散复算 δ²π − div
let maxErr = 0;
for (let k = 0; k < nz; k++) {
  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      const ip = i + 1 === nx ? 0 : i + 1;
      const im = i === 0 ? nx - 1 : i - 1;
      const jp = j + 1 === ny ? 0 : j + 1;
      const jm = j === 0 ? ny - 1 : j - 1;
      const n = i + nx * (j + ny * k);
      const nIp = ip + nx * (j + ny * k);
      const nIm = im + nx * (j + ny * k);
      const nJp = i + nx * (jp + ny * k);
      const nJm = i + nx * (jm + ny * k);

      let L = (pi[nIp] - 2 * pi[n] + pi[nIm]) / (dx * dx)
        + (pi[nJp] - 2 * pi[n] + pi[nJm]) / (dy * dy);

      if (k === 0) L += (pi[n + nXY] - pi[n]) / (dz * dz);
      else if (k === nz - 1) L += (pi[n - nXY] - pi[n]) / (dz * dz);
      else L += (pi[n + nXY] - 2 * pi[n] + pi[n - nXY]) / (dz * dz);

      const err = Math.abs(L - div0[n]);
      if (err > maxErr) maxErr = err;
    }
  }
}
console.log(`max|L(pi) - div| = ${maxErr.toExponential(3)}`);
console.log(maxErr < 1e-12 ? '判定: 通过' : '判定: 失败');
