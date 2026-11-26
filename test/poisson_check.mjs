import { Grid } from '../src/core/grid.js';
import { createPoissonSolver } from '../src/core/poisson.js';
import { makeBase } from '../src/model/base.js';

const g = new Grid({ nx: 16, nz: 8, Lx: 4000, H: 2000 });
const base = makeBase(g, {});
const solver = createPoissonSolver(g);
const { nx, nz, dx, dz } = g;

// 构造一个满足 w=0 边界、且总散度为零的任意速度场
const u = new Float64Array(nx * nz);
const w = new Float64Array(nx * (nz + 1));
for (let k = 0; k < nz; k++) {
  for (let i = 0; i < nx; i++) {
    u[i + nx * k] = Math.sin(2 * Math.PI * i / nx) * (1 + k);
  }
}
for (let k = 1; k < nz; k++) {
  for (let i = 0; i < nx; i++) {
    w[i + nx * k] = Math.cos(2 * Math.PI * i / nx) * (k / nz);
  }
}

const pi = new Float64Array(nx * nz);
const res = solver.project(u, w, 1.0, pi, base);
console.log('divBefore', res.divBefore.toExponential(3));
console.log('divAfter ', res.divAfter.toExponential(3));
console.log('max|pi|  ', Math.max(...Array.from(pi).map(Math.abs)).toExponential(3));

// 用离散拉普拉斯验证解：∇²π 是否等于 divBefore（逐点）
function lapAt(arr, i, k) {
  const ip = i + 1 === nx ? 0 : i + 1;
  const im = i === 0 ? nx - 1 : i - 1;
  let s = (arr[ip + nx * k] - 2 * arr[i + nx * k] + arr[im + nx * k]) / (dx * dx);
  const kp = k + 1 < nz ? k + 1 : k;
  const km = k > 0 ? k - 1 : k;
  if (k === 0) s += (arr[i + nx * 1] - arr[i]) / (dz * dz);
  else if (k === nz - 1) s += (arr[i + nx * (nz - 2)] - arr[i + nx * (nz - 1)]) / (dz * dz);
  else s += (arr[i + nx * kp] - 2 * arr[i + nx * k] + arr[i + nx * km]) / (dz * dz);
  return s;
}

// 重建投影前的散度
const divBefore = Float64Array.from(solver.divBeforeArr ?? []);
void divBefore;
console.log('检查：对 u,w 重新投影前先备份');

// 重新构造一次以取投影前的散度
const u2 = u.slice(), w2 = w.slice();
const solver2 = createPoissonSolver(g);
const divAll = new Float64Array(nx * nz);
solver2.computeDiv(u2, w2, base, divAll);

let maxErr = 0;
for (let k = 0; k < nz; k++) {
  for (let i = 0; i < nx; i++) {
    const n = i + nx * k;
    const lhs = lapAt(pi, i, k);
    const rhs = divAll[n] / 1.0;
    const err = Math.abs(lhs - rhs);
    if (err > maxErr) maxErr = err;
  }
}
console.log('max |∇²π - div/Δt| =', maxErr.toExponential(3));
console.log('max |div| =', Math.max(...Array.from(divAll).map(Math.abs)).toExponential(3));
