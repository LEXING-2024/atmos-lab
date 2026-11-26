/** 固壁投影求解器自检：散度消除精度、迭代次数、静力零流量 */
import { Grid3 } from '../src/three/core/grid3.js';
import { makeBase } from '../src/model/base.js';
import { makeTerrain, maskVelocities } from '../src/three/core/terrain3.js';
import { createPoissonSolver3IBM } from '../src/three/core/poisson3_ibm.js';

const g = new Grid3({ nx: 32, ny: 32, nz: 24, Lx: 12000, Ly: 12000, H: 9000 });
const base = makeBase(g, { TSfc: 299.5, RHsfc: 0.78 });
const terrain = makeTerrain(g, { kind: 'ridge', height: 1500, halfWidth: 2500 });
const solver = createPoissonSolver3IBM(g, terrain.mask);

console.log(`山脊掩膜：固壁 ${terrain.solidCount} / ${g.nCell} 格点（${(terrain.solidCount / g.nCell * 100).toFixed(1)}%）`);

function allocV() {
  const S = { u: new Float64Array(g.nCell), v: new Float64Array(g.nCell), w: new Float64Array(g.nW) };
  return S;
}

// 1) 静力：全零速度，投影后应仍为零
{
  const S = allocV();
  const res = solver.project(S.u, S.v, S.w, 6, new Float64Array(g.nCell), base);
  console.log(`静力：divBefore ${res.divBefore.toExponential(2)}  divAfter ${res.divAfter.toExponential(2)}  迭代 ${res.iterations}`);
}

// 2) 随机散度场：mask 掉固壁面后投影
{
  const S = allocV();
  let s = 12345;
  const rnd = () => ((s = (s * 1103515245 + 12345) >>> 0) / 4294967296 - 0.5);
  for (let n = 0; n < g.nCell; n++) S.u[n] = rnd() * 4;
  for (let n = 0; n < g.nCell; n++) S.v[n] = rnd() * 4;
  for (let n = 0; n < g.nW; n++) S.w[n] = rnd() * 2;
  maskVelocities(S, g, terrain);
  const res = solver.project(S.u, S.v, S.w, 6, new Float64Array(g.nCell), base);
  console.log(`随机：divBefore ${res.divBefore.toExponential(2)}  divAfter ${res.divAfter.toExponential(2)}  迭代 ${res.iterations}`);
  let inside = 0;
  for (let n = 0; n < g.nCell; n++) if (terrain.mask[n] && Math.abs(S.u[n]) + Math.abs(S.v[n]) > 0) inside++;
  console.log(`固壁内残余速度格点：${inside}`);
}
