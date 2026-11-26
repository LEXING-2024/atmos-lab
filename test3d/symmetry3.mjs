/** 诊断：三维场在 x 与 y 方向上的对称性 */
import { Grid3 } from '../src/three/core/grid3.js';
import { AnelasticModel3 } from '../src/three/model/anelastic3.js';

const g = new Grid3({ nx: 32, ny: 32, nz: 24, Lx: 12000, Ly: 12000, H: 9000 });
const m = new AnelasticModel3(g, { base: { TSfc: 299.5, RHsfc: 0.78 }, micro: { enable: true } });
m.initAtmosphere({ amplitude: 3, radius: 2000, zCenter: 1200, noise: 0.02 });
m.advance(1080, 0);

const S = m.S, d = m.diagnostics();
let maxv = 0, maxu = 0, maxw = 0;
for (let n = 0; n < g.nCell; n++) {
  maxu = Math.max(maxu, Math.abs(S.u[n]));
  maxv = Math.max(maxv, Math.abs(S.v[n]));
  maxw = Math.max(maxw, Math.abs(S.w[n]));
}

// x–z 切面（j = ny/2）与 y–z 切面（i = nx/2）：逐格比较 qc
const jc = g.ny >> 1, ic = g.nx >> 1;
let same = 0, diff = 0, dmax = 0;
for (let k = 0; k < g.nz; k++) {
  for (let x = 0; x < g.nx; x++) {
    const a = S.qc[x + g.nx * (jc + g.ny * k)];
    const b = S.qc[ic + g.nx * (x + g.ny * k)];
    const dd = Math.abs(a - b);
    dmax = Math.max(dmax, dd);
    if (dd < 1e-6) same++; else diff++;
  }
}
console.log(`t=${(d.time / 60).toFixed(0)} min  maxU(纬向 u)=${maxu.toFixed(2)}  maxV(经向 v)=${maxv.toFixed(2)}  maxW=${maxw.toFixed(2)}`);
console.log(`两面 qc 相同格点 ${same}，不同格点 ${diff}，最大差 ${dmax.toExponential(2)} kg/kg`);
console.log(`域中心对称性：qc[i,j] vs qc[j,i] 的最大偏差同上（若接近噪声量级即为近轴对称）`);
