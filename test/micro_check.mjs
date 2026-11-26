import { Grid } from '../src/core/grid.js';
import { makeBase } from '../src/model/base.js';
import { createMicrophysics } from '../src/physics/micro.js';

const g = new Grid({ nx: 32, nz: 48, Lx: 8000, H: 12000 });
const base = makeBase(g, {});
const micro = createMicrophysics(g, { k1: 1e-3, qc0: 0.6e-3, k2: 2.2, evapTau: 1200 });

const S = {};
S.th = new Float64Array(g.nCell); S.qv = new Float64Array(g.nCell);
S.qc = new Float64Array(g.nCell); S.qr = new Float64Array(g.nCell);
for (let k = 0; k < g.nz; k++) {
  for (let i = 0; i < g.nx; i++) {
    S.th[i + g.nx * k] = base.theta0[k];
    S.qv[i + g.nx * k] = base.qv0[k] * 0.5;
  }
}
// 在 2–5 km 放一层雨
for (let k = 8; k < 20; k++) {
  for (let i = 0; i < g.nx; i++) S.qr[i + g.nx * k] = 1e-3;
}

const total = () => {
  let s = 0;
  for (let k = 0; k < g.nz; k++) for (let i = 0; i < g.nx; i++) {
    s += base.rho0[k] * (S.qv[i + g.nx * k] + S.qc[i + g.nx * k] + S.qr[i + g.nx * k]);
  }
  return s * g.dx * g.dz;
};

const dt = 5;
let t0 = total();
let cumFlux = 0;
const steps = 200;
console.log('步  总水(×1e3 kg/m)  地表累计(×1e3)  两者和  残差');
for (let s = 1; s <= steps; s++) {
  const flux = micro.apply(S, base, dt);
  cumFlux += flux * dt * g.Lx / 1e3;
  if (s % 20 === 0) {
    const t = total() / 1e3;
    console.log(`${String(s).padStart(3)}  ${t.toFixed(3).padStart(12)}  ${cumFlux.toFixed(3).padStart(12)}  ${(t + cumFlux).toFixed(3)}  ${((t + cumFlux) / (t0 / 1e3) - 1).toExponential(2)}`);
  }
}
