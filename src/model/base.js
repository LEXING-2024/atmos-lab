/**
 * 基态（静力平衡、无扰动）。
 * 只依赖网格与参数，不依赖任何状态。
 */

import { CONST, qsat, exner } from '../physics/thermo.js';

export function makeBase(grid, opts = {}) {
  const g = CONST.g, Rd = CONST.Rd;

  const pSfc = opts.pSfc ?? 100000;
  const TSfc = opts.TSfc ?? 299.5;
  const lapse = opts.lapse ?? 0.007;
  const zTrop = opts.zTrop ?? 11000;
  const RHsfc = opts.RHsfc ?? 0.78;
  const RHmin = opts.RHmin ?? 0.03;
  const Hq = opts.Hq ?? 2500;

  const { nz, dz } = grid;

  const p0 = new Float64Array(nz);
  const T0 = new Float64Array(nz);
  const rho0 = new Float64Array(nz);
  const theta0 = new Float64Array(nz);
  const qv0 = new Float64Array(nz);
  const ex = new Float64Array(nz);
  const zc = new Float64Array(nz);

  for (let k = 0; k < nz; k++) {
    const z = (k + 0.5) * dz;
    zc[k] = z;

    let T, p;
    if (z < zTrop) {
      T = TSfc - lapse * z;
      p = pSfc * Math.pow(T / TSfc, g / (Rd * lapse));
    } else {
      T = TSfc - lapse * zTrop;
      const pT = pSfc * Math.pow(T / TSfc, g / (Rd * lapse));
      p = pT * Math.exp(-g * (z - zTrop) / (Rd * T));
    }

    T0[k] = T;
    p0[k] = p;
    rho0[k] = p / (Rd * T);
    ex[k] = exner(p);
    theta0[k] = T / ex[k];

    let rh = RHsfc * Math.exp(-z / Hq);
    if (rh < RHmin) rh = RHmin;
    if (rh > 1) rh = 1;
    qv0[k] = rh * qsat(T, p);
  }

  // 面密度：用于连续性方程与通量
  const rhoW = new Float64Array(nz + 1);
  rhoW[0] = rho0[0];
  for (let k = 1; k < nz; k++) rhoW[k] = 0.5 * (rho0[k - 1] + rho0[k]);
  rhoW[nz] = rho0[nz - 1];

  return { p0, T0, rho0, theta0, qv0, ex, zc, rhoW };
}
