/**
 * 压力投影求解器。
 *
 * 解  ∇²π = (1/Δt) ∇·(ρ₀ u*)
 * x 方向 FFT、z 方向三对角（上下 Neumann）。
 *
 * 关键：x 方向使用离散特征值 λm = 4 sin²(πm/nx)/Δx²，
 * 与离散散度/梯度算子严格一致，投影后剩余散度在机器精度量级。
 */

import { FFT1D } from './fft.js';

export function createPoissonSolver(grid) {
  const { nx, nz, dx, dz } = grid;
  const fft = new FFT1D(nx);
  const dz2 = dz * dz;
  const nCell = nx * nz;

  const specRe = new Float64Array(nCell);
  const specIm = new Float64Array(nCell);
  const rowRe = new Float64Array(nx);
  const rowIm = new Float64Array(nx);
  const cpRe = new Float64Array(nz);
  const cpIm = new Float64Array(nz);
  const dpRe = new Float64Array(nz);
  const dpIm = new Float64Array(nz);
  const xRe = new Float64Array(nz);
  const xIm = new Float64Array(nz);
  const div = new Float64Array(nCell);

  const lam = new Float64Array(nx);
  for (let m = 0; m < nx; m++) {
    const s = Math.sin(Math.PI * m / nx);
    lam[m] = 4 * s * s / (dx * dx);
  }

  /** 密度加权散度；返回最大绝对值 */
  function computeDiv(u, w, base, out) {
    const rho0 = base.rho0, rhoW = base.rhoW;
    let maxAbs = 0;
    for (let k = 0; k < nz; k++) {
      const rk = rho0[k];
      const off = nx * k;
      for (let i = 0; i < nx; i++) {
        const ip = i + 1 === nx ? 0 : i + 1;
        const d = (rk * (u[ip + off] - u[i + off])) / dx
          + (rhoW[k + 1] * w[i + off + nx] - rhoW[k] * w[i + off]) / dz;
        out[i + off] = d;
        const a = d < 0 ? -d : d;
        if (a > maxAbs) maxAbs = a;
      }
    }
    return maxAbs;
  }

  /**
   * @param {Float64Array} u
   * @param {Float64Array} w
   * @param {number} dt
   * @param {Float64Array} pi
   * @param {{rho0:Float64Array, rhoW:Float64Array}} base
   */
  function project(u, w, dt, pi, base) {
    const divBefore = computeDiv(u, w, base, div);

    for (let n = 0; n < nCell; n++) {
      specRe[n] = div[n] / dt;
      specIm[n] = 0;
    }

    for (let k = 0; k < nz; k++) {
      const off = nx * k;
      for (let i = 0; i < nx; i++) { rowRe[i] = specRe[off + i]; rowIm[i] = specIm[off + i]; }
      fft.transform(rowRe, rowIm, false);
      for (let i = 0; i < nx; i++) { specRe[off + i] = rowRe[i]; specIm[off + i] = rowIm[i]; }
    }

    const o1 = 1 / dz2;
    const dInterior = -2 * o1;
    const dBound = -1 * o1;

    for (let m = 0; m < nx; m++) {
      const l = lam[m];
      const isMean = m === 0;

      // k = 0：Neumann 面在 z=0，鬼点 π(-1)=π(0) → 对角 -1/dz²、上对角 1/dz²
      // m=0 时算子奇异，钉住 π(0)=0 消去零空间（兼容性条件由连续性保证）
      const d0 = isMean ? 1 : dBound - l;
      const c0 = isMean ? 0 : o1;
      cpRe[0] = c0 / d0;
      cpIm[0] = 0;
      dpRe[0] = (isMean ? 0 : specRe[m]) / d0;
      dpIm[0] = (isMean ? 0 : specIm[m]) / d0;

      for (let k = 1; k < nz; k++) {
        const isTop = k === nz - 1;
        const a = o1;
        const c = isTop ? 0 : o1;
        const diag = (isTop ? dBound : dInterior) - l;
        const sr = specRe[m + nx * k], si = specIm[m + nx * k];

        const dr = diag - a * cpRe[k - 1];
        const di = -a * cpIm[k - 1];
        const dmag = dr * dr + di * di;

        cpRe[k] = c * dr / dmag;
        cpIm[k] = -c * di / dmag;

        const nr = sr - a * dpRe[k - 1];
        const ni = si - a * dpIm[k - 1];
        dpRe[k] = (nr * dr + ni * di) / dmag;
        dpIm[k] = (ni * dr - nr * di) / dmag;
      }

      xRe[nz - 1] = dpRe[nz - 1];
      xIm[nz - 1] = dpIm[nz - 1];
      for (let k = nz - 2; k >= 0; k--) {
        const cr = cpRe[k], ci = cpIm[k];
        const xr = xRe[k + 1], xi = xIm[k + 1];
        xRe[k] = dpRe[k] - (cr * xr - ci * xi);
        xIm[k] = dpIm[k] - (cr * xi + ci * xr);
      }

      for (let k = 0; k < nz; k++) {
        specRe[m + nx * k] = xRe[k];
        specIm[m + nx * k] = xIm[k];
      }
    }

    for (let k = 0; k < nz; k++) {
      const off = nx * k;
      for (let i = 0; i < nx; i++) { rowRe[i] = specRe[off + i]; rowIm[i] = specIm[off + i]; }
      fft.transform(rowRe, rowIm, true);
      for (let i = 0; i < nx; i++) { pi[off + i] = rowRe[i]; }
    }

    const rho0 = base.rho0, rhoW = base.rhoW;

    for (let k = 0; k < nz; k++) {
      const rk = rho0[k], off = nx * k;
      for (let i = 0; i < nx; i++) {
        const im = i === 0 ? nx - 1 : i - 1;
        u[i + off] -= dt * (pi[i + off] - pi[im + off]) / (dx * rk);
      }
    }
    for (let k = 1; k < nz; k++) {
      const off = nx * k;
      for (let i = 0; i < nx; i++) {
        w[i + off] -= dt * (pi[i + off] - pi[i + off - nx]) / (dz * rhoW[k]);
      }
    }
    for (let i = 0; i < nx; i++) { w[i] = 0; w[i + nx * nz] = 0; }

    const divAfter = computeDiv(u, w, base, div);
    return { divBefore, divAfter };
  }

  return { fft, project, computeDiv };
}
