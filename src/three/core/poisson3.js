/**
 * 三维压力投影求解器。
 *
 * 解  δ²π = (1/Δt) ∇·(ρ₀ u*)
 * 水平两向 FFT（周期）、z 向三对角（上下 Neumann）。
 *
 * 与二维版同一套离散：水平特征值 λ = 4sin²(πi/nx)/Δx² + 4sin²(πj/ny)/Δy²，
 * 与离散散度/梯度严格一致（ρ 在散度与梯度中相互抵消），投影后剩余散度在机器精度量级。
 * 唯一奇异模 (i,j)=(0,0) 钉住 π(0)=0；其相容性由周期边界 + 上下 w=0 自动满足。
 */

import { createHorizontalFFT } from './fft3.js';

export function createPoissonSolver3(grid) {
  const { nx, ny, nz, dx, dy, dz } = grid;
  const nXY = nx * ny;
  const nCell = grid.nCell;

  const fftH = createHorizontalFFT(grid);
  const dz2 = dz * dz;
  const o1 = 1 / dz2;

  const specRe = new Float64Array(nCell);
  const specIm = new Float64Array(nCell);
  const modeRe = new Float64Array(nCell);   // (mode, k) 布局
  const modeIm = new Float64Array(nCell);
  const div = new Float64Array(nCell);

  const cpRe = new Float64Array(nz);
  const cpIm = new Float64Array(nz);
  const dpRe = new Float64Array(nz);
  const dpIm = new Float64Array(nz);
  const xRe = new Float64Array(nz);
  const xIm = new Float64Array(nz);

  const lamX = new Float64Array(nx);
  const lamY = new Float64Array(ny);
  for (let i = 0; i < nx; i++) { const s = Math.sin(Math.PI * i / nx); lamX[i] = 4 * s * s / (dx * dx); }
  for (let j = 0; j < ny; j++) { const s = Math.sin(Math.PI * j / ny); lamY[j] = 4 * s * s / (dy * dy); }

  /** 密度加权散度；返回最大绝对值 */
  function computeDiv(u, v, w, base, out) {
    const rho0 = base.rho0, rhoW = base.rhoW;
    let maxAbs = 0;
    for (let k = 0; k < nz; k++) {
      const rk = rho0[k];
      const off = nXY * k;
      for (let j = 0; j < ny; j++) {
        const jp = j + 1 === ny ? 0 : j + 1;
        const rowOff = off + nx * j;
        const rowOut = off + nx * j;
        for (let i = 0; i < nx; i++) {
          const ip = i + 1 === nx ? 0 : i + 1;
          const n = rowOff + i;
          const d = rk * ((u[rowOff + ip] - u[n]) / dx + (v[nx * jp + off + i] - v[n]) / dy)
            + (rhoW[k + 1] * w[n + nXY] - rhoW[k] * w[n]) / dz;
          out[rowOut + i] = d;
          const a = d < 0 ? -d : d;
          if (a > maxAbs) maxAbs = a;
        }
      }
    }
    return maxAbs;
  }

  /**
   * @param {Float64Array} u nCell
   * @param {Float64Array} v nCell
   * @param {Float64Array} w nW
   * @param {number} dt
   * @param {Float64Array} pi nCell（输出）
   */
  function project(u, v, w, dt, pi, base) {
    const divBefore = computeDiv(u, v, w, base, div);

    for (let n = 0; n < nCell; n++) { specRe[n] = div[n] / dt; specIm[n] = 0; }
    fftH.forward(specRe, specIm);

    // 转成 (mode,k) 布局：水平 modes 连续，z 步长 1
    for (let k = 0; k < nz; k++) {
      const off = nXY * k;
      for (let j = 0; j < ny; j++) {
        for (let i = 0; i < nx; i++) {
          const mode = j * nx + i;
          modeRe[mode * nz + k] = specRe[off + nx * j + i];
          modeIm[mode * nz + k] = specIm[off + nx * j + i];
        }
      }
    }

    for (let j = 0; j < ny; j++) {
      for (let i = 0; i < nx; i++) {
        const mode = j * nx + i;
        const bM = mode * nz;
        const l = lamX[i] + lamY[j];
        const isMean = i === 0 && j === 0;

        // k=0：Neumann 面在 z=0，鬼点 π(-1)=π(0) → 对角 -(o1+l)
        const d0 = isMean ? 1 : -(o1 + l);   // isMean: 钉住 π(0)=0
        const c0 = isMean ? 0 : o1;
        cpRe[0] = c0 / d0;
        cpIm[0] = 0;
        dpRe[0] = (isMean ? 0 : modeRe[bM]) / d0;
        dpIm[0] = (isMean ? 0 : modeIm[bM]) / d0;

        for (let k = 1; k < nz; k++) {
          const isTop = k === nz - 1;
          const a = o1;
          const c = isTop ? 0 : o1;
          const diag = -((isTop ? o1 : 2 * o1) + l);   // b_k = -(δ²对角) - λ
          const sr = modeRe[bM + k], si = modeIm[bM + k];

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
          modeRe[bM + k] = xRe[k];
          modeIm[bM + k] = xIm[k];
        }
      }
    }

    for (let k = 0; k < nz; k++) {
      const off = nXY * k;
      for (let j = 0; j < ny; j++) {
        for (let i = 0; i < nx; i++) {
          const mode = j * nx + i;
          specRe[off + nx * j + i] = modeRe[mode * nz + k];
          specIm[off + nx * j + i] = modeIm[mode * nz + k];
        }
      }
    }
    fftH.inverse(specRe, specIm);

    for (let n = 0; n < nCell; n++) pi[n] = specRe[n];

    const rho0 = base.rho0, rhoW = base.rhoW;

    for (let k = 0; k < nz; k++) {
      const rk = rho0[k], off = nXY * k;
      for (let j = 0; j < ny; j++) {
        const rowOff = off + nx * j;
        for (let i = 0; i < nx; i++) {
          const im = i === 0 ? nx - 1 : i - 1;
          u[rowOff + i] -= dt * (pi[rowOff + i] - pi[rowOff + im]) / (dx * rk);
        }
        for (let i = 0; i < nx; i++) {
          const jm = j === 0 ? ny - 1 : j - 1;
          v[rowOff + i] -= dt * (pi[rowOff + i] - pi[off + nx * jm + i]) / (dy * rk);
        }
      }
    }
    for (let k = 1; k < nz; k++) {
      const off = nXY * k;
      for (let n = 0; n < nXY; n++) {
        w[off + n] -= dt * (pi[off + n] - pi[off - nXY + n]) / (dz * rhoW[k]);
      }
    }
    for (let n = 0; n < nXY; n++) { w[n] = 0; w[nXY * nz + n] = 0; }

    const divAfter = computeDiv(u, v, w, base, div);
    return { divBefore, divAfter };
  }

  return { fftH, project, computeDiv };
}
