/**
 * 三维湍流：Smagorinsky（含三维应变张量）+ 稳定度修正 + 地表通量。
 * 接口契约：apply(state, base, tendencies) 直接累加到倾向数组。
 */

import { CONST } from '../../physics/thermo.js';

export function createTurbulence3(grid, cfg) {
  const { nx, ny, nz, dx, dy, dz } = grid;
  const nXY = nx * ny;
  const nCell = grid.nCell;

  const K = new Float64Array(nCell);
  const Kx = new Float64Array(nCell);
  const Ky = new Float64Array(nCell);
  const Kz = new Float64Array(nCell);

  const cs = cfg.cs ?? 0.15;
  const L = Math.cbrt(dx * dy * dz);
  const kMinCfg = cfg.kMin ?? 1.0;
  const kMaxCfg = cfg.kMax ?? 1500;
  const kMin = Math.min(kMinCfg, 0.01 * L * L);
  const kMax = Math.min(kMaxCfg, 0.05 * L * L);
  const cd = cfg.cd ?? 0.002;
  const dx2 = dx * dx, dy2 = dy * dy, dz2 = dz * dz;
  const g = CONST.g;

  function computeK(S, base) {
    const theta0 = base.theta0;
    for (let k = 0; k < nz; k++) {
      const off = nXY * k;
      const kp = k + 1 < nz ? k + 1 : k;
      const km = k > 0 ? k - 1 : k;
      const dzk = dz * (kp - km || 1);
      for (let j = 0; j < ny; j++) {
        const jp = j + 1 === ny ? 0 : j + 1;
        const jm = j === 0 ? ny - 1 : j - 1;
        const row = off + nx * j;
        const rowJp = off + nx * jp;
        const rowJm = off + nx * jm;
        for (let i = 0; i < nx; i++) {
          const ip = i + 1 === nx ? 0 : i + 1;
          const im = i === 0 ? nx - 1 : i - 1;
          const n = row + i;

          const du_dx = (S.u[row + ip] - S.u[n]) / dx;
          const dv_dy = (S.v[rowJp + i] - S.v[n]) / dy;
          const dw_dz = (S.w[n + nXY] - S.w[n]) / dz;
          const du_dy = (S.u[rowJp + i] - S.u[n]) / dy;
          const dv_dx = (S.v[row + ip] - S.v[n]) / dx;
          const du_dz = (S.u[i + nx * (j + ny * kp)] - S.u[i + nx * (j + ny * km)]) / dzk;
          const dv_dz = (S.v[i + nx * (j + ny * kp)] - S.v[i + nx * (j + ny * km)]) / dzk;
          const dw_dx = (S.w[row + ip] - S.w[n]) / dx;
          const dw_dy = (S.w[rowJp + i] - S.w[n]) / dy;

          const S2 = du_dx * du_dx + dv_dy * dv_dy + dw_dz * dw_dz
            + 0.5 * ((du_dy + dv_dx) * (du_dy + dv_dx)
              + (du_dz + dw_dx) * (du_dz + dw_dx)
              + (dv_dz + dw_dy) * (dv_dz + dw_dy));

          let N2 = (g / theta0[k]) * (S.th[i + nx * (j + ny * kp)] - S.th[i + nx * (j + ny * km)]) / dzk;
          if (N2 < 0) N2 = 0;

          const Sm = Math.sqrt(S2);
          let kv = cs * cs * L * L * Sm;
          if (N2 > 1e-8 && S2 > 1e-12) kv = kv / (1 + 10 * N2 / S2);
          if (kv < kMin) kv = kMin;
          if (kv > kMax) kv = kMax;
          K[n] = kv;
          void im;
        }
      }
    }

    // 面值
    for (let k = 0; k < nz; k++) {
      const off = nXY * k;
      const up = k + 1 < nz ? nXY : 0;
      for (let j = 0; j < ny; j++) {
        const jp = j + 1 === ny ? 0 : j + 1;
        const row = off + nx * j;
        const rowJp = off + nx * jp;
        for (let i = 0; i < nx; i++) {
          const ip = i + 1 === nx ? 0 : i + 1;
          const n = row + i;
          Kx[n] = 0.5 * (K[n] + K[row + ip]);
          Ky[n] = 0.5 * (K[n] + K[rowJp + i]);
          Kz[n] = up ? 0.5 * (K[n] + K[n + up]) : K[n];
        }
      }
    }
  }

  /** 二阶扩散；fb 为基态剖面（只扩散扰动，避免对基态产生虚假倾向） */
  function diffuse(f, tend, fb) {
    for (let k = 0; k < nz; k++) {
      const off = nXY * k;
      const fbHere = fb ? fb[k] : 0;
      const fbUp = fb ? (k + 1 < nz ? fb[k + 1] : fbHere) : 0;
      const fbDn = fb ? (k > 0 ? fb[k - 1] : fbHere) : 0;
      for (let j = 0; j < ny; j++) {
        const jp = j + 1 === ny ? 0 : j + 1;
        const jm = j === 0 ? ny - 1 : j - 1;
        const row = off + nx * j;
        const rowJp = off + nx * jp;
        const rowJm = off + nx * jm;
        for (let i = 0; i < nx; i++) {
          const ip = i + 1 === nx ? 0 : i + 1;
          const im = i === 0 ? nx - 1 : i - 1;
          const n = row + i;

          let lap = Kx[n] * (f[row + ip] - f[n]) / dx2;
          lap -= Kx[row + im] * (f[n] - f[row + im]) / dx2;
          lap += Ky[n] * (f[rowJp + i] - f[n]) / dy2;
          lap -= Ky[rowJm + i] * (f[n] - f[rowJm + i]) / dy2;

          if (k + 1 < nz) {
            const dUp = (f[n + nXY] - f[n]) - (fb ? fbUp - fbHere : 0);
            lap += Kz[n] * dUp / dz2;
          }
          if (k > 0) {
            const dDn = (f[n] - f[n - nXY]) - (fb ? fbHere - fbDn : 0);
            lap -= Kz[n - nXY] * dDn / dz2;
          }

          tend[n] += lap;
        }
      }
    }
  }

  // 地表通量的水平分布（构造时算一次）
  const sur = cfg.surface ?? {};
  const Hf = sur.sensible ?? 0;
  const LE = sur.latent ?? 0;
  const xc0 = sur.xCenter, yc0 = sur.yCenter;
  const xw = sur.xWidth ?? grid.Lx;
  const yw = sur.yWidth ?? grid.Ly;
  const facCol = new Float64Array(nXY);
  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      let fac = 1;
      if (xc0 !== undefined && xc0 !== null) {
        const dxr = (grid.xc(i) - xc0) / xw;
        const dyr = yc0 === undefined || yc0 === null ? 0 : (grid.yc(j) - yc0) / yw;
        fac = Math.exp(-(dxr * dxr + dyr * dyr));
      }
      facCol[i + nx * j] = fac;
    }
  }
  let waterRate = 0;
  for (let n = 0; n < nXY; n++) waterRate += (LE / CONST.Lv) * facCol[n] * dx * dy;

  function apply(S, base, T) {
    computeK(S, base);

    diffuse(S.u, T.tu, null);
    diffuse(S.v, T.tv, null);
    diffuse(S.th, T.tth, base.theta0);
    diffuse(S.qv, T.tqv, base.qv0);
    diffuse(S.qc, T.tqc, null);
    diffuse(S.qr, T.tqr, null);

    // w 只在中层；K 取上下两格的平均
    for (let k = 1; k < nz; k++) {
      const off = nXY * k, dn = -nXY;
      for (let j = 0; j < ny; j++) {
        const jp = j + 1 === ny ? 0 : j + 1;
        const jm = j === 0 ? ny - 1 : j - 1;
        const row = off + nx * j;
        for (let i = 0; i < nx; i++) {
          const ip = i + 1 === nx ? 0 : i + 1;
          const im = i === 0 ? nx - 1 : i - 1;
          const n = row + i;
          const kxw = 0.5 * (Kx[n] + Kx[n + dn]);
          const kyw = 0.5 * (Ky[n] + Ky[n + dn]);
          const kzw = Kz[n];
          const kzwn = Kz[n + dn];

          let lap = 0;
          lap += kxw * (S.w[row + ip] - S.w[n]) / dx2;
          lap -= 0.5 * (Kx[row + im] + Kx[row + im + dn]) * (S.w[n] - S.w[row + im]) / dx2;
          lap += kyw * (S.w[i + nx * (jp + ny * k)] - S.w[n]) / dy2;
          lap -= 0.5 * (Ky[i + nx * (jm + ny * k)] + Ky[i + nx * (jm + ny * (k - 1))]) * (S.w[n] - S.w[i + nx * (jm + ny * k)]) / dy2;
          if (k + 1 <= nz) lap += kzw * (S.w[n + nXY] - S.w[n]) / dz2;
          lap -= kzwn * (S.w[n] - S.w[n - nXY]) / dz2;

          T.tw[n] += lap;
        }
      }
    }

    // 地表：拖曳（u、v）+ 感热 + 潜热
    const rhoS = base.rho0[0];
    const tth = Hf / (rhoS * CONST.cp * dz);
    const tqv = LE / (rhoS * CONST.Lv * dz);

    for (let n = 0; n < nXY; n++) {
      const uS = S.u[n], vS = S.v[n];
      const spd = Math.sqrt(uS * uS + vS * vS);
      T.tu[n] -= cd * spd * uS / dz;
      T.tv[n] -= cd * spd * vS / dz;

      const fac = facCol[n];
      if (fac > 1e-3) {
        T.tth[n] += tth * fac;
        T.tqv[n] += tqv * fac;
      }
    }
  }

  return { apply, computeK, K, waterRate };
}
