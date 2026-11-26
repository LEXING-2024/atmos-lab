/**
 * 湍流：Smagorinsky + 稳定度修正 + 地表通量。
 * 接口契约：apply(state, base, tendencies)  直接累加到倾向数组。
 */

import { CONST } from './thermo.js';

export function createTurbulence(grid, cfg) {
  const { nx, nz, dx, dz } = grid;
  const nCell = grid.nCell;
  const K = new Float64Array(nCell);
  const Kx = new Float64Array(nCell);   // +x 面
  const Kz = new Float64Array(nCell);   // +z 面

  const cs = cfg.cs ?? 0.15;
  // 与网格自适应：粗网格取配置值，细网格(K 的物理量级随 Δ² 下降)自动收紧
  const kMinCfg = cfg.kMin ?? 1.0;
  const kMaxCfg = cfg.kMax ?? 1500;
  const kMin = Math.min(kMinCfg, 0.01 * dz * dz);
  const kMax = Math.min(kMaxCfg, 0.05 * dz * dz);
  const cd = cfg.cd ?? 0.002;
  const L = Math.sqrt(dx * dz);
  const dx2 = dx * dx, dz2 = dz * dz;
  const g = CONST.g, cp = CONST.cp, Lv = CONST.Lv;

  function computeK(S, base) {
    const theta0 = base.theta0;
    for (let k = 0; k < nz; k++) {
      const off = nx * k;
      const kp = k + 1 < nz ? k + 1 : k;
      const km = k > 0 ? k - 1 : k;
      for (let i = 0; i < nx; i++) {
        const n = off + i;
        const ip = i + 1 === nx ? 0 : i + 1;
        const du_dx = (S.u[ip + off] - S.u[n]) / dx;
        const dw_dz = (S.w[n + nx] - S.w[n]) / dz;
        const du_dz = (S.u[i + nx * kp] - S.u[i + nx * km]) / (dz * (kp - km || 1));
        const dw_dx = (S.w[ip + off] - S.w[n]) / dx;
        const S2 = du_dx * du_dx + dw_dz * dw_dz + 0.5 * (du_dz + dw_dx) * (du_dz + dw_dx);

        let N2 = (g / theta0[k]) * (S.th[i + nx * kp] - S.th[i + nx * km]) / (dz * (kp - km || 1));
        if (N2 < 0) N2 = 0;

        const Sm = Math.sqrt(S2);
        let kv = cs * cs * L * L * Sm;
        if (N2 > 1e-8 && S2 > 1e-12) kv = kv / (1 + 10 * N2 / S2);
        if (kv < kMin) kv = kMin;
        if (kv > kMax) kv = kMax;
        K[n] = kv;
      }
    }

    // 面值预先算好：6 个标量/动量方向共用，避免每个场重复插值
    for (let k = 0; k < nz; k++) {
      const off = nx * k;
      const up = k + 1 < nz ? nx : 0;
      for (let i = 0; i < nx; i++) {
        const n = off + i;
        const ip = i + 1 === nx ? 0 : i + 1;
        Kx[n] = 0.5 * (K[n] + K[ip + off]);
        Kz[n] = up ? 0.5 * (K[n] + K[n + up]) : K[n];
      }
    }
  }

  /**
   * 二阶扩散，K 取格点中心值半插值。
   * fb 为基态剖面：只对扰动量扩散，避免对基态产生虚假倾向。
   */
  function diffuse(f, tend, fb) {
    for (let k = 0; k < nz; k++) {
      const off = nx * k;
      const up = k + 1 < nz ? nx : 0;
      const dn = k > 0 ? -nx : 0;
      const fbHere = fb ? fb[k] : 0;
      const fbUp = fb ? (k + 1 < nz ? fb[k + 1] : fbHere) : 0;
      const fbDn = fb ? (k > 0 ? fb[k - 1] : fbHere) : 0;
      for (let i = 0; i < nx; i++) {
        const n = off + i;
        const ip = i + 1 === nx ? 0 : i + 1;
        const im = i === 0 ? nx - 1 : i - 1;

        let lap = 0;
        if (up) lap += Kz[n] * ((f[n + up] - f[n]) - (fb ? fbUp - fbHere : 0)) / dz2;
        if (dn) lap -= Kz[n + dn] * ((f[n] - f[n + dn]) - (fb ? fbHere - fbDn : 0)) / dz2;
        lap += Kx[n] * (f[ip + off] - f[n]) / dx2;
        lap -= Kx[im + off] * (f[n] - f[im + off]) / dx2;

        tend[n] += lap;
      }
    }
  }

  // 地表通量的水平分布（与状态无关，构造时算一次）
  const sur = cfg.surface ?? {};
  const H = sur.sensible ?? 0;
  const LE = sur.latent ?? 0;
  const xc0 = sur.xCenter;
  const xw = sur.xWidth ?? grid.Lx;
  const facCol = new Float64Array(nx);
  for (let i = 0; i < nx; i++) {
    let fac = 1;
    if (xc0 !== undefined && xc0 !== null) {
      const d = (grid.xc(i) - xc0) / xw;
      fac = Math.exp(-d * d);
    }
    facCol[i] = fac;
  }
  // 潜热通量换来的水汽输入，单位 kg/(m·s)（每单位宽度）
  let waterRate = 0;
  for (let i = 0; i < nx; i++) waterRate += (LE / Lv) * facCol[i] * dx;

  function apply(S, base, T) {
    computeK(S, base);

    diffuse(S.u, T.tu, null);
    diffuse(S.th, T.tth, base.theta0);
    diffuse(S.qv, T.tqv, base.qv0);
    diffuse(S.qc, T.tqc, null);
    diffuse(S.qr, T.tqr, null);

    // w 只在中层，边界固定为 0（顶/底取真实零梯度对应的面值）
    for (let k = 1; k < nz; k++) {
      const off = nx * k;
      const km = -nx;
      for (let i = 0; i < nx; i++) {
        const n = off + i;
        const ip = i + 1 === nx ? 0 : i + 1;
        const im = i === 0 ? nx - 1 : i - 1;
        const kzp = Kz[n];
        const kzm = Kz[n + km];
        const kxp = Kx[n];
        const kxm = Kx[im + off];

        let lap = 0;
        lap += kzp * (S.w[n + nx] - S.w[n]) / dz2;
        lap -= kzm * (S.w[n] - S.w[n + km]) / dz2;
        lap += kxp * (S.w[ip + off] - S.w[n]) / dx2;
        lap -= kxm * (S.w[n] - S.w[im + off]) / dx2;

        T.tw[n] += lap;
      }
    }

    // 地表：拖曳 + 感热 + 潜热
    const rhoS = base.rho0[0];
    const tth = H / (rhoS * cp * dz);
    const tqv = LE / (rhoS * Lv * dz);

    for (let i = 0; i < nx; i++) {
      const n = i;
      const uS = S.u[n];
      T.tu[n] -= cd * Math.abs(uS) * uS / dz;

      const fac = facCol[i];
      if (fac > 1e-3) {
        T.tth[n] += tth * fac;
        T.tqv[n] += tqv * fac;
      }
    }
  }

  return { apply, computeK, K, waterRate };
}
