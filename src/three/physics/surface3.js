/**
 * 地表能量平衡 + 海陆分布（v1.6）。
 *
 * 每个水平列一个地表温度 Tg（板式模型）：
 *   C_land · dTg/dt = SW_abs + LW_net − SH − LH
 *   陆地取 0.5 m 土壤热容，海洋取 50 m 混合层并回拉固定 SST。
 *
 * 感热/潜热用整体输运公式，作用在每列最低的流体层（与 terrain.surfaceK 一致）：
 *   SH = ρ cp Cd |V| (Tg − T_air)
 *   LH = ρ Lv Cd |V| β (q_sat(Tg) − q_air)
 *
 * 只做"地表给大气通量"这一件事，不含云、不含辐射（短波由 rad3 算好后传进来）。
 */

import { CONST, qsat } from '../../physics/thermo.js';

const SIGMA = 5.670374419e-8;

export function createSurface3(grid, cfg = {}, terrain = null) {
  const { nx, ny } = grid;
  const nXY = nx * ny;
  const mask = terrain ? terrain.mask : null;

  const kind = cfg.kind ?? 'land';
  const sst = cfg.sst ?? 301.0;
  const landAlbedo = cfg.albedoLand ?? 0.18;
  const seaAlbedo = cfg.albedoSea ?? 0.06;
  const cdLand = cfg.cdLand ?? 0.0025;
  const cdSea = cfg.cdSea ?? 0.0011;
  const betaLand = cfg.betaLand ?? 0.65;
  const betaSea = cfg.betaSea ?? 1.0;
  const Cland = cfg.heatCapacityLand ?? 2.0e6;    // J/(m²·K)，约 0.5 m 湿土壤
  const Csea = cfg.heatCapacitySea ?? 2.0e8;      // J/(m²·K)，约 50 m 混合层
  const tauSea = cfg.tauSea ?? 30 * 86400;

  // 海陆掩膜：land=1
  const land = new Uint8Array(nXY);
  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      const col = i + nx * j;
      if (kind === 'sea') land[col] = 0;
      else if (kind === 'land') land[col] = 1;
      else if (kind === 'coast') land[col] = grid.xc(i) < grid.Lx * 0.5 ? 1 : 0;
      else if (kind === 'island') {
        const r = Math.hypot(grid.xc(i) - grid.Lx * 0.5, grid.yc(j) - grid.Ly * 0.5) / (grid.Lx * 0.22);
        land[col] = r < 1 ? 1 : 0;
      } else if (kind === 'patch') {
        // 一方陆、一方海，交界面在 y 方向（用于海陆风）
        land[col] = grid.yc(j) < grid.Ly * 0.5 ? 1 : 0;
      }
    }
  }

  const Tg = new Float64Array(nXY);
  const shFlux = new Float64Array(nXY);
  const lhFlux = new Float64Array(nXY);
  const cdCol = new Float64Array(nXY);
  const albedo = new Float64Array(nXY);

  for (let col = 0; col < nXY; col++) {
    albedo[col] = land[col] ? landAlbedo : seaAlbedo;
    cdCol[col] = land[col] ? cdLand : cdSea;
  }

  function init(base) {
    for (let col = 0; col < nXY; col++) {
      const kS = terrain ? terrain.surfaceK[col] : 0;
      Tg[col] = land[col] ? base.T0[kS] : sst;
    }
  }

  /** 计算当前时刻的地表通量（W/m²，向上为正） */
  function computeFluxes(S, base) {
    for (let col = 0; col < nXY; col++) {
      const kS = terrain ? terrain.surfaceK[col] : 0;
      const n = col + nXY * kS;
      if (mask && mask[n]) { shFlux[col] = 0; lhFlux[col] = 0; continue; }

      const rho = base.rho0[kS];
      const Tair = S.th[n] * base.ex[kS];
      const qair = S.qv[n];
      const u = S.u[n], v = S.v[n];
      let spd = Math.sqrt(u * u + v * v + S.w[n] * S.w[n]);
      if (spd < 2.0) spd = 2.0;            // 静风时仍保留弱交换
      const cd = cdCol[col];

      shFlux[col] = rho * CONST.cp * cd * spd * (Tg[col] - Tair);

      const qsTg = qsat(Math.max(200, Tg[col]), base.p0[kS]);
      let dq = (qsTg - qair);
      const beta = land[col] ? betaLand : betaSea;
      lhFlux[col] = rho * CONST.Lv * cd * spd * beta * dq;
    }
  }

  /**
   * 推进地表温度一个时间步。
   * @param {Float64Array} swSurf 到达地表的短波（W/m²，已在 rad3 里按云衰减）
   */
  function advance(dt, S, base, swSurf) {
    computeFluxes(S, base);
    for (let col = 0; col < nXY; col++) {
      const kS = terrain ? terrain.surfaceK[col] : 0;
      const n = col + nXY * kS;
      if (mask && mask[n]) continue;

      const Tair = S.th[n] * base.ex[kS];
      const swAbs = (swSurf ? swSurf[col] : 0) * (1 - albedo[col]);
      const lwNet = 0.95 * SIGMA * (Tg[col] ** 4 - Tair ** 4);
      const net = swAbs - lwNet - shFlux[col] - lhFlux[col];

      if (land[col]) {
        Tg[col] += net / Cland * dt;
        if (Tg[col] < 220) Tg[col] = 220;
        if (Tg[col] > 340) Tg[col] = 340;
      } else {
        // 海洋：大热容 + 向 SST 的慢回拉，日变化被压到 0.1 K 量级
        const relax = (sst - Tg[col]) / tauSea;
        Tg[col] += (net / Csea + relax) * dt;
        if (Tg[col] < 271) Tg[col] = 271;
        if (Tg[col] > 305) Tg[col] = 305;
      }
    }
  }

  function stats() {
    let landMin = 1e9, landMax = -1e9, seaMin = 1e9, seaMax = -1e9;
    let nLand = 0, nSea = 0;
    for (let col = 0; col < nXY; col++) {
      if (land[col]) {
        nLand++;
        if (Tg[col] < landMin) landMin = Tg[col];
        if (Tg[col] > landMax) landMax = Tg[col];
      } else {
        nSea++;
        if (Tg[col] < seaMin) seaMin = Tg[col];
        if (Tg[col] > seaMax) seaMax = Tg[col];
      }
    }
    return {
      landMean: nLand ? (landMin + landMax) * 0.5 : null,
      landMin, landMax, seaMin, seaMax, nLand, nSea,
    };
  }

  return { init, computeFluxes, advance, stats, Tg, land, shFlux, lhFlux, cdCol, albedo, kind };
}
