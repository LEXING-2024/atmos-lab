/**
 * 三维辐射：日变化短波 + 牛顿冷却 + 顶部海绵层（v1.6）。
 *
 * 短波：太阳高度角由纬度、赤纬和模拟时间（小时）决定；
 *       沿每列自上而下积分，空气与凝结物各自吸收，剩余部分到达地表（swSurf）。
 * 长波：仍用牛顿冷却近似（明确不是辐射传输），云顶长波冷却不单独算。
 */

const DEG = Math.PI / 180;

export function createRadiation3(grid, cfg) {
  const { nx, ny, nz, dz, H, nXY } = grid;
  const tauRad = cfg.tauRad ?? 8 * 3600;
  const invTau = 1 / tauRad;
  const spongeMax = cfg.spongeMax ?? 1 / 300;
  const zSponge = H * (cfg.spongeDepth ?? 0.7);

  const sun = cfg.sun ?? {};
  const sunEnable = sun.enable === true;
  const lat = (sun.lat ?? 30) * DEG;
  const decl = (sun.decl ?? 8) * DEG;
  const hour0 = sun.hour0 ?? 6;
  const s0 = sun.s0 ?? 1361;
  const kAbsAir = sun.kAbsAir ?? 2.2e-5;      // 1/m，整层光学厚度 ≈ 0.2（9 km）
  const kExtCloud = sun.kExtCloud ?? 130;     // m²/kg，LWP 0.05 kg/m² → τ ≈ 6.5

  const swSurf = new Float64Array(nXY);

  const alpha = new Float64Array(nz + 1);
  for (let k = 0; k <= nz; k++) {
    const z = k * dz;
    if (z <= zSponge) { alpha[k] = 0; continue; }
    const s = (z - zSponge) / (H - zSponge || 1);
    alpha[k] = spongeMax * s * s;
  }

  /** 太阳天顶角余弦（模拟时间 time 秒） */
  function cosZenith(time) {
    if (!sunEnable) return 0;
    const hour = ((hour0 + time / 3600) % 24 + 24) % 24;
    const Ha = (hour - 12) * 15 * DEG;
    const mu = Math.sin(lat) * Math.sin(decl) + Math.cos(lat) * Math.cos(decl) * Math.cos(Ha);
    return mu > 0 ? mu : 0;
  }

  function apply(S, base, T, time = 0) {
    const theta0 = base.theta0;
    const cp = 1004.5;

    // 短波：逐列自上而下
    const mu = cosZenith(time);
    if (mu > 0) {
      const s0mu = s0 * mu;
      for (let col = 0; col < nXY; col++) {
        let F = s0mu;
        for (let k = nz - 1; k >= 0; k--) {
          const n = col + nXY * k;
          const qTot = S.qc[n] + S.qi[n] + S.qs[n] + S.qg[n];
          const dtau = dz * (kAbsAir + kExtCloud * base.rho0[k] * qTot);
          const abs = F * (1 - Math.exp(-dtau));
          const exk = base.ex[k];
          T.tth[n] += abs / (base.rho0[k] * cp * dz * exk);
          F -= abs;
        }
        swSurf[col] = F;
      }
    } else {
      swSurf.fill(0);
    }

    // 牛顿冷却 + 海绵
    for (let k = 0; k < nz; k++) {
      const off = nXY * k;
      const a = alpha[k];
      const th0 = theta0[k];
      for (let n = 0; n < nXY; n++) {
        const m = off + n;
        let d = -(S.th[m] - th0) * invTau;
        if (a > 0) {
          d -= a * (S.th[m] - th0);
          T.tu[m] -= a * S.u[m];
          T.tv[m] -= a * S.v[m];
        }
        T.tth[m] += d;
      }
    }

    for (let k = 1; k < nz; k++) {
      const off = nXY * k;
      const a = alpha[k];
      if (a <= 0) continue;
      for (let n = 0; n < nXY; n++) T.tw[off + n] -= a * S.w[off + n];
    }
  }

  return { apply, cosZenith, swSurf, sunEnable };
}
