/**
 * 辐射（牛顿冷却）+ 顶部海绵层。
 * 接口契约：apply(state, base, tendencies)
 */

export function createRadiation(grid, cfg) {
  const { nx, nz, dz, H } = grid;
  const tauRad = cfg.tauRad ?? 8 * 3600;
  const invTau = 1 / tauRad;
  const spongeMax = cfg.spongeMax ?? 1 / 300;

  const zSponge = H * (cfg.spongeDepth ?? 0.7);
  const alpha = new Float64Array(nz + 1);
  for (let k = 0; k <= nz; k++) {
    const z = k * dz;
    if (z <= zSponge) { alpha[k] = 0; continue; }
    const s = (z - zSponge) / (H - zSponge || 1);
    alpha[k] = spongeMax * s * s;
  }

  function apply(S, base, T) {
    const theta0 = base.theta0;

    for (let k = 0; k < nz; k++) {
      const off = nx * k;
      const a = alpha[k];
      const th0 = theta0[k];
      for (let i = 0; i < nx; i++) {
        const n = off + i;
        let d = -(S.th[n] - th0) * invTau;
        if (a > 0) {
          d -= a * (S.th[n] - th0);
          T.tu[n] -= a * S.u[n];
        }
        T.tth[n] += d;
      }
    }

    for (let k = 1; k < nz; k++) {
      const off = nx * k;
      const a = alpha[k];
      if (a <= 0) continue;
      for (let i = 0; i < nx; i++) {
        T.tw[off + i] -= a * S.w[off + i];
      }
    }
  }

  return { apply };
}
