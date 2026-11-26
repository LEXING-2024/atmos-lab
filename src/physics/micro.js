/**
 * Kessler 暖雨微物理 + 饱和调整。
 * 接口契约：apply(state, base, dt) → 地面降水通量 kg/(m²·s)
 */

import { CONST, qsat, satIncrement } from './thermo.js';

export function createMicrophysics(grid, cfg) {
  const { nx, nz, dz, dx, Lx } = grid;
  const fluxCol = new Float64Array(nx);

  const k1 = cfg.k1 ?? 1e-3;
  const qc0 = cfg.qc0 ?? 0.6e-3;
  const k2 = cfg.k2 ?? 2.2;
  const evapTau = cfg.evapTau ?? 1200;   // 雨蒸发时标：3 s 会让雨在几百米内蒸干，与 Kessler 量级不符
  const rhoCritNeg = 1e-12;

  function apply(S, base, dt) {
    if (cfg.enable === false) return 0;

    const { p0, ex, rho0, rhoW } = base;
    const Lv = CONST.Lv, cp = CONST.cp, Rv = CONST.Rv;

    for (let k = 0; k < nz; k++) {
      const p = p0[k], exk = ex[k];
      const off = nx * k;
      for (let i = 0; i < nx; i++) {
        const n = off + i;
        let th = S.th[n], qv = S.qv[n], qc = S.qc[n], qr = S.qr[n];

        // 1. 饱和调整
        if (qv > 0 || qc > 0) {
          const T = th * exk;
          const dq = satIncrement(T, p, qv, qc);
          qv -= dq;
          qc += dq;
          th += Lv * dq / (cp * exk);
        }

        // 2. 自动转换
        if (qc > qc0) {
          let d = k1 * (qc - qc0) * dt;
          if (d > qc) d = qc;
          qc -= d;
          qr += d;
        }

        // 3. 碰并
        if (qc > 0 && qr > 0) {
          let d = k2 * qc * qr * dt;
          if (d > qc) d = qc;
          qc -= d;
          qr += d;
        }

        // 4. 雨蒸发（未饱和层）
        if (qr > 0 && qv > 0) {
          const T = th * exk;
          const qs = qsat(T, p);
          if (qv < qs) {
            const f = 1 + (Lv * Lv * qs) / (cp * Rv * T * T);
            let d = (1 - qv / qs) * qr * dt / evapTau;
            const capEv = (qs - qv) / f;
            if (d > capEv) d = capEv;
            if (d > qr) d = qr;
            qr -= d;
            qv += d;
            th -= Lv * d / (cp * exk);
          }
        }

        S.th[n] = th;
        S.qv[n] = qv > rhoCritNeg ? qv : 0;
        S.qc[n] = qc > 0 ? qc : 0;
        S.qr[n] = qr > 0 ? qr : 0;
      }
    }

    // 5. 雨滴下落（逐列、自上而下）
    // 地表通量用「下落前后水量的实际差值」确定，保证收支精确闭合
    let qrPathBefore = 0;
    for (let k = 0; k < nz; k++) {
      const off = nx * k;
      for (let i = 0; i < nx; i++) qrPathBefore += rho0[k] * dz * S.qr[off + i];
    }
    fluxCol.fill(0);
    for (let k = nz - 1; k >= 0; k--) {
      const off = nx * k;
      const rk = rho0[k], rW = rhoW[k];
      for (let i = 0; i < nx; i++) {
        const n = off + i;
        const qr = S.qr[n];
        let vt = 0;
        if (qr > 1e-12) {
          vt = 36.34 * Math.pow(rk * qr, 0.1364);
          if (vt > 12) vt = 12;
          if (vt < 0.5) vt = 0.5;
        }
        const fOut = rW * vt * qr;
        const dq = (fluxCol[i] - fOut) / (rk * dz);
        const nq = qr + dt * dq;
        S.qr[n] = nq > 0 ? nq : 0;
        fluxCol[i] = fOut;
      }
    }

    let qrPathAfter = 0;
    for (let k = 0; k < nz; k++) {
      const off = nx * k;
      for (let i = 0; i < nx; i++) qrPathAfter += rho0[k] * dz * S.qr[off + i];
    }
    // kg/(m²·s)，域平均
    return (qrPathBefore - qrPathAfter) / (dt * nx);
  }

  return { apply };
}
