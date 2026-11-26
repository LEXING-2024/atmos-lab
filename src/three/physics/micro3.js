/**
 * 三维 Kessler 暖雨微物理 + 饱和调整。
 * 接口契约：apply(state, base, dt) → 地面降水通量 kg/(m²·s)（水平域平均）
 */

import { CONST, qsat, satIncrement } from '../../physics/thermo.js';

export function createMicrophysics3(grid, cfg) {
  const { nx, ny, nz, dx, dy, dz } = grid;
  const nXY = nx * ny;
  const fluxCol = new Float64Array(nXY);

  const k1 = cfg.k1 ?? 1e-3;
  const qc0 = cfg.qc0 ?? 0.6e-3;
  const k2 = cfg.k2 ?? 2.2;
  const evapTau = cfg.evapTau ?? 1200;

  function apply(S, base, dt) {
    if (cfg.enable === false) return 0;

    const { p0, ex, rho0, rhoW } = base;
    const Lv = CONST.Lv, cp = CONST.cp, Rv = CONST.Rv;

    for (let k = 0; k < nz; k++) {
      const p = p0[k], exk = ex[k];
      const off = nXY * k;
      for (let j = 0; j < ny; j++) {
        const row = off + nx * j;
        for (let i = 0; i < nx; i++) {
          const n = row + i;
          let th = S.th[n], qv = S.qv[n], qc = S.qc[n], qr = S.qr[n];

          if (qv > 0 || qc > 0) {
            const T = th * exk;
            const dq = satIncrement(T, p, qv, qc);
            qv -= dq;
            qc += dq;
            th += Lv * dq / (cp * exk);
          }

          if (qc > qc0) {
            let d = k1 * (qc - qc0) * dt;
            if (d > qc) d = qc;
            qc -= d;
            qr += d;
          }

          if (qc > 0 && qr > 0) {
            let d = k2 * qc * qr * dt;
            if (d > qc) d = qc;
            qc -= d;
            qr += d;
          }

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
          S.qv[n] = qv > 1e-12 ? qv : 0;
          S.qc[n] = qc > 0 ? qc : 0;
          S.qr[n] = qr > 0 ? qr : 0;
        }
      }
    }

    // 雨滴下落：逐列自上而下，地表通量按水路径实际差值确定
    let qrPathBefore = 0;
    for (let k = 0; k < nz; k++) {
      const off = nXY * k, rk = rho0[k];
      for (let n = 0; n < nXY; n++) qrPathBefore += rk * dz * S.qr[off + n];
    }
    fluxCol.fill(0);
    for (let k = nz - 1; k >= 0; k--) {
      const off = nXY * k;
      const rk = rho0[k], rW = rhoW[k];
      for (let n = 0; n < nXY; n++) {
        const m = off + n;
        const qr = S.qr[m];
        let vt = 0;
        if (qr > 1e-12) {
          vt = 36.34 * Math.pow(rk * qr, 0.1364);
          if (vt > 12) vt = 12;
          if (vt < 0.5) vt = 0.5;
        }
        const fOut = rW * vt * qr;
        const dq = (fluxCol[n] - fOut) / (rk * dz);
        const nq = qr + dt * dq;
        S.qr[m] = nq > 0 ? nq : 0;
        fluxCol[n] = fOut;
      }
    }

    let qrPathAfter = 0;
    for (let k = 0; k < nz; k++) {
      const off = nXY * k, rk = rho0[k];
      for (let n = 0; n < nXY; n++) qrPathAfter += rk * dz * S.qr[off + n];
    }
    return (qrPathBefore - qrPathAfter) / (dt * nXY);
  }

  return { apply };
}
