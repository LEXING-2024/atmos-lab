/**
 * 三维微物理：混合相（Kessler 暖雨 + 三类冰）。
 *
 * 物种：qc 云水、qr 雨、qi 云冰、qs 雪、qg 霰。
 * 接口契约：apply(state, base, dt) → 地面降水通量 kg/(m²·s)（水平域平均）
 *
 * 冰相不是"把暖雨换个名字"：饱和目标随冰相占比在水面/冰面之间过渡（WBF），
 * 冻结、冰→雪自动转换、雪淞附成霰、下落、融化各有一条速率链，
 * 每一环的潜热都回写到 θ，因此冰化砧会真的比液态砧冷。
 */

import { CONST, qsat } from '../../physics/thermo.js';

const Lf = 3.34e5;          // 融化潜热 J/kg
const T0 = 273.15;

/** 冰面饱和水汽压（Tetens 冰面系数） */
function esatIce(T) {
  return 611.2 * Math.exp(21.875 * (T - T0) / (T - 7.66));
}

/** 冰面饱和比湿 */
function qsatIce(T, p) {
  const es = esatIce(T);
  const d = p - es;
  return CONST.eps * es / (d > 1 ? d : 1);
}

export function createMicrophysics3(grid, cfg) {
  const { nx, ny, nz, dx, dy, dz } = grid;
  const nXY = nx * ny;

  const k1 = cfg.k1 ?? 1e-3;              // 云水→雨自动转换
  const qc0 = cfg.qc0 ?? 0.6e-3;
  const k2 = cfg.k2 ?? 2.2;               // 雨碰并云水
  const evapTau = cfg.evapTau ?? 1200;

  const iceEnable = cfg.ice !== false;   // 关掉即退回纯暖雨（老预设行为不变）
  const freezeRate = cfg.freezeRate ?? 0.012;   // 每 °C 每 s 的冻结率
  const autoIce = cfg.autoIce ?? 2e-3;          // 冰→雪自动转换
  const qi0 = cfg.qi0 ?? 1e-4;
  const rimK = cfg.rimK ?? 1.2e3;               // 雪淞附云水成霰
  const meltTau = cfg.meltTau ?? 90;
  const vtSnow = cfg.vtSnow ?? 1.4;
  const vtGraupel = cfg.vtGraupel ?? 3.0;

  const fluxR = new Float64Array(nXY);
  const fluxS = new Float64Array(nXY);
  const fluxG = new Float64Array(nXY);

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
          let qi = S.qi[n], qs = S.qs[n], qg = S.qg[n];

          const T = th * exk;
          const iceMass = qi + qs + qg;
          const condMass = qc + iceMass;

          // ---- 1. 饱和调整：冰相占比决定用冰面还是水面饱和 ----
          if (qv > 0 || condMass > 0) {
            const fIce = iceEnable ? (condMass > 1e-12 ? iceMass / condMass : (T < T0 ? 1 : 0)) : 0;
            const qsW = qsat(T, p);
            const qsI = qsatIce(T, p);
            const qsT = qsW * (1 - fIce) + qsI * fIce;
            const Leff = Lv * (1 - fIce) + (Lv + Lf) * fIce;
            const f = 1 + (Leff * Leff * qsT) / (cp * Rv * T * T);
            let dq = (qv - qsT) / f;

            if (dq > 0) {
              // 多余水汽：冷处直接凝华/凝露，按冰相占比分配
              dq = Math.min(dq, qv);
              qv -= dq;
              if (fIce > 0.5 && T < T0) { qi += dq; } else { qc += dq; }
              th += Leff * dq / (cp * exk);
            } else if (condMass > 0) {
              // 欠饱和：从现有凝结物蒸发/升华，优先液相（液相更快）
              dq = -dq;
              if (dq > condMass) dq = condMass;
              let take = dq;
              if (qc > 0) { const d = Math.min(qc, take); qc -= d; take -= d; }
              if (take > 0 && qi > 0) { const d = Math.min(qi, take); qi -= d; take -= d; }
              if (take > 0 && qs > 0) { const d = Math.min(qs, take); qs -= d; take -= d; }
              if (take > 0 && qg > 0) { const d = Math.min(qg, take); qg -= d; take -= d; }
              const real = dq - take;
              qv += real;
              th -= Leff * real / (cp * exk);
            }
          }

          if (iceEnable) {
            // ---- 2. 云水冻结（T<0，按温度线性加速）----
            if (qc > 0 && T < T0) {
              const fr = Math.min(1, freezeRate * (T0 - T) * dt);
              const d = qc * fr;
              qc -= d; qi += d;
              th += Lf * d / (cp * exk);          // 冻结放热
            }

            // ---- 3. 冰 → 雪（自动转换）----
            if (qi > qi0) {
              const d = Math.min(qi - qi0, autoIce * (qi - qi0) * dt);
              qi -= d; qs += d;
            }

            // ---- 4. 雪淞附云水 → 霰（冷于 -5°C）----
            if (qc > 1e-7 && qs > 1e-7 && T < 268.15) {
              let d = rimK * qc * qs * dt;
              if (d > qc) d = qc;
              if (d > qs) d = qs;
              qc -= d; qs -= d; qg += d;
              th += Lf * d / (cp * exk);
            }

            // ---- 5. 融化（T>0）----
            if (T > T0) {
              const rate = Math.min(1, (T - T0) / 5) * dt / meltTau;
              const dI = qi * Math.min(1, rate * 3);
              const dS = qs * Math.min(1, rate);
              const dG = qg * Math.min(1, rate * 1.5);
              qi -= dI; qc += dI;
              qs -= dS; qr += dS;
              qg -= dG; qr += dG;
              th -= Lf * (dI + dS + dG) / (cp * exk);
            }

            // ---- 6. 过冷雨冻结 ----
            if (qr > 1e-7 && T < 268.15) {
              const d = qr * Math.min(1, dt / 60);
              qr -= d; qg += d;
              th += Lf * d / (cp * exk);
            }
          }

          // ---- 7. 暖雨链（与冰相共存）----
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
            const qsW = qsat(T, p);
            if (qv < qsW) {
              const f = 1 + (Lv * Lv * qsW) / (cp * Rv * T * T);
              let d = (1 - qv / qsW) * qr * dt / evapTau;
              const capEv = (qsW - qv) / f;
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
          S.qi[n] = qi > 0 ? qi : 0;
          S.qs[n] = qs > 0 ? qs : 0;
          S.qg[n] = qg > 0 ? qg : 0;
        }
      }
    }

    // ---- 下落：雨、雪、霰各自一列自上而下，地面通量按水路径实际差值确定 ----
    let pathBefore = 0;
    for (let k = 0; k < nz; k++) {
      const off = nXY * k, rk = rho0[k];
      for (let n = 0; n < nXY; n++) {
        pathBefore += rk * dz * (S.qr[off + n] + S.qs[off + n] + S.qg[off + n]);
      }
    }

    fluxR.fill(0); fluxS.fill(0); fluxG.fill(0);
    for (let k = nz - 1; k >= 0; k--) {
      const off = nXY * k;
      const rk = rho0[k], rW = rhoW[k];
      for (let n = 0; n < nXY; n++) {
        const m = off + n;
        const qr = S.qr[m], qs = S.qs[m], qg = S.qg[m];

        let vR = 0;
        if (qr > 1e-12) {
          vR = 36.34 * Math.pow(rk * qr, 0.1364);
          if (vR > 12) vR = 12;
          if (vR < 0.5) vR = 0.5;
        }
        const vS = qs > 1e-12 ? vtSnow * Math.pow(rk / 1.0, 0.15) : 0;
        const vG = qg > 1e-12 ? vtGraupel * Math.pow(rk / 1.0, 0.15) : 0;

        const fR = rW * vR * qr;
        const fS = rW * vS * qs;
        const fG = rW * vG * qg;
        const sumIn = fluxR[n] + fluxS[n] + fluxG[n];
        const sumOut = fR + fS + fG;
        const dq = (sumIn - sumOut) / (rk * dz);
        const newR = qr + dt * dq * (fR > 0 ? fR / (sumOut || 1) : 0);
        const newS = qs + dt * dq * (fS > 0 ? fS / (sumOut || 1) : 0);
        const newG = qg + dt * dq * (fG > 0 ? fG / (sumOut || 1) : 0);
        S.qr[m] = newR > 0 ? newR : 0;
        S.qs[m] = newS > 0 ? newS : 0;
        S.qg[m] = newG > 0 ? newG : 0;
        fluxR[n] = fR; fluxS[n] = fS; fluxG[n] = fG;
      }
    }

    let pathAfter = 0;
    for (let k = 0; k < nz; k++) {
      const off = nXY * k, rk = rho0[k];
      for (let n = 0; n < nXY; n++) {
        pathAfter += rk * dz * (S.qr[off + n] + S.qs[off + n] + S.qg[off + n]);
      }
    }
    void dx; void dy;
    return (pathBefore - pathAfter) / (dt * nXY);
  }

  return { apply, iceEnable };
}
