/**
 * 带固壁掩膜的投影求解器（IBM / staircase cut）。
 *
 * 与 poisson3.js 的差别只有一处：水平不再是周期 FFT，而是把固壁面当成 Neumann
 * （法向通量 0），用预条件共轭梯度（Jacobi PCG）解同样离散的压力方程。
 *
 * 算子（与 anelastic3 的散度/梯度严格互为伴随）：
 *   A(π)_n = Σ_faces c_f (π_n − π_nb),  c_f = ρ_f / d²
 *   水平 c = ρ0_k/Δx² 或 ρ0_k/Δy²；垂直 c = ρW_k/Δz²
 *   k=0 的底面与 k=nz−1 的顶面取 Neumann（不连外界）
 *
 * 无固壁时请用 poisson3.js（快两个量级）。
 */

export function createPoissonSolver3IBM(grid, mask) {
  const { nx, ny, nz, dx, dy, dz } = grid;
  const nXY = nx * ny;
  const nCell = grid.nCell;

  const dx2 = dx * dx, dy2 = dy * dy, dz2 = dz * dz;

  // 每个格点 6 个方向的导纳（0 = 固壁面 / 域外）
  const cxp = new Float64Array(nCell);   // 东
  const cxm = new Float64Array(nCell);   // 西
  const cyp = new Float64Array(nCell);   // 北
  const cym = new Float64Array(nCell);   // 南
  const czb = new Float64Array(nCell);   // 下
  const czt = new Float64Array(nCell);   // 上
  const diag = new Float64Array(nCell);

  const r = new Float64Array(nCell);
  const z = new Float64Array(nCell);
  const p = new Float64Array(nCell);
  const q = new Float64Array(nCell);
  const div = new Float64Array(nCell);

  /** 依据 mask 重算导纳。注意：ρ 在散度与梯度中相互抵消，算子就是纯 Laplacian。 */
  function buildCoefficients() {
    const fluid = (n) => mask[n] === 0;
    const hx = 1 / dx2, hy = 1 / dy2, hz = 1 / dz2;
    for (let k = 0; k < nz; k++) {
      const off = nXY * k;
      for (let j = 0; j < ny; j++) {
        for (let i = 0; i < nx; i++) {
          const n = off + nx * j + i;
          if (!fluid(n)) {
            cxp[n] = cxm[n] = cyp[n] = cym[n] = czb[n] = czt[n] = 0;
            diag[n] = 1;
            continue;
          }
          const ip = i + 1 === nx ? 0 : i + 1;
          const im = i === 0 ? nx - 1 : i - 1;
          const jp = j + 1 === ny ? 0 : j + 1;
          const jm = j === 0 ? ny - 1 : j - 1;
          cxp[n] = fluid(off + nx * j + ip) ? hx : 0;
          cxm[n] = fluid(off + nx * j + im) ? hx : 0;
          cyp[n] = fluid(off + nx * jp + i) ? hy : 0;
          cym[n] = fluid(off + nx * jm + i) ? hy : 0;
          czb[n] = k > 0 && fluid(off - nXY + nx * j + i) ? hz : 0;
          czt[n] = k + 1 < nz && fluid(off + nXY + nx * j + i) ? hz : 0;
          const d = cxp[n] + cxm[n] + cyp[n] + cym[n] + czb[n] + czt[n];
          diag[n] = d > 0 ? d : 1;
        }
      }
    }
  }

  function applyA(src, dst) {
    for (let k = 0; k < nz; k++) {
      const off = nXY * k;
      for (let j = 0; j < ny; j++) {
        for (let i = 0; i < nx; i++) {
          const n = off + nx * j + i;
          if (mask[n]) { dst[n] = src[n]; continue; }
          let v = diag[n] * src[n];
          if (cxp[n]) v -= cxp[n] * src[off + nx * j + (i + 1 === nx ? 0 : i + 1)];
          if (cyp[n]) v -= cyp[n] * src[off + nx * (j + 1 === ny ? 0 : j + 1) + i];
          if (czb[n]) v -= czb[n] * src[off - nXY + nx * j + i];
          if (czt[n]) v -= czt[n] * src[off + nXY + nx * j + i];
          if (cxm[n]) v -= cxm[n] * src[off + nx * j + (i === 0 ? nx - 1 : i - 1)];
          if (cym[n]) v -= cym[n] * src[off + nx * (j === 0 ? ny - 1 : j - 1) + i];
          dst[n] = v;
        }
      }
    }
  }

  /** 密度加权散度（固壁内部置 0） */
  function computeDivMasked(u, v, w, base, out) {
    const rho0 = base.rho0, rhoW = base.rhoW;
    let maxAbs = 0;
    for (let k = 0; k < nz; k++) {
      const off = nXY * k;
      const rk = rho0[k];
      for (let j = 0; j < ny; j++) {
        for (let i = 0; i < nx; i++) {
          const n = off + nx * j + i;
          if (mask[n]) { out[n] = 0; continue; }
          const ip = i + 1 === nx ? 0 : i + 1;
          const jp = j + 1 === ny ? 0 : j + 1;
          const d = rk * ((u[off + nx * j + ip] - u[n]) / dx + (v[off + nx * jp + i] - v[n]) / dy)
            + (rhoW[k + 1] * w[n + nXY] - rhoW[k] * w[n]) / dz;
          out[n] = d;
          const a = d < 0 ? -d : d;
          if (a > maxAbs) maxAbs = a;
        }
      }
    }
    return maxAbs;
  }

  let lastIterations = 0;

  /**
   * @param {Float64Array} u v w 会被原地修正
   * @param {Float64Array} pi 输出
   * @returns {{divBefore:number, divAfter:number, iterations:number}}
   */
  function project(u, v, w, dt, pi, base) {
    if (!base) throw new Error('project 需要 base 提供密度剖面');
    buildCoefficients();

    const divBefore = computeDivMasked(u, v, w, base, div);

    // A π = b，b = −div/Δt
    let bsum = 0, nfluid = 0;
    for (let n = 0; n < nCell; n++) {
      if (mask[n]) { z[n] = 0; continue; }
      const b = -div[n] / dt;
      z[n] = b;
      bsum += b; nfluid++;
    }
    const bmean = nfluid > 0 ? bsum / nfluid : 0;   // Neumann 相容性：零均值
    let bn = 0;
    for (let n = 0; n < nCell; n++) {
      if (mask[n]) continue;
      z[n] -= bmean;
      bn += z[n] * z[n];
    }
    bn = Math.sqrt(bn);

    // 以 π 的上一帧解为初值（时间步之间变化很小，通常能省 3–5 倍迭代）
    for (let n = 0; n < nCell; n++) if (mask[n]) pi[n] = 0;
    applyA(pi, q);
    for (let n = 0; n < nCell; n++) r[n] = mask[n] ? 0 : z[n] - q[n];
    // 残差里若含常数分量（Neumann 算子的零空间），CG 无法消掉，必须先投影掉
    let rsum = 0;
    for (let n = 0; n < nCell; n++) rsum += r[n];
    const rmean = nfluid > 0 ? rsum / nfluid : 0;
    for (let n = 0; n < nCell; n++) if (!mask[n]) r[n] -= rmean;
    for (let n = 0; n < nCell; n++) q[n] = mask[n] ? 0 : r[n] / diag[n];
    for (let n = 0; n < nCell; n++) p[n] = q[n];

    let rz = 0;
    for (let n = 0; n < nCell; n++) rz += r[n] * q[n];

    // 散度已在机器精度：不需要求解，也不要拿旧的 π 去修速度
    if (bn < 1e-13) {
      lastIterations = 0;
      const divAfter0 = computeDivMasked(u, v, w, base, div);
      return { divBefore, divAfter: divAfter0, iterations: 0 };
    }

    const tol = 1e-4 * bn + 1e-13;
    const maxIter = 80;
    let it = 0;
    for (; it < maxIter; it++) {
      applyA(p, q);                     // q = A p
      let pq = 0;
      for (let n = 0; n < nCell; n++) if (!mask[n]) pq += p[n] * q[n];
      if (!(pq > 0)) break;
      const alpha = rz / pq;
      for (let n = 0; n < nCell; n++) { pi[n] += alpha * p[n]; r[n] -= alpha * q[n]; }
      let rn = 0;
      for (let n = 0; n < nCell; n++) rn += r[n] * r[n];
      if (Math.sqrt(rn) < tol) { it++; break; }

      for (let n = 0; n < nCell; n++) z[n] = mask[n] ? 0 : r[n] / diag[n];
      let rzNew = 0;
      for (let n = 0; n < nCell; n++) rzNew += r[n] * z[n];
      const beta = rzNew / rz;
      rz = rzNew;
      for (let n = 0; n < nCell; n++) p[n] = z[n] + beta * p[n];
      if (globalThis.__IBM_DEBUG && it % 10 === 9) {
        let rr = 0;
        for (let n = 0; n < nCell; n++) rr += r[n] * r[n];
        console.log(`   it ${it + 1} |r| ${Math.sqrt(rr).toExponential(2)} (tol ${tol.toExponential(2)}, bn ${bn.toExponential(2)})`);
      }
    }
    lastIterations = it;

    // 零均值（Neumann 下 π 只定到常数）
    let psum = 0, cnt = 0;
    for (let n = 0; n < nCell; n++) if (!mask[n]) { psum += pi[n]; cnt++; }
    const pmean = cnt > 0 ? psum / cnt : 0;
    for (let n = 0; n < nCell; n++) pi[n] -= pmean;

    const rho0 = base.rho0, rhoW = base.rhoW;
    for (let k = 0; k < nz; k++) {
      const off = nXY * k, rk = rho0[k];
      for (let j = 0; j < ny; j++) {
        for (let i = 0; i < nx; i++) {
          const n = off + nx * j + i;
          if (mask[n]) { u[n] = 0; v[n] = 0; continue; }
          const im = i === 0 ? nx - 1 : i - 1;
          const jm = j === 0 ? ny - 1 : j - 1;
          if (mask[off + nx * j + im]) u[n] = 0;
          else u[n] -= dt * (pi[n] - pi[off + nx * j + im]) / (dx * rk);
          if (mask[off + nx * jm + i]) v[n] = 0;
          else v[n] -= dt * (pi[n] - pi[off + nx * jm + i]) / (dy * rk);
        }
      }
    }
    for (let k = 1; k < nz; k++) {
      const off = nXY * k;
      for (let j = 0; j < ny; j++) {
        for (let i = 0; i < nx; i++) {
          const n = off + nx * j + i;
          if (mask[n] || mask[n - nXY]) w[n] = 0;
          else w[n] -= dt * (pi[n] - pi[n - nXY]) / (dz * rhoW[k]);
        }
      }
    }
    for (let n = 0; n < nXY; n++) { w[n] = 0; w[nXY * nz + n] = 0; }

    const divAfter = computeDivMasked(u, v, w, base, div);
    return { divBefore, divAfter, iterations: it };
  }

  return { project, computeDiv: computeDivMasked, mask, get iterations() { return lastIterations; } };
}
