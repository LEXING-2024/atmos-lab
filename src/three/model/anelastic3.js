/**
 * 三维非静力伪不可压大气模式（x–y–z）。
 *
 * 与二维版同一套离散：通量形式标量平流（单调限制器）、三阶迎风动量平流、
 * Δ⁴ 滤波后重新投影、RK3 时间推进、Kessler 暖雨、Smagorinsky 湍流。
 * 所有网格尺寸来自 Grid3，代码中不存在网格常量。
 */

import { makeBase } from '../../model/base.js';
import { allocState3, copyState3, combineState3, clampMoisture3, repairMoisture3, maxAbs3, FIELD_NAMES3 } from './state3.js';
import { createPoissonSolver3 } from '../core/poisson3.js';
import { createPoissonSolver3IBM } from '../core/poisson3_ibm.js';
import { makeTerrain, maskVelocities, resetSolidState } from '../core/terrain3.js';
import { createMicrophysics3 } from '../physics/micro3.js';
import { createTurbulence3 } from '../physics/turb3.js';
import { createRadiation3 } from '../physics/rad3.js';
import { CONST } from '../../physics/thermo.js';

const nowMs = () => (
  typeof performance !== 'undefined' && performance.now ? performance.now() : Date.now()
);

export const DEFAULT_CONFIG3 = {
  dtMax: 6,
  cfl: 0.4,
  filter: 0.02,
  micro: {
    enable: true, k1: 1e-3, qc0: 0.6e-3, k2: 2.2, evapTau: 1200,
    ice: true, freezeRate: 0.012, autoIce: 2e-3, qi0: 1e-4, rimK: 1.2e3, meltTau: 90,
    vtSnow: 1.4, vtGraupel: 3.0,
  },
  turb: {
    cs: 0.15, kMin: 1.0, kMax: 1500, cd: 0.002,
    surface: { sensible: 0, latent: 0, xCenter: null, yCenter: null, xWidth: 4000, yWidth: 4000 },
  },
  rad: { tauRad: 8 * 3600, spongeMax: 1 / 300, spongeDepth: 0.7 },
  // 环境风初值：u0/v0 为地面风，shearU/shearV 为线性风切变（1/s）
  // 缺省 0 = 静止大气（与二维一致）；给上非零切变是打破 x–y 对称最直接的办法
  // shearDepth：切变只加到这一层高度，之上保持常数，避免 9 km 域顶部出现 30 m/s 的假急流
  wind: { u0: 0, v0: 0, shearU: 0, shearV: 0, zRef: 0, shearDepth: 6000 },
};

export class AnelasticModel3 {
  constructor(grid, opts = {}) {
    this.grid = grid;
    this.opts = opts;
    this.cfg = {
      ...DEFAULT_CONFIG3,
      ...opts,
      micro: { ...DEFAULT_CONFIG3.micro, ...(opts.micro || {}) },
      turb: {
        ...DEFAULT_CONFIG3.turb,
        ...(opts.turb || {}),
        surface: { ...DEFAULT_CONFIG3.turb.surface, ...((opts.turb || {}).surface || {}) },
      },
      rad: { ...DEFAULT_CONFIG3.rad, ...(opts.rad || {}) },
      wind: { ...DEFAULT_CONFIG3.wind, ...(opts.wind || {}) },
    };

    this.base = makeBase(grid, opts.base || {});
    this.terrain = opts.terrain ? makeTerrain(grid, opts.terrain) : null;
    this.solver = this.terrain
      ? createPoissonSolver3IBM(grid, this.terrain.mask)
      : createPoissonSolver3(grid);
    this.micro = createMicrophysics3(grid, this.cfg.micro);
    this.turb = createTurbulence3(grid, this.cfg.turb, this.terrain);
    this.rad = createRadiation3(grid, this.cfg.rad);
    this.ibmIterations = 0;

    const nCell = grid.nCell, nW = grid.nW;
    this.S = allocState3(grid);
    this.B1 = allocState3(grid);
    this.B2 = allocState3(grid);
    this.S2 = allocState3(grid);

    this.T = {
      tu: new Float64Array(nCell),
      tv: new Float64Array(nCell),
      tw: new Float64Array(nW),
      tth: new Float64Array(nCell),
      tqv: new Float64Array(nCell),
      tqc: new Float64Array(nCell),
      tqr: new Float64Array(nCell),
      tqi: new Float64Array(nCell),
      tqs: new Float64Array(nCell),
      tqg: new Float64Array(nCell),
    };
    this.pi = new Float64Array(nCell);
    this.Fx = new Float64Array(grid.nx);
    this.Fy = new Float64Array(grid.ny);
    this.tmp3 = new Float64Array(Math.max(nCell, nW));   // w 比标量多一层

    this.time = 0;
    this.maxSpeed = 8;
    this.precipRate = 0;
    this.precipAccum = 0;
    this.waterInput = 0;
    this.divResidual = 0;
    this.cflNow = 0;
    this.lastDt = 0;
    this.initMass = 0;
    this.repairAmount = 0;
    this.repairLost = 0;
  }

  /* --------------------------- 初始化 --------------------------- */

  initAtmosphere(o = {}) {
    const g = this.grid, base = this.base, S = this.S;
    const amp = o.amplitude ?? 2.5;
    const radius = o.radius ?? 2000;
    const zc = o.zCenter ?? 1500;
    const xc = o.xCenter ?? g.Lx * 0.5;
    const yc = o.yCenter ?? g.Ly * 0.5;
    const noise = o.noise ?? 0.02;
    const seed = o.seed ?? 20261127;

    for (let k = 0; k < g.nz; k++) {
      const off = g.nXY * k;
      for (let n = 0; n < g.nXY; n++) {
        S.th[off + n] = base.theta0[k];
        S.qv[off + n] = base.qv0[k];
        S.qc[off + n] = 0;
        S.qr[off + n] = 0;
        S.qi[off + n] = 0;
        S.qs[off + n] = 0;
        S.qg[off + n] = 0;
      }
    }
    S.u.fill(0);
    S.v.fill(0);
    S.w.fill(0);
    this.pi.fill(0);

    // 环境风：u(z) = u0 + shearU·min(z, shearDepth)。只作为初值，之后由方程自己演化。
    const wind = { ...this.cfg.wind, ...(o.wind || {}) };
    if (wind.u0 || wind.v0 || wind.shearU || wind.shearV) {
      for (let k = 0; k < g.nz; k++) {
        const z = g.zc(k), off = g.nXY * k;
        const zs = Math.min(z - wind.zRef, wind.shearDepth);
        const uu = wind.u0 + wind.shearU * zs;
        const vv = wind.v0 + wind.shearV * zs;
        for (let n = 0; n < g.nXY; n++) {
          S.u[off + n] = uu;
          S.v[off + n] = vv;
        }
      }
    }

    for (let k = 0; k < g.nz; k++) {
      const z = g.zc(k), off = g.nXY * k;
      for (let j = 0; j < g.ny; j++) {
        const y = g.yc(j);
        for (let i = 0; i < g.nx; i++) {
          const x = g.xc(i);
          const r = Math.sqrt(((x - xc) / radius) ** 2 + ((y - yc) / radius) ** 2 + ((z - zc) / radius) ** 2);
          if (r < 1) S.th[off + g.nx * j + i] += amp * 0.5 * (1 + Math.cos(Math.PI * r));
        }
      }
    }

    let s = seed >>> 0;
    const rnd = () => {
      s = (s * 1664525 + 1013904223) >>> 0;
      return s / 4294967296 - 0.5;
    };
    for (let k = 0; k < g.nz && g.zc(k) < 3000; k++) {
      const off = g.nXY * k;
      for (let n = 0; n < g.nXY; n++) S.th[off + n] += noise * rnd();
    }

    this.time = 0;
    this.precipAccum = 0;
    this.precipRate = 0;
    this.waterInput = 0;
    if (this.terrain) {
      resetSolidState(this.S, this.base, g, this.terrain);
      maskVelocities(this.S, g, this.terrain);
    }
    this.initMass = this.totalWater();
  }

  /* --------------------------- 诊断 --------------------------- */

  totalWater() {
    const g = this.grid, base = this.base, S = this.S;
    const m0 = this.terrain ? this.terrain.mask : null;
    let sum = 0;
    for (let k = 0; k < g.nz; k++) {
      const off = g.nXY * k, rk = base.rho0[k];
      for (let n = 0; n < g.nXY; n++) {
        const m = off + n;
        if (m0 && m0[m]) continue;
        sum += rk * (S.qv[m] + S.qc[m] + S.qr[m] + S.qi[m] + S.qs[m] + S.qg[m]);
      }
    }
    return sum * g.dx * g.dy * g.dz;
  }

  /** 冰相总量（云冰+雪+霰），kg/m² 的域平均 */
  iceWater() {
    const g = this.grid, base = this.base, S = this.S;
    const m0 = this.terrain ? this.terrain.mask : null;
    let sum = 0;
    for (let k = 0; k < g.nz; k++) {
      const off = g.nXY * k, rk = base.rho0[k];
      for (let n = 0; n < g.nXY; n++) {
        const m = off + n;
        if (m0 && m0[m]) continue;
        sum += rk * (S.qi[m] + S.qs[m] + S.qg[m]) * g.dz;
      }
    }
    return sum / g.nXY;   // kg/m²（≈ 等价水柱 mm）
  }

  diagnostics() {
    const g = this.grid, S = this.S, nXY = g.nXY;
    const m0 = this.terrain ? this.terrain.mask : null;
    let maxW = 0, maxU = 0, cloudTop = 0;
    for (let k = 0; k < g.nz; k++) {
      const off = nXY * k;
      for (let n = 0; n < nXY; n++) {
        const m = off + n;
        if (m0 && m0[m]) continue;
        const w = Math.abs(S.w[m] + S.w[m + nXY]) * 0.5;
        if (w > maxW) maxW = w;
        const spd = Math.sqrt(S.u[m] * S.u[m] + S.v[m] * S.v[m]);
        if (spd > maxU) maxU = spd;
        if (S.qc[m] + S.qr[m] + S.qi[m] + S.qs[m] + S.qg[m] > 1e-5 && g.zc(k) > cloudTop) cloudTop = g.zc(k);
      }
    }
    const cfl = Math.max(
      (maxU * this.lastDt) / g.dx,
      (maxW * this.lastDt) / g.dz,
    );
    this.cflNow = cfl;
    return {
      time: this.time,
      maxW,
      maxU,
      cloudTop,
      precipRate: this.precipRate,
      precipAccum: this.precipAccum,
      cfl,
      divResidual: this.divResidual,
      totalWater: this.totalWater(),
      waterInput: this.waterInput,
      cells: g.nCell,
    };
  }

  /* --------------------------- 倾向 --------------------------- */

  computeTendencies(S, T) {
    const g = this.grid, base = this.base;
    const { nx, ny, nz, dx, dy, dz } = g;
    const nXY = nx * ny;

    T.tu.fill(0); T.tv.fill(0); T.tw.fill(0);
    T.tth.fill(0); T.tqv.fill(0); T.tqc.fill(0); T.tqr.fill(0);
    T.tqi.fill(0); T.tqs.fill(0); T.tqg.fill(0);

    const theta0 = base.theta0, qv0 = base.qv0;
    const rho0 = base.rho0, rhoW = base.rhoW;
    const gg = CONST.g;

    const scalars = [S.th, S.qv, S.qc, S.qr, S.qi, S.qs, S.qg];
    const tscalars = [T.tth, T.tqv, T.tqc, T.tqr, T.tqi, T.tqs, T.tqg];
    const Fx = this.Fx, Fy = this.Fy;

    for (let s = 0; s < 7; s++) {
      const f = scalars[s], tf = tscalars[s];

      // ---- x ----
      for (let k = 0; k < nz; k++) {
        const rk = rho0[k];
        for (let j = 0; j < ny; j++) {
          const off = nXY * k + nx * j;
          for (let i = 0; i < nx; i++) {
            const ip = i + 1 === nx ? 0 : i + 1;
            const im = i === 0 ? nx - 1 : i - 1;
            const im2 = i <= 1 ? (i - 2 + nx) % nx : i - 2;
            const vf = S.u[off + i];   // 面速度取本格 u，与散度算子 ρ(u_{i+1}-u_i)/dx 严格一致
            const a = f[off + im], b = f[off + i], c = f[off + ip], d = f[off + im2];
            let fq = vf > 0 ? (6 * a + 3 * b - d) / 8 : (6 * b + 3 * a - c) / 8;
            const lo = a < b ? a : b, hi = a < b ? b : a;
            if (fq < lo) fq = lo; else if (fq > hi) fq = hi;
            Fx[i] = rk * vf * fq;
          }
          for (let i = 0; i < nx; i++) {
            const ip = i + 1 === nx ? 0 : i + 1;
            tf[off + i] -= (Fx[ip] - Fx[i]) / (dx * rk);
          }
        }
      }

      // ---- y ----
      for (let k = 0; k < nz; k++) {
        const rk = rho0[k];
        for (let i = 0; i < nx; i++) {
          for (let j = 0; j < ny; j++) {
            const jp = j + 1 === ny ? 0 : j + 1;
            const jm = j === 0 ? ny - 1 : j - 1;
            const jm2 = j <= 1 ? (j - 2 + ny) % ny : j - 2;
            const n = nXY * k + nx * j + i;
            const vf = S.v[n];         // 同上：与 ρ(v_{j+1}-v_j)/dy 一致
            const a = f[nXY * k + nx * jm + i], b = f[n], c = f[nXY * k + nx * jp + i], d = f[nXY * k + nx * jm2 + i];
            let fq = vf > 0 ? (6 * a + 3 * b - d) / 8 : (6 * b + 3 * a - c) / 8;
            const lo = a < b ? a : b, hi = a < b ? b : a;
            if (fq < lo) fq = lo; else if (fq > hi) fq = hi;
            Fy[j] = rk * vf * fq;
          }
          for (let j = 0; j < ny; j++) {
            const jp = j + 1 === ny ? 0 : j + 1;
            const n = nXY * k + nx * j + i;
            tf[n] -= (Fy[jp] - Fy[j]) / (dy * rk);
          }
        }
      }

      // ---- z ----
      for (let k = 0; k < nz; k++) {
        const off = nXY * k;
        for (let j = 0; j < ny; j++) {
          const row = off + nx * j;
          for (let i = 0; i < nx; i++) {
            const n = row + i;
            const c0 = f[n];
            const km = k > 0 ? k - 1 : 0;
            const km2 = k > 1 ? k - 2 : 0;
            const kp = k + 1 < nz ? k + 1 : nz - 1;
            const kp2 = k + 2 < nz ? k + 2 : nz - 1;

            let fB = 0;
            if (k > 0) {
              const wB = S.w[n];
              const cD = f[n - nXY];
              const cD2 = f[n - nXY * (k - km2)];
              const cU = k < nz - 1 ? f[n + nXY] : c0;
              let fq = wB > 0 ? (6 * cD + 3 * c0 - cD2) / 8 : (6 * c0 + 3 * cD - cU) / 8;
              const lo = cD < c0 ? cD : c0, hi = cD < c0 ? c0 : cD;
              if (fq < lo) fq = lo; else if (fq > hi) fq = hi;
              fB = rhoW[k] * wB * fq;
            }

            let fT = 0;
            if (k < nz - 1) {
              const wT = S.w[n + nXY];
              const cU = f[n + nXY];
              const cU2 = f[n + nXY * (kp2 - k)];
              const cD = k > 0 ? f[n - nXY] : c0;
              let fq = wT > 0 ? (6 * c0 + 3 * cU - cD) / 8 : (6 * cU + 3 * c0 - cU2) / 8;
              const lo = c0 < cU ? c0 : cU, hi = c0 < cU ? cU : c0;
              if (fq < lo) fq = lo; else if (fq > hi) fq = hi;
              fT = rhoW[k + 1] * wT * fq;
            }

            tf[n] -= (fT - fB) / (dz * rho0[k]);
          }
        }
      }
    }

    // ---- u 动量 ----
    for (let k = 0; k < nz; k++) {
      const off = nXY * k;
      const km = k > 0 ? k - 1 : 0;
      const km2 = k > 1 ? k - 2 : 0;
      const kp = k + 1 < nz ? k + 1 : nz - 1;
      const kp2 = k + 2 < nz ? k + 2 : nz - 1;
      for (let j = 0; j < ny; j++) {
        const jp = j + 1 === ny ? 0 : j + 1;
        const jm = j === 0 ? ny - 1 : j - 1;
        const row = off + nx * j;
        for (let i = 0; i < nx; i++) {
          const n = row + i;
          const ua = S.u[n];
          const va = S.v[n];
          const wa = 0.5 * (S.w[n] + S.w[n + nXY]);
          T.tu[n] = -(ua * upX(S.u, n, i, ua, row, nx, dx)
            + va * upY(S.u, n, i, j, k, va, nx, ny, nXY, dy)
            + wa * upZ(S.u, n, k, wa, nXY, nz, dz, km, km2, kp, kp2));
          void jp; void jm;
        }
      }
    }

    // ---- v 动量 ----
    for (let k = 0; k < nz; k++) {
      const off = nXY * k;
      const km = k > 0 ? k - 1 : 0;
      const km2 = k > 1 ? k - 2 : 0;
      const kp = k + 1 < nz ? k + 1 : nz - 1;
      const kp2 = k + 2 < nz ? k + 2 : nz - 1;
      for (let j = 0; j < ny; j++) {
        const row = off + nx * j;
        for (let i = 0; i < nx; i++) {
          const n = row + i;
          const ua = S.u[n];
          const va = S.v[n];
          const wa = 0.5 * (S.w[n] + S.w[n + nXY]);
          T.tv[n] = -(ua * upX(S.v, n, i, ua, row, nx, dx)
            + va * upY(S.v, n, i, j, k, va, nx, ny, nXY, dy)
            + wa * upZ(S.v, n, k, wa, nXY, nz, dz, km, km2, kp, kp2));
        }
      }
    }

    // ---- w 动量 + 浮力（w 位于面层 k=1..nz-1）----
    for (let k = 1; k < nz; k++) {
      const off = nXY * k;
      const km = k - 1, km2 = k > 1 ? k - 2 : 0;
      const kp = k + 1 <= nz ? k + 1 : nz;
      const kp2 = k + 2 <= nz ? k + 2 : nz;
      for (let j = 0; j < ny; j++) {
        const row = off + nx * j;
        for (let i = 0; i < nx; i++) {
          const n = row + i;
          const wa = S.w[n];
          const ua = 0.5 * (S.u[n] + S.u[n - nXY]);
          const va = 0.5 * (S.v[n] + S.v[n - nXY]);

          const adv = -(ua * upX(S.w, n, i, ua, row, nx, dx)
            + va * upY(S.w, n, i, j, k, va, nx, ny, nXY, dy)
            + wa * upZ(S.w, n, k, wa, nXY, nz, dz, km, km2, kp, kp2));

          const thA = 0.5 * (S.th[n] + S.th[n - nXY]);
          const qvA = 0.5 * (S.qv[n] + S.qv[n - nXY]);
          const qcA = 0.5 * (S.qc[n] + S.qc[n - nXY]);
          const qrA = 0.5 * (S.qr[n] + S.qr[n - nXY]);
          const qiA = 0.5 * (S.qi[n] + S.qi[n - nXY]);
          const qsA = 0.5 * (S.qs[n] + S.qs[n - nXY]);
          const qgA = 0.5 * (S.qg[n] + S.qg[n - nXY]);
          const th0A = 0.5 * (theta0[k] + theta0[k - 1]);
          const qv0A = 0.5 * (qv0[k] + qv0[k - 1]);

          const b = gg * (thA / th0A - 1 + 0.61 * (qvA - qv0A) - qcA - qrA - qiA - qsA - qgA);
          T.tw[n] = adv + b;
        }
      }
    }

    this.turb.apply(S, base, T);
    this.rad.apply(S, base, T);
  }

  /* --------------------------- 数值滤波 --------------------------- */

  filterField3(arr, nLev, eps) {
    const g = this.grid;
    const { nx, ny } = g;
    const nXY = nx * ny;
    const tmp = this.tmp3;
    const nTot = nXY * nLev;

    for (let n = 0; n < nTot; n++) tmp[n] = arr[n];
    for (let k = 0; k < nLev; k++) {
      const off = nXY * k;
      for (let j = 0; j < ny; j++) {
        const row = off + nx * j;
        for (let i = 0; i < nx; i++) {
          const ip1 = i + 1 === nx ? 0 : i + 1;
          const ip2 = (i + 2) % nx;
          const im1 = i === 0 ? nx - 1 : i - 1;
          const im2 = (i - 2 + nx) % nx;
          const st = tmp[row + ip2] - 4 * tmp[row + ip1] + 6 * tmp[row + i]
            - 4 * tmp[row + im1] + tmp[row + im2];
          arr[row + i] = tmp[row + i] - eps * st / 16;
        }
      }
    }

    for (let n = 0; n < nTot; n++) tmp[n] = arr[n];
    for (let k = 0; k < nLev; k++) {
      const off = nXY * k;
      for (let j = 0; j < ny; j++) {
        const jp1 = j + 1 === ny ? 0 : j + 1;
        const jp2 = (j + 2) % ny;
        const jm1 = j === 0 ? ny - 1 : j - 1;
        const jm2 = (j - 2 + ny) % ny;
        const rowP1 = off + nx * jp1;
        const rowP2 = off + nx * jp2;
        const rowM1 = off + nx * jm1;
        const rowM2 = off + nx * jm2;
        const row = off + nx * j;
        for (let i = 0; i < nx; i++) {
          const st = tmp[rowP2 + i] - 4 * tmp[rowP1 + i] + 6 * tmp[row + i]
            - 4 * tmp[rowM1 + i] + tmp[rowM2 + i];
          arr[row + i] = tmp[row + i] - eps * st / 16;
        }
      }
    }

    for (let n = 0; n < nTot; n++) tmp[n] = arr[n];
    for (let k = 2; k < nLev - 2; k++) {
      const off = nXY * k;
      for (let j = 0; j < ny; j++) {
        const row = off + nx * j;
        for (let i = 0; i < nx; i++) {
          const st = tmp[row + nXY * 2 + i] - 4 * tmp[row + nXY + i] + 6 * tmp[row + i]
            - 4 * tmp[row - nXY + i] + tmp[row - nXY * 2 + i];
          arr[row + i] = tmp[row + i] - eps * st / 16;
        }
      }
    }
  }

  applyFilter(eps) {
    const g = this.grid;
    if (!(eps > 0)) return;
    this.filterField3(this.S.u, g.nz, eps);
    this.filterField3(this.S.v, g.nz, eps);
    this.filterField3(this.S.w, g.nz + 1, eps);
    for (let n = 0; n < g.nXY; n++) {
      this.S.w[n] = 0;
      this.S.w[n + g.nXY * g.nz] = 0;
    }
    if (this.terrain) {
      resetSolidState(this.S, this.base, g, this.terrain);
      maskVelocities(this.S, g, this.terrain);
    }
    const res = this.solver.project(this.S.u, this.S.v, this.S.w, this.lastDt > 0 ? this.lastDt : 1, this.pi, this.base);
    this.ibmIterations = res.iterations ?? this.ibmIterations;
    if (this.terrain) maskVelocities(this.S, g, this.terrain);
  }

  /* --------------------------- 时间推进 --------------------------- */

  stage(dt, src, dst) {
    this.computeTendencies(src, this.T);
    const g = this.grid;

    for (let n = 0; n < g.nCell; n++) {
      dst.u[n] = src.u[n] + dt * this.T.tu[n];
      dst.v[n] = src.v[n] + dt * this.T.tv[n];
      dst.th[n] = src.th[n] + dt * this.T.tth[n];
      dst.qv[n] = src.qv[n] + dt * this.T.tqv[n];
      dst.qc[n] = src.qc[n] + dt * this.T.tqc[n];
      dst.qr[n] = src.qr[n] + dt * this.T.tqr[n];
      dst.qi[n] = src.qi[n] + dt * this.T.tqi[n];
      dst.qs[n] = src.qs[n] + dt * this.T.tqs[n];
      dst.qg[n] = src.qg[n] + dt * this.T.tqg[n];
    }
    for (let n = 0; n < g.nW; n++) dst.w[n] = src.w[n] + dt * this.T.tw[n];

    if (this.terrain) {
      resetSolidState(dst, this.base, g, this.terrain);
      maskVelocities(dst, g, this.terrain);
    }
    const res = this.solver.project(dst.u, dst.v, dst.w, dt, this.pi, this.base);
    this.ibmIterations = res.iterations ?? this.ibmIterations;
    if (this.terrain) maskVelocities(dst, g, this.terrain);
    this.divResidual = res.divAfter;
  }

  step(dt) {
    const S0 = this.S, A = this.B1, B = this.B2, S2 = this.S2;

    this.stage(dt, S0, A);
    this.stage(dt, A, B);
    combineState3(S2, S0, 0.75, B, 0.25);
    this.stage(dt, S2, A);
    combineState3(S0, S0, 1 / 3, A, 2 / 3);
    this.applyFilter(this.cfg.filter);

    const rep = repairMoisture3(S0);
    this.repairAmount = rep.repaired;
    this.repairLost = rep.lost;

    if (this.cfg.micro.enable) {
      const flux = this.micro.apply(S0, this.base, dt);
      this.precipRate = flux * 3600;
      this.precipAccum += flux * dt;
    }
    this.waterInput += this.turb.waterRate * dt;
    clampMoisture3(S0);
    if (this.terrain) {
      resetSolidState(S0, this.base, this.grid, this.terrain);
      maskVelocities(S0, this.grid, this.terrain);
    }

    this.time += dt;
    this.lastDt = dt;
    let vmax = Math.max(maxAbs3(S0.u), maxAbs3(S0.v), maxAbs3(S0.w));
    if (vmax > 300) vmax = 300;
    this.maxSpeed = vmax > 6 ? vmax : 6;
  }

  suggestDt() {
    const g = this.grid;
    const dtAdv = this.cfg.cfl * Math.min(g.dx, g.dy, g.dz) / this.maxSpeed;
    const Lmin = Math.min(g.dx, g.dy, g.dz);
    const kEff = Math.min(this.cfg.turb.kMax, 0.05 * Lmin * Lmin);
    const invSum = 1 / (g.dx * g.dx) + 1 / (g.dy * g.dy) + 1 / (g.dz * g.dz);
    const dtDif = 0.5 / (kEff * invSum);
    let dt = Math.min(dtAdv, dtDif, this.cfg.dtMax);
    if (dt < 0.05) dt = 0.05;
    return dt;
  }

  advance(simSeconds, budgetMs = 0) {
    const t0 = nowMs();
    this.maxSpeed = Math.max(maxAbs3(this.S.u), maxAbs3(this.S.v), maxAbs3(this.S.w), 6);
    const dt = this.suggestDt();
    let remaining = simSeconds;
    let steps = 0;

    while (remaining > 1e-6) {
      if (budgetMs > 0 && steps > 0 && nowMs() - t0 >= budgetMs) break;
      const d = dt < remaining ? dt : remaining;
      this.step(d);
      remaining -= d;
      steps++;
      if (steps > 100000) break;
    }

    return { steps, dt, simAdvanced: simSeconds - remaining, ms: nowMs() - t0 };
  }

  /* --------------------------- 采样 --------------------------- */

  probe(x, y, z) {
    const g = this.grid, S = this.S;
    const fi = x / g.dx - 0.5;
    const fj = y / g.dy - 0.5;
    const fk = z / g.dz - 0.5;
    const i0 = Math.max(0, Math.min(g.nx - 1, Math.floor(fi)));
    const j0 = Math.max(0, Math.min(g.ny - 1, Math.floor(fj)));
    const k0 = Math.max(0, Math.min(g.nz - 1, Math.floor(fk)));
    const i1 = (i0 + 1) % g.nx;
    const j1 = (j0 + 1) % g.ny;
    const k1 = Math.min(g.nz - 1, k0 + 1);
    const tx = Math.max(0, Math.min(1, fi - i0));
    const ty = Math.max(0, Math.min(1, fj - j0));
    const tz = Math.max(0, Math.min(1, fk - k0));

    const lerp = (arr) => {
      const w00 = arr[i0 + g.nx * (j0 + g.ny * k0)], w10 = arr[i1 + g.nx * (j0 + g.ny * k0)];
      const w01 = arr[i0 + g.nx * (j1 + g.ny * k0)], w11 = arr[i1 + g.nx * (j1 + g.ny * k0)];
      const w0 = arr[i0 + g.nx * (j0 + g.ny * k1)], w1 = arr[i1 + g.nx * (j0 + g.ny * k1)];
      const w2 = arr[i0 + g.nx * (j1 + g.ny * k1)], w3 = arr[i1 + g.nx * (j1 + g.ny * k1)];
      const lo = (w00 * (1 - tx) + w10 * tx) * (1 - ty) + (w01 * (1 - tx) + w11 * tx) * ty;
      const hi = (w0 * (1 - tx) + w1 * tx) * (1 - ty) + (w2 * (1 - tx) + w3 * tx) * ty;
      return lo * (1 - tz) + hi * tz;
    };

    return {
      u: lerp(S.u), v: lerp(S.v), w: lerp(S.w), th: lerp(S.th),
      qv: lerp(S.qv), qc: lerp(S.qc), qr: lerp(S.qr), qi: lerp(S.qi), qs: lerp(S.qs), qg: lerp(S.qg),
      T: lerp(S.th) * this.base.ex[k0],
    };
  }

  addBubble({ x, y, z, amplitude = 2.5, radius = 2000 }) {
    const g = this.grid, S = this.S;
    const yc = y ?? g.Ly * 0.5;
    // 有地形时，把落在地下的热泡抬到该列最低的流体层之上
    let zc = z;
    if (this.terrain) {
      const col = (Math.rint(x / g.dx) + g.nx) % g.nx
        + g.nx * ((Math.rint(yc / g.dy) + g.ny) % g.ny);
      const kS = this.terrain.surfaceK[col];
      const zSurf = g.zc(kS);
      if (zc - radius < zSurf) zc = zSurf + radius * 0.6;
    }
    for (let k = 0; k < g.nz; k++) {
      const dzr = (g.zc(k) - z) / radius;
      if (dzr > 1 || dzr < -1) continue;
      const off = g.nXY * k;
      for (let j = 0; j < g.ny; j++) {
        const dyr = (g.yc(j) - yc) / radius;
        if (dyr > 1 || dyr < -1) continue;
        for (let i = 0; i < g.nx; i++) {
          let dxAbs = Math.abs(g.xc(i) - x);
          if (dxAbs > g.Lx * 0.5) dxAbs = g.Lx - dxAbs;
          const dxr = dxAbs / radius;
          const rr = Math.sqrt(dxr * dxr + dyr * dyr + dzr * dzr);
          if (rr < 1) S.th[off + g.nx * j + i] += amplitude * 0.5 * (1 + Math.cos(Math.PI * rr));
        }
      }
    }
  }

  snapshot() {
    const out = {};
    for (const f of FIELD_NAMES3) out[f] = this.S[f].slice();
    out.time = this.time;
    out.precipAccum = this.precipAccum;
    return out;
  }

  restore(s) {
    for (const f of FIELD_NAMES3) this.S[f].set(s[f]);
    this.time = s.time ?? 0;
    this.precipAccum = s.precipAccum ?? 0;
  }
}

/* -------------------- 三阶迎风差分（周期 x、y；z 钳位） -------------------- */

function upX(f, n, i, v, row, nx, dx) {
  const ip = i + 1 === nx ? 0 : i + 1;
  const im = i === 0 ? nx - 1 : i - 1;
  if (v > 0) {
    const im2 = i <= 1 ? (i - 2 + nx) % nx : i - 2;
    return (2 * f[row + ip] + 3 * f[n] - 6 * f[row + im] + f[row + im2]) / (6 * dx);
  }
  const ip2 = (i + 2) % nx;
  return -(2 * f[row + im] + 3 * f[n] - 6 * f[row + ip] + f[row + ip2]) / (6 * dx);
}

function upY(f, n, i, j, k, v, nx, ny, nXY, dy) {
  const layer = nXY * k + i;                 // 同一层、同一 i 的起点；行步长是 nx
  const jp = j + 1 === ny ? 0 : j + 1;
  const jm = j === 0 ? ny - 1 : j - 1;
  if (v > 0) {
    const jm2 = j <= 1 ? (j - 2 + ny) % ny : j - 2;
    return (2 * f[layer + nx * jp] + 3 * f[n] - 6 * f[layer + nx * jm] + f[layer + nx * jm2]) / (6 * dy);
  }
  const jp2 = (j + 2) % ny;
  return -(2 * f[layer + nx * jm] + 3 * f[n] - 6 * f[layer + nx * jp] + f[layer + nx * jp2]) / (6 * dy);
}

function upZ(f, n, k, v, nXY, nz, dz, km, km2, kp, kp2) {
  if (v > 0) {
    return (2 * f[n + nXY * (kp - k)] + 3 * f[n] - 6 * f[n + nXY * (km - k)] + f[n + nXY * (km2 - k)]) / (6 * dz);
  }
  return -(2 * f[n + nXY * (km - k)] + 3 * f[n] - 6 * f[n + nXY * (kp - k)] + f[n + nXY * (kp2 - k)]) / (6 * dz);
}
