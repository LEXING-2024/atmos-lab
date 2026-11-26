/**
 * 二维非静力伪不可压大气模式（x–z 剖面）。
 *
 * 这是全仓库唯一推进时间的地方。
 * 所有网格尺寸来自 Grid 对象，代码中不存在网格常量。
 */

import { makeBase } from './base.js';
import { allocState, copyState, combineState, clampMoisture, repairMoisture, maxAbs } from './state.js';
import { createPoissonSolver } from '../core/poisson.js';
import { createMicrophysics } from '../physics/micro.js';
import { createTurbulence } from '../physics/turb.js';
import { createRadiation } from '../physics/rad.js';
import { remap2D } from '../core/remap.js';
import { CONST } from '../physics/thermo.js';

const nowMs = () => (
  typeof performance !== 'undefined' && performance.now ? performance.now() : Date.now()
);

export const DEFAULT_CONFIG = {
  dtMax: 8,
  cfl: 0.4,
  filter: 0.02,
  micro: { enable: true, k1: 1e-3, qc0: 0.6e-3, k2: 2.2, evapTau: 1200 },
  turb: {
    cs: 0.15, kMin: 1.0, kMax: 1500, cd: 0.002,
    surface: { sensible: 0, latent: 0, xCenter: null, xWidth: 8000 },
  },
  rad: { tauRad: 8 * 3600, spongeMax: 1 / 300, spongeDepth: 0.7 },
};

/* ---------------------------- 三阶迎风差分 ---------------------------- */

function ddx(f, i, k, v, nx, dx) {
  const o = nx * k;
  if (v > 0) {
    const im1 = i === 0 ? nx - 1 : i - 1;
    const im2 = i <= 1 ? (i - 2 + nx) % nx : i - 2;
    const ip1 = i + 1 === nx ? 0 : i + 1;
    return (2 * f[ip1 + o] + 3 * f[i + o] - 6 * f[im1 + o] + f[im2 + o]) / (6 * dx);
  }
  const ip1 = i + 1 === nx ? 0 : i + 1;
  const ip2 = (i + 2) % nx;
  const im1 = i === 0 ? nx - 1 : i - 1;
  return -(2 * f[im1 + o] + 3 * f[i + o] - 6 * f[ip1 + o] + f[ip2 + o]) / (6 * dx);
}

function ddz(f, i, k, v, nx, kMax, dz) {
  const o = nx * k;
  if (v > 0) {
    const km1 = k > 0 ? k - 1 : 0;
    const km2 = k > 1 ? k - 2 : 0;
    const kp1 = k + 1 <= kMax ? k + 1 : kMax;
    return (2 * f[i + nx * kp1] + 3 * f[i + o] - 6 * f[i + nx * km1] + f[i + nx * km2]) / (6 * dz);
  }
  const kp1 = k + 1 <= kMax ? k + 1 : kMax;
  const kp2 = k + 2 <= kMax ? k + 2 : kMax;
  const km1 = k > 0 ? k - 1 : 0;
  return -(2 * f[i + nx * km1] + 3 * f[i + o] - 6 * f[i + nx * kp1] + f[i + nx * kp2]) / (6 * dz);
}

/* ------------------------------- 模型 ------------------------------- */

export class AnelasticModel {
 constructor(grid, opts = {}) {
    this.grid = grid;
    this.opts = opts;
    this.cfg = {
      ...DEFAULT_CONFIG,
      ...opts,
      micro: { ...DEFAULT_CONFIG.micro, ...(opts.micro || {}) },
      turb: {
        ...DEFAULT_CONFIG.turb,
        ...(opts.turb || {}),
        surface: { ...DEFAULT_CONFIG.turb.surface, ...((opts.turb || {}).surface || {}) },
      },
      rad: { ...DEFAULT_CONFIG.rad, ...(opts.rad || {}) },
    };

    this.base = makeBase(grid, opts.base || {});
    this.solver = createPoissonSolver(grid);
    this.micro = createMicrophysics(grid, this.cfg.micro);
    this.turb = createTurbulence(grid, this.cfg.turb);
    this.rad = createRadiation(grid, this.cfg.rad);

    const nCell = grid.nCell, nW = grid.nW;
    this.S = allocState(grid);
    this.B1 = allocState(grid);
    this.B2 = allocState(grid);
    this.S2 = allocState(grid);

    this.T = {
      tu: new Float64Array(nCell),
      tw: new Float64Array(nW),
      tth: new Float64Array(nCell),
      tqv: new Float64Array(nCell),
      tqc: new Float64Array(nCell),
      tqr: new Float64Array(nCell),
    };
    this.pi = new Float64Array(nCell);
    this.Fx = new Float64Array(grid.nx);
    this.tmpArr = new Float64Array(Math.max(grid.nCell, grid.nW));
    this.uc = new Float64Array(nCell);
    this.wc = new Float64Array(nCell);

    this.time = 0;
    this.maxSpeed = 8;
    this.precipRate = 0;
    this.precipAccum = 0;
    this.waterInput = 0;
    this.divResidual = 0;
    this.cflNow = 0;
    this.lastDt = 0;
    this.initMass = 0;
  }

  /* --------------------------- 初始化 --------------------------- */

  initAtmosphere(o = {}) {
    const g = this.grid, base = this.base;
    const amp = o.amplitude ?? 2.5;
    const radius = o.radius ?? 2000;
    const zc = o.zCenter ?? 1500;
    const xc = o.xCenter ?? g.Lx * 0.5;
    const noise = o.noise ?? 0.02;
    const seed = o.seed ?? 20261127;

    const S = this.S;
    for (let k = 0; k < g.nz; k++) {
      const off = g.nx * k;
      for (let i = 0; i < g.nx; i++) {
        S.th[off + i] = base.theta0[k];
        S.qv[off + i] = base.qv0[k];
        S.qc[off + i] = 0;
        S.qr[off + i] = 0;
      }
    }
    S.u.fill(0);
    S.w.fill(0);
    this.pi.fill(0);

    // 热泡（余弦钟）
    for (let k = 0; k < g.nz; k++) {
      const z = g.zc(k), off = g.nx * k;
      for (let i = 0; i < g.nx; i++) {
        const x = g.xc(i);
        const r = Math.sqrt(((x - xc) / radius) ** 2 + ((z - zc) / radius) ** 2);
        if (r < 1) S.th[off + i] += amp * 0.5 * (1 + Math.cos(Math.PI * r));
      }
    }

    // 低层微小噪声，打破对称（固定种子，保证可复现）
    let s = seed >>> 0;
    const rnd = () => {
      s = (s * 1664525 + 1013904223) >>> 0;
      return s / 4294967296 - 0.5;
    };
    for (let k = 0; k < g.nz && g.zc(k) < 3000; k++) {
      const off = g.nx * k;
      for (let i = 0; i < g.nx; i++) S.th[off + i] += noise * rnd();
    }

    this.time = 0;
    this.precipAccum = 0;
    this.precipRate = 0;
    this.waterInput = 0;
    this.initMass = this.totalWater();
  }

  /* --------------------------- 诊断 --------------------------- */

  totalWater() {
    const g = this.grid, base = this.base, S = this.S;
    let sum = 0;
    for (let k = 0; k < g.nz; k++) {
      const off = g.nx * k, rk = base.rho0[k];
      for (let i = 0; i < g.nx; i++) {
        const n = off + i;
        sum += rk * (S.qv[n] + S.qc[n] + S.qr[n]);
      }
    }
    return sum * g.dx * g.dz;
  }

  diagnostics() {
    const g = this.grid, S = this.S;
    let maxW = 0, maxU = 0, cloudTop = 0;
    for (let k = 0; k < g.nz; k++) {
      const off = g.nx * k;
      for (let i = 0; i < g.nx; i++) {
        const n = off + i;
        const w = Math.abs(S.w[n] + S.w[n + g.nx]) * 0.5;
        if (w > maxW) maxW = w;
        const u = Math.abs(S.u[n]);
        if (u > maxU) maxU = u;
        if (S.qc[n] + S.qr[n] > 1e-5 && g.zc(k) > cloudTop) cloudTop = g.zc(k);
      }
    }
    const cfl = Math.max(
      maxU * this.lastDt / g.dx,
      maxW * this.lastDt / g.dz,
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
    const { nx, nz, dx, dz } = g;
    const nCell = nx * nz;

    T.tu.fill(0); T.tw.fill(0);
    T.tth.fill(0); T.tqv.fill(0); T.tqc.fill(0); T.tqr.fill(0);

    const theta0 = base.theta0, qv0 = base.qv0;
    const rho0 = base.rho0, rhoW = base.rhoW;
    const gg = CONST.g;

    // 标量平流：通量形式 + 单调限制器（守恒、正定）
    const Fx = this.Fx;
    const scalars = [S.th, S.qv, S.qc, S.qr];
    const tscalars = [T.tth, T.tqv, T.tqc, T.tqr];

    for (let s = 0; s < 4; s++) {
      const f = scalars[s], tf = tscalars[s];

      // x 方向
      for (let k = 0; k < nz; k++) {
        const off = nx * k;
        for (let i = 0; i < nx; i++) {
          const n = off + i;
          const v = S.u[n];
          const im1 = i === 0 ? nx - 1 : i - 1;
          const im2 = i <= 1 ? (i - 2 + nx) % nx : i - 2;
          const ip1 = i + 1 === nx ? 0 : i + 1;

          const a = f[im1 + off], b = f[i + off];
          let fq = v > 0
            ? (6 * a + 3 * b - f[im2 + off]) / 8
            : (6 * b + 3 * a - f[ip1 + off]) / 8;
          const lo = a < b ? a : b, hi = a < b ? b : a;
          if (fq < lo) fq = lo; else if (fq > hi) fq = hi;
          Fx[i] = rho0[k] * v * fq;
        }
        const rk = rho0[k];
        for (let i = 0; i < nx; i++) {
          const ip = i + 1 === nx ? 0 : i + 1;
          tf[off + i] = -(Fx[ip] - Fx[i]) / (dx * rk);
        }
      }

      // z 方向（w 向上为正，底面/顶面通量分别取）
      for (let k = 0; k < nz; k++) {
        const off = nx * k;
        for (let i = 0; i < nx; i++) {
          const n = off + i;
          const c0 = f[n];
          let fB = 0;
          if (k > 0) {
            const wB = S.w[n];
            const cD = f[n - nx];
            const cD2 = k > 1 ? f[n - 2 * nx] : cD;
            const cU = k < nz - 1 ? f[n + nx] : c0;
            let fq = wB > 0
              ? (6 * cD + 3 * c0 - cD2) / 8
              : (6 * c0 + 3 * cD - cU) / 8;
            const lo = cD < c0 ? cD : c0, hi = cD < c0 ? c0 : cD;
            if (fq < lo) fq = lo; else if (fq > hi) fq = hi;
            fB = rhoW[k] * wB * fq;
          }
          let fT = 0;
          if (k < nz - 1) {
            const wT = S.w[n + nx];
            const cU = f[n + nx];
            const cU2 = k + 2 < nz ? f[n + 2 * nx] : cU;
            const cD = k > 0 ? f[n - nx] : c0;
            let fq = wT > 0
              ? (6 * c0 + 3 * cU - cD) / 8
              : (6 * cU + 3 * c0 - cU2) / 8;
            const lo = c0 < cU ? c0 : cU, hi = c0 < cU ? cU : c0;
            if (fq < lo) fq = lo; else if (fq > hi) fq = hi;
            fT = rhoW[k + 1] * wT * fq;
          }
          tf[n] -= (fT - fB) / (dz * rho0[k]);
        }
      }
    }

    // u 动量
    for (let k = 0; k < nz; k++) {
      const off = nx * k;
      for (let i = 0; i < nx; i++) {
        const n = off + i;
        const ua = S.u[n];
        const im = i === 0 ? nx - 1 : i - 1;
        const wa = 0.5 * (S.w[im + off] + S.w[n]);
        T.tu[n] = -(ua * ddx(S.u, i, k, ua, nx, dx) + wa * ddz(S.u, i, k, wa, nx, nz - 1, dz));
      }
    }

    // w 动量 + 浮力
    for (let k = 1; k < nz; k++) {
      const off = nx * k;
      for (let i = 0; i < nx; i++) {
        const n = off + i;
        const wa = S.w[n];
        const ua = 0.5 * (S.u[n] + S.u[n - nx]);

        const adv = -(ua * ddx(S.w, i, k, ua, nx, dx) + wa * ddz(S.w, i, k, wa, nx, nz, dz));

        const thA = 0.5 * (S.th[n] + S.th[n - nx]);
        const qvA = 0.5 * (S.qv[n] + S.qv[n - nx]);
        const qcA = 0.5 * (S.qc[n] + S.qc[n - nx]);
        const qrA = 0.5 * (S.qr[n] + S.qr[n - nx]);
        const th0A = 0.5 * (theta0[k] + theta0[k - 1]);
        const qv0A = 0.5 * (qv0[k] + qv0[k - 1]);

        const b = gg * (thA / th0A - 1 + 0.61 * (qvA - qv0A) - qcA - qrA);
        T.tw[n] = adv + b;
      }
    }

    // 参数化
    this.turb.apply(S, base, T);
    this.rad.apply(S, base, T);
    void nCell;
  }

  /* --------------------------- 数值滤波 --------------------------- */

  /**
   * 四阶（Δ⁴）滤波，抑制 2Δ 网格噪声。
   * 动量平流采用中心型离散时，近零速度下会出现奇偶失稳，此项为必需的稳定化手段。
   */
  filterField(arr, nLev, eps) {
    const nx = this.grid.nx;
    const tmp = this.tmpArr;

    for (let n = 0; n < nx * nLev; n++) tmp[n] = arr[n];

    for (let k = 0; k < nLev; k++) {
      const off = nx * k;
      for (let i = 0; i < nx; i++) {
        const ip1 = i + 1 === nx ? 0 : i + 1;
        const ip2 = (i + 2) % nx;
        const im1 = i === 0 ? nx - 1 : i - 1;
        const im2 = (i - 2 + nx) % nx;
        const st = tmp[ip2 + off] - 4 * tmp[ip1 + off] + 6 * tmp[i + off]
          - 4 * tmp[im1 + off] + tmp[im2 + off];
        arr[i + off] = tmp[i + off] - eps * st / 16;
      }
    }

    for (let n = 0; n < nx * nLev; n++) tmp[n] = arr[n];

    // z 向只作用于内部层：边界层钳位模板不满足 Σ系数=0（不守恒）且对线性剖面有偏置
    for (let k = 2; k < nLev - 2; k++) {
      const off = nx * k;
      for (let i = 0; i < nx; i++) {
        const st = tmp[i + nx * (k + 2)] - 4 * tmp[i + nx * (k + 1)] + 6 * tmp[i + off]
          - 4 * tmp[i + nx * (k - 1)] + tmp[i + nx * (k - 2)];
        arr[i + off] = tmp[i + off] - eps * st / 16;
      }
    }
  }

  applyFilter(eps) {
    const g = this.grid;
    const S = this.S;
    if (!(eps > 0)) return;
    // 只滤速度：标量由单调平流格式自行耗散。
    // 滤波会破坏 ∇·(ρ₀u)=0，必须重新投影，否则通量形式平流不再与连续方程一致，
    // 会逐步伪造位温倾向与水量收支。
    this.filterField(S.u, g.nz, eps);
    this.filterField(S.w, g.nz + 1, eps);
    for (let i = 0; i < g.nx; i++) { S.w[i] = 0; S.w[i + g.nx * g.nz] = 0; }
    this.solver.project(S.u, S.w, this.lastDt > 0 ? this.lastDt : 1, this.pi, this.base);
  }

  /* --------------------------- 时间推进 --------------------------- */

  stage(dt, src, dst) {
    this.computeTendencies(src, this.T);
    const g = this.grid;

    for (let n = 0; n < g.nCell; n++) {
      dst.u[n] = src.u[n] + dt * this.T.tu[n];
      dst.th[n] = src.th[n] + dt * this.T.tth[n];
      dst.qv[n] = src.qv[n] + dt * this.T.tqv[n];
      dst.qc[n] = src.qc[n] + dt * this.T.tqc[n];
      dst.qr[n] = src.qr[n] + dt * this.T.tqr[n];
    }
    for (let n = 0; n < g.nW; n++) dst.w[n] = src.w[n] + dt * this.T.tw[n];

    const res = this.solver.project(dst.u, dst.w, dt, this.pi, this.base);
    this.divResidual = res.divAfter;
  }

  step(dt) {
    const S0 = this.S, A = this.B1, B = this.B2, S2 = this.S2;

    this.stage(dt, S0, A);
    this.stage(dt, A, B);
    combineState(S2, S0, 0.75, B, 0.25);
    this.stage(dt, S2, A);
    combineState(S0, S0, 1 / 3, A, 2 / 3);
    this.applyFilter(this.cfg.filter);

    // 输运后做守恒式非负修补（负值置零会凭空造水）
    const rep = repairMoisture(S0, this.grid.nx, this.grid.nz);
    this.repairAmount = rep.repaired;
    this.repairLost = rep.lost;

    if (this.cfg.micro.enable) {
      const flux = this.micro.apply(S0, this.base, dt);   // kg/(m²·s)
      this.precipRate = flux * 3600;                      // mm/h
      this.precipAccum += flux * dt;                      // mm
    }
    this.waterInput += this.turb.waterRate * dt;           // kg/m（每单位宽度）
    clampMoisture(S0);

    this.time += dt;
    this.lastDt = dt;
    let vmax = Math.max(maxAbs(S0.u), maxAbs(S0.w));
    if (vmax > 300) vmax = 300;
    this.maxSpeed = vmax > 6 ? vmax : 6;
  }

  suggestDt() {
    const g = this.grid;
    const dtAdv = this.cfg.cfl * Math.min(g.dx, g.dz) / this.maxSpeed;
    const kEff = Math.min(this.cfg.turb.kMax, 0.05 * g.dz * g.dz);
    const dtDif = 0.5 * g.dz * g.dz / (2 * kEff);
    let dt = Math.min(dtAdv, dtDif, this.cfg.dtMax);
    if (dt < 0.05) dt = 0.05;
    return dt;
  }

  /**
   * 在墙钟预算内推进指定模拟时长。
   * 固定 dt（一次调用内不变），预算不足则提前返回。
   */
  advance(simSeconds, budgetMs = 0) {
    const t0 = nowMs();
    this.maxSpeed = Math.max(maxAbs(this.S.u), maxAbs(this.S.w), 6);
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

    return {
      steps,
      dt,
      simAdvanced: simSeconds - remaining,
      ms: nowMs() - t0,
    };
  }

  /* --------------------------- 采样 --------------------------- */

  probe(x, z) {
    const g = this.grid, S = this.S;
    const fi = x / g.dx - 0.5;
    const fk = z / g.dz - 0.5;
    const i0 = Math.max(0, Math.min(g.nx - 1, Math.floor(fi)));
    const k0 = Math.max(0, Math.min(g.nz - 1, Math.floor(fk)));
    const i1 = (i0 + 1) % g.nx;
    const k1 = Math.min(g.nz - 1, k0 + 1);
    const tx = Math.max(0, Math.min(1, fi - i0));
    const tz = Math.max(0, Math.min(1, fk - k0));

    const lerp = (arr) => {
      const a = arr[i0 + g.nx * k0], b = arr[i1 + g.nx * k0];
      const c = arr[i0 + g.nx * k1], d = arr[i1 + g.nx * k1];
      return (a * (1 - tx) + b * tx) * (1 - tz) + (c * (1 - tx) + d * tx) * tz;
    };

    return {
      u: lerp(S.u), w: lerp(S.w), th: lerp(S.th),
      qv: lerp(S.qv), qc: lerp(S.qc), qr: lerp(S.qr),
      T: lerp(S.th) * this.base.ex[k0],
    };
  }

  /** 在指定位置叠加一个暖泡（交互注入用） */
  addBubble({ x, z, amplitude = 2.5, radius = 2000 }) {
    const g = this.grid, S = this.S;
    const r2 = radius * radius;
    for (let k = 0; k < g.nz; k++) {
      const dz = (g.zc(k) - z) / radius;
      if (dz > 1 || dz < -1) continue;
      const off = g.nx * k;
      for (let i = 0; i < g.nx; i++) {
        let dxAbs = Math.abs(g.xc(i) - x);
        if (dxAbs > g.Lx * 0.5) dxAbs = g.Lx - dxAbs;      // 周期方向取最短
        const dxr = dxAbs / radius;
        const rr = Math.sqrt(dxr * dxr + dz * dz);
        if (rr < 1) S.th[off + i] += amplitude * 0.5 * (1 + Math.cos(Math.PI * rr));
      }
    }
    void r2;
  }

  snapshot() {
    const out = {};
    for (const f of ['u', 'w', 'th', 'qv', 'qc', 'qr']) out[f] = this.S[f].slice();
    out.time = this.time;
    out.precipAccum = this.precipAccum;
    return out;
  }

  restore(s) {
    for (const f of ['u', 'w', 'th', 'qv', 'qc', 'qr']) this.S[f].set(s[f]);
    this.time = s.time ?? 0;
    this.precipAccum = s.precipAccum ?? 0;
  }
}

/**
 * 分辨率切换。
 * mode='conservative' 保留场（守恒重映射）；mode='reset' 按新网格重建初值。
 */
export function resizeModel(oldModel, newGrid, o = {}) {
  const m = new AnelasticModel(newGrid, oldModel.opts);
  if (o.mode === 'conservative') {
    const g0 = oldModel.grid, g1 = newGrid;
    const b0 = oldModel.base, b1 = m.base;
    const tmp = new Float64Array(g1.nCell);
    const weighted = new Float64Array(g0.nCell);

    // 质量加权重映射：重映射 ρq，再除以新网格 ρ。
    // 直接重映射 q 只能守恒 ∫q dz，换网格后 ρ 剖面不同，总水量会跳变。
    for (const f of ['u', 'th', 'qv', 'qc', 'qr']) {
      const src = oldModel.S[f], dst = m.S[f];
      for (let k = 0; k < g0.nz; k++) {
        const off = g0.nx * k, rk = b0.rho0[k];
        for (let i = 0; i < g0.nx; i++) weighted[off + i] = rk * src[off + i];
      }
      remap2D(weighted, g0.nx, g0.nz, g1.nx, g1.nz, tmp);
      for (let k = 0; k < g1.nz; k++) {
        const off = g1.nx * k, rk = b1.rho0[k];
        for (let i = 0; i < g1.nx; i++) dst[off + i] = tmp[off + i] / rk;
      }
    }

    const wSrc = new Float64Array(g0.nW);
    const wOut = new Float64Array(g1.nW);
    for (let k = 0; k <= g0.nz; k++) {
      const off = g0.nx * k, rW = b0.rhoW[k];
      for (let i = 0; i < g0.nx; i++) wSrc[off + i] = rW * oldModel.S.w[off + i];
    }
    remap2D(wSrc, g0.nx, g0.nz + 1, g1.nx, g1.nz + 1, wOut);
    for (let k = 0; k <= g1.nz; k++) {
      const off = g1.nx * k, rW = b1.rhoW[k];
      for (let i = 0; i < g1.nx; i++) m.S.w[off + i] = wOut[off + i] / rW;
    }
    for (let i = 0; i < g1.nx; i++) { m.S.w[i] = 0; m.S.w[i + g1.nx * g1.nz] = 0; }

    m.time = oldModel.time;
    m.precipAccum = oldModel.precipAccum;
    m.waterInput = oldModel.waterInput;
    m.initMass = oldModel.initMass;
  } else {
    m.initAtmosphere(o.init || {});
  }
  return m;
}
