/**
 * 三维地形/障碍物：阶梯化掩膜（staircase mask）。
 *
 * 只做一件事：把每个格点标记成 fluid(0) 或 solid(1)，并给出每列的最低流体层。
 * 具体的边界条件由两个消费者实现：
 *   - core/poisson3_ibm.js  在固壁面上取 Neumann（法向速度=0）
 *   - model/anelastic3.js   每步把固壁面的速度清零（直接强迫 IBM）
 *
 * 支持的地形：
 *   flat      无地形（等价于不带 terrain，仍走 FFT 求解器）
 *   ridge     y 方向无限长的余弦山脊
 *   mountain  两个高斯峰
 *   city      规则排列的矩形街区
 */

/**
 * @param {Grid3} grid
 * @param {object} spec { kind, height, halfWidth, xCenter, yCenter, width, ...
 *                       blocks:{nx,ny,gap,size,height}, peaks:[{x,y,height,radius}] }
 */
export function makeTerrain(grid, spec = {}) {
  const { nx, ny, nz, nXY } = grid;
  const kind = spec.kind ?? 'flat';
  const height = new Float64Array(nXY);

  if (kind === 'ridge') {
    const H = spec.height ?? 1500;
    const halfWidth = spec.halfWidth ?? 2500;
    const xc = spec.xCenter ?? grid.Lx * 0.5;
    for (let j = 0; j < ny; j++) {
      for (let i = 0; i < nx; i++) {
        const dxr = Math.abs(grid.xc(i) - xc);
        if (dxr < halfWidth) height[i + nx * j] = 0.5 * H * (1 + Math.cos(Math.PI * dxr / halfWidth));
      }
    }
  } else if (kind === 'mountain') {
    const peaks = spec.peaks ?? [
      { x: grid.Lx * 0.32, y: grid.Ly * 0.5, height: 2200, radius: 1800 },
      { x: grid.Lx * 0.68, y: grid.Ly * 0.55, height: 1500, radius: 1400 },
    ];
    for (let j = 0; j < ny; j++) {
      for (let i = 0; i < nx; i++) {
        let h = 0;
        for (const p of peaks) {
          const r = Math.hypot(grid.xc(i) - p.x, grid.yc(j) - p.y) / p.radius;
          if (r < 1) h = Math.max(h, p.height * 0.5 * (1 + Math.cos(Math.PI * r)));
        }
        height[i + nx * j] = h;
      }
    }
  } else if (kind === 'city') {
    const b = spec.blocks ?? {};
    const nbx = b.nx ?? 5;
    const nby = b.ny ?? 5;
    const size = b.size ?? 800;
    const gapX = grid.Lx / nbx;
    const gapY = grid.Ly / nby;
    const H = b.height ?? 600;
    for (let j = 0; j < ny; j++) {
      for (let i = 0; i < nx; i++) {
        const x = grid.xc(i), y = grid.yc(j);
        const bi = Math.floor(x / gapX), bj = Math.floor(y / gapY);
        const lx = x - bi * gapX, ly = y - bj * gapY;
        if (lx < size && ly < size) height[i + nx * j] = H;
      }
    }
  }

  const mask = new Uint8Array(grid.nCell);
  const surfaceK = new Int32Array(nXY).fill(nz - 1);
  let solidCount = 0;
  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      const h = height[i + nx * j];
      let kTop = 0;
      for (let k = 0; k < nz; k++) {
        if (grid.zc(k) < h) {
          mask[i + nx * (j + ny * k)] = 1;
          solidCount++;
          kTop = k + 1;
        }
      }
      surfaceK[i + nx * j] = kTop >= nz ? nz - 1 : kTop;
    }
  }

  return { kind, height, mask, surfaceK, solidCount, nCell: grid.nCell };
}

/** 该列最低的流体层（地形之上的第一层） */
export function isSolid(terrain, n) {
  return terrain ? terrain.mask[n] === 1 : false;
}

/**
 * 把固壁面的法向速度清零（直接强迫 IBM 的一半；另一半在压力求解器里取 Neumann）。
 * 面的归属：u[n] 是格点 n 的西面，v[n] 是南面，w[n] 是底面。
 */
export function maskVelocities(S, grid, terrain) {
  if (!terrain) return;
  const { nx, ny, nz } = grid;
  const nXY = nx * ny;
  const m = terrain.mask;
  for (let k = 0; k < nz; k++) {
    const off = nXY * k;
    for (let j = 0; j < ny; j++) {
      const jm = j === 0 ? ny - 1 : j - 1;
      for (let i = 0; i < nx; i++) {
        const n = off + nx * j + i;
        if (m[n]) { S.u[n] = 0; S.v[n] = 0; S.w[n] = 0; S.w[n + nXY] = 0; continue; }
        const im = i === 0 ? nx - 1 : i - 1;
        if (m[off + nx * j + im]) S.u[n] = 0;
        if (m[off + nx * jm + i]) S.v[n] = 0;
        if (k > 0 && m[n - nXY]) S.w[n] = 0;
        if (k + 1 < nz && m[n + nXY]) S.w[n + nXY] = 0;
      }
    }
  }
  for (let n = 0; n < nXY; n++) { S.w[n] = 0; S.w[nXY * nz + n] = 0; }
}

/** 固壁内部填回基态，避免渲染时看到"地下的云" */
export function resetSolidState(S, base, grid, terrain) {
  if (!terrain) return;
  const { nx, ny, nz } = grid;
  const nXY = nx * ny;
  const m = terrain.mask;
  for (let k = 0; k < nz; k++) {
    const off = nXY * k;
    for (let n = 0; n < nXY; n++) {
      const idx = off + n;
      if (!m[idx]) continue;
      S.th[idx] = base.theta0[k];
      S.qv[idx] = base.qv0[k];
      S.qc[idx] = 0;
      S.qr[idx] = 0;
      S.qi[idx] = 0;
      S.qs[idx] = 0;
      S.qg[idx] = 0;
      S.u[idx] = 0; S.v[idx] = 0;
      S.w[idx] = 0; S.w[idx + nXY] = 0;
    }
  }
}
