/**
 * Grid3 —— 三维网格对象。
 * x、y 周期，z 上下刚性（w=0）。
 * 索引布局：n = i + nx·(j + ny·k)，x 最快，便于逐层做水平二维 FFT。
 */

export class Grid3 {
  /**
   * @param {{nx:number, ny:number, nz:number, Lx:number, Ly:number, H:number}} o
   */
  constructor(o) {
    const { nx, ny, nz, Lx, Ly, H } = o;
    const pow2 = (n) => Number.isInteger(n) && n >= 8 && (n & (n - 1)) === 0;
    if (!pow2(nx)) throw new Error('nx 必须为 ≥8 的 2 的幂：压力求解在 x 方向使用 FFT');
    if (!pow2(ny)) throw new Error('ny 必须为 ≥8 的 2 的幂：压力求解在 y 方向使用 FFT');
    if (!Number.isInteger(nz) || nz < 4) throw new Error('nz 非法');
    if (!(Lx > 0) || !(Ly > 0) || !(H > 0)) throw new Error('域尺寸非法');

    this.nx = nx;
    this.ny = ny;
    this.nz = nz;
    this.Lx = Lx;
    this.Ly = Ly;
    this.H = H;
    this.dx = Lx / nx;
    this.dy = Ly / ny;
    this.dz = H / nz;

    this.nCell = nx * ny * nz;          // 标量 / u / v
    this.nW = nx * ny * (nz + 1);       // w 多一层
    this.nXY = nx * ny;                 // 一个水平层
    this.volume = Lx * Ly * H;
    this.area = Lx * Ly;
  }

  /** 标量 / u / v 的偏移（k ∈ [0, nz-1]） */
  sIdx(i, j, k) {
    return i + this.nx * (j + this.ny * k);
  }

  /** w 的偏移（k ∈ [0, nz]） */
  wIdx(i, j, k) {
    return i + this.nx * (j + this.ny * k);
  }

  wrapX(i) {
    const nx = this.nx;
    return ((i % nx) + nx) % nx;
  }

  wrapY(j) {
    const ny = this.ny;
    return ((j % ny) + ny) % ny;
  }

  xc(i) { return (i + 0.5) * this.dx; }
  yc(j) { return (j + 0.5) * this.dy; }
  zc(k) { return (k + 0.5) * this.dz; }
}
