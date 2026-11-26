/**
 * Grid —— 运行时网格对象。
 * 全仓库唯一的几何来源：任何文件不得出现网格常量。
 */

export class Grid {
  /**
   * @param {{nx:number, nz:number, Lx:number, H:number}} o
   */
  constructor(o) {
    const { nx, nz, Lx, H } = o;
    if (!Number.isInteger(nx) || nx < 8) throw new Error('nx 非法');
    if ((nx & (nx - 1)) !== 0) {
      throw new Error('nx 必须为 2 的幂：压力求解在 x 方向使用 FFT');
    }
    if (!Number.isInteger(nz) || nz < 4) throw new Error('nz 非法');
    if (!(Lx > 0) || !(H > 0)) throw new Error('域尺寸非法');

    this.nx = nx;
    this.nz = nz;
    this.Lx = Lx;
    this.H = H;
    this.dx = Lx / nx;
    this.dz = H / nz;

    this.nCell = nx * nz;        // 标量与 u
    this.nW = nx * (nz + 1);     // w 多一层
    this.volume = Lx * H;        // 单位厚度
  }

  /** 标量 / u 的偏移（k ∈ [0, nz-1]） */
  sIdx(i, k) {
    return i + this.nx * k;
  }

  /** w 的偏移（k ∈ [0, nz]） */
  wIdx(i, k) {
    return i + this.nx * k;
  }

  /** 周期取模 */
  wrap(i) {
    const nx = this.nx;
    return ((i % nx) + nx) % nx;
  }

  /** 格点中心坐标 */
  xc(i) { return (i + 0.5) * this.dx; }
  zc(k) { return (k + 0.5) * this.dz; }
}
