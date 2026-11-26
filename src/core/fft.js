/**
 * 一维复数 FFT（radix-2，迭代版）。
 * 仅服务于压力求解，不含任何物理。
 */

export class FFT1D {
  constructor(n) {
    if ((n & (n - 1)) !== 0) throw new Error('FFT 长度必须为 2 的幂');
    this.n = n;

    const bits = Math.round(Math.log2(n));
    this.rev = new Uint32Array(n);
    for (let i = 0; i < n; i++) {
      let r = 0;
      for (let b = 0; b < bits; b++) r |= ((i >> b) & 1) << (bits - 1 - b);
      this.rev[i] = r;
    }

    this.cosTab = new Float64Array(n >> 1);
    this.sinTab = new Float64Array(n >> 1);
    for (let i = 0; i < (n >> 1); i++) {
      const a = -2 * Math.PI * i / n;
      this.cosTab[i] = Math.cos(a);
      this.sinTab[i] = Math.sin(a);
    }
  }

  /**
   * 原地变换。
   * @param {Float64Array} re
   * @param {Float64Array} im
   * @param {boolean} inverse
   */
  transform(re, im, inverse = false) {
    const n = this.n, rev = this.rev, ct = this.cosTab, st = this.sinTab;

    for (let i = 0; i < n; i++) {
      const j = rev[i];
      if (j > i) {
        let t = re[i]; re[i] = re[j]; re[j] = t;
        t = im[i]; im[i] = im[j]; im[j] = t;
      }
    }

    const sgn = inverse ? -1 : 1;
    for (let size = 2; size <= n; size <<= 1) {
      const half = size >> 1, step = n / size;
      for (let i = 0; i < n; i += size) {
        for (let j = 0; j < half; j++) {
          const k = j * step;
          const wr = ct[k], wi = st[k] * sgn;
          const a = i + j, b = a + half;
          const xr = re[b], xi = im[b];
          const tr = wr * xr - wi * xi;
          const ti = wr * xi + wi * xr;
          re[b] = re[a] - tr; im[b] = im[a] - ti;
          re[a] += tr; im[a] += ti;
        }
      }
    }

    if (inverse) {
      for (let i = 0; i < n; i++) { re[i] /= n; im[i] /= n; }
    }
  }
}
