/**
 * 守恒重映射（分辨率切换时使用）。
 *
 * 粗化：精确体积平均，严格守恒。
 * 细化：分段常数重建，守恒但有数值扩散。
 * 因此只在场景边界允许改 Δx，运行期不得改。
 */

/**
 * 一维均匀网格重映射。
 * @param {Float64Array} src 输入（nSrc 个格值）
 * @param {number} nSrc
 * @param {number} nDst
 * @param {Float64Array} out 输出（nDst 个格值）
 */
export function remapUniform(src, nSrc, nDst, out) {
  if (nSrc === nDst) {
    out.set(src.subarray(0, nDst));
    return;
  }
  const ratio = nSrc / nDst;
  for (let j = 0; j < nDst; j++) {
    const lo = j * ratio;
    const hi = (j + 1) * ratio;
    let sum = 0;
    const i0 = Math.floor(lo);
    const i1 = Math.min(nSrc, Math.ceil(hi));
    for (let i = i0; i < i1; i++) {
      const left = i > lo ? i : lo;
      const right = i + 1 < hi ? i + 1 : hi;
      const ov = right - left;
      if (ov > 0) sum += ov * src[i];
    }
    out[j] = sum / ratio;
  }
}

/**
 * 二维重映射（x 先、z 后）。
 * @param {Float64Array} src  按 sIdx(i,k)=i+nx0*k 排列
 */
export function remap2D(src, nx0, nz0, nx1, nz1, out) {
  const tmp = new Float64Array(nx1 * nz0);
  const row = new Float64Array(nx0);
  const rowOut = new Float64Array(nx1);
  const col = new Float64Array(nz0);
  const colOut = new Float64Array(nz1);

  for (let k = 0; k < nz0; k++) {
    const off = nx0 * k;
    for (let i = 0; i < nx0; i++) row[i] = src[off + i];
    remapUniform(row, nx0, nx1, rowOut);
    for (let i = 0; i < nx1; i++) tmp[i + nx1 * k] = rowOut[i];
  }

  for (let i = 0; i < nx1; i++) {
    for (let k = 0; k < nz0; k++) col[k] = tmp[i + nx1 * k];
    remapUniform(col, nz0, nz1, colOut);
    for (let k = 0; k < nz1; k++) out[i + nx1 * k] = colOut[k];
  }
}
