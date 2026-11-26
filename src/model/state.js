/**
 * 状态容器与线性运算。不含物理，不含时间。
 */

export const FIELD_NAMES = ['u', 'w', 'th', 'qv', 'qc', 'qr'];

export function allocState(grid) {
  return {
    u: new Float64Array(grid.nCell),
    w: new Float64Array(grid.nW),
    th: new Float64Array(grid.nCell),
    qv: new Float64Array(grid.nCell),
    qc: new Float64Array(grid.nCell),
    qr: new Float64Array(grid.nCell),
  };
}

export function copyState(dst, src) {
  for (const f of FIELD_NAMES) dst[f].set(src[f]);
}

/** dst = wa*a + wb*b */
export function combineState(dst, a, wa, b, wb) {
  for (const f of FIELD_NAMES) {
    const D = dst[f], A = a[f], B = b[f];
    for (let n = 0; n < D.length; n++) D[n] = wa * A[n] + wb * B[n];
  }
}

/** 水物质非负；同时给出被裁掉的量（用于守恒诊断） */
export function clampMoisture(S) {
  let clipped = 0;
  for (const key of ['qv', 'qc', 'qr']) {
    const arr = S[key];
    for (let n = 0; n < arr.length; n++) {
      if (arr[n] < 0) { clipped -= arr[n]; arr[n] = 0; }
    }
  }
  return clipped;
}

/**
 * 守恒式非负修补。
 * 输运后个别格点可能出现微小负值；直接把负值置零会凭空造水，
 * 因此把缺口从同格的水汽（其次云水、雨水）里扣掉，保持总水量不变。
 * 返回 { repaired, lost }（kg/kg，域内合计，供诊断）。
 */
export function repairMoisture(S, nx, nz) {
  const { qv, qc, qr } = S;
  let repaired = 0;
  let lost = 0;

  const cover = (n, deficit) => {
    // 先在本格内，从水量最大的分量扣
    let a = qv[n], c = qc[n], r = qr[n];
    if (deficit <= 0) return 0;
    const d0 = deficit;
    if (a >= c && a >= r) {
      const take = deficit < a ? deficit : a;
      a -= take; deficit -= take;
    } else if (c >= r) {
      const take = deficit < c ? deficit : c;
      c -= take; deficit -= take;
    } else {
      const take = deficit < r ? deficit : r;
      r -= take; deficit -= take;
    }
    qv[n] = a; qc[n] = c; qr[n] = r;
    repaired += d0 - deficit;
    return deficit;
  };

  for (let i = 0; i < nx; i++) {
    for (let k = 0; k < nz; k++) {
      const n = i + nx * k;
      let deficit = 0;
      if (qv[n] < 0) { deficit -= qv[n]; qv[n] = 0; }
      if (qc[n] < 0) { deficit -= qc[n]; qc[n] = 0; }
      if (qr[n] < 0) { deficit -= qr[n]; qr[n] = 0; }
      if (deficit <= 0) continue;

      deficit = cover(n, deficit);
      if (deficit > 1e-30) { lost += deficit; repaired += 1; }
    }
  }

  return { repaired, lost };
}

export function maxAbs(arr) {
  let m = 0;
  for (let n = 0; n < arr.length; n++) {
    const a = arr[n] < 0 ? -arr[n] : arr[n];
    if (a > m) m = a;
  }
  return m;
}
