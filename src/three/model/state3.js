/**
 * 三维状态容器与守恒式非负修补。不含物理，不含时间。
 */

export const FIELD_NAMES3 = ['u', 'v', 'w', 'th', 'qv', 'qc', 'qr'];

export function allocState3(grid) {
  return {
    u: new Float64Array(grid.nCell),
    v: new Float64Array(grid.nCell),
    w: new Float64Array(grid.nW),
    th: new Float64Array(grid.nCell),
    qv: new Float64Array(grid.nCell),
    qc: new Float64Array(grid.nCell),
    qr: new Float64Array(grid.nCell),
  };
}

export function copyState3(dst, src) {
  for (const f of FIELD_NAMES3) dst[f].set(src[f]);
}

/** dst = wa*a + wb*b */
export function combineState3(dst, a, wa, b, wb) {
  for (const f of FIELD_NAMES3) {
    const D = dst[f], A = a[f], B = b[f];
    for (let n = 0; n < D.length; n++) D[n] = wa * A[n] + wb * B[n];
  }
}

export function clampMoisture3(S) {
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
 * 守恒式非负修补（同二维版）：负值抹平后从同格水量最大的分量扣等质量。
 * 返回 { repaired, lost }。
 */
export function repairMoisture3(S) {
  const { qv, qc, qr } = S;
  const n = qv.length;
  let repaired = 0;
  let lost = 0;

  for (let m = 0; m < n; m++) {
    let deficit = 0;
    if (qv[m] < 0) { deficit -= qv[m]; qv[m] = 0; }
    if (qc[m] < 0) { deficit -= qc[m]; qc[m] = 0; }
    if (qr[m] < 0) { deficit -= qr[m]; qr[m] = 0; }
    if (deficit <= 0) continue;

    const d0 = deficit;
    let a = qv[m], c = qc[m], r = qr[m];
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
    qv[m] = a; qc[m] = c; qr[m] = r;

    repaired += d0 - deficit;
    if (deficit > 1e-30) { lost += deficit; repaired += 1; }
  }

  return { repaired, lost };
}

export function maxAbs3(arr) {
  let m = 0;
  for (let n = 0; n < arr.length; n++) {
    const a = arr[n] < 0 ? -arr[n] : arr[n];
    if (a > m) m = a;
  }
  return m;
}
