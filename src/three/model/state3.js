/**
 * 三维状态容器与守恒式非负修补。不含物理，不含时间。
 *
 * 湿物种：qv 水汽、qc 云水、qr 雨、qi 云冰、qs 雪、qg 霰。
 */

export const FIELD_NAMES3 = ['u', 'v', 'w', 'th', 'qv', 'qc', 'qr', 'qi', 'qs', 'qg'];

export const MOIST_NAMES3 = ['qv', 'qc', 'qr', 'qi', 'qs', 'qg'];

export function allocState3(grid) {
  return {
    u: new Float64Array(grid.nCell),
    v: new Float64Array(grid.nCell),
    w: new Float64Array(grid.nW),
    th: new Float64Array(grid.nCell),
    qv: new Float64Array(grid.nCell),
    qc: new Float64Array(grid.nCell),
    qr: new Float64Array(grid.nCell),
    qi: new Float64Array(grid.nCell),
    qs: new Float64Array(grid.nCell),
    qg: new Float64Array(grid.nCell),
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
  for (const key of MOIST_NAMES3) {
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
  const names = MOIST_NAMES3;
  const n = S.qv.length;
  let repaired = 0;
  let lost = 0;

  for (let m = 0; m < n; m++) {
    const vals = names.map((k) => S[k][m]);
    const defIdx = [];
    let deficit = 0;
    for (let s = 0; s < vals.length; s++) {
      if (vals[s] < 0) { deficit -= vals[s]; vals[s] = 0; defIdx.push(s); }
    }
    if (deficit <= 0) continue;

    const d0 = deficit;
    const maxIdx = vals.reduce((best, v, idx) => (v > vals[best] ? idx : best), 0);
    const take = Math.min(deficit, vals[maxIdx]);
    vals[maxIdx] -= take;
    deficit -= take;

    for (let s = 0; s < vals.length; s++) S[names[s]][m] = vals[s];

    repaired += d0 - deficit;
    if (deficit > 1e-30) { lost += deficit; repaired += 1; }
    void defIdx;
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
