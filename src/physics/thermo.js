/**
 * 热力学基础函数。纯函数，无状态，无 DOM。
 */

export const CONST = {
  g: 9.81,
  Rd: 287.05,
  Rv: 461.5,
  cp: 1004.5,
  p00: 100000,
  Lv: 2.5e6,
  eps: 0.622,
};

/** 饱和水汽压（Bolton 公式），Pa */
export function esat(T) {
  return 611.2 * Math.exp(17.67 * (T - 273.15) / (T - 29.65));
}

/** 饱和比湿 */
export function qsat(T, p) {
  const es = esat(T);
  const d = p - es;
  return CONST.eps * es / (d > 1 ? d : 1);
}

/** Exner 函数 */
export function exner(p) {
  return Math.pow(p / CONST.p00, CONST.Rd / CONST.cp);
}

/**
 * 饱和调整的隐式增量。
 * 返回应当凝结（正）或蒸发（负）的比湿增量。
 */
export function satIncrement(T, p, qv, qc) {
  const qs = qsat(T, p);
  const Lv = CONST.Lv, cp = CONST.cp, Rv = CONST.Rv;
  const f = 1 + (Lv * Lv * qs) / (cp * Rv * T * T);

  if (qv > qs) {
    let dq = (qv - qs) / f;
    if (dq > qv) dq = qv;
    return dq;                       // 凝结
  }
  if (qc > 0 && qv < qs) {
    let dq = (qs - qv) / f;
    if (dq > qc) dq = qc;
    return -dq;                      // 蒸发
  }
  return 0;
}
