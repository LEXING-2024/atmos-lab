import { Grid } from '../src/core/grid.js';
import { makeBase } from '../src/model/base.js';
import { CONST, qsat, exner } from '../src/physics/thermo.js';

const g = CONST.g, cp = CONST.cp, Rd = CONST.Rd, Rv = CONST.Rv, Lv = CONST.Lv;
const grid = new Grid({ nx: 64, nz: 64, Lx: 32000, H: 16000 });
const dzP = 25;

function parcelT(zStart, T0, p0) {
  let p = p0, T = T0;
  const out = [];
  for (let z = zStart; z < 16000; z += dzP) {
    out.push({ z, th: T / exner(p) });
    const qs = qsat(T, p);
    const num = 1 + Lv * qs / (Rd * T);
    const den = 1 + (Lv * Lv * qs) / (cp * Rv * T * T);
    T += -(g / cp) * num / den * dzP;
    p -= (p * g / (Rd * T)) * dzP;
  }
  return out;
}

function cape(base, zStart, label) {
  const dzGrid = base.zc[1] - base.zc[0];
  const kS = Math.round(zStart / dzGrid - 0.5);
  const par = parcelT(base.zc[kS], base.T0[kS], base.p0[kS]);
  const envAt = (z) => {
    const f = z / dzGrid - 0.5;
    const k0 = Math.max(0, Math.min(base.theta0.length - 2, Math.floor(f)));
    const t = Math.max(0, Math.min(1, f - k0));
    return base.theta0[k0] * (1 - t) + base.theta0[k0 + 1] * t;
  };
  let CAPE = 0, CIN = 0, LFC = null, EL = null, maxB = 0;
  for (const pt of par) {
    const thE = envAt(pt.z);
    const b = g * (pt.th - thE) / thE;
    if (b > 0) {
      if (LFC === null) LFC = pt.z;
      CAPE += b * dzP;
      if (b > maxB) maxB = b;
      EL = pt.z;
    } else if (LFC === null) {
      CIN += b * dzP;
    }
  }
  console.log(`${label.padEnd(26)} CAPE=${CAPE.toFixed(0).padStart(5)} J/kg  CIN=${CIN.toFixed(0).padStart(5)}  LFC=${LFC ? LFC.toFixed(0).padStart(5) : '   --'} m  EL=${EL ? EL.toFixed(0).padStart(5) : '   --'} m  maxB=${maxB.toFixed(2)}`);
}

const variants = [
  ['当前默认 RH0.78 lapse7', {}],
  ['RH0.90 lapse7', { RHsfc: 0.9 }],
  ['RH0.85 lapse6.5', { RHsfc: 0.85, lapse: 0.0065 }],
  ['RH0.88 lapse6.8 Hq4000', { RHsfc: 0.88, lapse: 0.0068, Hq: 4000 }],
  ['RH0.92 lapse6.5 Hq5000', { RHsfc: 0.92, lapse: 0.0065, Hq: 5000 }],
  ['RH0.90 lapse6.2 Hq6000', { RHsfc: 0.9, lapse: 0.0062, Hq: 6000 }],
];

for (const [label, o] of variants) cape(makeBase(grid, o), 500, label);
