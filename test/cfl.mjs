import { Grid } from '../src/core/grid.js';
import { AnelasticModel } from '../src/model/anelastic.js';

const stat = (a) => {
  let mx = 0, bad = 0;
  for (let n = 0; n < a.length; n++) {
    const v = a[n];
    if (!Number.isFinite(v)) { bad++; continue; }
    const abs = v < 0 ? -v : v;
    if (abs > mx) mx = abs;
  }
  return { mx, bad };
};

function run(label, extra) {
  const g = new Grid({ nx: 64, nz: 32, Lx: 16000, H: 8000 });
  const m = new AnelasticModel(g, { base: { TSfc: 299.5, RHsfc: 0.78 }, micro: { enable: false }, ...extra });
  m.initAtmosphere({ amplitude: 2.5, radius: 1500, noise: 0 });
  let last = null;
  for (let s = 0; s < 200; s++) {
    const dt = Math.min(m.suggestDt(), extra.fixedDt ?? m.suggestDt());
    m.step(dt);
    if (s % 20 === 0 || s === 199) {
      const w = stat(m.S.w), th = stat(m.S.th);
      last = `${label} step ${s}: maxW=${w.mx.toExponential(2)} maxTh=${th.mx.toExponential(2)} dt=${dt.toFixed(2)}`;
      if (w.bad || th.bad) { last += ` BAD=${w.bad + th.bad}`; break; }
    }
  }
  console.log(last);
}

run('dt=8   ', { fixedDt: 8, dtMax: 8 });
run('dt=2   ', { fixedDt: 2, dtMax: 2 });
run('dt=0.5 ', { fixedDt: 0.5, dtMax: 0.5 });
