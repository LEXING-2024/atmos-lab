import { Grid } from '../src/core/grid.js';
import { AnelasticModel } from '../src/model/anelastic.js';
import { maxAbs } from '../src/model/state.js';

const g = new Grid({ nx: 64, nz: 48, Lx: 32000, H: 12000 });
const m = new AnelasticModel(g, { micro: { enable: false }, turb: { surface: { sensible: 0, latent: 0 } } });
m.initAtmosphere({ amplitude: 0, noise: 0 });

console.log('初始 max|u|', maxAbs(m.S.u), 'max|w|', maxAbs(m.S.w));
console.log('初始 th 与基态最大差', (() => {
  let d = 0;
  for (let k = 0; k < g.nz; k++) for (let i = 0; i < g.nx; i++) {
    const v = Math.abs(m.S.th[i + g.nx * k] - m.base.theta0[k]);
    if (v > d) d = v;
  }
  return d;
})());

m.computeTendencies(m.S, m.T);
console.log('max|tu|', maxAbs(m.T.tu), 'max|tw|', maxAbs(m.T.tw), 'max|tth|', maxAbs(m.T.tth));

const dt = 8;
m.step(dt);
console.log('1步后 max|u|', maxAbs(m.S.u), 'max|w|', maxAbs(m.S.w), 'divResidual', m.divResidual);

for (let n = 0; n < 10; n++) m.step(dt);
console.log('11步后 max|u|', maxAbs(m.S.u), 'max|w|', maxAbs(m.S.w), 'divResidual', m.divResidual);

// 找出 w 最大的位置
let best = -1, bv = 0;
for (let n = 0; n < m.S.w.length; n++) if (Math.abs(m.S.w[n]) > bv) { bv = Math.abs(m.S.w[n]); best = n; }
console.log('max w 位置 i,k =', best % g.nx, Math.floor(best / g.nx), 'z=', Math.floor(best / g.nx) * g.dz);
