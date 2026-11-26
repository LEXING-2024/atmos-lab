/** 三维静力平衡：无扰动 60 步后不应出现任何自发运动 */
import { Grid3 } from '../src/three/core/grid3.js';
import { AnelasticModel3 } from '../src/three/model/anelastic3.js';
import { maxAbs3 } from '../src/three/model/state3.js';

const g = new Grid3({ nx: 32, ny: 32, nz: 24, Lx: 12000, Ly: 12000, H: 9000 });
const m = new AnelasticModel3(g, { micro: { enable: false } });
m.initAtmosphere({ amplitude: 0, noise: 0 });

let maxW = 0, divMax = 0;
const t0 = Date.now();
for (let s = 1; s <= 60; s++) {
  m.step(2);
  maxW = Math.max(maxW, maxAbs3(m.S.w));
  divMax = Math.max(divMax, m.divResidual);
}
const ms = Date.now() - t0;

console.log(`60 步 max|w| = ${maxW.toExponential(3)}   max 散度残差 = ${divMax.toExponential(3)}`);
console.log(`耗时 ${ms} ms，${(ms / 60).toFixed(1)} ms/步，格点 ${g.nCell}`);
console.log(maxW < 1e-10 ? '判定: 通过' : '判定: 失败');
