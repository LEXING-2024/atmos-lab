/** v1.3 地形自检：固壁无穿透、静力静止、绕流流场、耗时 */
import { Grid3 } from '../src/three/core/grid3.js';
import { AnelasticModel3 } from '../src/three/model/anelastic3.js';
import { makeBase } from '../src/model/base.js';

const g = new Grid3({ nx: 32, ny: 32, nz: 24, Lx: 12000, Ly: 12000, H: 9000 });
const base = makeBase(g, { TSfc: 299.5, RHsfc: 0.78 });

function checkSolid(m) {
  let maxIn = 0;
  for (let n = 0; n < g.nCell; n++) {
    if (!m.terrain.mask[n]) continue;
    const a = Math.abs(m.S.u[n]) + Math.abs(m.S.v[n]) + Math.abs(m.S.w[n]) + Math.abs(m.S.w[n + g.nXY]);
    if (a > maxIn) maxIn = a;
  }
  return maxIn;
}

// 1) 静力：无扰动、无风
{
  const m = new AnelasticModel3(g, {
    base: { TSfc: 299.5, RHsfc: 0.78 },
    terrain: { kind: 'ridge', height: 1500, halfWidth: 2500 },
    micro: { enable: false },
  });
  m.initAtmosphere({ amplitude: 0, radius: 2000, zCenter: 1200, noise: 0 });
  const t0 = Date.now();
  m.advance(300, 0);
  const d = m.diagnostics();
  console.log(`静力：${(d.time / 60).toFixed(0)} min  max|w| ${d.maxW.toExponential(2)}  max|u| ${d.maxU.toExponential(2)}  固壁内速度 ${checkSolid(m).toExponential(2)}`);
  console.log(`      散度残差 ${d.divResidual.toExponential(2)}  迭代 ${m.ibmIterations}  ${((Date.now() - t0) / Math.round(d.time / m.lastDt)).toFixed(0)} ms/步`);
}

// 2) 8 m/s 环境风越过山脊
{
  const m = new AnelasticModel3(g, {
    base: { TSfc: 299.5, RHsfc: 0.78 },
    terrain: { kind: 'ridge', height: 1500, halfWidth: 2500 },
    micro: { enable: false },
  });
  m.initAtmosphere({ amplitude: 0, radius: 1000, zCenter: 1000, noise: 0, wind: { u0: 8, shearU: 0 } });
  const t0 = Date.now();
  m.advance(900, 0);
  const d = m.diagnostics();
  console.log(`越山：${(d.time / 60).toFixed(0)} min  max|u| ${d.maxU.toFixed(2)}  max|w| ${d.maxW.toFixed(2)}  固壁内速度 ${checkSolid(m).toExponential(2)}`);
  console.log(`      散度残差 ${d.divResidual.toExponential(2)}  迭代 ${m.ibmIterations}  步数 ${Math.round(d.time / m.lastDt)}  耗时 ${Date.now() - t0} ms`);

  // 山脊上方 z=3 km 的 w 剖面（看迎风波/背风波）
  let row = '';
  const k = Math.round(3000 / g.dz);
  for (let i = 0; i < g.nx; i += 2) {
    const w = 0.5 * (m.S.w[i + g.nx * (16 + g.ny * k)] + m.S.w[i + g.nx * (16 + g.ny * (k + 1))]);
    row += `${w >= 0 ? '+' : ''}${w.toFixed(1)} `;
  }
  console.log(`      z≈3 km 的 w（每 2 格）：${row}`);
}
