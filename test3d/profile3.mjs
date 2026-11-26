/** 单步耗时分解：平流/倾向、压力投影、微物理、湍流、滤波各占多少 */
import { Grid3 } from '../src/three/core/grid3.js';
import { buildModel3 } from '../src/three/model/plan3.js';

const DOMAIN = { Lx: 12000, Ly: 12000, H: 9000 };
const CONFIGS = [{ nx: 32, ny: 32, nz: 24 }, { nx: 64, ny: 64, nz: 32 }];

for (const cfg of CONFIGS) {
  const grid = new Grid3({ ...DOMAIN, ...cfg });
  const m = buildModel3(grid, {
    Lx: DOMAIN.Lx, Ly: DOMAIN.Ly, H: DOMAIN.H, levels: [cfg], micro: true, rh: 0.5,
  }, { base: { RHsfc: 0.5 } });
  m.initAtmosphere({ amplitude: 0.5, radius: 1500, zCenter: 2000, noise: 0.01 });
  const dt = m.suggestDt();

  const acc = {};
  const wrap = (obj, name, key) => {
    const orig = obj[name].bind(obj);
    acc[key] = 0;
    obj[name] = (...args) => {
      const t0 = performance.now();
      const r = orig(...args);
      acc[key] += performance.now() - t0;
      return r;
    };
  };
  wrap(m, 'computeTendencies', 'tend');
  wrap(m.solver, 'project', 'proj');
  wrap(m.micro, 'apply', 'micro');
  wrap(m.turb, 'apply', 'turb');
  wrap(m, 'applyFilter', 'filter');

  m.step(dt);
  for (const k in acc) acc[k] = 0;

  const N = 3;
  const t0 = performance.now();
  for (let i = 0; i < N; i++) m.step(dt);
  const total = (performance.now() - t0) / N;
  const parts = Object.entries(acc).map(([k, v]) => [k, v / N]);
  parts.sort((a, b) => b[1] - a[1]);
  console.log(`\n${cfg.nx}×${cfg.ny}×${cfg.nz}  每步 ${total.toFixed(1)} ms  dt ${dt.toFixed(2)} s`);
  for (const [k, v] of parts) {
    const names = { tend: '倾向(含平流/浮力/辐射)', proj: '压力投影', micro: '微物理', turb: '湍流', filter: 'Δ⁴滤波' };
    console.log(`  ${names[k].padEnd(16)} ${v.toFixed(1).padStart(6)} ms  ${(v / total * 100).toFixed(0).padStart(3)}%`);
  }
  console.log(`  其它              ${(total - parts.reduce((s, p) => s + p[1], 0)).toFixed(1).padStart(6)} ms`);
}
