/**
 * 网格精度基准：同一域尺寸下比较不同 nx×ny×nz 的每步耗时与“每模拟分钟耗时”。
 * 用法：node test3d/gridbench.mjs
 */
import { Grid3 } from '../src/three/core/grid3.js';
import { buildModel3 } from '../src/three/model/plan3.js';

const DOMAIN = { Lx: 12000, Ly: 12000, H: 9000 };
const CONFIGS = [
  { nx: 32, ny: 32, nz: 24 },
  { nx: 64, ny: 64, nz: 32 },
  { nx: 128, ny: 128, nz: 48 },
  { nx: 128, ny: 128, nz: 64 },
];

const fmt = (v, d = 2) => Number(v).toFixed(d);

for (const cfg of CONFIGS) {
  const grid = new Grid3({ ...DOMAIN, ...cfg });
  const preset = {
    Lx: DOMAIN.Lx, Ly: DOMAIN.Ly, H: DOMAIN.H,
    levels: [cfg], micro: true, rh: 0.5,
  };
  const mem0 = process.memoryUsage().rss;
  const m = buildModel3(grid, preset, { base: { RHsfc: 0.5 } });
  m.initAtmosphere({ amplitude: 0.5, radius: 1500, zCenter: 2000, noise: 0.01 });
  const dt = m.suggestDt();

  m.step(dt); // 预热
  const t0 = performance.now();
  const N = 2;
  for (let i = 0; i < N; i++) m.step(dt);
  const ms = (performance.now() - t0) / N;
  const mem = (process.memoryUsage().rss - mem0) / 1e6;

  let nan = false, maxW = 0;
  for (const a of [m.S.u, m.S.v, m.S.w, m.S.th]) {
    for (let i = 0; i < a.length; i++) {
      if (!Number.isFinite(a[i])) { nan = true; break; }
      if (a === m.S.w && Math.abs(a[i]) > maxW) maxW = Math.abs(a[i]);
    }
    if (nan) break;
  }

  const msPerMin = ms * (60 / dt);
  console.log(
    `${String(cfg.nx).padStart(3)}×${String(cfg.ny).padStart(3)}×${String(cfg.nz).padStart(2)}  ` +
    `格点 ${String(grid.nCell).padStart(8)}  dt ${fmt(dt, 2).padStart(6)} s  ` +
    `每步 ${fmt(ms, 1).padStart(7)} ms  每模拟分钟 ${fmt(msPerMin / 1000, 2).padStart(6)} s  ` +
    `内存 +${fmt(mem, 0).padStart(4)} MB  maxW ${fmt(maxW, 1)}  ` +
    `${nan ? 'NaN!' : 'ok'}`
  );
}
