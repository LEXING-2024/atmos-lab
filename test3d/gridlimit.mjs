/** 单机可算网格上限探测：构造模型 → 量内存 → 跑一步测时间 */
import { Grid3 } from '../src/three/core/grid3.js';
import { buildModel3 } from '../src/three/model/plan3.js';

const DOMAIN = { Lx: 24000, Ly: 24000, H: 12000 };
const sizes = [
  [128, 128, 48],
  [256, 256, 64],
  [256, 256, 128],
  [512, 512, 64],
];

const mb = (b) => (b / 1048576).toFixed(0);
const maxCells = Number(process.argv[2] || Infinity);
const only = process.argv[3] || '';
const onlyCfg = /^\d+x\d+x\d+$/.test(only) ? only.split('x').map(Number) : null;
const list = onlyCfg ? [onlyCfg] : sizes;

for (const [nx, ny, nz] of list) {
  const cells = nx * ny * nz;
  if (cells > maxCells) continue;
  if (only && !onlyCfg && only !== `${nx}x${ny}x${nz}`) continue;
  const cfg = { nx, ny, nz };
  console.log(`\n=== ${nx}×${ny}×${nz}  ${(cells / 1e6).toFixed(2)} M 格点 ===`);
  const before = process.memoryUsage();
  let m;
  try {
    const grid = new Grid3({ ...DOMAIN, ...cfg });
    m = buildModel3(grid, {
      Lx: DOMAIN.Lx, Ly: DOMAIN.Ly, H: DOMAIN.H, levels: [cfg], micro: true, rh: 0.5,
    }, { base: { RHsfc: 0.5 } });
  } catch (e) {
    console.log('构造失败：', e.message);
    continue;
  }
  const after = process.memoryUsage();
  console.log(`内存：external ${mb(after.external)} MB（另有 rss ${mb(after.rss)} MB）  增量 external ${mb(after.external - before.external)} MB`);
  console.log(`理论：状态 5×80 B/格 + 倾向 80 + 求解器 40 + 杂项 16 ≈ 0.52 KB/格 → ${(cells * 0.52 / 1024).toFixed(0)} MB`);

  try {
    m.initAtmosphere({ amplitude: 0.3, radius: 1500, zCenter: 2000, noise: 0.01 });
    const dt = m.suggestDt();
    const t0 = performance.now();
    m.step(dt);
    const ms = performance.now() - t0;
    console.log(`首步 ${(ms / 1000).toFixed(2)} s（含一次冷启动）  dt ${dt.toFixed(2)} s  每模拟分钟 ${(ms * 60 / dt / 1000).toFixed(1)} s`);
    if (ms > 120000) console.log('超过 2 min/步，停止继续；内存与时间已足够外推');
    else {
      const t1 = performance.now();
      m.step(dt);
      console.log(`第二步 ${((performance.now() - t1) / 1000).toFixed(2)} s`);
    }
  } catch (e) {
    console.log('步进失败：', e.message);
  }
  m = null;
  await new Promise((r) => setTimeout(r, 300));
}
