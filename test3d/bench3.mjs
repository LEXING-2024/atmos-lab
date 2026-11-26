/** 稳定基准：固定步数计时（排除探测循环的抖动） */
import { SCALE_PRESETS3, buildGrid3, buildModel3 } from '../src/three/model/plan3.js';

const cases = [
  ['conv3d', 0, 40], ['conv3d', 1, 12], ['meso3d', 0, 40], ['les3d', 0, 40],
];

for (const [id, lv, steps] of cases) {
  const p = SCALE_PRESETS3.find((x) => x.id === id);
  const grid = buildGrid3(p, lv);
  const m = buildModel3(grid, p, {});
  m.initAtmosphere({ amplitude: 2, radius: Math.max(300, grid.Lx / 20), noise: 0.02 });
  const dt = m.suggestDt();
  m.step(dt); m.step(dt);                  // 热身
  const t0 = performance.now();
  for (let s = 0; s < steps; s++) m.step(dt);
  const ms = (performance.now() - t0) / steps;
  console.log(`${id} ${p.levels[lv].label}: ${ms.toFixed(1)} ms/步  (${grid.nCell} 格点, dt=${dt.toFixed(2)})`);
}
