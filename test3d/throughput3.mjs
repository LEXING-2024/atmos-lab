/** 三维吞吐实测：各预设各分辨率档（先热身一次，消除 JIT 冷启动） */
import { SCALE_PRESETS3, measureThroughput3 } from '../src/three/model/plan3.js';

measureThroughput3(SCALE_PRESETS3[0], 0, 150);   // 热身

for (const p of SCALE_PRESETS3) {
  for (let lv = 0; lv < p.levels.length; lv++) {
    const r = measureThroughput3(p, lv, 150);
    const g = p.levels[lv];
    console.log(`${p.id.padEnd(7)} ${g.label.padEnd(11)} ${(g.nx * g.ny * g.nz).toLocaleString().padStart(9)} 格点   ${r.msPerStep.toFixed(1)} ms/步   ${Math.round(r.cellStepsPerMs)} 格点·步/ms`);
  }
}
