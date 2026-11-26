import { Grid } from '../src/core/grid.js';
import { AnelasticModel } from '../src/model/anelastic.js';

/** 隔离测试：只保留平流，观察是否守恒、是否有虚假增长 */
function makeModel() {
  const g = new Grid({ nx: 64, nz: 32, Lx: 16000, H: 8000 });
  const m = new AnelasticModel(g, {
    base: { TSfc: 299.5, RHsfc: 0.05 },
    micro: { enable: false },
    turb: { cs: 0, kMin: 0, kMax: 0, surface: { sensible: 0, latent: 0 } },
    rad: { tauRad: 1e12, spongeMax: 0 },
  });
  m.initAtmosphere({ amplitude: 0.02, radius: 1500, noise: 0 });
  return m;
}

function maxOf(a) { let m = 0, bad = 0; for (let n = 0; n < a.length; n++) { const v = a[n]; if (!Number.isFinite(v)) { bad++; continue; } if (Math.abs(v) > m) m = Math.abs(v); } return { m, bad }; }

// A：水平均匀风 5 m/s
{
  const m = makeModel();
  m.S.u.fill(5);
  const t0 = m.diagnostics().totalWater;
  for (let s = 0; s < 400; s++) m.step(2);
  const th = maxOf(m.S.th), w = maxOf(m.S.w);
  console.log(`A 水平均匀风: maxTh=${th.m.toFixed(2)} maxW=${w.m.toExponential(2)} bad=${th.bad}`);
}

// B：垂直均匀风 2 m/s
{
  const m = makeModel();
  for (let i = 0; i < 32 + 1; i++) m.S.w.fill(0);
  for (let k = 1; k < 32; k++) for (let i = 0; i < 64; i++) m.S.w[i + 64 * k] = 2;
  for (let s = 0; s < 400; s++) m.step(2);
  const th = maxOf(m.S.th), w = maxOf(m.S.w);
  console.log(`B 垂直均匀风: maxTh=${th.m.toFixed(2)} maxW=${w.m.toExponential(2)} bad=${th.bad}`);
}

// C：无速度，只有浮力（检验重力波 CFL）
{
  const m = makeModel();
  for (let s = 0; s < 400; s++) m.step(2);
  const th = maxOf(m.S.th), w = maxOf(m.S.w);
  console.log(`C 纯浮力 dt=2: maxTh=${th.m.toFixed(2)} maxW=${w.m.toFixed(2)} bad=${th.bad}`);
}
