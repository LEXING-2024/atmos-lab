/**
 * 把模拟结果真的画成 PNG，用来肉眼检查渲染层与物理形态。
 * 复用 src/viz/render.js 的 colorize()（与浏览器里画的是同一套配色）。
 *
 * 用法：node test/render_shot.mjs [档位=cloud] [分钟=20]
 */

import { deflateSync } from 'node:zlib';
import { writeFileSync, mkdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';

import { buildGrid, presetById } from '../src/model/plan.js';
import { AnelasticModel } from '../src/model/anelastic.js';
import { colorize } from '../src/viz/render.js';

/* ---------- 极简 PNG 编码器 ---------- */
const CRC_TABLE = (() => {
  const t = new Int32Array(256);
  for (let n = 0; n < 256; n++) {
    let c = n;
    for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
    t[n] = c;
  }
  return t;
})();

function crc32(buf) {
  let c = -1;
  for (let i = 0; i < buf.length; i++) c = CRC_TABLE[(c ^ buf[i]) & 0xff] ^ (c >>> 8);
  return (c ^ -1) >>> 0;
}

function chunk(type, data) {
  const len = Buffer.alloc(4);
  len.writeUInt32BE(data.length);
  const body = Buffer.concat([Buffer.from(type, 'ascii'), data]);
  const crc = Buffer.alloc(4);
  crc.writeUInt32BE(crc32(body));
  return Buffer.concat([len, body, crc]);
}

function encodePNG(w, h, rgba) {
  const ihdr = Buffer.alloc(13);
  ihdr.writeUInt32BE(w, 0);
  ihdr.writeUInt32BE(h, 4);
  ihdr[8] = 8; ihdr[9] = 6; ihdr[10] = 0; ihdr[11] = 0; ihdr[12] = 0;

  const raw = Buffer.alloc((w * 4 + 1) * h);
  for (let y = 0; y < h; y++) {
    raw[y * (w * 4 + 1)] = 0;
    Buffer.from(rgba.buffer, y * w * 4, w * 4).copy(raw, y * (w * 4 + 1) + 1);
  }

  return Buffer.concat([
    Buffer.from([0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]),
    chunk('IHDR', ihdr),
    chunk('IDAT', deflateSync(raw, { level: 9 })),
    chunk('IEND', Buffer.alloc(0)),
  ]);
}

/** 最近邻放大 */
function upscale(src, nx, nz, fx) {
  const W = nx * fx, H = nz * fx;
  const out = new Uint8ClampedArray(W * H * 4);
  for (let y = 0; y < H; y++) {
    const sy = Math.floor(y / fx);
    for (let x = 0; x < W; x++) {
      const sx = Math.floor(x / fx);
      const s = (sy * nx + sx) * 4, d = (y * W + x) * 4;
      out[d] = src[s]; out[d + 1] = src[s + 1]; out[d + 2] = src[s + 2]; out[d + 3] = 255;
    }
  }
  return { buf: out, W, H };
}

/* ---------- 跑模型 ---------- */
const SCENARIOS = {
  // 云分辨：经典深对流热泡
  cloud: { preset: 'cloud', sens: 0, width: 0, init: [3, 4000, 2000], minutes: 35, field: 'cloud' },
  // 宏观：城市热岛 / 海陆风类的中尺度热力环流（地表感热条带）
  macro: { preset: 'macro', sens: 250, width: 60000, init: [1.0, 30000, 1500], minutes: 180, field: 'w' },
  // 中尺度：地表加热 + 弱热泡
  meso: { preset: 'meso', sens: 200, width: 30000, init: [2, 15000, 2000], minutes: 90, field: 'w' },
  // 微尺度：单个热泡/尘卷尺度
  les: { preset: 'les', sens: 0, width: 0, init: [1.5, 400, 250], minutes: 12, field: 'cloud' },
};

const name = process.argv[2] || 'cloud';
const sc = SCENARIOS[name];
if (!sc) { console.error(`未知场景 ${name}，可选：${Object.keys(SCENARIOS).join(' / ')}`); process.exit(1); }

const minutes = Number(process.argv[3] || sc.minutes);
const preset = presetById(sc.preset);
const grid = buildGrid(preset, 1);   // 中间档

const model = new AnelasticModel(grid, {
  base: { TSfc: 299.5, RHsfc: 0.78 },   // lapse/zTrop 用默认（K/m、m）
  micro: { enable: true },
  turb: {
    surface: {
      sensible: sc.sens,
      latent: 0,
      xCenter: sc.width ? grid.Lx / 2 : null,
      xWidth: sc.width || grid.Lx,
    },
  },
  filter: 0.02,
});

const [amp, radius, zc] = sc.init;
model.initAtmosphere({ amplitude: amp, radius, zCenter: zc, xCenter: grid.Lx / 2, noise: 0.02 });

const target = minutes * 60;
let t = 0;
while (t < target) {
  const r = model.advance(Math.min(30, target - t), 60000);
  if (r.simAdvanced <= 0) break;
  t += r.simAdvanced;
}

const d = model.diagnostics();
console.log(`${name}（${sc.preset} 档 ${grid.nx}×${grid.nz}，dx=${grid.dx.toFixed(0)}m dz=${grid.dz.toFixed(0)}m），积分 ${(t / 60).toFixed(1)} min`);
console.log(`maxW ${d.maxW.toFixed(1)} m/s，maxU ${d.maxU.toFixed(1)} m/s，云顶 ${(d.cloudTop / 1000).toFixed(1)} km`);
console.log(`雨强 ${d.precipRate.toFixed(2)} mm/h，累计 ${d.precipAccum.toFixed(2)} mm，CFL ${d.cfl.toFixed(3)}`);

/* ---------- 拼成四条（云 / w / u / RH / 雨） ---------- */
const fields = ['cloud', 'w', 'u', 'rh', 'rain'];
const FX = 4;
const viz = preset.viz;
const panels = fields.map((f) => {
  const c = colorize(model, f, viz);
  return upscale(c.data, c.nx, c.nz, FX);
});

const W = panels[0].W;
const GAP = 3;
const H = panels.reduce((s, p) => s + p.H + GAP, 0) - GAP;
const out = new Uint8ClampedArray(W * H * 4).fill(0);

let y0 = 0;
for (const p of panels) {
  for (let y = 0; y < p.H; y++) {
    const dst = ((y0 + y) * W) * 4;
    out.set(p.buf.subarray(y * W * 4, (y + 1) * W * 4), dst);
  }
  y0 += p.H;
  if (y0 < H) {
    for (let x = 0; x < W; x++) {
      const o = (y0 * W + x) * 4;
      out[o] = out[o + 1] = out[o + 2] = 255;
    }
    y0 += GAP;
  }
}

const outDir = fileURLToPath(new URL('./out/', import.meta.url));
mkdirSync(outDir, { recursive: true });
const file = `${outDir}shot_${name}.png`;
writeFileSync(file, encodePNG(W, H, out));
console.log(`已写出 ${file}（上→下：云与水 / 垂直速度 / 水平风 / 相对湿度 / 雨水）`);
