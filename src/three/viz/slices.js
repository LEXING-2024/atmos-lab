/**
 * 三维切片渲染。
 *
 * 只做"场 → 像素"的映射：任意一个坐标面上的二维切片，输出 w×h 的 RGBA。
 * 三个视图（x–z、y–z、x–y）共用这里的色标与场定义。
 */

import { qsat } from '../../physics/thermo.js';

export const FIELD_LIST3 = [
  { id: 'cloud', name: '云与水' },
  { id: 'ice', name: '冰相' },
  { id: 'w', name: '垂直速度' },
  { id: 'speed', name: '水平风速' },
  { id: 'temp', name: '位温扰动' },
  { id: 'rain', name: '雨水' },
  { id: 'rh', name: '相对湿度' },
];

export const DEFAULT_VIZ3 = { wMax: 25, uMax: 25, thMax: 6, qcMax: 2.2e-3, qrMax: 2.5e-3, iceMax: 1.5e-4 };

const clamp01 = (v) => (v < 0 ? 0 : v > 1 ? 1 : v);

/** 发散色标：蓝 → 白 → 红 */
function diverging(t, out) {
  const a = clamp01(t);
  if (a < 0.5) {
    const u = a * 2;
    out[0] = 30 + 225 * u; out[1] = 70 + 185 * u; out[2] = 200 + 55 * u;
  } else {
    const u = (a - 0.5) * 2;
    out[0] = 255; out[1] = 255 - 185 * u; out[2] = 255 - 200 * u;
  }
  return out;
}

/** 顺序色标：深蓝 → 青 → 黄 */
function sequential(t, out) {
  const a = clamp01(t);
  if (a < 0.5) {
    const u = a * 2;
    out[0] = 10 + 10 * u; out[1] = 30 + 170 * u; out[2] = 90 + 140 * u;
  } else {
    const u = (a - 0.5) * 2;
    out[0] = 20 + 235 * u; out[1] = 200 + 40 * u; out[2] = 230 - 180 * u;
  }
  return out;
}

/** 切片尺寸（像素） */
export function sliceSize(model, plane) {
  const g = model.grid;
  if (plane === 'xz') return { w: g.nx, h: g.nz };
  if (plane === 'yz') return { w: g.ny, h: g.nz };
  return { w: g.nx, h: g.ny };
}

/** 切片物理宽度与高度（米），供坐标轴使用 */
export function sliceExtent(model, plane) {
  const g = model.grid;
  if (plane === 'xz') return { lw: g.Lx, lh: g.H, labelW: 'x →  km', labelH: 'z ↑  km' };
  if (plane === 'yz') return { lw: g.Ly, lh: g.H, labelW: 'y →  km', labelH: 'z ↑  km' };
  return { lw: g.Lx, lh: g.Ly, labelW: 'x →  km', labelH: 'y ↑  km' };
}

/**
 * 生成一个坐标面的切片像素。
 * @param {number} index 该面的格点下标（xz 用 j，yz 用 i，xy 用 k）
 */
export function colorizeSlice(model, field = 'cloud', plane = 'xz', index = 0, viz = DEFAULT_VIZ3) {
  const g = model.grid, S = model.S, base = model.base;
  const terrain = model.terrain || null;
  const { w, h } = sliceSize(model, plane);
  const data = new Uint8ClampedArray(w * h * 4);
  const px = [0, 0, 0];
  const v = { ...DEFAULT_VIZ3, ...viz };
  void qsat;

  for (let hh = 0; hh < h; hh++) {
    for (let ww = 0; ww < w; ww++) {
      let i, j, k;
      if (plane === 'xz') { i = ww; k = hh; j = Math.min(g.ny - 1, Math.max(0, index)); }
      else if (plane === 'yz') { j = ww; k = hh; i = Math.min(g.nx - 1, Math.max(0, index)); }
      else { i = ww; j = hh; k = Math.min(g.nz - 1, Math.max(0, index)); }

      const n = i + g.nx * (j + g.ny * k);
      const o = ((h - 1 - hh) * w + ww) * 4;
      const zf = plane === 'xy' ? k / (g.nz - 1) : hh / (h - 1);

      if (terrain && terrain.mask[n]) {
        // 固壁：按高度给一点深浅，让山体和街区能看出来
        const shade = 0.55 + 0.45 * (plane === 'xy' ? 1 : zf);
        data[o] = 78 * shade;
        data[o + 1] = 66 * shade;
        data[o + 2] = 56 * shade;
        data[o + 3] = 255;
        continue;
      }

      if (field === 'cloud') {
        const skyR = 160 - 120 * zf, skyG = 200 - 110 * zf, skyB = 240 - 80 * zf;
        // 云（凝结物）用白，雨滴单独用偏灰蓝的半透明纱，避免雨把整根柱子刷白
        const aCloud = Math.pow(clamp01(S.qc[n] / v.qcMax), 0.65);
        const aIce = Math.pow(clamp01((S.qi[n] + S.qs[n] + S.qg[n]) / v.iceMax), 0.65);
        const aRain = Math.pow(clamp01(S.qr[n] / v.qrMax), 0.8) * 0.55;
        let r = skyR + (245 - skyR) * aCloud;
        let gg = skyG + (250 - skyG) * aCloud;
        let b = skyB + (255 - skyB) * aCloud;
        if (aIce > 0) {
          // 冰相偏冷白：云冰最白，雪次之，霰带一点灰
          const graupelShare = clamp01(S.qg[n] / Math.max(1e-12, S.qi[n] + S.qs[n] + S.qg[n]));
          const ir = 235 - 25 * graupelShare;
          const ig = 245 - 20 * graupelShare;
          const ib = 255 - 15 * graupelShare;
          r = r * (1 - aIce) + ir * aIce;
          gg = gg * (1 - aIce) + ig * aIce;
          b = b * (1 - aIce) + ib * aIce;
        }
        if (aRain > 0) {
          r = r * (1 - aRain) + 132 * aRain;
          gg = gg * (1 - aRain) + 146 * aRain;
          b = b * (1 - aRain) + 168 * aRain;
        }
        data[o] = r;
        data[o + 1] = gg;
        data[o + 2] = b;
      } else if (field === 'ice') {
        // 冰相单独成场：云冰白、雪亮青、霰灰紫，三种冰同一量纲叠色
        const qi = S.qi[n], qsv = S.qs[n], qg = S.qg[n];
        const tot = qi + qsv + qg;
        const a = Math.pow(clamp01(tot / v.iceMax), 0.6);
        const gf = tot > 1e-12 ? qg / tot : 0;
        const sf = tot > 1e-12 ? qsv / tot : 0;
        const r0 = 225 + 15 * gf;
        const g0 = 238 - 25 * gf + 5 * sf;
        const b0 = 255 - 40 * gf;
        data[o] = 12 + (r0 - 12) * a;
        data[o + 1] = 22 + (g0 - 22) * a;
        data[o + 2] = 45 + (b0 - 45) * a;
      } else if (field === 'rain') {
        const a = clamp01(S.qr[n] / v.qrMax);
        data[o] = 20 + 10 * a;
        data[o + 1] = 40 + 150 * a;
        data[o + 2] = 60 + 150 * a;
      } else if (field === 'temp') {
        diverging((S.th[n] - base.theta0[k]) / v.thMax * 0.5 + 0.5, px);
        data[o] = px[0]; data[o + 1] = px[1]; data[o + 2] = px[2];
      } else if (field === 'speed') {
        const spd = Math.sqrt(S.u[n] * S.u[n] + S.v[n] * S.v[n]);
        sequential(spd / v.uMax, px);
        data[o] = px[0]; data[o + 1] = px[1]; data[o + 2] = px[2];
      } else if (field === 'w') {
        const wv = 0.5 * (S.w[n] + S.w[n + g.nXY]);
        diverging(wv / v.wMax * 0.5 + 0.5, px);
        data[o] = px[0]; data[o + 1] = px[1]; data[o + 2] = px[2];
      } else {
        const T = S.th[n] * base.ex[k];
        const rh = clamp01(S.qv[n] / qsat(T, base.p0[k]));
        data[o] = 120 + 20 * rh;
        data[o + 1] = 90 + 130 * rh;
        data[o + 2] = 60 + 40 * rh;
      }
      data[o + 3] = 255;
    }
  }
  return { data, w, h };
}

/** 把一个切片画到画布上（含坐标轴） */
export function createSliceView(canvas) {
  const ctx = canvas.getContext('2d', { alpha: false });
  const off = document.createElement('canvas');
  const octx = off.getContext('2d', { willReadFrequently: true });
  let img = null, cw = 0, chh = 0;

  function ensure(w, h) {
    if (w === cw && h === chh && img) return;
    cw = w; chh = h;
    off.width = w; off.height = h;
    img = octx.createImageData(w, h);
  }

  function resize() {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const w = canvas.clientWidth || 400;
    const h = canvas.clientHeight || 300;
    const px = Math.round(w * dpr), py = Math.round(h * dpr);
    if (canvas.width !== px || canvas.height !== py) { canvas.width = px; canvas.height = py; }
    return { w: canvas.width, h: canvas.height, dpr };
  }

  function draw(model, opts = {}) {
    const plane = opts.plane || 'xz';
    const { data, w, h } = colorizeSlice(model, opts.field, plane, opts.index ?? 0, opts.viz);
    ensure(w, h);
    img.data.set(data);
    octx.putImageData(img, 0, 0);

    const { w: W, h: H, dpr } = resize();
    ctx.save();
    ctx.fillStyle = '#05070d';
    ctx.fillRect(0, 0, W, H);
    ctx.imageSmoothingEnabled = true;

    const padL = Math.round(38 * dpr), padB = Math.round(20 * dpr), padT = Math.round(16 * dpr);
    const plotW = Math.max(1, W - padL), plotH = Math.max(1, H - padB - padT);
    ctx.drawImage(off, padL, padT, plotW, plotH);

    const ext = sliceExtent(model, plane);
    const kmW = ext.lw / 1000, kmH = ext.lh / 1000;

    ctx.font = `${10 * dpr}px ui-monospace, Menlo, Consolas, monospace`;
    ctx.fillStyle = '#8fa3bd';
    ctx.strokeStyle = 'rgba(140,165,195,0.3)';
    ctx.lineWidth = 1;

    const vStep = niceStep(kmH, Math.max(2, Math.floor(plotH / (26 * dpr))));
    for (let z = 0; z <= kmH + 1e-9; z += vStep) {
      const y = padT + plotH - (z / kmH) * plotH;
      ctx.beginPath(); ctx.moveTo(padL - 4 * dpr, y); ctx.lineTo(padL, y); ctx.stroke();
      ctx.textAlign = 'right';
      ctx.fillText(`${z.toFixed(vStep < 1 ? 1 : 0)}`, padL - 6 * dpr, y + 3 * dpr);
    }
    ctx.textAlign = 'left';
    ctx.fillText('km', 2 * dpr, 10 * dpr);

    const hStep = niceStep(kmW, Math.max(2, Math.floor(plotW / (56 * dpr))));
    for (let x = 0; x <= kmW + 1e-9; x += hStep) {
      const px = padL + (x / kmW) * plotW;
      ctx.beginPath(); ctx.moveTo(px, padT + plotH); ctx.lineTo(px, padT + plotH + 4 * dpr); ctx.stroke();
      ctx.textAlign = x === 0 ? 'left' : (Math.abs(x - kmW) < 1e-9 ? 'right' : 'center');
      ctx.fillText(`${x.toFixed(hStep < 1 ? 1 : 0)}`, px, padT + plotH + 13 * dpr);
    }
    ctx.restore();
  }

  return { draw, resize };
}

function niceStep(span, maxTicks) {
  const raw = span / Math.max(1, maxTicks);
  const p = Math.pow(10, Math.floor(Math.log10(raw)));
  const m = raw / p;
  const f = m <= 1 ? 1 : m <= 2 ? 2 : m <= 5 ? 5 : 10;
  return f * p;
}
