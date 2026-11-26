/**
 * 渲染层。
 *
 * 只做"场 → 像素"的映射，不含任何数值逻辑。
 * 网格多大就渲染多大（离屏 nx×nz 的 ImageData），再上采样到画布，
 * 所以换分辨率时无需改这里。
 */

import { qsat } from '../physics/thermo.js';

export const FIELD_LIST = [
  { id: 'cloud', name: '云与水' },
  { id: 'w', name: '垂直速度' },
  { id: 'u', name: '水平风' },
  { id: 'temp', name: '位温扰动' },
  { id: 'rain', name: '雨水' },
  { id: 'rh', name: '相对湿度' },
];

const clamp01 = (v) => (v < 0 ? 0 : v > 1 ? 1 : v);

/** 默认色标（各尺度档会覆盖） */
export const DEFAULT_VIZ = { wMax: 15, uMax: 20, thMax: 6, qcMax: 2.2e-3, qrMax: 2.5e-3 };

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

/**
 * 把一个场映射成 nx×nz 的 RGBA 像素（第 0 行在顶部 = 模式顶）。
 * 与 canvas 无关，所以测试里可以直接调它导出 PNG 检查画面。
 */
export function colorize(model, field = 'cloud', viz = DEFAULT_VIZ) {
  const g = model.grid, S = model.S, base = model.base;
  const data = new Uint8ClampedArray(g.nx * g.nz * 4);
  const px = [0, 0, 0];

  const v = { ...DEFAULT_VIZ, ...viz };
  const dTemp = v.thMax;
  const dW = v.wMax;
  const dU = v.uMax;

  for (let k = 0; k < g.nz; k++) {
    const row = (g.nz - 1 - k) * g.nx;        // 第 0 行画在顶部
    const off = g.nx * k;
    const zf = k / (g.nz - 1);
    for (let i = 0; i < g.nx; i++) {
      const n = off + i;
      const o = (row + i) * 4;
      const qc = S.qc[n] + S.qr[n];

      if (field === 'cloud') {
        // 天空底色：近地面浅蓝，高空深蓝
        const skyR = 160 - 120 * zf, skyG = 200 - 110 * zf, skyB = 240 - 80 * zf;
        const a = clamp01(qc / v.qcMax);
        data[o] = skyR + (245 - skyR) * a;
        data[o + 1] = skyG + (250 - skyG) * a;
        data[o + 2] = skyB + (255 - skyB) * a;
      } else if (field === 'rain') {
        const a = clamp01(S.qr[n] / v.qrMax);
        data[o] = 20 + 10 * a;
        data[o + 1] = 40 + 150 * a;
        data[o + 2] = 60 + 150 * a;
      } else if (field === 'temp') {
        diverging((S.th[n] - base.theta0[k]) / dTemp * 0.5 + 0.5, px);
        data[o] = px[0]; data[o + 1] = px[1]; data[o + 2] = px[2];
      } else if (field === 'u') {
        const u = 0.5 * (S.u[n] + (i + 1 === g.nx ? S.u[off] : S.u[n + 1]));
        diverging(u / dU * 0.5 + 0.5, px);
        data[o] = px[0]; data[o + 1] = px[1]; data[o + 2] = px[2];
      } else if (field === 'w') {
        const w = 0.5 * (S.w[n] + S.w[n + g.nx]);
        diverging(w / dW * 0.5 + 0.5, px);
        data[o] = px[0]; data[o + 1] = px[1]; data[o + 2] = px[2];
      } else {
        const T = S.th[n] * base.ex[k];
        const rh = clamp01(S.qv[n] / qsat(T, base.p0[k]));
        // 干 → 棕，湿 → 绿
        data[o] = 120 + 20 * rh;
        data[o + 1] = 90 + 130 * rh;
        data[o + 2] = 60 + 40 * rh;
      }
      data[o + 3] = 255;
    }
  }
  return { data, nx: g.nx, nz: g.nz };
}

export function createRenderer(canvas) {
  const ctx = canvas.getContext('2d', { alpha: false });
  const off = document.createElement('canvas');
  const octx = off.getContext('2d', { willReadFrequently: true });

  let img = null;
  let nx = 0, nz = 0;
  const px = [0, 0, 0];

  function ensure(nxIn, nzIn) {
    if (nxIn === nx && nzIn === nz && img) return;
    nx = nxIn; nz = nzIn;
    off.width = nx; off.height = nz;
    img = octx.createImageData(nx, nz);
  }

  /** 画布尺寸随容器变化（含设备像素比） */
  function resize() {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const w = canvas.clientWidth || 800;
    const h = canvas.clientHeight || 480;
    const cw = Math.round(w * dpr), ch = Math.round(h * dpr);
    if (canvas.width !== cw || canvas.height !== ch) {
      canvas.width = cw; canvas.height = ch;
    }
    return { w: canvas.width, h: canvas.height, dpr };
  }

  /**
   * @param {import('../model/anelastic.js').AnelasticModel} model
   * @param {{field:string, scale:number, showAxes:boolean}} opts
   */
  function draw(model, opts = {}) {
    const g = model.grid;

    const field = opts.field || 'cloud';
    const viz = opts.viz;

    const { data, nx: cnx, nz: cnz } = colorize(model, field, viz);
    ensure(cnx, cnz);
    img.data.set(data);
    octx.putImageData(img, 0, 0);

    const { w: cw, h: ch } = resize();
    ctx.save();
    ctx.imageSmoothingEnabled = true;
    ctx.fillStyle = '#05070d';
    ctx.fillRect(0, 0, cw, ch);

    const padL = opts.showAxes === false ? 0 : Math.round(64 * (cw / 800));
    const padB = opts.showAxes === false ? 0 : Math.round(34 * (cw / 800));
    const plotW = cw - padL;
    const plotH = ch - padB;

    ctx.drawImage(off, padL, 0, plotW, plotH);

    if (opts.showAxes !== false) drawAxes(ctx, g, cw, ch, padL, padB, plotW, plotH);
    ctx.restore();
  }

  function drawAxes(ctx, g, cw, ch, padL, padB, plotW, plotH) {
    const dpr = cw / 800 > 1.6 ? 2 : 1;
    ctx.font = `${11 * dpr}px ui-monospace, Menlo, Consolas, monospace`;
    ctx.fillStyle = '#8fa3bd';
    ctx.strokeStyle = 'rgba(140,165,195,0.35)';
    ctx.lineWidth = 1;

    const kmX = g.Lx / 1000;
    const kmZ = g.H / 1000;

    // 高度刻度
    const zStep = kmZ > 10 ? 2 : kmZ > 3 ? 1 : 0.2;
    for (let z = 0; z <= kmZ + 1e-9; z += zStep) {
      const y = plotH - (z / kmZ) * plotH;
      ctx.beginPath(); ctx.moveTo(padL - 6 * dpr, y); ctx.lineTo(padL, y); ctx.stroke();
      ctx.textAlign = 'right';
      ctx.fillText(`${z.toFixed(zStep < 1 ? 1 : 0)}`, padL - 9 * dpr, y + 3.5 * dpr);
    }
    ctx.textAlign = 'left';
    ctx.fillText('km', 4 * dpr, 14 * dpr);

    // 水平刻度：按可用宽度取"整"步长，避免 0~60 km 每 1 km 一个标签糊成一片
    const maxTicks = Math.max(3, Math.floor(plotW / (58 * dpr)));
    const xStep = niceStep(kmX, maxTicks);
    for (let x = 0; x <= kmX + 1e-9; x += xStep) {
      const px = padL + (x / kmX) * plotW;
      ctx.beginPath(); ctx.moveTo(px, plotH); ctx.lineTo(px, plotH + 5 * dpr); ctx.stroke();
      ctx.textAlign = x === 0 ? 'left' : (Math.abs(x - kmX) < 1e-9 ? 'right' : 'center');
      ctx.fillText(`${x.toFixed(xStep < 1 ? 1 : 0)}`, px, plotH + 16 * dpr);
    }
  }

  /** 取 1/2/5×10ⁿ 的"整"步长，使刻度数不超过 maxTicks */
  function niceStep(span, maxTicks) {
    const raw = span / Math.max(1, maxTicks);
    const p = Math.pow(10, Math.floor(Math.log10(raw)));
    const m = raw / p;
    const f = m <= 1 ? 1 : m <= 2 ? 2 : m <= 5 ? 5 : 10;
    return f * p;
  }

  /** 画布坐标 → 模拟坐标（用于点击注入热泡） */
  function pick(model, clientX, clientY, showAxes = true) {
    const rect = canvas.getBoundingClientRect();
    const dpr = canvas.width / rect.width;
    const x = (clientX - rect.left) * dpr;
    const y = (clientY - rect.top) * dpr;
    const padL = showAxes ? Math.round(64 * (canvas.width / 800)) : 0;
    const padB = showAxes ? Math.round(34 * (canvas.width / 800)) : 0;
    const plotW = canvas.width - padL;
    const plotH = canvas.height - padB;
    const g = model.grid;
    return {
      x: Math.max(0, Math.min(g.Lx, ((x - padL) / plotW) * g.Lx)),
      z: Math.max(0, Math.min(g.H, (1 - y / plotH) * g.H)),
    };
  }

  return { draw, pick, resize, FIELD_LIST };
}
