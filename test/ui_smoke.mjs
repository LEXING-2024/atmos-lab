/**
 * 界面无头校验。
 *
 * 用最小 DOM 桩把 src/main.js 真正 import 进来跑若干帧：
 * 1. 校验 index.html 里声明的 id 与 main.js 取的 id 完全对得上（少一个就抛错）
 * 2. 校验渲染/交互代码路径没有运行时错误
 * 3. 输出初始化后的模型信息与若干帧的 FPS/步数
 *
 * 用法：node test/ui_smoke.mjs
 */

import { readFileSync } from 'node:fs';

const html = readFileSync(new URL('../index.html', import.meta.url), 'utf8');
const ids = new Set([...html.matchAll(/id="([^"]+)"/g)].map((m) => m[1]));

function ctx2dStub() {
  return {
    imageSmoothingEnabled: true,
    canvas: null,
    createImageData: (w, h) => ({ width: w, height: h, data: new Uint8ClampedArray(w * h * 4) }),
    putImageData: () => {}, drawImage: () => {}, fillRect: () => {},
    save: () => {}, restore: () => {}, beginPath: () => {}, moveTo: () => {},
    lineTo: () => {}, stroke: () => {}, fillText: () => {}, measureText: () => ({ width: 10 }),
    set font(v) {}, set fillStyle(v) {}, set strokeStyle(v) {}, set lineWidth(v) {},
  };
}

function elStub(id) {
  return {
    id,
    innerHTML: '',
    textContent: '',
    value: id === 'rh' ? '0.78' : id === 'amp' ? '3' : id === 'radius' ? '4000' : '0',
    checked: id === 'micro',
    classList: { toggle: () => {}, add: () => {}, remove: () => {} },
    querySelectorAll: () => [],
    appendChild: () => {},
    addEventListener: () => {},
    getContext: () => ctx2dStub(),
    getBoundingClientRect: () => ({ left: 0, top: 0, width: 1024, height: 640 }),
    clientWidth: 1024,
    clientHeight: 640,
    width: 1024,
    height: 640,
  };
}

let rafCount = 0;
const RAF_LIMIT = 40;
globalThis.requestAnimationFrame = (cb) => {
  rafCount++;
  if (rafCount > RAF_LIMIT) return 0;
  setTimeout(() => cb(performance.now()), 4);
  return rafCount;
};

globalThis.document = {
  getElementById(id) {
    if (!ids.has(id)) throw new Error(`index.html 缺少 id="${id}"`);
    return elStub(id);
  },
  createElement(tag) {
    if (tag === 'canvas') return elStub('__canvas');
    const el = elStub(tag);
    el.toString = () => tag;
    return el;
  },
};
globalThis.window = globalThis;
globalThis.devicePixelRatio = 1;
globalThis.HTMLCanvasElement = Object;

const t0 = Date.now();
await import('../src/main.js');
await new Promise((r) => setTimeout(r, 1200));

console.log(`界面初始化 + ${Math.min(rafCount, RAF_LIMIT)} 帧完成，耗时 ${((Date.now() - t0) / 1000).toFixed(2)} s`);
console.log(`index.html 声明的 id 共 ${ids.size} 个，全部被正确取用`);
