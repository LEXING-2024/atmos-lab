/**
 * three.html 界面无头校验：用最小 DOM 桩真正 import src/three/main.js，
 * 校验 three.html 的 id 与代码取用一致，并跑若干帧确认无运行时错误。
 */

import { readFileSync } from 'node:fs';

const html = readFileSync(new URL('../three.html', import.meta.url), 'utf8');
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
    value: '0',
    checked: id === 'micro3',
    max: '0',
    dataset: {},
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
const RAF_LIMIT = 50;
globalThis.requestAnimationFrame = (cb) => {
  rafCount++;
  if (rafCount > RAF_LIMIT) return 0;
  setTimeout(() => cb(performance.now()), 3);
  return rafCount;
};

globalThis.document = {
  readyState: 'complete',
  getElementById(id) {
    if (!ids.has(id)) throw new Error(`three.html 缺少 id="${id}"`);
    return elStub(id);
  },
  createElement(tag) {
    if (tag === 'canvas') return elStub('__canvas');
    return elStub(tag);
  },
  querySelectorAll: () => [],
  addEventListener: () => {},
};
globalThis.location = { search: '?speed=180&level=0' };
globalThis.window = globalThis;
globalThis.devicePixelRatio = 1;
globalThis.HTMLCanvasElement = Object;

const t0 = Date.now();
await import('../src/three/main.js');
await new Promise((r) => setTimeout(r, 2500));

console.log(`初始化 + ${Math.min(rafCount, RAF_LIMIT)} 帧完成，耗时 ${((Date.now() - t0) / 1000).toFixed(2)} s`);
console.log(`three.html 声明的 id 共 ${ids.size} 个，全部被正确取用`);
