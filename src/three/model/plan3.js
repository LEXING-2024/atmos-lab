/**
 * 三维场景档与性能定档。
 * 吞吐量定义与二维版一致：格点·步/毫秒。
 */

import { Grid3 } from '../core/grid3.js';
import { AnelasticModel3 } from './anelastic3.js';

const nowMs = () => (
  typeof performance !== 'undefined' && performance.now ? performance.now() : Date.now()
);

/** 三维预设：域水平为正方形，垂直层数按 H 取 */
export const SCALE_PRESETS3 = [
  {
    id: 'conv3d',
    name: '3D 深对流',
    domain: '12 × 12 × 9 km',
    Lx: 12000, Ly: 12000, H: 9000,
    levels: [
      { nx: 32, ny: 32, nz: 24, label: '32²×24' },
      { nx: 64, ny: 64, nz: 32, label: '64²×32' },
      { nx: 128, ny: 64, nz: 32, label: '128×64×32' },
    ],
    simSpeed: 20,
    note: '单个对流单体的三维结构：上升气流、云砧、降水',
    micro: true,
    bubble: { amplitude: 3, radius: 2000, zCenter: 1200 },
    surface: { sensible: 0, latent: 0, xWidthKm: 2, yWidthKm: 2 },
    viz: { wMax: 25, uMax: 25, thMax: 6, qcMax: 2.2e-3, qrMax: 2.5e-3 },
  },
  {
    id: 'meso3d',
    name: '3D 中尺度',
    domain: '120 × 120 × 16 km',
    Lx: 120000, Ly: 120000, H: 16000,
    levels: [
      { nx: 32, ny: 32, nz: 24, label: '32²×24' },
      { nx: 64, ny: 64, nz: 32, label: '64²×32' },
    ],
    simSpeed: 60,
    note: '海陆风、山谷环流的三维环流结构（无地形，用热力对比代理）',
    micro: true,
    bubble: { amplitude: 2, radius: 12000, zCenter: 3000 },
    surface: { sensible: 200, latent: 0, xWidthKm: 30, yWidthKm: 30 },
    viz: { wMax: 10, uMax: 20, thMax: 5, qcMax: 1.6e-3, qrMax: 1.6e-3 },
  },
  {
    id: 'les3d',
    name: '3D 微尺度',
    domain: '2.4 × 2.4 × 1.2 km',
    Lx: 2400, Ly: 2400, H: 1200,
    levels: [
      { nx: 32, ny: 32, nz: 24, label: '32²×24' },
      { nx: 64, ny: 64, nz: 32, label: '64²×32' },
    ],
    simSpeed: 60,
    note: '热泡羽流、对流边界层的三维湍流结构（LES）',
    micro: false,
    bubble: { amplitude: 0.5, radius: 300, zCenter: 400 },
    surface: { sensible: 120, latent: 0, xWidthKm: 0.2, yWidthKm: 0.2 },
    viz: { wMax: 4, uMax: 6, thMax: 2, qcMax: 8e-4, qrMax: 8e-4 },
  },
];

export function presetById3(id) {
  return SCALE_PRESETS3.find((p) => p.id === id) || SCALE_PRESETS3[0];
}

export function buildGrid3(preset, levelIndex) {
  const level = preset.levels[Math.max(0, Math.min(preset.levels.length - 1, levelIndex))];
  return new Grid3({ nx: level.nx, ny: level.ny, nz: level.nz, Lx: preset.Lx, Ly: preset.Ly, H: preset.H });
}

export function buildModel3(grid, preset, overrides = {}) {
  const base = { TSfc: 299.5, RHsfc: 0.78, ...(overrides.base || {}) };
  const cfg = {
    filter: 0.02,
    micro: { enable: preset.micro !== false, ...(overrides.micro || {}) },
    turb: { surface: { sensible: 0, latent: 0 }, ...(overrides.turb || {}) },
    rad: { ...(overrides.rad || {}) },
    ...(overrides.top || {}),
  };
  return new AnelasticModel3(grid, { ...cfg, base });
}

/** 实测吞吐量（格点·步/毫秒），用一次性模型，不污染正在运行的模型 */
export function measureThroughput3(preset, levelIndex, budgetMs = 120, overrides = {}) {
  const grid = buildGrid3(preset, levelIndex);
  const probe = buildModel3(grid, preset, overrides);
  probe.initAtmosphere({ amplitude: 0.5, radius: Math.max(300, grid.Lx / 20), noise: 0.01 });
  const dt = Math.max(probe.suggestDt(), 0.5);
  const t0 = nowMs();
  let steps = 0;
  while (nowMs() - t0 < budgetMs || steps < 2) {
    probe.step(dt);
    steps++;
    if (steps > 60) break;
  }
  const ms = nowMs() - t0;
  return {
    cellStepsPerMs: (grid.nCell * steps) / ms,
    msPerStep: ms / steps,
    steps,
    cells: grid.nCell,
    ms,
  };
}

/**
 * 自动定档：给定帧预算与倍速，选能跑得动的最大分辨率档。
 */
export function autoPlan3(probe, preset, want) {
  const frameMs = want.frameMs ?? 33;
  const simSpeed = want.simSpeed ?? preset.simSpeed;
  const dt = Math.max(want.dt ?? 1, 0.05);

  const simSecondsPerFrame = simSpeed * frameMs / 1000;
  const stepsPerFrame = Math.max(1, Math.ceil(simSecondsPerFrame / dt));
  const cellStepsPerFrame = probe.cellStepsPerMs * frameMs;
  const maxCells = cellStepsPerFrame / stepsPerFrame;

  const gridOf = (lv) => new Grid3({ nx: lv.nx, ny: lv.ny, nz: lv.nz, Lx: preset.Lx, Ly: preset.Ly, H: preset.H });
  const ranked = preset.levels
    .map((lv, i) => ({ i, lv, cells: lv.nx * lv.ny * lv.nz }))
    .sort((a, b) => a.cells - b.cells);

  let chosen = ranked[0];
  for (const r of ranked) if (r.cells <= maxCells) chosen = r;

  return {
    levelIndex: chosen.i,
    label: chosen.lv.label,
    cells: chosen.cells,
    maxCells: Math.floor(maxCells),
    stepsPerFrame,
    estMsPerFrame: (chosen.cells / probe.cellStepsPerMs) * stepsPerFrame,
    estMsPerStep: chosen.cells / probe.cellStepsPerMs,
    grid: gridOf(chosen.lv),
  };
}

export function frameCost3({ cells, cellStepsPerMs, dt, simSpeed, frameMs }) {
  const simSecondsPerFrame = simSpeed * frameMs / 1000;
  const stepsPerFrame = Math.max(1, Math.ceil(simSecondsPerFrame / dt));
  return { stepsPerFrame, msPerFrame: (cells / cellStepsPerMs) * stepsPerFrame };
}
