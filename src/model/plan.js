/**
 * 性能探测与自动定档。
 *
 * 把"这台机器能跑多大的网格"变成一个可测量的量：
 * 吞吐量 T = 格点·步/毫秒。
 * 一帧能负担的格点数 = T · 帧预算 / 每帧步数。
 * 每帧步数由时间倍速决定（模拟秒/帧 ÷ Δt）。
 */

import { Grid } from '../core/grid.js';
import { AnelasticModel } from './anelastic.js';

const nowMs = () => (
  typeof performance !== 'undefined' && performance.now ? performance.now() : Date.now()
);

/** 四档场景：域尺寸、目标分辨率、建议倍速 */
export const SCALE_PRESETS = [
  {
    id: 'macro',
    name: '宏观',
    domain: '约 600 km',
    Lx: 600000, H: 16000,
    dxTarget: 9400, dzTarget: 500,
    levels: [
      { nx: 64, nz: 32, label: '64×32' },
      { nx: 128, nz: 48, label: '128×48' },
      { nx: 256, nz: 64, label: '256×64' },
    ],
    simSpeed: 60,
    note: '大尺度波动与热力结构；本模式无地转，不含锋面气旋',
    micro: false,
    // 默认给一条 60 km 尺度的地表感热带：均匀加热只产生湍流，
    // 有水平对比才驱动热岛/海陆风那样的中尺度环流。
    surface: { sensible: 250, latent: 0, xWidthKm: 60 },
    // 各尺度的场幅不同，色标必须跟着尺度走，否则 0.5 m/s 的中尺度环流会被 ±15 m/s 的色标压成一片白
    viz: { wMax: 4, uMax: 12, thMax: 4, qcMax: 1.0e-3, qrMax: 1.0e-3 },
  },
  {
    id: 'meso',
    name: '中尺度',
    domain: '约 200 km',
    Lx: 200000, H: 16000,
    dxTarget: 1600, dzTarget: 333,
    levels: [
      { nx: 64, nz: 32, label: '64×32' },
      { nx: 128, nz: 48, label: '128×48' },
      { nx: 256, nz: 64, label: '256×64' },
    ],
    simSpeed: 20,
    note: '海陆风、山谷环流、中尺度对流系统',
    micro: true,
    surface: { sensible: 200, latent: 0, xWidthKm: 30 },
    viz: { wMax: 10, uMax: 20, thMax: 5, qcMax: 1.6e-3, qrMax: 1.6e-3 },
  },
  {
    id: 'cloud',
    name: '云分辨',
    domain: '约 60 km',
    Lx: 60000, H: 16000,
    dxTarget: 470, dzTarget: 333,
    levels: [
      { nx: 128, nz: 48, label: '128×48' },
      { nx: 256, nz: 64, label: '256×64' },
      { nx: 512, nz: 96, label: '512×96' },
    ],
    simSpeed: 20,
    note: '深对流、雷暴单体、地形降水',
    micro: true,
    surface: { sensible: 0, latent: 0, xWidthKm: 4 },
    viz: { wMax: 25, uMax: 25, thMax: 6, qcMax: 2.2e-3, qrMax: 2.5e-3 },
  },
  {
    id: 'les',
    name: '微尺度',
    domain: '约 2.4 km',
    Lx: 2400, H: 1200,
    dxTarget: 19, dzTarget: 25,
    levels: [
      { nx: 128, nz: 48, label: '128×48' },
      { nx: 256, nz: 64, label: '256×64' },
      { nx: 512, nz: 96, label: '512×96' },
    ],
    simSpeed: 60,
    note: '街区绕流、街谷环流、坡面湍流（LES）',
    micro: false,
    surface: { sensible: 0, latent: 0, xWidthKm: 0.3 },
    viz: { wMax: 4, uMax: 6, thMax: 2, qcMax: 8e-4, qrMax: 8e-4 },
  },
];

export function presetById(id) {
  return SCALE_PRESETS.find((p) => p.id === id) || SCALE_PRESETS[2];
}

/** 由场景档 + 分辨率档构造网格 */
export function buildGrid(preset, levelIndex) {
  const level = preset.levels[Math.max(0, Math.min(preset.levels.length - 1, levelIndex))];
  return new Grid({ nx: level.nx, nz: level.nz, Lx: preset.Lx, H: preset.H });
}

/** 用网格与配置造一个一次性模型，供探测/演示使用 */
export function buildModel(grid, preset, overrides = {}) {
  const base = { TSfc: 299.5, RHsfc: 0.78, ...(overrides.base || {}) };
  const cfg = {
    filter: 0.02,
    micro: { enable: preset.micro !== false, ...(overrides.micro || {}) },
    turb: { surface: { sensible: 0, latent: 0 }, ...(overrides.turb || {}) },
    rad: { ...(overrides.rad || {}) },
    ...(overrides.top || {}),
  };
  return new AnelasticModel(grid, { ...cfg, base });
}

/**
 * 实测吞吐量（格点·步/毫秒）。
 * 用一次性模型跑若干步，取平均；不污染正在运行的模型。
 * @returns {{cellStepsPerMs:number, msPerStep:number, steps:number}}
 */
export function measureThroughput(preset, levelIndex, budgetMs = 60, overrides = {}) {
  const grid = buildGrid(preset, levelIndex);
  const probe = buildModel(grid, preset, overrides);
  probe.initAtmosphere({ amplitude: 0.5, radius: Math.max(500, grid.Lx / 20), noise: 0.01 });
  const dt = Math.max(probe.suggestDt(), 0.25);
  const t0 = nowMs();
  let steps = 0;
  while (nowMs() - t0 < budgetMs || steps < 3) {
    probe.step(dt);
    steps++;
    if (steps > 400) break;
  }
  const ms = nowMs() - t0;
  const cells = grid.nCell;
  return {
    cellStepsPerMs: (cells * steps) / ms,
    msPerStep: ms / steps,
    steps,
    cells,
    ms,
  };
}

/**
 * 自动定档：在给定帧预算与倍速下，选能跑得动的最大分辨率档。
 * @param {{cellStepsPerMs:number}} probe 实测吞吐
 * @param {object} preset 场景档
 * @param {{frameMs:number, simSpeed:number, dt:number}} want 期望
 */
export function autoPlan(probe, preset, want) {
  const frameMs = want.frameMs ?? 16;
  const simSpeed = want.simSpeed ?? preset.simSpeed;
  const dt = Math.max(want.dt ?? 1, 0.05);

  const simSecondsPerFrame = simSpeed * frameMs / 1000;
  const stepsPerFrame = Math.max(1, Math.ceil(simSecondsPerFrame / dt));
  const cellStepsPerFrame = probe.cellStepsPerMs * frameMs;
  const maxCells = cellStepsPerFrame / stepsPerFrame;

  const gridOf = (lv) => new Grid({ nx: lv.nx, nz: lv.nz, Lx: preset.Lx, H: preset.H });
  const ranked = preset.levels
    .map((lv, i) => ({ i, lv, cells: lv.nx * lv.nz }))
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

/** 倍速能否维持：给定实测吞吐与当前网格，估算每帧耗时 */
export function frameCost({ cells, cellStepsPerMs, dt, simSpeed, frameMs }) {
  const simSecondsPerFrame = simSpeed * frameMs / 1000;
  const stepsPerFrame = Math.max(1, Math.ceil(simSecondsPerFrame / dt));
  return {
    stepsPerFrame,
    msPerFrame: (cells / cellStepsPerMs) * stepsPerFrame,
  };
}
