/**
 * v1.5 WebGPU 试点：三维标量平流核（WGSL）与 CPU 逐点比对。
 *
 * 只搬"一个标量的平流倾向"这一个算子，不接管主循环：
 *   CPU 版：anelastic3.js 里的单调限制三阶迎风（x/y 周期、z 钳位）
 *   GPU 版：同一离散写成 compute shader，逐格点算一个倾向值
 *
 * 比对口径：同一组随机速度场 + 光滑标量场，max|gpu − cpu| 与相对误差。
 */

const WGSL = /* wgsl */ `
struct Params {
  nx: u32, ny: u32, nz: u32,
  dx: f32, dy: f32, dz: f32,
};
@group(0) @binding(0) var<storage, read> f: array<f32>;
@group(0) @binding(1) var<storage, read> vu: array<f32>;
@group(0) @binding(2) var<storage, read> vv: array<f32>;
@group(0) @binding(3) var<storage, read> vw: array<f32>;
@group(0) @binding(4) var<storage, read> rho0: array<f32>;
@group(0) @binding(5) var<storage, read> rhoW: array<f32>;
@group(0) @binding(6) var<storage, read_write> outT: array<f32>;
@group(0) @binding(7) var<uniform> P: Params;

fn idx(i: u32, j: u32, k: u32) -> u32 {
  return i + P.nx * (j + P.ny * k);
}

fn lim(v: f32, a: f32, b: f32) -> f32 {
  let lo = min(a, b);
  let hi = max(a, b);
  return clamp(v, lo, hi);
}

@compute @workgroup_size(4, 4, 4)
fn main(@builtin(global_invocation_id) g: vec3<u32>) {
  let i = g.x;
  let j = g.y;
  let k = g.z;
  if (i >= P.nx || j >= P.ny || k >= P.nz) { return; }

  let n = idx(i, j, k);
  let nXY = P.nx * P.ny;
  let ip = select(i + 1u, 0u, i + 1u == P.nx);
  let im = select(i - 1u, P.nx - 1u, i == 0u);
  let jp = select(j + 1u, 0u, j + 1u == P.ny);
  let jm = select(j - 1u, P.ny - 1u, j == 0u);

  // i-2 / j-2：i<=1 时回绕，否则直接减 2
  var im2: u32;
  if (i <= 1u) { im2 = (i + P.nx - 2u) % P.nx; } else { im2 = i - 2u; }
  var jm2: u32;
  if (j <= 1u) { jm2 = (j + P.ny - 2u) % P.ny; } else { jm2 = j - 2u; }

  let rk = rho0[k];

  // ---- x ----
  let vx = vu[n];
  let a = f[idx(im, j, k)];
  let b = f[n];
  let c = f[idx(ip, j, k)];
  let d = f[idx(im2, j, k)];
  var fqx: f32;
  if (vx > 0.0) { fqx = (6.0 * a + 3.0 * b - d) / 8.0; } else { fqx = (6.0 * b + 3.0 * a - c) / 8.0; }
  fqx = lim(fqx, a, b);
  let fxHere = rk * vx * fqx;

  let vxP = vu[idx(ip, j, k)];
  let aP = f[n];
  let bP = c;
  let cP = f[idx(select(ip + 1u, 0u, ip + 1u == P.nx), j, k)];
  let dP = f[idx(im, j, k)];
  var fqxP: f32;
  if (vxP > 0.0) { fqxP = (6.0 * aP + 3.0 * bP - dP) / 8.0; } else { fqxP = (6.0 * bP + 3.0 * aP - cP) / 8.0; }
  fqxP = lim(fqxP, aP, bP);
  let fxNext = rk * vxP * fqxP;

  // ---- y ----
  let vy = vv[n];
  let aY = f[idx(i, jm, k)];
  let cY = f[idx(i, jp, k)];
  let dY = f[idx(i, jm2, k)];
  var fqy: f32;
  if (vy > 0.0) { fqy = (6.0 * aY + 3.0 * b - dY) / 8.0; } else { fqy = (6.0 * b + 3.0 * aY - cY) / 8.0; }
  fqy = lim(fqy, aY, b);
  let fyHere = rk * vy * fqy;

  let vyP = vv[idx(i, jp, k)];
  let aYP = f[n];
  let bYP = cY;
  let cYP = f[idx(i, select(jp + 1u, 0u, jp + 1u == P.ny), k)];
  let dYP = f[idx(i, jm, k)];
  var fqyP: f32;
  if (vyP > 0.0) { fqyP = (6.0 * aYP + 3.0 * bYP - dYP) / 8.0; } else { fqyP = (6.0 * bYP + 3.0 * aYP - cYP) / 8.0; }
  fqyP = lim(fqyP, aYP, bYP);
  let fyNext = rk * vyP * fqyP;

  // ---- z ----
  var fB: f32 = 0.0;
  if (k > 0u) {
    let wB = vw[n];
    let cD = f[n - nXY];
    var cD2: f32;
    if (k > 1u) { cD2 = f[n - nXY * 2u]; } else { cD2 = f[n - nXY]; }
    var cU: f32;
    if (k + 1u < P.nz) { cU = f[n + nXY]; } else { cU = b; }
    var fqB: f32;
    if (wB > 0.0) { fqB = (6.0 * cD + 3.0 * b - cD2) / 8.0; } else { fqB = (6.0 * b + 3.0 * cD - cU) / 8.0; }
    fqB = lim(fqB, cD, b);
    fB = rhoW[k] * wB * fqB;
  }

  var fT: f32 = 0.0;
  if (k + 1u < P.nz) {
    let wT = vw[n + nXY];
    let cU = f[n + nXY];
    var cU2: f32;
    if (k + 2u < P.nz) { cU2 = f[n + nXY * 2u]; } else { cU2 = f[n + nXY]; }
    var cD: f32;
    if (k > 0u) { cD = f[n - nXY]; } else { cD = b; }
    var fqT: f32;    if (wT > 0.0) { fqT = (6.0 * b + 3.0 * cU - cD) / 8.0; } else { fqT = (6.0 * cU + 3.0 * b - cU2) / 8.0; }
    fqT = lim(fqT, b, cU);
    fT = rhoW[k + 1u] * wT * fqT;
  }

  let tend = -((fxNext - fxHere) / (P.dx * rk))
    - ((fyNext - fyHere) / (P.dy * rk))
    - ((fT - fB) / (P.dz * rk));

  outT[n] = tend;
}
`;

/** CPU 参考：与 anelastic3 的标量平流同一个离散 */
export function advectTendencyCPU({ nx, ny, nz, dx, dy, dz }, f, u, v, w, rho0, rhoW, out) {
  const nXY = nx * ny;
  const im2Of = (i) => (i <= 1 ? (i - 2 + nx) % nx : i - 2);
  const jm2Of = (j) => (j <= 1 ? (j - 2 + ny) % ny : j - 2);
  const at = (i, j, k) => i + nx * (j + ny * k);

  for (let k = 0; k < nz; k++) {
    const off = nXY * k;
    const rk = rho0[k];
    for (let j = 0; j < ny; j++) {
      for (let i = 0; i < nx; i++) {
        const n = off + nx * j + i;

        // x
        const vx = u[n];
        const ip = i + 1 === nx ? 0 : i + 1;
        const im = i === 0 ? nx - 1 : i - 1;
        const im2 = im2Of(i);
        const a = f[at(im, j, k)];
        const b = f[n];
        const c = f[at(ip, j, k)];
        const d = f[at(im2, j, k)];
        let fq = vx > 0 ? (6 * a + 3 * b - d) / 8 : (6 * b + 3 * a - c) / 8;
        const lo = a < b ? a : b, hi = a < b ? b : a;
        if (fq < lo) fq = lo; else if (fq > hi) fq = hi;
        const fHere = rk * vx * fq;

        const vxP = u[at(ip, j, k)];
        const aP = f[n];
        const bP = c;
        const cP = f[at(ip + 1 === nx ? 0 : ip + 1, j, k)];
        const dP = f[at(im, j, k)];
        let fqP = vxP > 0 ? (6 * aP + 3 * bP - dP) / 8 : (6 * bP + 3 * aP - cP) / 8;
        const loP = aP < bP ? aP : bP, hiP = aP < bP ? bP : aP;
        if (fqP < loP) fqP = loP; else if (fqP > hiP) fqP = hiP;
        const fNext = rk * vxP * fqP;

        // y
        const vy = v[n];
        const jp = j + 1 === ny ? 0 : j + 1;
        const jm = j === 0 ? ny - 1 : j - 1;
        const jm2 = jm2Of(j);
        const ay = f[at(i, jm, k)];
        const cy = f[at(i, jp, k)];
        const dy2 = f[at(i, jm2, k)];
        let fqy = vy > 0 ? (6 * ay + 3 * b - dy2) / 8 : (6 * b + 3 * ay - cy) / 8;
        const loy = ay < b ? ay : b, hiy = ay < b ? b : ay;
        if (fqy < loy) fqy = loy; else if (fqy > hiy) fqy = hiy;
        const fyHere = rk * vy * fqy;

        const vyP = v[at(i, jp, k)];
        const aYP = f[n];
        const bYP = cy;
        const cYP = f[at(i, jp + 1 === ny ? 0 : jp + 1, k)];
        const dYP = f[at(i, jm, k)];
        let fqyP = vyP > 0 ? (6 * aYP + 3 * bYP - dYP) / 8 : (6 * bYP + 3 * aYP - cYP) / 8;
        const loYP = aYP < bYP ? aYP : bYP, hiYP = aYP < bYP ? bYP : aYP;
        if (fqyP < loYP) fqyP = loYP; else if (fqyP > hiYP) fqyP = hiYP;
        const fyNext = rk * vyP * fqyP;

        // z
        let fB = 0;
        if (k > 0) {
          const wB = w[n];
          const cD = f[n - nXY];
          const cD2 = f[k > 1 ? n - nXY * 2 : n - nXY];
          const cU = k + 1 < nz ? f[n + nXY] : f[n];
          let fqB = wB > 0 ? (6 * cD + 3 * b - cD2) / 8 : (6 * b + 3 * cD - cU) / 8;
          const loB = cD < b ? cD : b, hiB = cD < b ? b : cD;
          if (fqB < loB) fqB = loB; else if (fqB > hiB) fqB = hiB;
          fB = rhoW[k] * wB * fqB;
        }
        let fT = 0;
        if (k + 1 < nz) {
          const wT = w[n + nXY];
          const cU = f[n + nXY];
          const cU2 = f[k + 2 < nz ? n + nXY * 2 : n + nXY];
          const cD = k > 0 ? f[n - nXY] : f[n];
          let fqT = wT > 0 ? (6 * b + 3 * cU - cD) / 8 : (6 * cU + 3 * b - cU2) / 8;
          const loT = b < cU ? b : cU, hiT = b < cU ? cU : b;
          if (fqT < loT) fqT = loT; else if (fqT > hiT) fqT = hiT;
          fT = rhoW[k + 1] * wT * fqT;
        }

        out[n] = -((fNext - fHere) / (dx * rk))
          - ((fyNext - fyHere) / (dy * rk))
          - ((fT - fB) / (dz * rk));
      }
    }
  }
}

/** 生成测试场：光滑标量 + 随机但连续的速度 */
function makeCase(nx, ny, nz) {
  const nXY = nx * ny, nCell = nx * ny * nz;
  const f = new Float32Array(nCell);
  const u = new Float32Array(nCell);
  const v = new Float32Array(nCell);
  const w = new Float32Array(nCell + nXY);
  const rho0 = new Float32Array(nz);
  const rhoW = new Float32Array(nz + 1);
  let seed = 1234567;
  const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff; };

  for (let k = 0; k < nz; k++) rho0[k] = 1.2 * Math.exp(-k / 8);
  for (let k = 0; k <= nz; k++) rhoW[k] = 1.2 * Math.exp(-k / 8);

  for (let k = 0; k < nz; k++) {
    for (let j = 0; j < ny; j++) {
      for (let i = 0; i < nx; i++) {
        const n = i + nx * (j + ny * k);
        const sx = Math.sin(2 * Math.PI * i / nx) * Math.cos(2 * Math.PI * j / ny);
        f[n] = 1 + 0.5 * sx * Math.sin(Math.PI * (k + 0.5) / nz);
        u[n] = 2 + 1.5 * Math.sin(2 * Math.PI * (j / ny + k / nz));
        v[n] = 1 + 1.2 * Math.cos(2 * Math.PI * (i / nx - k / nz));
      }
    }
  }
  for (let k = 0; k <= nz; k++) {
    for (let j = 0; j < ny; j++) {
      for (let i = 0; i < nx; i++) {
        const n = i + nx * (j + ny * k);
        w[n] = 0.6 * Math.sin(2 * Math.PI * (i / nx + j / ny)) * Math.cos(Math.PI * k / nz);
      }
    }
  }
  return { f, u, v, w, rho0, rhoW, nCell, nW: nCell + nXY };
}

export async function runAdvectionCompare(opts = {}) {
  const nx = opts.nx ?? 32, ny = opts.ny ?? 32, nz = opts.nz ?? 24;
  const dx = opts.dx ?? 375, dy = opts.dy ?? 375, dz = opts.dz ?? 375;
  const c = makeCase(nx, ny, nz);

  const outCPU = new Float32Array(c.nCell);
  // 预热 2 次（规避 JIT 冷启动），再取 7 次的中位数
  for (let r = 0; r < 2; r++) advectTendencyCPU({ nx, ny, nz, dx, dy, dz }, c.f, c.u, c.v, c.w, c.rho0, c.rhoW, outCPU);
  const samples = [];
  for (let r = 0; r < 7; r++) {
    const t0 = performance.now();
    advectTendencyCPU({ nx, ny, nz, dx, dy, dz }, c.f, c.u, c.v, c.w, c.rho0, c.rhoW, outCPU);
    samples.push(performance.now() - t0);
  }
  samples.sort((a, b) => a - b);
  const cpuMs = samples[3];

  if (!navigator.gpu) return { ok: false, reason: 'navigator.gpu 不存在（浏览器未开启 WebGPU）', cpuMs };
  const adapter = await navigator.gpu.requestAdapter();
  if (!adapter) return { ok: false, reason: '没有可用的 GPU adapter', cpuMs };
  const device = await adapter.requestDevice();

  const usageIn = GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST;
  const mk = (arr, usage = usageIn) => {
    const buf = device.createBuffer({ size: arr.byteLength, usage, mappedAtCreation: true });
    new Float32Array(buf.getMappedRange()).set(arr);
    buf.unmap();
    return buf;
  };
  const bufF = mk(c.f), bufU = mk(c.u), bufV = mk(c.v), bufW = mk(c.w);
  const bufR0 = mk(c.rho0), bufRW = mk(c.rhoW);
  const bufOut = device.createBuffer({
    size: c.nCell * 4,
    usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC | GPUBufferUsage.COPY_DST,
  });
  const params = new Float32Array([nx, ny, nz, dx, dy, dz]);   // u32 与 f32 混排：前 3 个当 u32 读
  const bufP = device.createBuffer({ size: 24, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST });
  device.queue.writeBuffer(bufP, 0, new Uint32Array([nx, ny, nz]));
  device.queue.writeBuffer(bufP, 12, new Float32Array([dx, dy, dz]));

  const layout = device.createBindGroupLayout({
    entries: [0, 1, 2, 3, 4, 5, 6].map((i) => ({
      binding: i, visibility: GPUShaderStage.COMPUTE,
      buffer: { type: i === 6 ? 'storage' : 'read-only-storage' },
    })).concat([{ binding: 7, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } }]),
  });
  const mod = device.createShaderModule({ code: WGSL });
  const ci = await mod.getCompilationInfo();
  const errs = ci.messages.filter((m) => m.type === 'error');
  if (errs.length) {
    throw new Error('WGSL 编译失败：' + errs.map((m) => `L${m.lineNum}:${m.linePos} ${m.message}`).join('；'));
  }
  const pipeline = device.createComputePipeline({
    layout: device.createPipelineLayout({ bindGroupLayouts: [layout] }),
    compute: { module: mod, entryPoint: 'main' },
  });
  const bind = device.createBindGroup({
    layout,
    entries: [
      { binding: 0, resource: { buffer: bufF } },
      { binding: 1, resource: { buffer: bufU } },
      { binding: 2, resource: { buffer: bufV } },
      { binding: 3, resource: { buffer: bufW } },
      { binding: 4, resource: { buffer: bufR0 } },
      { binding: 5, resource: { buffer: bufRW } },
      { binding: 6, resource: { buffer: bufOut } },
      { binding: 7, resource: { buffer: bufP } },
    ],
  });

  const wgX = Math.ceil(nx / 4), wgY = Math.ceil(ny / 4), wgZ = Math.ceil(nz / 4);
  const reps = 30;
  // 预热一次，再计时（首次 dispatch 含管线/着色器编译）
  {
    const enc0 = device.createCommandEncoder();
    const p0 = enc0.beginComputePass();
    p0.setPipeline(pipeline); p0.setBindGroup(0, bind); p0.dispatchWorkgroups(wgX, wgY, wgZ); p0.end();
    device.queue.submit([enc0.finish()]);
    await device.queue.onSubmittedWorkDone();
  }
  const tGpu = performance.now();
  for (let r = 0; r < reps; r++) {
    const enc = device.createCommandEncoder();
    const pass = enc.beginComputePass();
    pass.setPipeline(pipeline);
    pass.setBindGroup(0, bind);
    pass.dispatchWorkgroups(wgX, wgY, wgZ);
    pass.end();
    device.queue.submit([enc.finish()]);
  }
  await device.queue.onSubmittedWorkDone();
  const gpuMs = (performance.now() - tGpu) / reps;

  const staging = device.createBuffer({ size: c.nCell * 4, usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ });
  const enc2 = device.createCommandEncoder();
  enc2.copyBufferToBuffer(bufOut, 0, staging, 0, c.nCell * 4);
  device.queue.submit([enc2.finish()]);
  await staging.mapAsync(GPUMapMode.READ);
  const outGPU = new Float32Array(staging.getMappedRange()).slice();
  staging.unmap();

  let maxErr = 0, ref = 0, nonFinite = 0;
  for (let n = 0; n < c.nCell; n++) {
    const cv = outCPU[n], gv = outGPU[n];
    if (!Number.isFinite(cv) || !Number.isFinite(gv)) { nonFinite++; continue; }
    const e = Math.abs(gv - cv);
    if (e > maxErr) maxErr = e;
    ref += Math.abs(cv);
  }
  const meanAbs = ref / c.nCell;
  const relErr = meanAbs > 0 ? maxErr / meanAbs : (maxErr === 0 ? 0 : Infinity);

  device.destroy();
  return {
    ok: maxErr < 1e-4 && nonFinite === 0,
    maxErr, relErr, nonFinite, cpuMs, gpuMs, nx, ny, nz,
    cells: c.nCell,
    adapter: adapter.info ? `${adapter.info.vendor} ${adapter.info.architecture}` : 'unknown',
  };
}
