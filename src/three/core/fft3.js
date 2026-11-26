/**
 * 水平二维 FFT（x、y 各做一次一维 FFT，按层执行）。
 *
 * 只服务压力求解：输入是 (i,j,k) 布局的实数组，输出同布局的复谱。
 * 归一化由两次逆变换各除以 n 完成，总计 1/(nx·ny)。
 */

import { FFT1D } from '../../core/fft.js';

export function createHorizontalFFT(grid) {
  const { nx, ny, nz } = grid;
  const fftX = new FFT1D(nx);
  const fftY = new FFT1D(ny);
  const nXY = nx * ny;
  const scratchRe = new Float64Array(Math.max(nx, ny));
  const scratchIm = new Float64Array(Math.max(nx, ny));

  function forward(re, im) {
    for (let k = 0; k < nz; k++) {
      const base = nXY * k;

      // 沿 x：每行连续
      for (let j = 0; j < ny; j++) {
        const off = base + nx * j;
        for (let i = 0; i < nx; i++) { scratchRe[i] = re[off + i]; scratchIm[i] = im[off + i]; }
        fftX.transform(scratchRe, scratchIm, false);
        for (let i = 0; i < nx; i++) { re[off + i] = scratchRe[i]; im[off + i] = scratchIm[i]; }
      }

      // 沿 y：步长为 nx
      for (let i = 0; i < nx; i++) {
        for (let j = 0; j < ny; j++) { scratchRe[j] = re[base + nx * j + i]; scratchIm[j] = im[base + nx * j + i]; }
        fftY.transform(scratchRe, scratchIm, false);
        for (let j = 0; j < ny; j++) { re[base + nx * j + i] = scratchRe[j]; im[base + nx * j + i] = scratchIm[j]; }
      }
    }
  }

  function inverse(re, im) {
    for (let k = 0; k < nz; k++) {
      const base = nXY * k;

      for (let i = 0; i < nx; i++) {
        for (let j = 0; j < ny; j++) { scratchRe[j] = re[base + nx * j + i]; scratchIm[j] = im[base + nx * j + i]; }
        fftY.transform(scratchRe, scratchIm, true);
        for (let j = 0; j < ny; j++) { re[base + nx * j + i] = scratchRe[j]; im[base + nx * j + i] = scratchIm[j]; }
      }

      for (let j = 0; j < ny; j++) {
        const off = base + nx * j;
        for (let i = 0; i < nx; i++) { scratchRe[i] = re[off + i]; scratchIm[i] = im[off + i]; }
        fftX.transform(scratchRe, scratchIm, true);
        for (let i = 0; i < nx; i++) { re[off + i] = scratchRe[i]; im[off + i] = scratchIm[i]; }
      }
    }
  }

  return { fftX, fftY, forward, inverse };
}
