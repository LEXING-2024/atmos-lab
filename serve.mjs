/**
 * 本地静态服务器（带 COOP/COEP 头）。
 *
 * 浏览器不允许从 file:// 加载 ES 模块，必须走 HTTP。
 * COOP/COEP 是为将来接 SharedArrayBuffer（多线程推进）预留的：
 * 人体工学的做法是现在就把头带上，以后加线程不用改部署方式。
 *
 * 用法：node serve.mjs [端口]
 */

import { createServer } from 'node:http';
import { readFile } from 'node:fs/promises';
import { extname, join, normalize, sep } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = fileURLToPath(new URL('.', import.meta.url));
const port = Number(process.argv[2] || process.env.PORT || 8787);

const MIME = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.mjs': 'text/javascript; charset=utf-8',
  '.css': 'text/css; charset=utf-8',
  '.json': 'application/json; charset=utf-8',
  '.md': 'text/markdown; charset=utf-8',
  '.png': 'image/png',
  '.svg': 'image/svg+xml',
  '.ico': 'image/x-icon',
};

const server = createServer(async (req, res) => {
  const url = new URL(req.url, 'http://localhost');
  let p = decodeURIComponent(url.pathname);
  if (p === '/' || p === '') p = '/index.html';

  // 只允许项目目录内的文件
  const target = normalize(join(root, p));
  if (!target.startsWith(root.endsWith(sep) ? root : root + sep)) {
    res.writeHead(403); res.end('403'); return;
  }

  try {
    const buf = await readFile(target);
    res.writeHead(200, {
      'Content-Type': MIME[extname(target).toLowerCase()] || 'application/octet-stream',
      'Cross-Origin-Opener-Policy': 'same-origin',
      'Cross-Origin-Embedder-Policy': 'require-corp',
      'Cross-Origin-Resource-Policy': 'same-origin',
      'Cache-Control': 'no-store',
    });
    res.end(buf);
  } catch {
    res.writeHead(404, { 'Content-Type': 'text/plain; charset=utf-8' });
    res.end('404 ' + p);
  }
});

server.listen(port, () => {
  console.log(`平面大气模拟框架已启动：http://localhost:${port}/`);
  console.log('按 Ctrl+C 停止。');
});
