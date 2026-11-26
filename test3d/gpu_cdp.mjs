/** 通过 CDP 读取 gpu3.html 的运行结果（绕过 headless 的虚拟时间不等待 GPU 异步） */
const port = process.argv[2] ? Number(process.argv[2]) : 9333;
const url = `http://127.0.0.1:${port}/json/list`;

const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

async function main() {
  let list;
  for (let i = 0; i < 30; i++) {
    try {
      list = await (await fetch(url)).json();
      if (list.length) break;
    } catch { /* retry */ }
    await sleep(500);
  }
  const page = list?.find((t) => t.type === 'page' && t.url.includes('gpu3'));
  if (!page) { console.log('找不到 gpu3 页面'); process.exit(1); }

  const ws = new WebSocket(page.webSocketDebuggerUrl);
  await new Promise((res, rej) => { ws.onopen = res; ws.onerror = rej; });

  let id = 0;
  const send = (method, params) => new Promise((res) => {
    const mid = ++id;
    const onMsg = (ev) => {
      const m = JSON.parse(ev.data);
      if (m.id === mid) { ws.removeEventListener('message', onMsg); res(m.result); }
    };
    ws.addEventListener('message', onMsg);
    ws.send(JSON.stringify({ id: mid, method, params }));
  });

  await send('Runtime.enable', {});
  for (let i = 0; i < 20; i++) {
    const r = await send('Runtime.evaluate', {
      expression: "document.getElementById('out') ? document.getElementById('out').textContent : 'no-pre'",
      returnByValue: true,
    });
    const txt = r?.result?.value ?? '';
    if (txt && txt !== '运行中…' && txt !== 'pending' && txt !== 'start') {
      console.log(txt);
      ws.close();
      return;
    }
    await sleep(1000);
  }
  console.log('超时：页面没有产出结果（可能 GPU 设备创建失败）');
  ws.close();
}

main().catch((e) => { console.error('CDP 失败:', e.message); process.exit(1); });
