#!/usr/bin/env node
"use strict";

const fs = require("fs/promises");
const path = require("path");
const http = require("http");

const PORT = Number(process.env.CHATGPT_ARCHIVER_DEBUG_PORT || 9222);
const OUTPUT = process.env.CHATGPT_ARCHIVE_OUTPUT || path.join(process.env.HOME, ".local/share/chatgpt-research-archive");
const MIN_DELAY = Number(process.env.CHATGPT_ARCHIVE_MIN_DELAY || 90);
const MAX_DELAY = Number(process.env.CHATGPT_ARCHIVE_MAX_DELAY || 180);
const MAX_ITEMS = Number(process.env.CHATGPT_ARCHIVE_MAX_ITEMS || 0);
const RATE_LIMIT_BACKOFF = Number(process.env.CHATGPT_ARCHIVE_RATE_LIMIT_BACKOFF || 3600);
const PROJECTS = [
  ["apery systems", "https://chatgpt.com/g/g-p-6a829c86de1081919010a9c156a218cf-apery-systems/project"],
  ["apery systems 4", "https://chatgpt.com/g/g-p-6a869e8cb0808191877658a38862861b-apery-systems-4/project"],
  ["apery systems 2", "https://chatgpt.com/g/g-p-6a82b0e7d91081918a3f6b61be2a9443-apery-systems-2/project"],
  ["apery systems p=3", "https://chatgpt.com/g/g-p-6a854b0b2a7081918bd6e18e8d13516e-apery-systems-p-3/project"],
  ["zeta", "https://chatgpt.com/g/g-p-6a801f369f748191955b34c90b30d548-zeta/project"]
];

const wait = (ms) => new Promise((resolve) => setTimeout(resolve, ms));
const slug = (value) => value.normalize("NFKD").replace(/[^a-zA-Z0-9]+/g, "-").replace(/^-|-$/g, "").slice(0, 80).toLowerCase() || "untitled";
const getJson = (route) => new Promise((resolve, reject) => http.get({ host: "127.0.0.1", port: PORT, path: route }, (response) => {
  let body = "";
  response.on("data", (chunk) => { body += chunk; });
  response.on("end", () => { try { resolve(JSON.parse(body)); } catch (error) { reject(error); } });
}).on("error", reject));

class Cdp {
  constructor(ws) { this.ws = ws; this.id = 0; this.pending = new Map(); }
  static async connect() {
    const targets = await getJson("/json/list");
    const target = targets.find((item) => item.type === "page" && item.url.startsWith("https://chatgpt.com"));
    if (!target) throw new Error("No authenticated ChatGPT page is open on the debug port");
    const ws = new WebSocket(target.webSocketDebuggerUrl);
    await new Promise((resolve, reject) => { ws.onopen = resolve; ws.onerror = reject; });
    const cdp = new Cdp(ws);
    ws.onmessage = (event) => {
      const message = JSON.parse(event.data);
      if (!message.id || !cdp.pending.has(message.id)) return;
      const { resolve, reject, timer } = cdp.pending.get(message.id);
      clearTimeout(timer); cdp.pending.delete(message.id);
      message.error ? reject(new Error(message.error.message)) : resolve(message.result);
    };
    return cdp;
  }
  call(method, params = {}, timeout = 30000) {
    const id = ++this.id;
    return new Promise((resolve, reject) => {
      const timer = setTimeout(() => { this.pending.delete(id); reject(new Error(`${method} timed out`)); }, timeout);
      this.pending.set(id, { resolve, reject, timer });
      this.ws.send(JSON.stringify({ id, method, params }));
    });
  }
  async evaluate(expression) {
    const result = await this.call("Runtime.evaluate", { expression, returnByValue: true, awaitPromise: true }, 60000);
    if (result.exceptionDetails) throw new Error(result.exceptionDetails.text || "Page evaluation failed");
    return result.result.value;
  }
  async navigate(url) {
    await this.call("Page.navigate", { url });
    for (let attempt = 0; attempt < 120; attempt++) {
      await wait(250);
      if ((await this.evaluate("document.readyState")) === "complete") break;
    }
    await wait(2000);
  }
  close() { this.ws.close(); }
}

async function loadState() {
  try { return JSON.parse(await fs.readFile(path.join(OUTPUT, "archive-state.json"), "utf8")); }
  catch { return { schemaVersion: 1, projects: [], conversations: [], createdAt: new Date().toISOString() }; }
}

async function saveState(state) {
  state.updatedAt = new Date().toISOString();
  await fs.mkdir(OUTPUT, { recursive: true });
  await fs.writeFile(path.join(OUTPUT, "archive-state.json"), JSON.stringify(state, null, 2));
}

async function discover(cdp, state) {
  const known = new Map(state.conversations.map((item) => [item.url, item]));
  for (const [projectTitle, projectUrl] of PROJECTS) {
    console.log(`Discovering: ${projectTitle}`);
    await cdp.navigate(projectUrl);
    const result = JSON.parse(await cdp.evaluate(`(async()=>{
      const wait=ms=>new Promise(r=>setTimeout(r,ms));
      const found=new Map();
      for(let page=0;page<200;page++){
        for(const a of document.querySelectorAll('main a[href*="/c/"]')){
          const url=new URL(a.getAttribute('href'),location.origin).href;
          const title=(a.innerText||a.getAttribute('aria-label')||'Untitled').trim().split('\\n')[0];
          found.set(url,{url,title});
        }
        const more=[...document.querySelectorAll('main button')].find(b=>b.innerText.trim()==='Load more conversations');
        if(!more||more.disabled)break;
        const before=found.size;more.click();
        for(let i=0;i<40;i++){await wait(250);for(const a of document.querySelectorAll('main a[href*="/c/"]'))found.set(new URL(a.getAttribute('href'),location.origin).href,{url:new URL(a.getAttribute('href'),location.origin).href,title:(a.innerText||a.getAttribute('aria-label')||'Untitled').trim().split('\\n')[0]});if(found.size>before)break;}
      }
      return JSON.stringify([...found.values()]);
    })()`));
    for (const chat of result) {
      const existing = known.get(chat.url);
      if (existing) { existing.project = projectTitle; existing.title ||= chat.title; }
      else { const item = { ...chat, project: projectTitle, status: "pending", addedAt: new Date().toISOString() }; state.conversations.push(item); known.set(chat.url, item); }
    }
    const project = state.projects.find((item) => item.url === projectUrl) || { title: projectTitle, url: projectUrl };
    project.discoveredAt = new Date().toISOString(); project.chatCount = result.length;
    if (!state.projects.includes(project)) state.projects.push(project);
    await saveState(state);
    console.log(`  ${result.length} conversations`);
  }
}

async function capture(cdp, item) {
  await cdp.navigate(item.url);
  for (let attempt = 0; attempt < 60; attempt++) {
    if (await cdp.evaluate(`document.querySelectorAll('[data-testid^="conversation-turn-"]').length`)) break;
    await wait(500);
  }
  const blocked = await cdp.evaluate(`(()=>{const t=document.body.innerText.toLowerCase();if(/too many requests|rate limit|unusual activity/.test(t))return 'rate-limit';if(/checking your browser|verify you are human|cloudflare/.test(t))return 'challenge';if(!document.querySelector('[data-testid^="conversation-turn-"]')&&/log in|sign up/.test(t))return 'login';return null})()`);
  if (blocked) throw new Error(`ChatGPT blocking state: ${blocked}`);
  await cdp.evaluate(`(async()=>{const wait=ms=>new Promise(r=>setTimeout(r,ms));let last='',stable=0;for(let i=0;i<80&&stable<5;i++){const turn=document.querySelector('[data-testid^="conversation-turn-"]');let s=turn?.parentElement;while(s&&s!==document.body){const y=getComputedStyle(s).overflowY;if(/auto|scroll/.test(y)&&s.scrollHeight>s.clientHeight)break;s=s.parentElement}s=s||document.scrollingElement;s.scrollTop=0;await wait(600);const sig=document.querySelectorAll('[data-testid^="conversation-turn-"]').length+':'+s.scrollHeight;stable=sig===last?stable+1:0;last=sig}})()`);
  const record = JSON.parse(await cdp.evaluate(`(()=>{const clean=e=>{const c=e.cloneNode(true);c.querySelectorAll('script,style,button,form,textarea,input,svg').forEach(n=>n.remove());for(const n of [c,...c.querySelectorAll('*')])for(const a of [...n.attributes])if(a.name.startsWith('on')||['nonce','integrity','srcset'].includes(a.name))n.removeAttribute(a.name);return c.innerHTML.trim()};const arts=[...document.querySelectorAll('[data-testid^="conversation-turn-"]')];const messages=arts.map((a,i)=>{const id=a.getAttribute('data-testid')||'';const role=id.includes('user')||a.querySelector('[data-message-author-role="user"]')?'user':'assistant';return{index:i+1,role,text:a.innerText.replace(/\\r/g,'').replace(/\\n{3,}/g,'\\n\\n').trim(),html:clean(a)}});return JSON.stringify({schemaVersion:1,url:location.href,title:document.title.replace(/\\s*[-–|]\\s*ChatGPT\\s*$/i,'').trim(),capturedAt:new Date().toISOString(),messages})})()`));
  if (!record.messages.length) throw new Error("No conversation turns found");
  return record;
}

async function writeRecord(item, record) {
  const id = item.url.match(/\/c\/([^/?#]+)/)?.[1] || "unknown";
  const directory = path.join(OUTPUT, "conversations", slug(item.project));
  const base = `${slug(item.title)}-${id.slice(0, 8)}`;
  await fs.mkdir(directory, { recursive: true });
  const markdown = [`# ${record.title}`, "", `- Project: ${item.project}`, `- Source: ${record.url}`, `- Captured: ${record.capturedAt}`, "", ...record.messages.flatMap((m) => [`## ${m.role === "user" ? "User" : "Assistant"} ${m.index}`, "", m.text, ""])].join("\n");
  const sections = record.messages.map((m) => `<section class="${m.role}"><h2>${m.role} ${m.index}</h2>${m.html}</section>`).join("\n");
  const html = `<!doctype html><html><head><meta charset="utf-8"><title>${record.title.replace(/[<>&]/g, "")}</title><style>body{font:16px/1.55 system-ui;max-width:900px;margin:auto;padding:2rem}section{padding:1rem;margin:1rem 0;background:#f4f4f4;border-radius:10px}.user{background:#e8f0fe}pre{overflow:auto;background:#222;color:#eee;padding:1rem}img{max-width:100%}</style></head><body><h1>${record.title.replace(/[<>&]/g, "")}</h1>${sections}</body></html>`;
  await Promise.all([
    fs.writeFile(path.join(directory, `${base}.json`), JSON.stringify(record, null, 2)),
    fs.writeFile(path.join(directory, `${base}.md`), markdown),
    fs.writeFile(path.join(directory, `${base}.html`), html)
  ]);
}

async function archive(cdp, state) {
  let processed = 0;
  for (const item of state.conversations.filter((entry) => entry.status !== "done")) {
    if (MAX_ITEMS && processed >= MAX_ITEMS) break;
    console.log(`Capturing [${item.project}]: ${item.title}`);
    item.status = "running"; item.attempts = (item.attempts || 0) + 1; await saveState(state);
    try {
      const record = await capture(cdp, item);
      await writeRecord(item, record);
      item.status = "done"; item.capturedAt = record.capturedAt; item.messageCount = record.messages.length; delete item.error;
      console.log(`  saved ${record.messages.length} messages`);
    } catch (error) {
      item.status = "error"; item.error = error.message; console.error(`  ERROR: ${error.message}`);
      if (/rate-limit/.test(error.message)) {
        item.status = "pending";
        await saveState(state);
        console.log(`  rate limit cooldown: waiting ${RATE_LIMIT_BACKOFF}s without requests`);
        await wait(RATE_LIMIT_BACKOFF * 1000);
        return archive(cdp, state);
      }
      if (/challenge|login/.test(error.message)) { await saveState(state); throw error; }
    }
    await saveState(state);
    processed++;
    if (MAX_ITEMS && processed >= MAX_ITEMS) break;
    const delay = MIN_DELAY + Math.floor(Math.random() * Math.max(1, MAX_DELAY - MIN_DELAY + 1));
    console.log(`  waiting ${delay}s`); await wait(delay * 1000);
  }
}

(async () => {
  await fs.mkdir(OUTPUT, { recursive: true });
  const state = await loadState();
  for (const item of state.conversations) if (item.status === "running") item.status = "pending";
  const cdp = await Cdp.connect();
  try {
    if (!state.projects.length || process.argv.includes("--rediscover")) await discover(cdp, state);
    if (!process.argv.includes("--discover-only")) await archive(cdp, state);
    await saveState(state);
    console.log(`Archive complete: ${OUTPUT}`);
  } finally { cdp.close(); }
})().catch((error) => { console.error(error.stack || error.message); process.exitCode = 1; });
