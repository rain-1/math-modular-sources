#!/usr/bin/env node
"use strict";

const fs = require("fs/promises");
const path = require("path");
const http = require("http");

const PORT = Number(process.env.CHATGPT_ARCHIVER_DEBUG_PORT || 9222);
const OUTPUT = process.env.CHATGPT_ARCHIVE_OUTPUT || path.join(process.env.HOME, ".local/share/chatgpt-research-archive");
const PROJECT = process.env.CHATGPT_AUDIT_PROJECT || "apery systems";
const MIN_DELAY = Number(process.env.CHATGPT_AUDIT_MIN_DELAY || 10);
const MAX_DELAY = Number(process.env.CHATGPT_AUDIT_MAX_DELAY || 20);
const wait = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

const getJson = (route) => new Promise((resolve, reject) => {
  const request = http.get({ host: "127.0.0.1", port: PORT, path: route }, (response) => {
    let body = "";
    response.on("data", (chunk) => { body += chunk; });
    response.on("end", () => { try { resolve(JSON.parse(body)); } catch (error) { reject(error); } });
  });
  request.on("error", reject);
});

class Cdp {
  constructor(ws) { this.ws = ws; this.id = 0; this.pending = new Map(); }
  static async connect() {
    const target = (await getJson("/json/list")).find((item) => item.type === "page" && item.url.startsWith("https://chatgpt.com"));
    if (!target) throw new Error("No ChatGPT debug page is available");
    const ws = new WebSocket(target.webSocketDebuggerUrl);
    await new Promise((resolve, reject) => { ws.onopen = resolve; ws.onerror = reject; });
    const cdp = new Cdp(ws);
    ws.onmessage = (event) => {
      const message = JSON.parse(event.data);
      const pending = cdp.pending.get(message.id);
      if (!pending) return;
      clearTimeout(pending.timer); cdp.pending.delete(message.id);
      message.error ? pending.reject(new Error(message.error.message)) : pending.resolve(message.result);
    };
    return cdp;
  }
  call(method, params = {}, timeout = 60000) {
    const id = ++this.id;
    return new Promise((resolve, reject) => {
      const timer = setTimeout(() => { this.pending.delete(id); reject(new Error(`${method} timed out`)); }, timeout);
      this.pending.set(id, { resolve, reject, timer });
      this.ws.send(JSON.stringify({ id, method, params }));
    });
  }
  async eval(expression) {
    const result = await this.call("Runtime.evaluate", { expression, returnByValue: true, awaitPromise: true });
    if (result.exceptionDetails) throw new Error(result.exceptionDetails.text || "Evaluation failed");
    return result.result.value;
  }
  async navigate(url) {
    await this.call("Page.navigate", { url });
    for (let attempt = 0; attempt < 120; attempt++) {
      await wait(250);
      if ((await this.eval("document.readyState")) === "complete") break;
    }
    for (let attempt = 0; attempt < 60; attempt++) {
      if (await this.eval(`document.querySelectorAll('[data-testid^="conversation-turn-"]').length`)) return;
      await wait(500);
    }
  }
}

async function listJson(directory) {
  const result = [];
  for (const entry of await fs.readdir(directory, { withFileTypes: true })) {
    const full = path.join(directory, entry.name);
    if (entry.isDirectory()) result.push(...await listJson(full));
    else if (entry.name.endsWith(".json")) result.push(full);
  }
  return result;
}

async function liveTurns(cdp) {
  await cdp.eval(`(async()=>{const wait=ms=>new Promise(r=>setTimeout(r,ms));let last='',stable=0;for(let i=0;i<80&&stable<5;i++){const turn=document.querySelector('[data-testid^="conversation-turn-"]');let s=turn?.parentElement;while(s&&s!==document.body){if(/auto|scroll/.test(getComputedStyle(s).overflowY)&&s.scrollHeight>s.clientHeight)break;s=s.parentElement}s=s||document.scrollingElement;s.scrollTop=0;await wait(600);const sig=[...document.querySelectorAll('[data-testid^="conversation-turn-"]')].map(x=>x.getAttribute('data-testid')).join('|')+':'+s.scrollHeight;stable=sig===last?stable+1:0;last=sig}})()`);
  return JSON.parse(await cdp.eval(`JSON.stringify({turns:[...document.querySelectorAll('[data-testid^="conversation-turn-"]')].map(a=>a.innerText.replace(/\\r/g,'').replace(/\\n{3,}/g,'\\n\\n').trim()),ids:[...document.querySelectorAll('[data-testid^="conversation-turn-"]')].map(a=>a.getAttribute('data-testid')),blocked:/too many requests|rate limit|unusual activity|verify you are human|checking your browser/i.test(document.body.innerText)})`));
}

(async () => {
  const state = JSON.parse(await fs.readFile(path.join(OUTPUT, "archive-state.json"), "utf8"));
  const items = state.conversations.filter((item) => item.project === PROJECT && item.status === "done");
  const records = new Map();
  for (const filename of await listJson(path.join(OUTPUT, "conversations"))) {
    const record = JSON.parse(await fs.readFile(filename, "utf8"));
    records.set(record.url, { filename, record });
  }
  const reportFile = path.join(OUTPUT, `audit-${PROJECT.replace(/[^a-z0-9]+/gi, "-").toLowerCase()}.json`);
  let report;
  try { report = JSON.parse(await fs.readFile(reportFile, "utf8")); }
  catch { report = { project: PROJECT, startedAt: new Date().toISOString(), results: [] }; }
  const completed = new Set(report.results.filter((item) => item.verified).map((item) => item.url));
  const cdp = await Cdp.connect();
  for (const [index, item] of items.entries()) {
    if (completed.has(item.url)) continue;
    console.log(`[${index + 1}/${items.length}] Auditing: ${item.title}`);
    const saved = records.get(item.url);
    if (!saved) throw new Error(`Saved JSON not found: ${item.url}`);
    await cdp.navigate(item.url);
    const live = await liveTurns(cdp);
    if (live.blocked) throw new Error("ChatGPT displayed a rate-limit or access warning; audit stopped");
    const stored = saved.record.messages.map((message) => message.text);
    const result = {
      url: item.url,
      title: item.title,
      savedFile: saved.filename,
      auditedAt: new Date().toISOString(),
      savedTurns: stored.length,
      liveTurns: live.turns.length,
      firstMatches: stored[0] === live.turns[0],
      lastMatches: stored.at(-1) === live.turns.at(-1),
      allTextMatches: stored.length === live.turns.length && stored.every((text, turn) => text === live.turns[turn]),
      liveTurnIds: live.ids
    };
    result.verified = result.savedTurns === result.liveTurns && result.firstMatches && result.lastMatches && result.allTextMatches;
    report.results = report.results.filter((entry) => entry.url !== item.url);
    report.results.push(result);
    report.updatedAt = new Date().toISOString();
    await fs.writeFile(reportFile, JSON.stringify(report, null, 2));
    console.log(`  ${result.verified ? "VERIFIED" : "MISMATCH"}: ${result.savedTurns} saved / ${result.liveTurns} live`);
    if (!result.verified) throw new Error(`Completeness mismatch: ${item.title}`);
    if (index < items.length - 1) await wait((MIN_DELAY + Math.floor(Math.random() * Math.max(1, MAX_DELAY - MIN_DELAY + 1))) * 1000);
  }
  report.completedAt = new Date().toISOString();
  report.verified = report.results.length === items.length && report.results.every((item) => item.verified);
  await fs.writeFile(reportFile, JSON.stringify(report, null, 2));
  console.log(`Audit complete: ${report.results.filter((item) => item.verified).length}/${items.length} verified`);
  cdp.ws.close();
})().catch((error) => { console.error(error.stack || error.message); process.exitCode = 1; });
