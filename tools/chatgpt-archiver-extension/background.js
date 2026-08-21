const DEFAULT_STATE = {
  queue: [],
  projects: [],
  running: false,
  workerTabId: null,
  currentUrl: null,
  currentProjectUrl: null,
  options: { delaySeconds: 90, formats: { markdown: true, html: true, json: true } }
};

async function getState() {
  const saved = await chrome.storage.local.get(DEFAULT_STATE);
  return { ...DEFAULT_STATE, ...saved, options: { ...DEFAULT_STATE.options, ...(saved.options || {}) } };
}

async function setState(changes) { await chrome.storage.local.set(changes); }

async function recoverInterruptedWork() {
  const state = await getState();
  for (const item of state.queue) {
    if (item.status === "running") item.status = "pending";
  }
  for (const project of state.projects) {
    if (project.status === "running") project.status = "pending";
  }
  await setState({ queue: state.queue, projects: state.projects, running: false, workerTabId: null, currentUrl: null, currentProjectUrl: null });
}

function normalizedConversationUrl(value) {
  try {
    const url = new URL(value);
    if (url.hostname !== "chatgpt.com") return null;
    const match = url.pathname.match(/\/c\/([a-zA-Z0-9-]+)/);
    return match ? `https://chatgpt.com/c/${match[1]}` : null;
  } catch { return null; }
}

function normalizedProjectUrl(value) {
  try {
    const url = new URL(value);
    if (url.hostname !== "chatgpt.com") return null;
    const isProject = /g-p-[^/]+/i.test(url.pathname) || /\/projects?\/[^/]+\/?$/i.test(url.pathname);
    if (!isProject) return null;
    url.search = "";
    url.hash = "";
    return url.href.replace(/\/$/, "");
  } catch { return null; }
}

async function addUrls(values) {
  const state = await getState();
  const existing = new Set(state.queue.map((item) => item.url));
  let added = 0;
  for (const value of values || []) {
    const url = normalizedConversationUrl(value);
    if (!url || existing.has(url)) continue;
    state.queue.push({ url, status: "pending", addedAt: new Date().toISOString(), attempts: 0 });
    existing.add(url);
    added++;
  }
  await setState({ queue: state.queue });
  return added;
}

async function addProjects(values) {
  const state = await getState();
  const existing = new Set(state.projects.map((item) => item.url));
  let added = 0;
  for (const value of values || []) {
    const url = normalizedProjectUrl(typeof value === "string" ? value : value?.url);
    if (!url || existing.has(url)) continue;
    state.projects.push({ url, title: typeof value === "object" ? value.title : undefined, status: "pending", addedAt: new Date().toISOString(), attempts: 0 });
    existing.add(url);
    added++;
  }
  await setState({ projects: state.projects });
  return added;
}

async function processNext() {
  const state = await getState();
  if (!state.running || state.currentUrl || state.currentProjectUrl) return;
  const project = state.projects.find((item) => item.status === "pending");
  if (project) {
    project.status = "running";
    project.attempts = (project.attempts || 0) + 1;
    let projectTabId = state.workerTabId;
    try {
      if (projectTabId) await chrome.tabs.update(projectTabId, { url: project.url, active: false });
      else projectTabId = (await chrome.tabs.create({ url: project.url, active: false })).id;
      await setState({ projects: state.projects, workerTabId: projectTabId, currentProjectUrl: project.url });
    } catch (error) {
      project.status = "error";
      project.error = error.message;
      await setState({ projects: state.projects, workerTabId: null, currentProjectUrl: null });
      scheduleNext(10);
    }
    return;
  }
  const next = state.queue.find((item) => item.status === "pending");
  if (!next) { await setState({ running: false, workerTabId: null }); return; }

  next.status = "running";
  next.attempts = (next.attempts || 0) + 1;
  next.startedAt = new Date().toISOString();
  let tabId = state.workerTabId;
  try {
    if (tabId) {
      await chrome.tabs.update(tabId, { url: next.url, active: false });
    } else {
      const tab = await chrome.tabs.create({ url: next.url, active: false });
      tabId = tab.id;
    }
    await setState({ queue: state.queue, workerTabId: tabId, currentUrl: next.url });
  } catch (error) {
    next.status = "error";
    next.error = error.message;
    await setState({ queue: state.queue, currentUrl: null, workerTabId: null });
    scheduleNext(10);
  }
}

function scheduleNext(seconds) {
  chrome.alarms.create("archive-next", { delayInMinutes: Math.max(seconds, 6) / 60 });
}

async function beginCapture(tabId) {
  const state = await getState();
  if (!state.running || tabId !== state.workerTabId || (!state.currentUrl && !state.currentProjectUrl)) return;
  try {
    if (state.currentProjectUrl) await chrome.tabs.sendMessage(tabId, { type: "CAPTURE_PROJECT", expectedUrl: state.currentProjectUrl });
    else await chrome.tabs.sendMessage(tabId, { type: "CAPTURE_CONVERSATION", expectedUrl: state.currentUrl });
  } catch {
    scheduleNext(8);
  }
}

chrome.tabs.onUpdated.addListener((tabId, changeInfo) => {
  if (changeInfo.status === "complete") beginCapture(tabId);
});

chrome.runtime.onInstalled.addListener(recoverInterruptedWork);
chrome.runtime.onStartup.addListener(recoverInterruptedWork);

chrome.alarms.onAlarm.addListener(async (alarm) => {
  if (alarm.name === "archive-next") {
    const state = await getState();
    if ((state.currentUrl || state.currentProjectUrl) && state.workerTabId) await beginCapture(state.workerTabId);
    else await processNext();
  }
});

chrome.tabs.onRemoved.addListener(async (tabId) => {
  const state = await getState();
  if (tabId !== state.workerTabId) return;
  const current = state.queue.find((item) => item.url === state.currentUrl) || state.projects.find((item) => item.url === state.currentProjectUrl);
  if (current?.status === "running") current.status = "pending";
  await setState({ queue: state.queue, projects: state.projects, workerTabId: null, currentUrl: null, currentProjectUrl: null });
  if (state.running) scheduleNext(10);
});

async function downloadText(filename, mime, body) {
  const dataUrl = `data:${mime};charset=utf-8,${encodeURIComponent(body)}`;
  await chrome.downloads.download({ url: dataUrl, filename: `chatgpt-archive/${filename}`, conflictAction: "uniquify", saveAs: false });
}

async function finishCapture(payload, sender) {
  const state = await getState();
  if (sender.tab?.id !== state.workerTabId || payload.url !== state.currentUrl) return;
  const item = state.queue.find((entry) => entry.url === state.currentUrl);
  if (!item) return;

  if (payload.blocked) {
    item.status = "error";
    item.error = payload.error;
    await setState({ queue: state.queue, running: false, currentUrl: null });
    return;
  }

  try {
    const formats = state.options.formats;
    if (formats.markdown) await downloadText(`${payload.slug}.md`, "text/markdown", payload.markdown);
    if (formats.html) await downloadText(`${payload.slug}.html`, "text/html", payload.html);
    if (formats.json) await downloadText(`${payload.slug}.json`, "application/json", JSON.stringify(payload.record, null, 2));
    item.status = "done";
    item.title = payload.record.title;
    item.messageCount = payload.record.messages.length;
    item.capturedAt = payload.record.capturedAt;
    delete item.error;
  } catch (error) {
    item.status = "error";
    item.error = `Download failed: ${error.message}`;
  }
  await setState({ queue: state.queue, currentUrl: null });
  if (state.running) {
    const base = Math.max(30, Number(state.options.delaySeconds) || 90);
    scheduleNext(base + Math.floor(Math.random() * base));
  }
}

async function finishProject(payload, sender) {
  const state = await getState();
  if (sender.tab?.id !== state.workerTabId || payload.url !== state.currentProjectUrl) return;
  const project = state.projects.find((entry) => entry.url === state.currentProjectUrl);
  if (!project) return;
  if (payload.blocked) {
    project.status = "error";
    project.error = payload.error;
    await setState({ projects: state.projects, running: false, currentProjectUrl: null });
    return;
  }
  const added = await addUrls(payload.urls);
  const refreshed = await getState();
  const refreshedProject = refreshed.projects.find((entry) => entry.url === payload.url);
  refreshedProject.status = "done";
  refreshedProject.title = payload.title;
  refreshedProject.chatCount = payload.urls.length;
  refreshedProject.newChatCount = added;
  refreshedProject.scannedAt = new Date().toISOString();
  delete refreshedProject.error;
  await setState({ projects: refreshed.projects, currentProjectUrl: null });
  scheduleNext(8);
}

chrome.runtime.onMessage.addListener((message, sender, respond) => {
  (async () => {
    if (message.type === "GET_STATE") return respond(await getState());
    if (message.type === "ADD_URLS") return respond({ added: await addUrls(message.urls) });
    if (message.type === "ADD_PROJECTS") return respond({ added: await addProjects(message.urls) });
    if (message.type === "PROJECTS_RESOLVED") return respond({ added: await addProjects(message.projects) });
    if (message.type === "PROJECT_RESOLUTION_ERROR") { console.error("Project resolution failed:", message.error); return respond({ ok: false }); }
    if (message.type === "SET_OPTIONS") { const state = await getState(); await setState({ options: { ...state.options, ...message.options } }); return respond({ ok: true }); }
    if (message.type === "START") { await setState({ running: true }); await processNext(); return respond({ ok: true }); }
    if (message.type === "PAUSE") { await setState({ running: false }); chrome.alarms.clear("archive-next"); return respond({ ok: true }); }
    if (message.type === "RETRY_ERRORS") {
      const state = await getState();
      for (const item of [...state.queue, ...state.projects]) if (item.status === "error") { item.status = "pending"; delete item.error; }
      await setState({ queue: state.queue, projects: state.projects }); return respond({ ok: true });
    }
    if (message.type === "CAPTURE_RESULT") { await finishCapture(message.payload, sender); return respond({ ok: true }); }
    if (message.type === "PROJECT_RESULT") { await finishProject(message.payload, sender); return respond({ ok: true }); }
    if (message.type === "EXPORT_INDEX") {
      const state = await getState();
      await downloadText("archive-index.json", "application/json", JSON.stringify({ exportedAt: new Date().toISOString(), projects: state.projects, queue: state.queue }, null, 2));
      return respond({ ok: true });
    }
  })().catch((error) => respond({ error: error.message }));
  return true;
});
