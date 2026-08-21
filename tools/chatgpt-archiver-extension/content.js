const wait = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

function conversationUrls(root = document) {
  return [...root.querySelectorAll('a[href*="/c/"]')]
    .map((anchor) => {
      try { return new URL(anchor.getAttribute("href"), location.origin).href; } catch { return null; }
    })
    .filter(Boolean);
}

function projectEntries(root = document) {
  const linked = [...root.querySelectorAll('a[href*="g-p-"], a[href*="/project/"], a[href*="/projects/"]')]
    .map((anchor) => {
      try {
        const url = new URL(anchor.getAttribute("href"), location.origin);
        const isProject = /\/g\/g-p-[^/]+(?:\/project)?\/?$/i.test(url.pathname) || /\/projects?\/[^/]+\/?$/i.test(url.pathname);
        return isProject ? { url: url.href, title: anchor.innerText.trim() || anchor.getAttribute("aria-label") || "Untitled Project" } : null;
      } catch { return null; }
    })
    .filter(Boolean);
  const buttonBased = [...root.querySelectorAll('button[aria-label="Open project home"]')].map((button, index) => ({
    url: null,
    title: button.closest("li")?.innerText.trim() || button.parentElement?.parentElement?.innerText.trim() || `Project ${index + 1}`,
    buttonBased: true
  }));
  const unique = new Map([...linked, ...buttonBased].map((project) => [project.url || `title:${project.title}`, project]));
  return [...unique.values()];
}

async function resolveButtonProjects(titles) {
  const startUrl = location.href;
  const resolved = [];
  for (const title of titles) {
    const row = [...document.querySelectorAll('button[aria-label="Open project home"]')]
      .map((button) => ({ button, title: button.closest("li")?.innerText.trim() || button.parentElement?.parentElement?.innerText.trim() || "" }))
      .find((entry) => entry.title === title);
    if (!row) continue;
    const before = location.href;
    row.button.click();
    for (let attempt = 0; attempt < 40 && location.href === before; attempt++) await wait(250);
    if (location.href !== before && /g-p-/i.test(location.pathname)) resolved.push({ url: location.href, title });
    history.back();
    for (let attempt = 0; attempt < 40 && location.href !== startUrl; attempt++) await wait(250);
    await wait(400);
  }
  await chrome.runtime.sendMessage({ type: "PROJECTS_RESOLVED", projects: resolved });
}

function sidebarScroller() {
  const links = [...document.querySelectorAll('nav a[href*="/c/"], aside a[href*="/c/"]')];
  for (const link of links) {
    let node = link.parentElement;
    while (node && node !== document.body) {
      const style = getComputedStyle(node);
      if (/(auto|scroll)/.test(style.overflowY) && node.scrollHeight > node.clientHeight * 1.2) return node;
      node = node.parentElement;
    }
  }
  return document.querySelector("nav") || document.querySelector("aside");
}

async function discoverConversations() {
  const found = new Set(conversationUrls());
  const projects = new Map(projectEntries().map((project) => [project.url, project]));
  if (/\/g\/g-p-[^/]+\/project\/?$/i.test(location.pathname) || /\/projects?\/[^/]+\/?$/i.test(location.pathname)) projects.set(location.href, { url: location.href, title: document.title.replace(/\s*[-–|]\s*ChatGPT\s*$/i, "").trim() || "Current Project" });
  const scroller = sidebarScroller();
  if (!scroller) return { urls: [...found], projects: [...projects.values()] };
  const original = scroller.scrollTop;
  let unchanged = 0;
  let previousSize = found.size;
  for (let i = 0; i < 500 && unchanged < 8; i++) {
    for (const url of conversationUrls(scroller)) found.add(url);
    for (const project of projectEntries(scroller)) projects.set(project.url, project);
    scroller.scrollTop = Math.min(scroller.scrollHeight, scroller.scrollTop + Math.max(250, scroller.clientHeight * .8));
    await wait(450);
    unchanged = found.size === previousSize && scroller.scrollTop + scroller.clientHeight >= scroller.scrollHeight - 5 ? unchanged + 1 : 0;
    previousSize = found.size;
  }
  scroller.scrollTop = original;
  return { urls: [...found], projects: [...projects.values()] };
}

async function discoverProjectChats(expectedUrl) {
  await wait(2500);
  const blocked = blockedReason();
  if (blocked) return { url: expectedUrl, blocked: true, error: blocked };
  const found = new Set(conversationUrls());
  for (let page = 0; page < 200; page++) {
    for (const url of conversationUrls(document.querySelector("main") || document)) found.add(url);
    const loadMore = [...document.querySelectorAll("main button")].find((button) => button.innerText.trim() === "Load more conversations");
    if (!loadMore || loadMore.disabled) break;
    const previousSize = found.size;
    loadMore.click();
    for (let attempt = 0; attempt < 40; attempt++) {
      await wait(250);
      for (const url of conversationUrls(document.querySelector("main") || document)) found.add(url);
      if (found.size > previousSize) break;
    }
  }
  const title = document.title.replace(/\s*[-–|]\s*ChatGPT\s*$/i, "").trim() || "ChatGPT Project";
  return { url: expectedUrl, title, urls: [...found] };
}

function findConversationScroller() {
  const turn = document.querySelector('[data-testid^="conversation-turn-"]');
  let node = turn?.parentElement;
  while (node && node !== document.body) {
    const style = getComputedStyle(node);
    if (/(auto|scroll)/.test(style.overflowY) && node.scrollHeight > node.clientHeight) return node;
    node = node.parentElement;
  }
  return document.scrollingElement;
}

async function loadAllTurns() {
  let stable = 0;
  let last = "";
  for (let i = 0; i < 80 && stable < 5; i++) {
    const scroller = findConversationScroller();
    if (!scroller) break;
    scroller.scrollTop = 0;
    await wait(700);
    const signature = `${document.querySelectorAll('[data-testid^="conversation-turn-"]').length}:${scroller.scrollHeight}`;
    stable = signature === last ? stable + 1 : 0;
    last = signature;
  }
}

function sanitizeFragment(element) {
  const clone = element.cloneNode(true);
  clone.querySelectorAll("script,style,button,form,textarea,input,svg").forEach((node) => node.remove());
  for (const node of [clone, ...clone.querySelectorAll("*")]) {
    for (const attr of [...node.attributes]) {
      if (attr.name.startsWith("on") || ["nonce", "integrity", "srcset"].includes(attr.name)) node.removeAttribute(attr.name);
    }
  }
  return clone.innerHTML.trim();
}

function roleFor(article, index) {
  const testId = article.getAttribute("data-testid") || "";
  if (testId.includes("user")) return "user";
  if (testId.includes("assistant")) return "assistant";
  if (article.querySelector('[data-message-author-role="user"]')) return "user";
  if (article.querySelector('[data-message-author-role="assistant"]')) return "assistant";
  const label = (article.getAttribute("aria-label") || "").toLowerCase();
  if (label.includes("you") || label.includes("user")) return "user";
  if (label.includes("assistant") || label.includes("chatgpt")) return "assistant";
  return index % 2 === 0 ? "user" : "assistant";
}

function markdownEscape(text) {
  return text.replace(/\r/g, "").replace(/\n{3,}/g, "\n\n").trim();
}

function safeSlug(title, id) {
  const base = title.normalize("NFKD").replace(/[^a-zA-Z0-9]+/g, "-").replace(/^-|-$/g, "").slice(0, 70).toLowerCase() || "conversation";
  return `${base}-${id.slice(0, 8)}`;
}

function blockedReason() {
  const text = document.body.innerText.toLowerCase();
  if (/too many requests|rate limit|unusual activity/.test(text)) return "ChatGPT displayed a rate-limit or unusual-activity message. Archiving was paused.";
  if (/checking your browser|verify you are human|cloudflare/.test(text)) return "An access challenge was detected. Archiving was paused.";
  if (/log in|sign up/.test(text) && !document.querySelector('[data-testid^="conversation-turn-"]')) return "The archive tab is not signed in to ChatGPT. Archiving was paused.";
  return null;
}

async function capture(expectedUrl) {
  await wait(2500);
  const blocked = blockedReason();
  if (blocked) return { url: expectedUrl, blocked: true, error: blocked };
  await loadAllTurns();
  const articles = [...document.querySelectorAll('[data-testid^="conversation-turn-"]')];
  if (!articles.length) return { url: expectedUrl, blocked: true, error: "No conversation messages were found. The page may still be loading or its layout may have changed." };

  const messages = articles.map((article, index) => ({
    index: index + 1,
    role: roleFor(article, index),
    text: markdownEscape(article.innerText),
    html: sanitizeFragment(article)
  }));
  const id = new URL(expectedUrl).pathname.split("/").filter(Boolean).pop();
  const title = document.title.replace(/\s*[-–|]\s*ChatGPT\s*$/i, "").trim() || `Conversation ${id}`;
  const capturedAt = new Date().toISOString();
  const record = { schemaVersion: 1, id, title, url: expectedUrl, capturedAt, note: "Captures the branch visible in ChatGPT at archive time.", messages };
  const slug = safeSlug(title, id);
  const markdown = [`# ${title}`, "", `- Source: ${expectedUrl}`, `- Captured: ${capturedAt}`, `- Conversation ID: ${id}`, "", ...messages.flatMap((m) => [`## ${m.role === "user" ? "User" : "Assistant"} ${m.index}`, "", m.text, ""])].join("\n");
  const sections = messages.map((m) => `<section class="message ${m.role}"><h2>${m.role === "user" ? "User" : "Assistant"} ${m.index}</h2><div>${m.html}</div></section>`).join("\n");
  const html = `<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>${title.replace(/[<>&]/g, "")}</title><style>body{font:16px/1.55 system-ui,sans-serif;max-width:900px;margin:auto;padding:2rem;color:#202124}.meta{color:#666}.message{padding:1rem 1.25rem;margin:1rem 0;border-radius:12px;background:#f5f5f5}.user{background:#e8f0fe}pre{overflow:auto;background:#222;color:#eee;padding:1rem}img{max-width:100%}table{border-collapse:collapse}td,th{border:1px solid #aaa;padding:.4rem}</style></head><body><h1>${title.replace(/[<>&]/g, "")}</h1><p class="meta">Source: <a href="${expectedUrl}">${expectedUrl}</a><br>Captured: ${capturedAt}</p>${sections}</body></html>`;
  return { url: expectedUrl, slug, record, markdown, html };
}

chrome.runtime.onMessage.addListener((message, _sender, respond) => {
  if (message.type === "DISCOVER_PROJECTS") {
    const unique = new Map(projectEntries(document).map((project) => [project.url, project]));
    if (/\/g\/g-p-[^/]+(?:\/project)?\/?$/i.test(location.pathname) || /\/projects?\/[^/]+\/?$/i.test(location.pathname)) {
      unique.set(location.href, { url: location.href, title: document.title.replace(/\s*[-–|]\s*ChatGPT\s*$/i, "").trim() || "Current Project" });
    }
    respond({ projects: [...unique.values()] });
    return;
  }
  if (message.type === "RESOLVE_PROJECTS") {
    resolveButtonProjects(message.titles).catch((error) => chrome.runtime.sendMessage({ type: "PROJECT_RESOLUTION_ERROR", error: error.message }));
    respond({ started: true });
    return;
  }
  if (message.type === "DISCOVER_CONVERSATIONS") {
    discoverConversations().then(respond).catch((error) => respond({ urls: [], projects: [], error: error.message }));
    return true;
  }
  if (message.type === "CAPTURE_PROJECT") {
    discoverProjectChats(message.expectedUrl).then((payload) => chrome.runtime.sendMessage({ type: "PROJECT_RESULT", payload }));
    respond({ started: true });
  }
  if (message.type === "CAPTURE_CONVERSATION") {
    capture(message.expectedUrl).then((payload) => chrome.runtime.sendMessage({ type: "CAPTURE_RESULT", payload }));
    respond({ started: true });
  }
});
