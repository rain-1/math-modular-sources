const $ = (id) => document.getElementById(id);
let discoveredProjects = [];

async function send(type, extra = {}) {
  return chrome.runtime.sendMessage({ type, ...extra });
}

async function activeTab() {
  const [tab] = await chrome.tabs.query({ active: true, currentWindow: true });
  return tab;
}

function notice(message, isError = false) {
  $("notice").textContent = message;
  $("notice").style.color = isError ? "#b42318" : "#287a3d";
}

async function refresh() {
  const state = await send("GET_STATE");
  const counts = { pending: 0, running: 0, done: 0, error: 0 };
  for (const item of state.queue || []) counts[item.status] = (counts[item.status] || 0) + 1;
  const projectCounts = { pending: 0, running: 0, done: 0, error: 0 };
  for (const item of state.projects || []) projectCounts[item.status] = (projectCounts[item.status] || 0) + 1;
  $("summary").textContent = `${state.running ? "Running" : "Paused"} · ${projectCounts.done}/${state.projects?.length || 0} projects scanned · ${counts.done} chats saved · ${counts.pending} pending · ${counts.error + projectCounts.error} errors`;
  $("markdown").checked = state.options?.formats?.markdown !== false;
  $("html").checked = state.options?.formats?.html !== false;
  $("json").checked = state.options?.formats?.json !== false;
  $("delay").value = String(state.options?.delaySeconds || 90);
  const queuedList = $("queued-project-list");
  queuedList.replaceChildren();
  for (const project of state.projects || []) {
    const item = document.createElement("li");
    item.textContent = project.title || project.url;
    queuedList.append(item);
  }
  $("queued-projects").hidden = !(state.projects || []).length;
}

async function saveOptions() {
  await send("SET_OPTIONS", {
    options: {
      formats: { markdown: $("markdown").checked, html: $("html").checked, json: $("json").checked },
      delaySeconds: Number($("delay").value)
    }
  });
}

$("add-current").addEventListener("click", async () => {
  const tab = await activeTab();
  const result = await send("ADD_PROJECTS", { urls: [{ url: tab?.url, title: tab?.title?.replace(/\s*[-–|]\s*ChatGPT\s*$/i, "") }] });
  notice(result.added ? "Added the current Project. Open the next Project and repeat." : "This page is not a Project, or it is already queued.", !result.added);
  refresh();
});

$("add-urls").addEventListener("click", async () => {
  const urls = $("urls").value.split(/\s+/).filter(Boolean);
  const result = await send("ADD_URLS", { urls });
  $("urls").value = "";
  notice(`Added ${result.added} conversation${result.added === 1 ? "" : "s"}.`);
  refresh();
});

$("discover").addEventListener("click", async () => {
  const tab = await activeTab();
  if (!tab?.url?.startsWith("https://chatgpt.com/")) return notice("Open chatgpt.com first.", true);
  notice("Reading the visible Project list…");
  try {
    const result = await chrome.tabs.sendMessage(tab.id, { type: "DISCOVER_PROJECTS" });
    discoveredProjects = result.projects;
    const list = $("project-list");
    list.replaceChildren();
    for (const [index, project] of discoveredProjects.entries()) {
      const label = document.createElement("label");
      const input = document.createElement("input");
      input.type = "checkbox";
      input.value = String(index);
      const text = document.createElement("span");
      text.textContent = project.title || project.url;
      label.append(input, text);
      list.append(label);
      input.addEventListener("change", () => {
        const count = document.querySelectorAll('#project-list input:checked').length;
        $("selection-count").textContent = `Selected ${count} of 5`;
      });
    }
    $("project-picker").hidden = false;
    $("selection-count").textContent = "Selected 0 of 5";
    notice(discoveredProjects.length ? `Found ${discoveredProjects.length} visible Projects. Tick five below.` : "No Project links were detected. Expand the Projects section in ChatGPT and try again.", !discoveredProjects.length);
  } catch (error) {
    notice(`Discovery failed: ${error.message}`, true);
  }
});

$("add-projects").addEventListener("click", async () => {
  const selected = [...document.querySelectorAll('#project-list input:checked')].map((input) => discoveredProjects[Number(input.value)]);
  if (selected.length !== 5) return notice(`Select exactly five Projects; currently selected: ${selected.length}.`, true);
  const buttonBased = selected.filter((project) => project.buttonBased || !project.url);
  const linked = selected.filter((project) => project.url);
  const result = linked.length ? await send("ADD_PROJECTS", { urls: linked }) : { added: 0 };
  if (buttonBased.length) {
    const tab = await activeTab();
    await chrome.tabs.sendMessage(tab.id, { type: "RESOLVE_PROJECTS", titles: buttonBased.map((project) => project.title) });
  }
  $("project-picker").hidden = true;
  notice(buttonBased.length ? "Resolving the five Project addresses now; the page will briefly open each Project." : `Queued ${result.added} new Projects.`);
  refresh();
});

$("start").addEventListener("click", async () => { await saveOptions(); await send("START"); notice("Archiving started."); refresh(); });
$("pause").addEventListener("click", async () => { await send("PAUSE"); notice("Paused after the current operation."); refresh(); });
$("retry").addEventListener("click", async () => { await send("RETRY_ERRORS"); notice("Errors returned to the queue."); refresh(); });
$("export-index").addEventListener("click", async () => { await send("EXPORT_INDEX"); notice("Index downloaded."); });
for (const id of ["markdown", "html", "json", "delay"]) $(id).addEventListener("change", saveOptions);

refresh().catch((error) => notice(error.message, true));
