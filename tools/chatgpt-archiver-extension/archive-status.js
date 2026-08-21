#!/usr/bin/env node
"use strict";

const fs = require("fs");
const path = require("path");

const output = process.env.CHATGPT_ARCHIVE_OUTPUT || path.join(process.env.HOME, ".local/share/chatgpt-research-archive");
const stateFile = path.join(output, "archive-state.json");
const watch = process.argv.includes("--watch") || process.argv.includes("-w");

function duration(seconds) {
  if (!Number.isFinite(seconds) || seconds < 0) return "unknown";
  const hours = Math.floor(seconds / 3600);
  const minutes = Math.floor((seconds % 3600) / 60);
  if (hours) return `${hours}h ${minutes}m`;
  return `${Math.max(1, minutes)}m`;
}

function render() {
  let state;
  try { state = JSON.parse(fs.readFileSync(stateFile, "utf8")); }
  catch (error) { return `Archive status unavailable: ${error.message}\nExpected: ${stateFile}`; }

  const conversations = state.conversations || [];
  const count = (status) => conversations.filter((item) => item.status === status).length;
  const done = count("done");
  const pending = count("pending");
  const running = count("running");
  const errors = count("error");
  const total = conversations.length;
  const percent = total ? (done / total * 100).toFixed(1) : "0.0";
  const completed = conversations.filter((item) => item.status === "done" && item.capturedAt);
  const first = completed.map((item) => Date.parse(item.capturedAt)).filter(Number.isFinite).sort()[0];
  const elapsedSeconds = first ? (Date.now() - first) / 1000 : NaN;
  const rate = done > 1 && elapsedSeconds > 0 ? (done - 1) / elapsedSeconds : 1 / 135;
  const remainingSeconds = rate > 0 ? (total - done) / rate : NaN;
  const eta = Number.isFinite(remainingSeconds) ? new Date(Date.now() + remainingSeconds * 1000).toLocaleString() : "unknown";
  const current = conversations.find((item) => item.status === "running") || conversations.find((item) => item.status === "pending");
  const projects = [...new Set(conversations.map((item) => item.project))];

  const lines = [
    "ChatGPT Research Archive",
    "=".repeat(48),
    `Progress:  ${done}/${total} (${percent}%)`,
    `Status:    ${running} running · ${pending} pending · ${errors} errors`,
    `Remaining: about ${duration(remainingSeconds)} · ETA ${eta}`,
    `Updated:   ${state.updatedAt || "unknown"}`,
    `Output:    ${output}`,
    "",
    "Projects"
  ];
  for (const project of projects) {
    const items = conversations.filter((item) => item.project === project);
    lines.push(`  ${project.padEnd(22)} ${String(items.filter((item) => item.status === "done").length).padStart(2)}/${String(items.length).padEnd(2)} saved`);
  }
  if (current) lines.push("", `${running ? "Currently saving" : "Next"}: [${current.project}] ${current.title}`);
  const failed = conversations.filter((item) => item.status === "error");
  if (failed.length) {
    lines.push("", "Errors");
    for (const item of failed.slice(0, 8)) lines.push(`  [${item.project}] ${item.title}: ${item.error}`);
  }
  if (watch) lines.push("", "Refreshing every 10 seconds. Press Ctrl+C to close this status view.");
  return lines.join("\n");
}

function display() {
  if (watch) process.stdout.write("\x1b[2J\x1b[H");
  process.stdout.write(`${render()}\n`);
}

display();
if (watch) setInterval(display, 10000);
