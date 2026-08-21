# ChatGPT Research Archiver

This unpacked Chrome extension slowly saves conversations from your own logged-in ChatGPT account. It creates Markdown, HTML, and JSON files in `Downloads/chatgpt-archive/` and keeps a resumable queue in Chrome extension storage.

Version 0.4.1 supports ChatGPT's button-based Project sidebar and its paginated Project chat lists. It lists Project rows by name, resolves the five selected hidden addresses, and repeatedly uses **Load more conversations** until each selected Project is fully enumerated.

It does not request or store your password, cookies, or session token. It saves rendered conversation content, links, and basic source metadata. It does not call undocumented ChatGPT APIs.

## Install

1. In Chrome, open `chrome://extensions`.
2. Turn on **Developer mode**.
3. Select **Load unpacked**.
4. Choose this `chatgpt-archiver-extension` directory.
5. Pin **ChatGPT Research Archiver** to the toolbar.
6. Sign in at <https://chatgpt.com/> normally.

## Discover and archive

1. Open ChatGPT with the conversation sidebar visible.
2. Expand the Projects section near the top of the ChatGPT sidebar, then select **Find projects**. Project links and names are read immediately without scrolling. Directly visible chats outside Projects are ignored.
3. Select exactly five relevant Projects and choose **Queue selected projects**.
   - If no checkboxes appear, open the first relevant Project, select **Add current project**, then repeat for the other four Projects.
4. Alternatively, add the current conversation or paste conversation URLs.
5. Choose the output formats and delay. Start with **90–180 seconds**.
6. Select **Start / resume**. The extension first opens and scans the five selected Projects, then archives the conversations it found. It uses one inactive ChatGPT tab and downloads each completed archive automatically.
7. Select **Download index** periodically to save a copy of queue status and source URLs.

Chrome may ask whether to allow multiple downloads. Choose **Allow** for this extension. Keep Chrome open and the computer awake. The queue survives browser and service-worker restarts, but after reopening Chrome you may need to select **Start / resume** again.

## Safety and limitations

- The extension archives only the branch currently displayed for each conversation URL. It does not automatically select alternate assistant responses inside a chat.
- Project and sidebar discovery depend on ChatGPT's current page structure. If a Project is missed, open that Project and run discovery again; its visible chat links will be added directly. You can also paste conversation URLs.
- Attached file links and rendered images may appear in the archive, but the underlying files are not downloaded by this version.
- The extension pauses when it detects a login page, access challenge, or rate-limit warning. Resolve the issue in the archive tab, wait, press **Retry errors**, and resume with a longer delay.
- Do not publish the resulting archive without reviewing it for personal or confidential information.

## Update after code changes

Return to `chrome://extensions` and press the extension's reload button. Existing queue data remains in extension storage unless the extension is removed.

## Dedicated remote-debugging browser

Run `./launch-debug-chrome.sh` to open a separate Chrome profile with DevTools available only at `127.0.0.1:9222`. Sign into ChatGPT yourself in that window. The profile is stored outside this repository at `~/.local/share/chatgpt-archiver-debug-profile` by default; it contains an authenticated browser session and must be treated as private.

## Watch archive progress

Run `./archive-status.sh` in a terminal. It displays overall and per-Project progress, errors, the current conversation, and an estimated completion time, refreshing every ten seconds. Press Ctrl+C to close the status view; this does not stop the archive.
