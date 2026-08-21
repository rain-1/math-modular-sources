#!/usr/bin/env bash
set -euo pipefail

debug_profile="${XDG_DATA_HOME:-$HOME/.local/share}/chatgpt-archiver-debug-profile"
debug_port="${CHATGPT_ARCHIVER_DEBUG_PORT:-9222}"
extension_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$debug_profile"

exec /usr/bin/google-chrome \
  --user-data-dir="$debug_profile" \
  --remote-debugging-address=127.0.0.1 \
  --remote-debugging-port="$debug_port" \
  --disable-extensions-except="$extension_dir" \
  --load-extension="$extension_dir" \
  --no-first-run \
  --no-default-browser-check \
  --new-window \
  https://chatgpt.com/
