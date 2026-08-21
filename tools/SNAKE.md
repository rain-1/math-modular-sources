# snake — dedicated compute box

- Login: `ssh snake` (key-based; config in `~/.ssh/config` on the dev box). Arch Linux, 8 cores, 15 GB RAM, RTX 2060 6 GB, PARI/GP 2.17.4, python3 + mpmath/sympy. No Sage/Lean.
- Working copy: `snake:~/math-modular-sources` (rsync'd from the dev box, excludes the chat archive, PDFs and `.git`). Sync with
  `rsync -az --exclude chatgpt-research-archive --exclude .git --exclude '*.pdf' ./ snake:~/math-modular-sources/`
  and pull results back with `rsync -az snake:~/math-modular-sources/lattice/ lattice/` (or a narrower path).
- Long jobs: run detached, one log per job:
  `ssh snake 'cd ~/math-modular-sources && nohup gp -q lattice/foo.gp > ~/jobs/logs/foo.log 2>&1 &'`
  Check: `ssh snake 'tail -n 5 ~/jobs/logs/foo.log; pgrep -a gp'`. Multi-hour runs are fine; up to 8 parallel gp processes.
- PARI memory: start scripts with `default(parisizemax, 8000000000)` or launch `gp -s 8G`.
- GPU: only useful for Python/CUDA experiments (batched LLL sweeps, etc.); PARI is CPU-only.
- Write results under the same relative paths as in the repo so rsync-back is trivial; stay inside `~/math-modular-sources` and `~/jobs`.
