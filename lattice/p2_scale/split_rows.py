#!/usr/bin/env python3
"""Split data/rows_scale.txt into NCH per-worker row files.

Two modes.

  contiguous (default)  -- NCH contiguous blocks of roughly equal scan cost
      (the per-n cost grows like a power of n, so the cumulative n^3 is cut
      evenly).  All workers then finish at about the same time, which is the
      fastest route to a *complete* range, but it leaves NCH gaps if the run
      is interrupted.

  interleaved (-i)      -- worker j takes n = NLO + j, NLO + j + NCH, ...
      All workers advance through n in lockstep, so at every moment the union
      of what has been written is contiguous up to a common front.  This is
      the right mode for an open-ended run that may be stopped by the clock.

Splitting at all is what keeps twelve gp processes from each holding the whole
row cache (2.5 GB at n = 10000) in memory.

usage: split_rows.py NCH NHI [NLO] [-i]
"""
import sys, os

args = [a for a in sys.argv[1:] if not a.startswith("-")]
INTER = "-i" in sys.argv
NCH = int(args[0]) if len(args) > 0 else 12
NHI = int(args[1]) if len(args) > 1 else 3000
NLO = int(args[2]) if len(args) > 2 else 4
here = os.path.dirname(os.path.abspath(__file__))
fh = [open(os.path.join(here, "data", "rows_c%02d.txt" % j), "w") for j in range(NCH)]

if INTER:
    rng = [(NLO + j, NHI) for j in range(NCH)]
    which = lambda n: (n - NLO) % NCH
else:
    bnd = [NLO] + [int(round((NLO ** 3 + (NHI ** 3 - NLO ** 3) * j / NCH) ** (1 / 3.0)))
                   for j in range(1, NCH + 1)]
    bnd[-1] = NHI
    rng = [((bnd[j] + 1 if j else NLO), bnd[j + 1]) for j in range(NCH)]

    def which(n):
        for j, (lo, hi) in enumerate(rng):
            if lo <= n <= hi:
                return j
        return None

for line in open(os.path.join(here, "data", "rows_scale.txt")):
    n = int(line.split(" ", 1)[0])
    if n < NLO or n > NHI:
        continue
    j = which(n)
    if j is not None:
        fh[j].write(line)
for f in fh:
    f.close()
with open(os.path.join(here, "data", "chunks.txt"), "w") as f:
    for j, (lo, hi) in enumerate(rng):
        f.write("%02d %d %d\n" % (j, lo, hi))
print("mode = %s" % ("interleaved" if INTER else "contiguous"))
print("\n".join("worker %2d : n = %5d .. %5d%s"
                % (j, lo, hi, "  step %d" % NCH if INTER else "")
                for j, (lo, hi) in enumerate(rng)))
