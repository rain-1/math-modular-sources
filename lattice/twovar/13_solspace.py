"""13: dimension of the joint power-series solution space (grid nullspace mod p,
several grid sizes to check stabilisation)."""
import sys, json, time
from libx import *
recs = json.load(open("recs.json"))
which = sys.argv[1:] if len(sys.argv) > 1 else list(recs)
p = BIGPRIMES[0]
for name in which:
    r = recs[name]
    reca = Rec(r["a"]["I"], r["a"]["mons"], r["a"]["coefs"], 0)
    recb = Rec(r["b"]["I"], r["b"]["mons"], r["b"]["coefs"], 1)
    line = "%-7s (order_a=%d, order_b=%d):  dim on grid " % (name, reca.I, recb.I)
    dims = []
    for N in (6, 8, 10, 12, 14, 16):
        t0 = time.time()
        ns, piv, free, cells, idx = grid_nullspace_modp(reca, recb, N, N, p)
        dims.append((N, len(ns)))
        line += " %d:%d" % (N, len(ns))
    print(line)
    ns, piv, free, cells, idx = grid_nullspace_modp(reca, recb, 12, 12, p)
    print("        free cells (N=12):", [cells[j] for j in free])
    sys.stdout.flush()
