#!/usr/bin/env python3
"""lattice/p2_scale/stats3.py -- the joint question of P2_SCALE.md sec.10:
is the parity independent of the size bias of a_{i(n)+1}?  Plus the drift of
the even-density with n, and the run-length distribution."""
import csv, glob, json, math, os, sys, collections
import numpy as np
from scipy import stats as ss
HERE = os.path.dirname(os.path.abspath(__file__))
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10**9
rows = [r for f in sorted(glob.glob(os.path.join(HERE, "data", "scan*_c*.csv")))
        for r in csv.DictReader(open(f)) if r.get("n") and int(r["n"]) <= NMAX]
KS = sorted({float(r["k"]) for r in rows})
D = {k: sorted([r for r in rows if float(r["k"]) == k], key=lambda r: int(r["n"])) for k in KS}
def A(k, c, cast=float): return np.array([cast(r[c]) for r in D[k]])
OUT = {}
print("=" * 78)
print("F. IS THE PARITY INDEPENDENT OF THE SIZE OF a_{i(n)+1} ?")
OUT["parity_vs_pq"] = {}
for k in KS:
    par = A(k, "parity", int); aq = A(k, "pqnext"); la = np.log(np.maximum(aq, 1))
    r, p = ss.pointbiserialr(par, la)
    m0, m1 = la[par == 0].mean(), la[par == 1].mean()
    u = ss.mannwhitneyu(la[par == 0], la[par == 1])
    print("  k=%-7.4f  mean log a_{i+1}: even %.4f, odd %.4f ; point-biserial "
          "r=%+.4f (p=%.3f) ; Mann-Whitney p=%.3f" % (k, m0, m1, r, p, u.pvalue))
    OUT["parity_vs_pq"][str(k)] = dict(mean_even=m0, mean_odd=m1, r=r, p=p, mw_p=u.pvalue)

print("\n" + "=" * 78)
print("G. DRIFT OF THE EVEN-DENSITY WITH n")
OUT["even_drift"] = {}
BINS = [(4, 1000), (1001, 2000), (2001, 3000), (3001, 4000), (4001, 5000),
        (5001, 6000), (6001, 7000), (7001, 8000), (8001, 9000), (9001, 10000)]
for k in KS:
    nn = A(k, "n", int); par = A(k, "parity", int); out = {}
    line = []
    for lo, hi in BINS:
        s = (nn >= lo) & (nn <= hi)
        if s.sum() < 50: continue
        f = float((par[s] == 0).mean()); out["%d-%d" % (lo, hi)] = f
        line.append("%d-%d: %.3f" % (lo, hi, f))
    b = np.polyfit(nn, (par == 0).astype(float), 1)
    print("  k=%-7.4f  " % k + "  ".join(line))
    print("            linear trend in n: %+.2e per unit" % b[0])
    OUT["even_drift"][str(k)] = dict(bins=out, slope=float(b[0]))

print("\n" + "=" * 78)
print("H. THE FULL RUN-LENGTH SPECTRUM OF ODD PARITY")
OUT["run_spectrum"] = {}
for k in KS:
    nn = A(k, "n", int); par = A(k, "parity", int)
    runs = []; cur = 0; prev = None
    for n, x in zip(nn, par):
        if x == 1 and cur > 0 and prev is not None and n == prev + 1: cur += 1
        elif x == 1: cur = 1
        else:
            if cur: runs.append(cur)
            cur = 0
        prev = n
    if cur: runs.append(cur)
    c = collections.Counter(runs)
    exp = {L: (len(nn) + 1) * 2.0 ** (-(L + 2)) for L in sorted(c)}
    print("  k=%-7.4f  runs of length L (observed / fair-coin expectation):" % k)
    print("            " + "  ".join("%d:%d/%.1f" % (L, c[L], exp[L]) for L in sorted(c)))
    OUT["run_spectrum"][str(k)] = dict(observed={str(L): c[L] for L in sorted(c)},
                                       expected={str(L): exp[L] for L in sorted(c)})
json.dump(OUT, open(os.path.join(HERE, "data", "summary3.json"), "w"), indent=1)
print("\nwritten data/summary3.json")
