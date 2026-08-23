#!/usr/bin/env python3
"""lattice/p2_structure/stats.py -- post-processing for consolidation/P2_STRUCTURE.md.

Reads data/struct_n120.csv and data/cfhist_n120.csv, checks the structural
propositions, fits the P2' regressions, and compares the population of
Catalan lattices with Haar-random unimodular planar lattices.
Usage:  python3 stats.py [datadir]
"""
import csv, math, sys, random, statistics as st
from collections import Counter

D = sys.argv[1] if len(sys.argv) > 1 else \
    "/home/ubuntu/code/math-modular-sources/lattice/p2_structure/data"
R2 = math.sqrt(2)

rows = [r for r in csv.DictReader(open(D + "/struct_n120.csv"))]
KS = ["22.4000", "23.0000", "23.9000"]

def f(r, c): return float(r[c])
def i(r, c): return int(r[c])

print("=" * 72); print("1. verification tallies (N = %d instances)" % len(rows))
print("   ok flag (integrality + form value + j-range stability): %d/%d"
      % (sum(i(r, 'ok') for r in rows), len(rows)))
print("   v1 is a CF convergent (cfkind=0): %d/%d"
      % (sum(1 for r in rows if r['cfkind'] == '0'), len(rows)))
c = Counter((i(r, 'cfidx') % 2, i(r, 'incone')) for r in rows)
print("   (parity of CF index, in-cone flag) joint counts:", dict(c))
print("   cone-vector kind (0=convergent,1=semiconvergent):",
      dict(Counter(r['cvkind'] for r in rows)))
print("   cone-vector index parity (must be even):",
      dict(Counter(i(r, 'cvidx') % 2 for r in rows)))
bad = [r for r in rows if i(r, 'incone') and (f(r, 'rho') > R2 + 1e-9 or abs(f(r, 'rho2') - 1) > 1e-9)]
print("   violations of  (v1 in cone) => rho2 = 1 and rho <= sqrt2 :", len(bad))
bad = [r for r in rows if not i(r, 'incone')
       and not (math.sqrt(3) / 2 * f(r, 'skew') - 1e-6 <= f(r, 'rho2') <= 1 + f(r, 'skew') + 1e-6)]
print("   violations of  (v1 not in cone) => (sqrt3/2) skew <= rho2 <= 1+skew :", len(bad))
q = [r for r in rows if f(r, 'rho') <= R2 and not i(r, 'incone')]
print("   rho <= sqrt2 with v1 outside the cone: %d, max skew there %.3f (bound 2sqrt2/sqrt3=%.3f)"
      % (len(q), max([f(r, 'skew') for r in q] or [0]), 2 * R2 / math.sqrt(3)))

print("=" * 72); print("2. P2' statistics by k   (n = 4..120, %d values each)" % (len(rows) // 3))
def linfit(xs, ys):
    n = len(xs); mx = sum(xs) / n; my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs); sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    b = sxy / sxx; a = my - b * mx
    res = [y - (a + b * x) for x, y in zip(xs, ys)]
    s2 = sum(t * t for t in res) / (n - 2); se = math.sqrt(s2 / sxx)
    return b, se, a, math.sqrt(s2)
for k in KS:
    rr = [r for r in rows if r['k'] == k]
    rho = [f(r, 'rho') for r in rr]; sk = [f(r, 'skew') for r in rr]
    ns = [f(r, 'n') for r in rr]
    b, se, a, sd = linfit(ns, [math.log(x) / n for x, n in zip(rho, ns)])
    b2, se2, a2, sd2 = linfit(ns, [math.log(x) for x in rho])
    print(" k=%s  in-cone %.3f  median rho %.4f  max rho %.1f  P(rho<=sqrt2) %.3f"
          % (k, sum(i(r, 'incone') for r in rr) / len(rr), st.median(rho), max(rho),
             sum(1 for x in rho if x <= R2) / len(rr)))
    print("        (1/n)log rho = (%+.5f +- %.5f) n %+.4f   [resid sd %.3f]" % (b, se, a, sd))
    print("        log rho      = (%+.5f +- %.5f) n %+.4f   [resid sd %.3f]" % (b2, se2, a2, sd2))
    print("        median skew %.3f  max skew %.1f  P(skew>20) %.3f"
          % (st.median(sk), max(sk), sum(1 for x in sk if x > 20) / len(rr)))
    blocks = [max(rr[j:j + 10], key=lambda r: -math.log(f(r, 'rho'))) for j in range(0, len(rr), 10)]
    print("        block minima of log rho (blocks of 10 consecutive n): "
          + " ".join("%.3f" % math.log(f(r, 'rho')) for r in blocks))

print("=" * 72); print("3. rho against the partial quotient at the balance index (parity odd)")
for k in KS:
    rr = [r for r in rows if r['k'] == k and not i(r, 'incone')]
    rat = [f(r, 'skew') / i(r, 'pqnext') for r in rr if i(r, 'pqnext') > 0]
    print(" k=%s  N=%d  skew/a_{i+1}: min %.3f med %.3f max %.3f ; corr(log skew, log a) = %.3f"
          % (k, len(rr), min(rat), st.median(rat), max(rat),
             corr := (lambda X, Y: (sum((x - st.mean(X)) * (y - st.mean(Y)) for x, y in zip(X, Y))
                                    / math.sqrt(sum((x - st.mean(X)) ** 2 for x in X)
                                                * sum((y - st.mean(Y)) ** 2 for y in Y))))(
                 [math.log(f(r, 'skew')) for r in rr], [math.log(i(r, 'pqnext')) for r in rr])))

print("=" * 72); print("4. Gauss-Kuzmin test on the partial quotients of h12/h11")
hr = [r for r in csv.DictReader(open(D + "/cfhist_n120.csv"))]
tot = [0] * 9; N = 0
for r in hr:
    for j in range(1, 9): tot[j - 1] += int(r["a%d" % j])
    tot[8] += int(r["a9plus"]); N += int(r["npq"])
print("   partial quotients counted: %d over %d instances" % (N, len(hr)))
gk = [math.log2(1 + 1 / (a * (a + 2))) for a in range(1, 9)]
gk.append(1 - sum(gk))
print("   a      :" + "".join("%8d" % a for a in range(1, 9)) + "     >=9")
print("   observed" + "".join("%8.4f" % (t / N) for t in tot))
print("   G-K    " + "".join("%8.4f" % g for g in gk))
chi = sum((tot[j] - N * gk[j]) ** 2 / (N * gk[j]) for j in range(9))
print("   chi^2 (8 df) = %.2f" % chi)

print("=" * 72); print("5. Haar-random unimodular planar lattices, same statistics")
def conedata(b1, b2):
    """exact cone minimum / first minimum for the basis (b1,b2), cone = first quadrant"""
    def gauss(u, v):
        if u[0] ** 2 + u[1] ** 2 > v[0] ** 2 + v[1] ** 2: u, v = v, u
        while True:
            mu = round((u[0] * v[0] + u[1] * v[1]) / (u[0] ** 2 + u[1] ** 2))
            v = (v[0] - mu * u[0], v[1] - mu * u[1])
            if v[0] ** 2 + v[1] ** 2 >= u[0] ** 2 + u[1] ** 2: return u, v
            u, v = v, u
    u, v = gauss(b1, b2)
    l1 = math.hypot(*u)
    best = None; bq = None
    for j in range(-4, 5):
        lo, hi = -1e18, 1e18; ok = True
        for r in range(2):
            aa = u[r]; bb = j * v[r]
            if aa == 0:
                if bb < 0: ok = False
            else:
                t = -bb / aa
                if aa > 0: lo = max(lo, t)
                else: hi = min(hi, t)
        if not ok: continue
        lo = math.ceil(lo); hi = math.floor(hi)
        for ii in set([lo, lo + 1, hi - 1, hi]):
            if ii < lo or ii > hi: continue
            x1 = ii * u[0] + j * v[0]; x2 = ii * u[1] + j * v[1]
            if x1 == 0 and x2 == 0: continue
            if x1 < -1e-12 or x2 < -1e-12: continue
            fv = x1 + x2
            if best is None or fv < best: best = fv
    incone = (u[0] >= 0 and u[1] >= 0) or (u[0] <= 0 and u[1] <= 0)
    return best / l1, math.hypot(*v) / l1, incone
random.seed(20260823)
M = 200000; rhos = []; sks = []; inc = 0
for _ in range(M):
    while True:
        x = random.uniform(-.5, .5); y = math.sqrt(3) / 2 / random.random()
        if x * x + y * y >= 1: break
    lam = 1 / math.sqrt(y); th = random.uniform(0, 2 * math.pi)
    ct, stn = math.cos(th), math.sin(th)
    b1 = (lam * ct, lam * stn); b2 = (lam * (x * ct - y * stn), lam * (x * stn + y * ct))
    rr, sk, ic = conedata(b1, b2)
    rhos.append(rr); sks.append(sk); inc += ic
print("   %d samples: in-cone %.4f  median rho %.4f  P(rho<=sqrt2) %.4f  P(rho>20) %.4f"
      % (M, inc / M, st.median(rhos), sum(1 for x in rhos if x <= R2) / M,
         sum(1 for x in rhos if x > 20) / M))
print("   median skew %.4f  P(skew>20) %.4f" % (st.median(sks), sum(1 for x in sks if x > 20) / M))
allr = [f(r, 'rho') for r in rows]; alls = [f(r, 'skew') for r in rows]
print("   Catalan (all %d): in-cone %.4f  median rho %.4f  P(rho<=sqrt2) %.4f  P(rho>20) %.4f"
      % (len(rows), sum(i(r, 'incone') for r in rows) / len(rows), st.median(allr),
         sum(1 for x in allr if x <= R2) / len(rows), sum(1 for x in allr if x > 20) / len(rows)))
print("   median skew %.4f  P(skew>20) %.4f" % (st.median(alls), sum(1 for x in alls if x > 20) / len(alls)))
srt = sorted(rhos)
def emp(x): 
    import bisect; return bisect.bisect_left(srt, x) / M
ks = max(abs(emp(x) - (j + 1) / len(allr)) for j, x in enumerate(sorted(allr)))
print("   Kolmogorov-Smirnov distance (rho, Catalan vs Haar) = %.4f ; 5%% level = %.4f"
      % (ks, 1.358 / math.sqrt(len(allr))))
