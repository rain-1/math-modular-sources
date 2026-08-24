#!/usr/bin/env python3
"""lattice/p2_scale/stats2.py -- the second-order questions of P2_SCALE.md:
robustness of the headline slope, the liminf statistic that P2' actually
asks about, the corrected midpoint law, the exact (24-k) n law, and the
archimedean rate."""
import csv, glob, json, math, os, sys
import numpy as np
from scipy import stats as ss
HERE = os.path.dirname(os.path.abspath(__file__))
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10**9
rows = []
for f in sorted(glob.glob(os.path.join(HERE, "data", "scan*_c*.csv"))):
    for r in csv.DictReader(open(f)):
        if r.get("n") and int(r["n"]) <= NMAX: rows.append(r)
KS = sorted({float(r["k"]) for r in rows})
D = {k: sorted([r for r in rows if float(r["k"]) == k], key=lambda r: int(r["n"])) for k in KS}
S2 = math.sqrt(2.0); OUT = {}
def A(k, c, cast=float): return np.array([cast(r[c]) for r in D[k]])
def ols(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float); nn = len(x)
    xb = x.mean(); sxx = ((x - xb) ** 2).sum()
    b = ((x - xb) * (y - y.mean())).sum() / sxx
    a = y.mean() - b * xb; res = y - (a + b * x)
    return b, a, math.sqrt((res ** 2).sum() / (nn - 2) / sxx), res

print("=" * 78)
print("A. ROBUSTNESS OF THE log rho SLOPE  (%d instances, n <= %d)"
      % (len(rows), max(int(r["n"]) for r in rows)))
OUT["slope_robust"] = {}
for k in KS:
    nn = A(k, "n", int); lr = np.log(A(k, "rho"))
    b, a, se, res = ols(nn, lr)
    ts = ss.theilslopes(lr, nn, 0.95)
    # winsorised at the 99th percentile, and the trend of P[rho<=sqrt2]
    w = np.minimum(lr, np.quantile(lr, 0.99))
    bw, aw, sew, _ = ols(nn, w)
    ind = (A(k, "rho") <= S2).astype(float)
    bi, ai, sei, _ = ols(nn, ind)
    half = len(nn) // 2
    m1, m2 = lr[:half].mean(), lr[half:].mean()
    s1 = lr[:half].std(ddof=1) / math.sqrt(half); s2 = lr[half:].std(ddof=1) / math.sqrt(len(lr) - half)
    print("  k=%-7.4f  OLS %+.2e +- %.2e | Theil-Sen %+.2e [%+.2e,%+.2e] | "
          "winsorised(99%%) %+.2e +- %.2e" % (k, b, se, ts[0], ts[2], ts[3], bw, sew))
    print("            mean log rho: first half %.4f +- %.4f, second half %.4f "
          "+- %.4f, difference %+.4f +- %.4f  (over %d in n)"
          % (m1, s1, m2, s2, m2 - m1, math.hypot(s1, s2), nn[half] - nn[0]))
    print("            P[rho<=sqrt2] trend per unit n: %+.2e +- %.2e" % (bi, sei))
    OUT["slope_robust"][str(k)] = dict(ols=b, ols_se=se, theilsen=ts[0],
        ts_lo=ts[2], ts_hi=ts[3], winsor=bw, winsor_se=sew,
        half1=m1, half2=m2, diff=m2 - m1, diff_se=math.hypot(s1, s2),
        indicator_slope=bi, indicator_se=sei)

print("\n" + "=" * 78)
print("B. WHAT P2' ACTUALLY ASKS: the liminf statistic")
print("   rho >= 1 always, so (1/n) log rho >= 0; P2' asks for infinitely many")
print("   n with (1/n) log rho -> 0.  Two exact finite-n statements:")
OUT["liminf"] = {}
for k in KS:
    nn = A(k, "n", int); rho = A(k, "rho"); lr = np.log(rho)
    worst = {}
    for W in (6, 8, 12, 16, 25, 50):
        # max over windows of W consecutive n of the window minimum of log rho
        mm = np.array([lr[i:i + W].min() for i in range(len(lr) - W + 1)])
        worst[W] = float(mm.max())
    dec = {}
    for lo, hi in ((4, 100), (101, 300), (301, 1000), (1001, 2000), (2001, 3000),
                   (3001, 4500), (4501, 6000), (6001, 8000), (8001, 10000)):
        sel = (nn >= lo) & (nn <= hi)
        if sel.sum() == 0: continue
        dec["%d-%d" % (lo, hi)] = float((lr[sel] / nn[sel]).max())
    print("  k=%-7.4f  max over windows of W consecutive n of (min log rho): "
          % k + "  ".join("W=%d: %.4f" % (W, v) for W, v in worst.items()))
    print("            max (1/n) log rho by range: "
          + "  ".join("%s: %.5f" % (r, v) for r, v in dec.items()))
    OUT["liminf"][str(k)] = dict(window_worst=worst, max_lograte=dec)

print("\n" + "=" * 78)
print("C. THE MIDPOINT LAW, CORRECTED FOR THE ABSORPTION CONTENT m_Z")
print("   i/L = 1/2 + (log m_Z - log h22 - (E2-E1)n) / (2 (log h11 - log m_Z))")
E21 = 1.2935564764
OUT["midpoint"] = {}
for k in KS:
    nn = A(k, "n", int); ir = A(k, "irat")
    lh11 = A(k, "logh11"); lh22 = A(k, "logh22"); lmz = A(k, "logmz")
    pred_old = 0.5 + E21 / (2 * lh11)
    pred_new = 0.5 + (lmz - lh22 - E21) / (2 * (lh11 - lmz))
    late = nn > 1000
    print("  k=%-7.4f  observed %.5f +- %.5f (s.e. %.5f) | old prediction "
          "%.5f | corrected prediction %.5f"
          % (k, ir[late].mean(), ir[late].std(ddof=1),
             ir[late].std(ddof=1) / math.sqrt(late.sum()),
             pred_old[late].mean(), pred_new[late].mean()))
    OUT["midpoint"][str(k)] = dict(obs=float(ir[late].mean()),
        sd=float(ir[late].std(ddof=1)),
        se=float(ir[late].std(ddof=1) / math.sqrt(late.sum())),
        pred_old=float(pred_old[late].mean()), pred_new=float(pred_new[late].mean()),
        logh11=float(lh11[late].mean()), logmz=float(lmz[late].mean()))

print("\n" + "=" * 78)
print("D. THE EXACT 2-ADIC LAW  v_2(q xi_2 - p) = (24 - k) n + O(1)")
OUT["v2law"] = {}
for k in KS:
    nn = A(k, "n", int); v2 = A(k, "v2qp")
    b, a, se, res = ols(nn, v2)
    r0 = v2 - (24 - k) * nn
    print("  k=%-7.4f  fitted slope %.6f (s.e. %.6f)  vs 24-k = %.6f  "
          "[difference %+.1e]" % (k, b, se, 24 - k, b - (24 - k)))
    print("            v2 - (24-k)n : mean %+.3f  sd %.3f  min %+.0f  max %+.0f"
          % (r0.mean(), r0.std(ddof=1), r0.min(), r0.max()))
    OUT["v2law"][str(k)] = dict(slope=b, se=se, pred=24 - k, diff=b - (24 - k),
        resid_mean=float(r0.mean()), resid_sd=float(r0.std(ddof=1)),
        resid_min=float(r0.min()), resid_max=float(r0.max()))

print("\n" + "=" * 78)
print("E. THE ARCHIMEDEAN RATE  (1/n) log |q_n G - p_n|  against F(k)")
kst = 22.3512905953; L2 = math.log(2)
OUT["Frate"] = {}
for k in KS:
    nn = A(k, "n", int); lc = A(k, "logcone"); l1 = A(k, "loglam1")
    late = nn > 4000 if (nn > 4000).sum() > 30 else nn > nn.max() // 2
    F = L2 / 2 * (kst - k)
    print("  k=%-7.4f  measured %.5f (n>%d)   F(k) = %.5f   half-log-covol "
          "%.5f   log q / n = %.4f"
          % (k, lc[late].mean(), int(nn[late].min()) - 1, F, l1[late].mean(),
             A(k, "logq")[late].mean()))
    OUT["Frate"][str(k)] = dict(measured=float(lc[late].mean()), F=F,
        loglam1=float(l1[late].mean()), logq=float(A(k, "logq")[late].mean()))

json.dump(OUT, open(os.path.join(HERE, "data", "summary2.json"), "w"), indent=1)
print("\nwritten data/summary2.json")
