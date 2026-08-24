#!/usr/bin/env python3
"""lattice/p2_scale/stats.py -- all statistics of consolidation/P2_SCALE.md.

Reads data/scan_c*.csv (one exact line per (k,n)) and writes
data/summary.json plus a human-readable report on stdout.

Sections
  1  parity of the balance index i(n): frequency, runs, autocorrelation,
     n mod m (m <= 60, Bonferroni), cross-k agreement, parity <-> cone
  2  log(cone-min / lambda_1) against n: the headline slope
  3  Gauss-Kuzmin at scale, and the partial quotient AT the balance index
  4  drift of i(n)/L(n) against 1/2 + (E2-E1)/(2 kappa) = 0.5238
  5  v_2(q_n xi_2 - p_n) against the (24.02-k) n law
"""
import csv, glob, json, math, os, sys
import numpy as np
from scipy import stats as ss

HERE = os.path.dirname(os.path.abspath(__file__))
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10**9

rows = []
for f in sorted(glob.glob(os.path.join(HERE, "data", "scan*_c*.csv"))):
    for r in csv.DictReader(open(f)):
        if r.get("n") and int(r["n"]) <= NMAX:
            rows.append(r)
rows.sort(key=lambda r: (float(r["k"]), int(r["n"])))
KS = sorted({float(r["k"]) for r in rows})
D = {k: [r for r in rows if float(r["k"]) == k] for k in KS}
OUT = {"n_instances": len(rows), "k_values": KS}
S2 = math.sqrt(2.0)

def arr(k, col, cast=float):
    return np.array([cast(r[col]) for r in D[k]])

def longest_run(nn, b):
    """longest run of True over CONSECUTIVE integers n (gaps break a run),
    and the n at which the record run starts"""
    best = cur = 0; at = start = None; prev = None
    for n, x in zip(nn, b):
        if x and cur > 0 and prev is not None and n == prev + 1:
            cur += 1
        elif x:
            cur = 1; start = n
        else:
            cur = 0; start = None
        if cur > best: best, at = cur, start
        prev = n
    return best, at

def runs_test(b):
    b = np.asarray(b, dtype=bool); nn = len(b)
    n1 = int(b.sum()); n0 = nn - n1
    if n1 == 0 or n0 == 0: return float("nan"), float("nan"), 0
    runs = 1 + int((b[1:] != b[:-1]).sum())
    mu = 2.0 * n0 * n1 / nn + 1
    var = 2.0 * n0 * n1 * (2.0 * n0 * n1 - nn) / (nn * nn * (nn - 1))
    z = (runs - mu) / math.sqrt(var)
    return z, 2 * ss.norm.sf(abs(z)), runs

def ols(x, y):
    """slope, intercept, HC1 (heteroscedasticity-consistent) s.e. of slope"""
    x = np.asarray(x, float); y = np.asarray(y, float); nn = len(x)
    xb = x.mean(); sxx = ((x - xb) ** 2).sum()
    b = ((x - xb) * (y - y.mean())).sum() / sxx
    a = y.mean() - b * xb
    res = y - (a + b * x)
    se_ols = math.sqrt((res ** 2).sum() / (nn - 2) / sxx)
    se_hc1 = math.sqrt(nn / (nn - 2.0) * (((x - xb) ** 2) * res ** 2).sum()) / sxx
    return b, a, se_ols, se_hc1

print("=" * 78)
print("P2 at scale -- %d instances, n <= %d, k = %s"
      % (len(rows), max(int(r["n"]) for r in rows), KS))
bad = [r for r in rows if r["ok"] != "1"]
print("instances failing the exactness certificate `ok`: %d" % len(bad))
if bad[:5]:
    for r in bad[:5]: print("   ", r["k"], r["n"], "brk=", r["brk"])
OUT["not_ok"] = len(bad)
red = sum(1 for r in rows if int(r["redst"]) > 1)
print("instances where the ladder pair was NOT already Gauss-reduced "
      "(Theorem 1 re-verified): %d" % red)
OUT["ladder_not_reduced"] = red

# ---------------------------------------------------------------- 1. parity
print("\n" + "=" * 78 + "\n1. PARITY OF THE BALANCE INDEX\n")
OUT["parity"] = {}
par = {}
for k in KS:
    p = arr(k, "parity", int); nn = arr(k, "n", int)
    inc = arr(k, "incone", int)
    par[k] = (nn, p)
    ev = int((p == 0).sum()); tot = len(p)
    binom = ss.binomtest(ev, tot, 0.5).pvalue
    z, pv, nruns = runs_test(p == 0)
    lr_odd, at_odd = longest_run(nn, p == 1)
    rho = arr(k, "rho")
    lr_bad, at_bad = longest_run(nn, rho > S2)
    agree = int((inc == (1 - p)).sum())
    print("k=%-7.4f  n=%d..%d  even %d/%d = %.4f  binom p=%.3f | runs z=%+.2f "
          "p=%.3f | longest odd run %d | longest run rho>sqrt2 %d | "
          "parity<->cone %d/%d" % (k, nn.min(), nn.max(), ev, tot, ev / tot,
                                   binom, z, pv, lr_odd, lr_bad, agree, tot))
    print("           record odd run starts at n=%s ; record rho>sqrt2 run "
          "starts at n=%s ; gaps in n: %d"
          % (at_odd, at_bad, int((np.diff(nn) != 1).sum())))
    OUT["parity"][str(k)] = dict(n_lo=int(nn.min()), n_hi=int(nn.max()),
        even=ev, tot=tot, freq=ev / tot, binom_p=binom, runs_z=z, runs_p=pv,
        longest_odd_run=lr_odd, longest_rho_gt_sqrt2_run=lr_bad, odd_run_at=int(at_odd) if at_odd else 0, rho_run_at=int(at_bad) if at_bad else 0, gaps=int((np.diff(nn)!=1).sum()),
        parity_cone_agree=agree)

print("\nautocorrelation of the parity sequence (lags 1..24, |z| shown):")
OUT["parity_acf"] = {}
for k in KS:
    nn, p = par[k]; x = 2.0 * (p == 0) - 1
    zs = []
    for L in range(1, 25):
        c = float((x[:-L] * x[L:]).mean()); zs.append(c * math.sqrt(len(x) - L))
    print("  k=%-7.4f  max|z| = %.2f at lag %d   (all: %s)"
          % (k, max(abs(z) for z in zs), int(np.argmax(np.abs(zs))) + 1,
             " ".join("%+.1f" % z for z in zs)))
    OUT["parity_acf"][str(k)] = zs

print("\nn mod m, m = 2..60 (chi-square of even-frequency across residues), "
      "Bonferroni over %d tests per k:" % 59)
OUT["parity_mod"] = {}
for k in KS:
    nn, p = par[k]; hits = []
    for m in range(2, 61):
        tab = np.zeros((m, 2), int)
        for a, b in zip(nn % m, p): tab[a, b] += 1
        keep = tab.sum(axis=1) >= 5
        if keep.sum() < 2: continue
        chi2, pv, _, _ = ss.chi2_contingency(tab[keep])
        hits.append((pv, m))
    hits.sort()
    bonf = 59
    sig = [(m, pv) for pv, m in hits if pv * bonf < 0.05]
    print("  k=%-7.4f  smallest p: %s   surviving Bonferroni: %s"
          % (k, ", ".join("m=%d p=%.4f" % (m, pv) for pv, m in hits[:4]),
             sig if sig else "none"))
    OUT["parity_mod"][str(k)] = dict(smallest=[(m, pv) for pv, m in hits[:6]],
                                     bonferroni_survivors=sig)

print("\ncross-k agreement of the parity (on the common n range):")
OUT["cross_k"] = {}
for i in range(len(KS)):
    for j in range(i + 1, len(KS)):
        ka, kb = KS[i], KS[j]
        na, pa = par[ka]; nb, pb = par[kb]
        common = np.intersect1d(na, nb)
        aa = pa[np.isin(na, common)]; bb = pb[np.isin(nb, common)]
        ag = float((aa == bb).mean()); m = len(common)
        z = (ag - 0.5) * 2 * math.sqrt(m)
        print("  k=%.4f vs k=%.4f : %.4f on %d values  (z=%+.2f)" % (ka, kb, ag, m, z))
        OUT["cross_k"]["%s|%s" % (ka, kb)] = dict(agree=ag, n=m, z=z)

# ------------------------------------------------------- 2. headline slope
print("\n" + "=" * 78 + "\n2. THE CONE-MIN / TRUE-MIN LOG-RATIO  (headline)\n")
OUT["rho"] = {}
WINS = [(4, 120), (4, 200), (4, 500), (4, 1000), (4, 2000), (4, 3000), (4, 6000),
        (1001, 10**9), (4, 10**9)]
for k in KS:
    nn = arr(k, "n", int); rho = arr(k, "rho"); lr = np.log(rho)
    print("  k=%.4f  median rho=%.4f  max rho=%.1f  P[rho<=sqrt2]=%.4f  "
          "P[rho>20]=%.4f  median lam2/lam1=%.4f"
          % (k, float(np.median(rho)), float(rho.max()),
             float((rho <= S2).mean()), float((rho > 20).mean()),
             float(np.median(arr(k, "skew")))))
    OUT["rho"][str(k)] = dict(median=float(np.median(rho)), max=float(rho.max()),
        p_le_sqrt2=float((rho <= S2).mean()), p_gt_20=float((rho > 20).mean()),
        median_skew=float(np.median(arr(k, "skew"))), fits={})
    for lo, hi in WINS:
        sel = (nn >= lo) & (nn <= hi)
        if sel.sum() < 30: continue
        b, a, se, se1 = ols(nn[sel], lr[sel])
        lab = "n in [%d,%d]" % (lo, min(hi, int(nn.max())))
        print("     log rho ~ n on %-18s : slope %+.6f +- %.6f (HC1 %.6f), "
              "intercept %+.4f, N=%d" % (lab, b, se, se1, a, int(sel.sum())))
        OUT["rho"][str(k)]["fits"][lab] = dict(slope=b, se=se, se_hc1=se1,
                                               intercept=a, N=int(sel.sum()))
    # the mis-specified model of POSITIVITY_PROGRAM sec.3.3, for comparison
    b, a, se, se1 = ols(nn, lr / nn)
    print("     (1/n) log rho ~ n (mis-specified)  : slope %+.8f +- %.8f" % (b, se, ))
    OUT["rho"][str(k)]["fits"]["over_n_misspecified"] = dict(slope=b, se=se)

print("\n  pooled over all k, log rho ~ n with k fixed effects:")
allb = []
for k in KS:
    nn = arr(k, "n", int); lr = np.log(arr(k, "rho"))
    b, a, se, se1 = ols(nn, lr); allb.append((b, se))
wb = sum(b / se ** 2 for b, se in allb) / sum(1 / se ** 2 for b, se in allb)
wse = 1 / math.sqrt(sum(1 / se ** 2 for b, se in allb))
print("     slope = %+.6f +- %.6f   (2-sigma band %+.6f .. %+.6f)"
      % (wb, wse, wb - 2 * wse, wb + 2 * wse))
OUT["rho_pooled_slope"] = dict(slope=wb, se=wse)

# --------------------------------------------------------- 3. Gauss-Kuzmin
print("\n" + "=" * 78 + "\n3. GAUSS-KUZMIN AND THE PARTIAL QUOTIENT AT THE BALANCE\n")
GK = [math.log2(1 + 1.0 / (a * (a + 2))) for a in range(1, 9)]
GK.append(1 - sum(GK))
tot = np.zeros(9)
for r in rows:
    for a in range(9): tot[a] += int(r["q%d" % (a + 1)])
obs = tot / tot.sum()
chi2 = float(((tot - tot.sum() * np.array(GK)) ** 2 / (tot.sum() * np.array(GK))).sum())
print("  %d partial quotients of h12/h11 in all" % int(tot.sum()))
print("  a         :" + "".join("%8d" % a for a in range(1, 9)) + "     >=9")
print("  observed  :" + "".join("%8.4f" % x for x in obs))
print("  Gauss-Kuzmin:" + "".join("%6.4f  " % x for x in GK))
print("  chi^2 = %.2f on 8 df   (p = %.3f)" % (chi2, ss.chi2.sf(chi2, 8)))
OUT["gauss_kuzmin"] = dict(total=int(tot.sum()), observed=list(obs),
                          predicted=GK, chi2=chi2, p=float(ss.chi2.sf(chi2, 8)))

print("\n  the partial quotient AT the balance index, a_{i(n)}, and the next, "
      "a_{i(n)+1}:")
OUT["pq_at"] = {}
for lab, col in (("a_{i}", "pqat"), ("a_{i+1}", "pqnext")):
    v = np.array([int(r[col]) for r in rows], float)
    v = v[v > 0]
    h = np.array([float((v == a).mean()) for a in range(1, 9)] + [float((v >= 9).mean())])
    c2 = float((((h - np.array(GK)) ** 2) / np.array(GK)).sum() * len(v))
    print("    %-8s N=%d  mean=%.3f  median=%.1f  max=%d  "
          "P[>=9]=%.4f (GK %.4f)  chi^2=%.2f (p=%.3f)"
          % (lab, len(v), v.mean(), float(np.median(v)), int(v.max()),
             h[8], GK[8], c2, ss.chi2.sf(c2, 8)))
    print("      histogram :" + "".join("%8.4f" % x for x in h))
    OUT["pq_at"][lab] = dict(N=len(v), mean=float(v.mean()),
        median=float(np.median(v)), max=int(v.max()), hist=list(h),
        chi2=c2, p=float(ss.chi2.sf(c2, 8)))
for col, lab in (("pqw5", "max a over i+-5"), ("pqw20", "max a over i+-20"),
                 ("pqmax", "max a over the whole CF")):
    v = np.array([int(r[col]) for r in rows], float)
    lg = np.log(v)
    nn = np.array([int(r["n"]) for r in rows], float)
    b, a, se, se1 = ols(nn, lg)
    print("    %-22s median=%8.1f  max=%10d   log ~ n slope %+.2e +- %.2e"
          % (lab, float(np.median(v)), int(v.max()), b, se))
    OUT["pq_at"][lab] = dict(median=float(np.median(v)), max=int(v.max()),
                             log_slope=b, log_slope_se=se)

# ------------------------------------------------- 1b. run-length null model
print("\n  run lengths against a fair coin (10^5 Monte-Carlo replicates of the "
      "same length and even-frequency):")
OUT["runs_null"] = {}
rng = np.random.default_rng(20260823)
for k in KS:
    nn = arr(k, "n", int); p = arr(k, "parity", int)
    obs, _ = longest_run(nn, p == 1)
    m = len(p); q = float((p == 1).mean())
    sim = rng.random((100000, m)) < q
    L = np.zeros(100000, int); cur = np.zeros(100000, int)
    for j in range(m):
        cur = np.where(sim[:, j], cur + 1, 0)
        L = np.maximum(L, cur)
    pv = float((L >= obs).mean())
    print("    k=%-7.4f  longest odd run %2d ; null mean %.2f, median %.0f, "
          "P[null >= observed] = %.3f ; log2(N) = %.1f"
          % (k, obs, L.mean(), float(np.median(L)), pv, math.log2(m)))
    OUT["runs_null"][str(k)] = dict(observed=obs, null_mean=float(L.mean()),
        null_median=float(np.median(L)), p=pv, log2N=math.log2(m))

# ------------------------------------------- 2b. tail of rho and Corollary D
print("\n  tail of rho against the Haar law P[rho>R] ~ 0.54/R, and Corollary D:")
OUT["rho_tail"] = {}
for k in KS:
    rho = arr(k, "rho"); skew = arr(k, "skew"); inc = arr(k, "incone", int)
    tail = {R: float((rho > R).mean()) for R in (5, 10, 20, 50, 100)}
    viol = int(((rho <= S2) & (inc == 0) & (skew > 2 * math.sqrt(2) / math.sqrt(3))).sum())
    dcheck = int(((rho <= S2) & (inc == 0)).sum())
    pc = int(((inc == 1) & (rho > S2)).sum())
    print("    k=%-7.4f  " % k + "  ".join("P[rho>%d]=%.4f (Haar %.4f)"
          % (R, v, 0.54 / R) for R, v in tail.items()))
    print("             rho<=sqrt2 with v1 outside the cone: %d (all have "
          "lam2/lam1<=1.633? %s) ; v1 in cone with rho>sqrt2: %d "
          "(Proposition C says 0)" % (dcheck, viol == 0, pc))
    OUT["rho_tail"][str(k)] = dict(tail=tail, corD=dcheck, corD_violations=viol,
                                   propC_violations=pc)

# --------------------------------------- 3b. lam2/lam1 against a_{i+1}
print("\n  the P2'_par equivalence: lam2/lam1 against a_{i(n)+1}")
OUT["skew_vs_pq"] = {}
for k in KS:
    sk = arr(k, "skew"); aq = arr(k, "pqnext"); nn = arr(k, "n", int)
    good = aq > 0
    rr = sk[good] / aq[good]
    c = float(np.corrcoef(np.log(sk[good]), np.log(aq[good]))[0, 1])
    b, a, se, se1 = ols(nn[good], np.log(aq[good]))
    print("    k=%-7.4f  lam2/lam1 in [%.3f,%.3f]*a_{i+1} ; corr of logs %.3f ; "
          "max (log a_{i+1})/n = %.4f ; log a_{i+1} ~ n slope %+.2e +- %.2e"
          % (k, rr.min(), rr.max(), c,
             float((np.log(aq[good]) / nn[good]).max()), b, se))
    OUT["skew_vs_pq"][str(k)] = dict(ratio_lo=float(rr.min()), ratio_hi=float(rr.max()),
        corr=c, max_log_a_over_n=float((np.log(aq[good]) / nn[good]).max()),
        slope=b, se=se)

# ------------------------------------------------------------- 4. i(n)/L(n)
print("\n" + "=" * 78 + "\n4. THE BALANCE POSITION i(n)/L(n)  (predicted 0.5238)\n")
OUT["irat"] = {}
for k in KS:
    nn = arr(k, "n", int); ir = arr(k, "irat")
    b, a, se, se1 = ols(nn, ir)
    hi = ir[nn > nn.max() // 2]
    print("  k=%-7.4f  mean %.5f +- %.5f (sd)   slope in n %+.3e +- %.3e   "
          "mean on the top half %.5f   L(n)/n = %.3f"
          % (k, ir.mean(), ir.std(ddof=1), b, se, hi.mean(),
             float((arr(k, "len") / nn).mean())))
    OUT["irat"][str(k)] = dict(mean=float(ir.mean()), sd=float(ir.std(ddof=1)),
        slope=b, se=se, mean_top_half=float(hi.mean()),
        L_over_n=float((arr(k, "len") / nn).mean()))

# ------------------------------------------------------------------ 5. 2-adic
print("\n" + "=" * 78 + "\n5. v_2(q_n xi_2 - p_n) AGAINST THE (24.02-k) n LAW\n")
OUT["adelic"] = {}
for k in KS:
    nn = arr(k, "n", int); v2 = arr(k, "v2qp")
    b, a, se, se1 = ols(nn, v2)
    late = nn > nn.max() // 2
    r = float((v2[late] / nn[late]).mean())
    print("  k=%-7.4f  v2/n on the top half = %.4f   (24.02-k = %.4f)   "
          "fit v2 = %.5f n %+.2f  (s.e. %.5f)   k+slope = %.4f"
          % (k, r, 24.02 - k, b, a, se, k + b))
    OUT["adelic"][str(k)] = dict(v2_over_n_top_half=r, pred=24.02 - k,
        slope=b, intercept=a, se=se, k_plus_slope=k + b,
        dinf=float(arr(k, "dinf")[late].mean()),
        d2=float(arr(k, "d2")[late].mean()))
    print("           delta_inf = %.4f  delta_2 = %.4f  sum = %.4f (top half)"
          % (OUT["adelic"][str(k)]["dinf"], OUT["adelic"][str(k)]["d2"],
             OUT["adelic"][str(k)]["dinf"] + OUT["adelic"][str(k)]["d2"]))

json.dump(OUT, open(os.path.join(HERE, "data", "summary.json"), "w"), indent=1)
print("\nwritten data/summary.json")
