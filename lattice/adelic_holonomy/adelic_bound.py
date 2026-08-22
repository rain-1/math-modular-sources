"""Adelic arithmetic holonomy bound: CDT (arXiv:2408.15403, Thm 6.0.2/7.0.1)
with a p-adic Liouville gain inserted.

CDT:   m <= BC(phi) / ( log|phi'(0)| - tau(b;e) ),  tau = tau^flat(b)+tau^sharp(e).
Here:  m <= BC(phi) / ( log|phi'(0)| - tau(b;e) + sum_p gamma_p )   with

   gamma_p = log p * [ s_min*(1-1/m) + (2/m) * int_0^1 delta^down(t)*max(m t -1,0) dt ]

where s_i = p-adic slope of f_i (v_p(coeff of x^n) >= s_i*n - o(n)),
s_min = min_i s_i, delta_i = s_i - s_min, delta^down the non-increasing
rearrangement viewed as a step function on [0,1] with steps of width 1/m.

Step profile (u functions of slope s, rest 0):  gamma = s*log p*(u-1)^2/m^2.
Uniform profile (all slope s):                  gamma = s*log p*(1-1/m).
Both are special cases of the formula above.  See consolidation/ADELIC_HOLONOMY.md.
"""
import math, sys, itertools
sys.path.insert(0, __file__.rsplit('/',1)[0].replace('adelic_holonomy','cdt_finder'))
from cdt_bound import tau_flat, tau_sharp, I_uvw
from fractions import Fraction

def gamma_p(m, slopes, p):
    """slopes: list of m p-adic slopes (reals, may be negative)."""
    assert len(slopes) == m
    smin = min(slopes)
    delta = sorted([s - smin for s in slopes], reverse=True)   # non-increasing
    # (2/m) * sum_i delta_i * int_{(i-1)/m}^{i/m} max(m t - 1,0) dt
    tot = 0.0
    for i, d in enumerate(delta, start=1):
        a, b = (i-1)/m, i/m
        lo = max(a, 1.0/m)
        if b > lo:
            tot += d*((m*(b*b-lo*lo)/2) - (b-lo))
    return math.log(p)*(smin*(1-1/m) + (2.0/m)*tot)

def adelic(m, cols, e, slopes_by_p, log_phiprime, BC):
    sm, tf = tau_flat(m, cols)
    ts, xi = tau_sharp(m, e)
    tau = float(tf) + ts
    g = sum(gamma_p(m, sl, p) for p, sl in slopes_by_p.items())
    entry = log_phiprime - tau + g
    return dict(sigma_m=float(sm), tau_flat=float(tf), tau_sharp=ts, tau=tau,
                gamma=g, entry=entry, BC=BC,
                bound=(BC/entry if entry > 0 else float('inf')),
                margin=m*entry - BC)

def fmt(name, r, m):
    print(f"{name:<46s} m={m:3d}  tau={r['tau']:.4f}  gamma={r['gamma']:+.4f}  "
          f"entry={r['entry']:+.4f}  bnd={r['bound']:8.3f}  margin={r['margin']:+9.3f}")

if __name__ == '__main__':
    # ---- sanity 1: gamma=0 reproduces CDT exactly -------------------------
    m = 14; cols=[(1,2),(3,2)]; e=[0,0,1,0,0,0,0,0,0,1,1,1,1,1]
    psi = 5448339453535586608000000000/8658833407565631122430056127
    lp = math.log(256*psi)
    r = adelic(m, cols, e, {}, lp, 11.845)
    print("CDT calibration (gamma=0):")
    fmt("  level 6, L(2,chi_-3)", r, m); print("  (CDT: bound 13.9938, margin +0.0053)\n")

    # ---- sanity 2: scale invariance --------------------------------------
    # rescaling the coordinate y -> y/p^v : log|phi'|-=v log p, BC-=v log p,
    # every slope += v.  Margin must be unchanged.
    m=21; cols=[(1,2),(3,2)]; e=[0]*21
    sl=[8]*7+[2]*7+[0]*7
    lp0, BC0 = math.log(64), 11.845+math.log(0.25)
    a = adelic(m, cols, e, {2:sl}, lp0, BC0)
    v=3; a2 = adelic(m, cols, e, {2:[s+v for s in sl]}, lp0-v*math.log(2), BC0-v*math.log(2))
    print(f"scale invariance: margin {a['margin']:.9f} vs {a2['margin']:.9f}  "
          f"(diff {a['margin']-a2['margin']:.2e})\n")

# ---------------------------------------------------------------------------
def host_table():
    """Level-8 Catalan host (Zagier E) and X_1(5) Sym^2, adelic margins."""
    import math
    print("="*100)
    print("LEVEL-8 CATALAN HOST (Zagier E, lambda_1=8, lambda_2=4, s=1/4, k=2)")
    L2 = math.log(2)
    ceil_ = math.log(256*0.25); real_ = math.log(256*0.25*0.6292232680)
    BC   = 11.845 + math.log(0.25)
    print(f"  ceiling log|phi'(0)| = log(256 s) = {ceil_:.6f};  realised (CDT loss .62922) = {real_:.6f}")
    print(f"  BC = 11.845 + log s = {BC:.6f}")
    print("  2-adic slopes in the symmetrised coordinate y (sigma_y >= max(0, 2 sigma_x - 2)):")
    print("     conditional orbit  sigma_x = 0  -> sigma_y = 0   (xi_2 = zeta_2(2)/2 != xi_infty)")
    print("     pure polylog orbit sigma_x = 2  -> sigma_y = 2   (Li_j(4x), v_2(4^n/n^j)=2n-O(log n))")
    print("     doubly-small orbit sigma_x = 3  -> sigma_y = 4   (measured on B_0=(1-4V_2)E companion)")
    print("     hypothetical best  sigma_x = 5  -> sigma_y = 8   (the Casoratian slope of B_E - xi_2 A)")
    for sy, lab in ((4,'measured sigma_y=4'), (8,'hypothetical sigma_y=8')):
        print(f"\n  --- {lab} ---")
        print(f"  {'p_dbl':>5} {'p_pure':>6} {'p_cond':>6} {'m':>4} {'tau':>7} {'gamma':>8} "
              f"{'entryC':>8} {'entryR':>8} {'marginC':>9} {'marginR':>9} {'cap m0':>7}")
        for pd in (0,2,3,5,7):
            for pp, pc in ((7,7),):
                m = pd+pp+pc
                if m == 0: continue
                cols = [(1,2),(3,2)]
                e = [0]*m
                for j in range(min(6,m)): e[-(j+1)] = 1     # CDT-proportional integration profile
                sl = [sy]*pd + [2]*pp + [0]*pc
                rC = adelic(m, cols, e, {2:sl}, ceil_, BC)
                rR = adelic(m, cols, e, {2:sl}, real_, BC)
                # unconditional cap on the doubly-small module alone
                cap = None
                for m0 in range(1, 40):
                    r0 = adelic(m0, [(0,2),(0,2)], [0]*m0, {2:[sy]*m0}, ceil_, BC)
                    if r0['entry'] > 0 and m0 > r0['bound']: cap = m0-1; break
                print(f"  {pd:>5} {pp:>6} {pc:>6} {m:>4} {rC['tau']:>7.4f} {rC['gamma']:>+8.4f} "
                      f"{rC['entry']:>+8.4f} {rR['entry']:>+8.4f} {rC['margin']:>+9.3f} {rR['margin']:>+9.3f} {str(cap):>7}")
    print("\n  reference (CDT inventory, gamma=0, m=14):")
    r = adelic(14, [(1,2),(3,2)], [0,0,1,0,0,0,0,0,0,1,1,1,1,1], {}, real_, BC)
    fmt("    Zagier E archimedean-only", r, 14)
    r = adelic(14, [(1,2),(3,2)], [0,0,1,0,0,0,0,0,0,1,1,1,1,1], {}, ceil_, BC)
    fmt("    Zagier E archimedean-only @ceiling", r, 14)

    print("\n"+"="*100)
    print("STRUCTURAL CEILING:  log(256/lambda_2) + sigma_p log p  <=  log(256 lambda_1)")
    for nm, l1, l2, c in (("Zagier E (Catalan)",8,4,32), ("Zagier C = CDT",9,1,9),
                          ("Zagier F",9,8,72), ("Apery zeta(3)",1,1,1),
                          ("X_1(5) Sym^2",11.09,1,-1), ("Domb",16,4,64)):
        best = max((v*math.log(p) for p in (2,3,5,7,11)
                    for v in [0] + [k for k in range(1,12) if abs(c) % p**k == 0]), default=0.0)
        print(f"  {nm:<20s} lambda_1={l1:>6}, lambda_2={l2:>4}, c={c:>4}: "
              f"ceil={math.log(256/l2):.4f}  max gain={best:.4f}  total={math.log(256/l2)+best:.4f} "
              f" (log(256*lambda_1)={math.log(256*l1):.4f})")
