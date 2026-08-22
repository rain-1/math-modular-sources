"""Two-prime functions in the adelic arithmetic holonomy bound.

Reuses lattice/cdt_finder/cdt_bound.py and lattice/adelic_holonomy/adelic_bound.py.
See consolidation/TWO_PRIME_HOLONOMY.md.
"""
import math, sys, os
R = '/home/ubuntu/code/math-modular-sources/lattice'
sys.path.insert(0, R + '/cdt_finder'); sys.path.insert(0, R + '/adelic_holonomy')
from cdt_bound import tau_flat, tau_sharp
from adelic_bound import gamma_p, adelic
from fractions import Fraction

L2, L3 = math.log(2), math.log(3)
CDT_E = [0,0,1,0,0,0,0,0,0,1,1,1,1,1]           # CDT's integration vector, m=14
def e_of(m):                                     # CDT-proportional: six 1's at the end
    e = [0]*m
    for j in range(min(6, m)): e[-(j+1)] = 1
    return e

HOSTS = {   # name: (s = 1/lambda_2, lambda_2, contour-loss factor)
 'CDT level 6   (Zagier C, lam2=1)': (1.0,  1, 0.6292232680),
 'Catalan lvl 8 (Zagier E, lam2=4)': (0.25, 4, 0.6292232680),
 'Catalan lvl16 (lam2=4, s=-1/4)  ': (0.25, 4, 0.6292232680),
}
def geom(s, loss):
    return math.log(256*s), math.log(256*s*loss), 11.845 + math.log(s)

def show(tag, m, cols, e, slopes, lphi, BC):
    r = adelic(m, cols, e, slopes, lphi, BC)
    print(f"  {tag:<52s} m={m:2d} tau={r['tau']:7.4f} gam={r['gamma']:+7.4f} "
          f"entry={r['entry']:+8.4f} margin={r['margin']:+9.3f}")
    return r

print("="*104)
print("1.  CALIBRATION AND PART (a): CDT's OWN 14 FUNCTIONS, MEASURED SLOPES")
print("="*104)
ce, cr, cB = geom(*HOSTS['CDT level 6   (Zagier C, lam2=1)'][::2][:1] + (1,), )  if False else geom(1.0, 0.6292232680)
print(f"  log|phi'(0)| ceiling = {ce:.6f}, CDT realised = {cr:.6f}, BC = {cB:.6f}")
show("gamma=0 (CDT's own numbers), realised contour", 14, [(1,2),(3,2)], CDT_E, {}, cr, cB)
show("MEASURED slopes: all 14 have varsigma_2=varsigma_3=0", 14, [(1,2),(3,2)], CDT_E,
     {2:[0.0]*14, 3:[0.0]*14}, cr, cB)
show("same, at the hard ceiling", 14, [(1,2),(3,2)], CDT_E, {2:[0.0]*14, 3:[0.0]*14}, ce, cB)

print("\n  counterfactual: if u of the 14 realised the host's 3-adic resource")
print("  (row C has c=9, sigma_3 = v_3(9) = 2 in x, hence varsigma_y = 4) -- NOT REALISABLE,")
print("  the overconvergent direction B - xi_3 A has 3-adically irrational coefficients:")
for u in (1,3,5,7,10,14):
    sl = [4.0]*u + [0.0]*(14-u)
    show(f"    u={u} functions with varsigma_3=4 (counterfactual)", 14, [(1,2),(3,2)], CDT_E,
         {3:sl}, cr, cB)

print("\n" + "="*104)
print("2.  THE PURE-MODULE MULTI-PRIME IDENTITY")
print("="*104)
print("  pure module Li_j(x/s), s=1/lambda_2:  varsigma_p = v_p(lambda_2), so")
print("  sum_p varsigma_p log p = log|lambda_2|  exactly, and")
print("  [archimedean ceiling] + [full multi-prime gain at weight 1] = log(256/|l2|) + log|l2| = log 256.")
print(f"  {'host':<26s} {'l2':>5s} {'ceil':>8s} {'sum v_p log p':>14s} {'total':>8s} "
      f"{'gam(u=7,m=14)':>14s} {'entry':>9s} {'margin':>9s}")
for nm, l2 in (("Apery/CDT/X1(5)", 1), ("Catalan lvl 8/16", 4), ("Zagier F", 8),
               ("Domb", 4), ("Cooper s_18", 12), ("sqrt(s_18)", 24), ("sqrt(s_10)", 8)):
    s = 1.0/l2
    prim = {p: 0 for p in (2,3,5)}
    n = l2
    for p in (2,3,5):
        while n % p == 0: n //= p; prim[p] += 1
    tot = sum(v*math.log(p) for p, v in prim.items())
    ceil_ = math.log(256*s)
    sl = {p: [float(v)]*7 + [0.0]*7 for p, v in prim.items() if v}
    r = adelic(14, [(1,2),(3,2)], CDT_E, sl, ceil_, 11.845 + math.log(s))
    print(f"  {nm:<26s} {l2:>5d} {ceil_:>8.4f} {tot:>14.4f} {ceil_+tot:>8.4f} "
          f"{r['gamma']:>+14.4f} {r['entry']:>+9.4f} {r['margin']:>+9.3f}")
print("  => the multi-prime total is log 256 = 5.5452 on every host; but the p-adic half enters")
print("     the entry condition at weight (u-1)^2/m^2 = %.5f and the archimedean half at weight 1,"
      % (36/196))
print("     so lambda_2 = 1 (the Apery-perfect hosts) is optimal and every prime factor of lambda_2")
print("     costs (1 - (u-1)^2/m^2) v_p log p = %.4f per v_p(lambda_2) log p." % (1-36/196))

print("\n" + "="*104)
print("3.  ADMITTING A WELL-POISED DECAYER'S GENERATING FUNCTION")
print("="*104)
# (name, alpha, kappa_2, kappa_3, nu(Q), logLambda)   [measured, see the .md]
DEC = [
 ("conductor 3,  alpha=2  ", 0.00, 3.000,  0.00,  6.9315),
 ("conductor 6,  alpha=1/4", 4.410, 0.000,  9.840, -10.198),
 ("conductor 6,  alpha=1/2", 4.950, 0.000,  8.340,  -7.371),
 ("conductor 6,  alpha=1  ", 5.830, 0.830,  5.480,  -2.062),
 ("conductor 6,  alpha=2  ", 7.820, 2.850,  0.660,  6.931472),   # log 1024, EXACT
 ("conductor 12, alpha=1/4", 2.000, 0.000, 29.501, -23.541),
 ("conductor 12, alpha=1/2", 4.000, 0.000, 27.283, -18.711),
 ("conductor 12, alpha=1  ", 8.000, 1.733, 20.101,  -8.214),
 ("conductor 12, alpha=2  ",16.000, 5.733,  1.896,  13.862944),   # log 2^20, EXACT
]
def add_decayer(hostname, s, loss, m0, cols0, slopes0, k2, k3, nu, logLam, tag):
    """Add one function of LCM-denominator rate nu and geometric p-parts kappa_p."""
    m = m0 + 1
    cols = list(cols0) + ([(m-1, nu)] if nu > 0 else [])
    e = e_of(m)
    ceil_, real_, BC = geom(s, loss)
    # optimal integral rescaling x -> x/N : N = smallest integer with N/Lambda >= s
    logN = max(0.0, logLam + math.log(s))
    k2p, k3p = k2 + logN/L2*0.0, k3          # log N distributed over primes; keep total
    extra = logN                              # nats of extra geometric denominator
    sl2 = [-(k2*L2 + extra)/L2] + [x for x in slopes0.get(2, [0.0]*m0)]
    sl3 = [-k3] + [x for x in slopes0.get(3, [0.0]*m0)]
    r = adelic(m, cols, e, {2: sl2, 3: sl3}, ceil_, BC)
    rr = adelic(m, cols, e, {2: sl2, 3: sl3}, real_, BC)
    g2 = gamma_p(m, sl2, 2); g3 = gamma_p(m, sl3, 3)
    print(f"  {tag:<24s} nu={nu:6.2f} k2={k2:6.2f} k3={k3:5.2f} logN={logN:6.3f} | "
          f"tau={r['tau']:7.4f} g2={g2:+7.4f} g3={g3:+7.4f} "
          f"entryC={r['entry']:+8.4f} marginC={r['margin']:+8.2f} entryR={rr['entry']:+8.4f}")
    return r, rr

for hostname, (s, l2, loss) in HOSTS.items():
    if 'lvl16' in hostname: continue
    print(f"\n  --- host {hostname},  s={s} ---")
    slopes0 = {2: ([2.0]*7 + [0.0]*7) if l2 == 4 else [0.0]*14, 3: [0.0]*14}
    ceil_, real_, BC = geom(s, loss)
    base = adelic(14, [(1,2),(3,2)], CDT_E, slopes0, ceil_, BC)
    baseR = adelic(14, [(1,2),(3,2)], CDT_E, slopes0, real_, BC)
    print(f"  {'BASELINE m=14 (no extra function)':<24s} {'':38s} | "
          f"tau={base['tau']:7.4f} g2={gamma_p(14,slopes0[2],2):+7.4f} g3={0.0:+7.4f} "
          f"entryC={base['entry']:+8.4f} marginC={base['margin']:+8.2f} entryR={baseR['entry']:+8.4f}")
    for tag, k2, k3, nu, lLam in DEC:
        add_decayer(hostname, s, loss, 14, [(1,2),(3,2)], slopes0, k2, k3, nu, lLam, tag)

print("\n  --- level-16 Catalan host (14 CDT functions + the doubly-small function, m=15) ---")
s, loss = 0.25, 0.6292232680
ceil_, real_, BC = geom(s, loss)
sl15 = {2: [1.0] + [2.0]*7 + [0.0]*7, 3: [0.0]*15}
b15 = adelic(15, [(1,2),(3,2)], e_of(15), sl15, ceil_, BC)
b15R = adelic(15, [(1,2),(3,2)], e_of(15), sl15, real_, BC)
print(f"  {'BASELINE m=15 (level 16)':<24s} {'':38s} | tau={b15['tau']:7.4f} "
      f"g2={gamma_p(15,sl15[2],2):+7.4f} entryC={b15['entry']:+8.4f} "
      f"marginC={b15['margin']:+8.2f} entryR={b15R['entry']:+8.4f} marginR={b15R['margin']:+8.2f}")
for tag, k2, k3, nu, lLam in DEC:
    add_decayer('lvl16', s, loss, 15, [(1,2),(3,2)], sl15, k2, k3, nu, lLam, tag)

print("\n" + "="*104)
print("4.  CONTOUR COST OF AN EXTRA SINGULAR POINT")
print("="*104)
print("  If the added function is singular at y = q (nontrivial monodromy) then phi factors")
print("  through the universal cover of P^1 minus the enlarged singular set, which is contained")
print("  in P^1 \\ {0, q, infinity}; hence  |phi'(0)| <= 16|q|  (256|q| if the Z/2 descent is")
print("  available with q as the outer point).  q = 1/Lambda for the decayer's own variable.")
print(f"  {'family':<24s} {'|q|=1/Lambda':>14s} {'log(16|q|)':>11s} {'log(256|q|)':>12s} "
      f"{'CDT tau':>8s} {'shortfall':>10s}")
for tag, k2, k3, nu, lLam in DEC:
    q = math.exp(-lLam)
    print(f"  {tag:<24s} {q:>14.4e} {math.log(16*q):>11.4f} {math.log(256*q):>12.4f} "
          f"{4.2355:>8.4f} {math.log(256*q)-4.2355:>10.4f}")
print("  (CDT's own configuration, q = s = 1: the 3-point bound gives log 256 = 5.5452 and CDT")
print("   realise 5.0819 -- the bound is not vacuous.)")

print("\n" + "="*104)
print("5.  BREAK-EVEN: what an extra 15th function would have to be")
print("="*104)
for hostname, (s, l2, loss) in (("CDT level 6", (1.0,1,0.6292232680)),
                                ("Catalan level 8", (0.25,4,0.6292232680))):
    ceil_, real_, BC = geom(s, loss)
    slopes0 = {2: ([2.0]*7 + [0.0]*7) if l2 == 4 else [0.0]*14}
    base = adelic(14, [(1,2),(3,2)], CDT_E, slopes0, ceil_, BC)
    lo, hi = 0.0, 200.0
    for _ in range(80):
        mid = (lo+hi)/2
        r = adelic(15, [(1,2),(3,2),(14,mid)], e_of(15),
                   {2: [0.0] + slopes0[2]}, ceil_, BC)
        if r['margin'] >= base['margin']: lo = mid
        else: hi = mid
    print(f"  {hostname}: a 15th function with p-adic slope 0 pays for itself iff its LCM")
    print(f"     denominator rate nu <= {lo:.3f}  (baseline margin {base['margin']:+.3f} at the ceiling).")
    lo2, hi2 = -60.0, 30.0
    for _ in range(120):
        mid = (lo2+hi2)/2
        r = adelic(15, [(1,2),(3,2)], e_of(15),
                   {2: [mid] + slopes0[2]}, ceil_, BC)
        if r['margin'] < base['margin']: lo2 = mid
        else: hi2 = mid
    print(f"     with nu = 0 it pays for itself iff its 2-adic slope varsigma_2 >= {hi2:+.3f}.")
    # how many extra functions of 2-adic slope vs to reach margin 0
    for vs in (1.0, 2.0, 3.0, 4.0, 8.0):
        d = None
        for k in range(0, 80):
            m = 14+k
            r = adelic(m, [(1,2),(3,2)], e_of(m),
                       {2: [vs]*k + slopes0[2]}, ceil_, BC)
            if r['margin'] >= 0: d = k; break
        print(f"     extra functions of 2-adic slope {vs:4.1f} needed to reach margin 0: {d}")
