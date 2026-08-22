"""12_adelic.py -- the (adelic) arithmetic holonomy bound on the Hadamard host
of (Zudilin at 3n, Nesterenko (4,7) at n).

All inputs are the MEASURED quantities of lattice/hadamard_host/ (see
consolidation/HADAMARD_HOST.md):
  singular set  {s1t1, s1t2, s2t1, s2t2} = {3.4875690e-7, 1.0971741, -0.64886152, -2041290.24}
  denominator type: [1..6n]^k with k measured per function (b_j = 6)
  2-adic slopes: W = -2, everything else = -26
Reuses lattice/adelic_holonomy/adelic_bound.py (calibrated on CDT).
"""
import math, sys, os
BASE = "/home/ubuntu/code/math-modular-sources/lattice"
sys.path.insert(0, BASE + "/adelic_holonomy")
sys.path.insert(0, BASE + "/cdt_finder")
from adelic_bound import adelic, gamma_p

# ---- measured singular points (06_operator.gp, exact roots of Pgf) ----------
s1t1 =  3.4875690311312326e-07
s1t2 =  1.0971741115134826
s2t1 = -0.6488615207279722
s2t2 = -2041290.2401220216

def rho_two_ray(a, b):
    """conformal radius at 0 of  C \\ ( (-inf,-a] u [b,+inf) ),  a,b > 0.
       phi = M^{-1} o Koebe;  M(w) = -(w+a)/(w-b);  rho = 4ab/(a+b).
       b = +inf gives the one-ray Koebe value 4a."""
    if b == float('inf'): return 4.0*a
    if a == float('inf'): return 4.0*b
    return 4.0*a*b/(a+b)

# domains
A_W    = abs(s2t1)          # nearest singularity of W on the left
B_W    = s1t2               # nearest on the right
A_C    = abs(s2t2)          # COND_Z: left singularity is only s2t2
B_C    = s1t2
A_D    = abs(s2t2)          # DBL: only s2t2 at all
B_D    = float('inf')
A_ALL  = abs(s2t1)          # everything, incl. A_Z(*)A_N
B_ALL  = s1t1

rho_W   = rho_two_ray(A_W, B_W)
rho_C   = rho_two_ray(A_C, B_C)
rho_D   = rho_two_ray(A_D, B_D)
rho_ALL = rho_two_ray(A_ALL, B_ALL)

print("="*104)
print("DOMAINS  (univalent architecture (K); BC = log|phi'(0)| exactly, Grunsky)")
print(f"  singular set of the host: {s1t1:.8e}, {s1t2:.8f}, {s2t1:.8f}, {s2t2:.6f}")
for nm, a, b, r in (("W and its theta-orbit", A_W, B_W, rho_W),
                    ("COND_Z = A_Z(*)(bB_N-aA_N)", A_C, B_C, rho_C),
                    ("DBL = (bB_Z-aA_Z)(*)(bB_N-aA_N)", A_D, B_D, rho_D),
                    ("everything (incl. A_Z(*)A_N)", A_ALL, B_ALL, rho_ALL)):
    print(f"  {nm:<34s} slits at -{a:.6g} and +{b:.6g}   rho = {r:.8g}   log rho = {math.log(r):+9.5f}"
          f"   Landau ceiling log(16 d_min) = {math.log(16*min(a,b)):+9.5f}")

# ---- inventories -----------------------------------------------------------
# cols = [(u_j, b_j)] one entry per LCM layer; e = integration profile (all 0 here)
INV = [
 # name,                              m, cols,                          slopes,             log|phi'|,  BC
 ("(K) {1,W,thW}",                    3, [(1,6),(1,6)],                 [0,-2,-2],          math.log(rho_W), math.log(rho_W)),
 ("(K) {1,W,thW,th2W}",               4, [(1,6)]*2,                     [0,-2,-2,-2],       math.log(rho_W), math.log(rho_W)),
 ("(K) {1,W,thW,th2W,th3W}",          5, [(1,6)]*2,                     [0,-2,-2,-2,-2],    math.log(rho_W), math.log(rho_W)),
 ("(K) {1,COND_Z,thCOND_Z}",          3, [(1,6),(1,6)],                 [0,-26,-26],        math.log(rho_C), math.log(rho_C)),
 ("(K) {1,DBL,thDBL}",                3, [(1,6)]*4,                     [0,-26,-26],        math.log(rho_D), math.log(rho_D)),
 ("(K) {1,W,thW,COND,thCOND}",        5, [(1,6),(1,6)],                 [0,-2,-2,-26,-26],  math.log(rho_W), math.log(rho_W)),
 ("(K) {1,W,thW,DBL,thDBL}",          5, [(1,6),(1,6),(3,6),(3,6)],     [0,-2,-2,-26,-26],  math.log(rho_W), math.log(rho_W)),
 ("(K) {1,AZ*AN,W,thW}",              4, [(1,6),(1,6)],                 [0,-26,-2,-2],      math.log(rho_ALL), math.log(rho_ALL)),
]
print()
print("="*104)
print("ADELIC BOUND, architecture (K) (univalent, rigorous end-to-end)")
print(f"{'inventory':<32s} {'m':>2} {'sigma_m':>7} {'tau^b':>8} {'gamma_2':>9} {'entry':>9} {'entry(no 2adic)':>16} {'margin':>10}")
for nm, m, cols, sl, lp, bc in INV:
    r  = adelic(m, cols, [0]*m, {2: sl}, lp, bc)
    r0 = adelic(m, cols, [0]*m, {}, lp, bc)
    print(f"{nm:<32s} {m:>2} {r['sigma_m']:>7.2f} {r['tau_flat']:>8.4f} {r['gamma']:>+9.4f} "
          f"{r['entry']:>+9.4f} {r0['entry']:>+16.4f} {r['margin']:>+10.3f}")

print()
print("="*104)
print("BEST CONCEIVABLE on this host: the Landau/Kodaira ceiling log(16 d_min), archimedean only,")
print("and with the 2-adic term.  (Not attainable: the covering map of a >3-punctured sphere is smaller.)")
for nm, m, cols, sl, d in (("{1,W,thW}", 3, [(1,6),(1,6)], [0,-2,-2], min(A_W,B_W)),
                           ("{1,COND_Z,thCOND_Z}", 3, [(1,6),(1,6)], [0,-26,-26], min(A_C,B_C)),
                           ("{1,DBL,thDBL}", 3, [(1,6)]*4, [0,-26,-26], min(A_D,B_D))):
    lp = math.log(16*d)
    r  = adelic(m, cols, [0]*m, {2: sl}, lp, lp)
    r0 = adelic(m, cols, [0]*m, {}, lp, lp)
    print(f"  {nm:<24s} log(16 d) = {lp:+8.4f}   tau = {r['tau']:8.4f}   entry_adelic = {r['entry']:+8.4f}"
          f"   entry_arch = {r0['entry']:+8.4f}")

# ---- architecture (D) where it is available (exactly 3 singular points) ----
print()
print("="*104)
print("Architecture (D) (Kodaira lambda-map, r*=0.46) -- available ONLY when the module has a")
print("single finite nonzero branch point.  On this host that is true for DBL alone (branch pt s2t2);")
print("W has three (s1t2,s2t1,s2t2) and COND_Z has two, so (D) is NOT available for them.")
for nm, m, cols, sl, c in (("{1,DBL,thDBL}", 3, [(1,6)]*4, [0,-26,-26], abs(s2t2)),):
    lpD = math.log(4*c) + 0.60977; bcD = lpD + 0.98175
    rD  = adelic(m, cols, [0]*m, {2: sl}, lpD, bcD)
    rD0 = adelic(m, cols, [0]*m, {}, lpD, bcD)
    print(f"  {nm:<22s} log|phi'| = {lpD:+9.4f}  tau = {rD['tau']:8.4f}  gamma2 = {rD['gamma']:+8.4f}"
          f"  entry = {rD['entry']:+8.4f}  (arch-only {rD0['entry']:+8.4f})  margin = {rD['margin']:+9.3f}")

# ---- POSITIVITY_PROGRAM adjacent pair (j0,j0+1), j0 = round(0.30 m), m = 3n ----
print()
print("="*104)
print("ADJACENT-PAIR HADAMARD HOST  (POSITIVITY_PROGRAM.md sec.4.3 fixed rule j0=round(0.30 m))")
print("  measured to n<=24 (14_pairs.gp): W_pair has LCM type [1..6n]^1 (sigma_m = 6, HALF the")
print("  Zudilin x Nesterenko value), 2-adic slope 0 (v_2(w_n) = O(1)), growth log|w_n|/n rising")
print("  through +0.019 at n=24; extrapolated limit ~ +0.26 (pre-asymptotic, uncertain in [0.02,0.40]).")
for g in (0.02, 0.26, 0.40):
    d = math.exp(-g)
    for lab, lp in (("Koebe ceiling log(4d)", math.log(4*d)), ("Landau ceiling log(16d)", math.log(16*d))):
        r = adelic(3, [(1,6)], [0]*3, {2:[0,0,0]}, lp, lp)
        print(f"  growth {g:+.2f} -> d = {d:.4f}   {lab:<24s} log|phi'| = {lp:+8.4f}  tau = {r['tau']:7.4f}"
              f"  entry <= {r['entry']:+8.4f}   margin <= {r['margin']:+9.3f}")

# ---- single-row baselines (un-Hadamarded), same measurement ---------------
print()
print("="*104)
print("BASELINES: the two source rows on their OWN hosts, m=3 = {1,H,thH}, H the conditional form")
print("  (K) univalent = Koebe 4|c|;  (D) Kodaira r*=0.46: log|phi'| = log(4|c|)+0.60977, BC = +0.98175")
PHI = (1+5**0.5)/2
BASE = [
  ("Zudilin, index m  (b=2,k=2, s2=-4)",      3, [(1,2),(1,2)], [0,-4,-4],   PHI**5),
  ("Zudilin 3-section (b=6,k=2, s2=-12)",     3, [(1,6),(1,6)], [0,-12,-12], PHI**15),
  ("Nesterenko (4,7)  (b=6,k=2, s2=-14)",     3, [(1,6),(1,6)], [0,-14,-14], 1496.5462924838050),
]
print(f"{'row':<38s} {'|c|':>12} {'tau':>8} {'gamma2':>9} {'(K)entry':>9} {'(K)marg':>9} {'(D)entry':>9} {'(D)marg':>9} {'(D)defic':>9}")
for nm, m, cols, sl, c in BASE:
    lpK = math.log(4*c)
    rK = adelic(m, cols, [0]*m, {2: sl}, lpK, lpK)
    lpD = lpK + 0.60977; bcD = lpD + 0.98175
    rD = adelic(m, cols, [0]*m, {2: sl}, lpD, bcD)
    print(f"{nm:<38s} {c:>12.5f} {rK['tau']:>8.4f} {rK['gamma']:>+9.4f} {rK['entry']:>+9.4f} "
          f"{rK['margin']:>+9.3f} {rD['entry']:>+9.4f} {rD['margin']:>+9.3f} {rD['margin']/(m-1):>+9.3f}")

print()
print("="*104)
print("COMPARISON with the modular hosts, in the per-function units of CDT_NONCONGRUENCE.md")
print("  (deficit = margin/(m-1) in (K)/(D); margin/13 in (S))")
for nm, m, cols, sl, lp, bc in INV[:1] + INV[3:5]:
    r = adelic(m, cols, [0]*m, {2: sl}, lp, bc)
    print(f"  {nm:<32s} margin {r['margin']:+9.3f}   deficit {r['margin']/(m-1):+9.3f}")
print("  Catalan level 8/16 (S)+adelic ceiling                      deficit   -0.613")
print("  Catalan level 8   (S)+measured slopes, ceiling (ADELIC sec.7.3)  margin -22.42 -> deficit  -1.725")
print("  Beukers row (D)                                            deficit   +0.978  (a proof)")
