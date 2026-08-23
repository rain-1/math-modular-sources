"""Optimise tau over the inventory on the EMN host and report the best entry."""
import sys, os, math, itertools
sys.path.insert(0, '/home/ubuntu/code/math-modular-sources/lattice/cdt_finder')
from fractions import Fraction as F
from cdt_bound import tau_flat, tau_sharp

CEIL = math.log(8.0)      # width law: log(16*|t2|) with t2=1/2 in Z=z/2

def entry_for(npoly, c, cond_layers=(1,2), maxe=None):
    """npoly lcm-free polylog rows Li_1..Li_npoly (e=1..npoly) plus the constant,
       plus c conditional rows of layers cond_layers and e=1."""
    e = [0] + list(range(1, npoly+1)) + [1]*c
    if maxe is not None:
        e = [min(x, maxe) for x in e]
    m = len(e)
    u = m - c
    cols = [(u, F(b)) for b in cond_layers]
    sm, tf = tau_flat(m, cols)
    ts, xi = tau_sharp(m, e)
    tau = float(tf) + ts
    return m, float(sm), float(tf), ts, tau, CEIL - tau

best = None
print("EMN host (Z=z/2), ceiling log 8 = %.4f ; conditional row n[1..n][1..2n]" % CEIL)
print(" npoly  c    m   sigma_m  tau^flat  tau^sharp    tau     ENTRY")
for c in (1,2,3,4,6,8):
    for npoly in range(1, 26):
        m, sm, tf, ts, tau, ent = entry_for(npoly, c)
        if best is None or ent > best[-1]:
            best = (npoly, c, m, sm, tf, ts, tau, ent)
        if npoly in (2,3,4,5,6,8,10,14,20) :
            print(f" {npoly:4d} {c:3d} {m:4d}  {sm:7.3f}  {tf:8.4f} {ts:9.4f} {tau:8.4f}  {ent:+8.4f}")
print("\nBEST:", f"npoly={best[0]} c={best[1]} m={best[2]} sigma_m={best[3]:.3f} "
      f"tau^f={best[4]:.4f} tau^#={best[5]:.4f} tau={best[6]:.4f} ENTRY={best[7]:+.4f}")

# same optimisation with a hypothetical conditional row of weight k (one layer)
print("\nsensitivity in the conditional weight k (best over npoly,c):")
for k in (2.0, 2.5, 2.95, 3.0, 3.64):
    bb = None
    for c in (1,2,3,4):
        for npoly in range(1, 26):
            m, sm, tf, ts, tau, ent = entry_for(npoly, c, cond_layers=(k,))
            if bb is None or ent > bb[-1]: bb = (npoly, c, m, tau, ent)
    print(f"   k={k:5.2f}:  best entry={bb[-1]:+.4f}  (npoly={bb[0]}, c={bb[1]}, m={bb[2]}, tau={bb[3]:.4f})")

# calibration: CDT's own host, conditional n[1..n]^2, ceiling log 16
print("\ncalibration -- CDT's own P^1-{0,1,oo} host, ceiling log16, cond [1..n]^2:")
C2 = math.log(16.0)
bb = None
for c in (1,2,3,4,7):
    for npoly in range(1, 26):
        e = [0]+list(range(1,npoly+1))+[1]*c
        m = len(e); u = m-c
        sm, tf = tau_flat(m, [(u, F(1)), (u, F(1))])
        ts, _ = tau_sharp(m, e)
        tau = float(tf)+ts
        if bb is None or C2-tau > bb[-1]: bb = (npoly, c, m, tau, C2-tau)
print(f"   best entry={bb[-1]:+.4f} (npoly={bb[0]}, c={bb[1]}, m={bb[2]}, tau={bb[3]:.4f});  CDT report tau=2.2756, entry +0.497")
