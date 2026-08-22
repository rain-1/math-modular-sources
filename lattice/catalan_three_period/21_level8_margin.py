"""Corrected ADELIC_HOLONOMY.md §4.2 table for the LEVEL-8 Catalan host.

Measured slopes in y = 4t^2/(4t-1) (20run.log, t-series to n=300, k<=147):
    Sym+/Sym- Li_j(4t)  (pure orbit)        +2      (confirms pure_2adic.py)
    Sym+ A                (host)            -2
    Sym+ B_E, Sym+(B_E + mu A) (conditional) -2     <-- assumed 0 in the old table
    Sym+ B_0 = (1-4V2)E companion            +3     (assumed 4; unusable anyway)
in the t-coordinate: A 0, B_E 0, Li_2(4t) 2, B_0 3 -- all as assumed.
The defect is the descent, not the x-slopes: sigma(t)=t/(4t-1) and
t = (y - sqrt(y^2-y))/2 carries 2-power denominators, so Sym^+ of an integral
t-series need not be integral in y and the max{0,.} floor of §3 does not apply.
"""
import math, sys
sys.path.insert(0, '/home/ubuntu/code/math-modular-sources/lattice/adelic_holonomy')
from adelic_bound import adelic, fmt

ceil_ = math.log(256*0.25)                  # 4.158883
real_ = ceil_ + math.log(0.6292232680)      # 3.695614
BC    = 11.845 + math.log(0.25)             # 10.458706
cols  = [(1,2),(3,2)]
e14   = [0,0,1,0,0,0,0,0,0,1,1,1,1,1]

def row(tag, m, sl, lp, e=None, c=None):
    r = adelic(m, c or cols, e if e is not None else e14, ({2:sl} if sl else {}), lp, BC)
    print(f"  {tag:<46s} m={m:3d}  tau={r['tau']:.6f}  gamma2={r['gamma']:+.4f}  "
          f"entry={r['entry']:+.4f}  margin={r['margin']:+8.3f}")
    return r

print(f"BC={BC:.6f}  ceiling={ceil_:.6f}  realised={real_:.6f}\n")
print("=== ADELIC_HOLONOMY.md §4.2, m=14 = 7 pure + 7 conditional ===")
print(" [SUPERSEDED: conditional slope assumed 0]")
row("archimedean only, ceiling",            14, None,             ceil_)
row("archimedean only, realised contour",   14, None,             real_)
row("adelic (pure 2, cond 0), ceiling",     14, [2]*7+[0]*7,      ceil_)
row("adelic (pure 2, cond 0), realised",    14, [2]*7+[0]*7,      real_)
print(" [CORRECTED: conditional slope measured -2]")
row("adelic (pure 2, cond -2), ceiling",    14, [2]*7+[-2]*7,     ceil_)
row("adelic (pure 2, cond -2), realised",   14, [2]*7+[-2]*7,     real_)
print()
print("=== sub-inventory optimisation (ceiling), pure slope 2 / conditional slope -2 ===")
best = bestu = None
for n2 in range(8):
    for nc in range(8):
        m = n2+nc
        if m < 2: continue
        c = [(1 if n2 else 0,2),(min(3,n2),2)]
        ne = min(3,n2) + (3 if nc==7 else max(0,nc-4))
        e  = [1]*ne + [0]*(m-ne)
        r  = adelic(m, c, e, {2:[2]*n2+[-2]*nc}, ceil_, BC)
        t  = (r['margin'], n2, nc)
        if nc >= 1 and (best  is None or t[0] > best[0]):  best  = t
        if nc == 0 and (bestu is None or t[0] > bestu[0]): bestu = t
for tag, t in (("BEST containing a conditional function", best), ("BEST unconditional (must be <=0)", bestu)):
    m = t[1]+t[2]
    c = [(1 if t[1] else 0,2),(min(3,t[1]),2)]
    ne = min(3,t[1]) + (3 if t[2]==7 else max(0,t[2]-4))
    row(f"{tag}  (n_pure,n_cond)=({t[1]},{t[2]})", m, [2]*t[1]+[-2]*t[2], ceil_,
        e=[1]*ne+[0]*(m-ne), c=c)
print()
print("=== consistency: the unconditional pure module alone (§2.4 check 3) ===")
for m in (7, 9, 10, 11):
    c = [(1,2),(3,2)]
    ne = min(3,m); e = [1]*ne+[0]*(m-ne)
    r = adelic(m, c, e, {2:[2]*m}, ceil_, BC)
    print(f"  pure only, m={m:2d}: tau={r['tau']:.4f} gamma2={r['gamma']:+.4f} "
          f"entry={r['entry']:+.4f} bound={r['bound']:.3f} margin={r['margin']:+8.3f}")

print()
print("=== the two extra level-8 rows under BOTH slope profiles (ceiling) ===")
for n2, nc in ((7,1),(7,0)):
    m = n2+nc
    c = [(1,2),(3,2)]
    ne = min(3,n2) + max(0,nc-4); e = [1]*ne+[0]*(m-ne)
    for lab, sc in (("P1 cond 0", 0), ("P2 cond -2", -2)):
        r = adelic(m, c, e, {2:[2]*n2+[sc]*nc}, ceil_, BC)
        print(f"  (n_pure,n_cond)=({n2},{nc}) {lab:11s} m={m:2d} tf={r['tau_flat']:.4f} "
              f"ts={r['tau_sharp']:.4f} tau={r['tau']:.4f} g2={r['gamma']:+.4f} "
              f"entry={r['entry']:+.4f} margin={r['margin']:+8.3f}")
