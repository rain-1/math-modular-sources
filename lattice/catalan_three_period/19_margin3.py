"""Adelic margins on the level-16 Catalan host with BOTH slope profiles:

  (P1) "assumed"  -- the profile used in CATALAN_TWO_CLASSES.md §6 and
       ADELIC_HOLONOMY.md §4.2: pure 2, doubly-small 1, conditional 0 in y,
       the 0 coming from the "max{0,.}" integrality floor.
  (P2) "measured" -- 18run.log: in y = 4x^2/(4x+1) the measured slopes are
       pure Li_j        +2
       inner  dbl-small +1
       outer  dbl-small  0
       host A, conditional -2      <-- the floor does NOT hold: y is not an
       integral coordinate (x(v)=(v^2+v sqrt(1+v^2))/2 carries 2-power
       denominators), so Sym^+ of an integral x-series need not be integral.

Only relative slopes matter: the bound is exactly covariant under y -> y/4,
which adds +2 to every slope and subtracts log 4 from log|phi'(0)| and from BC.
Both normalisations are computed below and agree.
"""
import math, sys
sys.path.insert(0, '/home/ubuntu/code/math-modular-sources/lattice/adelic_holonomy')
from adelic_bound import adelic, fmt

L4 = math.log(4)
ceil_ = math.log(256*0.25); real_ = ceil_ + math.log(0.6292232680)
rig_  = math.log(4.376879); BC = 11.845 + math.log(0.25)

def eprof(n2, n1, nd0, nc):
    """CDT-proportional integration profile: 3 of every 7 members carry n^1."""
    def k(n): return 3 if n >= 7 else max(0, n-4)
    ne = min(3, n2) + k(n1) + k(nd0) + k(nc)
    m = n2+n1+nd0+nc
    return [1]*ne + [0]*(m-ne)

def run(tag, n2, n1, nd0, nc, lp, profile, show=True):
    m = n2+n1+nd0+nc
    cols = [(1 if n2 >= 1 else 0, 2), (min(3, n2), 2)]
    sl = [profile['pure']]*n2 + [profile['din']]*n1 + [profile['dout']]*nd0 + [profile['cond']]*nc
    r = adelic(m, cols, eprof(n2, n1, nd0, nc), {2: sl}, lp, BC)
    if show:
        print(f"  {tag:<44s} m={m:3d} (n2,n1,nd0,nc)=({n2},{n1},{nd0},{nc})  "
              f"tf={r['tau_flat']:.4f} ts={r['tau_sharp']:.4f} tau={r['tau']:.4f} "
              f"g2={r['gamma']:+.4f} entry={r['entry']:+.4f} margin={r['margin']:+8.3f}")
    return r

P1 = dict(pure=2, din=1, dout=0, cond=0)     # assumed in the earlier documents
P2 = dict(pure=2, din=1, dout=0, cond=-2)    # measured, 18run.log

print(f"BC={BC:.6f}  ceiling={ceil_:.6f}  realised={real_:.6f}  rigorous={rig_:.6f}\n")

# covariance check: shift every slope by +2 and drop log|phi'| and BC by log 4
sl = [4]*7+[3]*7+[2]*7+[0]*7
r_a = adelic(28, [(1,2),(3,2)], eprof(7,7,7,7), {2:[s-2 for s in sl]}, ceil_, BC)
r_b = adelic(28, [(1,2),(3,2)], eprof(7,7,7,7), {2:sl}, ceil_-L4, BC-L4)
print(f"scale covariance check: margin {r_a['margin']:.9f} vs {r_b['margin']:.9f} "
      f"(diff {r_a['margin']-r_b['margin']:.2e})\n")

for name, P in (("P1 assumed (pure 2, din 1, dout 0, cond 0)", P1),
                ("P2 MEASURED (pure 2, din 1, dout 0, cond -2)", P2)):
    print(f"================ {name} ================")
    for lp, lab in ((ceil_, "ceiling"), (real_, "realised (transported)"), (rig_, "rigorous 3-punctured")):
        print(f" [{lab}]")
        run("previous: 7 pure + 7 conditional", 7,0,0,7, lp, P)
        run("previous: + 1 doubly-small",       7,1,0,7, lp, P)
        run("full (a)/(b): 7 pure + 3 orbits",  7,7,7,7, lp, P)
    print(" [ceiling, optimised over sub-inventories]")
    best = bestu = None
    for n2 in range(8):
        for n1 in range(8):
            for nd0 in range(8):
                for nc in range(8):
                    if n2+n1+nd0+nc < 2: continue
                    r = run("", n2, n1, nd0, nc, ceil_, P, show=False)
                    t = (r['margin'], n2, n1, nd0, nc)
                    if nc >= 1 and (best is None or t[0] > best[0]): best = t
                    if nc == 0 and (bestu is None or t[0] > bestu[0]): bestu = t
    run("BEST containing a conditional function", best[1], best[2], best[3], best[4], ceil_, P)
    run("BEST unconditional (must be <= 0)",      bestu[1], bestu[2], bestu[3], bestu[4], ceil_, P)
    print()
