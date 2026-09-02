"""10: minimal (order, degree) of the a- and b-recurrences, mod p, for every
decoupling.  All (a,b) sample points in 1<=a,b<=A are used and the relation is
verified on all of them."""
import sys, time
from libx import *
from lib2v import (z3_D1, z3_D2, z3_D3, z3_D4, z2_D1, z2_D2, z2_D3,
                   E_D1, s10_D1, s10_D2)

CASES = [
    ("z3_D1", z3_D1), ("z3_D2", z3_D2), ("z3_D3", z3_D3), ("z3_D4", z3_D4),
    ("z2_D1", z2_D1), ("z2_D2", z2_D2), ("z2_D3", z2_D3),
    ("E_D1", E_D1), ("s10_D1", s10_D1), ("s10_D2", s10_D2),
]
which = sys.argv[1:] if len(sys.argv) > 1 else [n for n, _ in CASES]
A = B = 28
for name, f in CASES:
    if name not in which:
        continue
    for axis, lab in ((0, "a"), (1, "b")):
        t0 = time.time()
        I, d, k, ns, mons, pts = scan_minrec(f, axis, A, B, maxI=5, maxd=8,
                                             amin=1, bmin=1)
        print("%-7s %s-rec:  order=%s  degree=%s  nullspace_dim=%s   (%.1fs, %d pts)"
              % (name, lab, I, d, k, time.time()-t0, len(pts) if pts else 0))
        sys.stdout.flush()
