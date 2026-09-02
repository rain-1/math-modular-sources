"""11: confirm no lower-order relation exists with higher degree (guards against
declaring order 3 when an order-2 relation of degree>8 exists)."""
import sys, time
from libx import *
from lib2v import (z3_D2, z3_D3, z3_D4, E_D1, s10_D1, s10_D2)

CASES = [("z3_D2", z3_D2), ("z3_D3", z3_D3), ("z3_D4", z3_D4),
         ("E_D1", E_D1), ("s10_D1", s10_D1), ("s10_D2", s10_D2)]
A = B = 30
p = PRIMES[0]
for name, f in CASES:
    for axis, lab in ((0, "a"), (1, "b")):
        tab = table_modp(f, A, B, p, axis, 2)
        hit = None
        for d in range(0, 13):
            mons = [(s, t) for s in range(d+1) for t in range(d+1) if s+t <= d]
            ncols = 3*len(mons)
            M, mons, pts = build_rows(tab, 2, d, axis, A, B, p, 1, 1)
            if M.shape[0] < ncols + 25:
                print("   (%s %s d=%d: not enough points %d < %d)" % (name, lab, d, M.shape[0], ncols+25))
                continue
            ns, piv, free = nullspace_modp_np(M, p)
            if ns:
                hit = (d, len(ns))
                break
        print("%-7s %s-rec order-2 relation up to degree 12: %s"
              % (name, lab, "NONE" if hit is None else "FOUND at d=%d dim=%d" % hit))
        sys.stdout.flush()
