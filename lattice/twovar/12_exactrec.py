"""12: exact (over Q, normalised primitive integer) two-variable recurrences."""
import sys, json, time
import sympy as sp
from libx import *
from lib2v import (z3_D1, z3_D2, z3_D3, z3_D4, z2_D1, z2_D2, z2_D3,
                   E_D1, s10_D1, s10_D2)

ORD = {  # (order_a, deg_a, order_b, deg_b) from 10/11
    "z3_D1": (2, 3, 2, 3), "z3_D2": (3, 6, 3, 6), "z3_D3": (3, 6, 3, 3),
    "z3_D4": (3, 6, 3, 3), "z2_D1": (2, 2, 2, 2), "z2_D2": (2, 2, 2, 2),
    "z2_D3": (2, 2, 2, 2), "E_D1": (3, 2, 3, 2),
    "s10_D1": (3, 6, 3, 6), "s10_D2": (3, 6, 3, 3),
}
FUN = dict(z3_D1=z3_D1, z3_D2=z3_D2, z3_D3=z3_D3, z3_D4=z3_D4, z2_D1=z2_D1,
           z2_D2=z2_D2, z2_D3=z2_D3, E_D1=E_D1, s10_D1=s10_D1, s10_D2=s10_D2)
which = sys.argv[1:] if len(sys.argv) > 1 else list(ORD)
A = B = 30
res = {}
try:
    res = json.load(open("recs.json"))
except Exception:
    pass
av, bv = sp.symbols('a b')
for name in which:
    f = FUN[name]
    Ia, da, Ib, db = ORD[name]
    res.setdefault(name, {})
    for axis, lab, (I, d) in ((0, "a", (Ia, da)), (1, "b", (Ib, db))):
        t0 = time.time()
        ints, mons, k = exact_rec(f, I, d, axis, A, B, nprimes=5, amin=1, bmin=1)
        if ints is None:
            print("%-7s %s-rec: FAILED (code %s)" % (name, lab, k)); sys.stdout.flush(); continue
        bad = verify_rec(f, ints, mons, I, axis, A, B, amin=0, bmin=0)
        nm = len(mons)
        polys = [poly_sympy(ints[i*nm:(i+1)*nm], mons) for i in range(I+1)]
        res[name][lab] = dict(I=I, d=d, mons=mons, coefs=ints)
        print("=" * 78)
        print("%s  %s-recurrence: order %d, coefficient degree %d   (%.1fs)"
              % (name, lab, I, d, time.time()-t0))
        print("   verified on 0<=a,b<=%d except %d points: %s"
              % (A, len(bad), bad[:8]))
        print("   LEADING  P_%d = %s" % (I, sp.factor(polys[I])))
        print("   TRAILING P_0 = %s" % sp.factor(polys[0]))
        other = bv if axis == 0 else av
        print("   leading depends on the OTHER variable (%s): %s"
              % (other, polys[I].has(other)))
        for i in range(1, I):
            print("   P_%d = %s" % (i, sp.factor(polys[i])))
        sys.stdout.flush()
json.dump(res, open("recs.json", "w"))
print("saved", list(res))
