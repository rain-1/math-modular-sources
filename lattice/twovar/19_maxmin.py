"""19: MIN- / MAX-type exponents with the bounded-excess (normalisation-robust)
criterion, for the row, the small direction and a companion direction."""
import sys, pickle
from fractions import Fraction as F
from math import gcd
from libx import lcmrange
import lib2v

CASES = [("z3_D1", [(1,1),(0,1)]), ("z2_D1", [(1,1),(0,1)]), ("z2_D2", [(1,1),(0,1)]),
         ("z2_D3", [(1,1),(1,0)]),
         ("z3_D2", [(2,1,1),(0,1,0),(1,0,0)]), ("z3_D3", [(3,1,1),(0,1,0),(1,0,0)]),
         ("z3_D4", [(1,1,1),(1,0,0),(0,1,0)]), ("E_D1", [(6,2,1),(2,1,0),(1,0,0)]),
         ("s10_D1", [(1,1,1),(3,2,0),(1,0,0)]), ("s10_D2", [(1,1,1),(2,1,0),(1,0,0)])]
NG = 22
for name, xs in CASES:
    Dd = pickle.load(open("sol_%s.pkl" % name, "rb"))
    Ts, NBIG = Dd["Ts"], Dd["NBIG"]
    LCM = [lcmrange(n) for n in range(NBIG+2)]
    CELLS = [(a, b) for a in range(2, NG+1) for b in range(2, NG+1)]
    BV = {ab: tuple(T[ab] for T in Ts) for ab in CELLS}
    print("-" * 74)
    print(name)
    for x in xs:
        out = []
        for mode in ("max", "min"):
            found = None
            for r in range(0, 8):
                acc = {N: 1 for N in (10, 16, NG)}
                for (a, b) in CELLS:
                    den = sum(F(xi)*v for xi, v in zip(x, BV[(a, b)]) if xi).denominator
                    if den == 1:
                        continue
                    m = max(a, b) if mode == "max" else min(a, b)
                    X = den // gcd(den, LCM[m]**r)
                    if X > 1:
                        for N in (10, 16, NG):
                            if a <= N and b <= N:
                                acc[N] = acc[N]*X//gcd(acc[N], X)
                if acc[10] == acc[16] == acc[NG]:
                    found = (r, acc[NG])
                    break
            out.append((mode, found))
        print("   x=%-12s  %s" % (str(x), "  ".join(
            "%s-type: r=%s (bounded excess %s)" % (m, f[0] if f else ">7", f[1] if f else "-")
            for m, f in out)))
