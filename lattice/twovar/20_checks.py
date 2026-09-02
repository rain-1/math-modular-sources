"""20: (a) growth of the excess for the z3_D2 / E_D1 companions (proof that they
are NOT of lcm type); (b) ray behaviour for lambda>1 in E_D1 and z3_D4;
(c) shape of the small solutions."""
import pickle
from fractions import Fraction as F
from math import gcd
from libx import lcmrange
import lib2v
import mpmath as mp
mp.mp.dps = 40

print("(a) excess of the companion over lcm(1..a)^4 lcm(1..b)^4 on nested grids")
for name, x in [("z3_D2", (1, 0, 0)), ("E_D1", (1, 0, 0)), ("z3_D3", (1, 0, 0)),
                ("z3_D1", (0, 1))]:
    Dd = pickle.load(open("sol_%s.pkl" % name, "rb"))
    Ts, NBIG = Dd["Ts"], Dd["NBIG"]
    LCM = [lcmrange(n) for n in range(NBIG+2)]
    row = []
    for N in (8, 12, 16, 20, 24):
        acc = 1
        for a in range(2, N+1):
            for b in range(2, N+1):
                den = sum(F(xi)*T[(a, b)] for xi, T in zip(x, Ts) if xi).denominator
                X = den // gcd(den, LCM[a]**4 * LCM[b]**4)
                if X > 1:
                    acc = acc*X//gcd(acc, X)
        row.append((N, acc))
    print("   %-7s x=%s : %s" % (name, str(x),
          ", ".join("N=%d: excess has %d digits (%s)" % (N, len(str(e)), str(e)[:18]) for N, e in row)))

print()
print("(b) ray behaviour t/c along a=lambda*b for lambda>1")
for name, x, lams in [("E_D1", (1, 0, 0), [F(1), F(2), F(3)]),
                      ("z3_D4", (0, 1, 0), [F(1), F(2), F(3), F(4), F(5)])]:
    Dd = pickle.load(open("sol_%s.pkl" % name, "rb"))
    Ts, NBIG = Dd["Ts"], Dd["NBIG"]
    cfun = getattr(lib2v, name)
    print("   %s  x=%s" % (name, str(x)))
    for lam in lams:
        vals = []
        for b in range(4, NBIG+1, max(1, NBIG//8)):
            a = lam*b
            if a.denominator != 1 or int(a) > NBIG:
                continue
            a = int(a)
            q = sum(F(xi)*T[(a, b)] for xi, T in zip(x, Ts) if xi)/F(cfun(a, b))
            vals.append((a, b, mp.nstr(mp.mpf(q.numerator)/mp.mpf(q.denominator), 12)))
        print("      lam=%-3s %s" % (lam, "  ".join("(%d,%d):%s" % v for v in vals[-5:])))

print()
print("(c) shape of the small (limit 0) solutions")
for name, x in [("z3_D2", (0, 1, 0)), ("z3_D3", (0, 1, 0)), ("z3_D4", (1, 0, 0)),
                ("z2_D3", (1, 0)), ("E_D1", (2, 1, 0)), ("s10_D1", (3, 2, 0)),
                ("s10_D2", (2, 1, 0))]:
    Dd = pickle.load(open("sol_%s.pkl" % name, "rb"))
    Ts = Dd["Ts"]
    def v(a, b):
        return sum(F(xi)*T[(a, b)] for xi, T in zip(x, Ts) if xi)
    print("   %-7s x=%s  values 0<=a,b<=5:" % (name, str(x)))
    for a in range(6):
        print("        a=%d : %s" % (a, [str(v(a, b)) for b in range(6)]))
    supp = [(a, b) for a in range(9) for b in range(9) if v(a, b) != 0]
    print("        support (a,b<=8): %s%s" % (supp[:14], "..." if len(supp) > 14 else ""))
