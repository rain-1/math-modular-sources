"""08: the full power-series solution space of the zeta(3) D1 two-variable system
(the pair of order-2 recurrences), with denominators and ray limits of each basis
solution.  This is the two-variable analogue of {row, companion}.

Initial data: t_{0,0}, t_{1,0}, t_{0,1}, t_{1,1}.
   a-recurrence -> t_{a,0}, t_{a,1} for all a;
   b-recurrence -> t_{a,b} for b>=2;
   consistency: the a-recurrence must still hold for b>=2.
"""
from lib2v import *
from fractions import Fraction as F
from math import comb
import mpmath as mp
mp.mp.dps = 40

NA, NB = 30, 30


def build(t00, t10, t01, t11, NA=NA, NB=NB):
    t = [[F(0)]*(NB+1) for _ in range(NA+1)]
    t[0][0], t[1][0], t[0][1], t[1][1] = map(F, (t00, t10, t01, t11))
    for b in (0, 1):
        B = b*(b+1)
        for n in range(1, NA):
            mid = 2*n**3 + 3*n**2 + (3+4*B)*n + 1 + 2*B
            t[n+1][b] = (F(mid)*t[n][b] - F(n**3)*t[n-1][b])/F((n+1)**3)
    for a in range(NA+1):
        for b in range(NB-1):
            Q1 = -(2*b**3 + 9*b**2 + (15+4*a+4*a*a)*b + 9 + 6*a + 6*a*a)
            t[a][b+2] = (-F(Q1)*t[a][b+1] - F((b+1)**3)*t[a][b])/F((b+2)**3)
    return t


def arec_ok(t):
    bad = []
    for a in range(NA-1):
        for b in range(NB+1):
            P0, P2 = (a+1)**3, (a+2)**3
            P1 = -(2*a**3 + 9*a**2 + (15+4*b+4*b*b)*a + 9 + 6*b + 6*b*b)
            if P2*t[a+2][b] + P1*t[a+1][b] + P0*t[a][b] != 0:
                bad.append((a, b))
    return bad


print("=" * 78)
print("A. which initial data (t00,t10,t01,t11) give a consistent solution?")
print("=" * 78)
basis = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
consistent = []
for v in basis:
    t = build(*v, NA=12, NB=12)
    bad = []
    for a in range(10):
        for b in range(10):
            P0, P2 = (a+1)**3, (a+2)**3
            P1 = -(2*a**3 + 9*a**2 + (15+4*b+4*b*b)*a + 9 + 6*b + 6*b*b)
            if P2*t[a+2][b] + P1*t[a+1][b] + P0*t[a][b] != 0:
                bad.append((a, b))
    print("   e_%s : a-recurrence violated at %d of 100 (a,b)" % (str(v), len(bad)))
    consistent.append((v, len(bad)))

# solve for the consistent subspace: impose the a-recurrence residuals as linear
# conditions on (t00,t10,t01,t11)
import itertools
rowsL = []
Ts = [build(*v, NA=12, NB=12) for v in basis]
for a in range(9):
    for b in range(2, 9):
        P0, P2 = (a+1)**3, (a+2)**3
        P1 = -(2*a**3 + 9*a**2 + (15+4*b+4*b*b)*a + 9 + 6*b + 6*b*b)
        rowsL.append([P2*T[a+2][b] + P1*T[a+1][b] + P0*T[a][b] for T in Ts])
import sympy as sp
M = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in r] for r in rowsL])
ns = M.nullspace()
print("   dimension of the consistent solution space:", len(ns))
for v in ns:
    v = v / max([abs(x) for x in v])
    print("     (t00,t10,t01,t11) =", [sp.nsimplify(x) for x in v])

print()
print("=" * 78)
print("B. the two basis solutions: row c and the 'symmetric companion'")
print("=" * 78)
sols = []
for v in ns:
    l = sp.lcm([sp.Rational(x).q for x in v])
    vv = [F(int(sp.Rational(x*l).p), int(sp.Rational(x*l).q)) for x in v]
    sols.append(vv)
    print("   solution with data", vv)

# the row
c = build(1, 1, 1, 5)
print("   row c matches z3_D1 :",
      all(c[a][b] == z3_D1(a, b) for a in range(NA+1) for b in range(NB+1)))


def report(t, name):
    print()
    print("   ---", name)
    print("     first rows: ", [[str(t[a][b]) for b in range(4)] for a in range(3)])
    print("     symmetric  :", all(t[a][b] == t[b][a] for a in range(NA+1)
                                   for b in range(min(NA, NB)+1)))
    print("     denominator type: minimal (r1,r2), r1+r2 minimal, over a grid")
    seen = {}
    for a in range(2, NA+1):
        for b in range(2, NB+1):
            den = t[a][b].denominator
            best = None
            for s in range(0, 9):
                for r1 in range(0, s+1):
                    r2 = s - r1
                    if (lcmrange(a)**r1 * lcmrange(b)**r2) % den == 0:
                        best = (r1, r2)
                        break
                if best:
                    break
            seen.setdefault(best, []).append((a, b))
    for k in sorted(seen, key=lambda z: (z[0]+z[1], z)):
        print("       (r1,r2)=%s : %d cells, e.g. %s" % (k, len(seen[k]), seen[k][:4]))
    print("     den | [1..max(a,b)]^3 everywhere :",
          all(lcmrange(max(a, b))**3 % t[a][b].denominator == 0
              for a in range(NA+1) for b in range(NB+1)))
    print("     den | [1..min(a,b)]^3 everywhere :",
          all(lcmrange(min(a, b))**3 % t[a][b].denominator == 0
              for a in range(NA+1) for b in range(NB+1)))
    print("     ray limits t/c :")
    for lam, (a, b) in [(1, (30, 30)), (2, (30, 15)), ('1/2', (15, 30)),
                        (3, (30, 10)), ('1/3', (10, 30)), (5, (30, 6)), ('1/5', (6, 30))]:
        r = F(t[a][b], c[a][b])
        print("       lambda=%-4s (a,b)=(%2d,%2d)  t/c = %s"
              % (lam, a, b, mp.nstr(mp.mpf(r.numerator)/mp.mpf(r.denominator), 15)))


for i, vv in enumerate(sols):
    t = build(*vv)
    report(t, "basis solution %d, data %s" % (i, vv))
print()
print("   zeta(3) =", mp.nstr(mp.zeta(3), 15))
