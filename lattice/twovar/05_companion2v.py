"""05: the two-variable companion for the zeta(3) D1 extension.

Three candidates:
  (T)  the triangular / pullback one:  d_{a,b} = sum_k g_k C(a+k,2k) C(b+k,2k),
       g_k determined by d_{n,n} = b_n  (Apery's companion);
  (R)  the "second solution in a" of the a-recurrence with d_{0,b}=0, d_{1,b}=gamma(b);
  (H)  harmonic-sum splittings of Apery's explicit formula.

For each: exact denominators, the two recurrences (homogeneous? rational
inhomogeneity?), the growth, and the ray limits d_{a,b}/c_{a,b}.
"""
from lib2v import *
from fractions import Fraction as F
from math import comb, log
import sys

NA = NB = 34


def triangular_solve(bseq, basis, N):
    g = []
    for n in range(N+1):
        s = F(0)
        for k in range(n):
            s += g[k]*basis(n, k)
        g.append((F(bseq[n]) - s)/F(basis(n, n)))
    return g


a3, b3 = apery3(NA)
g3 = triangular_solve(b3, lambda n, k: comb(n+k, 2*k)**2, NA)

# (T) the pullback companion
dT = [[sum(g3[k]*comb(a+k, 2*k)*comb(b+k, 2*k) for k in range(min(a, b)+1))
       for b in range(NB+1)] for a in range(NA+1)]
c = [[z3_D1(a, b) for b in range(NB+1)] for a in range(NA+1)]

print("=" * 78)
print("(T) pullback companion d_{a,b} = sum_k g_k C(a+k,2k)C(b+k,2k)")
print("=" * 78)
print("  diagonal check d_{n,n} == b_n :", all(dT[n][n] == b3[n] for n in range(NA+1)))
print("  symmetric :", all(dT[a][b] == dT[b][a] for a in range(NA+1) for b in range(NB+1)))


def arec(t, a, b):
    """P_2 c(a+2,b) + P_1 c(a+1,b) + P_0 c(a,b), the zeta3-D1 a-recurrence"""
    P0 = (a+1)**3
    P2 = (a+2)**3
    P1 = -(2*a**3 + 9*a**2 + (15 + 4*b + 4*b*b)*a + 9 + 6*b + 6*b*b)
    return P2*t[a+2][b] + P1*t[a+1][b] + P0*t[a][b]


print("  a-recurrence on the ROW      c: max |residual| =",
      max(abs(arec(c, a, b)) for a in range(NA-1) for b in range(NB+1)))
res = {(a, b): arec(dT, a, b) for a in range(NA-1) for b in range(NB+1)}
nz = [(k, v) for k, v in res.items() if v != 0]
print("  a-recurrence on companion dT : nonzero residuals:", len(nz), "of", len(res))
print("   sample residuals (a,b) -> value:")
for k in sorted(res)[:8]:
    print("     ", k, res[k])
for k in [(3, 3), (4, 4), (3, 5), (5, 3), (6, 2), (2, 6)]:
    print("     ", k, res[k])

print()
print("=" * 78)
print("denominators of the pullback companion dT")
print("=" * 78)
print("  test  den(d_{a,b}) | [1..min(a,b)]^3 :",
      all(lcmrange(min(a, b))**3 % dT[a][b].denominator == 0
          for a in range(NA+1) for b in range(NB+1)))
print("  test  den(d_{a,b}) | [1..a]^3  (all b) :",
      all(lcmrange(a)**3 % dT[a][b].denominator == 0
          for a in range(1, NA+1) for b in range(NB+1)))
print("  minimal r with den | [1..min(a,b)]^r , for a few (a,b):")
for (a, b) in [(5, 5), (5, 12), (12, 5), (9, 20), (20, 9), (17, 17), (13, 31), (31, 13)]:
    d = dT[a][b].denominator
    rs = [r for r in range(0, 7) if lcmrange(min(a, b))**r % d == 0]
    rs2 = [r for r in range(0, 7) if lcmrange(a)**r % d == 0]
    rs3 = [r for r in range(0, 7) if lcmrange(b)**r % d == 0]
    print("    (a,b)=(%2d,%2d)  min r for [1..min]^r = %s ; [1..a]^r = %s ; [1..b]^r = %s"
          % (a, b, rs[0] if rs else '>6', rs2[0] if rs2 else '>6',
             rs3[0] if rs3 else '>6'))

print()
print("=" * 78)
print("growth of dT and of the row c  (log|.| / (a+b))")
print("=" * 78)


def flog(q):
    if q == 0:
        return float('-inf')
    return log(abs(q.numerator)) - log(q.denominator) if q.denominator > 1 \
        else log(abs(q.numerator))


def flogint(n):
    return log(abs(n)) if n else float('-inf')


print("   (a,b)      log|c|/(a+b)   log|dT|/(a+b)   log(|dT|/|c|)")
for (a, b) in [(10, 10), (20, 20), (30, 30), (10, 20), (20, 10), (10, 30), (30, 10),
               (5, 30), (30, 5), (2, 32), (32, 2), (15, 25), (25, 15), (34, 34)]:
    if a > NA or b > NB:
        continue
    lc = flogint(c[a][b])
    ld = flog(F(dT[a][b]))
    print("   (%2d,%2d)     %8.4f       %8.4f        %8.3f"
          % (a, b, lc/(a+b), ld/(a+b), ld-lc))
