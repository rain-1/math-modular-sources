"""09: the exact denominator law of the two-variable companion s of the zeta(3)
D1 system, the speed of the ray limit, and the analytic domain of C and S.
"""
from lib2v import *
from fractions import Fraction as F
from math import comb, log
import mpmath as mp
mp.mp.dps = 40

NA = NB = 34


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


c = build(1, 1, 1, 5)
s = build(-1, 0, 0, 1)          # the companion: s/c -> zeta(3)-1
sc = build(0, 1, 1, 6)          # = c + s ;  sc/c -> zeta(3)

print("=" * 78)
print("A. exact denominator law of the two-variable companion s")
print("=" * 78)
print("  den(s_{a,b}) | [1..max(a,b)]^3 :",
      all(lcmrange(max(a, b))**3 % s[a][b].denominator == 0
          for a in range(NA+1) for b in range(NB+1)))
print("  den(s_{a,b}) | [1..max(a,b)]^2 :",
      all(lcmrange(max(a, b))**2 % s[a][b].denominator == 0
          for a in range(NA+1) for b in range(NB+1)))
print("  is [1..max]^3 / den bounded?  ratio [1..max(a,b)]^3/den for a>b:")
for (a, b) in [(6, 3), (10, 4), (14, 5), (20, 7), (26, 9), (30, 11), (34, 13)]:
    print("     (a,b)=(%2d,%2d)  [1..max]^3/den = %s"
          % (a, b, F(lcmrange(max(a, b))**3, s[a][b].denominator)))
print()
print("  minimal (r1,r2) needed uniformly on 2<=a,b<=%d:" % NA)
need1 = need2 = 0
for a in range(2, NA+1):
    for b in range(2, NB+1):
        den = s[a][b].denominator
        c1 = [r for r in range(9) if lcmrange(a)**r % den == 0]
        c2 = [r for r in range(9) if lcmrange(b)**r % den == 0]
        need1 = max(need1, c1[0] if c1 else 99)
        need2 = max(need2, c2[0] if c2 else 99)
print("    max over the grid of (min r with den|[1..a]^r) =", need1)
print("    max over the grid of (min r with den|[1..b]^r) =", need2)
print("    (99 means: no power of [1..a] alone clears the denominator, because")
print("     primes in (a,b] occur; so the only uniform split is (3,3))")
print("    (one-variable rate for this row is r = 3)")

print()
print("=" * 78)
print("B. speed of the ray limit  s_{a,b}/c_{a,b} -> zeta(3)-1")
print("=" * 78)
z = mp.zeta(3) - 1
print("   zeta(3)-1 =", mp.nstr(z, 25))
for (a, b) in [(10, 10), (20, 20), (30, 30), (34, 34), (30, 15), (15, 30),
               (30, 10), (10, 30), (30, 6), (6, 30), (34, 2), (2, 34)]:
    r = F(s[a][b], c[a][b])
    rv = mp.mpf(r.numerator)/mp.mpf(r.denominator)
    print("   (a,b)=(%2d,%2d)  s/c = %-24s  s/c - (zeta3-1) = %s"
          % (a, b, mp.nstr(rv, 18), mp.nstr(rv - z, 6)))

print()
print("=" * 78)
print("C. analytic domain: |c_{a,b}|^{1/(a+b)} and the polyradius constraint")
print("=" * 78)
print("   theory:  C = h(u)/((1-x)(1-y)),  u = xy/((1-x)^2(1-y)^2),")
print("            h = 2F1(1/2,1/2;1;16u) has radius 1/16, so the Reinhardt")
print("            domain is  16|x||y| < (1-|x|)^2 (1-|y|)^2.")
print("   diagonal:  16 r^2 < (1-r)^4  <=>  r < 3-2 sqrt2 = %.8f" % (3-2*mp.sqrt(2)))
print("   and (3-2 sqrt2)^2 = %.8f = (sqrt2-1)^4 = the one-variable Apery radius"
      % (3-2*mp.sqrt(2))**2)
print("   i.e. the two-variable polyradii obey  rho_1 rho_2 <= t_1  on the diagonal:")
print("        the product of the two radii is the one-variable radius.")
print()
print("   measured: log|c_{a,b}| / (a+b)")
for (a, b) in [(34, 34), (30, 30), (34, 17), (17, 34), (34, 8), (8, 34)]:
    print("     (a,b)=(%2d,%2d)  %.5f     [ -log(3-2sqrt2) = %.5f ]"
          % (a, b, log(c[a][b])/(a+b), -float(mp.log(3-2*mp.sqrt(2)))))
print("   log|s_{a,b}| / (a+b) (same growth, so S is analytic on the same domain)")
for (a, b) in [(34, 34), (30, 30), (34, 17), (17, 34), (34, 8), (8, 34)]:
    v = s[a][b]
    print("     (a,b)=(%2d,%2d)  %.5f"
          % (a, b, (log(abs(v.numerator))-log(v.denominator))/(a+b)))
