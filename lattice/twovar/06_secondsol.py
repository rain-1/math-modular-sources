"""06: the 'second solution in a' companion of the zeta(3) D1 family.

For each fixed b the row c_{.,b} satisfies the Apery-shaped recurrence

  (n+1)^3 u_{n+1} = (2n^3+3n^2+(3+4B)n+1+2B) u_n - n^3 u_{n-1},   B = b(b+1)

whose characteristic roots are both 1 (c_{a,b} is a polynomial in a of degree 2b).
d^{(R)}_{a,b} := second solution with d_{0,b}=0, d_{1,b}=1.
Denominator type: [1..a]^3 only (rates (3,0)).  Questions: (i) the fixed-b
Apery limit xi(b) = lim_a d_{a,b}/c_{a,b};  (ii) does d satisfy the b-recurrence?
"""
from lib2v import *
from fractions import Fraction as F
from math import comb, log
import mpmath as mp

mp.mp.dps = 40
NA = 4000
BMAX = 8

print("=" * 78)
print("A. fixed-b Apery limits of the second solution   (b = 0..8)")
print("=" * 78)
print("   zeta(3) = ", mp.zeta(3))
xis = {}
for b in range(BMAX+1):
    B = b*(b+1)
    c = [mp.mpf(z3_D1(0, b)), mp.mpf(z3_D1(1, b))]
    d = [mp.mpf(0), mp.mpf(1)]
    for n in range(1, NA):
        mid = 2*n**3 + 3*n**2 + (3+4*B)*n + 1 + 2*B
        c.append((mid*c[n] - n**3*c[n-1])/mp.mpf((n+1)**3))
        d.append((mid*d[n] - n**3*d[n-1])/mp.mpf((n+1)**3))
    xi = d[NA]/c[NA]
    xis[b] = xi
    print("   b=%d  B=%2d   xi(b) = %s     xi/zeta(3) = %s"
          % (b, B, mp.nstr(xi, 20), mp.nstr(xi/mp.zeta(3), 20)))

print()
print("   ratios xi(b)/xi(0):")
for b in range(BMAX+1):
    print("     b=%d  %s" % (b, mp.nstr(xis[b]/xis[0], 20)))

print()
print("=" * 78)
print("B. exact second solution; denominators; does it satisfy the b-recurrence?")
print("=" * 78)
NAE, NBE = 26, 12
d = [[F(0)]*(NBE+1) for _ in range(NAE+1)]
c = [[F(z3_D1(a, b)) for b in range(NBE+1)] for a in range(NAE+1)]
for b in range(NBE+1):
    B = b*(b+1)
    d[0][b] = F(0)
    d[1][b] = F(1)
    for n in range(1, NAE):
        mid = 2*n**3 + 3*n**2 + (3+4*B)*n + 1 + 2*B
        d[n+1][b] = (mid*d[n][b] - F(n**3)*d[n-1][b])/F((n+1)**3)

print("  denominators: den(d_{a,b}) | [1..a]^3 for all a<=%d, b<=%d : %s"
      % (NAE, NBE, all(lcmrange(a)**3 % d[a][b].denominator == 0
                       for a in range(1, NAE+1) for b in range(NBE+1))))
print("  d_{a,0} == H^(3)_a :",
      all(d[a][0] == sum(F(1, m**3) for m in range(1, a+1)) for a in range(NAE+1)))
print("  minimal r with den(d_{a,b}) | [1..a]^r  at a=20:",
      [min([r for r in range(7) if lcmrange(20)**r % d[20][b].denominator == 0] or [99])
       for b in range(NBE+1)])
print("  does den involve b at all?  den(d_{a,b})/den(d_{a,0}) for a=20:",
      [F(d[20][b].denominator, d[20][0].denominator) for b in range(min(6, NBE+1))])


def brec(t, a, b):
    """Q_2 t(a,b+2) + Q_1 t(a,b+1) + Q_0 t(a,b) (mirror of the a-recurrence)"""
    Q0 = (b+1)**3
    Q2 = (b+2)**3
    Q1 = -(2*b**3 + 9*b**2 + (15 + 4*a + 4*a*a)*b + 9 + 6*a + 6*a*a)
    return Q2*t[a][b+2] + Q1*t[a][b+1] + Q0*t[a][b]


print("  b-recurrence on the row c: max residual =",
      max(abs(brec(c, a, b)) for a in range(NAE+1) for b in range(NBE-1)))
print("  b-recurrence on d^(R): residuals (a,b) -> value")
for a in range(0, 9):
    print("     a=%d :" % a, [brec(d, a, b) for b in range(0, 5)])
