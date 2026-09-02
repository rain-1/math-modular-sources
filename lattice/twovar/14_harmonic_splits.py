"""14: explicit harmonic-sum splittings.

(a) The arithmetic of the SPLIT period is not the obstruction: the truncated
    MZV  Z(a,b) = sum_{m<=a} sum_{k<m, k<=b} 1/(m^2 k)  has denominators
    dividing [1..a]^2 [1..b] (split rates (2,1)) and tends to zeta(2,1)=zeta(3).
(b) But no such splitting can be a companion: every solution of the two-variable
    system is in span{c,s}, and s has max-type denominators.  We check that the
    natural harmonic-split candidates fail the a-recurrence.
"""
from lib2v import *
from fractions import Fraction as F
from math import comb
import mpmath as mp
mp.mp.dps = 30

print("=" * 78)
print("(a) the split truncated MZV Z(a,b) = sum_{a>=m>k>=1, k<=b} 1/(m^2 k)")
print("=" * 78)


def Z(a, b):
    return sum(F(1, m*m*k) for m in range(2, a+1) for k in range(1, min(m, b+1)))


print("   den(Z(a,b)) | [1..a]^2 [1..b] :",
      all(((lcmrange(a)**2)*lcmrange(b)) % Z(a, b).denominator == 0
          for a in range(1, 26) for b in range(1, 26)))
print("   den(Z(a,b)) | [1..a]^2 [1..b] but NOT [1..a]^2 in general:",
      any(lcmrange(a)**2 % Z(a, b).denominator != 0
          for a in range(1, 26) for b in range(1, 26)))
print("   limits (a=lambda b):")
for lam, (a, b) in [(1, (60, 60)), (2, (120, 60)), ('1/2', (60, 120)),
                    (5, (300, 60)), ('1/5', (60, 300))]:
    v = Z(a, b)
    print("     lambda=%-4s (a,b)=(%3d,%3d)  Z = %s   zeta(3) = %s"
          % (lam, a, b, mp.nstr(mp.mpf(v.numerator)/mp.mpf(v.denominator), 12),
             mp.nstr(mp.zeta(3), 12)))
print("   -> a split-denominator sequence with the right limit exists trivially;")
print("      the obstruction is holonomic, not arithmetic.")

print()
print("=" * 78)
print("(b) harmonic-split candidates versus the a-recurrence of the D1 system")
print("=" * 78)
c = [[F(z3_D1(a, b)) for b in range(23)] for a in range(23)]


def arec(t, a, b):
    P0, P2 = (a+1)**3, (a+2)**3
    P1 = -(2*a**3 + 9*a**2 + (15 + 4*b + 4*b*b)*a + 9 + 6*b + 6*b*b)
    return P2*t[a+2][b] + P1*t[a+1][b] + P0*t[a][b]


def H3(n):
    return sum(F(1, m**3) for m in range(1, n+1))


CAND = {
    "c * H3(a)                (rates (3,0))":
        [[c[a][b]*H3(a) for b in range(23)] for a in range(23)],
    "c * H3(b)                (rates (0,3))":
        [[c[a][b]*H3(b) for b in range(23)] for a in range(23)],
    "c * (H3(a)+H3(b))/2      (rates (3,3))":
        [[c[a][b]*(H3(a)+H3(b))/2 for b in range(23)] for a in range(23)],
    "c * Z(a,b)  (MZV split)  (rates (2,1))":
        [[c[a][b]*Z(a, b) for b in range(23)] for a in range(23)],
    "sum_k C(2k,k)^2 C(a+k,2k)C(b+k,2k) H3(k)  (rates (3,3) in min form)":
        [[sum(F(comb(2*k, k)**2*comb(a+k, 2*k)*comb(b+k, 2*k))*H3(k)
              for k in range(min(a, b)+1)) for b in range(23)] for a in range(23)],
}
for name, t in CAND.items():
    res = [arec(t, a, b) for a in range(3, 16) for b in range(3, 16)]
    nz = sum(1 for v in res if v != 0)
    print("   %-62s nonzero residuals %d/%d" % (name, nz, len(res)))
print()
print("   (the genuine companion s of the system has data (t00,t10,t01,t11)")
print("    = (-1,0,0,1); no harmonic splitting reproduces it because the")
print("    solution space is only 2-dimensional and s is forced.)")
