"""07: two-variable Apery limits along rays a = lambda*b, for the candidate
companions of the zeta(3) D1 extension, and the corrected candidate
 dtilde = d^(R) + H^(3)_b * c   (fixed-b limit zeta(3) for EVERY b).
"""
from lib2v import *
from fractions import Fraction as F
from math import comb
import mpmath as mp

mp.mp.dps = 60


def dR_row(b, AMAX):
    """second solution in a for fixed b, exact"""
    B = b*(b+1)
    d = [F(0), F(1)]
    for n in range(1, AMAX):
        mid = 2*n**3 + 3*n**2 + (3+4*B)*n + 1 + 2*B
        d.append((F(mid)*d[n] - F(n**3)*d[n-1])/F((n+1)**3))
    return d


def c_row(b, AMAX):
    B = b*(b+1)
    c = [F(z3_D1(0, b)), F(z3_D1(1, b))]
    for n in range(1, AMAX):
        mid = 2*n**3 + 3*n**2 + (3+4*B)*n + 1 + 2*B
        c.append(F(mid)*c[n] - F(n**3)*c[n-1])
    return c


def H3(b):
    return sum(F(1, m**3) for m in range(1, b+1))


print("=" * 78)
print("A. ray limits of d^(R)  (rate profile (3,0)) and of dtilde = d^(R)+H3(b)c")
print("=" * 78)
print("   zeta(3) =", mp.nstr(mp.zeta(3), 25))
print()
print("  lambda   b      d^(R)/c        dtilde/c        (zeta3-H3(b))")
for lam in [F(1), F(2), F(1, 2), F(4), F(1, 4), F(10), F(1, 10)]:
    for b in [20, 40, 80, 160]:
        a = int(lam*b)
        AMAX = a+2
        if AMAX < 3:
            continue
        d = dR_row(b, AMAX)
        c = c_row(b, AMAX)
        def tomp(q):
            q = F(q)
            return mp.mpf(q.numerator)/mp.mpf(q.denominator)
        r1 = tomp(d[a]/c[a])
        dt = d[a] + H3(b)*c[a]
        r2 = tomp(dt/c[a])
        tail = mp.zeta(3) - tomp(H3(b))
        print("  %6s  %3d   %s   %s   %s"
              % (lam, b, mp.nstr(r1, 12), mp.nstr(r2, 12), mp.nstr(tail, 12)))
    print()

print("=" * 78)
print("B. does dtilde satisfy the b-recurrence?  and what are its denominators?")
print("=" * 78)
NAE, NBE = 24, 14
c = [[F(z3_D1(a, b)) for b in range(NBE+1)] for a in range(NAE+1)]
dR = [[F(0)]*(NBE+1) for _ in range(NAE+1)]
for b in range(NBE+1):
    col = dR_row(b, NAE+1)
    for a in range(NAE+1):
        dR[a][b] = col[a]
dt = [[dR[a][b] + H3(b)*c[a][b] for b in range(NBE+1)] for a in range(NAE+1)]


def brec(t, a, b):
    Q0 = (b+1)**3
    Q2 = (b+2)**3
    Q1 = -(2*b**3 + 9*b**2 + (15 + 4*a + 4*a*a)*b + 9 + 6*a + 6*a*a)
    return Q2*t[a][b+2] + Q1*t[a][b+1] + Q0*t[a][b]


def arec(t, a, b):
    P0 = (a+1)**3
    P2 = (a+2)**3
    P1 = -(2*a**3 + 9*a**2 + (15 + 4*b + 4*b*b)*a + 9 + 6*b + 6*b*b)
    return P2*t[a+2][b] + P1*t[a+1][b] + P0*t[a][b]


print("  a-recurrence residual of dtilde (should be H3(b)*0 + brec-free):")
for a in range(0, 6):
    print("    a=%d :" % a, [arec(dt, a, b) for b in range(0, 5)])
print("  b-recurrence residual of dtilde:")
for a in range(0, 6):
    print("    a=%d :" % a, [brec(dt, a, b) for b in range(0, 5)])
print()
print("  denominators of dtilde: minimal (r1,r2) with den | [1..a]^r1 [1..b]^r2")
for (a, b) in [(6, 6), (10, 4), (4, 10), (13, 13), (20, 7), (7, 20), (24, 14), (14, 24)]:
    if a > NAE or b > NBE:
        continue
    den = dt[a][b].denominator
    best = None
    for r1 in range(0, 5):
        for r2 in range(0, 5):
            if (lcmrange(a)**r1 * lcmrange(b)**r2) % den == 0:
                if best is None or r1+r2 < best[0]+best[1]:
                    best = (r1, r2)
    rmin = [r for r in range(6) if lcmrange(min(a, b))**r % den == 0]
    print("    (a,b)=(%2d,%2d)  best (r1,r2) = %s   ;  [1..min]^r : r = %s"
          % (a, b, best, rmin[0] if rmin else '>5'))
