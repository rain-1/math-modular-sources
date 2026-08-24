#!/usr/bin/env python3
"""04_denom.py -- exact denominator array b, the integration vector e, the 2-adic
coefficient slopes, and the search for pure (unconditional, denominator-free)
functions on the six mixed-exponent four-term hosts.

  * A = sum a_n x^n   with a_n in Z   (checked to n = NMAX)
  * B = sum b_n x^n   the companion (b_0 = 0, b_1 = 1) of the SAME four-term
    recurrence; the sharp CDT denominator type is  b_n = (integer)/([1..k n]^...)
    -- we measure the least k with d_n^k b_n in Z, d_n = lcm(1..n), and the least
    pair (k,e) with n^e d_n^k b_n in Z.
  * 2-adic slopes  alpha_p = lim v_p(a_n)/n  (this is the R_p of ADELIC_HOLONOMY
    Theorem A applied to the PURE function A: |a_n|_p <= R_p^{-n}).
  * pure functions: sqrt(Rc_quad(t)) where Rc_quad = 1 - sigma t + pi t^2 is the
    reciprocal of the quadratic factor of the characteristic cubic -- an algebraic
    function branched exactly at the two order-2 orbifold points; and the
    polylogarithm module Li_j(r x) = sum r^n x^n / n^j on P^1 - {0, 1/r, oo}.
"""
import os, sys, json
from fractions import Fraction
from math import gcd, log

HERE = os.path.dirname(os.path.abspath(__file__))
NMAX = int(os.environ.get('NMAX', '160'))

ROWS = [   # label, r, a, c, d, f, C
    ("R1", 8,  16,  8,  48,  0, -128),
    ("R2", 2,  14,  8,  28,  4,    8),
    ("R3", -2,  6,  4, -32, -8,   32),
    ("R4", 8,  16,  8,  68,  8,   32),
    ("R5", 1,  17, 10,  32,  8,   16),
    ("R6", -1, 13,  8, -13, -1,   -1),
    ("R7", 4,   8,  4,  32,  8,   64),
]
RHO_P, RHO_R = Fraction(-1, 2), Fraction(0)


def coeffs(r, a, c, d, f, C):
    g = C
    pi = g // r
    sig = a - r
    b = int((1 - RHO_R) * r + (1 - RHO_P) * sig)
    e = int(-RHO_P * (2 * pi + r * sig) - RHO_R * r * sig)
    h = int(-(1 + 2 * RHO_P + RHO_R) * Fraction(g))
    return dict(a=a, b=b, c=c, d=d, e=e, f=f, g=g, h=h, j=0, r=r, sigma=sig, pi=pi)


def seq(co, N, u0, u1m1=0, u1=None, inhom=False):
    """(n+1)^2 u_{n+1} = P(n) u_n - Q(n) u_{n-1} + R(n) u_{n-2}."""
    P = lambda n: co['a']*n*n + co['b']*n + co['c']
    Q = lambda n: co['d']*n*n + co['e']*n + co['f']
    R = lambda n: co['g']*n*n + co['h']*n + co['j']
    u = {-2: Fraction(0), -1: Fraction(0), 0: Fraction(u0)}
    if u1 is not None:
        u[1] = Fraction(u1)
        start = 1
    else:
        start = 0
    for n in range(start, N):
        u[n+1] = (P(n)*u[n] - Q(n)*u[n-1] + R(n)*u[n-2]) / Fraction((n+1)**2)
    return [u[n] for n in range(N+1)]


def lcms(N):
    out = [1]
    cur = 1
    for n in range(1, N+1):
        cur = cur*n//gcd(cur, n)
        out.append(cur)
    return out


def vp(x, p):
    if x == 0:
        return None
    n, d = x.numerator, x.denominator
    v = 0
    while n % p == 0:
        n //= p; v += 1
    while d % p == 0:
        d //= p; v -= 1
    return v


def sqrt_series(poly, N):
    """power series sqrt of 1 + poly (poly a list of coefficients, poly[0]=1)."""
    s = [Fraction(1)] + [Fraction(0)]*N
    for n in range(1, N+1):
        acc = Fraction(poly[n]) if n < len(poly) else Fraction(0)
        for k in range(1, n):
            acc -= s[k]*s[n-k]
        s[n] = acc/2
    return s


def main():
    res = []
    L = lcms(NMAX+2)
    for lab, r, a, c, d, f, C in ROWS:
        co = coeffs(r, a, c, d, f, C)
        A = seq(co, NMAX, 1)
        B = seq(co, NMAX, 0, u1=1)
        aint = all(x.denominator == 1 for x in A)
        # sharp k: least k with d_n^k b_n in Z for all n<=NMAX
        kmin = 0
        for n in range(1, NMAX+1):
            den = B[n].denominator
            k = 0
            while den > 1:
                gg = gcd(den, L[n])
                if gg == 1:
                    break
                den //= gg
                k += 1
            if den > 1:
                kmin = 99
                break
            kmin = max(kmin, k)
        # is n^e d_n^{k} enough with smaller k?
        best = None
        for k in range(0, 5):
            for e in range(0, 4):
                ok = all((Fraction(n)**e * L[n]**k * B[n]).denominator == 1
                         for n in range(1, NMAX+1))
                if ok:
                    best = (k, e); break
            if best:
                break
        # 2-adic slope of a_n and of b_n
        sl = {}
        for p in (2, 3, 5, 7):
            vs = [vp(A[n], p) for n in range(NMAX//2, NMAX+1)]
            vs = [v for v in vs if v is not None]
            if vs:
                ns = list(range(NMAX//2, NMAX+1))[:len(vs)]
                # least squares slope
                m = len(vs); sx = sum(ns); sy = sum(vs)
                sxx = sum(x*x for x in ns); sxy = sum(x*y for x, y in zip(ns, vs))
                sl[p] = (m*sxy - sx*sy)/(m*sxx - sx*sx)
        # pure algebraic function: sqrt of the reciprocal quadratic factor
        quad = [Fraction(1), Fraction(-co['sigma']), Fraction(co['pi'])]
        S = sqrt_series(quad, 60)
        sqint = all(x.denominator == 1 for x in S)
        # sqrt of the full symbol Rc = 1 - a t + d t^2 - g t^3
        full = [Fraction(1), Fraction(-co['a']), Fraction(co['d']), Fraction(-co['g'])]
        SF = sqrt_series(full, 60)
        sqfull = all(x.denominator == 1 for x in SF)
        # sqrt(1 - r t) (the rational factor alone)
        SR = sqrt_series([Fraction(1), Fraction(-r)], 60)
        sqr = all(x.denominator == 1 for x in SR)
        rec = dict(label=lab, a_integral=aint, k_sharp=kmin, k_e=best,
                   slopes={str(k): float(v) for k, v in sl.items()},
                   sqrt_quad_integral=sqint, sqrt_full_integral=sqfull,
                   sqrt_lin_integral=sqr,
                   quad="1 - %d t + %d t^2" % (co['sigma'], co['pi']),
                   sqrt_quad_head=[str(x) for x in S[:8]],
                   r=r, v2_a=[vp(A[n], 2) for n in (20, 40, 80, NMAX)])
        res.append(rec)
        print("%s  r=%3d  A integral: %s   sharp k = %d   (k,e) = %s" %
              (lab, r, aint, kmin, best))
        print("     2-adic slope of a_n: %s   (other primes: %s)" %
              (sl.get(2), {k: round(v, 4) for k, v in sl.items() if k != 2}))
        print("     sqrt(%s) integral: %s     head %s" %
              (rec['quad'], sqint, rec['sqrt_quad_head'][:6]))
        print("     sqrt(full symbol) integral: %s ; sqrt(1-%dt) integral: %s"
              % (sqfull, r, sqr))
    json.dump(res, open(os.path.join(HERE, 'out', 'denom.json'), 'w'), indent=1)
    print("\nwrote out/denom.json  (NMAX = %d)" % NMAX)


if __name__ == '__main__':
    main()
