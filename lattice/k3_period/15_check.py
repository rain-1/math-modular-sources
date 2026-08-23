#!/usr/bin/env python3
"""15_check.py -- (i) independent recurrence-side confirmation of the fold
constant, (ii) p-adic slopes v_p(b_n/a_n) at the bad primes 2,3."""
from fractions import Fraction as Fr
from mpmath import mp, mpf, mpc, sqrt, nstr

def seqs(N):
    P = lambda n: 11*n*n + 11*n + 4
    Q = lambda n: 37*n*n + 3
    R = lambda n: 27*n*n - 27*n + 6
    a = [Fr(0), Fr(0), Fr(1)]
    for n in range(0, N):
        a.append((P(n)*a[n+2] - Q(n)*a[n+1] + R(n)*a[n])/Fr((n+1)**2))
    b = [Fr(0), Fr(0), Fr(0), Fr(1)]
    for n in range(1, N):
        b.append((P(n)*b[n+2] - Q(n)*b[n+1] + R(n)*b[n])/Fr((n+1)**2))
    return a[2:], b[2:]

N = 600
A, B = seqs(N)
print("a_n:", [int(x) for x in A[:10]])
print("b_n:", [str(x) for x in B[:8]])
print("integrality of a_n up to", N, ":", all(x.denominator == 1 for x in A))

# ---- (i) the complex direction ------------------------------------------
mp.dps = 60
Xc = mpf('0.52874945559576683377595806385986587134115392973305253672438')
Yc = mpf('0.37388232560046202502756426769073809774397550237286934171496')
mu = mpc(5, 1)*0  + mpc(5, float(sqrt(2)))
mu = mpc(5, 0) + mpc(0, 1)*sqrt(2)
print()
print("mu = 5 + i sqrt2 =", nstr(mu, 25))
for n in [100, 200, 300, 400, 500, 590]:
    an  = mpf(A[n].numerator)/mpf(A[n].denominator)
    an1 = mpf(A[n+1].numerator)/mpf(A[n+1].denominator)
    rn  = mpf((B[n]-0).numerator)/mpf(B[n].denominator) - Xc*an
    rn1 = mpf(B[n+1].numerator)/mpf(B[n+1].denominator) - Xc*an1
    z  = an  + mpc(0,1)*(rn /Yc)
    z1 = an1 + mpc(0,1)*(rn1/Yc)
    print("n=%4d   z_{n+1}/z_n - mu = %s" % (n, nstr(abs(z1/z - mu), 8)))
    z  = an  - mpc(0,1)*(rn /Yc)
    z1 = an1 - mpc(0,1)*(rn1/Yc)
    print("           (other sign)   = %s" % nstr(abs(z1/z - mu), 8))

# ---- (ii) p-adic slopes --------------------------------------------------
def vp(x, p):
    if x == 0: return None
    n, d = x.numerator, x.denominator
    v = 0
    while n % p == 0: n //= p; v += 1
    while d % p == 0: d //= p; v -= 1
    return v

print()
for p in [2, 3, 5, 7]:
    print("p =", p)
    for n in [10, 20, 40, 60, 80, 120, 160, 200, 300, 400, 500]:
        va = vp(A[n], p); vb = vp(B[n], p)
        print("   n=%4d  v(a_n)=%4s  v(b_n)=%5s  v(b_n/a_n)=%5s   ratio/n=%s"
              % (n, va, vb, (None if (va is None or vb is None) else vb-va),
                 (None if (va is None or vb is None) else round((vb-va)/n, 4))))
