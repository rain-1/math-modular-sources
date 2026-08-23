#!/usr/bin/env python3
"""18_padic2.py -- the 2-adic limit of b_n/a_n along n = 2^s."""
from fractions import Fraction as Fr
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
def vp(x, p):
    if x == 0: return None
    n, d = x.numerator, x.denominator
    v = 0
    while n % p == 0: n //= p; v += 1
    while d % p == 0: d //= p; v -= 1
    return v
def resid(x, p, k):
    """unit part of x mod p^k (x = p^v * u)"""
    if x == 0: return 0
    v = vp(x, p)
    y = x / Fr(p)**v
    n, d = y.numerator, y.denominator
    m = p**k
    return (n % m) * pow(d % m, -1, m) % m

N = 5000
A, B = seqs(N)
print("s   n=2^s   v2(r_n)   r_n mod 2^14")
prev = None
for s in range(0, 13):
    n = 2**s
    if n >= N-1: break
    r = B[n]/A[n]
    print("%2d  %6d   %4s      %s" % (s, n, vp(r,2), bin(resid(r,2,14))))
print()
print("differences v2(r_{2^{s+1}} - r_{2^s}):")
for s in range(0, 12):
    n0, n1 = 2**s, 2**(s+1)
    if n1 >= N-1: break
    print("  s=%2d : %s" % (s, vp(B[n1]/A[n1] - B[n0]/A[n0], 2)))
# also along n = 3*2^s and 5*2^s to test independence of the limit
for a0 in [3,5,7]:
    print("\n a0=%d : v2(r_{a0 2^{s+1}} - r_{a0 2^s}) and r mod 2^10" % a0)
    s = 0
    while a0*2**(s+1) < N-1:
        n0, n1 = a0*2**s, a0*2**(s+1)
        print("   s=%2d  v2diff=%s   r_%d mod 2^10 = %s" % (s, vp(B[n1]/A[n1]-B[n0]/A[n0],2), n1, resid(B[n1]/A[n1],2,10)))
        s += 1
