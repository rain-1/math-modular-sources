#!/usr/bin/env python3
"""16_padic.py -- p-adic Cauchy test for an Apery limit of b_n/a_n along
n = a*p^s (the Euler-criterion protocol of consolidation/EULER_CRITERION.md)."""
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
    if x == 0: return 10**6
    n, d = x.numerator, x.denominator
    v = 0
    while n % p == 0: n //= p; v += 1
    while d % p == 0: d //= p; v -= 1
    return v

N = 2200
A, B = seqs(N)
print("n, p:  v_p(r_{n} - r_{n/p})   with r_n = b_n/a_n")
for p in [2,3,5,7,11,13]:
    for a0 in [1,2,4]:
        ns = []
        s = 0
        while a0*p**s < N-1:
            ns.append(a0*p**s); s += 1
        row = []
        for i in range(1, len(ns)):
            n1, n0 = ns[i], ns[i-1]
            if A[n1] == 0 or A[n0] == 0: row.append('a=0'); continue
            row.append(vp(B[n1]/A[n1] - B[n0]/A[n0], p))
        print("p=%2d a=%d  n=%s  ->  %s" % (p, a0, ns, row))
