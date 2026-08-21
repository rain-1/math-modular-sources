"""Integer/rational q-series utilities + eta quotients + Ligozat cusp orders."""
from fractions import Fraction
from math import gcd

def divisors(n):
    return [d for d in range(1, n+1) if n % d == 0]

def phi(n):
    r = n; p = 2; m = n
    while p*p <= m:
        if m % p == 0:
            while m % p == 0: m //= p
            r -= r//p
        p += 1
    if m > 1: r -= r//m
    return r

def eta_factor(d, P):
    """prod_{n>=1}(1-q^{dn}) to O(q^P)."""
    f = [0]*P; f[0] = 1
    n = 1
    while d*n < P:
        m = d*n
        for i in range(P-1, m-1, -1):
            f[i] -= f[i-m]
        n += 1
    return f

def smul(a, b, P):
    r = [0]*P
    for i, ai in enumerate(a):
        if ai == 0: continue
        for j in range(min(len(b), P-i)):
            if b[j]: r[i+j] += ai*b[j]
    return r

def sinv(a, P):
    assert a[0] == 1
    r = [0]*P; r[0] = 1
    for n in range(1, P):
        s = 0
        for k in range(1, min(n, len(a)-1)+1):
            if a[k]: s += a[k]*r[n-k]
        r[n] = -s
    return r

def spow(a, e, P):
    if e == 0:
        r = [0]*P; r[0] = 1; return r
    neg = e < 0; e = abs(e)
    base = a[:]; res = [0]*P; res[0] = 1
    while e:
        if e & 1: res = smul(res, base, P)
        e >>= 1
        if e: base = smul(base, base, P)
    return sinv(res, P) if neg else res

def etaquo_series(D, r, P):
    """q-expansion of prod eta(d tau)^{r_d} divided by q^{ord}, i.e. the series
    part; returns (ord24, series) where ord24 = sum d*r_d (so q-power = ord24/24)."""
    res = [0]*P; res[0] = 1
    for d, rd in zip(D, r):
        if rd == 0: continue
        res = smul(res, spow(eta_factor(d, P), rd, P), P)
    return sum(d*rd for d, rd in zip(D, r)), res

def ligozat_matrix(N):
    """rows indexed by cusp denominators c|N, cols by d|N; ord_c(prod eta(d)^{r_d})."""
    D = divisors(N)
    M = []
    for c in D:
        row = [Fraction(N, 24*gcd(c*c, N)) * Fraction(gcd(c, d)**2, d) for d in D]
        M.append(row)
    return D, M

def cusp_counts(N):
    return {c: phi(gcd(c, N//c)) for c in divisors(N)}
