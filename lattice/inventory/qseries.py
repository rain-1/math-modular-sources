"""Exact power-series utilities over Q, and the descent x -> y = x^2/(x-s)."""
from fractions import Fraction as F
from math import gcd

def mul(a, b, P):
    r = [F(0)]*P
    for i, ai in enumerate(a[:P]):
        if ai:
            for j, bj in enumerate(b[:P-i]):
                if bj: r[i+j] += ai*bj
    return r

def compose(f, g, P):
    """f(g(x)); requires g[0]==0."""
    assert not g[0]
    r = [F(0)]*P; pw = [F(0)]*P; pw[0] = F(1)
    for n in range(P):
        if n < len(f) and f[n]:
            for i in range(P): r[i] += f[n]*pw[i]
        if n < P-1: pw = mul(pw, g, P)
    return r

def w_series(s, P):
    """w(x) = s x/(x-s) = -x/(1-x/s) = -sum_{n>=1} x^n/s^{n-1}."""
    r = [F(0)]*P
    for n in range(1, P): r[n] = -F(1)/F(s)**(n-1)
    return r

def y_series(s, P):
    """y = x + w(x) = x^2/(x-s) = -sum_{n>=2} x^n / s^{n-1}."""
    r = [F(0)]*P
    for n in range(2, P): r[n] = -F(1)/F(s)**(n-1)
    return r

def to_y(gx, s, Py):
    """Given a symmetric g(x) (x-series), return its y-expansion, Py terms.
    y starts at x^2 so we peel off coefficients of x^{2n}."""
    Px = len(gx)
    Y = y_series(s, Px)
    Yp = [F(0)]*Px; Yp[0] = F(1)
    cur = list(gx)
    out = [F(0)]*Py
    for n in range(Py):
        if 2*n >= Px: break
        out[n] = cur[2*n]/Yp[2*n] if Yp[2*n] else (F(0) if not cur[2*n] else None)
        if out[n] is None: raise RuntimeError('non-symmetric?')
        if out[n]:
            cur = [c - out[n]*t for c, t in zip(cur, Yp)]
        Yp = mul(Yp, Y, Px)
    return out

# ---------------------------------------------------------------- lcm's
_L = [1]
def lcm_upto(n):
    while len(_L) <= n:
        k = len(_L); g = _L[-1]; _L.append(g*k//gcd(g, k))
    return _L[n]
