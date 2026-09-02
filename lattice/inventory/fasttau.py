"""Fast float version of I_u^v(w) and tau^sharp (cross-checked against tau.py)."""
import math
from functools import lru_cache

def I_uvw(u, v, w):
    u, v, w = float(u), float(v), float(w)
    mw = max(1.0, w)
    a, b = min(u, 1.0), 1.0
    t1 = 0.0
    if b > a:
        lo = max(a, w)
        if b > lo: t1 = 0.5*(b-lo)**2
    A = max(u, 1.0)
    t2 = 0.0; h = 0; H = 0.0
    while 1.0 + h*mw < v:
        c0, c1 = max(A, 1.0+h*mw), min(v, 1.0+(h+1)*mw)
        if c1 > c0: t2 += H*(c1-c0)
        h += 1; H += 1.0/h
    t3 = 0.0; s = max(0.0, w-1.0); k = 1
    while k*mw - s < v:
        c0, c1 = max(A, k*mw-s), min(v, (k+1)*mw-s)
        if c1 > c0:
            lo = max(c0, k*w)
            if c1 > lo: t3 += (c1*c1-lo*lo)/(2.0*k) - w*(c1-lo)
        k += 1
    return t1+t2+t3

@lru_cache(maxsize=None)
def tau_sharp_f(m, se, me):
    if se == 0 or me == 0: return 0.0, 0.0
    f = lambda xi: (2.0/m**2)*(xi*se + me*I_uvw(xi, m, xi))
    hi = min(float(m), 30.0)
    N = 60; best = None; xb = 0.0
    for i in range(N+1):
        xi = hi*i/N; v = f(xi)
        if best is None or v < best: best, xb = v, xi
    step = hi/N
    for _ in range(28):
        step /= 2
        for xi in (xb-step, xb+step):
            if 0 <= xi <= m:
                v = f(xi)
                if v < best: best, xb = v, xi
    return best, xb
