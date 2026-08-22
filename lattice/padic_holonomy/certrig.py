"""FULLY A-PRIORI certificate: no sampled quantity enters the error bound.

  Lip(g) <= 2 pi * M(r) * (1 + sum_k k|c_k|),
  M(r)   := 1 + 6 sum_n n ( r^n/(1-r^n) + 5 r^{5n}/(1-r^{5n}) )   >= max_{|q|<=r}|Lam(q)|
            (triangle inequality on the Lambert series, r = max|q| on the contour)
  |RE - RE_N| <= Lip(g)/N.
"""
import os, sys as _sys
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in _sys.path: _sys.path.insert(0, HERE)
import json, math, warnings
import numpy as np
warnings.filterwarnings('ignore')
import haupt, outer
from targets import TARGETS

idx = int(_sys.argv[1]); c = np.array(json.load(open(_sys.argv[2]))['c'], float)
T = TARGETS[idx]; H, m, L, tau = T['H'], T['m'], T['L'], T['tau']
Sk = float(np.sum(np.arange(1, len(c)) * np.abs(c[1:])))


def build(N):
    t = np.arange(N) / N
    u = np.full(N, c[0])
    for k in range(1, len(c)):
        u += c[k] * np.cos(2 * np.pi * k * t)
    lg = outer.conj_fn(u)
    q = np.exp(2j * np.pi * t + lg)
    return q, u


def Mmaj(r):
    s = 0.0; n = 1
    while n < 20000:
        t = n * (r ** n / (1 - r ** n) + 5 * r ** (5 * n) / (1 - r ** (5 * n)))
        s += t; n += 1
        if t < 1e-20:
            break
    return 1 + 6 * s


qs, us = build(1 << 20)
r = float(np.abs(qs).max())
Lip = 2 * math.pi * Mmaj(r) * (1 + Sk)
print("max|q| = %.9f   M(r) = %.4f   sum k|c_k| = %.6f   Lip_apriori = %.2f" % (r, Mmaj(r), Sk, Lip), flush=True)
den = c[0] + L - tau
for e in (20, 21, 22, 23, 24):
    N = 1 << e
    q, u = build(N)
    rr = float(np.abs(q).max())
    g = haupt.logabs_x_vec(H, q)
    RE = outer.rearr(g)
    eps = Lip / N
    print("N=2^%d  RE_N=%.12f  eps=%.3e   bound <= %.9f   margin >= %+.9f   (max|q|=%.9f)"
          % (e, RE, eps, (RE + eps + L) / den, m * den - (RE + eps + L), rr), flush=True)
    del q, u, g
