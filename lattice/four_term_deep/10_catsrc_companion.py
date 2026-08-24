#!/usr/bin/env python3
"""10_catsrc_companion.py

The genuine companion of the four-term signed-binomial transforms of Zagier E.

Discovery driving this script: T(b) (the signed binomial transform of the
three-term companion) is NOT a solution of the four-term row -- the second
solution of the ODE picks up an extra  -A_new(t) log(1+nu t)  under the gauge
(1+nu t)^{-1} y(t/(1+nu t)), because log(t/(1+nu t)) = log t - log(1+nu t).

So we compute the row's own companion in the project's convention
(05_report.py seq(...,start=1)):  u_{-2}=u_{-1}=u_0=0, u_1=1, recurrence
applied from n=1 onwards; and identify its Apery limit.
"""
import sys
from fractions import Fraction as Fr
from math import gcd, comb
import mpmath as mp
mp.mp.dps = 220

A0, C0, D0 = 12, 4, 32
N = 700
NK = 150

def sharp_k(seq, upto):
    L = 1; kmax = 0
    for n in range(1, upto+1):
        L = L*n//gcd(L, n)
        den = seq[n].denominator; kk = 0
        while den > 1:
            g = gcd(den, L)
            if g == 1: kk = 99; break
            den //= g; kk += 1
        kmax = max(kmax, kk)
    return kmax

def four_term(P, Q, R, N, start):
    """start=0: u_0=1 analytic (recurrence from n=0, u_{-1}=u_{-2}=0)
       start=1: companion u_0=0,u_1=1 (recurrence from n=1)"""
    pad = [Fr(0), Fr(0)]
    if start == 0:
        u = pad + [Fr(1)]; n0 = 0
    else:
        u = pad + [Fr(0), Fr(1)]; n0 = 1
    off = 2
    for n in range(n0, N):
        nx = (Fr(P[0]*n*n+P[1]*n+P[2])*u[n+off]
              - Fr(Q[0]*n*n+Q[1]*n+Q[2])*u[n-1+off]
              + Fr(R[0]*n*n+R[1]*n+R[2])*u[n-2+off]) / Fr((n+1)**2)
        u.append(nx)
    return u[off:]

GCAT = mp.catalan
lines = []
for nu in (1, 2, -1, -2):
    P = (A0-3*nu, A0-3*nu, C0-nu)
    Q = (D0-2*A0*nu+3*nu*nu, 0, 0)
    CC = -nu*(D0-A0*nu+nu*nu)
    R = (CC, -CC, 0)
    v = four_term(P, Q, R, N, 0)
    w = four_term(P, Q, R, N, 1)
    assert all(x.denominator == 1 for x in v)
    kk = sharp_k(w, NK)
    l1, l2 = sorted([abs(8-nu), abs(4-nu), abs(-nu)], reverse=True)[:2]
    print("="*70)
    print("nu = %+d :  P=%s Q=%s R=%s" % (nu, P, Q, R))
    print("  v_n[0..7] =", [int(x) for x in v[:8]])
    print("  companion w_n[0..7] =", [str(x) for x in w[:8]])
    print("  sharp k of companion (n<=%d) = %d ; score = %.4f"
          % (NK, kk, float(-mp.log(l2)-kk)))
    with mp.workdps(240):
        def rat(x): return mp.mpf(x.numerator)/mp.mpf(x.denominator)
        r1 = rat(w[N-1])/rat(v[N-1]); r2 = rat(w[N-2])/rat(v[N-2])
        conv = abs(r1-r2)
        xi = r1
    print("  xi = %s" % mp.nstr(xi, 70))
    print("  stability |xi_n - xi_{n-1}| = 10^%.1f" % float(mp.log10(conv)))
    # candidate identification
    lg = mp.log(mp.mpf(8)/mp.mpf(8-nu))
    rel = mp.pslq([mp.mpf(1), xi, GCAT, lg], maxcoeff=10**12, maxsteps=10**5)
    print("  pslq[1, xi, G, log(8/(8-nu))] = %s" % (rel,))
    if rel and rel[1] != 0:
        rec = -(rel[0] + rel[2]*GCAT + rel[3]*lg)/rel[1]
        print("     => xi = %s ; agree %s digits"
              % (mp.nstr(rec, 50), int(-mp.log10(abs(rec-xi)/abs(xi)))))
    lines.append("nu%+d %s" % (nu, mp.nstr(xi, 80)))

with open("/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/catsrc_limits.txt", "w") as fh:
    fh.write("\n".join(lines) + "\n")
print("\nwrote out/catsrc_limits.txt")
