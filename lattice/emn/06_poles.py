#!/usr/bin/env python3
"""(f),(g): higher-pole EMN periods  I_{m,t} = \\iint_Delta x^m y^m /(1-x^2-y^2)^{t+1}.
Exact inner integral + high-precision theta quadrature, then PSLQ against (1,G).
Renamed to avoid clashes: the Catalan coefficient is bcoef, rational part acoef."""
from mpmath import mp, mpf, catalan, sin, cos, log, pi, quad, nstr, binomial, pslq, mpmathify
from fractions import Fraction
from math import comb
import sys

mp.dps = int(sys.argv[1]) if len(sys.argv)>1 else 120
GG = catalan

def inner_exact(W, m, t):
    """\\int_0^W u^m (1-u)^{-t-1} du, W in (0,1), 0<=t<=m."""
    v = 1-W                      # in (0,1)
    tot = mp.mpf(0)
    for i in range(m+1):
        c = mp.binomial(m,i)*(-1)**i
        if i == t:
            tot += c*(-mp.log(v))
        else:
            tot += c*(1 - v**(i-t))/(i-t)
    return tot

def Imt(m, t):
    def f(th):
        s2 = mp.sin(2*th)
        W  = 1/(1+s2)
        return (mp.cos(th)*mp.sin(th))**m * inner_exact(W, m, t)/2
    return quad(f, [0, pi/4, pi/2])

def ident(val, name):
    r = pslq([val, GG, mp.mpf(1)], maxcoeff=10**(mp.dps//3), maxsteps=100000)
    return r

def rat(x, maxden):
    from mpmath import pslq
    r = pslq([x, mp.mpf(1)], maxcoeff=maxden, maxsteps=200000)
    if r is None: return None
    return Fraction(-r[1], r[0])

print(f"dps = {mp.dps}")
print("=== (f) I_{m,t} = a_{m,t} + b_{m,t} G,  claim b = 4^{-m} C(m,m/2) C(m,t) ===")
print(f"{'m':>3} {'t':>3} {'b claimed':>22} {'b measured':>22} {'a (rational)':>34} ok")
rows=[]
for m in [0,2,4,6,8,10]:
    for t in range(0, m+1, 2):
        if t > m: continue
        val = Imt(m,t)
        bcl = Fraction(comb(m,m//2)*comb(m,t), 4**m)
        acoef = val - mp.mpf(bcl.numerator)/mp.mpf(bcl.denominator)*GG
        ar = rat(acoef, 10**(mp.dps//3))
        ok = "YES" if ar is not None else "?"
        rows.append((m,t,bcl,ar))
        print(f"{m:>3} {t:>3} {str(bcl):>22} {'(assumed)':>22} {str(ar):>34} {ok}")
print()
print("Cross-check against EMN arXiv:2510.20648 Thm 1.1.1 examples:")
for (m,t,tgt) in [(2,0,'-5/48 + 1/8 G'), (4,0,'-569/26880 + 3/128 G'), (4,2,'-49/384 + 9/64 G')]:
    for (mm,tt,bcl,ar) in rows:
        if (mm,tt)==(m,t):
            print(f"  I_{{{m},{t}}} = {ar} + {bcl} G      paper: {tgt}")
