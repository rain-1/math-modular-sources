#!/usr/bin/env python3
"""10_catsrc_identify.py -- closed form for the four-term companion's Apery limit.

Derivation (see log): with L_0 the Zagier-E operator,  L_0 B = t  for the
three-term companion and L' W = t for the four-term one, while the gauge
transform sends  L_0 y = t  to  L' (Ty) = t/(1+nu t).  Hence
   W = T(B) + T(y),   (n+1)^2 y_{n+1} = P_0 y_n - Q_0 y_{n-1} + nu^n (n>=1), y_0=y_1=0,
and by the Casoratian  (m+1)^2 (a_m b_{m+1}-a_{m+1}b_m) = 32^m,
   lim y_n/a_n = sum_{m>=1} (nu/32)^m (xi a_m - b_m).
Therefore
   XI_NEW = xi * A(nu/32) - B(nu/32),      xi = G/2,
with A(t)=sum a_n t^n, B(t)=sum b_n t^n the Zagier-E row and companion series.
"""
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 120
A0, C0, D0 = 12, 4, 32
N = 400
def tt(N, start):
    P = lambda n: A0*(n*n+n)+C0; Q = lambda n: D0*n*n
    u = [Fr(1), Fr(P(0))] if start == 0 else [Fr(0), Fr(1)]
    for n in range(1, N): u.append((P(n)*u[n]-Q(n)*u[n-1])/Fr((n+1)**2))
    return u
a = tt(N, 0); b = tt(N, 1)
XI = mp.catalan/2
def ser(seq, x):
    s = mp.mpf(0); p = mp.mpf(1)
    for n in range(N):
        s += (mp.mpf(seq[n].numerator)/mp.mpf(seq[n].denominator))*p
        p *= x
    return s
meas = {1: "0.4863068456612505674910777191229431581349968112070476667871014530573948",
        2: "0.5203715170059685771322550428796926054562462673777764908045843387208894",
        -1: "0.433916123619711430683647037867754208338259942801335546203939452940266",
        -2: "0.4131171653710453594220235058204506785721254797785678398499641881505225"}
for nu in (1, 2, -1, -2):
    x = mp.mpf(nu)/32
    pred = XI*ser(a, x) - ser(b, x)
    m = mp.mpf(meas[nu])
    print("nu=%+d  A(nu/32)=%s" % (nu, mp.nstr(ser(a, x), 40)))
    print("        B(nu/32)=%s" % mp.nstr(ser(b, x), 40))
    print("        predicted xi_new = %s" % mp.nstr(pred, 60))
    print("        measured  xi_new = %s" % mp.nstr(m, 60))
    print("        agreement: %d digits" % int(-mp.log10(abs(pred-m)/abs(m))))
