"""18: row E (Catalan) - is the min-type companion really divergent?
Extend the hidden-variable coefficients g_k to k=180 and fit
log|g_k| = gamma k log k + delta k + ... ."""
from lib2v import *
from fractions import Fraction as F
from math import comb, log
N = 180
aE, bE = rowE(N)
g = []
for n in range(N+1):
    s = F(0)
    for k in range(n):
        s += g[k]*comb(n, k)*comb(2*n-2*k, n-k)
    g.append((F(bE[n]) - s))
print("row E: b_n = sum_k g_k C(n,k) C(2n-2k,n-k)")
print("   den(g_k) | [1..k]^2 for k<=60 :",
      all(lcmrange(k)**2 % g[k].denominator == 0 for k in range(1, 61)))
prev = None
print("     k   log|g_k|/k   d(log|g_k|/k)/d(log k)   [gamma; 0 => convergent]")
for k in [20, 40, 60, 80, 100, 120, 140, 160, 180]:
    L = (log(abs(g[k].numerator)) - log(g[k].denominator))/k
    if prev:
        k0, L0 = prev
        gam = (L - L0)/(log(k) - log(k0))
        print("   %4d   %9.4f            %7.4f" % (k, L, gam))
    else:
        print("   %4d   %9.4f              --" % (k, L))
    prev = (k, L)
print("   (gamma ~ 1 means |g_k| ~ (k!)^1 C^k : Gevrey-1, radius 0)")
