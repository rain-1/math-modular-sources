import sys, math
from fractions import Fraction as F
sys.path.insert(0, '/home/ubuntu/code/math-modular-sources/lattice/two_worlds')
from chi3_forms import *
import mpmath as mp
mp.mp.dps = 60

L3 = mp.nsum(lambda n: (mp.mpf(1)/(3*n+1)**2 - mp.mpf(1)/(3*n+2)**2), [0, mp.inf])
Z2 = mp.zeta(2)

def numer_beukers(m):
    """N = (2t+m+1) * prod_{j=1..m}(t-j+1) * prod_{j=1..m}(t+m+j)  -- sigma-invariant x odd"""
    N = linear(2, m+1)
    for j in range(1, m+1):
        N = pmul(N, linear(1, -j+1))
    for j in range(1, m+1):
        N = pmul(N, linear(1, m+j))
    return N

def numer_plain(m):
    return linear(2, m+1)

def numer_half(m):
    """N = (2t+m+1)*prod_{j=1..m}(t-j+1)(t+m+j) but only half the reflected factors"""
    N = linear(2, m+1)
    for j in range(1, m+1):
        N = pmul(N, pmul(linear(1, -j+1), linear(1, m+j)))
    return N

def check(name, numf, mmax=8, verbose=True):
    print("=== %s ===" % name)
    rows = []
    for m in range(0, mmax+1):
        P = poles_thirds(m)
        N = numf(m)
        if len(N)-1 > 2*len(P)-2:
            print(" m=%d numerator degree %d too big (need <=%d)" % (m, len(N)-1, 2*len(P)-2)); continue
        d = linform(P, N)
        q = d['L']; p = -d['Crat']
        val = mp.mpf(q.numerator)/q.denominator*L3 - mp.mpf(p.numerator)/p.denominator \
              + mp.mpf(d['zeta2'].numerator)/d['zeta2'].denominator*Z2
        rows.append((m, q, p, d, val))
        if verbose:
            print(" m=%2d  zeta2coef=%s  Lcoef=%s  pi/sqrt3=%s  |form|=%s" %
                  (m, d['zeta2'], d['L'], d['pi_over_sqrt3'], mp.nstr(val, 10)))
    return rows

check("plain N=(2t+m+1)", numer_plain, 6)
check("beukers-like", numer_beukers, 6)
check("half", numer_half, 5)
