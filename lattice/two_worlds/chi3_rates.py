import sys, math
from fractions import Fraction as F
sys.path.insert(0,'/home/ubuntu/code/math-modular-sources/lattice/two_worlds')
from chi3_forms import *
from chi3_family import numerator, row
import mpmath as mp
mp.mp.dps = 200
L3 = mp.nsum(lambda n: (mp.mpf(1)/(3*n+1)**2 - mp.mpf(1)/(3*n+2)**2), [0, mp.inf])

def v3(x):
    if x == 0: return None
    n, d = x.numerator, x.denominator
    v = 0
    while n % 3 == 0: n//=3; v+=1
    while d % 3 == 0: d//=3; v-=1
    return v

def cofactor3(x):
    """the part of den(x) prime to 3"""
    d = x.denominator
    while d % 3 == 0: d//=3
    return d

def report(cname, af, ms):
    print("=== %s ===" % cname)
    print("  m    log|Q|/m   log|form|/m   v3(Q)/m   den'(Q)  v3(P)  den'(P) log(den'(P))/m")
    for m in ms:
        a = af(m)
        Q, P = row(m, a)
        form = mp.mpf(Q.numerator)/Q.denominator*L3 - mp.mpf(P.numerator)/P.denominator
        lq = mp.log(abs(mp.mpf(Q.numerator)/Q.denominator))/m
        lf = mp.log(abs(form))/m
        dq = cofactor3(Q); dp = cofactor3(P)
        print("  %3d  %9s  %11s  %7s  %s  %5s  %s  %s" % (
            m, mp.nstr(lq,7), mp.nstr(lf,7), mp.nstr(mp.mpf(v3(Q))/m,6),
            ("1" if dq==1 else str(dq)[:20]+"..."), v3(P),
            ("1" if dp==1 else str(dp)[:14]+"..."), mp.nstr(mp.log(dp)/m,6) if dp>1 else "0"))

ms = list(range(4, 31, 2))
report("a=2m", lambda m: 2*m, ms)
report("a=7m//4", lambda m: (7*m)//4, ms)
report("a=3m//2", lambda m: (3*m)//2, ms)
