import sys, math
from fractions import Fraction as F
sys.path.insert(0,'/home/ubuntu/code/math-modular-sources/lattice/two_worlds')
from chi3_row import row
from chi3_asym import maxG, maxH
import mpmath as mp
mp.mp.dps = 40
L3 = mp.nsum(lambda n: (mp.mpf(1)/(3*n+1)**2 - mp.mpf(1)/(3*n+2)**2), [0, mp.inf])

def split3(x):
    n,d = x.numerator, x.denominator; v=0
    while n%3==0: n//=3; v+=1
    while d%3==0: d//=3; v-=1
    return v, d          # v3 and the prime-to-3 part of the denominator

def fl(x): return mp.mpf(x.numerator)/x.denominator

ALPHAS = [(1,2),(3,4),(1,1),(5,4),(3,2),(7,4),(15,8),(2,1)]
print("%-6s %-4s %10s %10s %10s %8s %8s %8s %8s" %
      ("alpha","m","log|Q|/m","log|f|/m","logLam*","-v3Q/m","nu(Q)/m","kd(P)/m","-v3P/m"))
for (p,q) in ALPHAS:
    al = mp.mpf(p)/q
    xg,lg = maxG(al,1); yh,lh = maxH(al,1)
    poly = 2*al-4        # log-m exponent
    ms = [q*k for k in range(4, 25, 4)]
    for m in ms:
        a = p*m//q
        Q,P = row(m,a)
        v3Q, dQ = split3(Q); v3P, dP = split3(P)
        form = fl(Q)*L3 - fl(P)
        print("%-6s %-4d %10s %10s %10s %8s %8s %8s %8s" % (
            "%d/%d"%(p,q), m,
            mp.nstr(mp.log(abs(fl(Q)))/m,7), mp.nstr(mp.log(abs(form))/m,7),
            mp.nstr(lh,7),
            mp.nstr(mp.mpf(-v3Q)/m,5), mp.nstr(mp.log(dQ)/m if dQ>1 else 0,5),
            mp.nstr(mp.log(dP)/m if dP>1 else 0,5), mp.nstr(mp.mpf(-v3P)/m,5)))
    print("   [asymptotic:  logLambda=%s  loglambda=%s  (times m^%s)]" % (
        mp.nstr(lh,8), mp.nstr(lg,8), mp.nstr(poly,4)))
    print()
