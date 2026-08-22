import sys
from fractions import Fraction as F
sys.path.insert(0,'/home/ubuntu/code/math-modular-sources/lattice/two_worlds')
from chi3_row import row
from chi3_asym import maxG, maxH
import mpmath as mp

def L3(dps):
    mp.mp.dps = dps
    return mp.nsum(lambda n: (mp.mpf(1)/(3*n+1)**2 - mp.mpf(1)/(3*n+2)**2), [0, mp.inf])

def split3(x):
    n,d = x.numerator, x.denominator; v=0
    while n%3==0: n//=3; v+=1
    while d%3==0: d//=3; v-=1
    return v, d

ALPHAS = [(1,2),(3,4),(1,1),(5,4),(3,2),(7,4),(2,1)]
print("%-5s %-5s %11s %11s %11s %11s %8s %8s %8s" %
      ("alpha","m","log|Q|/m","logLam_pred","log|f|/m","loglam_pred","kap3","nu","kd"))
for (p,q) in ALPHAS:
    al = mp.mpf(p)/q
    xg,lg = maxG(al,1); yh,lh = maxH(al,1)
    poly = float(2*al-4)
    for m in [q*k for k in (6,12,18,24,30)]:
        a = p*m//q
        dps = 40 + int(4*m)
        Lv = L3(dps)
        Q,P = row(m,a)
        v3Q,dQ = split3(Q); v3P,dP = split3(P)
        fq = mp.mpf(Q.numerator)/Q.denominator
        fp = mp.mpf(P.numerator)/P.denominator
        form = fq*Lv - fp
        # remove the m^poly factor from the measured rates for comparison
        corr = poly*mp.log(m)/m
        print("%-5s %-5d %11s %11s %11s %11s %8s %8s %8s" % (
            "%d/%d"%(p,q), m,
            mp.nstr(mp.log(abs(fq))/m - corr, 7), mp.nstr(lh,7),
            mp.nstr(mp.log(abs(form))/m - corr, 7), mp.nstr(lg,7),
            mp.nstr(mp.mpf(-v3Q)/m,5),
            mp.nstr(mp.log(dQ)/m if dQ>1 else 0,5),
            mp.nstr(mp.log(dP)/m if dP>1 else 0,5)))
    print()
