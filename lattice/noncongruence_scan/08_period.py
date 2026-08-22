#!/usr/bin/env python3
"""Apery limit xi of a row in the class
   (n+1)^2 a_{n+1} = (al n^2+be n+ga) a_n - (de n^2+ep n+ze) a_{n-1}
via the exact Casoratian series, then identification attempts."""
from fractions import Fraction as Fr
import sys, mpmath as mp

def rows(al,be,ga,de,ep,ze,N):
    a=[Fr(1),Fr(ga)]; b=[Fr(0),Fr(1)]
    for n in range(1,N):
        P=al*n*n+be*n+ga; Q=de*n*n+ep*n+ze
        a.append((P*a[n]-Q*a[n-1])/Fr((n+1)**2))
        b.append((P*b[n]-Q*b[n-1])/Fr((n+1)**2))
    return a,b

def xi(al,be,ga,de,ep,ze,N=420,dps=200):
    mp.mp.dps=dps
    a,b=rows(al,be,ga,de,ep,ze,N)
    W=[a[m]*b[m+1]-a[m+1]*b[m] for m in range(N)]
    S=Fr(0)
    for m in range(N-1): S+=W[m]/(a[m]*a[m+1])
    return mp.mpf(S.numerator)/mp.mpf(S.denominator)

if __name__=="__main__":
    cases=[]
    for line in sys.stdin:
        p=line.split()
        if len(p)<6: continue
        cases.append(tuple(int(x) for x in p[:6]))
    mp.mp.dps=120
    for c in cases:
        al,be,ga,de,ep,ze=c
        v=xi(al,be,ga,de,ep,ze)
        print("(al,be,ga,de,ep,ze)=%s\n   xi = %s"%(str(c), mp.nstr(v,60)))
        # structured identification
        basis=[("1",mp.mpf(1)),("xi",v)]
        cands={}
        for m in range(2,60):
            cands["log(%d/%d)"%(m+1,m-1)]=mp.log(mp.mpf(m+1)/(m-1))
        cands["pi^2"]=mp.pi**2; cands["zeta(3)"]=mp.zeta(3); cands["log2^2"]=mp.log(2)**2
        hit=None
        for name,val in cands.items():
            r=mp.pslq([v,val,mp.mpf(1)],maxcoeff=10**12,maxsteps=8000)
            if r and r[0]!=0 and max(abs(x) for x in r)<10**10:
                hit=(name,r); break
        print("   ->", hit if hit else "no 3-term hit in the log/pi^2 battery")
