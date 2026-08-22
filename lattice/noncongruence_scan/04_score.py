#!/usr/bin/env python3
"""For every integral row found by 01_class_scan, compute
   lam1, lam2, the companion b_n (census convention b_0=0,b_1=1),
   k = min{ k : d_n^k b_n in Z for n<=NK }, and score = log(1/|lam2|) - k."""
import sys, math
from fractions import Fraction
from math import gcd

NK = 40
D = [1]*(NK+2)
for n in range(1, NK+2): D[n] = D[n-1]*n//gcd(D[n-1], n)

def row(al,be,ga,de,ep,ze,N,a0,a1):
    a=[Fraction(a0),Fraction(a1)]
    for n in range(1,N):
        P=al*n*n+be*n+ga; Q=de*n*n+ep*n+ze
        a.append((P*a[n]-Q*a[n-1])/Fraction((n+1)**2))
    return a

def main(path, out):
    seen=set(); res=[]
    for line in open(path):
        t=tuple(map(int,line.split()))
        if t[2]==0: continue
        if t in seen: continue
        seen.add(t)
        al,be,ga,de,ep,ze=t
        Dd=al*al-4*de
        if Dd<0:                       # complex conjugate pair
            l1=l2=math.sqrt(abs(de)); cplx=1
        else:
            r1=(al+math.sqrt(Dd))/2; r2=(al-math.sqrt(Dd))/2
            l1,l2=(r1,r2) if abs(r1)>=abs(r2) else (r2,r1); cplx=0
        if abs(l2)<1e-13: continue
        # Apery ingredient 3: Casoratian W_n = Q(n) W_{n-1}/(n+1)^2 must never vanish
        deg_=False
        for n in range(1,200):
            if de*n*n+ep*n+ze==0: deg_=True; break
        if deg_: continue
        a=row(al,be,ga,de,ep,ze,NK,1,ga)
        if any(x.denominator!=1 for x in a): continue     # re-verify integrality
        b=row(al,be,ga,de,ep,ze,NK,0,1)
        k=0
        while k<=8:
            if all((D[n]**k*b[n]).denominator==1 for n in range(NK+1)): break
            k+=1
        sc=math.log(1/abs(l2))-k
        res.append((sc,al,ga,de,ze,l1,l2,k,cplx))
    res.sort(reverse=True)
    with open(out,'w') as f:
        for r in res:
            f.write("%.6f %d %d %d %d %.10g %.12g %d %d\n"%r)
    print("rows",len(res))
    for r in res[:40]:
        print("score %+.5f k=%d (al,ga,de,ze)=(%d,%d,%d,%d) lam1=%.6g lam2=%.8g%s"%(
            r[0],r[7],r[1],r[2],r[3],r[4],r[5],r[6]," CPLX" if r[8] else ""))

main(sys.argv[1], sys.argv[2])
