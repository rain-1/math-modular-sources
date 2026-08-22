#!/usr/bin/env python3
"""Deep analysis of every non-degenerate hit with |lambda_2| < 1 (the only ones
that can have a positive score)."""
import glob, os, math, sys
from fractions import Fraction as F
from math import gcd
NK=60
D=[1]*(NK+2)
for n in range(1,NK+2): D[n]=D[n-1]*n//gcd(D[n-1],n)
def seq(A,be,B,M,j1,j2,C,N,u0,u1):
    de=C*M*M; ep=-C*M*(j1+j2); ze=C*j1*j2
    u=[F(u0),F(u1)]
    for n in range(1,N):
        P=A*n*n+be*n+B; Qn=de*n*n+ep*n+ze
        u.append((P*u[n]-Qn*u[n-1])/F((n+1)**2))
    return u
cands=[]
for p in sorted(glob.glob('out/c_*.txt'))+['out/a0.txt']:
    if not os.path.exists(p): continue
    for line in open(p):
        f=line.split()
        if len(f)!=6: continue
        M,j1,j2,A,B,C=map(int,f)
        if (j1%M==0 and j1//M>=1) or (j2%M==0 and j2//M>=1): continue
        Dd=C*M*M; disc=A*A-4*Dd
        if disc<=0: continue
        r1=(A+math.sqrt(disc))/2.0; r2=(A-math.sqrt(disc))/2.0
        l2=min(abs(r1),abs(r2)); l1=max(abs(r1),abs(r2))
        if l2>=1.0: continue
        cands.append((l2,l1,M,j1,j2,A,B,C))
cands.sort()
print("non-degenerate real-root hits with |lam2|<1 :",len(cands))
res=[]
for (l2,l1,M,j1,j2,A,B,C) in cands:
    S=2*M-j1-j2; be=A*S//(2*M)
    b=seq(A,be,B,M,j1,j2,C,NK,0,1)
    k=0
    while k<=6:
        if all((D[n]**k*b[n]).denominator==1 for n in range(NK+1)): break
        k+=1
    sc=math.log(1/l2)-k
    res.append((sc,k,l2,l1,M,j1,j2,A,B,C))
res.sort(reverse=True)
print("top 25 by score:")
for r in res[:25]:
    print("  score=%+.4f k=%d |lam2|=%.6f lam1=%.4f class(%d;%d,%d) A=%d B=%d C=%d"%r)
print("number with score>0 and k>=2:", len([r for r in res if r[0]>0 and r[1]>=2]))
for r in res:
    if r[0]>0 and r[1]>=2: print("   POSITIVE k>=2:", r)
