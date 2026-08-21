import csv
from math import gcd, log, log2
from fractions import Fraction
from mpmath import mp, mpf, catalan
exec(open('/tmp/bridge.py').read().split("Z={int(r")[0])   # reuse helpers

Z={int(r['n']):(int(r['X_n']),int(r['Y_n'])) for r in csv.DictReader(open('zudilin_rows.csv'))}
N={int(r['n']):(int(r['V_n']),int(r['U_n'])) for r in csv.DictReader(open('nesterenko_rows.csv'))}
ns=sorted(set(Z)&set(N))
D={}; cur=1
for i in range(1,6*max(ns)+2): cur=cur*i//gcd(cur,i); D[i]=cur
mp.dps=4000; G=+catalan

def delta_of(q,p):
    if q==0: return None
    lf=mpf(q)*G-mpf(p)
    if lf==0: return None
    lq=log(abs(q)) if abs(q)<10**300 else float(mp.log(abs(mpf(q))))
    lq=float(mp.log(abs(mpf(q)))); le=float(mp.log(abs(lf)))
    return lq, le, 1-le/lq

print("  n | v2(h)/n | log2 M/n | log2 covol/n | H/n | F/n | delta   || baseline delta")
for n in ns:
    if n not in (20,40,60,80,90,98): continue
    X,Y=Z[n]; V,U=N[n]; S=D[6*n]**2
    a1,a2=X//S,V//S
    h=a1*U-a2*Y; vh=v2(h); T=1<<vh; M=S*T
    b1,b2=cong_basis(a1,a2,T); b1,b2=intersect(b1,b2,Y,U,M); b1,b2=gauss(b1,b2)
    covol=abs(b1[0]*b2[1]-b1[1]*b2[0])
    best=None
    for (c1,c2) in [b1,b2,(b1[0]+b2[0],b1[1]+b2[1]),(b1[0]-b2[0],b1[1]-b2[1])]:
        nq,np_=c1*X+c2*V, c1*Y+c2*U
        if nq%M or np_%M: continue
        r=delta_of(nq//M, np_//M)
        if r and (best is None or r[2]>best[2]): best=r
    bb1,bb2=gauss(*cong_basis(Y,U,S)); bb=None
    for (c1,c2) in [bb1,bb2]:
        nq,np_=c1*X+c2*V, c1*Y+c2*U
        if nq%S: continue
        r=delta_of(nq//S,np_//S)
        if r and (bb is None or r[2]>bb[2]): bb=r
    print(f" {n:3d} | {vh/n:7.3f} | {log2(M)/n:8.3f} | {log2(covol)/n:12.3f} |"
          f" {best[0]/n if best else 0:6.3f} | {best[1]/n if best else 0:7.3f} | {best[2] if best else 0:.5f} || {bb[2] if bb else 0:.5f}")
