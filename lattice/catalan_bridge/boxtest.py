import csv
from math import gcd, log, log2, exp
from mpmath import mp, mpf, catalan
def v2(x):
    x=abs(x); return None if x==0 else (x&-x).bit_length()-1
def gauss(v1,v2):
    v1,v2=list(v1),list(v2)
    dot=lambda a,b:a[0]*b[0]+a[1]*b[1]
    while True:
        if dot(v2,v2)<dot(v1,v1): v1,v2=v2,v1
        n1=dot(v1,v1)
        if n1==0: return tuple(v1),tuple(v2)
        mu=(2*dot(v1,v2)+n1)//(2*n1) if dot(v1,v2)>=0 else -((-2*dot(v1,v2)+n1)//(2*n1))
        if mu==0: return tuple(v1),tuple(v2)
        v2=[v2[0]-mu*v1[0],v2[1]-mu*v1[1]]
exec(open('/tmp/bridge.py').read().split("def gauss")[0].split("Z={int(r")[0])
Z={int(r['n']):(int(r['X_n']),int(r['Y_n'])) for r in csv.DictReader(open('zudilin_rows.csv'))}
N={int(r['n']):(int(r['V_n']),int(r['U_n'])) for r in csv.DictReader(open('nesterenko_rows.csv'))}
D={}; cur=1
for i in range(1,600): cur=cur*i//gcd(cur,i); D[i]=cur
mp.dps=2500; G=+catalan
print("BOX-CONSTRAINED (anisotropic) selection -- the decisive test")
print("  n | log|q|/n | log|qG-p|/n | delta   | predicted F/n")
for n in (20,40,60,80,98):
    X,Y=Z[n]; V,U=N[n]; S=D[6*n]**2
    a1,a2=X//S,V//S
    h=a1*U-a2*Y; T=1<<v2(h); M=S*T
    b1,b2=cong_basis(a1,a2,T); b1,b2=intersect(b1,b2,Y,U,M)
    sig=log(M)/n
    E1=float(mp.log(abs(mpf(X)*G-mpf(Y))))/n
    E2=float(mp.log(abs(mpf(V)*G-mpf(U))))/n
    x=(sig+E2-E1)/2
    w=int(mp.floor(mp.e**((2*x-sig)*n)))     # weight on c2
    if w<1: w=1
    s1=(b1[0], b1[1]*w); s2=(b2[0], b2[1]*w)
    r1,r2=gauss(s1,s2)
    best=None
    for r in (r1,r2,(r1[0]+r2[0],r1[1]+r2[1]),(r1[0]-r2[0],r1[1]-r2[1])):
        c1,c2=r[0], r[1]//w
        if (c1,c2)==(0,0): continue
        nq,np_=c1*X+c2*V, c1*Y+c2*U
        if nq%M or np_%M: continue
        q,p=nq//M,np_//M
        if q==0: continue
        e=mpf(q)*G-mpf(p)
        if e==0: continue
        lq=float(mp.log(abs(mpf(q)))); le=float(mp.log(abs(e)))
        if best is None or le<best[1]: best=(lq,le)
    if best: print("  %3d | %8.3f | %11.3f | %7.4f | %+.4f"%(n,best[0]/n,best[1]/n,1-best[1]/best[0],(E1+E2-sig)/2))
    else: print("  %3d | no valid vector"%n)
