import csv
from math import gcd, log
from mpmath import mp, mpf, catalan
exec(open('/tmp/bridge.py').read().split("Z={int(r")[0])
def gauss2(v1,v2):
    v1,v2=list(v1),list(v2)
    dot=lambda a,b:a[0]*b[0]+a[1]*b[1]
    while True:
        if dot(v2,v2)<dot(v1,v1): v1,v2=v2,v1
        n1=dot(v1,v1)
        if n1==0: return tuple(v1),tuple(v2)
        mu=(2*dot(v1,v2)+n1)//(2*n1) if dot(v1,v2)>=0 else -((-2*dot(v1,v2)+n1)//(2*n1))
        if mu==0: return tuple(v1),tuple(v2)
        v2=[v2[0]-mu*v1[0],v2[1]-mu*v1[1]]
def v2f(x):
    x=abs(x); return None if x==0 else (x&-x).bit_length()-1
Z={int(r['n']):(int(r['X_n']),int(r['Y_n'])) for r in csv.DictReader(open('zudilin_rows.csv'))}
N={int(r['n']):(int(r['V_n']),int(r['U_n'])) for r in csv.DictReader(open('nesterenko_rows.csv'))}
D={}; cur=1
for i in range(1,600): cur=cur*i//gcd(cur,i); D[i]=cur
mp.dps=5000; G=+catalan
print(" n  v2(h)/n  logS/n  sigma  covolK/n  M?  H_n     F_n    log|l|  log|l|/n")
for n in list(range(2,99,4))+[98]:
    X,Y=Z[n]; V,U=N[n]; S=D[6*n]**2
    a1,a2=X//S,V//S
    h=a1*U-a2*Y; vh=v2f(h); T=1<<vh; M=S*T
    b1,b2=cong_basis(a1,a2,T); b1,b2=intersect(b1,b2,Y,U,M)
    covolK=abs(b1[0]*b2[1]-b1[1]*b2[0])
    E1=float(mp.log(abs(mpf(X)*G-mpf(Y))))/n; E2=float(mp.log(abs(mpf(V)*G-mpf(U))))/n
    sig=log(M)/n; x=(sig+E2-E1)/2
    w=max(1,int(mp.floor(mp.e**((2*x-sig)*n))))
    r1,r2=gauss2((b1[0],b1[1]*w),(b2[0],b2[1]*w))
    best=None
    for r in (r1,r2):
        c1,c2=r[0], r[1]//w
        nq,np_=c1*X+c2*V, c1*Y+c2*U
        if nq%M or np_%M: continue
        q,p=nq//M,np_//M
        if q==0: continue
        e=mpf(q)*G-mpf(p)
        if e==0: continue
        lq=float(mp.log(abs(mpf(q)))); le=float(mp.log(abs(e)))
        if best is None or le<best[1]: best=(lq,le)
    if best is None: print(n,"none"); continue
    print(f"{n:3d} {vh/n:8.3f} {log(S)/n:7.3f} {sig:7.3f} {log(covolK)/n:8.3f} {covolK==M} {best[0]/n:7.3f} {(E1+E2-sig)/2:7.3f} {best[1]:10.2f} {best[1]/n:8.4f}")
