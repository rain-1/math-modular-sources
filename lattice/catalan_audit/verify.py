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
mp.dps=3000; G=+catalan
for n in (20,40):
    X,Y=Z[n]; V,U=N[n]; S=D[6*n]**2
    a1,a2=X//S,V//S
    print(f"n={n} v2(X)={v2f(X)} v2(Y)={v2f(Y)} v2(V)={v2f(V)} v2(U)={v2f(U)} v2(a1)={v2f(a1)} v2(a2)={v2f(a2)}")
    h=a1*U-a2*Y
    print("  v2(h)=",v2f(h)," min(v2(a1)+v2(U),v2(a2)+v2(Y))=",min(v2f(a1)+v2f(U),v2f(a2)+v2f(Y)))
    T=1<<v2f(h); M=S*T
    b1,b2=cong_basis(a1,a2,T); b1,b2=intersect(b1,b2,Y,U,M)
    # check basis really in K
    for b in (b1,b2):
        print("   basis",("OK" if (a1*b[0]+a2*b[1])%T==0 and (Y*b[0]+U*b[1])%M==0 else "BAD"))
    E1=float(mp.log(abs(mpf(X)*G-mpf(Y))))/n; E2=float(mp.log(abs(mpf(V)*G-mpf(U))))/n
    sig=log(M)/n; x=(sig+E2-E1)/2
    w=int(mp.floor(mp.e**((2*x-sig)*n)));  w=max(w,1)
    s1=(b1[0],b1[1]*w); s2=(b2[0],b2[1]*w)
    r1,r2=gauss2(s1,s2)
    for r in (r1,r2,(r1[0]+r2[0],r1[1]+r2[1]),(r1[0]-r2[0],r1[1]-r2[1])):
        c1,c2=r[0], r[1]//w
        exact = (r[1]%w==0)
        nq,np_=c1*X+c2*V, c1*Y+c2*U
        div = (nq%M==0 and np_%M==0)
        if not div: 
            print(f"   cand exact_lattice={exact} divisible=False  -> skipped"); continue
        q,p=nq//M,np_//M
        if q==0: print("   q=0"); continue
        e=mpf(q)*G-mpf(p)
        print(f"   cand exact_lattice={exact} divisible=True q!=0  log|q|/n={float(mp.log(abs(mpf(q))))/n:.3f} log|qG-p|/n={float(mp.log(abs(e)))/n:.3f} gcd(q,p)={gcd(q,p)} v2(q)={v2f(q)}")
        print(f"      recheck cong: aC%T={ (a1*c1+a2*c2)%T ==0 }  YC%M={ (Y*c1+U*c2)%M==0 }")
        print(f"      |c1|={abs(c1)} bits={abs(c1).bit_length()} |c2| bits={abs(c2).bit_length()}")
        # direct: is q*G-p really small? compare p/q to G
        print("      p/q vs G first digits:", mp.nstr(mpf(p)/mpf(q),30), mp.nstr(G,30))
