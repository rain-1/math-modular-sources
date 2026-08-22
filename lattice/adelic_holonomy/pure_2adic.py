"""2-adic slope of the pure polylogarithm module on the level-8 Catalan host.
Host: outer singularity s = 1/lambda_2 = 1/4;  w(x)=s x/(x-s);  y = x + w(x).
Pure module: Li_j(x/s) = sum_n (4x)^n/n^j  (integral coefficients 4^n/n^j).
We expand Sym^{+/-} Li_j in the invariant coordinate y and read off v_2."""
from fractions import Fraction as F
N = 40
def mul(a,b):
    c=[F(0)]*N
    for i,ai in enumerate(a):
        if ai:
            for k in range(0,N-i):
                if b[k]: c[i+k]+=ai*b[k]
    return c
def comp(f,g):                       # f(g(x)), g[0]==0
    r=[F(0)]*N; p=[F(0)]*N; p[0]=F(1)
    for n in range(N):
        if f[n]:
            for k in range(N): r[k]+=f[n]*p[k]
        p=mul(p,g)
    return r
s=F(1,4)
x=[F(0)]*N; x[1]=F(1)
# w = s x/(x-s) = -x/(1-x/s) * ... compute as series:  w = s*x/(x-s)
# 1/(x-s) = -1/s * 1/(1-x/s) = -1/s * sum (x/s)^n
inv=[F(0)]*N
for n in range(N): inv[n]=-(F(1)/s)*(F(1)/s)**n
w=mul([F(0),F(1)]+[F(0)]*(N-2),[s*c for c in inv])
y=[w[i]+x[i] for i in range(N)]
def li(j):
    f=[F(0)]*N
    for n in range(1,N): f[n]=F(4)**n/F(n)**j
    return f
def tosy(g):
    """write a symmetric series g(x) as sum c_k y^k (y has x-valuation 2)."""
    c=[]; h=g[:]; yp=[F(0)]*N; yp[0]=F(1)
    for k in range(0, N//2):
        # leading x-order of y^k is 2k
        v=2*k
        ck=h[v]/(yp[v] if yp[v]!=0 else F(1))
        c.append(ck)
        h=[h[i]-ck*yp[i] for i in range(N)]
        yp=mul(yp,y)
    return c
def v2(q):
    if q==0: return None
    n,d=q.numerator,q.denominator; v=0
    while n%2==0: n//=2; v+=1
    while d%2==0: d//=2; v-=1
    return v
print("host s=1/4 (lambda_2=4);  y = x + w(x) = %s ..."%([str(y[i]) for i in range(2,6)]))
for lab,sg in (("Sym+",1),("Sym-",-1)):
    for j in range(1,5):
        f=li(j); fw=comp(f,w)
        g=[f[i]+sg*fw[i] for i in range(N)]
        c=tosy(g)
        vs=[(k,v2(c[k])) for k in range(1,13)]
        sl=[ (v/k if v is not None else None) for k,v in vs]
        print(f"  {lab} Li_{j}: v2(c_k), k=1..12 = {[v for _,v in vs]}")
        nz=[(k,v) for k,v in vs if v is not None]
        if nz: print(f"        min v2(c_k)/k over k<=12 = {min(v/k for k,v in nz):.3f}")
