from math import gcd, log, log2
from fractions import Fraction
from mpmath import mp, mpf, zeta
def v2(x):
    x=abs(x); return None if x==0 else (x&-x).bit_length()-1
def ext_gcd(a,b):
    o_r,r=a,b; o_s,s=1,0; o_t,t=0,1
    while r: q=o_r//r; o_r,r=r,o_r-q*r; o_s,s=s,o_s-q*s; o_t,t=t,o_t-q*t
    return o_r,o_s,o_t
def cong_basis(a,b,m):
    if m==1: return (1,0),(0,1)
    g1,p,q=ext_gcd(a,b)
    if g1==0: return (1,0),(0,1)
    ap,bp=a//g1,b//g1; g=gcd(g1,m); m2=m//g
    return (m2*p,m2*q),(bp,-ap)
def intersect(b1,b2,a,b,m):
    A=a*b1[0]+b*b1[1]; B=a*b2[0]+b*b2[1]
    w1,w2=cong_basis(A,B,m)
    c=lambda w:(w[0]*b1[0]+w[1]*b2[0], w[0]*b1[1]+w[1]*b2[1])
    return c(w1),c(w2)
def gauss(v1,v2):
    v1,v2=list(v1),list(v2); dot=lambda a,b:a[0]*b[0]+a[1]*b[1]
    while True:
        if dot(v2,v2)<dot(v1,v1): v1,v2=v2,v1
        n1=dot(v1,v1)
        if n1==0: return tuple(v1),tuple(v2)
        d=dot(v1,v2); mu=(2*d+n1)//(2*n1) if d>=0 else -((-2*d+n1)//(2*n1))
        if mu==0: return tuple(v1),tuple(v2)
        v2=[v2[0]-mu*v1[0],v2[1]-mu*v1[1]]
NM=700
def rows(a,b,C):
    A=[0]*(NM+2); c=[0]*(NM+2); A[0],A[1]=1,b; c[0],c[1]=0,1
    for n in range(1,NM):
        A[n+1]=((2*n+1)*(a*n*n+a*n+b)*A[n]-C*n**3*A[n-1])//(n+1)**3
        c[n+1]=(2*n+1)*(a*n*n+a*n+b)*c[n]-C*n**6*c[n-1]
    return A,c
AD,CD=rows(10,4,64); AT,CT=rows(12,4,16)
fact=[1]*(NM+2)
for i in range(1,NM+2): fact[i]=fact[i-1]*i
D={}; cur=1
for i in range(1,3*100+2): cur=cur*i//gcd(cur,i); D[i]=cur
mp.dps=2000; TH=+zeta(3)
print("CONTROL: zeta(3) Domb+T through the SAME pipeline.  Theory: F=+1.163, delta=0.9010")
print("  n | log|q|/n | log|qZ3-p|/n | delta   | predicted F/n | v2(h)/n")
for n in (10,20,30,40,50):
    m,k=2*n,3*n
    S=D[3*n]**3
    Y1f=Fraction(S*24*CD[m], fact[m]**3); Y2f=Fraction(S*32*CT[k], fact[k]**3)
    if Y1f.denominator!=1 or Y2f.denominator!=1: print(n,"non-integral Y"); continue
    Y1,Y2=int(Y1f),int(Y2f)
    a1,a2=7*AD[m],7*AT[k]
    X1,X2=S*a1,S*a2
    h=a1*Y2-a2*Y1; T=1<<v2(h); M=S*T
    b1,b2=cong_basis(a1,a2,T); b1,b2=intersect(b1,b2,Y1,Y2,M)
    sig=log(M)/n
    E1=float(mp.log(abs(mpf(X1)*TH-mpf(Y1))))/n
    E2=float(mp.log(abs(mpf(X2)*TH-mpf(Y2))))/n
    x=(sig+E2-E1)/2
    w=int(mp.floor(mp.e**((2*x-sig)*n))) or 1
    r1,r2=gauss((b1[0],b1[1]*w),(b2[0],b2[1]*w))
    best=None
    for r in (r1,r2,(r1[0]+r2[0],r1[1]+r2[1]),(r1[0]-r2[0],r1[1]-r2[1])):
        c1,c2=r[0],r[1]//w
        if (c1,c2)==(0,0): continue
        nq,np_=c1*X1+c2*X2, c1*Y1+c2*Y2
        if nq%M or np_%M: continue
        q,p=nq//M,np_//M
        if q==0: continue
        e=mpf(q)*TH-mpf(p)
        if e==0: continue
        lq=float(mp.log(abs(mpf(q)))); le=float(mp.log(abs(e)))
        if best is None or le<best[1]: best=(lq,le)
    if best: print("  %3d | %8.3f | %12.3f | %7.4f | %+13.4f | %7.3f"%(n,best[0]/n,best[1]/n,1-best[1]/best[0],(E1+E2-sig)/2,v2(h)/n))
