import csv
from math import gcd, log, log2
from fractions import Fraction

def ext_gcd(a,b):
    old_r,r=a,b; old_s,s=1,0; old_t,t=0,1
    while r: q=old_r//r; old_r,r=r,old_r-q*r; old_s,s=s,old_s-q*s; old_t,t=t,old_t-q*t
    return old_r,old_s,old_t

def cong_basis(a,b,m):
    """basis of {(c1,c2) in Z^2 : a c1 + b c2 = 0 mod m}"""
    if m==1: return (1,0),(0,1)
    g1,p,q=ext_gcd(a,b)
    if g1==0: return (1,0),(0,1)
    ap,bp=a//g1,b//g1
    g=gcd(g1,m); m2=m//g
    return (m2*p,m2*q),(bp,-ap)

def intersect(b1,b2,a,b,m):
    """{x*b1+y*b2 : a*(..)+b*(..) = 0 mod m}"""
    A=a*b1[0]+b*b1[1]; B=a*b2[0]+b*b2[1]
    w1,w2=cong_basis(A,B,m)
    def comb(w): return (w[0]*b1[0]+w[1]*b2[0], w[0]*b1[1]+w[1]*b2[1])
    return comb(w1),comb(w2)

def gauss(v1,v2):
    v1,v2=list(v1),list(v2)
    dot=lambda x,y:x[0]*y[0]+x[1]*y[1]
    while True:
        if dot(v2,v2)<dot(v1,v1): v1,v2=v2,v1
        n1=dot(v1,v1)
        if n1==0: return tuple(v1),tuple(v2)
        mu=round(Fraction(dot(v1,v2),n1))
        if mu==0: return tuple(v1),tuple(v2)
        v2=[v2[0]-mu*v1[0],v2[1]-mu*v1[1]]

def v2(x):
    x=abs(x)
    if x==0: return None
    return (x & -x).bit_length()-1

Z={int(r['n']):(int(r['X_n']),int(r['Y_n'])) for r in csv.DictReader(open('zudilin_rows.csv'))}
N={int(r['n']):(int(r['V_n']),int(r['U_n'])) for r in csv.DictReader(open('nesterenko_rows.csv'))}
ns=sorted(set(Z)&set(N))

# lcm table
D={}; cur=1
for i in range(1,6*max(ns)+2):
    cur=cur*i//gcd(cur,i); D[i]=cur

print(" n | v2(h)/n | log2(S)/n | log2(M)/n | sigma nats | delta_noT | delta_withT | q integral")
for n in ns:
    if n not in (10,20,30,40,50,60,70,80,90,98): continue
    X,Y=Z[n]; V,U=N[n]
    S=D[6*n]**2
    if X%S or V%S:
        print(f" {n}: S does not divide rows -- abort"); break
    a1,a2=X//S,V//S
    h=a1*U-a2*Y
    vh=v2(h); T=1<<vh; M=S*T
    # K = {c: a.c=0 mod T, (Y,U).c=0 mod M}
    b1,b2=cong_basis(a1,a2,T)
    b1,b2=intersect(b1,b2,Y,U,M)
    b1,b2=gauss(b1,b2)
    best=None; okint=True
    for (c1,c2) in [b1,b2,(b1[0]+b2[0],b1[1]+b2[1]),(b1[0]-b2[0],b1[1]-b2[1])]:
        nq=c1*X+c2*V; np_=c1*Y+c2*U
        if nq%M or np_%M: okint=False; continue
        q,p=nq//M,np_//M
        if q==0: continue
        best=(q,p) if best is None or abs(q)<abs(best[0]) else best
    # baseline: no T
    bb1,bb2=gauss(*cong_basis(Y,U,S))
    bq=None
    for (c1,c2) in [bb1,bb2]:
        nq=c1*X+c2*V; np_=c1*Y+c2*U
        if nq%S: continue
        q,p=nq//S,np_//S
        if q: bq=(q,p) if bq is None or abs(q)<abs(bq[0]) else bq
    print(f" {n:3d} | {vh/n:7.3f} | {log2(S)/n:9.3f} | {log2(M)/n:9.3f} | {log(M)/n:10.3f} |"
          f" {'y' if bq else '-':^9} | {'y' if best else '-':^11} | {okint}")
