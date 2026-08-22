#!/usr/bin/env python3
"""Post-process the Herfurtner-class scan.

For each hit (M,j1,j2,A,B,C):
   P(n) = A n^2 + A(2M-j1-j2)/(2M) n + B,  Q(n) = C(Mn-j1)(Mn-j2)
compute  u_n (exact, n <= NVER),  lambda_{1,2} = roots of x^2 - A x + C M^2,
the companion b_n (b_0=0, b_1=1), k = min{k : d_n^k b_n in Z},
score = log(1/|lambda_2|) - k, budget = log|lambda_1| - k,
and flags: Casoratian degeneracy (Q(n)=0 for some n>=1), double root,
lambda_2 = 0 (order drops).
Dedupe: (A,B,C) ~ (cA,cB,c^2 C) for integer c>=1 (t -> t/c); keep the minimal one.
"""
import sys, os, math, json
from fractions import Fraction
from math import gcd

NVER = int(os.environ.get('NVER','40'))
NK   = int(os.environ.get('NK','40'))

def dnl(N):
    D=[1]*(N+2)
    for n in range(1,N+2): D[n]=D[n-1]*n//gcd(D[n-1],n)
    return D

def seq(A,be,B,C,M,j1,j2,N,u0,u1):
    de=C*M*M; ep=-C*M*(j1+j2); ze=C*j1*j2
    u=[Fraction(u0),Fraction(u1)]
    for n in range(1,N):
        P=A*n*n+be*n+B; Q=de*n*n+ep*n+ze
        u.append((P*u[n]-Q*u[n-1])/Fraction((n+1)**2))
    return u

def analyse(M,j1,j2,A,B,C,Nver=NVER,Nk=NK):
    S=2*M-j1-j2
    assert (A*S)%(2*M)==0
    be=A*S//(2*M)
    D=C*M*M
    r={'M':M,'j1':j1,'j2':j2,'A':A,'B':B,'C':C,'be':be,'D':D}
    # Casoratian degeneracy: Q(n)=0 for some integer n>=1  <=> j1 or j2 = M*n
    r['casdeg'] = (j1%M==0 and j1//M>=1) or (j2%M==0 and j2//M>=1)
    disc=A*A-4*D; r['disc']=disc
    if disc<0:
        r['cplx']=1; l1=l2=math.sqrt(abs(D))
    else:
        r['cplx']=0
        x1=(A+math.sqrt(disc))/2; x2=(A-math.sqrt(disc))/2
        l1,l2=(x1,x2) if abs(x1)>=abs(x2) else (x2,x1)
    r['lam1'],r['lam2']=l1,l2
    r['dblroot']= (disc==0)
    u=seq(A,be,B,C,M,j1,j2,Nver,1,B)
    r['integral']=all(x.denominator==1 for x in u)
    r['u']=[int(x) for x in u[:8]] if r['integral'] else None
    if not r['integral']: return r
    b=seq(A,be,B,C,M,j1,j2,Nk,0,1)
    Dn=dnl(Nk)
    k=0
    while k<=8:
        if all((Dn[n]**k*b[n]).denominator==1 for n in range(Nk+1)): break
        k+=1
    r['k']=k
    r['score']= (math.log(1/abs(l2))-k) if abs(l2)>1e-13 and not r['cplx'] else None
    r['budget']= math.log(abs(l1))-k if abs(l1)>0 else None
    return r

def canon(A,B,C,M):
    """smallest (A,B,C) in the ray (cA,cB,c^2C); returns (A,B,C,c)."""
    best=(A,B,C,1)
    g=gcd(abs(A),abs(B))
    c=1
    for d in range(2,int(abs(g))+1):
        if g%d: continue
        if C%(d*d)==0:
            best2=(A//d,B//d,C//(d*d),d)
            if abs(best2[0])<abs(best[0]): best=best2
    return best

if __name__=='__main__':
    rows={}
    for path in sys.argv[1:]:
        for line in open(path):
            p=line.split()
            if len(p)!=6: continue
            M,j1,j2,A,B,C=map(int,p)
            rows[(M,j1,j2,A,B,C)]=None
    out=[]
    for key in sorted(rows):
        M,j1,j2,A,B,C=key
        r=analyse(M,j1,j2,A,B,C)
        if not r['integral']: continue
        out.append(r)
    json.dump(out,open('out/rows.json','w'))
    print("verified integral rows:",len(out))
