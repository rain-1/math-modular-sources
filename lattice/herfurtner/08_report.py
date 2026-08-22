#!/usr/bin/env python3
"""Full post-processing: scan hits -> verified rows -> Herfurtner match -> table.

Filters applied, in order:
  (F0) integrality of u_n re-verified exactly for n <= NVER (default 300);
  (F1) Casoratian non-degeneracy: Q(n) != 0 for all n >= 1
       (otherwise b_n is a rational multiple of a_n and the row is vacuous);
  (F2) lam_1 != lam_2 recorded (double root -> no second growth rate);
  (F3) deduplication modulo t -> t/c, i.e. (A,B,C) ~ (cA,cB,c^2C), and modulo
       t -> -t, i.e. (A,B,C) ~ (-A,-B,C) with u_n -> (-1)^n u_n; a row is kept
       once, in its smallest representative.
"""
import sys, os, glob, json, math
from fractions import Fraction
from math import gcd
from sympy import Rational as Q, simplify, sqrt, N as num
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from importlib import import_module
H = import_module('04_herfurtner')

NVER = int(os.environ.get('NVER','300'))
NK   = int(os.environ.get('NK','60'))

def dnl(N):
    D=[1]*(N+2)
    for n in range(1,N+2): D[n]=D[n-1]*n//gcd(D[n-1],n)
    return D
DN = dnl(max(NVER,NK))

def seq(A,be,B,M,j1,j2,C,N,u0,u1):
    de=C*M*M; ep=-C*M*(j1+j2); ze=C*j1*j2
    u=[Fraction(u0),Fraction(u1)]
    for n in range(1,N):
        P=A*n*n+be*n+B; Qn=de*n*n+ep*n+ze
        u.append((P*u[n]-Qn*u[n-1])/Fraction((n+1)**2))
    return u

def canon(M,j1,j2,A,B,C):
    """smallest representative in the orbit (A,B,C)~(cA,cB,c^2C), c>=1, and
       (A,B,C)~(-A,-B,C)."""
    best=(A,B,C)
    g=gcd(abs(A),abs(B)) if (A or B) else 0
    if g:
        d=2
        while d*d<=abs(C) or d<=g:
            if d>g: break
            if g%d==0 and C%(d*d)==0:
                A,B,C = A//d,B//d,C//(d*d); g=gcd(abs(A),abs(B))
                d=2; continue
            d+=1
    if A<0 or (A==0 and B<0): A,B=-A,-B
    return (A,B,C)

def main():
    files=sys.argv[1:] or (glob.glob('out/c_*.txt')+['out/a0.txt'])
    raw=set()
    for p in files:
        if not os.path.exists(p): continue
        for line in open(p):
            f=line.split()
            if len(f)!=6: continue
            raw.add(tuple(map(int,f)))
    print("# raw hits:",len(raw))
    keep = {k:0 for k in raw}
    T=H.table()
    out=[]
    verified=set()
    for (M,j1,j2,A,B,C) in sorted(keep):
        S=2*M-j1-j2
        if (A*S)%(2*M): continue
        be=A*S//(2*M)
        uu=seq(A,be,B,M,j1,j2,C,NVER,1,B)
        if all(x.denominator==1 for x in uu): verified.add((M,j1,j2,A,B,C))
    print("# integral to n=%d: %d"%(NVER,len(verified)))
    # a row is primitive if it is not c*(a smaller verified row)
    prim=set()
    for (M,j1,j2,A,B,C) in verified:
        g=gcd(abs(A),abs(B)) if (A or B) else 0
        red=False
        for d in range(2, (abs(g) or 1)+1):
            if g%d: continue
            if C%(d*d): continue
            if (M,j1,j2,A//d,B//d,C//(d*d)) in verified: red=True; break
        if not red: prim.add((M,j1,j2,A,B,C))
    print("# primitive (not a t->t/c rescaling of another): %d"%len(prim))
    for (M,j1,j2,A,B,C) in sorted(prim):
        S=2*M-j1-j2
        be=A*S//(2*M); D=C*M*M
        # F1 Casoratian
        deg = (j1%M==0 and j1//M>=1) or (j2%M==0 and j2//M>=1)
        u=seq(A,be,B,M,j1,j2,C,NVER,1,B)
        if any(x.denominator!=1 for x in u): continue
        b=seq(A,be,B,M,j1,j2,C,NK,0,1)
        k=0
        while k<=8:
            if all((DN[n]**k*b[n]).denominator==1 for n in range(NK+1)): break
            k+=1
        disc=A*A-4*D
        if disc<0: l1=l2=math.sqrt(abs(D)); cplx=1
        else:
            cplx=0
            x1=(A+math.sqrt(disc))/2; x2=(A-math.sqrt(disc))/2
            l1,l2=(x1,x2) if abs(x1)>=abs(x2) else (x2,x1)
        score = (math.log(1/abs(l2))-k) if (not cplx and abs(l2)>1e-12) else None
        Iv,(tt,ti),hits = H.match_row(M,j1,j2,A,C,T) if D else (None,(None,None),[])
        out.append(dict(M=M,j1=j1,j2=j2,A=A,B=B,C=C,be=be,D=D,disc=disc,
                        casdeg=deg,dbl=(disc==0),cplx=cplx,lam1=l1,lam2=l2,k=k,
                        score=score,I=str(Iv),types=[str(tt),str(ti)],
                        u=[int(x) for x in u[:9]],
                        herf=[dict(row=r['row'],fib=" ".join(r['fib']),degJ=r['degJ'],
                                   at0=r['at0'],atinf=r['atinf']) for r in hits]))
    json.dump(out,open('out/rows_full.json','w'),indent=0)
    print("# verified integral rows (n <= %d): %d"%(NVER,len(out)))
    nd=[r for r in out if not r['casdeg'] and not r['dbl']]
    print("# non-degenerate (Casoratian != 0, lam_1 != lam_2): %d"%len(nd))
    print("# with a rigid Herfurtner match: %d"%len([r for r in nd if r['herf']]))

main()
