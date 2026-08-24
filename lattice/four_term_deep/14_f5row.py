#!/usr/bin/env python3
"""14_f5row.py -- full invariants of one five-term row: exact integrality, the
quartic, exponents, companion, sharp k, Apery limit."""
import sys
from fractions import Fraction as Fr
from math import gcd
import mpmath as mp
mp.mp.dps=140

def run(RN,RD,M,J1,J2,a,c,d,f,g,j,C, NEXACT=120, NK=60):
    rho=Fr(RN,RD)
    b=int(Fr(a)*(1-rho)); e=int(-2*rho*Fr(d)); h=int(-(1+3*rho)*Fr(g))
    k=C*M*M; l=-C*M*(J1+J2); m=C*J1*J2
    P=lambda n:a*n*n+b*n+c; Q=lambda n:d*n*n+e*n+f
    R=lambda n:g*n*n+h*n+j;  T=lambda n:k*n*n+l*n+m
    print("P(n) = %d n^2 + %d n + %d" % (a,b,c))
    print("Q(n) = %d n^2 + %d n + %d" % (d,e,f))
    print("R(n) = %d n^2 + %d n + %d" % (g,h,j))
    print("T(n) = %d n^2 + %d n + %d" % (k,l,m))
    def seq(start,N):
        off=3
        v=[Fr(0)]*3+[Fr(1)] if start==0 else [Fr(0)]*4+[Fr(1)]
        n0=0 if start==0 else 1
        for n in range(n0,N):
            while len(v)<n+off+1: v.append(Fr(0))
            v.append((P(n)*v[n+off]-Q(n)*v[n-1+off]+R(n)*v[n-2+off]-T(n)*v[n-3+off])/Fr((n+1)**2))
        return v[off:]
    u=seq(0,NEXACT)
    print("u_n =",[int(x) if x.denominator==1 else str(x) for x in u[:12]])
    print("integral to n=%d :"%NEXACT, all(x.denominator==1 for x in u))
    lam=mp.polyroots([1,-a,d,-g,k],maxsteps=500,extraprec=500)
    lam=sorted(lam,key=lambda z:-abs(z))
    print("lambda =",[mp.nstr(z,14) for z in lam])
    print("|lambda| =",[mp.nstr(abs(z),10) for z in lam])
    bb=seq(1,max(NK,60))
    L=1; ks=[]
    for n in range(1,NK+1):
        L=L*n//gcd(L,n); den=bb[n].denominator; kk=0
        while den>1:
            gg=gcd(den,L)
            if gg==1: kk=99; break
            den//=gg; kk+=1
        ks.append(kk)
    kd=max(ks); print("sharp k =",kd," score =",float(-mp.log(abs(lam[1]))-kd))
    l1,l2=abs(lam[0]),abs(lam[1])
    if l2<l1*(1-1e-12):
        rate=mp.log10(l1/l2); N=int(min(60000,max(400,150/rate+200)))
        with mp.workdps(mp.mp.dps+80):
            A=[mp.mpf(0)]*3+[mp.mpf(1)]; B=[mp.mpf(0)]*4+[mp.mpf(1)]
            for n in range(0,N):
                pp,qq,rr,tt=mp.mpf(P(n)),mp.mpf(Q(n)),mp.mpf(R(n)),mp.mpf(T(n))
                A.append((pp*A[n+3]-qq*A[n+2]+rr*A[n+1]-tt*A[n])/mp.mpf((n+1)**2))
                if n>=1: B.append((pp*B[n+3]-qq*B[n+2]+rr*B[n+1]-tt*B[n])/mp.mpf((n+1)**2))
            x1=B[N+3]/A[N+3]; x2=B[N+2]/A[N+2]
            print("n =",N," xi =",mp.nstr(+x1,110))
            print("convergence log10|xi_n - xi_{n-1}| =",float(mp.log10(abs(x1-x2))) if x1!=x2 else 'exact')
            return mp.nstr(+x1,100)
    else:
        print("no archimedean Apery limit (|lam_2| = |lam_1|)")
    return None

if __name__=='__main__':
    v=list(map(int,sys.argv[1:13]))
    xi=run(*v)
    if xi and len(sys.argv)>13:
        open(sys.argv[13],'a').write("f5_%s %s\n"%("_".join(map(str,v)).replace("-","m"),xi))
