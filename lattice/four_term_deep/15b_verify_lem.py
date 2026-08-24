#!/usr/bin/env python3
"""15b_verify_lem.py -- independent verification of the gauge-partner (rho_p=+1/2)
rows whose Apery limits hit Gamma(1/4)^4/(64 Pi)."""
from fractions import Fraction as Fr
from math import gcd
import mpmath as mp, sys
mp.mp.dps = 240
ROWS = [
 (Fr(1,2), Fr(0), 1,1,1,  8, 16,  4,  48,  0, -128),
 (Fr(1,2), Fr(0), 1,1,1,  2, 14,  2,  28,  4,    8),
 (Fr(1,2), Fr(0), 1,1,1, -2,  6,  0, -32, -8,   32),
]
NEX=200; NK=80
def coeffs(rp, rr, M, J1, J2, r, a, c, d, f, C):
    g=C*M*M; assert r!=0 and g%r==0
    p=g//r; s=a-r; assert d==s*r+p, (d,s*r+p)
    b=(1-rr)*r+(1-rp)*s
    e=-rp*(2*p+r*s)-rr*r*s
    h=-(1+2*rp+rr)*Fr(g); j=Fr(C*J1*J2)
    for x in (b,e,h,j): assert x.denominator==1
    return a,int(b),c,d,int(e),f,g,int(h),int(j)
def rho_at(co, lam):
    a,b,c,d,e,f,g,h,j = co
    t=1/lam
    Rp=-a+2*d*t-3*g*t*t
    Tv=1-b*t+(d+e)*t*t-(2*g+h)*t**3
    return -Tv/(t*Rp)
out=open('out/lem_rows_xi.txt','w')
for (rp,rr,M,J1,J2,r,a,c,d,f,C) in ROWS:
    co=coeffs(rp,rr,M,J1,J2,r,a,c,d,f,C)
    a_,b_,c_,d_,e_,f_,g_,h_,j_=co
    P=lambda n:a_*n*n+b_*n+c_; Q=lambda n:d_*n*n+e_*n+f_; R=lambda n:g_*n*n+h_*n+j_
    print("="*72); print("row (a,c,d,f,C)=(%d,%d,%d,%d,%d) r=%d"%(a,c,d,f,C,r))
    print("  P=%dn^2+%dn+%d  Q=%dn^2+%dn+%d  R=%dn^2+%dn+%d"%(a_,b_,c_,d_,e_,f_,g_,h_,j_))
    v=[Fr(0),Fr(0),Fr(1)]
    for n in range(NEX): v.append((P(n)*v[n+2]-Q(n)*v[n+1]+R(n)*v[n])/Fr((n+1)**2))
    u=v[2:]
    print("  u:",[int(x) for x in u[:8]],"  integral to n=%d: %s"%(NEX,all(x.denominator==1 for x in u)))
    lam=sorted(mp.polyroots([1,-a_,d_,-g_],maxsteps=400,extraprec=400),key=lambda z:-abs(z))
    print("  lambda:",[mp.nstr(z,16) for z in lam])
    print("  rho_i :",[mp.nstr(rho_at(co,z),10) for z in lam])
    w=[Fr(0),Fr(0),Fr(0),Fr(1)]
    for n in range(1,NK+1):
        while len(w)<n+3: w.append(Fr(0))
        w.append((P(n)*w[n+2]-Q(n)*w[n+1]+R(n)*w[n])/Fr((n+1)**2))
    bb=w[2:]; L=1; ks=[]
    for n in range(1,NK+1):
        L=L*n//gcd(L,n); den=bb[n].denominator; kk=0
        while den>1:
            gg=gcd(den,L)
            if gg==1: kk=99;break
            den//=gg; kk+=1
        ks.append(kk)
    print("  sharp k =",max(ks))
    rate=mp.log10(abs(lam[0])/abs(lam[1])); N=int(min(200000,max(600,165/rate+400)))
    with mp.workdps(mp.mp.dps+120):
        A=[mp.mpf(0),mp.mpf(0),mp.mpf(1)]; B=[mp.mpf(0),mp.mpf(0),mp.mpf(0),mp.mpf(1)]
        for n in range(N):
            pp,qq,rr2=mp.mpf(P(n)),mp.mpf(Q(n)),mp.mpf(R(n))
            A.append((pp*A[n+2]-qq*A[n+1]+rr2*A[n])/mp.mpf((n+1)**2))
            if n>=1: B.append((pp*B[n+2]-qq*B[n+1]+rr2*B[n])/mp.mpf((n+1)**2))
        xi=+(B[N+2]/A[N+2])
    tgt=mp.gamma(mp.mpf(1)/4)**4/(64*mp.pi)
    print("  n=%d  xi = %s"%(N, mp.nstr(xi,130)))
    print("  Gamma(1/4)^4/(64 Pi) = %s"%mp.nstr(tgt,130))
    print("  |xi - target| = %s   agreeing digits = %d"%(mp.nstr(abs(xi-tgt),8), int(-mp.log10(abs(xi-tgt)/abs(xi)))))
    out.write("lem_%d_%d_%d_%d_%d_r%s %s\n"%(a,c,d,f,C,str(r).replace('-','m'),mp.nstr(xi,130)))
    sys.stdout.flush()
out.close()
