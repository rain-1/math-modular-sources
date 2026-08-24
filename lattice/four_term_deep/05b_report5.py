#!/usr/bin/env python3
"""Exact verification + invariants for FIVE-term (six singular point) hits.
row line: RN RD M J1 J2 a c d f g j C
  P=a n^2+(1-rho)a n+c, Q=d n^2-2 rho d n+f, R=g n^2-(1+3rho)g n+j,
  T=C(Mn-J1)(Mn-J2)   [k=CM^2, l=-CM(J1+J2), m=C J1 J2, (J1+J2)/M=2+4rho]
"""
import sys, json
from fractions import Fraction as Fr
from math import gcd
import mpmath as mp
mp.mp.dps = 80

def coeffs(RN,RD,M,J1,J2,a,c,d,f,g,j,C):
    rho=Fr(RN,RD)
    b=Fr(a)*(1-rho); e=-2*rho*Fr(d); h=-(1+3*rho)*Fr(g)
    k=C*M*M; l=-C*M*(J1+J2); m=C*J1*J2
    assert all(x.denominator==1 for x in (b,e,h)), (a,d,g,rho)
    return dict(a=a,b=int(b),c=c,d=d,e=int(e),f=f,g=g,h=int(h),j=j,k=k,l=l,m=m,rho=rho,M=M,J1=J1,J2=J2,C=C)

def seq(co,N,start=0):
    a,b,c,d,e,f,g,h,j,k,l,m=(co[x] for x in ['a','b','c','d','e','f','g','h','j','k','l','m'])
    P=lambda n:a*n*n+b*n+c; Q=lambda n:d*n*n+e*n+f; R=lambda n:g*n*n+h*n+j; T=lambda n:k*n*n+l*n+m
    off=3
    if start==0: v=[Fr(0)]*3+[Fr(1)]
    else:        v=[Fr(0)]*4+[Fr(1)]
    n0=0 if start==0 else 1
    for n in range(n0,N):
        while len(v) < n+off+1: v.append(Fr(0))
        nx=(P(n)*v[n+off]-Q(n)*v[n-1+off]+R(n)*v[n-2+off]-T(n)*v[n-3+off])/Fr((n+1)**2)
        v.append(nx)
    return v[off:]

def taylor_shift(poly,t0):
    work=[mp.mpmathify(x) for x in poly]; res=[]
    n=len(poly)-1
    for _ in range(n+1):
        r=work[-1]; newc=[work[-1]]
        for i in range(len(work)-2,-1,-1):
            r=work[i]+r*t0; newc.append(r)
        newc.reverse(); res.append(newc[0]); work=newc[1:]
        if not work: break
    while len(res)<n+1: res.append(mp.mpf(0))
    return res

def frob_obstruction(A,B,Cc,nu,off,mstar):
    cs=[mp.mpmathify(1)]
    for mm in range(1,mstar+1):
        r=mm+off; S=mp.mpmathify(0); I=mp.mpmathify(0)
        for jj,Aj in enumerate(A):
            kk=r+2-jj
            if kk<0: continue
            w=Aj*(kk+nu)*(kk+nu-1)
            if kk==mm: I+=w
            elif kk<mm: S+=w*cs[kk]
        for jj,Bj in enumerate(B):
            kk=r+1-jj
            if kk<0: continue
            w=Bj*(kk+nu)
            if kk==mm: I+=w
            elif kk<mm: S+=w*cs[kk]
        for jj,Cj in enumerate(Cc):
            kk=r-jj
            if kk<0: continue
            if kk==mm: I+=Cj
            elif kk<mm: S+=Cj*cs[kk]
        if mm==mstar: return S,I
        cs.append(-S/I)
    return None,None

def analyse(RN,RD,M,J1,J2,a,c,d,f,g,j,C, NEXACT=60):
    co=coeffs(RN,RD,M,J1,J2,a,c,d,f,g,j,C)
    u=seq(co,NEXACT)
    if not all(x.denominator==1 for x in u): return None
    a_,b_,c_,d_,e_,f_,g_,h_,j_,k_,l_,m_=(co[x] for x in ['a','b','c','d','e','f','g','h','j','k','l','m'])
    rho=co['rho']
    lams=mp.polyroots([1,-a_,d_,-g_,k_],maxsteps=400,extraprec=400)
    lams=sorted(lams,key=lambda z:-abs(z))
    s1,s2=Fr(J1,M),Fr(J2,M); delta=abs(s2-s1)
    # Rc = 1 - a t + d t^2 - g t^3 + k t^4 ; Sc ; Vc
    Apoly=[0,0,1,-a_,d_,-g_,k_]                        # t^2 Rc
    Bpoly=[0,1,-(a_+b_),3*d_+e_,-(5*g_+h_),7*k_+l_]    # t Sc
    Cpoly=[0,-c_,d_+e_+f_,-(4*g_+2*h_+j_),9*k_+3*l_+m_]# t Vc
    apparent=[]
    if rho.denominator==1 and rho!=0:
        mm=abs(int(rho)); nu=min(0,int(rho))
        for L in lams:
            t0=1/L
            A=taylor_shift(Apoly,t0);B=taylor_shift(Bpoly,t0);Cc=taylor_shift(Cpoly,t0)
            S,I=frob_obstruction(A,B,Cc,nu,-1,mm)
            sc=max(abs(x) for x in A+B+Cc if x!=0)
            apparent.append(bool(abs(S)<mp.mpf(10)**(-30)*max(mp.mpf(1),sc)))
    # discriminant of the quartic: repeated root?
    dmin=min(abs(lams[i]-lams[jx]) for i in range(4) for jx in range(i+1,4))
    return dict(cls=[RN,RD,M,J1,J2],row=[a,c,d,f,g,j,C],rho=str(rho),delta=str(delta),
                lam=[[float(mp.re(z)),float(mp.im(z))] for z in lams],
                sep=float(dmin), apparent=apparent, u=[int(x) for x in u[:9]])

if __name__=='__main__':
    for line in sys.stdin:
        p=line.split()
        if not p or p[0].startswith('#'): continue
        v=list(map(int,p))
        try: r=analyse(*v)
        except Exception as ex:
            print(json.dumps({'row':v,'ERR':str(ex)})); continue
        print(json.dumps(r) if r else json.dumps({'row':v,'FAIL':'integrality'}))
