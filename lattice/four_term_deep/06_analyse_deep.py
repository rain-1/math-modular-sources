#!/usr/bin/env python3
"""Deduplication / invariants for the DEEP four-term scan (equal-rho and mixed-rho).

Input lines:
  equal-rho (10 fields):  RN RD M J1 J2 a c d f C
  mixed-rho (13 fields):  RPN RPD RRN RRD M J1 J2 r a c d f C
Mixed-rho coefficients (see 03_fmix.c / 01_mixcheck.py):
  s = a - r, p = g/r, g = C M^2,
  b = (1-rho_r) r + (1-rho_p) s,
  e = -rho_p (2p + r s) - rho_r r s,
  h = -(1 + 2 rho_p + rho_r) g,  j = C J1 J2.

Verdicts: NOT-INTEGRAL / REPEATED (disc 0) / DISGUISED (apparent singularity =
four-point projective local system) / RESCALED / CANDIDATE.
"""
import sys, json, glob, os
from fractions import Fraction as Fr
from math import gcd
import mpmath as mp
mp.mp.dps = 120

def coeffs_eq(RN,RD,M,j1,j2,a,c,d,f,C):
    rho=Fr(RN,RD); b=Fr(a)*(1-rho); e=-2*rho*Fr(d)
    assert b.denominator==1 and e.denominator==1
    g=C*M*M; h=-C*M*(j1+j2); jj=C*j1*j2
    return dict(a=a,b=int(b),c=c,d=d,e=int(e),f=f,g=g,h=h,j=jj,
                rhos=[rho,rho,rho],delta=abs(Fr(j2-j1,M)),M=M,j1=j1,j2=j2,C=C)

def coeffs_mix(RPN,RPD,RRN,RRD,M,j1,j2,r,a,c,d,f,C):
    rp=Fr(RPN,RPD); rr=Fr(RRN,RRD)
    g=C*M*M
    assert r!=0 and g % r == 0
    p=g//r; s=a-r
    b=(1-rr)*r+(1-rp)*s
    e=-rp*(2*p+r*s)-rr*r*s
    h=-(1+2*rp+rr)*Fr(g)
    jj=Fr(C*j1*j2)
    assert d == s*r+p, (d, s*r+p)
    for x in (b,e,h,jj): assert x.denominator==1, (b,e,h,jj)
    return dict(a=a,b=int(b),c=c,d=d,e=int(e),f=f,g=g,h=int(h),j=int(jj),
                rhos=[rr,rp,rp],delta=abs(Fr(j2-j1,M)),M=M,j1=j1,j2=j2,C=C,r=r)

def seq(co,N,start=0):
    a,b,c,d,e,f,g,h,jj=(co[k] for k in 'abcdefghj')
    P=lambda n:a*n*n+b*n+c; Q=lambda n:d*n*n+e*n+f; R=lambda n:g*n*n+h*n+jj
    off=2
    v=[Fr(0),Fr(0),Fr(1)] if start==0 else [Fr(0),Fr(0),Fr(0),Fr(1)]
    n0=0 if start==0 else 1
    for n in range(n0,N):
        v.append((P(n)*v[n+off]-Q(n)*v[n-1+off]+R(n)*v[n-2+off])/Fr((n+1)**2))
    return v[off:]

def taylor_shift(poly,t0):
    work=[mp.mpmathify(x) for x in poly]; res=[]; n=len(poly)-1
    for _ in range(n+1):
        rr=work[-1]; newc=[work[-1]]
        for i in range(len(work)-2,-1,-1):
            rr=work[i]+rr*t0; newc.append(rr)
        newc.reverse(); res.append(newc[0]); work=newc[1:]
        if not work: break
    while len(res)<n+1: res.append(mp.mpf(0))
    return res

def frob_obstruction(A,B,Cc,nu,off,mstar):
    cs=[mp.mpmathify(1)]
    for m in range(1,mstar+1):
        rr=m+off; S=mp.mpmathify(0); I=mp.mpmathify(0)
        for jx,Aj in enumerate(A):
            k=rr+2-jx
            if k<0: continue
            w=Aj*(k+nu)*(k+nu-1)
            if k==m: I+=w
            elif k<m: S+=w*cs[k]
        for jx,Bj in enumerate(B):
            k=rr+1-jx
            if k<0: continue
            w=Bj*(k+nu)
            if k==m: I+=w
            elif k<m: S+=w*cs[k]
        for jx,Cj in enumerate(Cc):
            k=rr-jx
            if k<0: continue
            if k==m: I+=Cj
            elif k<m: S+=Cj*cs[k]
        if m==mstar: return S,I
        cs.append(-S/I)
    return None,None

def analyse(co, NEXACT=120, NK=60, full=True):
    u=seq(co,NEXACT)
    if not all(x.denominator==1 for x in u): return None
    a_,b_,c_,d_,e_,f_,g_,h_,j_=(co[k] for k in 'abcdefghj')
    lams=mp.polyroots([1,-a_,d_,-g_],maxsteps=300,extraprec=300)
    lams=sorted(lams,key=lambda z:-abs(z))
    disc=18*a_*d_*g_-4*a_**3*g_+a_*a_*d_*d_-4*d_**3-27*g_*g_
    delta=co['delta']
    Apoly=[0,0,1,-a_,d_,-g_]; Bpoly=[0,1,-(a_+b_),3*d_+e_,-(5*g_+h_)]
    Cpoly=[0,-c_,d_+e_+f_,-(4*g_+2*h_+j_)]
    # apparent test at each t_i whose rho is a non-zero integer.  Pair each root
    # with its rho: for mixed classes the rational root r carries rho_r.
    apparent=[]
    rhos=co['rhos']
    if 'r' in co:
        rr=co['r']
        idx=min(range(3),key=lambda i:abs(lams[i]-rr))
        rho_of=[rhos[1]]*3; rho_of[idx]=rhos[0]
    else:
        rho_of=[rhos[0]]*3
    for i,L in enumerate(lams):
        rh=rho_of[i]
        if rh.denominator==1 and rh!=0:
            m=abs(int(rh)); nu=min(0,int(rh)); t0=1/L
            A=taylor_shift(Apoly,t0);B=taylor_shift(Bpoly,t0);Cc=taylor_shift(Cpoly,t0)
            S,I=frob_obstruction(A,B,Cc,nu,-1,m)
            sc=max(abs(x) for x in A+B+Cc if x!=0)
            apparent.append(bool(abs(S)<mp.mpf(10)**(-40)*max(mp.mpf(1),sc)))
        else:
            apparent.append(False)
    inf_apparent=None
    s1,s2=Fr(co['j1'],co['M']),Fr(co['j2'],co['M'])
    if delta.denominator==1 and delta!=0:
        At=[0,0,-g_,d_,-a_,1]; St=[-(5*g_+h_),3*d_+e_,-(a_+b_),1]
        Bt=[0]+[2*At[i+2]-St[i] for i in range(len(St))]
        Vt=[-(4*g_+2*h_+j_),d_+e_+f_,-c_]
        nu=mp.mpmathify(float(min(2-s1,2-s2)))
        S,I=frob_obstruction([mp.mpmathify(x) for x in At],[mp.mpmathify(x) for x in Bt],
                             [mp.mpmathify(x) for x in Vt],nu,0,int(delta))
        sc=max(abs(mp.mpmathify(x)) for x in At+Bt+Vt if x!=0)
        inf_apparent=bool(abs(S)<mp.mpf(10)**(-40)*max(mp.mpf(1),sc))
    out=dict(lam=[[float(mp.re(z)),float(mp.im(z))] for z in lams],disc=int(disc),
             apparent=apparent,inf_apparent=inf_apparent,
             rhos=[str(x) for x in rhos],delta=str(delta),u=[int(x) for x in u[:9]])
    if not full: return out
    bb=seq(co,max(NK,300),start=1)
    L=1; ks=[]
    for n in range(1,NK+1):
        L=L*n//gcd(L,n); den=bb[n].denominator; kk=0
        while den>1:
            gg=gcd(den,L)
            if gg==1: kk=99; break
            den//=gg; kk+=1
        ks.append(kk)
    k=max(ks); out['k']=k
    l1,l2=abs(lams[0]),abs(lams[1])
    out['score']=float(-mp.log(l2)-k) if l2>0 else None
    out['real_l2']=bool(abs(mp.im(lams[1]))<1e-30)
    out['xi']=None; out['nxi']=None
    if l2<l1*(1-1e-12) and l1>0:
        rate=mp.log10(l1/l2); nxi=int(min(30000,max(300,95/rate+80)))
        with mp.workdps(mp.mp.dps+60):
            A=[mp.mpf(0),mp.mpf(0),mp.mpf(1)]; B=[mp.mpf(0),mp.mpf(0),mp.mpf(0),mp.mpf(1)]
            for n in range(0,nxi):
                P=a_*n*n+b_*n+c_; Q=d_*n*n+e_*n+f_; R=g_*n*n+h_*n+j_
                A.append((P*A[n+2]-Q*A[n+1]+R*A[n])/mp.mpf((n+1)**2))
                if n>=1: B.append((P*B[n+2]-Q*B[n+1]+R*B[n])/mp.mpf((n+1)**2))
            xi=B[nxi+2]/A[nxi+2]; xi2=B[nxi+1]/A[nxi+1]
            out['xi']=mp.nstr(+xi,72); out['nxi']=nxi
            out['conv']=float(mp.log10(abs(xi-xi2))) if xi!=xi2 else None
    return out

def parse(line):
    v=line.split()
    if len(v)==10:
        v=list(map(int,v)); return ('eq',v[:5],v[5:],coeffs_eq(*v))
    if len(v)==13:
        v=list(map(int,v)); return ('mix',v[:7]+[v[7]],v[8:],coeffs_mix(*v))
    return None

def reduce_scale(kind,cls,row,co):
    a,c,d,f,C=row[-5:] if kind=='eq' else row[-5:]
    best=1
    for mu in range(2,60):
        if a%mu or c%mu or d%mu**2 or f%mu**2 or C%mu**3: continue
        r2=[a//mu,c//mu,d//mu**2,f//mu**2,C//mu**3]
        if r2[4]==0: continue
        try:
            if kind=='eq': co2=coeffs_eq(*cls,*r2)
            else:
                rr=co['r']
                if rr%mu: continue
                co2=coeffs_mix(*cls[:7],rr//mu,*r2)
        except Exception: continue
        u=seq(co2,60)
        if all(x.denominator==1 for x in u): best=mu
    return best

def canon(v):
    """t -> -t : equal-rho (a,c,d,f,C) -> (-a,-c,d,f,-C); mixed also r -> -r.
    Normalise to a > 0 (or, when a = 0, to C > 0)."""
    v=list(v)
    if len(v)==10:
        a,c,d,f,C=v[5:]
        if a<0 or (a==0 and C<0): v[5:]=[-a,-c,d,f,-C]
    elif len(v)==13:
        r=v[7]; a,c,d,f,C=v[8:]
        if a<0 or (a==0 and C<0): v[7]=-r; v[8:]=[-a,-c,d,f,-C]
    return tuple(v)

def main(paths,out_json,full=False):
    seen=set(); recs=[]
    for p in paths:
        for line in open(p):
            if not line.strip() or line.startswith('#'): continue
            key=canon(tuple(int(x) for x in line.split()))
            if key in seen: continue
            seen.add(key); recs.append(" ".join(map(str,key)))
    sys.stderr.write("# %d distinct hits\n"%len(recs))
    out=[]
    for line in recs:
        pr=parse(line)
        if pr is None: continue
        kind,cls,row,co=pr
        try: r=analyse(co,full=full)
        except Exception as ex:
            out.append(dict(kind=kind,cls=cls,row=row,verdict='ERR',err=str(ex))); continue
        if r is None:
            out.append(dict(kind=kind,cls=cls,row=row,verdict='NOT-INTEGRAL')); continue
        r['kind']=kind; r['cls']=cls; r['row']=row
        if r['disc']==0: r['verdict']='REPEATED'
        elif any(r['apparent']) or r['inf_apparent']: r['verdict']='DISGUISED'
        else:
            mu=reduce_scale(kind,cls,row,co)
            r['verdict']='RESCALED' if mu>1 else 'CANDIDATE'; r['mu']=mu
        out.append(r)
    json.dump(out,open(out_json,'w'))
    from collections import Counter
    cnt=Counter((tuple(r['cls'][:-1]) if r['kind']=='mix' else tuple(r['cls']),r['verdict']) for r in out)
    for k in sorted(cnt,key=str): sys.stderr.write("%s %d\n"%(k,cnt[k]))
    sys.stderr.write("TOTAL %d  CANDIDATE %d\n"%(len(out),sum(1 for r in out if r['verdict']=='CANDIDATE')))
    return out

if __name__=='__main__':
    full = '--full' in sys.argv
    args=[x for x in sys.argv[1:] if x!='--full']
    main(args[1:],args[0],full=full)
