#!/usr/bin/env python3
"""Enumerate eta-quotient pairs (t,F) on Gamma_0(N), N<=NMAX:

  t : weight 0, ord_inf(t) = 1, polar degree <= DEGMAX   (a Hauptmodul of
      Gamma_0(N) or of an Atkin-Lehner quotient, or a degree<=DEGMAX function)
  F : weight 2, HOLOMORPHIC, ord_inf(F) = 0             (so F = 1 + O(q))

Both are eta quotients, so div is supported on the cusps and F has no zero on H:
hypothesis (H3) of ROOT_ROWS Thm R4 holds automatically for g = sqrt(F).

Ligozat: for f = prod_{d|N} eta(d tau)^{r_d},
  ord_{a/c}(f) = N/(24 gcd(c^2,N)) * sum_d gcd(c,d)^2 r_d / d      (c | N)
and #cusps with denominator c is phi(gcd(c,N/c)).  Total divisor degree = k*mu/12.
"""
import sys, json
from fractions import Fraction
from math import gcd

def divisors(n): return [d for d in range(1,n+1) if n%d==0]
def phi(n):
    r=n; p=2; m=n
    while p*p<=m:
        if m%p==0:
            while m%p==0: m//=p
            r-=r//p
        p+=1
    if m>1: r-=r//m
    return r
def index0(N):
    mu=N; m=N; p=2
    ps=[]
    while p*p<=m:
        if m%p==0:
            ps.append(p)
            while m%p==0: m//=p
        p+=1
    if m>1: ps.append(m)
    for p in ps: mu=mu*(p+1)//p
    return mu

def ligozat(N):
    ds=divisors(N)
    k=len(ds)
    A=[[Fraction(N,24*gcd(c*c,N))*Fraction(gcd(c,d)**2,d) for d in ds] for c in ds]
    return ds,A

def matinv(A):
    n=len(A)
    M=[row[:]+[Fraction(1 if i==j else 0) for j in range(n)] for i,row in enumerate(A)]
    for i in range(n):
        piv=None
        for r in range(i,n):
            if M[r][i]!=0: piv=r;break
        M[i],M[piv]=M[piv],M[i]
        pv=M[i][i]
        M[i]=[x/pv for x in M[i]]
        for r in range(n):
            if r!=i and M[r][i]!=0:
                f=M[r][i]
                M[r]=[a-f*b for a,b in zip(M[r],M[i])]
    return [row[n:] for row in M]

def enum_ord(weights, total, free_idx, cap):
    """nonneg integer vectors x over free_idx with sum weights[i]*x[i] == total, x[i]<=cap"""
    res=[]
    m=len(free_idx)
    def rec(i, rem, cur):
        if i==m:
            if rem==0: res.append(cur[:])
            return
        w=weights[free_idx[i]]
        # remaining slots can contribute at most cap*sum(weights of rest)
        maxrest=sum(weights[free_idx[j]] for j in range(i+1,m))*cap
        v=0
        while v*w<=rem and v<=cap:
            if rem-v*w<=maxrest:
                cur.append(v); rec(i+1, rem-v*w, cur); cur.pop()
            v+=1
    rec(0,total,[])
    return res

def enum_t(weights, ds, Ainv, N, degmax):
    """ord vectors for weight-0 functions with ord_inf=1, polar degree<=degmax"""
    k=len(ds); inf=k-1   # cusp denominator c=N is the cusp infinity
    out=[]
    # positive part (excluding the forced 1 at infinity) and negative part
    # weighted positive total P, weighted negative total P (degree 0), P<=degmax
    idx=[i for i in range(k) if i!=inf]
    for P in range(1, degmax+1):
        # ord_inf = 1 contributes weight w_inf*1 to positive part
        winf=weights[inf]
        if winf>P: continue
        pos=enum_ord(weights, P-winf, idx, degmax)
        neg=enum_ord(weights, P, idx, degmax)
        for pv in pos:
            for nv in neg:
                ordv=[0]*k; ordv[inf]=1
                bad=False
                for j,i in enumerate(idx):
                    ordv[i]=pv[j]-nv[j]
                    if pv[j]>0 and nv[j]>0: bad=True;break
                if bad: continue
                out.append(ordv)
    return out

def r_from_ord(Ainv, ordv):
    k=len(ordv)
    r=[sum(Ainv[i][j]*ordv[j] for j in range(k)) for i in range(k)]
    if any(x.denominator!=1 for x in r): return None
    return [int(x) for x in r]

def run(NMAX=60, DEGMAX=4, CAP=8, MAXF=400000):
    data={}
    for N in range(2,NMAX+1):
        ds,A=ligozat(N); k=len(ds); Ainv=matinv(A)
        weights=[phi(gcd(c,N//c)) for c in ds]
        mu=index0(N)
        if (2*mu)%12: 
            Ftot=None
        else:
            Ftot=(2*mu)//12
        # ---- forms F: weight 2, holomorphic, ord_inf=0
        Fs=[]
        if Ftot is not None:
            inf=k-1
            idx=[i for i in range(k) if i!=inf]
            cap=min(CAP, Ftot)
            cand=enum_ord(weights, Ftot, idx, cap)
            if len(cand)>MAXF:
                cand=enum_ord(weights, Ftot, idx, min(3,cap))
            for cv in cand:
                ordv=[0]*k
                for j,i in enumerate(idx): ordv[i]=cv[j]
                r=r_from_ord(Ainv,ordv)
                if r is None: continue
                if sum(r)!=4: continue
                Fs.append(r)
        # ---- parameters t
        Ts=[]
        for ordv in enum_t(weights, ds, Ainv, N, DEGMAX):
            r=r_from_ord(Ainv,ordv)
            if r is None: continue
            if sum(r)!=0: continue
            Ts.append((r,ordv))
        # dedup
        Fs=[list(x) for x in dict.fromkeys(tuple(f) for f in Fs)]
        Tsd={}
        for r,o in Ts: Tsd[tuple(r)]=o
        Ts=[(list(r),o) for r,o in Tsd.items()]
        data[N]={"divs":ds,"weights":weights,"mu":mu,"F":Fs,
                 "t":[{"r":r,"ord":o,"deg":sum(weights[i]*max(o[i],0) for i in range(k))} for r,o in Ts]}
        print(f"N={N:3d} divs={k:2d} mu={mu:4d}  #t={len(Ts):5d}  #F={len(Fs):6d}  pairs={len(Ts)*len(Fs)}",flush=True)
    json.dump(data, open("eta_pairs.json","w"))
    tot=sum(len(v["t"])*len(v["F"]) for v in data.values())
    print("TOTAL PAIRS", tot)

if __name__=="__main__":
    import os
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    run(int(sys.argv[1]) if len(sys.argv)>1 else 60)
