#!/usr/bin/env python3
"""Exact verification and invariants for four-term scan hits.

per row (rho=RN/RD ; M, j1, j2 ; a, c, d, f, C):
   P = a n^2 + a(1-rho) n + c ,  Q = d n^2 - 2 rho d n + f ,  R = C(Mn-j1)(Mn-j2)
   * exact integrality of u_n to n = NEXACT
   * cubic 1 - a t + d t^2 - g t^3  (g = C M^2): roots t_i, char roots lambda_i = 1/t_i
   * apparent-singularity (no-log) test at each t_i with rho in Z and at infinity if delta in Z
   * companion b_n (b_0=0,b_1=1), sharp k, Apery limit xi, score log(1/|lambda_2|) - k
   * primitivity under u_n -> mu^n u_n
"""
import sys, json, math
from fractions import Fraction as Fr
from math import gcd
import mpmath as mp

mp.mp.dps = 120

def coeffs(RN,RD,M,j1,j2,a,c,d,f,C):
    rho = Fr(RN,RD)
    b = Fr(a)*(1-rho); e = -2*rho*Fr(d)
    assert b.denominator==1 and e.denominator==1, (a,d,rho)
    g = C*M*M; h = -C*M*(j1+j2); jj = C*j1*j2
    return dict(a=a,b=int(b),c=c,d=d,e=int(e),f=f,g=g,h=h,j=jj,rho=rho,M=M,j1=j1,j2=j2,C=C)

def seq(co, N, u0=1, u1=None, start=0):
    """u_0=1 analytic solution (start=0) or companion b_0=0,b_1=1 (start=1)."""
    a,b,c,d,e,f,g,h,jj = (co[k] for k in 'abcdefghj')
    P=lambda n:a*n*n+b*n+c; Q=lambda n:d*n*n+e*n+f; R=lambda n:g*n*n+h*n+jj
    if start==0:
        u=[Fr(1)]; um1=Fr(0); um2=Fr(0)
        seqv=[Fr(0),Fr(0),Fr(1)]          # u_{-2},u_{-1},u_0
    else:
        seqv=[Fr(0),Fr(0),Fr(0),Fr(1)]    # u_{-2},u_{-1},u_0=0,u_1=1
    off = 2
    v = list(seqv)
    n0 = 0 if start==0 else 1
    for n in range(n0, N):
        nx = (P(n)*v[n+off] - Q(n)*v[n-1+off] + R(n)*v[n-2+off])/Fr((n+1)**2)
        v.append(nx)
    return v[off:]

def integral(u):
    return all(x.denominator==1 for x in u)

def taylor_shift(poly, t0):
    """poly given low->high in t; return coefficients at x = t - t0 (mpmath)."""
    n = len(poly)-1
    out = [mp.mpf(0)]*(n+1)
    cur = list(poly)
    # Horner-based synthetic division
    res = []
    cur = [mp.mpmathify(x) for x in poly]
    work = cur[:]
    for k in range(n+1):
        # evaluate remainder
        r = work[-1]
        newc = [work[-1]]
        for i in range(len(work)-2, -1, -1):
            r = work[i] + r*t0
            newc.append(r)
        newc.reverse()          # newc[0] = value at t0 = coefficient
        res.append(newc[0])
        work = newc[1:]
        if not work: break
    while len(res) < n+1: res.append(mp.mpf(0))
    return res

def frob_obstruction(A,B,C,nu,off,mstar):
    """Return the obstruction S at m = mstar (0 => apparent, no log)."""
    cs = [mp.mpmathify(1)]
    for m in range(1, mstar+1):
        r = m + off
        S = mp.mpmathify(0); I = mp.mpmathify(0)
        for j,Aj in enumerate(A):
            k = r+2-j
            if k < 0: continue
            w = Aj*(k+nu)*(k+nu-1)
            if k == m: I += w
            elif k < m: S += w*cs[k]
        for j,Bj in enumerate(B):
            k = r+1-j
            if k < 0: continue
            w = Bj*(k+nu)
            if k == m: I += w
            elif k < m: S += w*cs[k]
        for j,Cj in enumerate(C):
            k = r-j
            if k < 0: continue
            if k == m: I += Cj
            elif k < m: S += Cj*cs[k]
        if m == mstar:
            return S, I
        cs.append(-S/I)
    return None, None

def analyse(RN,RD,M,j1,j2,a,c,d,f,C, NEXACT=120, NK=60, NXI=600, full=True):
    co = coeffs(RN,RD,M,j1,j2,a,c,d,f,C)
    u = seq(co, NEXACT)
    if not integral(u): return None
    a_,b_,c_,d_,e_,f_,g_,h_,j_ = (co[k] for k in 'abcdefghj')
    rho = co['rho']
    # cubic  Rc(t) = 1 - a t + d t^2 - g t^3   <->  lam^3 - a lam^2 + d lam - g
    lams = mp.polyroots([1,-a_,d_,-g_], maxsteps=200, extraprec=200)
    lams = sorted(lams, key=lambda z: -abs(z))
    disc = 18*a_*d_*g_ - 4*a_**3*g_ + a_*a_*d_*d_ - 4*d_**3 - 27*g_*g_
    # exponents
    s1, s2 = Fr(j1,M), Fr(j2,M)
    delta = abs(s2-s1)
    # ---- apparent-singularity tests -------------------------------------
    Apoly = [0,0,1,-a_,d_,-g_]                       # t^2 Rc(t)
    Bpoly = [0,1,-(a_+b_),3*d_+e_,-(5*g_+h_)]        # t Sc(t)
    Cpoly = [0,-c_,d_+e_+f_,-(4*g_+2*h_+j_)]         # t Vc(t)
    apparent = []
    if rho.denominator == 1 and rho != 0:
        m = abs(int(rho)); nu = min(0, int(rho))
        for t0 in [1/L for L in lams]:
            A = taylor_shift(Apoly, t0); B = taylor_shift(Bpoly, t0); Cc = taylor_shift(Cpoly, t0)
            S,I = frob_obstruction(A,B,Cc,nu,-1,m)
            sc = max(abs(x) for x in A+B+Cc if x!=0)
            apparent.append(bool(abs(S) < mp.mpf(10)**(-40)*max(mp.mpf(1),sc)))
    inf_apparent = None
    if delta.denominator == 1 and delta != 0:
        At = [0,0,-g_,d_,-a_,1]                      # w^2 * Atil(w), Atil = w^3-a w^2+d w-g
        St = [-(5*g_+h_), 3*d_+e_, -(a_+b_), 1]
        Bt = [0] + [2*At[i+2]-St[i] for i in range(len(St))]   # w(2Atil - Stil)
        Vt = [-(4*g_+2*h_+j_), d_+e_+f_, -c_]
        nu = min(2-s1, 2-s2); nu = mp.mpmathify(float(nu))
        S,I = frob_obstruction([mp.mpmathify(x) for x in At],[mp.mpmathify(x) for x in Bt],
                               [mp.mpmathify(x) for x in Vt], nu, 0, int(delta))
        sc = max(abs(mp.mpmathify(x)) for x in At+Bt+Vt if x!=0)
        inf_apparent = bool(abs(S) < mp.mpf(10)**(-40)*max(mp.mpf(1),sc))
    if not full:
        return dict(cls=[RN,RD,M,j1,j2], row=[a,c,d,f,C], rho=str(rho), delta=str(delta),
                    lam=[[float(mp.re(z)),float(mp.im(z))] for z in lams], disc=int(disc),
                    apparent=apparent, inf_apparent=inf_apparent,
                    u=[int(x) for x in u[:8]])
    # ---- companion, k, xi ------------------------------------------------
    bb = seq(co, max(NK,NXI), start=1)
    dn = 1; k = 0
    for n in range(1, NK+1):
        dn = dn*(n//gcd(dn,n)) if False else dn
    # lcm(1..n)
    L = 1; ks = []
    for n in range(1, NK+1):
        L = L*n//gcd(L,n)
        den = bb[n].denominator
        kk = 0
        while den > 1:
            gg = gcd(den, L)
            if gg == 1: kk = 99; break
            den //= gg; kk += 1
        ks.append(kk)
    k = max(ks)
    # ---- Apery limit xi by high-precision forward recursion ---------------
    xi = None; conv = None; nxi = None
    l1, l2m = abs(lams[0]), abs(lams[1])
    if l2m < l1*(1-1e-12) and l1 > 0:
        rate = mp.log10(l1/l2m)
        nxi = int(min(20000, max(300, 90/rate + 60)))
        with mp.workdps(int(mp.mp.dps)+40):
            A = [mp.mpf(0), mp.mpf(0), mp.mpf(1)]     # a_{-2},a_{-1},a_0
            B = [mp.mpf(0), mp.mpf(0), mp.mpf(0), mp.mpf(1)]
            aP = [mp.mpf(x) for x in (a_,b_,c_)]; aQ=[mp.mpf(x) for x in (d_,e_,f_)]
            aR = [mp.mpf(x) for x in (g_,h_,j_)]
            for n in range(0, nxi):
                P = aP[0]*n*n+aP[1]*n+aP[2]; Q = aQ[0]*n*n+aQ[1]*n+aQ[2]; Rr = aR[0]*n*n+aR[1]*n+aR[2]
                A.append((P*A[n+2]-Q*A[n+1]+Rr*A[n])/mp.mpf((n+1)**2))
                if n >= 1:
                    B.append((P*B[n+2]-Q*B[n+1]+Rr*B[n])/mp.mpf((n+1)**2))
            r1 = B[nxi+2]/A[nxi+2]; r2 = B[nxi+1]/A[nxi+1]
            xi = +r1; conv = abs(r1-r2)
    l2 = abs(lams[1])
    score = (-mp.log(l2) - k) if l2>0 else None
    return dict(cls=[RN,RD,M,j1,j2], row=[a,c,d,f,C], rho=str(rho), delta=str(delta),
                lam=[[float(mp.re(z)),float(mp.im(z))] for z in lams], disc=int(disc),
                apparent=apparent, inf_apparent=inf_apparent,
                k=k, score=(float(score) if score is not None else None),
                xi=(mp.nstr(xi,70) if xi is not None else None), nxi=nxi,
                conv=(float(mp.log10(conv)) if conv and conv>0 else None),
                u=[int(x) for x in u[:8]])

if __name__ == '__main__':
    full = ('full' in sys.argv)
    if not full: mp.mp.dps = 60
    for line in sys.stdin:
        p = line.split()
        if not p or p[0].startswith('#'): continue
        v = list(map(int,p))
        try:
            r = analyse(*v, full=full)
        except Exception as ex:
            print(json.dumps({'row':v,'ERR':str(ex)})); continue
        print(json.dumps(r) if r else json.dumps({'row':v,'FAIL':'integrality'}))
