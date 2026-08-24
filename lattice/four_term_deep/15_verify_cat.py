#!/usr/bin/env python3
"""15_verify_cat.py -- independent, from-scratch verification of the mixed-class
rows whose Apery limits hit the Catalan battery.  Nothing is imported from the
scan pipeline: coefficients are rebuilt from Theorem D3, integrality is exact,
and the limit is recomputed by high-precision forward recursion."""
from fractions import Fraction as Fr
from math import gcd
import mpmath as mp, sys
mp.mp.dps = 260

ROWS = [
 # (rho_p, rho_r, M, J1, J2, r, a, c, d, f, C)
 (Fr(-1,2), Fr(0), 1,0,0,  8, 16,  8,  48,  0, -128),
 (Fr(-1,2), Fr(0), 1,0,0,  2, 14,  8,  28,  4,    8),
 (Fr(-1,2), Fr(0), 1,0,0, -2,  6,  4, -32, -8,   32),
 (Fr(-1,2), Fr(0), 1,0,0,  8, 16,  8,  68,  8,   32),
 (Fr(-1,2), Fr(0), 1,0,0,  1, 17, 10,  32,  8,   16),
 (Fr(-1,2), Fr(0), 1,0,0, -1, 13,  8, -13, -1,   -1),
]
NEX = 200; NK = 80

def coeffs(rp, rr, M, J1, J2, r, a, c, d, f, C):
    g = C*M*M
    assert r != 0 and g % r == 0
    p = g//r; s = a - r
    assert d == s*r + p, (d, s*r+p)
    b = (1-rr)*r + (1-rp)*s
    e = -rp*(2*p + r*s) - rr*r*s
    h = -(1 + 2*rp + rr)*Fr(g)
    j = Fr(C*J1*J2)
    for x in (b,e,h,j): assert x.denominator == 1, (b,e,h,j)
    return a,int(b),c,d,int(e),f,g,int(h),int(j)

def rho_at(co, lam):
    """rho_i = -T(t_i)/(t_i Rc'(t_i)), T = Sc - t Rc' -- the Theorem F1 dictionary."""
    a,b,c,d,e,f,g,h,j = co
    t = 1/lam
    Rp = -a + 2*d*t - 3*g*t*t
    Tv = 1 - b*t + (d+e)*t*t - (2*g+h)*t**3
    return -Tv/(t*Rp)

def frob_has_log(co, lam):
    """at a point with coincident exponents (0,0) a logarithm is forced; here we
    only record which points have integral exponent difference at all."""
    return abs(rho_at(co, lam) - mp.nint(rho_at(co, lam).real)) < mp.mpf(10)**(-40)

for (rp,rr,M,J1,J2,r,a,c,d,f,C) in ROWS:
    co = coeffs(rp,rr,M,J1,J2,r,a,c,d,f,C)
    a_,b_,c_,d_,e_,f_,g_,h_,j_ = co
    P=lambda n:a_*n*n+b_*n+c_; Q=lambda n:d_*n*n+e_*n+f_; R=lambda n:g_*n*n+h_*n+j_
    print("="*78)
    print("row (a,c,d,f,C)=(%d,%d,%d,%d,%d)  rational root r=%d" % (a,c,d,f,C,r))
    print("  P(n) = %d n^2 + %d n + %d" % (a_,b_,c_))
    print("  Q(n) = %d n^2 + %d n + %d" % (d_,e_,f_))
    print("  R(n) = %d n^2 + %d n + %d" % (g_,h_,j_))
    # exact integrality
    v=[Fr(0),Fr(0),Fr(1)]
    for n in range(NEX):
        v.append((P(n)*v[n+2]-Q(n)*v[n+1]+R(n)*v[n])/Fr((n+1)**2))
    u=v[2:]
    print("  u_n =", [int(x) for x in u[:8]], "...")
    print("  integral to n=%d : %s   (u_%d has %d digits)" %
          (NEX, all(x.denominator==1 for x in u), NEX, len(str(abs(int(u[NEX]))))))
    # companion, sharp k
    w=[Fr(0),Fr(0),Fr(0),Fr(1)]
    for n in range(1,NK+1):
        while len(w) < n+3: w.append(Fr(0))
        w.append((P(n)*w[n+2]-Q(n)*w[n+1]+R(n)*w[n])/Fr((n+1)**2))
    bb=w[2:]
    L=1; ks=[]
    for n in range(1,NK+1):
        L=L*n//gcd(L,n); den=bb[n].denominator; kk=0
        while den>1:
            gg=gcd(den,L)
            if gg==1: kk=99; break
            den//=gg; kk+=1
        ks.append(kk)
    kd=max(ks)
    lam=mp.polyroots([1,-a_,d_,-g_],maxsteps=400,extraprec=400)
    lam=sorted(lam,key=lambda z:-abs(z))
    disc = 18*a_*d_*g_-4*a_**3*g_+a_*a_*d_*d_-4*d_**3-27*g_*g_
    print("  lambda =", [mp.nstr(z,18) for z in lam], "  disc =", disc)
    print("  rho_i  =", [mp.nstr(rho_at(co,z),12) for z in lam])
    print("  sharp k =", kd, "  score = %.4f" % float(-mp.log(abs(lam[1]))-kd))
    # high-precision limit
    rate = mp.log10(abs(lam[0])/abs(lam[1]))
    N = int(min(200000, max(600, 175/rate + 400)))
    with mp.workdps(mp.mp.dps+120):
        A=[mp.mpf(0),mp.mpf(0),mp.mpf(1)]; B=[mp.mpf(0),mp.mpf(0),mp.mpf(0),mp.mpf(1)]
        for n in range(N):
            pp,qq,rr2 = mp.mpf(P(n)),mp.mpf(Q(n)),mp.mpf(R(n))
            A.append((pp*A[n+2]-qq*A[n+1]+rr2*A[n])/mp.mpf((n+1)**2))
            if n>=1: B.append((pp*B[n+2]-qq*B[n+1]+rr2*B[n])/mp.mpf((n+1)**2))
        x1=B[N+2]/A[N+2]; x2=B[N+1]/A[N+1]
        conv = float(mp.log10(abs(x1-x2))) if x1!=x2 else None
        xi = +x1
    print("  n = %d   log10|xi_n - xi_{n-1}| = %s" % (N, conv))
    print("  xi =", mp.nstr(xi, 150))
    sys.stdout.flush()
    open('out/cat_rows_xi.txt','a').write("cat_%d_%d_%d_%d_%d_r%s %s\n" %
        (a,c,d,f,C,str(r).replace('-','m'), mp.nstr(xi, 140)))
