#!/usr/bin/env python3
"""Verify the mixed-exponent normalisation forms of 03_fmix.c against the
Theorem F1 dictionary rho_i = -T(t_i)/(t_i Rc'(t_i)),  T = Sc - t Rc'."""
from fractions import Fraction as Fr
import mpmath as mp, itertools, random
mp.mp.dps = 40

def build(rp, rr, M, J1, J2, r, s, p, c, f):
    a = r+s; d = s*r+p; g = r*p
    assert Fr(g, M*M).denominator == 1
    C = Fr(g, M*M)
    b = (1-rr)*r + (1-rp)*s
    e = -rp*(2*p + r*s) - rr*r*s
    h = -(1 + 2*rp + rr)*g
    j = C*J1*J2
    return dict(a=a,b=b,c=c,d=d,e=e,f=f,g=g,h=h,j=j,C=C)

def rhos(co):
    a,b,c,d,e,f,g,h,j = (co[k] for k in 'abcdefghj')
    Rc = [1,-a,d,-g]           # low->high
    T  = [1,-b,d+e,-(2*g+h)]
    lam = mp.polyroots([1,-mp.mpf(a),mp.mpf(d),-mp.mpf(g)], maxsteps=300, extraprec=300)
    out=[]
    for L in lam:
        t = 1/L
        Rp = -a + 2*d*t - 3*g*t*t
        Tv = 1 - b*t + (d+e)*t*t - (2*g+h)*t**3
        out.append(-Tv/(t*Rp))
    return lam, out

random.seed(1)
bad=0
CLS = [
 (Fr(1,2), Fr(0), 4, 3, 5), (Fr(1,2), Fr(0), 6, 5, 7), (Fr(1,2), Fr(0), 3, 2, 4),
 (Fr(-1,2), Fr(0), 1, 0, 0), (Fr(-1,2), Fr(0), 4,-1, 1), (Fr(-1,2), Fr(0), 3,-1, 1),
 (Fr(0), Fr(1,2), 4, 3, 3), (Fr(0), Fr(1,2), 12, 7, 11), (Fr(0), Fr(1,2), 2, 0, 3),
 (Fr(0), Fr(-1,2), 4, 1, 1), (Fr(0), Fr(-1,2), 2, 0, 1), (Fr(0), Fr(-1,2), 12, 1, 5),
 (Fr(1,2), Fr(-1,2), 4, 3, 3),
 (Fr(0), Fr(1), 4, 3, 5), (Fr(0), Fr(-1), 4,-1, 1),
]
for (rp, rr, M, J1, J2) in CLS:
    assert Fr(J1+J2, M) == 1 + 2*rp + rr, (rp,rr,M,J1,J2, Fr(J1+J2,M))
    for _ in range(6):
        r = random.choice([x for x in range(-6,7) if x])
        s = random.randint(-9,9); p = random.choice([x for x in range(-9,10) if x])
        if Fr(r*p, M*M).denominator != 1: continue
        co = build(rp,rr,M,J1,J2,r,s,p,random.randint(-5,5),random.randint(-5,5))
        if all(x.denominator==1 for x in [Fr(co['b']),Fr(co['e']),Fr(co['h']),Fr(co['j'])]) is False: continue
        disc = 18*co['a']*co['d']*co['g']-4*co['a']**3*co['g']+co['a']**2*co['d']**2-4*co['d']**3-27*co['g']**2
        if disc == 0: continue
        lam, rr_ = rhos(co)
        # the root closest to r should carry rho_r; the other two rho_p
        idx = min(range(3), key=lambda i: abs(lam[i]-r))
        ok = abs(rr_[idx]-mp.mpf(float(rr))) < 1e-20
        for i in range(3):
            if i!=idx: ok = ok and abs(rr_[i]-mp.mpf(float(rp))) < 1e-20
        # exponents at infinity from R(n)=g n^2 + h n + j
        gg,hh,jj = mp.mpf(co['g']), mp.mpf(float(co['h'])), mp.mpf(float(co['j']))
        d_inf = mp.sqrt(abs(hh*hh-4*gg*jj))/abs(gg)
        ok2 = abs(d_inf - abs(mp.mpf(J2-J1)/M)) < 1e-20
        if not (ok and ok2):
            bad += 1
            print("BAD", rp,rr,M,J1,J2, r,s,p, [mp.nstr(x,12) for x in rr_], mp.nstr(d_inf,12))
print("checked classes:", len(CLS), " failures:", bad)
