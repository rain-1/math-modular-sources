"""Exact evaluation of  T(p,q,t) = \\iint_Delta x^p y^q (1-x^2-y^2)^{-(t+1)} dx dy
for even p,q and any integer t <= (p+q)/2.  Extends emnexact.py to p != q.

Reduction (see EMN_VERIFICATION.md sec.3):
  polar + u = r^2 + the substitution p_ang = sin(2 theta):
     T(p,q,t) = int_0^{pi/2} cos^p sin^q (theta) * Inner_{mm}(theta) d theta ,  mm = (p+q)/2,
     Inner_mm(theta) = (1/2) int_0^{W} u^{mm}(1-u)^{-t-1} du ,  W = 1/(1+sin 2theta).
  For an x<->y symmetric numerator the angular weight cos^p sin^q + cos^q sin^p is a
  polynomial in (cos sin) = sin(2theta)/2, and int_0^{pi/2} f(sin 2theta) dtheta
  = int_0^{pi/2} f(sin psi) dpsi.  So everything reduces to atoms with weight sin^k(psi).

Module basis: (1, pi, log2, pi*log2, G) -- 5-tuples of Fractions.
"""
from fractions import Fraction as Fr
from math import comb
from functools import lru_cache

NB = 5
ZERO = (Fr(0),)*NB
def vadd(a,b): return tuple(x+y for x,y in zip(a,b))
def vsub(a,b): return tuple(x-y for x,y in zip(a,b))
def vmul(c,a):
    c = Fr(c); return tuple(c*x for x in a)

@lru_cache(maxsize=None)
def Cw(k):
    """int_0^{pi/2} sin^k."""
    if k == 0: return (Fr(0), Fr(1,2), Fr(0), Fr(0), Fr(0))
    if k == 1: return (Fr(1), Fr(0), Fr(0), Fr(0), Fr(0))
    return vmul(Fr(k-1,k), Cw(k-2))

@lru_cache(maxsize=None)
def Pw(r):
    """int_0^{pi/2}(1+sin)^r, r in Z."""
    if r >= 0:
        out = ZERO
        for l in range(r+1): out = vadd(out, vmul(comb(r,l), Cw(l)))
        return out
    q = -r
    poly = {0: Fr(1)}; base = {2: Fr(1), 1: Fr(-2), 0: Fr(2)}
    for _ in range(q-1):
        new = {}
        for d1,c1 in poly.items():
            for d2,c2 in base.items(): new[d1+d2] = new.get(d1+d2, Fr(0)) + c1*c2
        poly = new
    tot = Fr(0)
    for d,c in poly.items():
        e = d - 2*q
        assert e != -1
        tot += c*(Fr(2)**(e+1) - 1)/(e+1)
    return (2*tot, Fr(0), Fr(0), Fr(0), Fr(0))

@lru_cache(maxsize=None)
def Ew(k,i):
    out = ZERO
    for j in range(k+1): out = vadd(out, vmul(comb(k,j)*(-1)**(k-j), Pw(j-i)))
    return out

@lru_cache(maxsize=None)
def Lw(m):
    """int_0^{pi/2} sin^m log sin = Cw(m)*(digamma((m+1)/2)-digamma(m/2+1))/2."""
    c = Cw(m)
    if m % 2 == 0:
        k = m//2
        rp = sum(Fr(1,2*j-1) for j in range(1,k+1)) - sum(Fr(1,j) for j in range(1,k+1))/2
        r = c[1]; assert c[0]==0
        return (Fr(0), r*rp, Fr(0), -r, Fr(0))
    else:
        k = (m-1)//2
        rp = sum(Fr(1,j) for j in range(1,k+1))/2 - sum(Fr(1,2*j-1) for j in range(1,k+2))
        r = c[0]; assert c[1]==0
        return (r*rp, Fr(0), r, Fr(0), Fr(0))

@lru_cache(maxsize=None)
def Sw(m):
    """int_0^{pi/2} sin^m log(1+sin);  m*Sw(m) = (m-1)Sw(m-2)+Cw(m-1)-Cw(m)."""
    if m == 0: return (Fr(0), Fr(0), Fr(0), Fr(-1,2), Fr(2))    # 2G - (pi/2)log2
    if m == 1: return vsub(Cw(0), Cw(1))                        # pi/2 - 1
    return vmul(Fr(1,m), vadd(vmul(m-1, Sw(m-2)), vsub(Cw(m-1), Cw(m))))

@lru_cache(maxsize=None)
def genT(k, mm, t):
    """int_0^{pi/2} (cos sin)^k(theta) * Inner_mm(theta) dtheta
       = 2^{-(k+1)} int_0^{pi/2} sin^k(psi) [int_0^W u^mm (1-u)^{-t-1} du] dpsi."""
    out = ZERO
    for i in range(mm+1):
        if i == t: continue
        out = vadd(out, vmul(Fr(comb(mm,i)*(-1)**i, i-t), vsub(Cw(k), Ew(k+i-t, i-t))))
    if 0 <= t <= mm:
        out = vadd(out, vmul(comb(mm,t)*(-1)**t, vsub(Sw(k), Lw(k))))
    return vmul(Fr(1, 2**(k+1)), out)

@lru_cache(maxsize=None)
def symweight(p, q):
    """cos^p sin^q + cos^q sin^p (p<q, q-p even) as {k: coef} in (cos sin)^k;
       for p==q returns {p:1}."""
    if p == q: return ((p, Fr(1)),)
    if p > q: p, q = q, p
    e = (q-p)//2
    # P_r = s^{2r}+c^{2r} : P_0=2, P_1=1, P_r = P_{r-1} - X*P_{r-2}, X=(cs)^2
    P0 = {0: Fr(2)}; P1 = {0: Fr(1)}
    if e == 0: Pe = P0
    elif e == 1: Pe = P1
    else:
        for _ in range(2, e+1):
            new = dict(P1)
            for d,c in P0.items(): new[d+1] = new.get(d+1, Fr(0)) - c
            P0, P1 = P1, new
        Pe = P1
    return tuple(sorted((p + 2*d, c) for d,c in Pe.items()))

def Tmono(p, q, t):
    """T(p,q,t)+T(q,p,t) if p!=q, else T(p,p,t).  Returns the 5-tuple."""
    mm = (p+q)//2
    out = ZERO
    for k, c in symweight(p,q):
        out = vadd(out, vmul(c, genT(k, mm, t)))
    return out

def linform_gen(v):
    assert v[1]==0 and v[2]==0 and v[3]==0, f"pi/log2 parts did not cancel: {v}"
    return v[0], v[4]
