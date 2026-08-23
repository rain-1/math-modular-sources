"""Exact evaluation of the EMN higher-pole periods
      I_{m,t} = \\iint_Delta x^m y^m (1-x^2-y^2)^{-(t+1)} dx dy ,   m even, t in Z.

Derivation (all steps elementary; see EMN_VERIFICATION.md sec.3):
  polar/half-angle reduction with p = sin(2*theta) gives
      I_{m,t} = 2^{-(m+1)} \\int_0^{pi/2} sin^m(psi) * [ \\int_0^{W} u^m (1-u)^{-t-1} du ] dpsi,
      W = 1/(1+sin psi),  1-W = sin psi/(1+sin psi).
  Expanding u^m = ((1-(1-u))^m and integrating in u:
      inner = sum_{i != t} C(m,i)(-1)^i (1 - v^{i-t})/(i-t)  +  [i=t] C(m,t)(-1)^t (-log v),
      v = 1-W = s/(1+s).
Atoms (exact, in the module  Q + Q*pi + Q*pi*log2 + Q*G):
  Cw(k)   = int_0^{pi/2} sin^k                     (Wallis)
  Pw(r)   = int_0^{pi/2} (1+sin)^r,  r in Z        (rational for r<0, Q+Q pi for r>=0)
  Ew(k,i) = int_0^{pi/2} sin^k/(1+sin)^i           = sum_j C(k,j)(-1)^{k-j} Pw(j-i)
  Lw(m)   = int_0^{pi/2} sin^m log sin             (= Cw(m)*(psi-difference)/2)
  Sw(m)   = int_0^{pi/2} sin^m log(1+sin),  m*Sw(m) = (m-1)Sw(m-2) + Cw(m-1) - Cw(m),
            Sw(0) = 2G - (pi/2) log 2

Module element = 4-tuple of Fractions (q, qpi, qpilog2, qG).
Names deliberately avoid gp/py builtins (no psi, M, Phi, S, cmp).
"""
from fractions import Fraction as Fr
from math import comb
from functools import lru_cache

ZERO = (Fr(0),)*4
def vadd(a,b): return tuple(x+y for x,y in zip(a,b))
def vsub(a,b): return tuple(x-y for x,y in zip(a,b))
def vmul(c,a):
    c = Fr(c); return tuple(c*x for x in a)

@lru_cache(maxsize=None)
def Cw(k):
    """int_0^{pi/2} sin^k psi dpsi."""
    if k == 0: return (Fr(0), Fr(1,2), Fr(0), Fr(0))     # pi/2
    if k == 1: return (Fr(1), Fr(0), Fr(0), Fr(0))
    return vmul(Fr(k-1,k), Cw(k-2))

@lru_cache(maxsize=None)
def Pw(r):
    """int_0^{pi/2} (1+sin psi)^r dpsi, r in Z."""
    if r >= 0:
        out = ZERO
        for l in range(r+1):
            out = vadd(out, vmul(comb(r,l), Cw(l)))
        return out
    q = -r
    # 2 * int_0^1 (1+t^2)^{q-1}/(1+t)^{2q} dt, t = tan(psi/2); set u = 1+t
    # (1+t^2) = u^2-2u+2
    poly = {0: Fr(1)}                                    # (u^2-2u+2)^{q-1} as dict deg->coef
    base = {2: Fr(1), 1: Fr(-2), 0: Fr(2)}
    for _ in range(q-1):
        new = {}
        for d1,c1 in poly.items():
            for d2,c2 in base.items():
                new[d1+d2] = new.get(d1+d2, Fr(0)) + c1*c2
        poly = new
    tot = Fr(0)
    for d,c in poly.items():
        e = d - 2*q                                      # integrate u^e from 1 to 2
        assert e != -1, "unexpected log term"
        tot += c * (Fr(2)**(e+1) - 1)/(e+1)
    return (2*tot, Fr(0), Fr(0), Fr(0))

@lru_cache(maxsize=None)
def Ew(k, i):
    """int_0^{pi/2} sin^k/(1+sin)^i dpsi, k>=0, i in Z."""
    out = ZERO
    for j in range(k+1):
        out = vadd(out, vmul(comb(k,j)*(-1)**(k-j), Pw(j-i)))
    return out

@lru_cache(maxsize=None)
def Lw(m):
    """int_0^{pi/2} sin^m log(sin) dpsi = Cw(m)*(digamma((m+1)/2)-digamma(m/2+1))/2."""
    assert m % 2 == 0
    k = m//2
    ratpart = sum(Fr(1, 2*j-1) for j in range(1, k+1)) - sum(Fr(1,j) for j in range(1,k+1))/2
    c = Cw(m)                                            # = (0, r, 0, 0)
    assert c[0]==0 and c[2]==0 and c[3]==0
    r = c[1]
    return (Fr(0), r*ratpart, -r, Fr(0))                 # r*pi*ratpart - r*pi*log2

@lru_cache(maxsize=None)
def Sw(m):
    """int_0^{pi/2} sin^m log(1+sin) dpsi."""
    if m == 0: return (Fr(0), Fr(0), Fr(-1,2), Fr(2))    # 2G - (pi/2) log 2
    if m == 1:
        # 1*S_1 = 0*S_{-1} + C_0 - C_1 ... recurrence invalid at m=1; do directly:
        # int_0^{pi/2} sin*log(1+sin) = 1 - pi/2 + log2 ... not needed (m even only)
        raise NotImplementedError
    return vmul(Fr(1,m), vadd(vmul(m-1, Sw(m-2)), vsub(Cw(m-1), Cw(m))))

def Imt(m, t):
    """I_{m,t}; returns the 4-tuple.  m even >= 0, t any integer <= m."""
    out = ZERO
    for i in range(m+1):
        if i == t: continue
        c = Fr(comb(m,i)*(-1)**i, i-t)
        out = vadd(out, vmul(c, vsub(Cw(m), Ew(m+i-t, i-t))))
    if 0 <= t <= m:
        out = vadd(out, vmul(comb(m,t)*(-1)**t, vsub(Sw(m), Lw(m))))
    return vmul(Fr(1, 2**(m+1)), out)

def linform(m, t):
    """Return (a, b) with I_{m,t} = a + b*G, asserting the pi and pi*log2 parts vanish."""
    v = Imt(m, t)
    assert v[1] == 0 and v[2] == 0, f"pi-part did not cancel at (m,t)=({m},{t}): {v}"
    return v[0], v[3]
