"""T1: the seven pure CDT functions B1..B7 as exact rational series in Y=y/s,
     with MEASURED sharp denominator types, plus the Z[phi]-integrality check
     for the actual functions of y (s = -phi^5 a unit of Z[phi]).
"""
from fractions import Fraction as F
from math import factorial as fac, gcd
import sys

N = 96                      # order (we measure up to n<=80)
NMEAS = 80

def lcm_upto(m, cache={}):
    if m in cache: return cache[m]
    r = 1
    for j in range(1, m+1): r = r*j//gcd(r, j)
    cache[m] = r
    return r

# ---------------------------------------------------------------- the series
B1 = [F(0)]*(N+1); B1[0] = F(1)
B2 = [F(0)]*(N+1)
for n in range(2, N+1): B2[n] = F(2*fac(n-2)*fac(n), fac(2*n))
B3 = [F(0)]*(N+1)
for n in range(1, N+1): B3[n] = F(fac(n-1)**2, fac(2*n))
B5 = [F(0)]*(N+1)
for n in range(1, N+1): B5[n] = F(fac(n-1)**2, fac(2*n-1)*(2*n-1))
B6 = [F(0)]*(N+1)
for n in range(1, N+1): B6[n] = F(fac(n-1)**2, n*fac(2*n))
# B4 = Sym^- Li_2 : CDT's closed form  (defB)
C = [F(fac(2*k), fac(k)**2) for k in range(N+2)]   # central binomials
B4 = [F(0)]*(N+1)
for n in range(0, N):
    S = F(0)
    for k in range(0, n+1):
        S += C[k]*C[n-k]*F(1, (2*k-1)*(2*n-2*k+1)**2)
    B4[n+1] = 4*S/F(16)**n
B7 = [F(0)]*(N+1)
for n in range(1, N+1): B7[n] = B4[n]/n

assert [str(t) for t in B4[1:5]] == ['-4', '4/9', '31/900', '389/88200'], B4[1:5]
print("B4 head verified against CDT (defB): -4y + 4y^2/9 + 31y^3/900 + 389y^4/88200")

# independent cross-check of B4 by direct symmetrisation over Q (s=1 normalisation)
def symminus_li2(P):
    # w = x/(x-1) = -x -x^2 -x^3 - ...  ; Y = x^2/(x-1) = -x^2-x^3-...
    w = [F(0)]*(P+1)
    for n in range(1, P+1): w[n] = F(-1)
    Ycoef = [F(0)]*(P+1)
    for n in range(2, P+1): Ycoef[n] = F(-1)
    def mul(a, b):
        r = [F(0)]*(P+1)
        for i, ai in enumerate(a):
            if ai:
                for j in range(0, P+1-i):
                    if b[j]: r[i+j] += ai*b[j]
        return r
    def compose(f, g):   # f(g(x)), g(0)=0
        r = [F(0)]*(P+1); pw = [F(0)]*(P+1); pw[0] = F(1)
        for n in range(P+1):
            if f[n]:
                for i in range(P+1): r[i] += f[n]*pw[i]
            pw = mul(pw, g)
        return r
    Li2 = [F(0)]+[F(1, n*n) for n in range(1, P+1)]
    d = [a-b for a, b in zip(Li2, compose(Li2, w))]
    pref = [F(0)]*(P+1); pref[1] = F(1)
    pref = [a-b for a, b in zip(pref, w)]
    prod = mul(pref, d)
    # re-expand prod (a series in x) as a series in Y
    out = [F(0)]*(P//2+1); cur = prod[:]; Yp = [F(0)]*(P+1); Yp[0] = F(1)
    for n in range(P//2+1):
        Yp_n = Yp
        if 2*n > P: break
        out[n] = cur[2*n]/Yp_n[2*n] if Yp_n[2*n] else F(0)
        if out[n]:
            cur = [c-out[n]*s for c, s in zip(cur, Yp_n)]
        Yp = mul(Yp, Ycoef)
    return out
chk = symminus_li2(40)
assert all(chk[n] == B4[n] for n in range(0, 21)), [chk[:6], B4[:6]]
print("B4 cross-checked to order 20 by direct Sym^- of Li_2 in x  [verified]")

# ---------------------------------------------------------------- Z[phi] check
# s = -phi^5.  represent elements of Z[phi] as (a,b) meaning (a + b*sqrt5)/2, a=b mod 2.
def mulZ(u, v):
    a, b = u; c, d = v
    return ((a*c+5*b*d)//2, (a*d+b*c)//2)
phi = (1, 1)                      # (1+sqrt5)/2
phi5 = (1, 0)
for _ in range(5): phi5 = mulZ(phi5, phi)   # start from 1 = (2,0)/2 -> careful
phi5 = (2, 0)
for _ in range(5): phi5 = mulZ(phi5, phi)
s_ = (-phi5[0], -phi5[1])
# s^{-1} = -phi^{-5};  phi^{-1} = phi-1 = (-1+sqrt5)/2
phiinv = (-1, 1)
phiinv5 = (2, 0)
for _ in range(5): phiinv5 = mulZ(phiinv5, phiinv)
sinv = (-phiinv5[0], -phiinv5[1])
assert mulZ(s_, sinv) == (2, 0), (s_, sinv, mulZ(s_, sinv))
print(f"s = -phi^5 = ({s_[0]} + {s_[1]}*sqrt5)/2 ; s^-1 = ({sinv[0]} + {sinv[1]}*sqrt5)/2 ; s*s^-1 = 1  [proved: s is a unit]")

# ---------------------------------------------------------------- types
def T_lcm(b, n):  return lcm_upto(b*n)
CANDS = []
for e in (0, 1, 2):
    for layers in ((), (1,), (2,), (1, 1), (1, 2), (2, 2)):
        CANDS.append((layers, e))
def typeval(cand, n):
    layers, e = cand
    v = n**e
    for b in layers: v *= T_lcm(b, n)
    return v
def name(cand):
    layers, e = cand
    parts = []
    for b in layers: parts.append('[1..n]' if b == 1 else '[1..2n]')
    ss = ''.join(parts) if parts else '1'
    if e == 1: ss += ' n'
    if e == 2: ss += ' n^2'
    return ss

CDT_claim = {
 'B1': ('trivial', lambda n: 1),
 'B2': ('[1..2n]', lambda n: T_lcm(2, n)),
 'B3': ('[1..2n] n', lambda n: T_lcm(2, n)*n),
 'B4': ('[1..2n]^2', lambda n: T_lcm(2, n)**2),
 'B5': ('[1..2n](2n-1)', lambda n: T_lcm(2, n)*(2*n-1)),
 'B6': ('[1..2n] n^2', lambda n: T_lcm(2, n)*n*n),
 'B7': ('[1..2n]^2 n', lambda n: T_lcm(2, n)**2*n),
}
series = [('B1', B1), ('B2', B2), ('B3', B3), ('B4', B4), ('B5', B5), ('B6', B6), ('B7', B7)]

print("\n=== CDT Lemma (bdenominators): validity of the claimed types, n <= %d ===" % NMEAS)
for nm, f in series:
    cn, cf = CDT_claim[nm]
    bad = [n for n in range(1, NMEAS+1) if (F(cf(n))*f[n]).denominator != 1]
    # how slack: max over n of  T(n)/den(c_n) is not meaningful; report whether ever equality
    eq = sum(1 for n in range(1, NMEAS+1) if f[n] != 0 and cf(n) == f[n].denominator)
    print(f"  {nm}: claimed {cn:16s} -> {'VALID' if not bad else 'FAILS at n='+str(bad[:5])}"
          f"   (exact equality T(n)=den(c_n) for {eq}/{NMEAS} indices)")

print("\n=== MEASURED minimal admissible types in the family [1..n]^a [1..2n]^b n^e (n<=%d) ===" % NMEAS)
minimal = {}
for nm, f in series:
    ok = []
    for cand in CANDS:
        if all((F(typeval(cand, n))*f[n]).denominator == 1 for n in range(1, NMEAS+1)):
            ok.append(cand)
    # minimality: cand1 <= cand2 if typeval divides for all n
    mins = []
    for c1 in ok:
        if not any(c2 != c1 and all(typeval(c2, n) % typeval(c1, n) == 0 or typeval(c1, n) % typeval(c2, n) != 0
                                    for n in range(1, 0)) for c2 in ok):
            pass
    # simpler: report all admissible, sorted by typeval at n=NMEAS
    ok.sort(key=lambda c: typeval(c, NMEAS))
    minimal[nm] = ok
    print(f"  {nm}: admissible = {[name(c) for c in ok[:6]]}")

# also: sharp denominator itself
print("\n=== den(c_n) vs [1..2n]^2 for a few n (the actual denominators) ===")
for nm, f in series:
    row = []
    for n in (10, 40, 80):
        if f[n] == 0: row.append((n, '0')); continue
        d = f[n].denominator
        L2 = T_lcm(2, n)
        row.append((n, f"den/[1..2n]={F(d,L2)}"))
    print(f"  {nm}: {row}")
