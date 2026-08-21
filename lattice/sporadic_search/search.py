#!/usr/bin/env python3
"""
Systematic search for sporadic Apery-like sequences from eta-quotient data
on genus-zero Gamma_0(N).

Construction (see consolidation/ZETA3_TWO_LATTICE.md 7.4 item 3):
  t = prod_{d|N} eta(d tau)^{r_d},  sum r_d = 0, sum d r_d = 24   (t = q + O(q^2))
  F = prod_{d|N} eta(d tau)^{s_d},  sum s_d = 2w, sum d s_d = 0    (F = 1 + O(q))
  A(t) = F  =>  a_n
Test integrality and a 3-term recurrence with poly coefficients of degree <= 3.
"""
import sys, itertools
from fractions import Fraction

PREC = 46          # scan precision (number of q-coefficients)
PREC_HI = 66       # verification precision
GENUS0 = [2,3,4,5,6,7,8,9,10,12,13,16,18,25]

# ---------------- power series utilities (integer coefficient lists) -------
def divisors(n):
    return [d for d in range(1, n+1) if n % d == 0]

def eta_factor(d, P):
    """prod_{n>=1} (1 - q^{d n}) truncated to P terms (index 0..P-1)."""
    f = [0]*P; f[0] = 1
    n = 1
    while d*n < P:
        # multiply by (1 - q^{dn})
        m = d*n
        for i in range(P-1, m-1, -1):
            f[i] -= f[i-m]
        n += 1
    return f

def smul(a, b, P):
    r = [0]*P
    for i, ai in enumerate(a):
        if ai == 0: continue
        lim = P - i
        for j in range(min(len(b), lim)):
            bj = b[j]
            if bj: r[i+j] += ai*bj
    return r

def sinv(a, P):
    """inverse of series with a[0] = 1"""
    assert a[0] == 1
    r = [0]*P; r[0] = 1
    for n in range(1, P):
        s = 0
        for k in range(1, min(n, len(a)-1)+1):
            if a[k]: s += a[k]*r[n-k]
        r[n] = -s
    return r

def spow(a, e, P):
    if e == 0:
        r = [0]*P; r[0] = 1; return r
    neg = e < 0
    e = abs(e)
    base = a[:]; res = [0]*P; res[0] = 1
    while e:
        if e & 1: res = smul(res, base, P)
        e >>= 1
        if e: base = smul(base, base, P)
    return sinv(res, P) if neg else res

# ---------------- Ligozat cusp orders -------------------------------------
def cusp_orders(N, D, r):
    """order of prod eta(d tau)^{r_d} at cusp c (c | N), in local uniformiser.
       ord_c = N/(24 gcd(c^2,N)) * sum_d gcd(c,d)^2 r_d / d
       returned as Fractions, keyed by c."""
    from math import gcd
    out = {}
    for c in D:
        s = Fraction(0)
        for d, rd in zip(D, r):
            if rd: s += Fraction(gcd(c, d)**2 * rd, d)
        out[c] = Fraction(N, 24*gcd(c*c, N)) * s
    return out

def is_rational_square(D, r):
    """prod d^{r_d} a rational square?"""
    from sympy import factorint
    exps = {}
    for d, rd in zip(D, r):
        if d == 1 or rd == 0: continue
        for p, e in factorint(d).items():
            exps[p] = exps.get(p, 0) + e*rd
    return all(e % 2 == 0 for e in exps.values())

def ligozat_ok(N, D, r):
    if sum(d*rd for d, rd in zip(D, r)) % 24 != 0: return False
    if sum((N//d)*rd for d, rd in zip(D, r)) % 24 != 0: return False
    return True  # nebentypus allowed: only the mod-24 Ligozat conditions imposed

# ---------------- enumeration of exponent vectors -------------------------
def enum_vectors(D, S0, S1, B):
    """all integer r with sum r = S0, sum d r = S1, |r_d| <= B."""
    k = len(D)
    if k == 1:
        r0 = S0
        if D[0]*r0 == S1 and abs(r0) <= B: yield (r0,)
        return
    dA, dB = D[-2], D[-1]
    det = dB - dA
    rng = range(-B, B+1)
    for free in itertools.product(rng, repeat=k-2):
        s = sum(free); t = sum(d*x for d, x in zip(D[:-2], free))
        u = S0 - s; v = S1 - t
        # rA + rB = u ; dA rA + dB rB = v
        num = v - dA*u
        if num % det: continue
        rB = num//det; rA = u - rB
        if abs(rA) > B or abs(rB) > B: continue
        yield tuple(free) + (rA, rB)

# ---------------- recurrence fitting --------------------------------------
def fit_recurrence(a, deg=3, nfit=None):
    """find p2(n)a[n+1] + p1(n)a[n] + p0(n)a[n-1] = 0, deg <= deg.
       returns list of coefficient vectors spanning the solution space."""
    from sympy import Matrix, Rational
    L = len(a)
    if nfit is None: nfit = min(L-2, 22)
    rows = []
    for n in range(1, nfit+1):
        row = []
        for shift, val in ((1, a[n+1]), (0, a[n]), (-1, a[n-1])):
            for j in range(deg+1):
                row.append(val * n**j)
        rows.append(row)
    M = Matrix(rows)
    return M.nullspace()

def poly_eval(c, n):
    return sum(ci * n**i for i, ci in enumerate(c))

def check_recurrence(vec, a, deg=3):
    p2 = vec[0:deg+1]; p1 = vec[deg+1:2*deg+2]; p0 = vec[2*deg+2:3*deg+3]
    for n in range(1, len(a)-1):
        if poly_eval(p2, n)*a[n+1] + poly_eval(p1, n)*a[n] + poly_eval(p0, n)*a[n-1] != 0:
            return False, None
    return True, (p2, p1, p0)

# ---------------- classification into Zagier normalisation ----------------
def ptrim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0: p.pop()
    return p

def pdivexact(p, q):
    """exact division of poly p by poly q (coeff lists, low->high). None if not exact."""
    p = ptrim(p); q = ptrim(q)
    if len(p) < len(q): return None if any(p) else [Fraction(0)]
    out = [Fraction(0)]*(len(p)-len(q)+1)
    p = [Fraction(x) for x in p]
    for i in range(len(p)-len(q), -1, -1):
        cf = Fraction(p[i+len(q)-1], q[-1])
        out[i] = cf
        for j in range(len(q)):
            p[i+j] -= cf*q[j]
    if any(x != 0 for x in p): return None
    return out

def pmul(a, b):
    r = [Fraction(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): r[i+j] += x*y
    return r

def classify(p2, p1, p0):
    """Return (order, a, b, c) in Zagier / Almkvist-Zudilin normalisation, else None."""
    P2 = [Fraction(x) for x in p2]; P1 = [Fraction(x) for x in p1]; P0 = [Fraction(x) for x in p0]
    for k in (2, 3):
        tgt = [Fraction(1)]
        for _ in range(k): tgt = pmul(tgt, [Fraction(1), Fraction(1)])
        q = pdivexact(P2, tgt)
        if q is None or len(ptrim(q)) != 1 or q[0] == 0: continue
        s = q[0]
        A1 = ptrim([-x/s for x in P1]); A0 = ptrim([x/s for x in P0])
        if k == 2:
            ca = A1 + [Fraction(0)]*4; cb = A0 + [Fraction(0)]*4
            if len(ptrim(A1)) > 3 or len(ptrim(A0)) > 3: continue
            a, b = ca[2], ca[0]
            if ca[1] != a: continue
            if cb[0] or cb[1]: continue
            return (2, a, b, cb[2])
        else:
            quo = pdivexact(A1, [Fraction(1), Fraction(2)])
            if quo is None: continue
            cq = ptrim(quo) + [Fraction(0)]*4
            if len(ptrim(quo)) > 3: continue
            a, b = cq[2], cq[0]
            if cq[1] != a: continue
            cb = ptrim(A0) + [Fraction(0)]*5
            if len(ptrim(A0)) > 4: continue
            if cb[0] or cb[1] or cb[2]: continue
            return (3, a, b, cb[3])
    return None

def roots_of(order, a, c):
    import math
    A = float(a) if order == 2 else 2*float(a)
    C = float(c)
    disc = A*A - 4*C
    if disc >= 0:
        s = math.sqrt(disc)
        return ((A+s)/2, (A-s)/2)
    return None
