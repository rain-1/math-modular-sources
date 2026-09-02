"""04: exact (over Q) two-variable recurrences for the decoupled rows, via
mod-p nullspace at two primes + CRT + rational reconstruction."""
from lib2v import *
from math import comb
from fractions import Fraction as F
import sys

PRIMES = [1000003, 999983, 999979]


def nullspace_modp(rows, ncols, p):
    M = [r[:] for r in rows]
    nrows = len(M)
    piv, r = [], 0
    for c in range(ncols):
        pr = None
        for i in range(r, nrows):
            if M[i][c] % p:
                pr = i
                break
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        inv = pow(M[r][c], p-2, p)
        M[r] = [(v*inv) % p for v in M[r]]
        for i in range(nrows):
            if i != r and M[i][c] % p:
                f = M[i][c]
                M[i] = [(M[i][j] - f*M[r][j]) % p for j in range(ncols)]
        piv.append(c)
        r += 1
        if r == nrows:
            break
    free = [c for c in range(ncols) if c not in piv]
    basis = []
    for fc in free:
        v = [0]*ncols
        v[fc] = 1
        for i, pc in enumerate(piv):
            v[pc] = (-M[i][fc]) % p
        basis.append(v)
    return basis, piv, free


def ratrec(a, m):
    """rational reconstruction of a mod m"""
    import math
    bound = int(math.isqrt(m // 2))
    r0, r1 = m, a % m
    s0, s1 = 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q*r1
        s0, s1 = s1, s0 - q*s1
    if s1 == 0 or math.gcd(r1, s1) != 1:
        return None
    return F(r1 if s1 > 0 else -r1, abs(s1))


def find_exact(cfun, I, d, axis, A, B, shift=1, normalise_index=None):
    mons = [(s, t) for s in range(d+1) for t in range(d+1) if s+t <= d]
    ncols = (I+1)*len(mons)
    tabs = {}
    vecs = []
    for p in PRIMES:
        tab = [[cfun(a, b) % p for b in range(B+1)] for a in range(A+1)]
        rows = []
        amax = A - (I if axis == 0 else 0)
        bmax = B - (I if axis == 1 else 0)
        for a in range(shift, amax+1):
            for b in range(shift, bmax+1):
                row = []
                for i in range(I+1):
                    val = tab[a+i][b] if axis == 0 else tab[a][b+i]
                    for (s, t) in mons:
                        row.append((val*pow(a, s, p)*pow(b, t, p)) % p)
                rows.append(row)
        ns, piv, free = nullspace_modp(rows, ncols, p)
        if len(ns) != 1:
            return None, mons, len(ns)
        vecs.append((p, ns[0], free[0]))
    # CRT
    M = 1
    for p, _, _ in vecs:
        M *= p
    out = []
    for j in range(ncols):
        x = 0
        for p, v, _ in vecs:
            Mp = M//p
            x = (x + v[j]*Mp*pow(Mp, -1, p)) % M
        out.append(ratrec(x, M))
    return out, mons, 1


def show(vec, mons, I, label, var="a"):
    nm = len(mons)
    print("   " + label)
    for idx in range(I+1):
        seg = vec[idx*nm:(idx+1)*nm]
        terms = []
        for (c, (s, t)) in zip(seg, mons):
            if c is None:
                terms.append("??")
                continue
            if c == 0:
                continue
            terms.append("%+s a^%d b^%d" % (str(c), s, t))
        print("     P_%d = %s" % (idx, " ".join(terms) if terms else "0"))


print("=" * 78)
print("zeta3 D1: a-recurrence, order 2, degree 3")
v, mons, k = find_exact(z3_D1, 2, 3, 0, 30, 30)
show(v, mons, 2, "sum_i P_i(a,b) c(a+i,b) = 0")

print()
print("=" * 78)
print("zeta2 D1: a-recurrence, order 2, degree 2")
v, mons, k = find_exact(z2_D1, 2, 2, 0, 30, 30)
show(v, mons, 2, "sum_i P_i(a,b) c(a+i,b) = 0")
print("zeta2 D1: b-recurrence, order 2, degree 2")
v, mons, k = find_exact(z2_D1, 2, 2, 1, 30, 30)
show(v, mons, 2, "sum_j Q_j(a,b) c(a,b+j) = 0")

print()
print("=" * 78)
print("rowE D1: a-recurrence, order 3, degree 2")
v, mons, k = find_exact(E_D1, 3, 2, 0, 30, 30)
show(v, mons, 3, "sum_i P_i(a,b) c(a+i,b) = 0")
