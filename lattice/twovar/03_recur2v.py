"""03: two-variable holonomic systems of the decoupled rows.

For each decoupling c(a,b) we search mod p for recurrences

    sum_{i=0..I} P_i(a,b) c(a+i,b) = 0      (the "a-recurrence")
    sum_{j=0..J} Q_j(a,b) c(a,b+j) = 0      (the "b-recurrence")

with P_i, Q_j of total degree <= d, by nullspace over F_p.  The leading
coefficient P_I(a,b) is the (recurrence-side) singular locus; its differential
analogue is the singular curve of the D-module.

Structure test:
  * TENSOR       c(a,b) = f(a) g(b)  -- recurrences of order 1 in each variable
  * PULLBACK     c = R * h(u): a first-order "fibration" operator exists,
                  detected as: the a-recurrence and b-recurrence have the SAME
                  order = rank of the pulled-back system.
  * COUPLED      neither.
"""
from lib2v import *
from math import comb
import sys

P = 1000003


def build(cfun, A, B):
    return [[cfun(a, b) % P for b in range(B+1)] for a in range(A+1)]


def nullspace_modp(rows, ncols, p=P):
    """Return a basis of the nullspace of the matrix `rows` over F_p."""
    M = [r[:] for r in rows]
    nrows = len(M)
    piv = []
    r = 0
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
    return basis


def find_rec(tab, I, d, axis, A, B, shift=1):
    """axis=0: recurrence in a (shifts c(a+i,b)); axis=1: in b.
    monomials a^s b^t with s+t<=d."""
    mons = [(s, t) for s in range(d+1) for t in range(d+1) if s+t <= d]
    ncols = (I+1)*len(mons)
    rows = []
    amax = A - (I if axis == 0 else 0)
    bmax = B - (I if axis == 1 else 0)
    for a in range(shift, amax+1):
        for b in range(shift, bmax+1):
            row = []
            for i in range(I+1):
                val = tab[a+i][b] if axis == 0 else tab[a][b+i]
                for (s, t) in mons:
                    row.append((val * pow(a, s, P) * pow(b, t, P)) % P)
            rows.append(row)
    if len(rows) < ncols + 20:
        return None, mons
    return nullspace_modp(rows, ncols), mons


def minimal_rec(tab, axis, A, B, maxI=4, maxd=7, name=""):
    for I in range(1, maxI+1):
        for d in range(0, maxd+1):
            ns, mons = find_rec(tab, I, d, axis, A, B)
            if ns:
                return I, d, len(ns), ns, mons
    return None, None, None, None, None


def poly_str(vec, mons, I, idx):
    """print the coefficient polynomial of shift idx"""
    nm = len(mons)
    seg = vec[idx*nm:(idx+1)*nm]
    terms = []
    for (c, (s, t)) in zip(seg, mons):
        if c % P == 0:
            continue
        cc = c if c < P//2 else c - P
        terms.append("%+d*a^%d*b^%d" % (cc, s, t))
    return " ".join(terms) if terms else "0"


CASES = [
    ("zeta3 D1  sum C(a,k)C(b,k)C(a+k,k)C(b+k,k)", z3_D1),
    ("zeta3 D2  sum C(a,k)^2 C(b+k,k)^2", z3_D2),
    ("zeta3 D3  sum C(a,k)^2 C(a+k,k) C(b+k,k)", z3_D3),
    ("zeta3 D4  sum C(a,k)C(b,k)C(a+k,k)^2", z3_D4),
    ("zeta2 D1  sum C(a,k)C(2k,k)C(b+k,2k)", z2_D1),
    ("zeta2 D2  sum C(a,k)C(b,k)C(a+k,k)", z2_D2),
    ("zeta2 D3  sum C(a,k)^2 C(b+k,k)", z2_D3),
    ("rowE  D1  sum C(a,k)C(2k,k)C(2b-2k,b-k)", E_D1),
    ("s10   D1  sum C(a,k)^2 C(b,k)^2", s10_D1),
    ("s10   D2  sum C(a,k)^3 C(b,k)", s10_D2),
]

A = B = 30
for name, f in CASES:
    tab = build(f, A, B)
    Ia, da, na, nsa, monsa = minimal_rec(tab, 0, A, B, name=name)
    Ib, db, nb, nsb, monsb = minimal_rec(tab, 1, A, B, name=name)
    print("-"*78)
    print(name)
    print("   a-recurrence: order %s, min coeff degree %s, nullspace dim %s"
          % (Ia, da, na))
    if nsa:
        print("     leading coeff P_%d(a,b) = %s" % (Ia, poly_str(nsa[0], monsa, Ia, Ia)))
        print("     trailing coeff P_0(a,b) = %s" % poly_str(nsa[0], monsa, Ia, 0))
    print("   b-recurrence: order %s, min coeff degree %s, nullspace dim %s"
          % (Ib, db, nb))
    if nsb:
        print("     leading coeff Q_%d(a,b) = %s" % (Ib, poly_str(nsb[0], monsb, Ib, Ib)))
    sys.stdout.flush()
