"""Helper library for the 10_* series: fast mod-p nullspaces (numpy), exact
two-variable recurrence fitting (CRT + rational reconstruction), and
propagation / denominator utilities."""
import numpy as np
from fractions import Fraction as F
from math import comb, gcd, isqrt
import sys

PRIMES = [1000003, 999983, 999979, 999961, 999959]

# ---------------- mod p linear algebra ----------------

def nullspace_modp_np(M, p):
    """M: numpy int64 array (nrows x ncols) already reduced mod p.
    Returns list of nullspace basis vectors (python lists of ints)."""
    M = M.astype(np.int64) % p
    nrows, ncols = M.shape
    piv = []
    r = 0
    for c in range(ncols):
        if r >= nrows:
            break
        col = M[r:, c]
        nz = np.nonzero(col)[0]
        if nz.size == 0:
            continue
        pr = r + int(nz[0])
        if pr != r:
            M[[r, pr]] = M[[pr, r]]
        inv = pow(int(M[r, c]), p-2, p)
        M[r] = (M[r] * inv) % p
        colvals = M[:, c].copy()
        colvals[r] = 0
        nzr = np.nonzero(colvals)[0]
        if nzr.size:
            M[nzr] = (M[nzr] - np.outer(colvals[nzr], M[r])) % p
        piv.append(c)
        r += 1
    pivset = set(piv)
    free = [c for c in range(ncols) if c not in pivset]
    basis = []
    for fc in free:
        v = [0]*ncols
        v[fc] = 1
        for i, pc in enumerate(piv):
            v[pc] = int((-M[i, fc]) % p)
        basis.append(v)
    return basis, piv, free


def build_rows(tab, I, d, axis, A, B, p, amin=0, bmin=0):
    """Rows of the ansatz  sum_i P_i(a,b) c(a+i,b) = 0  (axis=0)
    with P_i = sum_{s+t<=d} coef * a^s b^t.  Returns (numpy matrix, mons, pts)."""
    mons = [(s, t) for s in range(d+1) for t in range(d+1) if s+t <= d]
    nm = len(mons)
    ncols = (I+1)*nm
    amax = A - (I if axis == 0 else 0)
    bmax = B - (I if axis == 1 else 0)
    pts = [(a, b) for a in range(amin, amax+1) for b in range(bmin, bmax+1)]
    M = np.zeros((len(pts), ncols), dtype=np.int64)
    for ri, (a, b) in enumerate(pts):
        apw = [pow(a, s, p) for s in range(d+1)]
        bpw = [pow(b, t, p) for t in range(d+1)]
        mv = [(apw[s]*bpw[t]) % p for (s, t) in mons]
        for i in range(I+1):
            val = tab[a+i][b] if axis == 0 else tab[a][b+i]
            base = i*nm
            for j in range(nm):
                M[ri, base+j] = (val*mv[j]) % p
    return M, mons, pts


def table_modp(cfun, A, B, p, axis=0, I=0):
    AA = A + (I if axis == 0 else 0)
    BB = B + (I if axis == 1 else 0)
    return [[cfun(a, b) % p for b in range(BB+1)] for a in range(AA+1)]


def scan_minrec(cfun, axis, A, B, maxI=5, maxd=8, p=PRIMES[0], amin=0, bmin=0,
                need_margin=25):
    """Find the minimal (order, degree) with a nonzero relation, verified on ALL
    sample points."""
    for I in range(1, maxI+1):
        tab = table_modp(cfun, A, B, p, axis, I)
        for d in range(0, maxd+1):
            mons = [(s, t) for s in range(d+1) for t in range(d+1) if s+t <= d]
            ncols = (I+1)*len(mons)
            M, mons, pts = build_rows(tab, I, d, axis, A, B, p, amin, bmin)
            if M.shape[0] < ncols + need_margin:
                continue
            ns, piv, free = nullspace_modp_np(M, p)
            if ns:
                # verify on all points (they already are all points)
                good = []
                for v in ns:
                    vv = np.array(v, dtype=np.int64) % p
                    res = (M.dot(vv)) % p
                    if not res.any():
                        good.append(v)
                if good:
                    return I, d, len(good), good, mons, pts
    return None, None, 0, None, None, None


# ---------------- exact reconstruction ----------------

def ratrec(a, m):
    bound = isqrt(m // 2)
    r0, r1 = m, a % m
    s0, s1 = 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q*r1
        s0, s1 = s1, s0 - q*s1
    if s1 == 0 or gcd(r1, abs(s1)) != 1:
        return None
    return F(r1 if s1 > 0 else -r1, abs(s1))


def exact_rec(cfun, I, d, axis, A, B, nprimes=3, amin=0, bmin=0):
    """Exact (over Q) recurrence coefficients; returns (list of Fractions, mons)
    normalised to primitive integer vector."""
    mons = [(s, t) for s in range(d+1) for t in range(d+1) if s+t <= d]
    ncols = (I+1)*len(mons)
    vecs = []
    freeidx = None
    for p in PRIMES[:nprimes]:
        tab = table_modp(cfun, A, B, p, axis, I)
        M, _, pts = build_rows(tab, I, d, axis, A, B, p, amin, bmin)
        ns, piv, free = nullspace_modp_np(M, p)
        if len(ns) != 1:
            return None, mons, len(ns)
        if freeidx is None:
            freeidx = free[0]
        elif freeidx != free[0]:
            return None, mons, -1
        vecs.append((p, ns[0]))
    Mod = 1
    for p, _ in vecs:
        Mod *= p
    out = []
    for j in range(ncols):
        x = 0
        for p, v in vecs:
            Mp = Mod//p
            x = (x + v[j]*Mp*pow(Mp, -1, p)) % Mod
        rr = ratrec(x, Mod)
        if rr is None:
            return None, mons, -2
        out.append(rr)
    # normalise: multiply by lcm of denominators, divide by gcd of numerators
    L = 1
    for x in out:
        L = L*x.denominator//gcd(L, x.denominator)
    ints = [int(x*L) for x in out]
    g = 0
    for x in ints:
        g = gcd(g, abs(x))
    if g:
        ints = [x//g for x in ints]
    # sign: make leading coefficient's first nonzero positive
    for x in reversed(ints):
        if x:
            if x < 0:
                ints = [-y for y in ints]
            break
    return ints, mons, 1


def polyval(seg, mons, a, b):
    return sum(c*a**s*b**t for c, (s, t) in zip(seg, mons) if c)


def poly_sympy(seg, mons):
    import sympy as sp
    a, b = sp.symbols('a b')
    e = sum(sp.Integer(c)*a**s*b**t for c, (s, t) in zip(seg, mons) if c)
    return sp.expand(e)


def verify_rec(cfun, ints, mons, I, axis, A, B, amin=0, bmin=0):
    nm = len(mons)
    bad = []
    amax = A - (I if axis == 0 else 0)
    bmax = B - (I if axis == 1 else 0)
    for a in range(amin, amax+1):
        for b in range(bmin, bmax+1):
            s = 0
            for i in range(I+1):
                val = cfun(a+i, b) if axis == 0 else cfun(a, b+i)
                s += polyval(ints[i*nm:(i+1)*nm], mons, a, b)*val
            if s != 0:
                bad.append((a, b))
    return bad


def lcmrange(n):
    r = 1
    for i in range(1, n+1):
        r = r*i//gcd(r, i)
    return r


# ---------------- grid solution space of a pair of recurrences ----------------

BIGPRIMES = [2147483647, 2147483629, 2147483587, 2147483579, 2147483563,
             2147483549, 2147483543, 2147483497, 2147483489, 2147483477,
             2147483423, 2147483399]


class Rec:
    """one recurrence:  sum_i P_i(a,b) t[a+i][b] = 0   (axis 0)
                        sum_j Q_j(a,b) t[a][b+j] = 0   (axis 1)"""
    def __init__(self, I, mons, coefs, axis):
        self.I, self.mons, self.axis = I, [tuple(m) for m in mons], axis
        nm = len(self.mons)
        self.segs = [coefs[i*nm:(i+1)*nm] for i in range(I+1)]

    def P(self, i, a, b):
        return sum(c*a**s*b**t for c, (s, t) in zip(self.segs[i], self.mons) if c)

    def cells(self, a, b):
        if self.axis == 0:
            return [(a+i, b) for i in range(self.I+1)]
        return [(a, b+i) for i in range(self.I+1)]


def grid_index(NA, NB):
    cells = [(a, b) for a in range(NA+1) for b in range(NB+1)]
    cells.sort(key=lambda ab: (-(ab[0]+ab[1]), -ab[0]))
    return cells, {c: i for i, c in enumerate(cells)}


def grid_equations(reca, recb, NA, NB):
    eqs = []
    for a in range(NA - reca.I + 1):
        for b in range(NB+1):
            eqs.append([((a+i, b), reca.P(i, a, b)) for i in range(reca.I+1)])
    for a in range(NA+1):
        for b in range(NB - recb.I + 1):
            eqs.append([((a, b+j), recb.P(j, a, b)) for j in range(recb.I+1)])
    return eqs


def grid_nullspace_modp(reca, recb, NA, NB, p):
    cells, idx = grid_index(NA, NB)
    eqs = grid_equations(reca, recb, NA, NB)
    M = np.zeros((len(eqs), len(cells)), dtype=np.int64)
    for r, eq in enumerate(eqs):
        for cell, co in eq:
            M[r, idx[cell]] = co % p
    return nullspace_modp_np(M, p) + (cells, idx)


def grid_nullspace_exact(reca, recb, NA, NB, nprimes=6):
    """Exact rational nullspace basis (rref form), via CRT of the mod-p rrefs."""
    cells, idx = grid_index(NA, NB)
    eqs = grid_equations(reca, recb, NA, NB)
    rows = []
    for eq in eqs:
        rows.append([(idx[cell], co) for cell, co in eq])
    ncols = len(cells)
    vecs, freeref = [], None
    for p in BIGPRIMES[:nprimes]:
        M = np.zeros((len(rows), ncols), dtype=np.int64)
        for r, eq in enumerate(rows):
            for j, co in eq:
                M[r, j] = co % p
        ns, piv, free = nullspace_modp_np(M, p)
        if freeref is None:
            freeref = free
        elif freeref != free:
            return None, None, None, "free-set mismatch at p=%d" % p
        vecs.append((p, ns))
    Mod = 1
    for p, _ in vecs:
        Mod *= p
    D = len(freeref)
    out = []
    for k in range(D):
        v = []
        for j in range(ncols):
            x = 0
            for p, ns in vecs:
                Mp = Mod//p
                x = (x + ns[k][j]*Mp*pow(Mp, -1, p)) % Mod
            rr = ratrec(x, Mod)
            if rr is None:
                return None, None, None, "ratrec failed at basis %d col %d" % (k, j)
            v.append(rr)
        out.append(v)
    return out, cells, freeref, None


def check_solution(reca, recb, T, NA, NB):
    """T: dict (a,b)->Fraction. Returns list of violated equations."""
    bad = []
    for rec in (reca, recb):
        for a in range(NA+1):
            for b in range(NB+1):
                cs = rec.cells(a, b)
                if max(c[0] for c in cs) > NA or max(c[1] for c in cs) > NB:
                    continue
                s = sum(rec.P(i, a, b)*T[c] for i, c in enumerate(cs))
                if s != 0:
                    bad.append((rec.axis, a, b))
    return bad


def propagate(reca, recb, T, NA, NB, maxsweep=200):
    """Fill in unknown cells of dict T using either recurrence wherever exactly
    one cell of a relation is unknown and its coefficient is nonzero."""
    for _ in range(maxsweep):
        changed = 0
        for rec in (reca, recb):
            for a in range(NA+1):
                for b in range(NB+1):
                    cs = rec.cells(a, b)
                    if max(c[0] for c in cs) > NA or max(c[1] for c in cs) > NB:
                        continue
                    unk = [k for k, c in enumerate(cs) if c not in T]
                    if len(unk) != 1:
                        continue
                    k = unk[0]
                    lead = rec.P(k, a, b)
                    if lead == 0:
                        continue
                    s = sum(rec.P(i, a, b)*T[c] for i, c in enumerate(cs) if i != k)
                    T[cs[k]] = F(-s, 1)/F(lead)
                    changed += 1
        if not changed:
            break
    missing = [(a, b) for a in range(NA+1) for b in range(NB+1) if (a, b) not in T]
    return T, missing
