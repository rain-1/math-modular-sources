#!/usr/bin/env python3
"""
Exhaustive (box-free) scan.

For a parameter t and a weight w, the candidate rows are a_n(F) with
F = F0 + sum_j m_j g_j ranging over a coset of the integral lattice of weight-w
forms on Gamma_1(N).  a_n(F) is LINEAR in F, so the condition

    sum_{i<=r, d<=D} c_{i,d} n^d a_{n+i}(F) = 0     (a recurrence of length r+1,
                                                     polynomial degree D)

is BILINEAR in (c, m).  Linearising with u_{i,d,j} = c_{i,d} m_j (m_0 := 1) turns
it into a single linear system; its nullspace is computed mod a prime, and the
rank-one points of the nullspace (which are exactly the genuine (c,m) pairs) are
extracted with a matrix-pencil computation.  No search box is involved: for every
(t,w,r,D) the answer is complete over the sublattice actually used.
"""
import sys, os, json, time, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from qser import etaquo_series

Q = 99999989                # scan prime (< 2.1e8 keeps int64 convolutions exact)
NSDEF = 205

# ---------------------------------------------------------------- linear algebra
def rref_nullspace(M, q):
    """M: (rows, cols) int64 mod q. Returns basis of the right nullspace as (s, cols)."""
    M = M % q
    rows, cols = M.shape
    piv = []
    r = 0
    for c in range(cols):
        if r >= rows: break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0: continue
        i = r + nz[0]
        if i != r: M[[r, i]] = M[[i, r]]
        inv = pow(int(M[r, c]), q-2, q)
        M[r] = (M[r]*inv) % q
        col = M[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr):
            M[nzr] = (M[nzr] - np.outer(col[nzr], M[r])) % q
        piv.append(c); r += 1
    free = [c for c in range(cols) if c not in piv]
    ns = np.zeros((len(free), cols), dtype=np.int64)
    for k, fc in enumerate(free):
        ns[k, fc] = 1
        for i, pc in enumerate(piv):
            ns[k, pc] = (-M[i, fc]) % q
    return ns


# ------------------------------------------------- polynomial roots mod a prime
def pstrip(a):
    while a and a[-1] == 0: a.pop()
    return a

def pdivmod(a, b, q):
    a = a[:]; db = len(b)-1
    if db < 0: raise ZeroDivisionError
    inv = pow(b[-1], q-2, q)
    out = [0]*max(0, len(a)-db)
    for i in range(len(a)-1, db-1, -1):
        c = a[i]*inv % q
        out[i-db] = c
        if c:
            for j in range(db+1):
                a[i-db+j] = (a[i-db+j] - c*b[j]) % q
    return out, pstrip(a)

def pmul(a, b, q):
    if not a or not b: return []
    r = [0]*(len(a)+len(b)-1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                r[i+j] = (r[i+j] + ai*bj) % q
    return pstrip(r)

def pmulmod(a, b, f, q):
    return pdivmod(pmul(a, b, q), f, q)[1]

def ppowmod(a, e, f, q):
    r = [1]; a = pdivmod(a, f, q)[1]
    while e:
        if e & 1: r = pmulmod(r, a, f, q)
        a = pmulmod(a, a, f, q); e >>= 1
    return r

def pgcd(a, b, q):
    a = pstrip(a[:]); b = pstrip(b[:])
    while b:
        a, b = b, pdivmod(a, b, q)[1]
    if a:
        inv = pow(a[-1], q-2, q)
        a = [x*inv % q for x in a]
    return a

def psub(a, b, q):
    n = max(len(a), len(b)); r = [0]*n
    for i in range(n):
        r[i] = ((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % q
    return pstrip(r)

def poly_roots(coeffs, q):
    f = pstrip([c % q for c in coeffs])
    if len(f) <= 1: return []
    g = pgcd(f, psub(ppowmod([0, 1], q, f, q), [0, 1], q), q)
    res = []
    import random
    rnd = random.Random(7)
    def split(h):
        if len(h) <= 1: return
        if len(h) == 2:
            res.append((-h[0]) * pow(h[1], q-2, q) % q); return
        for _ in range(60):
            a = rnd.randrange(q)
            t = psub(ppowmod([a, 1], (q-1)//2, h, q), [1], q)
            d = pgcd(h, t, q)
            if 0 < len(d)-1 < len(h)-1:
                split(d); split(pdivmod(h, d, q)[0]); return
        return
    split(g)
    return res

def ratrecon(a, q, bound=None):
    """rational reconstruction of a mod q with |num|,|den| <= bound."""
    if bound is None: bound = int((q//2)**0.5)
    r0, r1 = q, a % q
    s0, s1 = 0, 1
    while r1 > bound:
        qq = r0//r1
        r0, r1 = r1, r0 - qq*r1
        s0, s1 = s1, s0 - qq*s1
    if s1 == 0 or abs(s1) > bound: return None
    if s1 < 0: r1, s1 = -r1, -s1
    return (r1, s1)

_LAMS = None
def all_lams(q):
    global _LAMS
    if _LAMS is None: _LAMS = np.arange(q, dtype=np.int64)
    return _LAMS

def poly_roots_table(coeffs, q):
    return poly_roots([int(c) for c in coeffs], q)

def rank_mod(M, q):
    M = M % q
    rows, cols = M.shape
    r = 0
    for c in range(cols):
        if r >= rows: break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0: continue
        i = r + nz[0]
        if i != r: M[[r, i]] = M[[i, r]]
        inv = pow(int(M[r, c]), q-2, q)
        M[r] = (M[r]*inv) % q
        col = M[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr): M[nzr] = (M[nzr] - np.outer(col[nzr], M[r])) % q
        r += 1
    return r

def recon_vec(vals, q):
    out = []
    for v in vals:
        rr = ratrecon(v, q, bound=6000)
        if rr is None: return None
        out.append(rr)
    return tuple(out)

def rank1_points(NS, U0, K1, q):
    """NS: (s, U0*K1) nullspace basis. Return list of m-vectors (length K1-1) for
    every rank-one point X = c (x) (1,m) in the span."""
    s = NS.shape[0]
    if s == 0: return []
    X = NS.reshape(s, U0, K1)
    out = []
    if K1 == 1:
        # no free directions: any nullspace vector works, m = ()
        return [()] if s > 0 else []
    if s == 1:
        M = X[0]
        if rank_mod(M.copy(), q) == 1:
            col0 = M[:, 0]
            nz = np.nonzero(col0)[0]
            if len(nz) == 0: return []
            k = nz[0]; inv = pow(int(col0[k]), q-2, q)
            m = recon_vec([int(M[k, j]*inv % q) for j in range(1, K1)], q)
            return [m] if m is not None else []
        return []
    # pencil method: need x with columns of X(x) all proportional.
    A0 = X[:, :, 0].T % q            # (U0, s)
    A1 = X[:, :, 1].T % q            # (U0, s)
    # find lambda with rank(A1 - lambda A0) < s  via an s x s minor determinant
    rng = np.random.default_rng(12345)
    cand = None
    for _ in range(6):
        idx = rng.choice(U0, size=min(s, U0), replace=False)
        B0 = A0[idx]; B1 = A1[idx]
        if B0.shape[0] < s: return []
        # det(B1 - l B0) as a polynomial in l : interpolate at s+1 points
        pts = list(range(s+1))
        vals = []
        for l in pts:
            Mm = (B1 - l*B0) % q
            vals.append(det_mod(Mm, q))
        co = lagrange_coeffs(pts, vals, q)
        if any(c % q for c in co):
            cand = co; break
    if cand is None: return []
    roots = poly_roots_table(cand, q)
    for lam in roots:
        Mm = (A1 - int(lam)*A0) % q
        ns2 = rref_nullspace(Mm.copy(), q)
        for x in ns2:
            Xx = np.tensordot(x, X, axes=(0, 0)) % q
            if rank_mod(Xx.copy(), q) != 1: continue
            col0 = Xx[:, 0]
            nz = np.nonzero(col0)[0]
            if len(nz) == 0: continue
            k = nz[0]; inv = pow(int(col0[k]), q-2, q)
            m = recon_vec([int(Xx[k, j]*inv % q) for j in range(1, K1)], q)
            if m is not None and m not in out: out.append(m)
    return out

def det_mod(M, q):
    M = M.copy() % q
    n = M.shape[0]
    det = 1
    for c in range(n):
        nz = np.nonzero(M[c:, c])[0]
        if len(nz) == 0: return 0
        i = c + nz[0]
        if i != c:
            M[[c, i]] = M[[i, c]]; det = (-det) % q
        det = det*int(M[c, c]) % q
        inv = pow(int(M[c, c]), q-2, q)
        M[c] = (M[c]*inv) % q
        if c+1 < n:
            col = M[c+1:, c]
            nzr = np.nonzero(col)[0]
            if len(nzr): M[c+1:][nzr] = (M[c+1:][nzr] - np.outer(col[nzr], M[c])) % q
    return det % q

def lagrange_coeffs(xs, ys, q):
    n = len(xs)
    co = [0]*n
    for i in range(n):
        # basis poly prod_{j!=i}(x - xj)/(xi - xj)
        num = [1]
        den = 1
        for j in range(n):
            if j == i: continue
            num = poly_mul(num, [(-xs[j]) % q, 1], q)
            den = den*((xs[i]-xs[j]) % q) % q
        inv = pow(den, q-2, q)
        for k in range(len(num)):
            co[k] = (co[k] + ys[i]*num[k] % q * inv) % q
    return co

def poly_mul(a, b, q):
    r = [0]*(len(a)+len(b)-1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                r[i+j] = (r[i+j] + ai*bj) % q
    return r

# ---------------------------------------------------------------- series
def tseries(D, r, ns):
    o24, s = etaquo_series(D, r, ns)
    assert o24 == 24
    T = [0]*ns
    for n in range(1, ns): T[n] = s[n-1]
    return T

def apeel(T, Bas, ns, q):
    """Bas: (nb, ns) integer q-expansions. Return a-vectors (nb, ns) mod q."""
    Tq = np.array([x % q for x in T], dtype=np.int64)
    Pm = np.zeros((ns, ns), dtype=np.int64)
    cur = np.zeros(ns, dtype=np.int64); cur[0] = 1
    Pm[0] = cur
    for n in range(1, ns):
        cur = np.convolve(cur, Tq)[:ns] % q
        Pm[n] = cur
    B = Bas % q
    A = np.zeros_like(B)
    for n in range(ns):
        v = B[:, n].copy()
        if n:
            v = (v - A[:, :n] @ Pm[:n, n]) % q
        A[:, n] = v % q
    return A, Pm

PREPDIR = "prep"

def load_prep(M, w, ch):
    fn = os.path.join(HERE, PREPDIR, f"pr_{M}_{w}_{ch}.txt")
    if not os.path.exists(fn): return None
    with open(fn) as f:
        lines = [l.strip() for l in f if l.strip()]
    h = lines[0].split()
    if int(h[4]) < 0: return None
    out = []
    for l in lines[1:]:
        try:
            v = json.loads(l.replace(" ", ""))
        except Exception:
            break
        if not isinstance(v, list): break
        out.append(v)
    return out or None

def space_index():
    """map level -> list of (w, chi) with a prep file"""
    idx = {}
    d = os.path.join(HERE, PREPDIR)
    for fn in os.listdir(d):
        if not fn.startswith("pr_"): continue
        try:
            _, M, w, ch = fn[:-4].split("_")
        except ValueError:
            continue
        idx.setdefault(int(M), []).append((int(w), int(ch)))
    for k in idx: idx[k].sort()
    return idx

SHAPES = {1: [(2,2),(2,3),(3,2),(3,3),(4,3)],
          2: [(2,3),(2,4),(3,3),(3,4),(4,4)],
          3: [(2,4),(2,5),(3,4),(3,5),(4,5)]}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--degmin", type=int, default=1)
    ap.add_argument("--degmax", type=int, default=2)
    ap.add_argument("--weights", type=str, default="1,2,3")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--ns", type=int, default=NSDEF)
    ap.add_argument("--kcap", type=int, default=24)
    ap.add_argument("--levcap", type=int, default=6)   # F-level M <= levcap * t-level
    ap.add_argument("--out", type=str, default="hits2")
    ap.add_argument("--prepdir", type=str, default="prep")
    args = ap.parse_args()
    global PREPDIR
    PREPDIR = args.prepdir

    tl = json.load(open(os.path.join(HERE, "t_list.json")))
    tl = [(i, x) for i, x in enumerate(tl) if args.degmin <= x["deg"] <= args.degmax]
    weights = set(int(x) for x in args.weights.split(","))
    SIDX = space_index()
    fout = open(os.path.join(HERE, f"{args.out}_{args.shard}.jsonl"), "a")
    ns = args.ns
    t0 = time.time(); nsolve = 0; nhit = 0
    prep_cache = {}
    for cnt, (ti, tt) in enumerate(tl):
        if cnt % args.nshard != args.shard: continue
        N = tt["N"]
        levels = [M for M in SIDX if M % N == 0 and M <= args.levcap*N]
        if not levels: continue
        T = tseries(tt["D"], tt["r"], ns)
        Aq = None
        for M in sorted(levels):
            for (w, ch) in SIDX[M]:
                if w not in weights: continue
                key = (M, w, ch)
                if key not in prep_cache: prep_cache[key] = load_prep(M, w, ch)
                pr = prep_cache[key]
                if pr is None: continue
                Bas = np.array([[int(x) % Q for x in row[:ns]] for row in pr], dtype=np.int64)
                A, _ = apeel(T, Bas, ns, Q)
                found = {}
                for (r, D) in SHAPES[w]:
                    U0 = (r+1)*(D+1)
                    K1max = (ns - r - 10)//U0
                    K1 = min(A.shape[0], K1max, args.kcap+1)
                    if K1 < 1: continue
                    NEQ = U0*K1 + 6
                    if NEQ + r >= ns: NEQ = ns - r - 1
                    if NEQ < U0*K1: continue
                    nn = np.arange(1, NEQ+1, dtype=np.int64)
                    npow = np.ones((NEQ, D+1), dtype=np.int64)
                    for d in range(1, D+1): npow[:, d] = (npow[:, d-1]*nn) % Q
                    sh = np.stack([A[:K1, 1+i:1+i+NEQ] for i in range(r+1)], axis=0)
                    Mx = np.zeros((NEQ, U0*K1), dtype=np.int64)
                    for i in range(r+1):
                        for d in range(D+1):
                            blk = (sh[i].T * npow[:, d:d+1]) % Q
                            Mx[:, (i*(D+1)+d)*K1:(i*(D+1)+d)*K1+K1] = blk
                    NSb = rref_nullspace(Mx, Q)
                    nsolve += 1
                    if NSb.shape[0] == 0: continue
                    for m in rank1_points(NSb, U0, K1, Q):
                        if m in found: continue
                        found[m] = (r, D)
                        fout.write(json.dumps({"ti": ti, "N": N, "M": M, "chi": ch,
                                               "rq": tt["r"], "Dq": tt["D"],
                                               "deg": tt["deg"], "w": w, "r": r, "Dg": D,
                                               "K1": K1, "m": [[int(x[0]), int(x[1])] for x in m]}) + "\n")
                        nhit += 1
                    fout.flush()
        if cnt % 25 == 0:
            print(f"[{time.time()-t0:.0f}s] {cnt}/{len(tl)} N={N} solves={nsolve} hits={nhit}", flush=True)
    print(f"DONE solves={nsolve} hits={nhit} time={time.time()-t0:.0f}s", flush=True)

if __name__ == "__main__":
    main()
