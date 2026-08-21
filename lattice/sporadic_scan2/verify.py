#!/usr/bin/env python3
"""Exact verification + invariants for the scan hits.

For each (t, w, F):
  * recompute a_n exactly (big integers) to n = NV
  * fit the minimal (r,D) recurrence exactly over Q and verify it for all n <= 200
  * characteristic roots lambda_i (roots of the leading symbol), lambda_1, lambda_2, c
  * companion  B = F * D^{-(w+1)}(F * Dt),  b_n,  sharp denominator exponent k
  * p-adic slopes sigma_p, Apery limit to high precision
"""
import sys, os, json, math, argparse
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qser import etaquo_series, divisors

NV = 210

def tser(D, r, P):
    o24, s = etaquo_series(D, r, P)
    assert o24 == 24
    T = [0]*P
    for n in range(1, P): T[n] = s[n-1]
    return T

def peel(T, F, P):
    """a_n with F = sum a_n t^n (exact integers)."""
    tp = [[0]*P for _ in range(P)]
    tp[0][0] = 1
    for n in range(1, P):
        prev = tp[n-1]; cur = tp[n]
        for i in range(P):
            pi = prev[i]
            if pi:
                for j in range(1, P-i):
                    if T[j]: cur[i+j] += pi*T[j]
    rem = list(F); a = []
    for n in range(P):
        an = rem[n]; a.append(an)
        if an:
            row = tp[n]
            for i in range(n, P):
                if row[i]: rem[i] -= an*row[i]
    return a, tp

def expand_in_t(G, tp, P):
    rem = list(G); b = []
    for n in range(P):
        bn = rem[n]; b.append(bn)
        if bn:
            row = tp[n]
            for i in range(n, P):
                if row[i]: rem[i] -= bn*row[i]
    return b

def smul(a, b, P):
    r = [0]*P
    for i, ai in enumerate(a):
        if ai:
            for j in range(min(len(b), P-i)):
                if b[j]: r[i+j] += ai*b[j]
    return r

def fit_recurrence(a, r, D, nmax):
    """solve sum_{i<=r,d<=D} c_{i,d} n^d a_{n+i} = 0 exactly; return coefficient
    matrix c (list of lists) or None."""
    U = (r+1)*(D+1)
    rows = []
    for n in range(1, min(nmax, len(a)-r-1)+1):
        row = []
        for i in range(r+1):
            an = a[n+i]
            v = 1
            for d in range(D+1):
                row.append(Fraction(v*an)); v *= n
        rows.append(row)
        if len(rows) >= U + 8: break
    # gaussian elimination over Q
    M = [row[:] for row in rows]
    nr = len(M)
    piv = []
    rr = 0
    for c in range(U):
        p = None
        for i in range(rr, nr):
            if M[i][c] != 0: p = i; break
        if p is None: continue
        M[rr], M[p] = M[p], M[rr]
        pv = M[rr][c]
        M[rr] = [x/pv for x in M[rr]]
        for i in range(nr):
            if i != rr and M[i][c] != 0:
                f = M[i][c]
                M[i] = [x - f*y for x, y in zip(M[i], M[rr])]
        piv.append(c); rr += 1
        if rr == U: break
    if rr == U: return None
    free = [c for c in range(U) if c not in piv]
    sols = []
    for fc in free:
        v = [Fraction(0)]*U
        v[fc] = Fraction(1)
        for i, pc in enumerate(piv):
            v[pc] = -M[i][fc]
        sols.append(v)
    return sols

def check_rec(a, c, r, D, nmax):
    for n in range(1, nmax+1):
        if n + r >= len(a): break
        s = Fraction(0)
        for i in range(r+1):
            v = 1
            for d in range(D+1):
                s += c[i*(D+1)+d]*v*a[n+i]
                v *= n
        if s != 0: return n-1
    return nmax

def charroots(c, r, D):
    lead = [c[i*(D+1)+D] for i in range(r+1)]
    while lead and lead[-1] == 0: lead.pop()
    if len(lead) < 2: return None
    import numpy as np
    co = [float(x) for x in lead][::-1]
    rts = np.roots(co)
    return sorted([complex(z) for z in rts], key=lambda z: -abs(z))

def vp(x, p):
    if x == 0: return 10**9
    n = 0
    if isinstance(x, Fraction):
        return vp(x.numerator, p) - vp(x.denominator, p)
    x = abs(x)
    while x % p == 0:
        x //= p; n += 1
    return n

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", nargs="+", required=True)
    ap.add_argument("--out", default="rows.jsonl")
    ap.add_argument("--nv", type=int, default=NV)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0)
    args = ap.parse_args()
    HERE = os.path.dirname(os.path.abspath(__file__))
    P = args.nv
    tl = json.load(open(os.path.join(HERE, "t_list.json")))
    prep = {}
    seen = set()
    recs = []
    for fn in args.hits:
        for line in open(fn):
            h = json.loads(line)
            if any(x[1] != 1 for x in h["m"]): continue     # F must be integral
            key = (h["ti"], h["w"], tuple(x[0] for x in h["m"]))
            if key in seen: continue
            seen.add(key)
            recs.append(h)
    print(f"{len(recs)} integral candidates", flush=True)
    fout = open(os.path.join(HERE, args.out), "a")
    for idx, h in enumerate(recs):
        if idx % args.nshard != args.shard: continue
        N, w, ti = h["N"], h["w"], h["ti"]
        tt = tl[ti]
        key = (N, w)
        if key not in prep:
            pf = os.path.join(HERE, "prep", f"pr_{N}_{w}.txt")
            lines = [l.strip() for l in open(pf) if l.strip()]
            prep[key] = [json.loads(l.replace(" ", "")) for l in lines[1:]]
        vecs = prep[key]
        F = [int(x) for x in vecs[0][:P]]
        for j, mv in enumerate(h["m"]):
            if mv[0]:
                kv = vecs[1+j]
                for i in range(P): F[i] += mv[0]*kv[i]
        T = tser(tt["D"], tt["r"], P)
        a, tp = peel(T, F, P)
        # minimal recurrence
        best = None
        for D in range(1, w+3):
            for r in range(2, 6):
                if (r+1)*(D+1) + 10 > P - r: continue
                sols = fit_recurrence(a, r, D, P-r-2)
                if not sols: continue
                c = sols[0]
                nok = check_rec(a, c, r, D, min(200, P-r-2))
                if nok >= min(200, P-r-2):
                    best = (r, D, c, len(sols), nok); break
            if best: break
        if not best:
            continue
        r, D, c, nsol, nok = best
        rts = charroots(c, r, D)
        if rts is None: continue
        lam1 = abs(rts[0]); lam2 = abs(rts[1]) if len(rts) > 1 else 0.0
        lead = [c[i*(D+1)+D] for i in range(r+1)]
        while lead and lead[-1] == 0: lead.pop()
        cc = float(lead[0]/lead[-1])*((-1)**(len(lead)-1))
        # companion
        Dt = [n*T[n] for n in range(P)]
        Phi = smul(F, Dt, P)
        G = Phi[:]
        for _ in range(w+1):
            G = [Fraction(G[n], n) if n else Fraction(0) for n in range(P)]
        Bq = [sum(Fraction(F[i])*G[n-i] for i in range(n+1)) for n in range(P)]
        b = expand_in_t(Bq, tp, P)
        # denominator exponent
        dn = [1]*(P)
        for n in range(1, P): dn[n] = dn[n-1]*n//math.gcd(dn[n-1], n)
        kk = None
        for ktry in range(0, 7):
            ok = all((b[n]*dn[n]**ktry).denominator == 1 for n in range(1, min(P, 201)))
            if ok: kk = ktry; break
        # slopes
        slopes = {}
        for p in (2, 3, 5, 7):
            vals = []
            for n in (100, 150, 200):
                if n+1 >= P: continue
                x = Fraction(a[n])*b[n+1] - Fraction(a[n+1])*b[n]
                if x == 0: vals.append(None); continue
                vals.append(vp(x, p))
            if len(vals) >= 2 and all(v is not None for v in vals):
                slopes[p] = round((vals[-1]-vals[0])/(200-100), 2)
        # apery limit (numeric)
        lim = None
        try:
            import mpmath as mp
            mp.mp.dps = 120
            n0 = min(P-2, 200)
            lim = mp.mpf(b[n0].numerator)/mp.mpf(b[n0].denominator)/mp.mpf(a[n0]) if a[n0] else None
            lim = mp.nstr(lim, 80) if lim is not None else None
        except Exception as e:
            lim = None
        out = {"ti": ti, "N": N, "w": w, "m": [x[0] for x in h["m"]],
               "etaq_t": tt["r"], "divs": tt["D"], "tdeg": tt["deg"],
               "r": r, "D": D, "nsol": nsol, "nok": nok,
               "lam1": lam1, "lam2": lam2, "c": cc, "k": kk,
               "slopes": slopes, "a": [str(x) for x in a[:12]],
               "rec": [str(x) for x in c], "limit": lim,
               "score": (math.log(1/lam2) - kk) if (kk is not None and lam2 > 0) else None,
               "budget": (math.log(lam1) - kk) if kk is not None else None}
        fout.write(json.dumps(out) + "\n"); fout.flush()
        if idx % 10 == 0: print(f"{idx}/{len(recs)} N={N} w={w} r={r} D={D} lam1={lam1:.4f}", flush=True)
    print("done", flush=True)

if __name__ == "__main__":
    main()
