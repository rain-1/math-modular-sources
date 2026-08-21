#!/usr/bin/env python3
"""Fast full scan. See SPORADIC_SEARCH.md for method."""
import sys, os, math, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from search import *
from run_scan import hauptmodul_candidates, form_candidates, etaq
import run_scan
run_scan.STRICT_HAUPTMODUL = False

PS = 28          # scan precision
PH = 66          # verify precision
DEG = 3
NU = 3*(DEG+1)   # 12 unknowns
PRIME = (1 << 61) - 1

def peel(tpow, F, P):
    """a_n with F = sum a_n t^n; tpow[n] = t^n (integer series)."""
    rem = list(F); a = []
    for n in range(P):
        an = rem[n]; a.append(an)
        if an:
            tn = tpow[n]
            for i in range(n, P):
                if tn[i]: rem[i] -= an*tn[i]
    return a

def rank_deficient_modp(a, nrows):
    p = PRIME
    M = []
    for n in range(1, nrows+1):
        row = []
        for val in (a[n+1], a[n], a[n-1]):
            v = val % p; nk = 1
            for j in range(DEG+1):
                row.append(v*nk % p); nk = nk*n % p
        M.append(row)
    # gaussian elimination
    rank = 0
    for col in range(NU):
        piv = None
        for i in range(rank, len(M)):
            if M[i][col] % p: piv = i; break
        if piv is None: continue
        M[rank], M[piv] = M[piv], M[rank]
        inv = pow(M[rank][col], p-2, p)
        for i in range(rank+1, len(M)):
            f = M[i][col]
            if f:
                f = f*inv % p
                Mi = M[i]; Mr = M[rank]
                for c in range(col, NU): Mi[c] = (Mi[c] - f*Mr[c]) % p
        rank += 1
        if rank == NU: return False
    return rank < NU

def main(levels):
    results = []; seen = {}
    t0 = time.time(); ntot = 0; ncand = 0
    for N in levels:
        D = divisors(N); B = 24 if len(D) <= 4 else 16
        hs = hauptmodul_candidates(N, B)[1]
        fs = {w: form_candidates(N, w, B)[1] for w in (1, 2)}
        print(f"N={N} t-cands={len(hs)} F(w=1)={len(fs[1])} F(w=2)={len(fs[2])}", flush=True)
        cache = {}
        Fs = {w: [(s, etaq(D, s, PS, cache)) for s in fs[w]] for w in (1, 2)}
        for (r, cpole) in hs:
            t = [0] + etaq(D, r, PS, cache)[:PS-1]
            tpow = [[0]*PS for _ in range(PS)]; tpow[0][0] = 1
            cur = [0]*PS; cur[0] = 1
            for n in range(1, PS):
                cur = smul(cur, t, PS); tpow[n] = cur
            for w in (1, 2):
                for s, F in Fs[w]:
                    ntot += 1
                    a = peel(tpow, F, PS)
                    if not rank_deficient_modp(a, 20): continue
                    ncand += 1
                    key = tuple(a[:14])
                    if key in seen:
                        seen[key].append((N, list(r), list(s), w)); continue
                    seen[key] = [(N, list(r), list(s), w)]
                    # exact confirm at high precision
                    cH = {}
                    tH = [0] + etaq(D, r, PH, cH)[:PH-1]
                    tpH = [[0]*PH for _ in range(PH)]; tpH[0][0] = 1
                    cu = [0]*PH; cu[0] = 1
                    for n in range(1, PH):
                        cu = smul(cu, tH, PH); tpH[n] = cu
                    aH = peel(tpH, etaq(D, s, PH, cH), PH)
                    ns = fit_recurrence(aH[:26], DEG)
                    good = None
                    for v in ns:
                        den = 1
                        for x in v: den = den*x.q//math.gcd(den, x.q)
                        vv = [int(x*den) for x in v]
                        ok, pp = check_recurrence(vv, aH)
                        if ok: good = pp; break
                    if good is None: continue
                    cl = classify(*good); flip = False; aa = aH
                    if cl is None or cl[1] < 0:
                        af = [((-1)**i)*x for i, x in enumerate(aH)]
                        for v in fit_recurrence(af[:26], DEG):
                            den = 1
                            for x in v: den = den*x.q//math.gcd(den, x.q)
                            ok, pp = check_recurrence([int(x*den) for x in v], af)
                            if ok:
                                c2 = classify(*pp)
                                if c2 is not None and (cl is None or c2[1] > 0):
                                    cl, good, aa, flip = c2, pp, af, True
                                break
                    rec = dict(N=N, r=list(r), s=list(s), w=w, poles=str(cpole), flip=flip,
                               a=aa[:14], cls=[str(x) for x in cl] if cl else None,
                               poly=[[str(x) for x in q] for q in good], key=key)
                    results.append(rec)
                    print("  HIT", N, r, s, "w", w, "cls", cl, aa[:7], flush=True)
    print("scanned", ntot, "rank-deficient", ncand, "elapsed", time.time()-t0)
    for rec in results:
        rec['also'] = seen[rec.pop('key')]
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hits.json')
    json.dump(results, open(out, 'w'), indent=1)
    print("wrote", out, len(results))

if __name__ == '__main__':
    lv = [int(x) for x in sys.argv[1:]] or GENUS0
    main(lv)
