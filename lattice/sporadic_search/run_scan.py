#!/usr/bin/env python3
import sys, os, math, json, time
STRICT_HAUPTMODUL = False
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from search import *
from fractions import Fraction

ZAG2 = {(7,2,-8),(9,3,27),(10,3,9),(11,3,-1),(12,4,32),(17,6,72)}
ZAG3 = {(7,3,81),(9,3,-27),(10,4,64),(11,5,125),(12,4,16),(17,5,1)}

def hauptmodul_candidates(N, B):
    D = divisors(N); out = []
    for r in enum_vectors(D, 0, 24, B):
        if not ligozat_ok(N, D, r): continue
        co = cusp_orders(N, D, r)
        if any(x.denominator != 1 for x in co.values()): continue
        if co[N] != 1: continue
        neg = [c for c in D if co[c] < 0]
        deg = -sum(co[c] for c in neg)
        if STRICT_HAUPTMODUL and (len(neg) != 1 or co[neg[0]] != -1
                                  or any(co[c] != 0 for c in D if c != N and c != neg[0])):
            continue
        out.append((r, tuple(sorted((int(c), int(co[c])) for c in neg))))
    return D, out

def form_candidates(N, w, B):
    D = divisors(N); out = []
    for s in enum_vectors(D, 2*w, 0, B):
        if not ligozat_ok(N, D, s): continue
        co = cusp_orders(N, D, s)
        if any(x.denominator != 1 for x in co.values()): continue
        if co[N] != 0: continue
        if any(x < 0 for x in co.values()): continue
        out.append(s)
    return D, out

def etaq(D, exps, P, cache):
    res = [0]*P; res[0] = 1
    for d, e in zip(D, exps):
        if e == 0: continue
        if (d, P) not in cache: cache[(d, P)] = eta_factor(d, P)
        res = smul(res, spow(cache[(d, P)], e, P), P)
    return res

def revert_coeffs(tser, Fser, P):
    """a_n with F = sum a_n t^n, t = q + O(q^2). Returns a[0..P-1] as Fractions/ints."""
    a = []
    rem = [Fraction(x) for x in Fser]
    tp = [Fraction(0)]*P; tp[0] = Fraction(1)   # t^0
    tpow = [tp]
    cur = [Fraction(x) for x in tser]
    for n in range(1, P):
        tpow.append(cur)
        cur = smul(cur, [Fraction(x) for x in tser], P)
    for n in range(P):
        an = rem[n]
        a.append(an)
        if an:
            tn = tpow[n]
            for i in range(n, P):
                rem[i] -= an*tn[i]
    return a

def main():
    B_small, B_big = 24, 16
    results = []
    seen = {}
    t0 = time.time()
    for N in GENUS0:
        D = divisors(N)
        B = B_small if len(D) <= 4 else B_big
        _, hs = hauptmodul_candidates(N, B)
        fs = {}
        for w in (1, 2):
            _, fs[w] = form_candidates(N, w, B)
        print(f"N={N} divisors={D} hauptmoduls={len(hs)} forms w=1:{len(fs[1])} w=2:{len(fs[2])}", flush=True)
        cache = {}
        P = PREC
        for (r, cpole) in hs:
            tser = [0] + etaq(D, r, P, cache)[:P-1]
            for w in (1, 2):
                for s in fs[w]:
                    Fser = etaq(D, s, P, cache)
                    a = revert_coeffs(tser, Fser, P)
                    if any(x.denominator != 1 for x in a): continue
                    ai = [int(x) for x in a]
                    if ai[1] == 0 and all(x == 0 for x in ai[1:6]): continue
                    key = tuple(ai[:14])
                    if key in seen:
                        seen[key].append((N, r, s, w)); continue
                    ns = fit_recurrence(ai[:26], 3)
                    if len(ns) != 1: 
                        seen[key] = [(N, r, s, w)]
                        if len(ns) == 0: continue
                    else:
                        seen[key] = [(N, r, s, w)]
                    # verify at high precision
                    tser2 = [0] + etaq(D, r, PREC_HI, cache)[:PREC_HI-1]
                    Fser2 = etaq(D, s, PREC_HI, cache)
                    a2 = revert_coeffs(tser2, Fser2, PREC_HI)
                    if any(x.denominator != 1 for x in a2): continue
                    ai2 = [int(x) for x in a2]
                    good = None
                    for v in ns:
                        den = 1
                        for x in v: den = den*x.q//math.gcd(den, x.q)
                        vv = [int(x*den) for x in v]
                        ok, pp = check_recurrence(vv, ai2)
                        if ok: good = pp; break
                    if good is None: continue
                    cl = classify(*good)
                    flip = False
                    if cl is None or cl[1] < 0:
                        af = [((-1)**i)*x for i, x in enumerate(ai2)]
                        ns2 = fit_recurrence(af[:26], 3)
                        for v in ns2:
                            den = 1
                            for x in v: den = den*x.q//math.gcd(den, x.q)
                            vv = [int(x*den) for x in v]
                            ok, pp = check_recurrence(vv, af)
                            if ok:
                                c2 = classify(*pp)
                                if c2 is not None and c2[1] > 0:
                                    cl, good, ai2, flip = c2, pp, af, True
                                break
                    results.append(dict(N=N, r=list(r), s=list(s), w=w, cusp=str(cpole), flip=flip,
                                        a=ai2[:12], cls=[str(x) for x in cl] if cl else None,
                                        poly=[[str(x) for x in p] for p in good]))
                    print("   HIT", N, r, s, "w=", w, "cls=", cl, "a=", ai2[:8], flush=True)
    print("elapsed", time.time()-t0)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hits.json'), 'w') as f:
        json.dump(results, f, indent=1)

if __name__ == '__main__':
    main()
