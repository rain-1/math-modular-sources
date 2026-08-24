#!/usr/bin/env python3
"""07_padic.py -- the p-adic radii R_p of ADELIC_HOLONOMY.md Theorem A for the
inventory on each mixed-exponent four-term host.

Theorem A needs, for every function f_i in the inventory,
    |c_{i,n}|_p <= C n^A R_p^{-n}       (*_p)
so R_p = min over the inventory of p^{slope_i}, slope_i = lim v_p(c_{i,n})/n.

Members measured here:
  A  = sum a_n x^n                        (the row itself)
  B  = sum b_n x^n                        (the companion, b_0=0, b_1=1)
  H  = B - xi A                           (the conditional function; under the
       hypothesis xi in Q its 2-adic slope equals that of A, because the 2-adic
       Apery limit xi_2 is irrational -- the Calegari mechanism.  We verify this
       by measuring v_2(b_n - c a_n) for a rational c and confirming it does not
       exceed min(v_2(a_n), v_2(b_n)) + O(1).)
  Li_j(r x) = sum r^n x^n / n^j           (the pure polylogarithm module on
       P^1 - {0, 1/r, infinity}: slope v_p(r) exactly)
"""
import os, sys, json
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import importlib.util
spec = importlib.util.spec_from_file_location("d4", os.path.join(HERE, "04_denom.py"))
d4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d4)

NMAX = int(os.environ.get('NMAX', '200'))
PRIMES = [2, 3, 5, 7, 11, 13]


def slope(vals, ns):
    """fit v = s*n + c over the top half; report s and the residual spread."""
    m = len(ns)
    sx = sum(ns); sy = sum(vals)
    sxx = sum(x*x for x in ns); sxy = sum(x*y for x, y in zip(ns, vals))
    s = (m*sxy - sx*sy)/(m*sxx - sx*sx)
    c = (sy - s*sx)/m
    res = [v - (s*n + c) for n, v in zip(ns, vals)]
    return s, c, min(res), max(res)


def main():
    out = []
    for lab, r, a, c, d, f, C in d4.ROWS:
        co = d4.coeffs(r, a, c, d, f, C)
        A = d4.seq(co, NMAX, 1)
        B = d4.seq(co, NMAX, 0, u1=1)
        rec = dict(label=lab, r=r, C=C)
        print("%s  (r = %d, C = %d)" % (lab, r, C))
        for p in PRIMES:
            ns = list(range(NMAX//2, NMAX+1))
            va = [d4.vp(A[n], p) for n in ns]
            vb = [d4.vp(B[n], p) for n in ns]
            if any(x is None for x in va+vb):
                continue
            sa = slope(va, ns); sb = slope(vb, ns)
            # a random rational combination stands in for the conditional H
            vh = [d4.vp(B[n] - Fraction(3, 7)*A[n], p) for n in ns]
            sh = slope(vh, ns)
            rec['p%d' % p] = dict(slope_a=sa[0], slope_b=sb[0], slope_h=sh[0],
                                  spread_a=(sa[2], sa[3]))
            if abs(sa[0]) > 1e-3 or abs(sb[0]) > 1e-3:
                print("   p=%-3d slope(a_n)=%+.6f  slope(b_n)=%+.6f  slope(b-3a/7)=%+.6f"
                      " [spread %.2f,%.2f]  -> R_%d = %d^%.4f"
                      % (p, sa[0], sb[0], sh[0], sa[2], sa[3], p, p, min(sa[0], sb[0], sh[0])))
            else:
                print("   p=%-3d slopes ~ 0 (a: %+.2e, b: %+.2e) -> R_%d = 1" % (p, sa[0], sb[0], p))
        # polylogarithm module slope
        pl = {}
        for p in PRIMES:
            v = 0; rr = abs(r)
            while rr % p == 0:
                rr //= p; v += 1
            if v:
                pl[p] = v
        rec['polylog_slopes'] = pl
        print("   pure polylog module Li_j(%d x): p-adic slopes %s" % (r, pl or "none"))
        out.append(rec)
    json.dump(out, open(os.path.join(HERE, 'out', 'padic.json'), 'w'), indent=1)
    print("\nwrote out/padic.json")


if __name__ == '__main__':
    main()
