#!/usr/bin/env python3
"""Enumerate eta-quotient parameters t on Gamma_0(N) with

    weight 0,  ord_inf(t) = 1,  divisor supported on cusps,  deg(t) <= DEGMAX.

deg(t)=1  <=> t is a Hauptmodul of Gamma_0(N)  (genus 0)
deg(t)=2  <=> t is (the pullback of) a Hauptmodul of an Atkin-Lehner quotient
              Gamma_0(N)+W with |W|=2 (this is how Apery's and Domb's parameters arise)
deg(t)=4  <=> quotient by a Klein four-group of Atkin-Lehner involutions.

Enumeration is over the *divisor* (cusp-order vector) rather than the exponent
vector r, which makes it exhaustive and tiny; r is recovered by inverting the
Ligozat matrix and kept only when integral.
"""
import sys, os, json, itertools
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qser import divisors, phi, ligozat_matrix, cusp_counts, etaquo_series
from math import gcd

DEGMAX = int(os.environ.get("DEGMAX", "4"))
NMAX   = int(os.environ.get("NMAX", "120"))
PREC   = int(os.environ.get("TPREC", "80"))

def matinv(M):
    n = len(M)
    A = [row[:] + [Fraction(1 if i == j else 0) for j in range(n)] for i, row in enumerate(M)]
    for col in range(n):
        piv = next((i for i in range(col, n) if A[i][col] != 0), None)
        if piv is None: return None
        A[col], A[piv] = A[piv], A[col]
        pv = A[col][col]
        A[col] = [x / pv for x in A[col]]
        for i in range(n):
            if i != col and A[i][col] != 0:
                f = A[i][col]
                A[i] = [x - f*y for x, y in zip(A[i], A[col])]
    return [row[n:] for row in A]

def parts_with_weight(weights, total):
    """all dicts idx->mult>=1 with sum weights[idx]*mult == total."""
    out = []
    n = len(weights)
    def rec(i, rem, cur):
        if rem == 0:
            out.append(dict(cur)); return
        if i == n: return
        rec(i+1, rem, cur)
        w = weights[i]
        m = 1
        while w*m <= rem:
            cur[i] = m
            rec(i+1, rem - w*m, cur)
            del cur[i]
            m += 1
    rec(0, total, {})
    return out

def gen_level(N):
    D, M = ligozat_matrix(N)
    Minv = matinv(M)
    if Minv is None: return []
    tau = len(D)
    cc = cusp_counts(N)
    w = [cc[c] for c in D]            # multiplicity of each cusp class
    iinf = D.index(N)                 # cusp infinity  (denominator N)
    # columns of Minv scaled to integers: r = sum_c v_c * col_c
    L = 1
    for i in range(tau):
        for j in range(tau):
            L = L*Minv[i][j].denominator // gcd(L, Minv[i][j].denominator)
    cols = [[int(Minv[i][c]*L) for i in range(tau)] for c in range(tau)]  # cols[c][i]
    others = [c for c in range(tau) if c != iinf]
    wo = [w[c] for c in others]

    res = []
    for deg in range(1, DEGMAX+1):
        extra = deg - w[iinf]
        if extra < 0: continue
        posparts = parts_with_weight(wo, extra) if extra > 0 else [{}]
        negparts = parts_with_weight(wo, deg)
        # residues mod L
        table = {}
        for pp in posparts:
            vec = [0]*tau
            for i in range(tau):
                vec[i] = cols[iinf][i]
                for k, m in pp.items():
                    vec[i] += m*cols[others[k]][i]
            table.setdefault(tuple(x % L for x in vec), []).append(pp)
        for np_ in negparts:
            vec = [0]*tau
            for i in range(tau):
                s = 0
                for k, m in np_.items():
                    s -= m*cols[others[k]][i]
                vec[i] = s
            key = tuple((-x) % L for x in vec)
            for pp in table.get(key, []):
                if set(pp) & set(np_): continue
                r = [0]*tau
                for i in range(tau):
                    s = cols[iinf][i]
                    for k, m in pp.items(): s += m*cols[others[k]][i]
                    for k, m in np_.items(): s -= m*cols[others[k]][i]
                    assert s % L == 0
                    r[i] = s // L
                if sum(r) != 0: continue
                ordv = {D[iinf]: 1}
                for k, m in pp.items(): ordv[D[others[k]]] = m
                for k, m in np_.items(): ordv[D[others[k]]] = -m
                res.append((r, deg, ordv))
    return D, res

def main():
    out = []
    seen = {}
    for N in range(1, NMAX+1):
        g = gen_level(N)
        if not g: continue
        D, cands = g
        for r, deg, ordv in cands:
            o24, ser = etaquo_series(D, r, PREC)
            assert o24 == 24, (N, r, o24)
            key = tuple(ser[:40])
            if key in seen: continue
            seen[key] = (N, r)
            out.append({"N": N, "D": D, "r": r, "deg": deg,
                        "ord": {str(k): v for k, v in ordv.items()},
                        "ser": ser[:PREC]})
        print(f"N={N}: {len(cands)} raw", flush=True)
    print("total distinct t:", len(out))
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "t_list.json"), "w") as f:
        json.dump(out, f)

if __name__ == "__main__":
    main()
