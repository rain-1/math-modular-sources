"""Measure the sharp denominator type of a rational power series.

CDT (6.0.9): f = sum a_n y^n/(n^e prod_j [1..b_j n]),  a_n in Z.
So the test is: denom(c_n) divides n^e * prod_j [1..b_j n] for every n.
We search the minimal (e, multiset of b's) over a small grid.
"""
from fractions import Fraction as F
from qseries import lcm_upto

def denom_profile(c, N0=1, N=200):
    return [(n, c[n].denominator) for n in range(N0, min(N, len(c))) if n < len(c)]

def divides_type(c, e, bs, N0=1, N=None):
    """Does denom(c_n) | n^e prod_{b in bs} [1..b n] for all n in [N0,N)?"""
    N = len(c) if N is None else N
    for n in range(max(1, N0), N):
        d = c[n].denominator
        if d == 1: continue
        T = n**e
        for b in bs:
            T *= lcm_upto(int(b*n))
        if T % d: return False, n
    return True, None

def minimal_type(c, N0=1, N=None, max_e=6, max_layers=4, b_choices=(1, 2)):
    """Smallest total 'cost' (#layers, then e) type that works."""
    N = len(c) if N is None else N
    best = None
    for nlay in range(0, max_layers+1):
        from itertools import combinations_with_replacement as cwr
        for bs in cwr(b_choices, nlay):
            for e in range(0, max_e+1):
                ok, _ = divides_type(c, e, bs, N0, N)
                if ok:
                    cand = (sum(bs), e, bs)
                    if best is None or (cand[0], cand[1]) < (best[0], best[1]):
                        best = cand
                    break
    return best
