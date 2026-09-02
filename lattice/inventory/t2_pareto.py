"""Pareto frontier of (number of [1..2n] layers, exponent e) for the arcsine module."""
from fractions import Fraction as F
from math import comb
from qseries import lcm_upto
from dtypes import divides_type

N = 201
def pareto(c, maxlay=3, maxe=8, N=N):
    out = []
    for L in range(0, maxlay+1):
        bs = (2,)*L
        for e in range(0, maxe+1):
            ok, _ = divides_type(c, e, bs, 1, N)
            if ok:
                out.append((L, e)); break
        else:
            out.append((L, None))
    return out

print("Arcsine module F_e = sum lam^n y^n/(n^e binom(2n,n)); Pareto (layers of [1..2n], min e), n<=200")
for lam in (1, 2, 4, 8):
    print(f"-- lam_2 = {lam}")
    for e in range(0, 11):
        c = [F(0)] + [F(lam**n, n**e*comb(2*n, n)) for n in range(1, N)]
        p = pareto(c)
        print(f"   F_{e:2d}: " + "   ".join(f"{L} layer(s): e={ee}" for L, ee in p if ee is not None))
