"""02: the hidden-variable companion series h_comp(u) = sum g_k u^k.
Growth, holonomy, inhomogeneity, and the polydisc of convergence of
D(x,y) = R(x,y) h_comp(u(x,y)).
"""
from lib2v import *
from fractions import Fraction as F
from math import comb, gcd, log
import sympy as sp

N = 80


def triangular_solve(bseq, basis, N):
    g = []
    for n in range(N+1):
        s = F(0)
        for k in range(n):
            s += g[k] * basis(n, k)
        g.append((F(bseq[n]) - s) / F(basis(n, n)))
    return g


def flog(q):
    """log |q| for a Fraction"""
    if q == 0:
        return None
    n, d = abs(q.numerator), q.denominator
    return len(str(n)) * log(10) - len(str(d)) * log(10) + \
        log(float(str(n)[:15]) / 10**14) - log(float(str(d)[:15]) / 10**14)


a3, b3 = apery3(N)
g3 = triangular_solve(b3, lambda n, k: comb(n+k, 2*k)**2, N)

print("=" * 78)
print("A. growth of g_k   (zeta(3), basis C(n+k,2k)^2)")
print("=" * 78)
print("   k   log|g_k|/k    log|g_k/g_{k-1}|      log|g_k|/(k log k)")
for k in list(range(2, 20)) + list(range(20, N+1, 10)):
    L = flog(g3[k])
    Lm = flog(g3[k-1])
    print("  %3d   %9.4f      %10.4f            %8.4f"
          % (k, L/k, L-Lm, L/(k*log(k))))
print("   [ log 16 = %.4f, log 34 = %.4f ]" % (log(16), log(34)))

print()
print("=" * 78)
print("B. holonomy test for g_k: linear recurrence with polynomial coefficients?")
print("=" * 78)
# search for sum_{i=0..I} p_i(k) g_{k+i} = 0, deg p_i <= d, over Q
for I in range(1, 5):
    for d in range(0, 6):
        nunk = (I+1)*(d+1)
        rows = []
        for k in range(1, min(nunk + 12, N - I)):
            rows.append([F(k)**j * g3[k+i] for i in range(I+1) for j in range(d+1)])
        if len(rows) < nunk + 4:
            continue
        M = sp.Matrix([[sp.Rational(v.numerator, v.denominator) for v in r] for r in rows])
        ns = M.nullspace()
        if ns:
            print("   FOUND recurrence: order %d, degree %d, nullspace dim %d" % (I, d, len(ns)))
            break
    else:
        continue
    break
else:
    print("   no recurrence of order<=4, degree<=5 found  ==> g_k is (very likely)"
          " NOT holonomic")

print()
print("=" * 78)
print("C. the inhomogeneity: apply L_u = theta^2 - 4u(2 theta+1)^2 to h_comp")
print("   (L_u annihilates h(u) = sum C(2k,k)^2 u^k)")
print("=" * 78)
eps = [F(0)]
for k in range(1, N+1):
    eps.append(F(k*k)*g3[k] - F(4*(2*k-1)**2)*g3[k-1])
print("   eps_k for k=1..12:")
for k in range(1, 13):
    print("     k=%2d  %s" % (k, eps[k]))
print("   log|eps_k|/k :", ["%.3f" % (flog(eps[k])/k) for k in range(5, 30, 3)])
