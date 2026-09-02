"""13: growth of the hidden-variable coefficients g_k of the MIN-type
(pullback) two-variable companion, for all four rows.  If |g_k|^{1/k} stays
bounded the min-type companion is an analytic function on a polydisc and the
split-denominator two-variable count is live; if log|g_k| ~ C k log k it is a
divergent formal series and the candidate is inadmissible.
"""
from lib2v import *
from fractions import Fraction as F
from math import comb, log

N = 90


def triangular_solve(bseq, basis, N):
    g = []
    for n in range(N+1):
        s = F(0)
        for k in range(n):
            s += g[k]*basis(n, k)
        g.append((F(bseq[n]) - s)/F(basis(n, n)))
    return g


def flog(q):
    if q == 0:
        return None
    return log(abs(q.numerator)) - log(q.denominator)


CASES = []
a3, b3 = apery3(N)
CASES.append(("zeta(3) D1  basis C(n+k,2k)^2      (r=3)",
              b3, lambda n, k: comb(n+k, 2*k)**2))
a2, b2 = apery2(N)
CASES.append(("zeta(2) D1  basis C(n,k)C(n+k,2k)  (r=2)",
              b2, lambda n, k: comb(n, k)*comb(n+k, 2*k)))
aE, bE = rowE(N)
CASES.append(("row E   D1  basis C(n,k)C(2n-2k,n-k) (r=2)",
              bE, lambda n, k: comb(n, k)*comb(2*n-2*k, n-k)))
a10, b10 = cooper10(N)
CASES.append(("s_10  'D1'  basis C(n,k)^2          (r=2)",
              b10, lambda n, k: comb(n, k)**2))

for name, bs, basis in CASES:
    g = triangular_solve(bs, basis, N)
    print("=" * 78)
    print(name)
    print("=" * 78)
    print("   den(g_k) | [1..k]^r ?  minimal r for k<=40:",
          end=" ")
    rs = []
    for k in range(1, 41):
        d = g[k].denominator
        L = lcmrange(k)
        cand = [r for r in range(0, 9) if L**r % d == 0]
        rs.append(cand[0] if cand else 99)
    print(max(rs))
    print("     k    log|g_k|/k    log|g_k|/(k log k)")
    for k in [5, 10, 20, 30, 40, 50, 60, 70, 80, 90]:
        if k > N:
            continue
        L = flog(g[k])
        if L is None:
            print("   %3d    (zero)" % k)
            continue
        print("   %3d    %9.4f      %9.4f" % (k, L/k, L/(k*log(k))))
    print()
