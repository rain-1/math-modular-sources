"""01: closed forms of the separated two-variable extensions (verified by the
elementary coefficient identities, no CAS), and the triangular expansion of the
one-variable companion in the pullback basis."""
from lib2v import *
from fractions import Fraction as F
from math import comb, gcd

N = 60

print("=" * 78)
print("A. pullback closed forms  C(x,y) = R(x,y) * h(u(x,y)), verified coefficientwise")
print("=" * 78)
# [x^a] x^k/(1-x)^{2k+1} = C(a+k,2k);  [x^a] x^k/(1-x)^{k+1} = C(a,k);
# [y^b] y^k/sqrt(1-4y)   = C(2b-2k,b-k)
ok = all(z3_D1(a, b) == sum(comb(2*k, k)**2 * comb(a+k, 2*k) * comb(b+k, 2*k)
                            for k in range(min(a, b)+1))
         for a in range(26) for b in range(26))
print("  zeta3 D1 : C = 1/((1-x)(1-y)) h(u),  u = xy/((1-x)^2(1-y)^2),"
      " h = sum C(2k,k)^2 u^k  ->", ok)
ok = all(z2_D1(a, b) == sum(comb(2*k, k) * comb(a, k) * comb(b+k, 2*k)
                            for k in range(min(a, b)+1))
         for a in range(26) for b in range(26))
print("  zeta2 D1 : C = 1/((1-x)(1-y)) h(v),  v = xy/((1-x)(1-y)^2),"
      "  h = sum C(2k,k) v^k   ->", ok)
ok = all(E_D1(a, b) == sum(comb(2*k, k) * comb(a, k) * comb(2*b-2*k, b-k)
                           for k in range(min(a, b)+1))
         for a in range(26) for b in range(26))
print("  rowE  D1 : C = 1/((1-x)sqrt(1-4y)) h(w), w = xy/(1-x),"
      "       h = sum C(2k,k) w^k   ->", ok)
print("             i.e. C = 1/((1-x) sqrt(1-4y) sqrt(1-4xy/(1-x)))   -- ALGEBRAIC")

print()
print("=" * 78)
print("B. triangular expansion of the ONE-VARIABLE companion in the pullback basis")
print("=" * 78)


def triangular_solve(bseq, basis, N):
    g = []
    for n in range(N+1):
        s = F(0)
        for k in range(n):
            s += g[k] * basis(n, k)
        g.append((F(bseq[n]) - s) / F(basis(n, n)))
    return g


def den_report(g, name, kmax=30, show=14):
    print("  " + name)
    minr = {}
    for k in range(1, min(len(g), kmax+1)):
        d = g[k].denominator
        L = lcmrange(k)
        r = None
        for rr in range(0, 9):
            if L**rr % d == 0:
                r = rr
                break
        minr[k] = r
        if k <= show:
            print("    k=%2d  g_k = %-38s den=%-16d min r with den|[1..k]^r : %s"
                  % (k, str(g[k])[:38], d, r))
    vals = [v for v in minr.values() if v is not None]
    bad = [k for k, v in minr.items() if v is None]
    print("    ==> max over k<=%d of minimal r : %s   (no r<=8 at k=%s)"
          % (kmax, max(vals) if vals else None, bad))
    return minr


a3, b3 = apery3(N)
g3 = triangular_solve(b3, lambda n, k: comb(n+k, 2*k)**2, N)
den_report(g3, "zeta(3) D1:  b_n = sum_k g_k C(n+k,2k)^2   [row: g_k = C(2k,k)^2]")
print("    ratio g_k / C(2k,k)^2:")
for k in range(1, 15):
    print("       k=%2d  %s" % (k, g3[k] / F(comb(2*k, k)**2)))

a2, b2 = apery2(N)
g2 = triangular_solve(b2, lambda n, k: comb(n, k)*comb(n+k, 2*k), N)
den_report(g2, "zeta(2) D1:  b_n = sum_k g_k C(n,k)C(n+k,2k)  [row: g_k = C(2k,k)]")
print("    ratio g_k / C(2k,k):")
for k in range(1, 15):
    print("       k=%2d  %s" % (k, g2[k] / F(comb(2*k, k))))

aE, bE = rowE(N)
gE = triangular_solve(bE, lambda n, k: comb(n, k)*comb(2*n-2*k, n-k), N)
den_report(gE, "row E D1:  b_n = sum_k g_k C(n,k)C(2n-2k,n-k)  [row: g_k = C(2k,k)]")
print("    ratio g_k / C(2k,k):")
for k in range(1, 15):
    print("       k=%2d  %s" % (k, gE[k] / F(comb(2*k, k))))
