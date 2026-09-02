"""Sharpness diagnostics for CDT Lemma (bdenominators): is the claimed type
improvable by a constant factor, and are the individual ingredients necessary?"""
from fractions import Fraction as F
from math import factorial as fac, gcd
N = 80
def L(m, c={}):
    if m in c: return c[m]
    r = 1
    for j in range(1, m+1): r = r*j//gcd(r, j)
    c[m] = r; return r
B = {}
B['B1'] = [F(1)]+[F(0)]*N
B['B2'] = [F(0)]+[F(2*fac(n-2)*fac(n), fac(2*n)) if n >= 2 else F(0) for n in range(1, N+1)]
B['B3'] = [F(0)]+[F(fac(n-1)**2, fac(2*n)) for n in range(1, N+1)]
B['B5'] = [F(0)]+[F(fac(n-1)**2, fac(2*n-1)*(2*n-1)) for n in range(1, N+1)]
B['B6'] = [F(0)]+[F(fac(n-1)**2, n*fac(2*n)) for n in range(1, N+1)]
C = [F(fac(2*k), fac(k)**2) for k in range(N+2)]
b4 = [F(0)]*(N+1)
for n in range(N): b4[n+1] = 4*sum(C[k]*C[n-k]*F(1, (2*k-1)*(2*n-2*k+1)**2) for k in range(n+1))/F(16)**n
B['B4'] = b4
B['B7'] = [F(0)]+[b4[n]/n for n in range(1, N+1)]
CL = {'B1': ('trivial', lambda n: 1), 'B2': ('[1..2n]', lambda n: L(2*n)),
      'B3': ('[1..2n] n', lambda n: L(2*n)*n), 'B4': ('[1..2n]^2', lambda n: L(2*n)**2),
      'B5': ('[1..2n](2n-1)', lambda n: L(2*n)*(2*n-1)),
      'B6': ('[1..2n] n^2', lambda n: L(2*n)*n*n), 'B7': ('[1..2n]^2 n', lambda n: L(2*n)**2*n)}
print("For each B_i: g = gcd_{2<=n<=%d} ( T_claimed(n) / den(c_n) ).  g=1 means no constant" % N)
print("factor can be stripped from CDT's claimed type.")
for k in ['B1','B2','B3','B4','B5','B6','B7']:
    nm, T = CL[k]
    rs = []
    for n in range(2, N+1):
        if B[k][n] == 0: continue
        r = F(T(n), B[k][n].denominator)
        assert r.denominator == 1
        rs.append(r.numerator)
    g = 0
    for r in rs: g = gcd(g, r)
    print(f"  {k}: {nm:16s} g = {g}   (min ratio {min(rs) if rs else '-'} at some n)")
# necessity of each ingredient
print("\nNecessity checks (n<=%d):" % N)
def valid(k, T): return all((F(T(n))*B[k][n]).denominator == 1 for n in range(1, N+1))
print("  B3: [1..2n] alone           ->", valid('B3', lambda n: L(2*n)))
print("  B5: [1..2n] alone           ->", valid('B5', lambda n: L(2*n)))
print("  B5: [1..2n]*n^2             ->", valid('B5', lambda n: L(2*n)*n*n))
print("  B5: [1..2n]*(2n-1) (CDT)    ->", valid('B5', lambda n: L(2*n)*(2*n-1)))
print("  B6: [1..2n]*n               ->", valid('B6', lambda n: L(2*n)*n))
print("  B6: [1..n][1..2n]           ->", valid('B6', lambda n: L(n)*L(2*n)))
print("  B7: [1..2n]^2 alone         ->", valid('B7', lambda n: L(2*n)**2))
print("  B4: [1..n][1..2n]           ->", valid('B4', lambda n: L(n)*L(2*n)))
print("  B4: [1..2n]*n^2             ->", valid('B4', lambda n: L(2*n)*n*n))
print("\nCDT Remark (central binomial): [1..2n] -> n(n-1)binom(2n,n) in the claimed types")
for k in ['B2','B3','B4','B5','B6','B7']:
    nm, T = CL[k]
    def Tp(n, k=k, T=T):
        cb = n*(n-1)*(fac(2*n)//fac(n)**2) if n >= 2 else 1
        return T(n)//L(2*n)*cb if k != 'B4' and k != 'B7' else T(n)//L(2*n)**2*cb*cb
    print(f"  {k}: relaxed type valid ->", all((F(Tp(n))*B[k][n]).denominator == 1 for n in range(2, N+1)))
