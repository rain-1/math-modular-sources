"""CDT Remark (central binomial): test which relaxations [1..2n] -> n(n-1)C(2n,n) hold."""
from fractions import Fraction as F
from math import factorial as fac, gcd
N = 60
def L(m, c={}):
    if m in c: return c[m]
    r = 1
    for j in range(1, m+1): r = r*j//gcd(r, j)
    c[m] = r; return r
def cb(n): return n*(n-1)*(fac(2*n)//fac(n)**2)
C = [F(fac(2*k), fac(k)**2) for k in range(N+2)]
b4 = [F(0)]*(N+1)
for n in range(N): b4[n+1] = 4*sum(C[k]*C[n-k]*F(1, (2*k-1)*(2*n-2*k+1)**2) for k in range(n+1))/F(16)**n
S = {'B2': [F(0)]+[F(2*fac(n-2)*fac(n), fac(2*n)) if n >= 2 else F(0) for n in range(1, N+1)],
     'B3': [F(0)]+[F(fac(n-1)**2, fac(2*n)) for n in range(1, N+1)],
     'B4': b4,
     'B5': [F(0)]+[F(fac(n-1)**2, fac(2*n-1)*(2*n-1)) for n in range(1, N+1)],
     'B6': [F(0)]+[F(fac(n-1)**2, n*fac(2*n)) for n in range(1, N+1)],
     'B7': [F(0)]+[b4[n]/n for n in range(1, N+1)]}
def ok(k, T): return all((F(T(n))*S[k][n]).denominator == 1 for n in range(2, N+1))
V = {'B2': [('[1..2n]', lambda n: L(2*n)), ('cb', cb)],
     'B3': [('[1..2n]n', lambda n: L(2*n)*n), ('cb*n', lambda n: cb(n)*n)],
     'B4': [('[1..2n]^2', lambda n: L(2*n)**2), ('cb^2', lambda n: cb(n)**2),
            ('cb*[1..2n]', lambda n: cb(n)*L(2*n)), ('cb^2*n^2', lambda n: cb(n)**2*n*n),
            ('cb^2*[1..n]', lambda n: cb(n)**2*L(n)), ('cb^2*n^4', lambda n: cb(n)**2*n**4)],
     'B5': [('[1..2n](2n-1)', lambda n: L(2*n)*(2*n-1)), ('cb*(2n-1)', lambda n: cb(n)*(2*n-1))],
     'B6': [('[1..2n]n^2', lambda n: L(2*n)*n*n), ('cb*n^2', lambda n: cb(n)*n*n)],
     'B7': [('[1..2n]^2 n', lambda n: L(2*n)**2*n), ('cb^2*n', lambda n: cb(n)**2*n),
            ('cb*[1..2n]*n', lambda n: cb(n)*L(2*n)*n), ('cb^2*n^3', lambda n: cb(n)**3)]}
for k in ['B2','B3','B4','B5','B6','B7']:
    print(f"  {k}: " + ",  ".join(f"{nm}: {ok(k,T)}" for nm, T in V[k]))
