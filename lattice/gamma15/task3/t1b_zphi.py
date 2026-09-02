"""T1 (part b): the actual functions of y.  B_i(y/s) = sum_n c_n s^{-n} y^n with
s = -phi^5 a unit of Z[phi].  Verify D_n * c_n * s^{-n} in Z[phi] for the CDT types."""
from fractions import Fraction as F
from math import factorial as fac, gcd
NMEAS = 80
def lcm_upto(m, c={}):
    if m in c: return c[m]
    r = 1
    for j in range(1, m+1): r = r*j//gcd(r, j)
    c[m] = r; return r
N = NMEAS
B = {}
B['B1'] = [F(1)]+[F(0)]*N
B['B2'] = [F(0)]+[F(2*fac(n-2)*fac(n), fac(2*n)) if n >= 2 else F(0) for n in range(1, N+1)]
B['B3'] = [F(0)]+[F(fac(n-1)**2, fac(2*n)) for n in range(1, N+1)]
B['B5'] = [F(0)]+[F(fac(n-1)**2, fac(2*n-1)*(2*n-1)) for n in range(1, N+1)]
B['B6'] = [F(0)]+[F(fac(n-1)**2, n*fac(2*n)) for n in range(1, N+1)]
C = [F(fac(2*k), fac(k)**2) for k in range(N+2)]
b4 = [F(0)]*(N+1)
for n in range(0, N):
    b4[n+1] = 4*sum(C[k]*C[n-k]*F(1, (2*k-1)*(2*n-2*k+1)**2) for k in range(n+1))/F(16)**n
B['B4'] = b4
B['B7'] = [F(0)]+[b4[n]/n for n in range(1, N+1)]
TY = {'B1': lambda n: 1, 'B2': lambda n: lcm_upto(2*n), 'B3': lambda n: lcm_upto(2*n)*n,
      'B4': lambda n: lcm_upto(2*n)**2, 'B5': lambda n: lcm_upto(2*n)*(2*n-1),
      'B6': lambda n: lcm_upto(2*n)*n*n, 'B7': lambda n: lcm_upto(2*n)**2*n}
# Z[phi] as pairs (a,b) <-> (a+b sqrt5)/2 with a=b mod 2;  general K element: Fractions.
def mulK(u, v):
    a, b = u; c, d = v; return (F(a*c+5*b*d, 2), F(a*d+b*c, 2))
sinv = (F(11), F(-5))              # s^{-1} = -phi^{-5} = (11-5 sqrt5)/2
s_   = (F(-11), F(-5))
assert mulK(s_, sinv) == (F(2), F(0))
def isZphi(u):
    a, b = u
    return a.denominator == 1 and b.denominator == 1 and (a-b) % 2 == 0
pw = [(F(2), F(0))]
for n in range(1, N+1): pw.append(mulK(pw[-1], sinv))
print("Z[phi]-integrality of the y-coefficients  D_n * c_n * s^{-n}  (n <= %d):" % NMEAS)
for nm in ['B1','B2','B3','B4','B5','B6','B7']:
    bad = []
    for n in range(1, NMEAS+1):
        c = B[nm][n]*TY[nm](n)
        assert c.denominator == 1
        u = (c*pw[n][0], c*pw[n][1])
        if not isZphi(u): bad.append(n)
    print(f"  {nm}: {'ALL IN Z[phi]' if not bad else 'FAILS '+str(bad[:5])}")
print("\n(and s^{-n} in Z[phi]^* for all n, so this is automatic; verified numerically anyway)")
