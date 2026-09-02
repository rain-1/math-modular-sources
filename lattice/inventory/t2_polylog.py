"""Task 2: y-expansions and sharp denominator types of the symmetrised
polylogarithms on the descent orbifold, and of CDT's factorial series."""
from fractions import Fraction as F
from math import comb, factorial as fac
from qseries import mul, compose, w_series, to_y
from dtypes import divides_type

def pareto(c, maxlay=4, maxe=10, N=None):
    N = len(c) if N is None else N
    out = []
    for L in range(0, maxlay+1):
        bs = (2,)*L
        for e in range(0, maxe+1):
            ok, _ = divides_type(c, e, bs, 1, N)
            if ok:
                out.append((L, e)); break
    return out

def show(name, c, N=None):
    p = pareto(c, N=N)
    txt = "   ".join(f"[1..2n]^{L}: e={e}" for L, e in p) or "no type found"
    lead = " ".join(str(v) for v in c[1:5])
    print(f"  {name:26s} {txt}")
    return p

Px, Py = 170, 80
for s in (F(1), F(1, 4)):
    lam = 1/s
    print(f"\n=== descent host s = {s} (lambda_2 = {lam}), branch y = {4*s}, y = x^2/(x-s) ===")
    W = w_series(s, Px)
    # Li_j(x/s) = sum lam^n x^n/n^j    (integral coefficients over the lcm-free type n^j in x)
    for j in (1, 2, 3, 4):
        Li = [F(0)]+[F(int(lam)**n, n**j) for n in range(1, Px)]
        LiW = compose(Li, W, Px)
        sym = [a+b for a, b in zip(Li, LiW)]
        asym = [a-b for a, b in zip(Li, LiW)]
        xmw = [a-b for a, b in zip([F(0), F(1)]+[F(0)]*(Px-2), W)]   # x - w(x)
        symy = to_y(sym, s, Py)
        prod = to_y(mul(xmw, asym, Px), s, Py)
        show(f"Sym+ Li_{j}", symy)
        show(f"(x-w)*Sym- Li_{j}", prod)

print("\n=== CDT's factorial series B_2,B_3,B_5,B_6 (their (10.1.2),(10.1.4),(10.2.1)) ===")
N = 160
B2 = [F(0), F(0)]+[F(2*fac(n-2)*fac(n), fac(2*n)) for n in range(2, N)]
B3 = [F(0)]+[F(fac(n-1)**2, fac(2*n)) for n in range(1, N)]
B5 = [F(0)]+[F(fac(n-1)**2, fac(2*n-1)*(2*n-1)) for n in range(1, N)]
B6 = [F(0)]+[F(fac(n-1)**2, n*fac(2*n)) for n in range(1, N)]
for nm, c in (('B2 = 2(n-2)!n!/(2n)!', B2), ('B3 = ((n-1)!)^2/(2n)!  [= F_2]', B3),
              ('B5 = ((n-1)!)^2/((2n-1)!(2n-1))', B5), ('B6 = ((n-1)!)^2/(n(2n)!) [= F_3]', B6)):
    show(nm, c)
print("\n  identities:  B3_n = 1/(n^2 C(2n,n)) = F_2 ;  B6_n = 1/(n^3 C(2n,n)) = F_3 ;"
      "  B2_n = 2/(n(n-1)C(2n,n))")
print("  check B3==F_2:", all(B3[n] == F(1, n**2*comb(2*n, n)) for n in range(1, N)))
print("  check B6==F_3:", all(B6[n] == F(1, n**3*comb(2*n, n)) for n in range(1, N)))
