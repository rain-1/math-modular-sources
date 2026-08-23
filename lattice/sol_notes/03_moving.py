"""SOL note 1 sec 5 / note 2 sec 8-9: the naive moving EMN period.
   Phi(z) = int_Delta dxdy/(1-z(x^2+y^2)),  H = z*Phi.
Checks: moment formula, the 3-term recurrence, Phi(1)=G, the inhomogeneous
Picard-Fuchs equation, H'(z) closed form, the log branch at z=1, and the
"2*int_0^q artanh(u)/(1+u^2) du" identity."""
from mpmath import mp, mpf, binomial, catalan, quad, atanh, asin, sqrt, log, pi, nsum, inf
import sympy as sp

mp.dps = 60
print("== A1. moment formula a_n = int_Delta (x^2+y^2)^n  vs  2^{-(n+1)}/(n+1) * sum C(n,k)/(2k+1)")
z = sp.symbols('z')
def exact_moment(n):
    # int_Delta x^a y^b = a! b!/(a+b+2)!
    return sum(sp.binomial(n,j)*sp.factorial(2*j)*sp.factorial(2*n-2*j)/sp.factorial(2*n+2) for j in range(n+1))
for n in range(0,9):
    exact = sp.nsimplify(exact_moment(n))
    claim = sp.Rational(1,2**(n+1)*(n+1))*sum(sp.binomial(n,k)/sp.Integer(2*k+1) for k in range(n+1))
    print(f"  n={n}: exact={exact}  claim={claim}  equal={sp.simplify(exact-claim)==0}")

print("\n== A2. recurrence (n+3)(2n+5)a_{n+2} - (n+2)(3n+5)a_{n+1} + (n+1)^2 a_n = 0")
def a(n):
    return sp.Rational(1,2**(n+1)*(n+1))*sum(sp.binomial(n,k)/sp.Integer(2*k+1) for k in range(n+1))
A=[a(n) for n in range(0,30)]
res=[sp.simplify((n+3)*(2*n+5)*A[n+2]-(n+2)*(3*n+5)*A[n+1]+(n+1)**2*A[n]) for n in range(0,28)]
print("  residuals n=0..27 all zero:", all(r==0 for r in res), " first few:", res[:4])

print("\n== A3. Phi(1) = sum a_n  vs  G")
mp.dps=40
def af(n):
    return mpf(1)/(mpf(2)**(n+1)*(n+1))*sum(binomial(n,k)/mpf(2*k+1) for k in range(n+1))
Ssum = nsum(lambda n: af(int(n)), [0, inf])
print("  sum a_n =", Ssum)
print("  G       =", catalan)
print("  diff    =", Ssum-catalan)

print("\n== A4. inhomogeneous PF:  z(z-2)H'' + (z-1)H' = 1/(2(z-1)),  H=z*Phi")
N=22
Hs = sum(A[n]*z**(n+1) for n in range(N))
lhs = sp.expand(z*(z-2)*sp.diff(Hs,z,2) + (z-1)*sp.diff(Hs,z))
rhs = sp.series(sp.Rational(1,2)/(z-1), z, 0, N).removeO()
d = sp.expand(lhs-rhs)
print("  lhs-rhs series (should be O(z^%d)):"%(N-2), sp.Poly(d,z).all_coeffs()[-(N-2):][::-1][:8], "...")
print("  lowest surviving order:", min([m for m in range(0,N) if sp.expand(d).coeff(z,m)!=0], default=None))

print("\n== A5. H'(z) = artanh(sqrt(z/(2-z)))/sqrt(z(2-z))")
mp.dps=40
def Hp_series(zz, N=400):
    return sum(af(n)*(n+1)*zz**n for n in range(N))
for zz in [mpf('0.3'), mpf('0.7'), mpf('0.95')]:
    lhs = Hp_series(zz)
    rhs = atanh(sqrt(zz/(2-zz)))/sqrt(zz*(2-zz))
    print(f"  z={zz}:  series={lhs}  closed={rhs}  diff={lhs-rhs}")

print("\n== A6. G = 2 int_0^1 artanh(u)/(1+u^2) du ; and H(z(q)) = 2 int_0^q ...")
val = 2*quad(lambda t: atanh(t)/(1+t**2), [0,1])
print("  2*int =", val, "   G =", catalan, "   diff =", val-catalan)
for qq in [mpf('0.4'), mpf('0.8')]:
    zz = 2*qq**2/(1+qq**2)
    lhs = zz*sum(af(n)*zz**n for n in range(400))
    rhs = 2*quad(lambda t: atanh(t)/(1+t**2), [0,qq])
    print(f"  q={qq} (z={zz}): H={lhs}  2*int={rhs}  diff={lhs-rhs}")

print("\n== A7. branch at z=1:  H(z) = G + ((1-z)/2) log(1-z) + O(1-z)?")
mp.dps=40
def Hval(zz,N=4000):
    return zz*sum(af(n)*zz**n for n in range(N))
for eps in [mpf('1e-2'), mpf('1e-3'), mpf('1e-4')]:
    zz = 1-eps
    lhs = Hval(zz)
    pred = catalan + (eps/2)*log(eps)
    print(f"  1-z={eps}: H={lhs}  G+(eps/2)log eps={pred}   (H-G)/(eps*log eps)={(lhs-catalan)/(eps*log(eps))}")

print("\n== A8. homogeneous solutions 1 and 2*arcsin(sqrt(z/2))")
zz=mpf('0.6')
f = lambda t: 2*asin(sqrt(t/2))
h2 = mp.diff(f,zz,2); h1=mp.diff(f,zz,1)
print("  z(z-2)f''+(z-1)f' at z=0.6 :", zz*(zz-2)*h2+(zz-1)*h1)
