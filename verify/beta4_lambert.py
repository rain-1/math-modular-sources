"""beta(4) Lambert-series identity (book, beta(4) section, pure CM value at tau=i/sqrt(24)).
Run: python3 verify/beta4_lambert.py"""
from mpmath import mp, mpf, exp, cosh, pi, sqrt, dirichlet
mp.dps = 40
chi = lambda n: 0 if n % 2 == 0 else (1 if n % 4 == 1 else -1)
N = 300
J = lambda y: sum(chi(n) / (mpf(n)**4 * (exp(n*y) - 1)) for n in range(1, N))
K = lambda y: sum(1 / (mpf(n)**4 * cosh(n*y)) for n in range(1, N))
y = pi / sqrt(6)
rhs = (6*J(y) - 21*J(2*y) + 14*J(3*y) - J(6*y)
       + sqrt(6) * (3*K(y)/16 - 21*K(2*y)/64 + 7*K(3*y)/48 - K(6*y)/192))
err = abs(rhs - dirichlet(4, [0, 1, 0, -1]))
print(("PASS" if err < mpf(10)**-35 else "FAIL"), "beta(4) Lambert identity, error =", err)
