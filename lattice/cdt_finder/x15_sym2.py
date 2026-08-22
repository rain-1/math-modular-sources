"""The X_1(5) Sym^2 host (Zagier's D parametrisation) -- the Beukers 1987 Thm 4
configuration -- as a candidate for a CDT-style amplification.

Weight-1 layer  = Zagier D row (a,b,c) = (11,3,-1):
      (n+1)^2 a_{n+1} = (11n^2+11n+3) a_n + n^2 a_{n-1},  a_0=1,a_1=3.
Char. roots  lam = (11 +- 5 sqrt5)/2 = phi^5, -phi^{-5};  singular t-values
      t_1 = 1/phi^5 = phi^{-5} = 0.0901699...,  t_2 = -phi^5 = -11.0901699...,
      t_1 t_2 = -1   (Apery-perfect fold; N_{Q(sqrt5)/Q}(t_2) = -1, a unit).

Sym^2 layer: A_n := [t^n] F(t)^2 with F(t) = sum a_n t^n.  We fit the minimal
recurrence, extract its characteristic roots, build the census companion
(b_0=0,b_1=1) and measure the sharp denominator exponent k with d_n^k B_n in Z.
"""
from fractions import Fraction
from math import log, sqrt, gcd
import sympy as sp

N = 260

# ---- weight-1 row -----------------------------------------------------------
a = [Fraction(1), Fraction(3)]
for n in range(1, N):
    a.append((Fraction(11*n*n + 11*n + 3)*a[n] + Fraction(n*n)*a[n-1]) / Fraction((n+1)**2))
assert all(x.denominator == 1 for x in a[:60]), 'weight-1 row not integral'
print('Zagier D a_n :', [int(x) for x in a[:8]])
phi = (1+sqrt(5))/2
print(f'  char roots  {(11+5*sqrt(5))/2:.10f}, {(11-5*sqrt(5))/2:.10f}   phi^5={phi**5:.10f}')
print(f'  singular t  t1={2/(11+5*sqrt(5)):.10f}  t2={-(11+5*sqrt(5))/2:.10f}   t1*t2=-1')

# ---- Sym^2 row: A_n = [t^n] F^2 ---------------------------------------------
A = [sum(a[i]*a[n-i] for i in range(n+1)) for n in range(N)]
assert all(x.denominator == 1 for x in A[:80]), 'Sym^2 row not integral'
print('Sym^2 A_n    :', [int(x) for x in A[:8]])
print('  (n+1) | A_n ?', [int(A[n]) % (n+1) == 0 for n in range(8)],
      ' -> free integration' if all(int(A[n]) % (n+1) == 0 for n in range(40)) else ' -> NO free integration')

# ---- fit the minimal recurrence  sum_{j=0}^{J} P_j(n) A_{n+j} = 0 -----------
def fit(seq, J, deg):
    n = sp.symbols('n')
    unk = sp.symbols(f'c0:{(J+1)*(deg+1)}')
    P = [sum(unk[j*(deg+1)+d]*n**d for d in range(deg+1)) for j in range(J+1)]
    eqs = []
    for nn in range(0, min(len(seq)-J-1, 6*(J+1)*(deg+1))):
        eqs.append(sum(P[j].subs(n, nn)*sp.Rational(seq[nn+j]) for j in range(J+1)))
    sol = sp.linsolve(eqs, unk)
    sol = list(sol)[0]
    free = [s for s in sol if s.free_symbols]
    if all(s == 0 for s in sol):
        return None
    return sol, P, unk

rec = None
for J in (2, 3):
    for deg in (2, 3, 4):
        r = fit([int(x) for x in A], J, deg)
        if r is None:
            continue
        sol, P, unk = r
        fs = sorted(set().union(*[s.free_symbols for s in sol if s.free_symbols]) or set(), key=str)
        if not fs:
            continue
        sub = {f: (1 if i == 0 else 0) for i, f in enumerate(fs)}
        coeffs = [sp.simplify(sp.expand(P[j].subs(dict(zip(unk, sol))).subs(sub))) for j in range(J+1)]
        if all(c == 0 for c in coeffs):
            continue
        # verify
        nn = sp.symbols('n')
        ok = all(sum(coeffs[j].subs(nn, t)*sp.Rational(int(A[t+j])) for j in range(J+1)) == 0
                 for t in range(0, 150))
        if ok:
            rec = (J, deg, coeffs)
            break
    if rec:
        break

J, deg, coeffs = rec
nn = sp.symbols('n')
print(f'\nSym^2 minimal recurrence: order {J}, poly degree {deg}, verified n<=150')
for j, c in enumerate(coeffs):
    print(f'   P_{j}(n) = {sp.factor(c)}')
lead = [sp.LT(sp.Poly(c, nn)) if c != 0 else 0 for c in coeffs]
charpoly = sum(sp.Poly(coeffs[j], nn).LC()*sp.symbols('L')**j for j in range(J+1))
roots = sp.nroots(sp.Poly(charpoly, sp.symbols('L')))
print('   characteristic roots:', [complex(r) for r in roots])
print('   -> singular t-values 1/root:', [complex(1/r) for r in roots if abs(r) > 1e-12])

# ---- companion B_n (b_0=0, b_1=1) and sharp denominator exponent ------------
B = [Fraction(0), Fraction(1)] + [Fraction(0)]*(N-2)
Pn = [sp.lambdify(nn, c, 'math') for c in coeffs]
for t in range(0, N-J-1):
    lc = coeffs[J].subs(nn, t)
    if lc == 0:
        break
    s = sum(sp.Rational(coeffs[j].subs(nn, t))*B[t+j] for j in range(J))
    if t+J < N:
        B[t+J] = Fraction(-sp.Rational(s)/sp.Rational(lc))
lim = float(B[200]/A[200]) if A[200] else None
print(f'\ncompanion limit B_n/A_n at n=200 : {lim!r}')

def dn(n):
    from math import gcd
    L = 1
    for i in range(1, n+1):
        L = L*i//gcd(L, i)
    return L

worst = 0
for n in range(2, 121):
    d = dn(n)
    for k in range(0, 8):
        if (Fraction(d)**k*B[n]).denominator == 1:
            worst = max(worst, k)
            break
print(f'sharp denominator exponent k (d_n^k B_n in Z, n<=120): k = {worst}')
