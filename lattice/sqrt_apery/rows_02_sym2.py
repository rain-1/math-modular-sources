"""
rows_02_sym2.py -- EXACT symbolic proof that Sym^2(L_1) = L_parent for the eight
non-Apery sporadic rows.  (Adapted from 01_sym2.py, which does the Apery case.)

Setup.  Let the parent satisfy  (n+1)^3 A_{n+1} = q1(n) A_n - q0(n) A_{n-1};
its generating function F(t) is annihilated by
      L_par = theta^3 - t*q1(theta) + t^2*q0(theta+1),      theta = t d/dt.
Let the square root row satisfy (n+1)^2 a_{n+1} = p1(n) a_n - p0(n) a_{n-1};
since a_n = lambda^n [t^n]sqrt(F),  the series  y(u) = sum a_n u^n = sqrt(F(lambda u))
and y is annihilated by
      L_1 = theta^2 - u*p1(theta) + u^2*p0(theta+1),        theta = u d/du.
Put both in the variable u (t = lambda*u, theta_t = theta_u), write L_1 in normal
form y'' + p y' + q y = 0, and form its classical symmetric square
      w''' + 3p w'' + (2p^2 + p' + 4q) w' + (4pq + 2q') w = 0.
If that operator equals the monic-normalised L_par (in u) then y^2 solves L_par;
both L_1 and L_par have exponent 0 at u=0 with a UNIQUE analytic solution normalised
to value 1 (indicial roots (0,0) resp. (0,0,0)), so y^2 = F(lambda u) identically,
i.e. sqrt(F) is exactly the a-row.  [PROOF, not verification.]
"""
import sympy as sp

u = sp.symbols('u')
x = sp.symbols('x')
Y = sp.Function('Y')(u)

def theta_apply(expr, k):
    for _ in range(k):
        expr = sp.expand(u*sp.diff(expr, u))
    return expr

def apply_poly_theta(poly_in_x, f):
    """apply P(theta) to f, P a sympy Poly/expr in x"""
    c = sp.Poly(sp.expand(poly_in_x), x).all_coeffs()[::-1]
    return sp.expand(sum(c[k]*theta_apply(f, k) for k in range(len(c))))

def op_coeffs(expr, order):
    """extract coefficients of Y^{(order)}, ..., Y in a linear differential expression"""
    e = sp.expand(expr)
    co = []
    for k in range(order, 0, -1):
        ck = e.coeff(sp.Derivative(Y, (u, k)) if k > 1 else sp.Derivative(Y, u))
        co.append(sp.simplify(ck))
        e = sp.expand(e - ck*(sp.Derivative(Y, (u, k)) if k > 1 else sp.Derivative(Y, u)))
    co.append(sp.simplify(sp.expand(e).coeff(Y)))
    return co   # [c_order, ..., c_1, c_0]

# ---------------------------------------------------------------- row data
# (name, lambda, q1(x), q0(x)  [parent], p1(x), p0(x) [sqrt row])
def AZ(a, b, c):
    return ((2*x+1)*(a*x**2 + a*x + b), c*x**3)

ROWS = [
 ("Domb (10,4,64)",  1, AZ(10,4,64),
      20*x**2+10*x+2, 16*(2*x-1)**2),
 ("T (12,4,16)",     1, AZ(12,4,16),
      24*x**2+12*x+2, 4*(2*x-1)**2),
 ("AZ(9,3,-27)",     4, AZ(9,3,-27),
      72*x**2+36*x+6, -108*(2*x-1)**2),
 ("AZ(11,5,125)",    4, AZ(11,5,125),
      88*x**2+44*x+10, 500*(2*x-1)**2),
 ("AZ(7,3,81)",      4, AZ(7,3,81),
      56*x**2+28*x+6, 324*(2*x-1)**2),
 ("Cooper s7",       1, ((2*x+1)*(13*x**2+13*x+4), -3*x*(9*x**2-1)),
      26*x**2+13*x+2, -3*(3*x-1)*(3*x-2)),
 ("Cooper s10",      2, (2*(2*x+1)*(3*x**2+3*x+1), -4*x*(16*x**2-1)),
      24*x**2+12*x+2, -4*(8*x-3)*(8*x-5)),
 ("Cooper s18",      2, (2*(2*x+1)*(7*x**2+7*x+3), 12*x*(16*x**2-1)),
      56*x**2+28*x+6, 12*(8*x-3)*(8*x-5)),
]

allok = True
for name, lam, (q1, q0), p1, p0 in ROWS:
    t = lam*u
    # parent operator in u
    Lpar = sp.expand(theta_apply(Y,3)
                     - t*apply_poly_theta(q1, Y)
                     + t**2*apply_poly_theta(q0.subs(x, x+1), Y))
    # sqrt operator in u
    L1 = sp.expand(theta_apply(Y,2)
                   - u*apply_poly_theta(p1, Y)
                   + u**2*apply_poly_theta(p0.subs(x, x+1), Y))

    c2, c1, c0 = op_coeffs(L1, 2)
    p = sp.cancel(c1/c2); q = sp.cancel(c0/c2)
    P2 = sp.cancel(3*p)
    P1 = sp.cancel(2*p**2 + sp.diff(p,u) + 4*q)
    P0 = sp.cancel(4*p*q + 2*sp.diff(q,u))

    d3, d2, d1, d0 = op_coeffs(Lpar, 3)
    D2 = sp.cancel(d2/d3); D1 = sp.cancel(d1/d3); D0 = sp.cancel(d0/d3)

    diffs = [sp.simplify(sp.cancel(P2-D2)), sp.simplify(sp.cancel(P1-D1)),
             sp.simplify(sp.cancel(P0-D0))]
    ok = all(d == 0 for d in diffs)
    allok &= ok
    print(f"{name:18s} lambda={lam}   Sym^2(L_1) == L_parent ?  {'YES (exact, over Q(u))' if ok else 'NO'}")
    if not ok:
        print("   residuals:", diffs)
        print("   L_1 normal form: p =", sp.factor(p), "  q =", sp.factor(q))
print()
print("ALL EIGHT EXACT:" , allok)
