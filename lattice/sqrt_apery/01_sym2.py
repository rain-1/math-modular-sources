"""Exact symbolic proof-check: Sym^2(L_1) = Apery's operator, in u = t/4."""
import sympy as sp

u = sp.symbols('u')
th = sp.symbols('theta')          # for printing only

def theta_apply(expr, k):
    """apply theta = u d/du, k times"""
    for _ in range(k):
        expr = sp.expand(u*sp.diff(expr, u))
    return expr

y = sp.Function('y')

# --- Apery operator in t:  L_Ap = th^3 - t(2th+1)(17th^2+17th+5) + t^2 (th+1)^3
# rewrite in u with t = 4u  (theta_t = theta_u)
def L_Ap(f):
    t = 4*u
    T = lambda g,k: theta_apply(g,k)
    p = sp.expand((2*sp.Symbol('x')+1)*(17*sp.Symbol('x')**2+17*sp.Symbol('x')+5))
    # (2th+1)(17th^2+17th+5) = 34 th^3 + 51 th^2 + 27 th + 5
    c = sp.Poly(p, sp.Symbol('x')).all_coeffs()[::-1]   # [5,27,51,34]
    part = sum(c[k]*T(f,k) for k in range(4))
    q = sp.expand((sp.Symbol('x')+1)**3)
    d = sp.Poly(q, sp.Symbol('x')).all_coeffs()[::-1]
    part2 = sum(d[k]*T(f,k) for k in range(4))
    return sp.expand(T(f,3) - t*part + t**2*part2)

# --- L_1 = th^2 - u(136 th^2 + 68 th + 10) + 16 u^2 (th+1/2)^2
def L_1(f):
    T = lambda g,k: theta_apply(g,k)
    A = T(f,2)
    B = 136*T(f,2) + 68*T(f,1) + 10*f
    # (th+1/2)^2 = th^2 + th + 1/4
    C = T(f,2) + T(f,1) + sp.Rational(1,4)*f
    return sp.expand(A - u*B + 16*u**2*C)

# put L_1 into normal form y'' + p y' + q y
yy = sp.Function('Y')(u)
e = sp.expand(L_1(yy))
e = sp.collect(sp.expand(e), [sp.Derivative(yy,u,2), sp.Derivative(yy,u)])
c2 = sp.simplify(e.coeff(sp.Derivative(yy,u,2)))
c1 = sp.simplify(e.coeff(sp.Derivative(yy,u,1)))
c0 = sp.simplify(sp.expand(e - c2*sp.Derivative(yy,u,2) - c1*sp.Derivative(yy,u)).coeff(yy))
print("L_1 in d/du :  c2 =", sp.factor(c2), "  c1 =", sp.factor(c1), "  c0 =", sp.factor(c0))
p = sp.cancel(c1/c2); q = sp.cancel(c0/c2)
print("p =", sp.factor(p)); print("q =", sp.factor(q))

# classical symmetric square of y''+p y'+q y :
#   w''' + 3p w'' + (2p^2 + p' + 4q) w' + (4 p q + 2 q') w = 0
P3 = sp.Integer(1)
P2 = 3*p
P1 = sp.together(2*p**2 + sp.diff(p,u) + 4*q)
P0 = sp.together(4*p*q + 2*sp.diff(q,u))

# Apery operator in normal (d/du) form
e2 = sp.expand(L_Ap(yy))
d3 = sp.simplify(e2.coeff(sp.Derivative(yy,u,3)))
d2 = sp.simplify(e2.coeff(sp.Derivative(yy,u,2)))
d1 = sp.simplify(e2.coeff(sp.Derivative(yy,u,1)))
d0 = sp.simplify(sp.expand(e2 - d3*sp.Derivative(yy,u,3) - d2*sp.Derivative(yy,u,2) - d1*sp.Derivative(yy,u)).coeff(yy))
print("L_Ap in d/du: ", sp.factor(d3), sp.factor(d2), sp.factor(d1), sp.factor(d0))

print("\nDIFFERENCES (should all be 0):")
for name, X, Y in [("w''", P2, sp.cancel(d2/d3)), ("w'", P1, sp.cancel(d1/d3)), ("w", P0, sp.cancel(d0/d3))]:
    print("  ", name, sp.simplify(sp.cancel(X - Y)))
