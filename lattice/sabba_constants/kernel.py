"""Moment kernel: y_n = int_Gamma t^n w(t) dt solves y_n = n y_{n-1} + a(-1)^n y_{n-2}
   iff  w' + w = a w(-t)/t^2.  With w = A + B (A even, B odd):
        B' = A(a/t^2 - 1),  A' = -B(1 + a/t^2),
   whence A'' + 2A'/(t(1+a t^2)) - (1 - t^-4) A = 0.
   Question: are the regular singular points t^2 = -1/a apparent (no log)?"""
import sympy as sp
t, a = sp.symbols('t a')
A = sp.Function('A')
ode = sp.diff(A(t),t,2) + 2*sp.diff(A(t),t)/(t*(1+a*t**2)) - (1 - t**(-4))*A(t)
# work at a = 1, t0 = I  (t0^2 = -1/a)
for aval, t0 in ((1, sp.I), (-1, sp.Integer(1))):
    u = sp.symbols('u')
    P = sp.simplify((2/(t*(1+aval*t**2))).subs(t, t0+u))
    Qc = sp.simplify((-(1 - t**(-4))).subs(t, t0+u))
    # Frobenius: A = sum c_k u^(k+rho), exponents 0 and 2 -> test the k=2 obstruction for rho=0
    c = sp.symbols('c0:6')
    ser = sum(c[k]*u**k for k in range(6))
    expr = sp.diff(ser,u,2) + P*sp.diff(ser,u) + Qc*ser
    expr = sp.series(sp.expand(sp.simplify(expr)), u, 0, 2).removeO()
    p = sp.Poly(sp.expand(expr), u)
    e0 = sp.simplify(p.coeff_monomial(u**-1)) if p.coeff_monomial(u**-1) else 0
    lead = sp.simplify(sp.limit(sp.expand(expr)*u, u, 0))
    print("a = %s, t0 = %s"%(aval, t0))
    print("   coefficient of u^-1 in L[sum c_k u^k] (obstruction at rho=0, k=2):")
    print("   ", sp.simplify(sp.expand(lead)))
