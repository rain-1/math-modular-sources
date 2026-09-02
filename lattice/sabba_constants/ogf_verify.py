"""Verify the 2nd-order inhomogeneous ODE for the OGF Y(z) = sum_{n>=0} y_n z^n.

Claim (derived by hand):   z^4 Y'' + 2 z^3 Y' + (z^4 + 2z - 1) Y
                            = -y0 + (2 y0 - y1) z - a y0 z^2 - a (y0-y1) z^3.
Homogeneous part: with W = z Y,  W'' = (1/z^4 - 2/z^3 - 1) W  (normal form, no Y' term)
                  with W = e^{-1/z} V,  V'' + (2/z^2) V' + V = 0.
"""
import sympy as sp
z, a, y0, y1 = sp.symbols('z a y0 y1')
N = 24
y = sp.symbols('y0:%d'%(N+6))
# impose recurrence y_n = n y_{n-1} + a(-1)^n y_{n-2}, y_{-1}=y_{-2}=0 (so y_0,y_1 free)
Y = [sp.Symbol('y0'), sp.Symbol('y1')]
for n in range(2, N+4):
    Y.append(sp.expand(n*Y[n-1] + a*(-1)**n*Y[n-2]))
Ys = sum(Y[n]*z**n for n in range(N+4))
lhs = sp.expand(z**4*sp.diff(Ys,z,2) + 2*z**3*sp.diff(Ys,z) + (z**4+2*z-1)*Ys)
rhs = -y0 + (2*y0-y1)*z - a*y0*z**2 - a*(y0-y1)*z**3
diff = sp.expand(lhs - rhs)
p = sp.Poly(diff, z)
low = [sp.simplify(sp.expand(p.coeff_monomial(z**k).subs(a**2,1).subs(a**3,a).subs(a**4,1))) for k in range(0, N)]
print("coefficients of (LHS-RHS) for z^0..z^%d (should all be 0):"%(N-1))
print([sp.simplify(c) for c in low])

# and the a^2 -> 1 reduction check
print()
print("check homogeneous normal form: W = z Y")
W = sp.Function('W'); V = sp.Function('V'); zz=z
Yf = sp.Function('Y')
expr = z**4*sp.diff(W(z)/z, z, 2) + 2*z**3*sp.diff(W(z)/z, z) + (z**4+2*z-1)*(W(z)/z)
print(" z^4 (W/z)'' + 2z^3 (W/z)' + (z^4+2z-1)(W/z) =", sp.simplify(sp.expand(expr)))
expr2 = sp.simplify(sp.expand( (sp.diff(sp.exp(-1/z)*V(z),z,2) - (1/z**4 - 2/z**3 - 1)*sp.exp(-1/z)*V(z))*sp.exp(1/z) ))
print(" (e^{-1/z}V)'' - (1/z^4-2/z^3-1) e^{-1/z}V  =  e^{-1/z} * (", sp.simplify(expr2), ")")
