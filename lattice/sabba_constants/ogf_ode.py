"""Eliminate Y(-z) from the OGF functional-differential equation, keeping inhomogeneous terms."""
import sympy as sp

z, a, y0, y1 = sp.symbols('z a y0 y1')
Y = sp.Function('Y'); V = sp.Function('V')

# Eq1:  z^2 Y'(z) + (z-1) Y(z) + a z^2 V(z) = p(z),  V(z)=Y(-z), p = -y0 + (y0-y1) z
p  = -y0 + (y0-y1)*z
pm = p.subs(z, -z)
E1 = z**2*sp.diff(Y(z),z) + (z-1)*Y(z) + a*z**2*V(z) - p
# Eq2:  a z^2 Y(z) = z^2 V'(z) + (z+1) V(z) + p(-z)
E2 = a*z**2*Y(z) - z**2*sp.diff(V(z),z) - (z+1)*V(z) - pm

Vexpr = sp.solve(E1, V(z))[0]
E2s = E2.subs(V(z), Vexpr).doit()
E2s = E2s.subs(sp.Derivative(Vexpr, z), sp.diff(Vexpr, z))
E2s = sp.simplify(sp.expand(E2s.replace(sp.Derivative(V(z), z), sp.diff(Vexpr, z))))
# safer: build explicitly
Vex = Vexpr
E2b = a*z**2*Y(z) - z**2*sp.diff(Vex, z) - (z+1)*Vex - pm
E2b = sp.simplify(sp.expand(E2b*a*z**2))   # clear the 1/(a z^2)
E2b = sp.expand(E2b.subs(a**2, 1))
# substitute a^2 -> 1 thoroughly
E2b = sp.simplify(sp.Poly(sp.expand(E2b), a).as_expr())
print("raw:", sp.simplify(E2b))
