import sympy as sp
t, a, x, z = sp.symbols('t a x z')
N=30
Y = [sp.Symbol('y0'), sp.Symbol('y1')]
for n in range(2, N+8): Y.append(sp.expand(n*Y[n-1] + a*(-1)**n*Y[n-2]))
fs = sum(Y[n]*t**n/sp.factorial(n) for n in range(N+8))
def red(e):
    """reduce modulo a^2 = 1"""
    e = sp.expand(e)
    p = sp.Poly(e, a)
    out = 0
    for (k,), c in p.terms():
        out += c * (a if k % 2 else 1)
    return sp.simplify(sp.expand(out))
def test(expr, name, kmax=22):
    p = sp.Poly(sp.expand(expr), t)
    cs = [red(p.coeff_monomial(t**k)) for k in range(kmax)]
    bad = [(k,c) for k,c in enumerate(cs) if c != 0]
    print("%-52s -> %s"%(name, "ALL ZERO (t^0..t^%d)"%(kmax-1) if not bad else bad[:3]))
test(sp.diff((1-t)*fs, t, 2) - a*fs.subs(t,-t), "((1-t)f)'' - a f(-t)")
test((1-t**2)*sp.diff(fs,t,4) - (6*t+2)*sp.diff(fs,t,3) - 6*sp.diff(fs,t,2) - fs, "(1-t^2)f'''' - (6t+2)f''' - 6f'' - f")
F = (1-t)*fs
test((1-t)*sp.diff((1+t)*sp.diff(F,t,2), t, 2) - F, "(1-t)((1+t)F'')'' - F,  F=(1-t)f")
Fx = sp.expand(F.subs(t, 2*x-1))
test(sp.expand((x*(1-x)*sp.diff(Fx,x,4) + 2*(1-x)*sp.diff(Fx,x,3) - 4*Fx)).subs(x,(t+1)/2), "x(1-x)F''''+2(1-x)F'''-4F  [GPT-5.6 L_1]")
# also the a=-1 / L_a form with 4a^2: identical since a^2=1
Ys = sum(Y[n]*z**n for n in range(N+8))
R = -sp.Symbol('y0') + (2*sp.Symbol('y0')-sp.Symbol('y1'))*z - a*sp.Symbol('y0')*z**2 - a*(sp.Symbol('y0')-sp.Symbol('y1'))*z**3
test((z**4*sp.diff(Ys,z,2) + 2*z**3*sp.diff(Ys,z) + (z**4+2*z-1)*Ys - R).subs(z,t), "Y-operator - RHS  (OGF)", 26)
