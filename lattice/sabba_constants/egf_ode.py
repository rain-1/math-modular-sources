import sympy as sp
t, a = sp.symbols('t a')
N=30
Y = [sp.Symbol('y0'), sp.Symbol('y1')]
for n in range(2, N+8): Y.append(sp.expand(n*Y[n-1] + a*(-1)**n*Y[n-2]))
fs = sum(Y[n]*t**n/sp.factorial(n) for n in range(N+8))
red = lambda e: sp.simplify(sp.expand(e).subs(a**6,1).subs(a**5,a).subs(a**4,1).subs(a**3,a).subs(a**2,1))
def test(expr, name, kmax=18):
    p = sp.Poly(sp.expand(expr), t)
    cs = [red(p.coeff_monomial(t**k)) for k in range(kmax)]
    print("%-52s -> %s"%(name, "ALL ZERO" if all(c==0 for c in cs) else cs[:6]))
# functional equation
test(sp.diff((1-t)*fs, t, 2) - a*fs.subs(t,-t), "((1-t)f)'' - a f(-t)")
# candidate 4th order operators
test((1-t**2)*sp.diff(fs,t,4) - (6*t+2)*sp.diff(fs,t,3) - 6*sp.diff(fs,t,2) - fs, "(1-t^2)f'''' - (6t+2)f''' - 6f'' - f")
F = (1-t)*fs
test((1-t)*sp.diff((1+t)*sp.diff(F,t,2), t, 2) - F, "(1-t)((1+t)F'')'' - F,  F=(1-t)f")
x = sp.symbols('x')
Fx = F.subs(t, 2*x-1)
test((x*(1-x)*sp.diff(Fx,x,4) + 2*(1-x)*sp.diff(Fx,x,3) - 4*Fx).subs(x,(t+1)/2), "x(1-x)F'''' + 2(1-x)F''' - 4F   [GPT-5.6 L_1]")
# Laplace dual bookkeeping: Fourier-Laplace of the f-operator
print()
z = sp.symbols('z')
Ys = sum(Y[n]*z**n for n in range(N+8))
# operator M := z^4 D^2 + 2 z^3 D + (z^4+2z-1) is the Y-operator; check the inhomogeneity again
R = -sp.Symbol('y0') + (2*sp.Symbol('y0')-sp.Symbol('y1'))*z - a*sp.Symbol('y0')*z**2 - a*(sp.Symbol('y0')-sp.Symbol('y1'))*z**3
test((z**4*sp.diff(Ys,z,2) + 2*z**3*sp.diff(Ys,z) + (z**4+2*z-1)*Ys - R).subs(z,t), "Y-operator minus RHS", 20)
