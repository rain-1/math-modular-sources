"""Structural checks:
 (a) which operator annihilates Y (Laplace dual bookkeeping),
 (b) Kovacic case-1 Riccati search for W'' = (1/z^4 - 2/z^3 - 1) W,
 (c) bilinear invariant coupling a=+1 and a=-1,
 (d) values f(+-1), f'(+-1)."""
import sympy as sp
z, a, y0, y1, t = sp.symbols('z a y0 y1 t')

# ---- (a) 4th-order EGF operator, verified on the series
N=26
Y = [sp.Symbol('y0'), sp.Symbol('y1')]
for n in range(2, N+6): Y.append(sp.expand(n*Y[n-1] + a*(-1)**n*Y[n-2]))
fs = sum(Y[n]*t**n/sp.factorial(n) for n in range(N+6))
Lf = sp.expand((1-t**2)*sp.diff(fs,t,4) - (2*t+6)*sp.diff(fs,t,3) + 6*sp.diff(fs,t,2) - fs)
pf = sp.Poly(Lf, t)
red = lambda e: sp.simplify(sp.expand(e).subs(a**4,1).subs(a**3,a).subs(a**2,1))
print("(a) EGF operator (1-t^2)f'''' - (2t+6)f''' + 6f'' - f  applied to the series:")
print("    coeffs of t^0..t^15:", [red(pf.coeff_monomial(t**k)) for k in range(16)])

# Laplace-dual candidate from Fourier-Laplace of L_f
Ys = sum(Y[n]*z**n for n in range(N+6))
for k in range(-2,4):
    expr = sp.expand(z**4*sp.diff(z**k*Ys,z,2) + 6*z**3*sp.diff(z**k*Ys,z) + (-z**4+6*z**3+1)*(z**k*Ys))
    pp = sp.Poly(sp.expand(expr*z**2), z)
    cs = [red(pp.coeff_monomial(z**m)) for m in range(4, 18)]
    print("    FL-operator on z^%d Y : coeffs z^4..z^17 all zero? %s"%(k, all(c==0 for c in cs)))

# ---- (b) Kovacic case 1: search rational Riccati solutions w' + w^2 = Q
print()
print("(b) Kovacic / rational Riccati search for W'' = (1/z^4 - 2/z^3 - 1) W")
Q = 1/z**4 - 2/z**3 - 1
srq = sp.series(sp.sqrt(sp.together(Q)), z, 0, 3)
print("    sqrt(Q) at 0 :", sp.simplify(sp.expand(sp.series(sp.sqrt(1-2*z-z**4), z, 0, 5).removeO()/z**2)))
# [sqrtQ]_0 = 1/z^2 - 1/z ; ([sqrtQ]_0)^2 ; b = coeff of z^-3 in Q - ([sqrtQ]_0)^2
S0 = 1/z**2 - 1/z
print("    Q - ([sqrtQ]_0)^2 =", sp.simplify(Q - S0**2), "  -> b (coeff of z^-3) = 0, v=2, a=1, alpha_0^pm = (0+2)/2 = 1")
print("    at infinity Q = -1 - 2/z^3 + 1/z^4, so [sqrtQ]_inf = i, b_inf = coeff of 1/z = 0, alpha_inf^pm = 0")
print("    d = alpha_inf - alpha_0 = -1 < 0 for every sign choice  =>  Kovacic case 1 fails.")
# brute force confirmation: w = eps/z^2 + b1/z + c + P'/P with deg P <= 6
zz = sp.symbols('zz')
found=[]
for eps in (1,-1):
    for dP in range(0,7):
        coeffs = sp.symbols('p0:%d'%(dP+1))
        P = sum(coeffs[i]*z**i for i in range(dP)) + z**dP
        b1, c = sp.symbols('b1 c')
        w = eps/z**2 + b1/z + c + sp.diff(P,z)/P
        eq = sp.together(sp.simplify(sp.diff(w,z) + w**2 - Q))
        num = sp.expand(sp.numer(sp.cancel(eq)))
        polyn = sp.Poly(num, z)
        sols = sp.solve(polyn.all_coeffs(), list(coeffs[:dP])+[b1,c], dict=True)
        if sols: found.append((eps,dP,sols))
print("    brute-force rational Riccati solutions (eps,degP,sols):", found if found else "NONE up to deg 6")
