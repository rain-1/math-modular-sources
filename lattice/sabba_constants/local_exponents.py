import sympy as sp
x, X, a, r = sp.symbols('x X a rho')
# L_a[y] = x(1-x)y'''' + 2(1-x)y''' - 4a^2 y ;  a^2 = 1
def indicial(sub, var):
    y = var**r
    e = sp.expand(sub(y, var))
    e = sp.simplify(e/var**(r-3))
    return sp.simplify(sp.limit(e, var, 0))
Lx = lambda y,v: v*(1-v)*sp.diff(y,v,4) + 2*(1-v)*sp.diff(y,v,3) - 4*y
print("indicial at x=0 :", sp.factor(indicial(Lx, x)), " -> exponents", sp.solve(indicial(Lx,x), r))
# at x=1 : X=1-x
LX = lambda y,v: (1-v)*v*sp.diff(y,v,4) - 2*v*sp.diff(y,v,3) - 4*y
print("indicial at x=1 :", sp.factor(indicial(LX, X)), " -> exponents", sp.solve(indicial(LX,X), r))
# obstruction at X^3 for the exponent-0 solution
c = sp.symbols('c0:8'); Xs = sp.Symbol('X')
ser = sum(c[k]*Xs**k for k in range(8))
res = sp.expand(LX(ser, Xs))
p = sp.Poly(res, Xs)
print()
print("Frobenius at x=1, rho=0, coefficients of X^-3..X^1 of L[sum c_k X^k]:")
for k in range(-3, 2):
    print("   X^%2d : %s"%(k, sp.simplify(p.coeff_monomial(Xs**k) if k>=0 else 0)))
print("   (lowest surviving order is X^0 :", sp.simplify(p.coeff_monomial(Xs**0)), ") -> forces c0 = 0 unless a log is added")
print()
# same at x=0 (repeated exponent 1)
d = sp.symbols('d0:8'); xs = sp.Symbol('x')
ser0 = sum(d[k]*xs**k for k in range(8))
p0 = sp.Poly(sp.expand(Lx(ser0, xs)), xs)
print("Frobenius at x=0, rho=0, coefficients of x^0..x^2:")
for k in range(0,3):
    print("   x^%d : %s"%(k, sp.simplify(p0.coeff_monomial(xs**k))))
