#!/usr/bin/env python3
"""01_general_move.py -- symbolic derivation of the cusp move for a GENERAL
second-order Apery-like row

      (n+1)^2 u_{n+1} = P(n) u_n - Q(n) u_{n-1},   deg P, deg Q <= 2,
      P(n) = p2 n^2 + p1 n + p0,  Q(n) = q2 n^2 + q1 n + q0,
      p2 = lam + mu,  q2 = lam*mu  (lam, mu = characteristic roots).

Operator:  L = th^2 - t P(th) + t^2 Q(th+1),  th = t d/dt,  i.e.
           L y = R2 t^2 y'' + (R2+R1) t y' + R0 y.
Move:      s = t/(1-lam t),  z(s) = (1-lam t)^al y(t).

Output: the exact operator identity  L# z = (1-lam t)^(al-1) L y,
        the transformed row (P#, Q#), and the constraint on al.

Claude (Opus 5), 2026-08-22.
"""
import sympy as sp

t, s, lam, mu, al, n = sp.symbols('t s lam mu al n')
p1, p0, q1, q0 = sp.symbols('p1 p0 q1 q0')
p2 = lam + mu
q2 = lam*mu

P = p2*n**2 + p1*n + p0
Q = q2*n**2 + q1*n + q0

Qsh = sp.expand(Q.subs(n, n+1))
R2 = 1 - p2*t + q2*t**2
R1 = -sp.Poly(P, n).coeff_monomial(n)*t + sp.Poly(Qsh, n).coeff_monomial(n)*t**2
R0 = -sp.Poly(P, n).coeff_monomial(1)*t + sp.Poly(Qsh, n).coeff_monomial(1)*t**2

z0s, z1s, z2s = sp.symbols('z0 z1 z2')
u  = t/(1-lam*t)
up, upp = sp.diff(u, t), sp.diff(u, t, 2)
g   = (1-lam*t)**(-al)
gp, gpp = sp.diff(g, t), sp.diff(g, t, 2)
Z0, Z1, Z2 = z0s, z1s*up, z2s*up**2 + z1s*upp
yv, yp, ypp = g*Z0, gp*Z0 + g*Z1, gpp*Z0 + 2*gp*Z1 + g*Z2

expr = sp.expand(R2*t**2*ypp + (R2+R1)*t*yp + R0*yv)
expr = sp.powdenest(sp.simplify(expr.subs(t, s/(1+lam*s))), force=True)
pol  = sp.Poly(sp.expand(expr), z0s, z1s, z2s)
red  = lambda e: sp.cancel(sp.powsimp(sp.expand(sp.powdenest(sp.simplify(e), force=True))))
c2, c1, c0 = (red(pol.coeff_monomial(v)) for v in (z2s, z1s, z0s))

rho = (1+lam*s)**(al-1)
n2 = sp.cancel(sp.simplify(sp.powdenest(c2/(rho*s**2), force=True)))
n1 = sp.cancel(sp.simplify(sp.powdenest(c1/(rho*s), force=True)))
n0 = sp.cancel(sp.simplify(sp.powdenest(c0/rho, force=True)))

print("=== L# z = rho^{-1} L y  with rho = (1+lam s)^(al-1) = (1-lam t)^(1-al) ===")
print("R2#     =", sp.factor(sp.expand(n2)))
r1, r2 = sp.symbols('r1 r2')
# Q(n) = lam*mu*(n-r1)(n-r2)  =>  q1 = -lam*mu*(r1+r2), q0 = lam*mu*r1*r2
alsub = {q1: -lam*mu*(r1+r2), q0: lam*mu*r1*r2, al: 1-r1}

def polyify(e):
    e = sp.cancel(sp.together(sp.expand(e)))
    num, den = sp.fraction(e)
    qt, rm = sp.div(sp.expand(num), sp.expand(den), s)
    if sp.simplify(rm) != 0:
        print("   [residue before imposing the alpha-condition:]", sp.factor(sp.simplify(rm)))
        e2 = sp.cancel(sp.together(sp.expand(e.subs(alsub))))
        num, den = sp.fraction(e2)
        qt, rm = sp.div(sp.expand(num), sp.expand(den), s)
        assert sp.simplify(rm) == 0, ("still non-polynomial", sp.simplify(rm))
    return sp.expand(qt)
R2n, R1n, R0n = polyify(n2), polyify(n1 - n2), polyify(n0)
for nm, e in (("R2#", R2n), ("R1#", R1n), ("R0#", R0n)):
    print(nm, "=", sp.collect(e, s), "   [deg_s =", sp.Poly(e, s).degree(), "]")

R2n = sp.expand(R2n.subs(alsub)); R1n = sp.expand(R1n.subs(alsub)); R0n = sp.expand(R0n.subs(alsub))
P2n = -R2n.coeff(s, 1); Q2n = R2n.coeff(s, 2)
P1n = -R1n.coeff(s, 1); QQ1 = R1n.coeff(s, 2)
P0n = -R0n.coeff(s, 1); QQ0 = R0n.coeff(s, 2)
Q1n = sp.simplify(QQ1 - 2*Q2n)
Q0n = sp.simplify(QQ0 - Q2n - Q1n)

print()
print("=== the transformed row ===")
print("P#(n) =", sp.collect(sp.expand(P2n*n**2 + P1n*n + P0n), n))
print("Q#(n) =", sp.collect(sp.expand(Q2n*n**2 + Q1n*n + Q0n), n))
print()
print("consistency: R2# const =", sp.simplify(R2n.coeff(s,0)),
      "  R1# const =", sp.simplify(R1n.coeff(s,0)),
      "  R0# const =", sp.simplify(R0n.coeff(s,0)))

print()
print("=== ALPHA CONDITION ===")
print("the move gives a 3-term row of degree <=2  iff  (al-1+r1)(al-1+r2) = 0,")
print("i.e. al = 1-r1 or al = 1-r2, the two exponents of L at t = infinity.")

print()
print("=== Zagier class  P = a(n^2+n)+b,  Q = d n^2  (r1=r2=0, al=1) ===")
zg = {p1: lam+mu, p0: sp.Symbol('b'), r1: 0, r2: 0}
print("  P#(n) =", sp.collect(sp.expand((P2n*n**2+P1n*n+P0n).subs(zg)), n))
print("  Q#(n) =", sp.collect(sp.expand((Q2n*n**2+Q1n*n+Q0n).subs(zg)), n))
print("  ACF Theorem 1: (a#,b#,d#) = (mu-2lam, b-lam, lam^2-lam*mu)")

print()
print("=== root-row class  Q = d(n-1/2)^2, al = 1/2 ===")
rr = {r1: sp.Rational(1,2), r2: sp.Rational(1,2)}
print("  P#(n) =", sp.collect(sp.expand((P2n*n**2+P1n*n+P0n).subs(rr)), n))
print("  Q#(n) =", sp.collect(sp.expand((Q2n*n**2+Q1n*n+Q0n).subs(rr)), n))
