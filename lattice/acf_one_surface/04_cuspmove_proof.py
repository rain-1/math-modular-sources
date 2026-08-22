#!/usr/bin/env python3
"""Symbolic proof of the cusp-move lemma (Theorem 1 of ACF_ONE_SURFACE.md).

Row: (n+1)^2 u_{n+1} = P(n) u_n - Q(n) u_{n-1}, P = a(n^2+n)+b, Q = d n^2.
Picard-Fuchs operator in D = d/dt:
      t^2 R y'' + t Sf y' + (-b t + d t^2) y = 0,
      R = 1 - a t + d t^2,  Sf = 1 - 2a t + 3d t^2.
Write a = lam + mu, d = lam*mu (lam, mu the characteristic roots).
Claim: with  s = t/(1 - lam t)  (so t = s/(1 + lam s))  and  z(s) = (1 - lam t) y(t),
z satisfies the same shape of equation with
      a# = mu - 2 lam,   b# = b - lam,   d# = lam^2 - lam*mu = -lam(mu-lam).
"""
import sympy as sp

t, s, lam, mu, b = sp.symbols('t s lam mu b')

def op(y, x, A, B, D):
    R  = 1 - A*x + D*x**2
    Sf = 1 - 2*A*x + 3*D*x**2
    return x**2*R*sp.diff(y, x, 2) + x*Sf*sp.diff(y, x) + (-B*x + D*x**2)*y

a  = lam + mu
d  = lam*mu
aS = mu - 2*lam
dS = lam**2 - lam*mu
bS = b - lam

z = sp.Function('z')
# t as a function of s, and y(t) = (1 + lam s) z(s)
t_of_s = s/(1 + lam*s)
y_expr = (1 + lam*s)*z(s)

# original operator, everything pushed to the variable s
R  = 1 - a*t_of_s + d*t_of_s**2
Sf = 1 - 2*a*t_of_s + 3*d*t_of_s**2
# d/dt = (ds/dt) d/ds ; ds/dt = 1/(1-lam t)^2 = (1+lam s)^2
dsdt = (1 + lam*s)**2
y1 = sp.diff(y_expr, s)*dsdt
y2 = sp.diff(y1, s)*dsdt
orig = t_of_s**2*R*y2 + t_of_s*Sf*y1 + (-b*t_of_s + d*t_of_s**2)*y_expr

new = op(z(s), s, aS, bS, dS)

ratio = sp.simplify(sp.expand(sp.simplify(orig)) / sp.expand(new))
print("orig / new =", sp.simplify(ratio))
print("difference (should be 0):", sp.simplify(sp.expand(orig - ratio*new)))

# also: the exponents / singular points
print()
print("original char. roots :", sp.solve(sp.Symbol('L')**2 - a*sp.Symbol('L') + d, sp.Symbol('L')))
print("new  char. roots     :", sp.solve(sp.Symbol('L')**2 - aS*sp.Symbol('L') + dS, sp.Symbol('L')))

# instances
print()
for (aa,bb,dd,name) in [(10,3,9,'C'),(7,2,-8,'A'),(17,6,72,'F'),(12,4,32,'E'),(9,3,27,'B'),(11,3,-1,'D')]:
    L = sp.symbols('L')
    rts = sp.solve(L**2 - aa*L + dd, L)
    out = []
    for r in rts:
        m = aa - r
        out.append((sp.simplify(m - 2*r), sp.simplify(bb - r), sp.simplify(r**2 - r*m)))
    print(f"{name} (a,b,d)=({aa},{bb},{dd}) roots={rts} -> cusp-moved rows {out}")
