#!/usr/bin/env python3
"""Exponent dictionary for FOUR-term rank-2 rows with FIVE singular points.

Row:  (n+1)^2 u_{n+1} = P(n) u_n - Q(n) u_{n-1} + R(n) u_{n-2},  u_0=1, u_{-1}=u_{-2}=0
      P = a n^2 + b n + c,  Q = d n^2 + e n + f,  R = g n^2 + h n + j.

Operator  L = theta^2 - t P(theta) + t^2 Q(theta+1) - t^3 R(theta+2),  theta = t d/dt.
"""
import sympy as sp

t, n, mu, rho = sp.symbols('t n mu rho')
a,b,c,d,e,f,g,h,j = sp.symbols('a b c d e f g h j')

# --- 1. check the operator convention on the series level -------------------
def series_check():
    N = 12
    P = lambda k: a*k**2 + b*k + c
    Q = lambda k: d*k**2 + e*k + f
    R = lambda k: g*k**2 + h*k + j
    u = {-2: 0, -1: 0, 0: 1}
    for k in range(0, N):
        u[k+1] = sp.expand((P(k)*u[k] - Q(k)*u[k-1] + R(k)*u[k-2])/sp.Integer(k+1)**2)
    y = sum(u[k]*t**k for k in range(N))
    th = lambda F: sp.expand(t*sp.diff(F, t))
    L = th(th(y)) - t*(a*th(th(y)) + b*th(y) + c*y) \
        + t**2*(d*th(th(y)) + (2*d+e)*th(y) + (d+e+f)*y) \
        - t**3*(g*th(th(y)) + (4*g+h)*th(y) + (4*g+2*h+j)*y)
    L = sp.expand(L)
    bad = [k for k in range(N-2) if sp.simplify(L.coeff(t, k)) != 0]
    return bad

# --- 2. the D-form coefficients --------------------------------------------
# L = t^2 Rc(t) D^2 + t Sc(t) D + t Vc(t)
Rc = 1 - a*t + d*t**2 - g*t**3
Sc = 1 - (a+b)*t + (3*d+e)*t**2 - (5*g+h)*t**3
Vc = -c + (d+e+f)*t - (4*g+2*h+j)*t**2

def dform_check():
    y = sp.Function('y')(t)
    th = lambda F: t*sp.diff(F, t)
    L1 = th(th(y)) - t*(a*th(th(y)) + b*th(y) + c*y) \
         + t**2*(d*th(th(y)) + (2*d+e)*th(y) + (d+e+f)*y) \
         - t**3*(g*th(th(y)) + (4*g+h)*th(y) + (4*g+2*h+j)*y)
    L2 = t**2*Rc*sp.diff(y,t,2) + t*Sc*sp.diff(y,t) + t*Vc*y
    return sp.simplify(sp.expand(L1 - L2))

T = sp.expand(Sc - t*sp.diff(Rc, t))     # = 1 - b t + (d+e) t^2 - (2g+h) t^3

if __name__ == '__main__':
    print("series identity Ly=0, nonzero coefficients:", series_check())
    print("D-form difference:", dform_check())
    print("T(t) = S - t R' =", sp.factor(T))
    # exponent at a root t_i of Rc:   rho_i = - T(t_i)/(t_i Rc'(t_i))
    # sum of rho_i via residues:
    ti = sp.symbols('t1 t2 t3')
    # generic cubic with roots t1,t2,t3 and Rc(0)=1:
    Rgen = sp.expand((1-t/ti[0])*(1-t/ti[1])*(1-t/ti[2]))
    print("\nRc in terms of roots:", sp.collect(Rgen, t))
    # sum rho_i  should be  s1+s2-1 = -h/g - 1
    tot = 0
    for i in range(3):
        others = [ti[k] for k in range(3) if k != i]
        Rp = sp.diff(Rgen, t).subs(t, ti[i])
        Tt = T.subs(t, ti[i])
        tot += -Tt/(ti[i]*Rp)
    # substitute a,d,g by the elementary symmetric functions of 1/t_i
    l = [1/x for x in ti]
    subs = {a: sum(l), d: l[0]*l[1]+l[0]*l[2]+l[1]*l[2], g: l[0]*l[1]*l[2]}
    tot = sp.simplify(sp.expand(sp.simplify(tot.subs(subs))))
    print("sum rho_i =", tot, "   (should be -1 - h/g)")
    print("check:", sp.simplify(tot - (-1 - h/g)))
    # ---- equal-exponent class -------------------------------------------
    cond = sp.expand(T + rho*t*sp.diff(Rc, t) - Rc)
    print("\nT + rho t R' - R  =", sp.collect(cond, t))
    sol = sp.solve([cond.coeff(t,1), cond.coeff(t,2), cond.coeff(t,3)], [b,e,h], dict=True)
    print("all rho_i = rho  <=>", sol)
