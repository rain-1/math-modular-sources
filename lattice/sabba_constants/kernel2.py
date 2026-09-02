"""Kernel equation  A'' + 2A'/(t(1+a t^2)) - (1 - t^-4) A = 0.
Regular singular points at t^2 = -1/a.  Exponents 0,2.  Claim: APPARENT.
Reason: at such t0, t0^4 = 1/a^2 = 1, so the zeroth-order coefficient Q = -(1-t^-4)
vanishes there; the u^{-1} relation forces c1=0 and the only obstruction (at u^0)
is Q(t0) c0 = 0.  Verified below symbolically and by numerical monodromy."""
import sympy as sp
from mpmath import mp, mpf, mpc, odefun, nstr, exp, pi, j as I
u = sp.Symbol('u')
for aval, t0 in ((1, sp.I), (1, -sp.I), (-1, sp.Integer(1)), (-1, sp.Integer(-1))):
    t = t0 + u
    P = sp.series(sp.simplify(2/(t*(1+aval*t**2))), u, 0, 3).removeO()
    Q = sp.series(sp.simplify(-(1 - t**(-4))), u, 0, 3).removeO()
    pm1 = sp.simplify(sp.limit(sp.simplify(2/(t*(1+aval*t**2)))*u, u, 0))
    q0  = sp.simplify(Q.subs(u,0))
    print("a=%+d  t0=%-3s : residue of P = %s (indicial rho(rho-1)+rho*p_-1 = rho(rho+p_-1-1))"%(aval,t0,pm1),
          " exponents:", sp.solve(sp.Symbol('r')*(sp.Symbol('r')-1)+pm1*sp.Symbol('r'), sp.Symbol('r')),
          " Q(t0) =", q0, "-> obstruction Q(t0)*c0 =", q0)
print()
print("numerical monodromy check (a=-1, loop around t0=1, radius 0.3):")
mp.dps = 30
a = -1
def rhs(x, y):
    tt = mpc(1) + mpf('0.3')*exp(2j*pi*x)      # circle around t0=1
    dt = 2j*pi*mpf('0.3')*exp(2j*pi*x)
    A, Ap = y
    App = -(2*Ap/(tt*(1+a*tt**2))) + (1 - tt**(-4))*A
    return [Ap*dt, App*dt]
for ic in ([mpc(1),mpc(0)], [mpc(0),mpc(1)]):
    f = odefun(rhs, 0, ic, tol=mpf(10)**-24)
    v = f(1)
    print("   start", [nstr(c,4) for c in ic], "-> after one loop", [nstr(c,20) for c in v],
          " |deviation| =", nstr(max(abs(v[0]-ic[0]), abs(v[1]-ic[1])), 6))
