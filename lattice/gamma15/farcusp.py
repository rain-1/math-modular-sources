# -*- coding: utf-8 -*-
"""The FAR-cusp periods, and what the second real place actually demands.

At place v2 the pure module is B_i(y/sigma(s)) = B_i(y/t1), so x = t1 = sigma(s) is where
singularities are allowed (it maps to y = infinity) and **x = t2 is the point that must be
removed** -- the v2-fold.  The conditional function at v2 is
    sigma(H) = sigma(a) A + sigma(b) B_D + sigma(c) B'_new,   B'_new = sigma(B_new),
and the requirement is that it be REGULAR AT t2.

Since the inhomogeneities are single-valued rational functions, for any f in
{A, B_D, B'_new} the monodromy difference Delta_f := (M_{t2} - 1) f solves the homogeneous
equation and is fixed by M_{t2}, hence lies in the line <u_{t2}> spanned by the homogeneous
solution holomorphic at t2.  Writing Delta_f = kappa_f * u_{t2}, the regularity condition is
    lambda kappa_A + mu kappa_{B_D} + nu kappa_{B'} = 0,
so the FAR-CUSP PERIODS are  pi_D := kappa_{B_D}/kappa_A  and  pi' := kappa_{B'}/kappa_A.
(The near-cusp analogues are the Apery limits zeta(2)/5 and xi.)

This script computes them by numerical analytic continuation around t2.
"""
import numpy as np, math
from fractions import Fraction as Fr
from scipy.integrate import solve_ivp

s5 = 5**0.5
t1 = (-11+5*s5)/2
t2 = (-11-5*s5)/2

def load(fn='farcusp_data.txt'):
    A, BD, U, V = [], [], [], []
    for line in open(fn):
        if line.startswith('#'): continue
        p = line.split()
        A.append(Fr(p[1])); BD.append(Fr(p[2])); U.append(Fr(p[3])); V.append(Fr(p[4]))
    return A, BD, U, V

def evalser(c, x):
    v = 0.0; vp = 0.0
    for n in range(len(c)-1, -1, -1):
        v = v*x + float(c[n])
    for n in range(len(c)-1, 0, -1):
        vp = vp*x + n*float(c[n])
    return v, vp

# operator L y = x(1-11x-x^2) y'' + (1-22x-3x^2) y' - (3+x) y
P2 = lambda x: x*(1-11*x-x*x)
P1 = lambda x: 1-22*x-3*x*x
P0 = lambda x: -(3+x)

def check_rhs(c, name, N=40):
    """Apply L to the power series and print the first coefficients of the result."""
    out = []
    for n in range(N):
        # [x^n] of L y
        v = 0.0
        def a(k): return float(c[k]) if 0 <= k < len(c) else 0.0
        # x*y'' : coeff n  -> (n)(n-1) a_n ... careful: x*y'' = sum n(n+1) a_{n+1} x^n
        v += (n+1)*n*a(n+1)
        v += -11*(n)*(n-1)*a(n)          # -11x^2 y''
        v += -1*(n-1)*(n-2)*a(n-1)       # -x^3 y''
        v += (n+1)*a(n+1)                # y'
        v += -22*n*a(n)                  # -22x y'
        v += -3*(n-1)*a(n-1)             # -3x^2 y'
        v += -3*a(n)                     # -3y
        v += -1*a(n-1)                   # -x y
        out.append(v)
    print('   L(%s) = %s' % (name, ['%.6g' % t for t in out[:8]]))
    return out

def continue_loop(rhs, x0, Y0, loop, rtol=1e-12, atol=1e-14):
    """continue (y,y') along the polyline `loop` for  P2 y'' + P1 y' + P0 y = rhs(x)."""
    Y = list(Y0); cur = x0
    for nxt in loop:
        a, b = cur, nxt
        def f(t, Y, a=a, b=b):
            z = a + t*(b-a); dz = b-a
            y = Y[0]+1j*Y[1]; yp = Y[2]+1j*Y[3]
            ypp = (rhs(z) - P1(z)*yp - P0(z)*y)/P2(z)
            dy, dyp = yp*dz, ypp*dz
            return [dy.real, dy.imag, dyp.real, dyp.imag]
        sol = solve_ivp(f, [0, 1], Y, rtol=rtol, atol=atol, method='DOP853')
        Y = [sol.y[i][-1] for i in range(4)]
        cur = b
    return (Y[0]+1j*Y[1], Y[2]+1j*Y[3])

if __name__ == '__main__':
    A, BD, U, V = load()
    print('series loaded: %d coefficients' % len(A))
    print('right-hand sides of the operator (should be 0, 1, and 2/(1 - x/t2) resp. 2/(1 - x/t1)):')
    check_rhs(A, 'A')
    check_rhs(BD, 'B_D')
    Bn  = [U[i] + V[i]*0 for i in range(len(U))]        # placeholder
    # B_new = U + V sqrt5 ; B'_new = U - V sqrt5 -- as REAL series
    Bnew  = [float(U[i]) + float(V[i])*s5 for i in range(len(U))]
    Bnewp = [float(U[i]) - float(V[i])*s5 for i in range(len(U))]
    print('   (numeric) L(B_new)  first coeffs: %s   vs 2*(1/t2)^n = %s'
          % (['%.6f' % t for t in check_rhs(Bnew, 'B_new', 8)[:5]],
             ['%.6f' % (2*(1/t2)**n) for n in range(5)]))
    print('   (numeric) L(B\'_new) first coeffs: %s   vs 2*(1/t1)^n = %s'
          % (['%.6f' % t for t in check_rhs(Bnewp, "B'_new", 8)[:5]],
             ['%.6f' % (2*(1/t1)**n) for n in range(5)]))

    x0 = 0.05
    loop = [0.05+2j, -6.0+2j, -6.0+0j]
    th0 = np.angle(-6.0-t2)
    loop += [t2 + 5.0*np.exp(1j*(th0+2*np.pi*k/96)) for k in range(0, 97)]
    loop += [-6.0+2j, 0.05+2j, 0.05+0j]

    print()
    print('monodromy differences around t2 = %.10f (base point x0 = %.3f):' % (t2, x0))
    res = {}
    for name, c, rhs in (('A', [float(t) for t in A], lambda x: 0.0),
                         ('B_D', [float(t) for t in BD], lambda x: 1.0),
                         ("B'_new", Bnewp, lambda x: 2.0/(1 - x/t1))):
        y0, yp0 = evalser(c, x0)
        y1, yp1 = continue_loop(rhs, x0, [y0, 0.0, yp0, 0.0], loop)
        d = (y1-y0, yp1-yp0)
        res[name] = d
        print('   %-8s  f(x0) = %+.12f   Delta = (%+.10e, %+.10e)'
              % (name, y0, d[0].real if abs(d[0].imag) < 1e-6*max(1,abs(d[0])) else d[0], d[1].real))
    # the three Delta's must be proportional (all in <u_t2>)
    dA, dD, dN = res['A'], res['B_D'], res["B'_new"]
    print()
    print('   proportionality check (all Delta must span the same line <u_t2>):')
    for nm, d in (('B_D/A', dD), ("B'_new/A", dN)):
        r0 = d[0]/dA[0]; r1 = d[1]/dA[1]
        print('      %-10s  Delta_0 ratio = %+.12f   Delta_1 ratio = %+.12f   |diff| = %.2e'
              % (nm, r0.real, r1.real, abs(r0-r1)))
    piD = (dD[0]/dA[0]).real; pip = (dN[0]/dA[0]).real
    print()
    print('   FAR-CUSP PERIODS:  pi_D = %.12f      pi\' = %.12f' % (piD, pip))
    import mpmath as mp
    mp.mp.dps = 20
    z2 = float(mp.zeta(2))
    ReL = 0.9587161227168831553919364293311785264159715307582960686724044479103954560598657
    ImL = 0.1455658767850895904617045118119864537208051468890990138342266111925591477167064
    xi  = 0.6556341888406567663309814138723994024111384591367020284631349608915842669423314
    xip = -(5**0.5*0+1)*ReL - ((-11+5*5**0.5)/2)*0 - ImL*((11-5*5**0.5)/2)  # -ReL - phi^-5 ImL
    xip = -ReL - ((5*5**0.5-11)/2)*ImL
    print('   candidates:  zeta(2)/5 = %.12f   xi = %.12f   xi\' = -ReL - phi^-5 ImL = %.12f'
          % (z2/5, xi, xip))
    for nm, val in (('pi_D', piD), ("pi'", pip)):
        print('   lindep(%s): %s' % (nm, mp.pslq([mp.mpf(val), 1, mp.mpf(z2), mp.mpf(ReL),
                                                 mp.mpf(ImL), mp.sqrt(5), mp.sqrt(5)*mp.mpf(ReL),
                                                 mp.sqrt(5)*mp.mpf(ImL)], maxcoeff=10**6, maxsteps=10**5)))

def near_cusp():
    """Calibration: the SAME construction at the NEAR cusp t1 must return the Apery limits."""
    A, BD, U, V = load()
    Bnew  = [float(U[i]) + float(V[i])*s5 for i in range(len(U))]
    Bnewp = [float(U[i]) - float(V[i])*s5 for i in range(len(U))]
    x0 = 0.05
    loop = [0.05+0.06j] + [t1 + 0.05*np.exp(1j*(np.angle(0.05+0.06j-t1)+2*np.pi*k/96))
                           for k in range(0, 97)] + [0.05+0.06j, 0.05+0j]
    res = {}
    for name, c, rhs in (('A', [float(t) for t in A], lambda x: 0.0),
                         ('B_D', [float(t) for t in BD], lambda x: 1.0),
                         ('B_new', Bnew,  lambda x: 2.0/(1 - x/t2)),
                         ("B'_new", Bnewp, lambda x: 2.0/(1 - x/t1))):
        y0, yp0 = evalser(c, x0)
        y1, yp1 = continue_loop(rhs, x0, [y0, 0.0, yp0, 0.0], loop)
        res[name] = (y1-y0, yp1-yp0)
    dA = res['A']
    print()
    print('CALIBRATION at the NEAR cusp t1 = %.10f:' % t1)
    for nm in ('B_D', 'B_new', "B'_new"):
        d = res[nm]
        r0 = (d[0]/dA[0]).real; r1 = (d[1]/dA[1]).real
        print('   period of %-8s at t1 = %+.12f   (consistency %.1e)' % (nm, r0, abs(r0-r1)))
    print('   zeta(2)/5 = %+.12f    xi = %+.12f' % (float(__import__("mpmath").zeta(2))/5,
                                                     0.6556341888406568))
