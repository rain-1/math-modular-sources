# -*- coding: utf-8 -*-
"""Closing the one loophole: are the (2,+-1) preimages really singular for the lift?

Set-up.  The lift of the conditional function to the disc is regular at the deck translate
gamma.z0 iff rho(gamma)H is regular at the fold t1.  The two deepest bad preimages -- at
CDT's host (|z| = 0.401921) and at BOTH of ours (|z| = 0.536032 at v1, 0.213693 at v2) --
have bottom rows (2,+-1), i.e. they are the images of z0 under the parabolic generator
[[1,0],[+-2,1]] of the stabiliser of the cusp tau = 0 of Gamma_0(2), which is the cusp
y = infinity.  The fibre of y = x^2/(x-s) over y = infinity is {s, infinity}, two points,
so the loop lifts to the x-line as a loop around the OUTER cusp t2 = s.  Hence

   rho(gamma)H = H + Delta,   Delta = (M_{t2} - 1)H .

Delta solves the HOMOGENEOUS equation (H and rho(gamma)H satisfy the same inhomogeneous one
with single-valued rational right-hand side) and is killed by M_{t2} (unipotent), so
Delta lies in the line ker(M_{t2}-1) = <u_{t2}>.  Delta != 0 because H is singular at t2
(the radius of convergence of the fold-regular linear form is exactly |t2|).  Therefore

   rho(gamma)H is regular at t1  <==>  u_{t2} is regular at t1  <==>  <u_{t1}> = <u_{t2}>.

This script computes the monodromy matrices M_{t1}, M_{t2} of the Picard-Fuchs operator
numerically, on CDT's host and on Zagier's row D, and tests whether the two parabolic fixed
lines coincide.  (They do not: the monodromy group is the image of a congruence subgroup and
M_{t1}, M_{t2} are parabolic generators at DISTINCT cusps, whose fixed points on the boundary
of H -- hence whose fixed lines -- are distinct.)
"""
import numpy as np
from scipy.integrate import solve_ivp

def monodromy(P2, P1, P0, x0, loop, rtol=1e-12, atol=1e-14):
    """Monodromy matrix along the closed polyline `loop` (list of complex points, starting
    and ending at x0), in the basis (y,y')(x0) = (1,0), (0,1)."""
    def seg(a, b, Y):
        def f(t, Y):
            z = a + t*(b-a); dz = b-a
            y = Y[0]+1j*Y[1]; yp = Y[2]+1j*Y[3]
            ypp = -(P1(z)*yp + P0(z)*y)/P2(z)
            dy, dyp = yp*dz, ypp*dz
            return [dy.real, dy.imag, dyp.real, dyp.imag]
        s = solve_ivp(f, [0, 1], Y, rtol=rtol, atol=atol, method='DOP853')
        return [s.y[i][-1] for i in range(4)]
    cols = []
    for init in ((1, 0), (0, 1)):
        Y = [float(init[0]), 0.0, float(init[1]), 0.0]
        cur = x0
        for nxt in loop:
            Y = seg(cur, nxt, Y); cur = nxt
        cols.append([Y[0]+1j*Y[1], Y[2]+1j*Y[3]])
    return np.array(cols).T          # columns = images of the basis vectors

def circle(c, r, x0, n=64, sign=+1):
    """polyline: x0 -> c+r*e^{i th0} -> around the circle -> back to x0."""
    th0 = np.angle(x0-c)
    pts = [c + r*np.exp(1j*th0)]
    for k in range(1, n+1):
        pts.append(c + r*np.exp(1j*(th0 + sign*2*np.pi*k/n)))
    pts.append(x0)
    return pts

def fixed_line(M):
    """the (1-dimensional) kernel of M - I, as a normalised vector."""
    A = M - np.eye(2)
    u, s, vh = np.linalg.svd(A)
    v = vh[-1].conj()
    return v/np.linalg.norm(v), s

def report(name, P2, P1, P0, x0, loops):
    print('=== %s ===' % name)
    Ms = {}
    for tag, loop in loops.items():
        M = monodromy(P2, P1, P0, x0, loop)
        Ms[tag] = M
        ev = np.linalg.eigvals(M)
        print('  M_%s  = [[%.6f%+.6fj, %.6f%+.6fj], [%.6f%+.6fj, %.6f%+.6fj]]'
              % (tag, M[0,0].real, M[0,0].imag, M[0,1].real, M[0,1].imag,
                 M[1,0].real, M[1,0].imag, M[1,1].real, M[1,1].imag))
        print('        det = %.10f%+.10fj   trace = %.8f%+.8fj   eigs = %s'
              % (np.linalg.det(M).real, np.linalg.det(M).imag,
                 np.trace(M).real, np.trace(M).imag,
                 ', '.join('%.6f%+.6fj' % (e.real, e.imag) for e in ev)))
        v, s = fixed_line(M)
        print('        parabolic? |trace-2| = %.2e ; smallest singular value of M-I = %.2e'
              % (abs(np.trace(M)-2), s[-1]))
        print('        fixed line u_%s ~ (%.8f%+.8fj, %.8f%+.8fj)'
              % (tag, v[0].real, v[0].imag, v[1].real, v[1].imag))
        Ms[tag+'_v'] = v
    v1, v2 = Ms['t1_v'], Ms['t2_v']
    cross = abs(v1[0]*v2[1] - v1[1]*v2[0])
    print('  |u_t1 ^ u_t2| = %.10f   ->  the two parabolic fixed lines are %s'
          % (cross, 'DISTINCT (loophole closed)' if cross > 1e-6 else 'EQUAL (loophole open!)'))
    print()
    return Ms

if __name__ == '__main__':
    # ---- CDT's host: x(1-x)(1-9x) y'' + (1-20x+27x^2) y' + 3(3x-1) y = 0  (Prop 11.1.4)
    P2 = lambda x: x*(1-x)*(1-9*x)
    P1 = lambda x: 1-20*x+27*x*x
    P0 = lambda x: 3*(3*x-1)
    x0 = 0.5+0j
    loops = {'t1': circle(1/9, 0.08, x0), 't2': circle(1.0, 0.4, x0)}
    report("CDT's host (Zagier C on Gamma_0(6)): t1 = 1/9, t2 = 1", P2, P1, P0, x0, loops)

    # ---- Zagier D on Gamma_1(5):  x(1-11x-x^2) y'' + (1-22x-3x^2) y' - (3+x) y = 0
    Q2 = lambda x: x*(1-11*x-x*x)
    Q1 = lambda x: 1-22*x-3*x*x
    Q0 = lambda x: -(3+x)
    t1 = (-11+5*5**0.5)/2; t2 = (-11-5*5**0.5)/2
    print('  Zagier D singular points: t1 = %.10f, t2 = %.10f' % (t1, t2))
    x0 = 0.5+0j
    L1 = circle(t1, 0.045, x0)
    # a loop around t2 = -11.09 that avoids x = 0: go up, across, down, circle, back
    big = [0.5+3j, -6.0+3j, -6.0+0j]
    L2 = big + circle(t2, 5.0, -6.0+0j) + [-6.0+3j, 0.5+3j, x0]
    report("Zagier D (Gamma_1(5)): t1 = phi^-5, t2 = -phi^5", Q2, Q1, Q0, x0,
           {'t1': L1, 't2': L2})
