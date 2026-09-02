# -*- coding: utf-8 -*-
"""Design and evaluate a contour for one place of a CDT-type configuration."""
import numpy as np, math, cmath, itertools, contour as C, bcfast as BF
from region2 import Reg2

def bad_points(Y, R, cdmax=90):
    pre = C.preimages_h(Y, rmax=0.9999999, cdmax=cdmax)
    return pre[0][0], [z for (z, cd, r) in pre[1:] if abs(z) < R]

def build(R, lunes, pts, order, over=1.0, rho=1.0):
    reg = Reg2(R, rho=rho)
    for (th, c) in lunes: reg.add_lune_at(th, c)
    for i in order: reg.add_slit_to(pts[i], over=over)
    return reg

def best_design(R, lunes, pts, over=1.0, rho=1.0, maxperm=5040, ntrial=400):
    n = len(pts)
    if n == 0: return build(R, lunes, pts, (), over, rho), ()
    if math.factorial(n) <= maxperm:
        bv, bo = -1, None
        for o in itertools.permutations(range(n)):
            v = build(R, lunes, pts, o, over, rho).rad()
            if v > bv: bv, bo = v, o
        return build(R, lunes, pts, bo, over, rho), bo
    best = (-1, None); rng = np.random.default_rng(0)
    for trial in range(ntrial):
        if trial == 0:
            reg = Reg2(R, rho=rho)
            for (th, c) in lunes: reg.add_lune_at(th, c)
            rem = list(range(n)); o = []
            while rem:
                vals = [(abs(complex(reg._inv_raw(np.array([pts[j]]))[0])), j) for j in rem]
                vals.sort(reverse=True); j = vals[0][1]
                reg.add_slit_to(pts[j], over=over); o.append(j); rem.remove(j)
            v = reg.rad()
        else:
            o = tuple(rng.permutation(n)); v = build(R, lunes, pts, o, over, rho).rad()
        if v > best[0]: best = (v, tuple(o))
    return build(R, lunes, pts, best[1], over, rho), best[1]

def evaluate(R, lunes, Y, s_abs, m, tau, over=1.0, rho=1.0, N=8192, cdmax=30, verify=True):
    principal, pts = bad_points(Y, R)
    reg, order = best_design(R, lunes, pts, over, rho)
    rad = reg.rad()
    lp = math.log(256*s_abs*rad)
    bcv, shape, av, r0, r1 = BF.BC_fast(reg, lp, N=N, cdmax=cdmax, rmax=R*(1+1e-9),
                                        return_parts=True)
    out = dict(R=R, lunes=lunes, nbad=len(pts), nslit=sum(1 for o in reg.ops if o[0] == 's'),
               psi=rad, logphip=lp, BC=bcv, shape=shape, J=m*lp-bcv, entry=lp-tau,
               reg=reg, order=order, pts=pts, principal=principal)
    if verify:
        ok = True; worst = 2.0
        for p in pts:
            e, v = reg.excluded(p); ok = ok and e; worst = min(worst, v)
        pe, pv = reg.excluded(principal)
        out['ok'] = ok; out['worst_excl'] = worst; out['principal_inside'] = (not pe)
    return out
