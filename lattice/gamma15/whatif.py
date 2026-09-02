# -*- coding: utf-8 -*-
"""How much would the ONE open loophole be worth?

We exclude every non-principal h-preimage of Y_v (the conservative reading of CDT's
Lemma A.4.4).  In fact only those preimages at which the analytically continued conditional
function actually acquires its logarithm need be excluded.  This script quantifies the gain
if the k deepest conjugate pairs turned out to be harmless at each place."""
import numpy as np, math, json, contour as C, bcfast as BF, opt2, places
from region2 import Reg2
import itertools

m = 14; TAU = places.TAU

def eval_drop(nm, sa, Y, R, c, drop):
    """drop = number of deepest bad preimages (individually, not pairs) to ignore."""
    principal, pts = opt2.bad_points(Y, R)
    pts = sorted(pts, key=abs)[drop:]
    reg, order = opt2.best_design(R, [(0.0, c)], pts, 1.0, 1.0)
    rad = reg.rad(); lp = math.log(256*sa*rad)
    bc = BF.BC_fast(reg, lp, N=8192, cdmax=30, rmax=R*(1+1e-9))
    return rad, lp, bc, m*lp-bc, len(pts)

if __name__ == '__main__':
    print('Baseline (nothing dropped): J1 = 91.17577 (v1, R=0.78, c=6.5), J2 = 24.53386 (v2, R=0.76, c=7.5)')
    print('margin = (J1+J2)/2 - 14 tau = %+.5f' % ((91.17577+24.53386)/2 - 14*TAU))
    print()
    best = {}
    for (nm, sa, Y) in places.PLACES:
        print('=== %s ===' % nm)
        for drop in (0, 2, 4):
            b = None
            for R in (0.70, 0.76, 0.80, 0.84, 0.88, 0.92):
                for c in (6.5, 7.5, 10.0):
                    rad, lp, bc, J, npts = eval_drop(nm, sa, Y, R, c, drop)
                    if b is None or J > b[3]: b = (R, c, rad, J, lp, bc, npts)
            print('  drop %d deepest: best R=%.2f c=%.1f  #cut=%d  |psi\'|=%.6f  L=%.6f  BC=%.5f  14L-BC=%.5f'
                  % (drop, b[0], b[1], b[6], b[2], b[4], b[5], b[3]))
            best[(nm, drop)] = b[3]
        print()
    print('resulting averaged margins:')
    for d1 in (0, 2, 4):
        for d2 in (0, 2, 4):
            J1 = best[('v1', d1)]; J2 = best[('v2', d2)]
            print('  drop %d at v1, %d at v2:  margin = %+.5f' % (d1, d2, (J1+J2)/2 - 14*TAU))
