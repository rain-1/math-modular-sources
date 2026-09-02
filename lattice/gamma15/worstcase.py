# -*- coding: utf-8 -*-
"""Does the verdict depend on any UNPROVED exclusion?

`monodromy.py` proves that the pair of preimages with bottom rows (2,+-1) must be excluded
(the parabolic fixed lines u_{t1} and u_{t2} are distinct, so continuing the conditional
function around the cusp y = infinity destroys its regularity at the fold).  Every other
non-principal preimage is excluded here only by the conservative reading of CDT's Lemma A.4.4.

This script therefore RELAXES the problem to its extreme: only the (2,+-1) pair is REQUIRED
to be excluded; every other bad preimage may be slit or not, whichever is better.  Extra slits
are always admissible (they only shrink Omega), so the search is over
   base radius R, lune parameter c, and the set of optional points with |z| < T (T swept).
If the margin is still negative, the verdict does not rest on any unproved exclusion.
"""
import numpy as np, math, contour as C, bcfast as BF, places
from region2 import Reg2
import itertools

m = 14; TAU = places.TAU

def run():
    best = {}
    for (nm, sa, Y) in places.PLACES:
        pre = C.preimages_h(Y, rmax=0.999, cdmax=60)
        req = [z for (z, cd, r) in pre[1:] if cd in ((2, 1), (2, -1))]
        opt = [z for (z, cd, r) in pre[1:] if cd not in ((2, 1), (2, -1))]
        print('  %s: required exclusions |z| = %s' % (nm, [round(abs(z), 6) for z in req]))
        b = None
        for R in (0.60, 0.68, 0.74, 0.78, 0.82, 0.86, 0.90, 0.94, 0.97):
            for c in (5.0, 6.5, 7.5, 10.0, 20.0):
                for T in (0.0, 0.30, 0.55, 0.70, 0.80, 0.90, 1.0):
                    pts = [z for z in req if abs(z) < R] + \
                          [z for z in opt if abs(z) < min(R, T)]
                    if len(pts) > 9: continue
                    n = len(pts)
                    orders = list(itertools.permutations(range(n))) if math.factorial(max(n,1)) <= 720 \
                             else [tuple(np.argsort([-abs(p) for p in pts])), tuple(np.argsort([abs(p) for p in pts]))]
                    for o in orders:
                        rg = Reg2(R); rg.add_lune_at(0.0, c)
                        for i in o: rg.add_slit_to(pts[i])
                        rad = rg.rad()
                        lp = math.log(256*sa*rad)
                        bc = BF.BC_fast(rg, lp, N=4096, cdmax=30, rmax=R*(1+1e-9))
                        J = m*lp - bc
                        if b is None or J > b[0]: b = (J, R, c, T, rad, lp, bc, n)
        print('     best: R=%.2f c=%.1f T=%.2f  #slits=%d  |psi\'|=%.6f  L=%.6f  BC=%.5f  14L-BC=%.5f'
              % (b[1], b[2], b[3], b[7], b[4], b[5], b[6], b[0]))
        best[nm] = b[0]
    print()
    print('  RELAXED margin (only the proved exclusions enforced) = %+.5f'
          % ((best['v1']+best['v2'])/2 - 14*TAU))
    return best

if __name__ == '__main__':
    print('Relaxed problem: only the (2,+-1) pair is required to be excluded.')
    run()
