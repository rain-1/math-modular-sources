# -*- coding: utf-8 -*-
"""RELAXED problem: only the PROVED exclusions (bottom rows (2,+-1), whose singularity is
established in monodromy.py) are enforced; every other non-principal preimage is assumed
harmless.  Extra slits are always admissible, so the search is over base radius R, lune
parameter c, and the set of optional points with |z| < T.
If the margin is still negative, the verdict rests on no unproved exclusion."""
import numpy as np, math, itertools, contour as C, bcfast as BF, places
from region2 import Reg2

m = 14; TAU = places.TAU
best = {}
for (nm, sa, Y) in places.PLACES:
    pre = C.preimages_h(Y, rmax=0.999, cdmax=60)
    req = [z for (z, cd, r) in pre[1:] if cd in ((2, 1), (2, -1))]
    opt = [z for (z, cd, r) in pre[1:] if cd not in ((2, 1), (2, -1))]
    print(' ', nm, ': required |z| =', [round(abs(z), 6) for z in req], flush=True)
    b = None
    for R in (0.60, 0.70, 0.80, 0.88, 0.94, 0.97, 0.99):
        for c in (5.0, 7.5, 12.0):
            for T in (0.0, 0.55, 0.75, 0.85, 1.0):
                pts = [z for z in req if abs(z) < R] + [z for z in opt if abs(z) < min(R, T)]
                n = len(pts)
                if n > 8: continue
                ords = list(itertools.permutations(range(n))) if math.factorial(max(n, 1)) <= 720 \
                       else [tuple(np.argsort([-abs(p) for p in pts]))]
                bo, br = None, -1.0
                for o in ords:
                    rg = Reg2(R); rg.add_lune_at(0.0, c)
                    for i in o: rg.add_slit_to(pts[i])
                    if rg.rad() > br: br, bo = rg.rad(), o
                rg = Reg2(R); rg.add_lune_at(0.0, c)
                for i in bo: rg.add_slit_to(pts[i])
                lp = math.log(256*sa*rg.rad())
                bc = BF.BC_fast(rg, lp, N=4096, cdmax=30, rmax=R*(1+1e-9))
                J = m*lp - bc
                if b is None or J > b[0]: b = (J, R, c, T, rg.rad(), lp, bc, n)
    print('    best: R=%.2f c=%.1f T=%.2f #slits=%d |psi|=%.6f L=%.6f BC=%.5f 14L-BC=%.5f'
          % (b[1], b[2], b[3], b[7], b[4], b[5], b[6], b[0]), flush=True)
    best[nm] = b[0]
print()
print('  RELAXED margin = %+.5f   (fully constrained: -1.4416)'
      % ((best['v1']+best['v2'])/2 - 14*TAU), flush=True)
