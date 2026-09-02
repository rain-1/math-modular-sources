# -*- coding: utf-8 -*-
"""The variant that needs only ONE hypothesis.

The far-cusp periods (farcusp.py) show that the conditional function at the second place,
sigma(H) = sigma(a) A + sigma(b) B_D + sigma(c) B'_new, is regular at the v2-fold t2 iff
    (ii)   sigma(a) - (11/5) sigma(b) zeta(2) + sigma(c) xi' = 0,     xi' = -Re L - phi^-5 Im L,
whereas the target hypothesis is
    (i)    a + (1/5) b zeta(2) + c xi = 0,                            xi  = -Re L + phi^5  Im L.
(ii) is NOT sigma of (i): they differ exactly in the zeta(2) coefficient, +1/5 against -11/5.
So one K-relation gives a conditional function at ONE place only.

There is, however, a variant that survives on the single hypothesis: at v2 one may simply
require Omega_2 to avoid EVERY preimage of Y_2 -- the principal one included.  Then
sigma(G) pulls back holomorphically whatever its behaviour at Y_2, and no second relation is
needed.  The price is that the principal preimage, at |z| = 0.0742, is by far the deepest
point to cut.  This script prices it.
"""
import numpy as np, math, contour as C, bcfast as BF, opt2, places
from region2 import Reg2
import itertools

m = 14; TAU = places.TAU

def run():
    nm1, sa1, Y1 = places.PLACES[0]
    nm2, sa2, Y2 = places.PLACES[1]
    # place v1 unchanged (its conditional function IS regular there, by the hypothesis)
    J1 = 91.17577
    print('  v1 (unchanged, hypothesis (i) holds there):  14L-BC = %.5f' % J1)
    print('  v2 with the principal preimage ALSO excluded:')
    r0 = C.preimages_h(Y2, rmax=0.999, cdmax=60)
    print('     principal preimage at |z| = %.6f  ->  hard cap |psi_2\'(0)| <= %.6f'
          % (abs(r0[0][0]), 4*abs(r0[0][0])/(1+abs(r0[0][0]))**2))
    best = None
    for R in (0.30, 0.40, 0.50, 0.60, 0.70, 0.78, 0.86):
        pre = [z for (z, cd, r) in r0 if abs(z) < R]          # principal INCLUDED
        n = len(pre)
        if n > 8: continue
        for c in (5.0, 7.5, 10.0):
            orders = list(itertools.permutations(range(n))) if math.factorial(max(n,1)) <= 5040 \
                     else [tuple(np.argsort([-abs(p) for p in pre]))]
            bo, br = None, -1.0
            for o in orders:
                rg = Reg2(R); rg.add_lune_at(0.0, c)
                for i in o: rg.add_slit_to(pre[i])
                if rg.rad() > br: br, bo = rg.rad(), o
            rg = Reg2(R); rg.add_lune_at(0.0, c)
            for i in bo: rg.add_slit_to(pre[i])
            lp = math.log(256*sa2*rg.rad())
            bc = BF.BC_fast(rg, lp, N=8192, cdmax=30, rmax=R*(1+1e-9))
            J = m*lp - bc
            print('     R=%.2f c=%.1f #cut=%d  |psi_2\'|=%.6f  L=%.6f  BC=%.5f  14L-BC=%.5f'
                  % (R, c, n, rg.rad(), lp, bc, J))
            if best is None or J > best[0]: best = (J, R, c, rg.rad(), lp, bc)
    print('     best 14L-BC at v2 = %.5f' % best[0])
    print()
    print('  MARGIN on the single hypothesis = (%.5f + %.5f)/2 - 14 tau = %+.5f'
          % (J1, best[0], (J1+best[0])/2 - 14*TAU))
    print('  (compare: with the doubled hypothesis, margin = -1.4416)')

if __name__ == '__main__':
    run()
