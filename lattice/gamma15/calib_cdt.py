# -*- coding: utf-8 -*-
"""CALIBRATION on CDT's own L(2,chi_-3) proof (their Appendix A).

(1) Reproduce |psi'(0)| = 0.6292232680 from the exact product formula.
(2) Enumerate the h-preimages of -1/72 and check CDT's Lemma A.4.4 for their published
    parameters (it fails: six bad preimages, four slits, five of the six inside Omega).
(3) Compute BC of their published phi (12.1207, not 11.845).
(4) Rebuild a VALID contour of the same family (base radius R, one lune at the cusp z=1,
    one curved slit per bad preimage) and show that CDT's published pair
    (|psi'(0)|, BC) = (0.6292232680, 11.845) lies exactly on its curve at their own
    R = 77/100, c = 75/10, shrink = 995/1000.
"""
import numpy as np, math, cmath, contour as C, bcfast as BF, bcdirect as BD, opt2
from wind import winding

TARGET = -1/72.0
TAU = 16603/3920

class CDTmap:
    """CDT's published GGG, with its exact inverse."""
    R = float(C.CDT_R)
    def __init__(self):
        self.G = C.make_G(C.CDT_R, C.CDT_c, C.CDT_r, C.CDT_th, C.CDT_shrink)
        self.sh = float(C.CDT_shrink); self.c = float(C.CDT_c)
        self.rs = [float(t) for t in C.CDT_r]; self.ths = [float(t) for t in C.CDT_th]
    def fwd(self, z): return self.G(z)
    def inv(self, w):
        w = np.asarray(w, dtype=complex)
        u = C.lune_inv(-w/float(C.CDT_R), self.c); u = -u*np.exp(-2j*np.pi*self.ths[0])
        u = C.Slit_inv(u, self.rs[3])*np.exp(-2j*np.pi*self.ths[1])
        u = C.Slit_inv(u, self.rs[2])*np.exp(-2j*np.pi*self.ths[2])
        u = C.Slit_inv(u, self.rs[1])*np.exp(-2j*np.pi*self.ths[3])
        u = C.Slit_inv(u, self.rs[0]); return u/self.sh

if __name__ == '__main__':
    print('=== 1. CDT conformal radius, exactly ===')
    ex = C.G_rad(C.CDT_R, C.CDT_c, C.CDT_r, C.CDT_shrink)
    print('   |GGG\'(0)| =', ex)
    print('             = %.10f     CDT print: 0.6292232680' % float(ex))

    print()
    print('=== 2. h-preimages of -1/72 and CDT Lemma A.4.4 ===')
    pre = C.preimages_h(TARGET, rmax=0.9, cdmax=60)
    print('   principal preimage z0 = %.12g' % pre[0][0].real)
    for z, cd, r in pre[1:11]:
        print('     (2c,d)=%-9s |z|=%.6f  arg/2pi=%+.6f  h(z)+1/72 = %.2e'
              % (str(cd), r, cmath.phase(z)/(2*math.pi), abs(complex(C.h(z, 2000))+1/72.)))
    for R in (0.77,):
        N = 200000; z = R*np.exp(2j*np.pi*np.arange(N)/N)
        f = C.h(z, 1200) + 1/72.0
        a = np.angle(f); d = np.diff(np.concatenate([a, a[:1]])); d = (d+np.pi) % (2*np.pi)-np.pi
        print('   argument principle: #{h = -1/72 in |z| < %.2f} = %.4f' % (R, d.sum()/(2*np.pi)))
    m = CDTmap()
    bd = np.array(m.fwd(np.exp(2j*np.pi*np.arange(200000)/200000)))
    print('   winding numbers of CDT\'s published contour about the bad preimages:')
    for z, cd, r in pre[1:7]:
        print('     |z|=%.6f -> winding %+0.4f  (%s)' % (r, winding(bd, z),
              'INSIDE (violates Lemma A.4.4)' if abs(winding(bd, z)-1) < .3 else 'outside'))

    print()
    print('=== 3. BC of CDT\'s published phi ===')
    lp = math.log(256*float(ex))
    for N in (4096, 16384):
        v, sh, av, r0, r1 = BF.BC_fast(m, lp, N=N, cdmax=25, rmax=0.7700001, return_parts=True)
        print('   Jensen  N=%6d : BC = %.6f   (mean #extra preimages %.3f)' % (N, v, av))
    print('   direct  N= 4096 : BC = %.6f' % BD.BC_direct(m, N=4096))
    print('   CDT print: 11.845')

    print()
    print('=== 4. a VALID contour of the same family ===')
    print('   R      over   rho     |psi\'(0)|      BC       entry      bound     margin')
    for R, over, rho in ((0.77, 0.999, 0.995), (0.77, 1.0, 1.0), (0.80, 1.0, 1.0),
                         (0.83, 1.0, 1.0), (0.86, 1.0, 1.0)):
        r = opt2.evaluate(R, [(0.0, 7.5)], TARGET, 1.0, 14, TAU, over=over, rho=rho,
                          N=8192, cdmax=25)
        print('   %.2f   %.4f %.4f  %.7f  %8.4f  %8.6f  %8.4f  %+8.4f'
              % (R, over, rho, r['psi'], r['BC'], r['entry'], r['BC']/r['entry'],
                 14*r['entry']-r['BC']))
    print('   (CDT: |psi\'(0)| = 0.6292232680, BC = 11.845, bound 13.9938, margin +0.0053)')
