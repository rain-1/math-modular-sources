# -*- coding: utf-8 -*-
"""The contours at the two real places of K = Q(sqrt5), and the averaged CDT bound
   m <= BCbar/(Lbar - tau)   [CDT Remark BCboundK; NUMBER_FIELD_HOLONOMY.md Thm 2.2].
"""
import numpy as np, math, cmath, json, contour as C, opt2, places
from wind import winding

m = 14
TAU = places.TAU

def geometry():
    print('=== the two places ===')
    print('  t1 = phi^-5 = %.15f   t2 = -phi^5 = %.15f   t1 t2 = %.15f'
          % (places.t1, places.t2, places.t1*places.t2))
    s5 = places.s5
    print('  place v1: s = t2, |s| = %.12f, extra point y1 = t1^2/(t1-t2) = %.15e,'
          % (abs(places.s_1), places.y_1))
    print('            Y1 = y1/s = %.15e   (exact: (682 sqrt5 - 1525)/25 = %.15e)'
          % (places.Y_1, (682*s5-1525)/25))
    print('  place v2: s = t1, |s| = %.12f, extra point y2 = t2^2/(t2-t1) = %.15f,'
          % (abs(places.s_2), places.y_2))
    print('            Y2 = y2/s = %.15f   (exact: -(682 sqrt5 + 1525)/25 = %.15f)'
          % (places.Y_2, -(682*s5+1525)/25))
    print('  4 s1 = %.6f   4 s2 = %.6f   |Y2|/4 = %.4f   log(256|s1|) = %.6f  log(256|s2|) = %.6f'
          % (4*places.s_1, 4*places.s_2, abs(places.Y_2)/4,
             math.log(256*abs(places.s_1)), math.log(256*abs(places.s_2))))
    for nm, sa, Y in places.PLACES:
        pre = C.preimages_h(Y, rmax=0.95, cdmax=40)
        print('  %s: principal h-preimage of Y is z0 = %.12g;  bad preimages:' % (nm, pre[0][0].real))
        for z, cd, r in pre[1:9]:
            print('      (2c,d)=%-9s |z| = %.6f   arg/2pi = %+.6f'
                  % (str(cd), r, cmath.phase(z)/(2*math.pi)))

def curves(fname='curves.json'):
    Rs = [0.40, 0.45, 0.50, 0.55, 0.60, 0.64, 0.68, 0.70, 0.72, 0.74, 0.76, 0.78,
          0.80, 0.82, 0.84, 0.86, 0.88, 0.90, 0.92, 0.94]
    out = {}
    for nm, sa, Y in places.PLACES:
        rows = []
        print('  --- %s ---' % nm)
        for R in Rs:
            r = opt2.evaluate(R, [(0.0, 7.5)], Y, sa, m, TAU, N=8192, cdmax=30, verify=False)
            rows.append((R, r['nbad'], r['psi'], r['logphip'], r['BC']))
            print('    R=%.2f nbad=%3d |psi\'|=%.6f  L=%.6f  BC=%.5f  14L-BC=%.4f'
                  % (R, r['nbad'], r['psi'], r['logphip'], r['BC'], r['J']))
        out[nm] = rows
    json.dump(out, open(fname, 'w'))
    return out

def optimum():
    print('=== optimised contours (scan over R and the lune parameter c) ===')
    best = {}
    for nm, sa, Y in places.PLACES:
        b = None
        for R in (0.72, 0.74, 0.76, 0.78, 0.80, 0.82):
            for c in (5.0, 6.5, 7.5, 9.0, 12.0):
                r = opt2.evaluate(R, [(0.0, c)], Y, sa, m, TAU, N=4096, cdmax=30, verify=False)
                if b is None or r['J'] > b['J']: b = r
        best[nm] = (b['R'], b['lunes'][0][1])
        print('  %s: best R = %.2f, c = %.1f  ->  |psi\'| = %.6f, L = %.6f, BC = %.5f, 14L-BC = %.5f'
              % (nm, b['R'], b['lunes'][0][1], b['psi'], b['logphip'], b['BC'], b['J']))
    return best

def final(best, over=0.9995, rho=0.999):
    print('=== final certified contours and the averaged bound ===')
    res = {}
    for tag, ov, rh in (('ideal', 1.0, 1.0), ('certified', over, rho)):
        vals = []
        for nm, sa, Y in places.PLACES:
            R, c = best[nm]
            r = opt2.evaluate(R, [(0.0, c)], Y, sa, m, TAU, over=ov, rho=rh,
                              N=16384, cdmax=30, verify=True)
            vals.append(r)
            print('  [%s] %s R=%.2f c=%.1f  nbad=%d nslit=%d  |psi\'(0)|=%.9f  L=%.7f  BC=%.6f  shape=%.6f'
                  % (tag, nm, R, c, r['nbad'], r['nslit'], r['psi'], r['logphip'], r['BC'], r['shape']))
            print('        all bad preimages excluded: %s (worst |Phi^-1(p)| = %.8f); principal inside: %s'
                  % (r['ok'], r['worst_excl'], r['principal_inside']))
            print('        slit radii r_j = %s' % [round(o[2], 6) for o in r['reg'].ops if o[0] == 's'])
        L = (vals[0]['logphip']+vals[1]['logphip'])/2
        B = (vals[0]['BC']+vals[1]['BC'])/2
        print('  [%s] Lbar = %.7f   BCbar = %.6f   tau = %.7f' % (tag, L, B, TAU))
        print('  [%s] entry = Lbar - tau = %+.7f      bound m <= %.5f      margin = 14*entry - BCbar = %+.5f'
              % (tag, L-TAU, B/(L-TAU), m*(L-TAU)-B))
        res[tag] = (L, B, vals)
    return res

if __name__ == '__main__':
    geometry(); print()
    print('=== trade-off curves L_v(R), BC_v(R) ==='); curves(); print()
    best = optimum(); print()
    final(best)
