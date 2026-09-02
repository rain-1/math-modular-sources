# -*- coding: utf-8 -*-
"""Cross-validation, attribution of the shortfall, hard caps, inventory optimisation."""
import numpy as np, math, json, contour as C, bcfast as BF, bcdirect as BD, opt2, places
from cdt_bound import tau_flat, tau_sharp

m = 14; TAU = places.TAU

def crossvalidate():
    print('=== BC: Jensen reduction vs the model-free double sum ===')
    print('  (a) concentric discs  phi = -h on |Q| < r  (host C of lattice/catalan_mu4)')
    class Disc:
        def __init__(s, R): s.R = float(R)
        def fwd(s, z): return s.R*np.asarray(z, dtype=complex)
        def inv(s, w): return np.asarray(w, dtype=complex)/s.R
    for R in (0.05, 0.2, 0.3, 0.5):
        d = Disc(R); lp = math.log(256*R)
        print('     r=%.2f: log|phi\'(0)|=%.6f  Jensen=%.6f  direct(2048)=%.6f'
              % (R, lp, BF.BC_fast(d, lp, N=4096, cdmax=25, rmax=R*1.0000001), BD.BC_direct(d, N=2048)))
    print('  (b) the place-v2 designs (deep slits, boundary near the h -> 0 cusps)')
    nm, sa, Y = places.PLACES[1]
    for R in (0.75, 0.83, 0.90):
        d, princ, pts = None, None, None
        r = opt2.evaluate(R, [(0.0, 7.5)], Y, sa, m, TAU, over=1.0, rho=0.999,
                          N=16384, cdmax=30, verify=False)
        reg = r['reg']
        print('     R=%.2f: log|phi\'(0)|=%.6f  Jensen=%.5f  direct 2048/4096 = %.5f %.5f'
              % (R, r['logphip'], r['BC'], BD.BC_direct(reg, N=2048, s=sa),
                 BD.BC_direct(reg, N=4096, s=sa)))

def caps():
    print('=== hard cap on |psi_v\'(0)| from the single deepest bad preimage ===')
    print('  (classical: among simply connected Omega in D with 0 in Omega and p not in Omega,')
    print('   the radial-slit domain is extremal, conformal radius 4|p|/(1+|p|)^2)')
    for nm, sa, Y in places.PLACES:
        pre = C.preimages_h(Y, rmax=0.99, cdmax=40)
        r0 = pre[1][2]
        cap = 4*r0/(1+r0)**2
        print('  %s: deepest bad preimage |p| = %.6f -> |psi\'(0)| <= %.6f, log(256|s|.cap) = %.6f'
              % (nm, r0, cap, math.log(256*sa*cap)))

def attribution():
    print('=== attribution of the shortfall ===')
    L1, B1 = 7.5999530, 15.223765          # certified, place v1
    L2, B2 = 2.2100944, 6.408897           # certified, place v2
    lp1 = L1 - math.log(256*abs(places.s_1)); lp2 = L2 - math.log(256*abs(places.s_2))
    sh1, sh2 = B1-L1, B2-L2
    print('  margin = 13 log256 + 13 * mean(log|psi_v\'(0)|) - mean(shape_v) - 14 tau')
    print('    13 log 256                 = %+.5f' % (13*math.log(256)))
    print('    13 * mean log|psi\'|        = %+.5f   (log|psi_1\'| = %.6f, log|psi_2\'| = %.6f)'
          % (13*(lp1+lp2)/2, lp1, lp2))
    print('    - mean shape               = %+.5f   (shape_1 = %.5f, shape_2 = %.5f)'
          % (-(sh1+sh2)/2, sh1, sh2))
    print('    - 14 tau                   = %+.5f' % (-14*TAU))
    tot = 13*math.log(256) + 13*(lp1+lp2)/2 - (sh1+sh2)/2 - 14*TAU
    print('    ------------------------------------')
    print('    margin                     = %+.5f' % tot)
    # CDT's own host, same accounting
    lpC, shC = math.log(0.6461355), 12.2007 - 5.107874 - math.log(0.6461355/0.6458)*0
    lpC = math.log(0.6461355); LC = math.log(256*0.6461355); shC = 12.2007 - LC
    print('  CDT\'s own host, same accounting (their optimum in this family, R = 0.80):')
    print('    13 log256 + 13 log|psi\'| - shape - 14 tau = %+.5f'
          % (13*math.log(256) + 13*lpC - shC - 14*TAU))
    print('  per-place score  s_v := 13 log|psi_v\'| - shape_v  (margin = 13 log256 + mean(s_v) - 14 tau):')
    print('    v1 : %+.5f      v2 : %+.5f      CDT (single place) : %+.5f'
          % (13*lp1-sh1, 13*lp2-sh2, 13*lpC-shC))
    print('    v1 is BETTER than CDT by %+.5f nats, v2 is WORSE by %+.5f nats;'
          % ((13*lp1-sh1)-(13*lpC-shC), (13*lpC-shC)-(13*lp2-sh2)))
    print('    halved and summed this is %+.5f nats of margin relative to CDT\'s own +%.4f.'
          % (((13*lp1-sh1)+(13*lp2-sh2))/2 - (13*lpC-shC),
             13*math.log(256) + 13*lpC - shC - 14*TAU))
    print()
    print('  what place v2 would have to deliver for margin 0 (v1 held at its optimum):')
    J1 = 14*L1 - B1
    need = 2*14*TAU - J1
    print('    need 14 L_2 - BC_2 >= %.4f ; achieved %.4f ; deficit %.4f' % (need, 14*L2-B2, need-(14*L2-B2)))
    print('    at fixed BC_2 that is |psi_2\'(0)| = %.6f (achieved %.6f, hard cap %.6f)'
          % (math.exp((need-(14*L2-B2))/14)*0.394939, 0.394939,
             (lambda r: 4*r/(1+r)**2)(0.213693)))

def inventory():
    print('=== inventory optimisation (u_1 = 1 forced by INVENTORY_BOUND Cor. 2.2) ===')
    C_ = json.load(open('curves.json'))
    def best_contour(mm):
        out = []
        for nm in ('v1', 'v2'):
            b = None
            for (R, nbad, psi, L, BC) in C_[nm]:
                v = mm*L - BC
                if b is None or v > b[0]: b = (v, R, L, BC)
            out.append(b)
        return (out[0][2]+out[1][2])/2, (out[0][3]+out[1][3])/2
    base = [0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1]     # CDT's fourteen
    print('  adding k two-layer functions of exponent 2 to CDT\'s fourteen:')
    print('   k    m       tau      entry     bound     margin')
    for k in range(0, 8):
        e = base + [2]*k; mm = len(e)
        sm, tf = tau_flat(mm, [(1, 2), (3, 2)]); ts, xi = tau_sharp(mm, e)
        tau = float(tf)+ts
        L, B = best_contour(mm); entry = L-tau
        if entry <= 0:
            print('  %3d %4d  %9.5f   ENTRY FAILS' % (k, mm, tau)); continue
        print('  %3d %4d  %9.5f  %8.5f  %8.4f  %+8.4f' % (k, mm, tau, entry, B/entry, mm*entry-B))
    print('  dropping conditional or pure functions (m < 14) only lowers m*entry:')
    for mm in (8, 10, 12, 13):
        e = base[:mm]
        sm, tf = tau_flat(mm, [(1, 2), (3, 2)]); ts, xi = tau_sharp(mm, e); tau = float(tf)+ts
        L, B = best_contour(mm); entry = L-tau
        print('   m=%2d tau=%.5f entry=%.5f bound=%.4f margin=%+.4f' % (mm, tau, entry, B/entry, mm*entry-B))

def convexity():
    print('=== CDT\'s convexity improvements, transported ===')
    print('  On their own host CDT improve the bound 13.9938 -> 13.730, 13.7206, 13.678, 13.621')
    print('  (four radii), i.e. a factor 13.621/13.9938 = %.6f on the bound.' % (13.621/13.9938))
    L, B = 4.9050237, 10.816331
    print('  Our bound 16.15428 * %.6f = %.4f, still > 14; equivalent margin %+.4f.'
          % (13.621/13.9938, 16.15428*13.621/13.9938, 14*(L-TAU) - B*13.621/13.9938))

if __name__ == '__main__':
    crossvalidate(); print()
    caps(); print()
    attribution(); print()
    inventory(); print()
    convexity()
