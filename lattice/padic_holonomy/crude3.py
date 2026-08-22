"""Crude (sup) multi-place bound for zeta_5(3), with a MAXIMUM-PRINCIPLE guard.

The crude form needs  S = sup_{Dbar} |x o psi|  (equivalently the H^infty norm),
NOT merely the sampled boundary values.  A contour internally tangent to |q|=1 at
an image of the cusp 0 has bounded a.e. boundary values but x o psi is UNBOUNDED
inside (singular-inner-function behaviour), so it is inadmissible.  We therefore
verify  max_boundary |x o psi|  >=  max over a grid of interior points, and reject
otherwise.
"""
import os, sys as _sys
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in _sys.path: _sys.path.insert(0, HERE)
import math, warnings, json
import numpy as np
warnings.filterwarnings('ignore')
from scipy.optimize import minimize
import haupt, family
from targets import TARGETS

T = TARGETS[5]; H, m, L, d = T['H'], T['m'], T['L'], T['d']
BUD = (m - 1) * L - m * d


def ev(specs, N=20000, guard=True):
    psi, dr = family.compose(specs)
    z = np.exp(2j * np.pi * (np.arange(N) + 0.5) / N)
    q = np.array([psi(zz) for zz in z], dtype=complex)
    if np.abs(q).max() > 1.0:
        return None
    g = haupt.logabs_x_vec(H, q)
    if not np.all(np.isfinite(g)):
        return None
    S = float(g.max())
    if guard:                      # maximum-principle check on interior points
        rr = np.array([0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999])
        tt = np.arange(600) / 600.
        zi = (rr[:, None] * np.exp(2j * np.pi * tt[None, :])).ravel()
        qi = np.array([psi(zz) for zz in zi], dtype=complex)
        gi = haupt.logabs_x_vec(H, qi)
        if gi.max() > S + 1e-6:
            return dict(reject=True, Sbd=S, Sint=float(gi.max()), logR=math.log(dr))
    return dict(reject=False, cost=max(S, 0.0) - m * math.log(dr), logR=math.log(dr),
                logS=S, maxq=float(np.abs(q).max()))


print("zeta_5(3) crude budget = (m-1)L - m*tau = %.7f    (need log^+S - 4 log R < this)" % BUD)
e = ev([('OFF', (0.0, 0.17, 0.83))])
print("\n[rejected example] D(0.17,0.83), tangent at q=+1 (an image of the CUSP 0):")
print("   boundary max log|x| = %.4f  but interior max = %.4f  ->  REJECTED (max principle)"
      % (e['Sbd'], e['Sint']))

best = (1e9, None)
def scan(tag, mk, grid, N=20000):
    global best
    loc = (1e9, None, None)
    nrej = 0
    for prm in grid:
        r = ev(mk(prm), N)
        if r is None: continue
        if r.get('reject'): nrej += 1; continue
        if r['cost'] < loc[0]: loc = (r['cost'], prm, r)
    if loc[1] is not None:
        print("   %-34s best %-22s logR=%+.6f logS=%+.6f cost=%9.6f  (%d rejected)"
              % (tag, str(np.round(np.atleast_1d(loc[1]), 4)), loc[2]['logR'], loc[2]['logS'], loc[0], nrej))
        if loc[0] < best[0]: best = (loc[0], tag, loc[1], loc[2])

print()
scan("centred disc |q|<=r", lambda r: [('SCALE', (r,))], np.arange(0.05, 0.75, 0.005))
scan("OFF D(-a,b), a+b<=0.999", lambda ab: [('OFF', (math.pi, ab[0], ab[1]))],
     [(a, b) for a in np.arange(0.01, 0.5, 0.01) for b in np.arange(0.05, 1.0, 0.01) if a < b <= 0.999 - a])
scan("OFF D(+a,b), a+b<=0.999", lambda ab: [('OFF', (0.0, ab[0], ab[1]))],
     [(a, b) for a in np.arange(0.01, 0.5, 0.01) for b in np.arange(0.05, 1.0, 0.01) if a < b <= 0.999 - a])
scan("TANG at zeta_5^j (cusp inf)", lambda p: [('TANG', (2 * math.pi * p[0] / 5., p[1]))],
     [(j, b) for j in (1, 2) for b in np.arange(0.505, 0.999, 0.005)])
scan("lune: SCALE(r) o BITE(0,c)", lambda p: [('SCALE', (p[0],)), ('BITE', (0.0, p[1]))],
     [(r, c) for r in np.arange(0.15, 1.0, 0.025) for c in (1.1, 1.2, 1.35, 1.5, 1.75, 2., 2.5, 3., 4., 6., 10., 20.)])
scan("gobble: r o BITE(0,c1) o BITE(pi,c2)",
     lambda p: [('SCALE', (p[0],)), ('BITE', (0.0, p[1])), ('BITE', (math.pi, p[2]))],
     [(r, c1, c2) for r in np.arange(0.2, 1.0, 0.05) for c1 in (1.2, 1.5, 2., 3., 5.) for c2 in (1.5, 2., 3., 6., 20.)])
print("\nBEST: cost = %.6f   vs budget %.6f   ->   %s"
      % (best[0], BUD, "CONTRADICTION" if best[0] < BUD else "FAILS"))
print("  attained by %s %s : logR=%+.6f (R=%.5f)  logS=%+.6f (S=%.5f)"
      % (best[1], np.round(np.atleast_1d(best[2]), 5), best[3]['logR'], math.exp(best[3]['logR']),
         best[3]['logS'], math.exp(best[3]['logS'])))
print("  crude bound = (log^+S + L)/(logR + L - tau) = %.6f  (m = 4)"
      % ((max(best[3]['logS'], 0) + L) / (best[3]['logR'] + L - d)))
json.dump(dict(cost=best[0], tag=best[1], prm=[float(x) for x in np.atleast_1d(best[2])],
               logR=best[3]['logR'], logS=best[3]['logS'], budget=BUD), open('crude_final.json', 'w'))
