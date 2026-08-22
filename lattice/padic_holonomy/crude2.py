"""Crude (sup-norm) multi-place bound: calibration + optimisation for zeta_5(3)."""
import os, sys as _sys
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in _sys.path: _sys.path.insert(0, HERE)
import math, warnings, json, time
import numpy as np
warnings.filterwarnings('ignore')
from scipy.optimize import minimize
import haupt, outer, family
from targets import TARGETS


def ev(H, specs, m, N=20000):
    psi, dr = family.compose(specs)
    z = np.exp(2j * np.pi * (np.arange(N) + 0.5) / N)
    q = np.array([psi(zz) for zz in z], dtype=complex)
    if np.abs(q).max() > 1.0:
        return None
    g = haupt.logabs_x_vec(H, q)
    if not np.all(np.isfinite(g)):
        return None
    S = float(g.max()); R = math.log(dr)
    return dict(cost=max(S, 0.0) - m * R, logR=R, logS=S, maxq=float(np.abs(q).max()))


# ---------------- calibration: Dimitrov's zeta_2(5) ------------------------
X2 = haupt.X0p(2); L2, m2, d2 = 12 * math.log(2), 6, 5
print("== calibration zeta_2(5): tau=5 m=6 L=12log2=%.6f ; crude budget=%.6f" % (L2, 5 * L2 - 6 * 5))
r = ev(X2, [('OFF', (math.pi, 3 / 16., 5 / 16.))], m2, 40000)
print("   Dimitrov B={|q+3/16|<=5/16}: logR=%.6f (log 1/5=%.6f)  logS=%.6f  S=%.4f (his |x(1/8)|=3.2316)"
      % (r['logR'], math.log(.2), r['logS'], math.exp(r['logS'])))
print("   crude bound = %.6f   (his 9.5/1.7 < 5.58; exact %.4f)  cost=%.6f budget=%.6f -> %s"
      % ((max(r['logS'], 0) + L2) / (r['logR'] + L2 - d2),
         (math.log(3.2316) + L2) / (math.log(.2) + L2 - 5), r['cost'], 5 * L2 - 30,
         "CONTRADICTION" if r['cost'] < 5 * L2 - 30 else "fails"))

# ---------------- zeta_5(3) ------------------------------------------------
T = TARGETS[5]; H, m, L, d = T['H'], T['m'], T['L'], T['d']
BUD = (m - 1) * L - m * d
print("\n== zeta_5(3): tau=3 m=4 L=3log5=%.7f ; crude budget=(m-1)L-m*tau=%.7f" % (L, BUD))
best = (1e9, None)


def scan(tag, mk, grid):
    global best
    loc = (1e9, None)
    for prm in grid:
        e = ev(H, mk(prm), m)
        if e and e['cost'] < loc[0]:
            loc = (e['cost'], prm, e)
    if loc[1] is not None:
        print("   %-28s best %s : logR=%+.6f logS=%+.6f cost=%.6f" % (tag, loc[1], loc[2]['logR'], loc[2]['logS'], loc[0]))
        if loc[0] < best[0]:
            best = (loc[0], tag, loc[1], loc[2])


scan("centred disc r", lambda r: [('SCALE', (r,))], np.arange(0.05, 0.62, 0.01))
scan("OFF D(a,b) tangent at +1", lambda ab: [('OFF', (0.0, ab[0], 1 - ab[0]))],
     [(a,) * 1 for a in np.arange(0.02, 0.49, 0.01)])
scan("OFF D(-a,b) tangent at -1", lambda ab: [('OFF', (math.pi, ab[0], 1 - ab[0]))],
     [(a,) * 1 for a in np.arange(0.02, 0.49, 0.01)])
scan("TANG at zeta_5^j", lambda p: [('TANG', (2 * math.pi * p[0] / 5., p[1]))],
     [(j, b) for j in (1, 2) for b in np.arange(0.51, 0.99, 0.01)])
scan("OFF general (a,b) th=0", lambda ab: [('OFF', (0.0, ab[0], ab[1]))],
     [(a, b) for a in np.arange(0.02, 0.5, 0.02) for b in np.arange(0.05, 1.0, 0.02) if a < b <= 1 - a + 1e-12])
scan("lune r,c at +1", lambda p: [('SCALE', (p[0],)), ('BITE', (0.0, p[1]))],
     [(r, c) for r in np.arange(0.2, 1.0, 0.05) for c in (1.2, 1.5, 2., 2.5, 3., 4., 6., 10.)])
print("\n   BEST simple family: cost=%.6f  (%s %s)  budget=%.6f  ->  %s"
      % (best[0], best[1], best[2], BUD, "CONTRADICTION" if best[0] < BUD else "FAILS"))
json.dump(dict(cost=best[0], tag=best[1], prm=[float(x) for x in np.atleast_1d(best[2])],
               logR=best[3]['logR'], logS=best[3]['logS']), open('crude_best_simple.json', 'w'))
