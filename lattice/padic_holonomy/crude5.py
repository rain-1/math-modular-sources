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
sg = lambda x: 1 / (1 + math.exp(-max(-40, min(40, x))))


def ev(specs, N=40000):
    psi, dr = family.compose(specs)
    z = np.exp(2j * np.pi * (np.arange(N) + 0.5) / N)
    q = np.array([psi(zz) for zz in z], dtype=complex)
    if np.abs(q).max() > 0.9999: return None
    g = haupt.logabs_x_vec(H, q)
    if not np.all(np.isfinite(g)): return None
    S = float(g.max())
    return max(S, 0.0) - m * math.log(dr), math.log(dr), S, float(np.abs(q).max())


def mk2(v):
    return [('SCALE', (0.02 + 0.977 * sg(v[0]),)), ('BITE', (0.0, 1 + math.exp(min(9, v[1]))))]
def mk3(v):
    return mk2(v) + [('BITE', (2 * math.pi * v[2], 1 + math.exp(min(9, v[3]))))]
def mk4(v):
    return mk3(v) + [('BITE', (2 * math.pi * v[4], 1 + math.exp(min(9, v[5]))))]

for nm, mk, n0 in (('2-par lune', mk2, 2), ('3-par (2 bites)', mk3, 4), ('4-par (3 bites)', mk4, 6)):
    rng = np.random.default_rng(7); best = (1e9, None)
    f = lambda v: (ev(mk(list(v)), 12000) or (1e6,))[0]
    for s in range(8):
        v0 = rng.normal(0, 1.0, n0); v0[0] = rng.normal(0.7, 0.6)
        r = minimize(f, v0, method='Nelder-Mead', options=dict(maxfev=3000, xatol=1e-7, fatol=1e-10))
        r = minimize(f, r.x, method='Powell', options=dict(maxfev=3000, xtol=1e-8, ftol=1e-11))
        if r.fun < best[0]: best = (r.fun, list(r.x))
    e = ev(mk(best[1]), 200000)
    sp = mk(best[1])
    print("%-16s cost=%.6f  logR=%+.6f (R=%.5f)  logS=%+.6f (S=%.5f)  maxq=%.5f  bound=%.6f"
          % (nm, e[0], e[1], math.exp(e[1]), e[2], math.exp(e[2]), e[3], (max(e[2], 0) + L) / (e[1] + L - d)), flush=True)
    print("    specs:", [(k, tuple(round(float(x), 6) for x in p)) for k, p in sp], flush=True)
    json.dump(dict(name=nm, cost=e[0], logR=e[1], logS=e[2],
                   specs=[[k, [float(x) for x in p]] for k, p in sp]), open('crude_%s.json' % nm.split()[0], 'w'))
print("budget = %.7f  -> all FAIL" % BUD)
