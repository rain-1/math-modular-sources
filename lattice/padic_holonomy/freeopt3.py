"""Free-template optimisation allowing a NON-EVEN boundary profile
u(t) = -exp( a_0 + sum_k (a_k cos 2pi k t + b_k sin 2pi k t) ) - umin .
usage: python3 freeopt3.py <idx> <K> <N> <restarts> <tag>"""
import os, sys as _sys
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in _sys.path: _sys.path.insert(0, HERE)
import json, time, math, warnings
import numpy as np
warnings.filterwarnings('ignore')
from scipy.optimize import minimize
import haupt, outer
from targets import TARGETS
UMIN = 1e-3


def u_of(v, N, K):
    t = np.arange(N) / N
    w = np.full(N, v[0])
    for k in range(1, K + 1):
        w = w + v[2 * k - 1] * np.cos(2 * np.pi * k * t) + v[2 * k] * np.sin(2 * np.pi * k * t)
    return -np.exp(np.clip(w, -30, 30)) - UMIN


def cost(v, H, m, N, K):
    c = -1e9
    for NN in (N, 2 * N):
        u = u_of(v, NN, K)
        lg = outer.conj_fn(u)
        q = np.exp(2j * np.pi * np.arange(NN) / NN + lg)
        if np.abs(q).max() >= 1.0:
            return 1e7
        g = haupt.logabs_x_vec(H, q)
        if not np.all(np.isfinite(g)):
            return 1e7
        c = max(c, outer.rearr(g) - m * float(np.mean(u)))
    return c


idx = int(_sys.argv[1]); K = int(_sys.argv[2]); N = int(_sys.argv[3]); RS = int(_sys.argv[4]); tag = _sys.argv[5]
T = TARGETS[idx]; H, m = T['H'], T['m']
print("### %s NON-EVEN K=%d N=%d budget=%.6f" % (T['key'], K, N, T['budget']), flush=True)
rng = np.random.default_rng(4242 + K)
best = (1e9, None); t0 = time.time()
for s in range(RS):
    v0 = np.zeros(2 * K + 1); v0[0] = rng.normal(-0.5, 0.7); v0[1:] = rng.normal(0, 0.3, 2 * K)
    f = lambda v: cost(v, H, m, N, K)
    for meth in ('Powell', 'Nelder-Mead'):
        r = minimize(f, v0, method=meth, options=dict(maxfev=12000, xtol=1e-7, ftol=1e-10)
                     if meth == 'Powell' else dict(maxfev=12000, xatol=1e-7, fatol=1e-10))
        v0 = r.x
    if r.fun < best[0]: best = (r.fun, list(r.x))
    print("  restart %d cost=%.7f margin=%+.7f [%.0fs]" % (s, r.fun, T['budget'] - r.fun, time.time() - t0), flush=True)
json.dump(dict(key=T['key'], idx=idx, v=best[1], K=K, kind='noneven'), open(tag + '.json', 'w'), indent=1)
print("best %.7f  margin %+.7f" % (best[0], T['budget'] - best[0]), flush=True)
