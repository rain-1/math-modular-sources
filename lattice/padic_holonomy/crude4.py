"""Free-template optimisation of the CRUDE (sup-norm) bound for zeta_5(3).
u <= -UMIN keeps Omega compactly inside the unit q-disc, so x o psi is
holomorphic on the closed disc and sup over the boundary = sup over Dbar."""
import os, sys as _sys
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in _sys.path: _sys.path.insert(0, HERE)
import math, warnings, json, time
import numpy as np
warnings.filterwarnings('ignore')
from scipy.optimize import minimize
import haupt, outer
from targets import TARGETS
T = TARGETS[5]; H, m, L, d = T['H'], T['m'], T['L'], T['d']
BUD = (m - 1) * L - m * d
UMIN = float(_sys.argv[3]) if len(_sys.argv) > 3 else 0.005
K = int(_sys.argv[1]); N = int(_sys.argv[2])


def u_of(a, N):
    t = np.arange(N) / N
    w = np.full(N, a[0])
    for k in range(1, len(a)):
        w = w + a[k] * np.cos(2 * np.pi * k * t)
    return -np.exp(np.clip(w, -30, 30)) - UMIN


def cost(a, N):
    c = -1e9
    for NN in (N, 2 * N):
        u = u_of(a, NN)
        lg = outer.conj_fn(u)
        q = np.exp(2j * np.pi * np.arange(NN) / NN + lg)
        if np.abs(q).max() >= 1.0: return 1e7
        g = haupt.logabs_x_vec(H, q)
        if not np.all(np.isfinite(g)): return 1e7
        c = max(c, max(float(g.max()), 0.0) - m * float(np.mean(u)))
    return c


rng = np.random.default_rng(31337 + K)
best = (1e9, None); t0 = time.time()
print("crude budget = %.7f ; UMIN=%g K=%d N=%d" % (BUD, UMIN, K, N), flush=True)
for s in range(6):
    a0 = np.zeros(K + 1); a0[0] = rng.normal(-0.8, 0.6); a0[1:] = rng.normal(0, 0.3, K)
    f = lambda a: cost(a, N)
    for meth in ('Powell', 'Nelder-Mead', 'Powell'):
        r = minimize(f, a0, method=meth, options=dict(maxfev=6000, xtol=1e-7, ftol=1e-10)
                     if meth == 'Powell' else dict(maxfev=6000, xatol=1e-7, fatol=1e-10))
        a0 = r.x
    if r.fun < best[0]: best = (r.fun, list(r.x))
    print("  restart %d cost=%.6f  margin=%+.6f  [%.0fs]" % (s, r.fun, BUD - r.fun, time.time() - t0), flush=True)
a = np.array(best[1])
for NN in (20000, 80000, 320000):
    u = u_of(a, NN); lg = outer.conj_fn(u)
    q = np.exp(2j * np.pi * np.arange(NN) / NN + lg)
    g = haupt.logabs_x_vec(H, q)
    S = float(g.max()); R = float(np.mean(u))
    print("  N=%7d  logR=%+.6f logS=%+.6f  cost=%.6f  bound=%.6f  maxq=%.6f"
          % (NN, R, S, max(S, 0) - m * R, (max(S, 0) + L) / (R + L - d), float(np.abs(q).max())), flush=True)
json.dump(dict(a=list(a), cost=best[0], UMIN=UMIN, K=K), open('crude_free_K%d.json' % K, 'w'))
print("done", flush=True)
