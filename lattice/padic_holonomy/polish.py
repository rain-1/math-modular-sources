"""Local polish of a template directly in the polynomial coefficients c_0..c_K,
using the fast vectorised evaluator at high N (no aliasing possible).
usage: python3 polish.py <idx> <json with 'c' or 'a'> <K> <N> <iters> <tag>
"""
import os, sys as _sys
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in _sys.path: _sys.path.insert(0, HERE)
import json, math, time, warnings
import numpy as np
warnings.filterwarnings('ignore')
from scipy.optimize import minimize
import haupt, outer, certify
from targets import TARGETS


def make(c, N):
    t = np.arange(N) / N
    u = np.full(N, c[0])
    for k in range(1, len(c)):
        u = u + c[k] * np.cos(2 * np.pi * k * t)
    return u


def cost_c(c, H, m, N):
    u = make(c, N)
    if u.max() > -1e-6:
        return 1e6 + 1e3 * u.max()
    lg = outer.conj_fn(u)
    q = np.exp(2j * np.pi * np.arange(N) / N + lg)
    if np.abs(q).max() >= 1.0:
        return 1e7
    g = haupt.logabs_x_vec(H, q)
    if not np.all(np.isfinite(g)):
        return 1e7
    return outer.rearr(g) - m * c[0]


def main():
    idx = int(_sys.argv[1]); jf = _sys.argv[2]; K = int(_sys.argv[3])
    N = int(_sys.argv[4]); iters = int(_sys.argv[5]); tag = _sys.argv[6]
    T = TARGETS[idx]; H, m = T['H'], T['m']
    d = json.load(open(jf))
    if 'c' in d:
        c0 = np.array(d['c'], float)
    else:
        import freeopt
        u = freeopt.u_of(np.array(d['a']), 8192)
        c0 = certify.coeffs_from_u(u, K)
    c = np.zeros(K + 1)
    c[:min(len(c0), K + 1)] = c0[:K + 1]
    mx = certify.fine_max_u(c)
    if mx > -1e-9:
        c[0] -= mx + 1e-6
    print("### %s  start cost=%.9f (budget %.6f)  K=%d N=%d" %
          (T['key'], cost_c(c, H, m, N), T['budget'], K, N), flush=True)
    t0 = time.time()
    f = lambda v: cost_c(v, H, m, N)
    for it in range(iters):
        r = minimize(f, c, method='Powell', options=dict(maxfev=20000, xtol=1e-9, ftol=1e-12))
        c = r.x
        r = minimize(f, c, method='Nelder-Mead',
                     options=dict(maxfev=20000, xatol=1e-8, fatol=1e-11))
        c = r.x
        print("  iter %d cost=%.9f margin=%+.9f [%.0fs]" % (it, r.fun, T['budget'] - r.fun, time.time() - t0), flush=True)
    json.dump(dict(idx=idx, key=T['key'], c=list(map(float, c))), open(tag + '.json', 'w'), indent=1)
    print("done", flush=True)


if __name__ == '__main__':
    main()
