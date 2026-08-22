"""Fast free-template optimisation (vectorised log|x|), with anti-aliasing guard.

usage: python3 freeopt2.py <idx> <K> <N> <restarts> <umin> <tag>
The objective is evaluated at N AND at 2N and the WORSE (larger) cost is used, so
under-resolved spikes near |q|=1 can no longer be exploited.
"""
import os, sys as _sys
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in _sys.path: _sys.path.insert(0, HERE)
import json, time, math, warnings
import numpy as np
warnings.filterwarnings('ignore')
from scipy.optimize import minimize
import haupt, outer
from targets import TARGETS


def u_of(a, N, umin):
    t = np.arange(N) / N
    w = np.full(N, a[0])
    for k in range(1, len(a)):
        w = w + a[k] * np.cos(2 * np.pi * k * t)
    return -np.exp(np.clip(w, -30, 30)) - umin


def re_of(H, u):
    lg = outer.conj_fn(u)
    q = np.exp(2j * np.pi * np.arange(len(u)) / len(u) + lg)
    if np.abs(q).max() >= 1.0:
        return None, None
    g = haupt.logabs_x_vec(H, q)
    if not np.all(np.isfinite(g)):
        return None, None
    return outer.rearr(g), float(np.abs(q).max())


def cost(a, H, m, N, umin):
    c = -1e9
    for NN in (N, 2 * N):
        u = u_of(a, NN, umin)
        RE, mq = re_of(H, u)
        if RE is None:
            return 1e7
        c = max(c, RE - m * float(np.mean(u)))
    return c


def main():
    idx = int(_sys.argv[1]); K = int(_sys.argv[2]); N = int(_sys.argv[3])
    RS = int(_sys.argv[4]); umin = float(_sys.argv[5]); tag = _sys.argv[6]
    T = TARGETS[idx]; H, m = T['H'], T['m']
    print("### %s m=%d budget=%.6f K=%d N=%d umin=%g" % (T['key'], m, T['budget'], K, N, umin), flush=True)
    rng = np.random.default_rng(9000 + idx + 7 * K + int(1e4 * umin))
    t0 = time.time(); best = (1e9, None)
    for s in range(RS):
        a0 = np.zeros(K + 1)
        a0[0] = rng.normal(-0.5, 0.7)
        a0[1:] = rng.normal(0, 0.35, K)
        f = lambda a: cost(a, H, m, N, umin)
        for meth in ('Powell', 'Nelder-Mead', 'Powell'):
            r = minimize(f, a0, method=meth, options=dict(maxfev=9000, xtol=1e-6, ftol=1e-9)
                         if meth == 'Powell' else dict(maxfev=9000, xatol=1e-6, fatol=1e-9))
            a0 = r.x
        if r.fun < best[0]:
            best = (r.fun, list(r.x))
        print("  restart %d: cost=%.7f margin=%+.7f [%.0fs]" % (s, r.fun, T['budget'] - r.fun, time.time() - t0), flush=True)
    a = np.array(best[1])
    for NN in (8192, 16384, 65536, 262144):
        u = u_of(a, NN, umin)
        RE, mq = re_of(H, u)
        lg = float(np.mean(u))
        den = lg + T['L'] - T['tau']
        print("  N=%7d RE=%.9f logdr=%+.9f cost=%.9f margin=%+.9f bound=%.9f maxq=%.6f"
              % (NN, RE, lg, RE - m * lg, m * den - (RE + T['L']), (RE + T['L']) / den, mq), flush=True)
    json.dump(dict(key=T['key'], idx=idx, a=list(a), K=K, umin=umin), open(tag + '.json', 'w'), indent=1)
    print("done %.0fs" % (time.time() - t0), flush=True)


if __name__ == '__main__':
    main()
