"""Free optimisation of the archimedean template over ALL admissible psi.

psi ranges over holomorphic self-maps of D fixing 0; the extremal ones for a
given boundary modulus are  psi(z) = z * (outer function with |G|=e^u).
u <= 0 is the only degree of freedom.  |x| is invariant under q -> conj(q), so
we restrict to even u (cosine series).
        u(t) = -exp( a_0 + sum_{k=1..K} a_k cos(2 pi k t) )
usage:  python3 freeopt.py <target-index> [K] [N] [restarts]
"""
import os, sys as _sys
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in _sys.path: _sys.path.insert(0, HERE)
import sys, json, time, math, warnings
import numpy as np
warnings.filterwarnings('ignore')
from scipy.optimize import minimize
import haupt, outer
from targets import TARGETS

import os
UMIN = -float(os.environ.get('UMIN', '1e-3'))   # |q| <= 1-UMIN : contour strictly interior


def u_of(a, N):
    t = np.arange(N) / N
    w = np.full(N, a[0])
    for k in range(1, len(a)):
        w = w + a[k] * np.cos(2 * np.pi * k * t)
    return -np.exp(np.clip(w, -30, 30)) + UMIN


def cost(a, H, m, N):
    u = u_of(a, N)
    r = outer.evaluate(H, u, False)
    if r is None:
        return 1e7
    return r['RE'] - m * r['logdr']


def main():
    idx = int(sys.argv[1]); K = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    N = int(sys.argv[3]) if len(sys.argv) > 3 else 1024
    RS = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    T = TARGETS[idx]; H, m = T['H'], T['m']
    print("### %s  m=%d budget=%.6f   K=%d N=%d" % (T['key'], m, T['budget'], K, N), flush=True)
    rng = np.random.default_rng(1234 + idx)
    t0 = time.time(); best = (1e9, None)
    for s in range(RS):
        a0 = np.zeros(K + 1)
        a0[0] = rng.normal(-0.5, 0.6)
        a0[1:] = rng.normal(0, 0.4, K)
        f = lambda a: cost(a, H, m, N)
        for meth, mx in (('Powell', 6000), ('Nelder-Mead', 6000)):
            res = minimize(f, a0, method=meth, options=dict(maxfev=mx, xtol=1e-5, ftol=1e-7)
                           if meth == 'Powell' else dict(maxfev=mx, xatol=1e-5, fatol=1e-7))
            a0 = res.x
        if res.fun < best[0]:
            best = (res.fun, list(res.x))
        print("  restart %d: cost=%.6f  margin=%+.6f  [%.0fs]" % (s, res.fun, T['budget'] - res.fun, time.time() - t0), flush=True)
    # final high-resolution evaluation
    a = np.array(best[1])
    outrec = dict(key=T['key'], m=m, d=T['d'], L=T['L'], tau=T['tau'], budget=T['budget'], a=list(a), K=K)
    for NN in (2048, 4096, 8192):
        u = u_of(a, NN)
        r = outer.evaluate(H, u, True)
        for w in ('RE', 'BC'):
            if w not in r: continue
            V = r[w]; lg = r['logdr']
            den = lg + T['L'] - T['tau']; num = V + T['L']
            print("  N=%5d %-2s  val=%10.6f logdr=%+9.6f cost=%10.6f margin=%+10.6f bound=%s maxq=%.6f gmax=%.3f" %
                  (NN, w, V, lg, V - m * lg, m * den - num,
                   ("%9.6f" % (num / den)) if den > 1e-12 else " n/a ", r['maxq'], r['gmax']), flush=True)
            outrec['N%d_%s' % (NN, w)] = dict(val=V, logdr=lg, cost=V - m * lg, margin=m * den - num,
                                              bound=(num / den if den > 1e-12 else None))
    json.dump(outrec, open(os.environ.get('TAG','free_%02d' % idx)+'.json', 'w'), indent=1)
    print("done %.0fs" % (time.time() - t0), flush=True)


if __name__ == '__main__':
    main()
