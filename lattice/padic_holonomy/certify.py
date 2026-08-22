"""Turn an optimised boundary profile into an EXPLICIT admissible template

        psi(z) = z * exp( c_0 + c_1 z + ... + c_K z^K ),   c_k real

for which
  * psi is holomorphic on all of C, psi(0)=0, |psi'(0)| = e^{c_0};
  * log|psi(e^{2 pi i t})| = c_0 + sum_k c_k cos(2 pi k t)  =: u(t);
  * psi(D) subset D  as soon as  max_t u(t) <= 0   (checked on a fine grid, with
    a rigorous slack:  max_t u <= sum_k |c_k| + c_0 is a cheap sufficient test,
    and we also report the fine-grid max).
Then phi = x o psi is admissible for the CDT holonomy bound.
usage: python3 certify.py <idx> <json-with-'a'> [K] 
"""
import os, sys as _sys
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in _sys.path: _sys.path.insert(0, HERE)
import sys, json, math, warnings
import numpy as np
warnings.filterwarnings('ignore')
import haupt, outer
from targets import TARGETS


def coeffs_from_u(u, K):
    N = len(u)
    U = np.fft.fft(u) / N
    c = np.zeros(K + 1)
    c[0] = U[0].real
    for k in range(1, K + 1):
        c[k] = 2 * U[k].real
    return c


def u_from_c(c, N):
    t = np.arange(N) / N
    u = np.full(N, c[0])
    for k in range(1, len(c)):
        u = u + c[k] * np.cos(2 * np.pi * k * t)
    return u


def evaluate_c(T, c, N, want_bc=True):
    u = u_from_c(c, N)
    r = outer.evaluate(T['H'], u, want_bc)
    if r is None:
        return None
    V = {}
    for w in ('RE', 'BC'):
        if w not in r:
            continue
        val = r[w]; lg = c[0]
        den = lg + T['L'] - T['tau']; num = val + T['L']
        V[w] = dict(val=val, logdr=lg, cost=val - T['m'] * lg,
                    margin=T['m'] * den - num, bound=(num / den if den > 1e-12 else None))
    V['maxq'] = r['maxq']; V['gmax'] = r['gmax']; V['gmin'] = r['gmin']
    return V


def fine_max_u(c, M=200000):
    t = np.arange(M) / M
    u = np.full(M, c[0])
    for k in range(1, len(c)):
        u = u + c[k] * np.cos(2 * np.pi * k * t)
    return float(u.max())


if __name__ == '__main__':
    idx = int(sys.argv[1]); jf = sys.argv[2]; K = int(sys.argv[3]) if len(sys.argv) > 3 else 48
    T = TARGETS[idx]
    d = json.load(open(jf))
    a = np.array(d['a'])
    import freeopt
    u = freeopt.u_of(a, 8192)
    c = coeffs_from_u(u, K)
    mx = fine_max_u(c)
    if mx > 0:
        c[0] -= mx + 1e-12
        mx = fine_max_u(c)
    print("K=%d  c0=%.10f  max_t u = %.3e  sum|c_k|(k>=1)=%.6f" % (K, c[0], mx, np.abs(c[1:]).sum()), flush=True)
    print("c =", ", ".join("%.8f" % v for v in c), flush=True)
    for N in (4096, 8192, 16384, 32768):
        V = evaluate_c(T, c, N)
        s = "N=%6d maxq=%.6f gmax=%8.4f gmin=%9.4f" % (N, V['maxq'], V['gmax'], V['gmin'])
        for w in ('RE', 'BC'):
            if w in V:
                s += "  %s: val=%.7f cost=%.7f margin=%+.7f bound=%.7f" % (
                    w, V[w]['val'], V[w]['cost'], V[w]['margin'], V[w]['bound'])
        print(s, flush=True)
    json.dump(dict(idx=idx, key=T['key'], c=list(c)), open('cert_%02d.json' % idx, 'w'), indent=1)
