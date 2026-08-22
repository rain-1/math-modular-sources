"""High-precision 2-parameter optimisation SCALE(r) o BITE(0,c) for X_0(5)."""
import os, sys as _sys
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in _sys.path: _sys.path.insert(0, HERE)
import numpy as np, math, warnings, itertools, sys, time
warnings.filterwarnings('ignore')
import haupt, family, outer
from targets import TARGETS

T = [t for t in TARGETS if t['key'].startswith('zeta_5(3)')][0]
H, m, L, tau, bud = T['H'], T['m'], T['L'], T['tau'], T['budget']


def ev(r, c, N, want_bc=True):
    psi, dr = family.compose([('SCALE', (r,)), ('BITE', (0.0, c))])
    z = np.exp(2j * np.pi * (np.arange(N) + 0.5) / N)
    q = np.array([psi(zz) for zz in z], dtype=complex)
    g = haupt.logabs_x(H, q)
    RE = outer.rearr(g)
    out = dict(RE=RE, dr=dr, gmax=g.max())
    if want_bc:
        ph = H.values_fast(q, 0.80)
        out['BC'] = outer.bc(z, ph)
    return out


def line(tag, r, c, N):
    e = ev(r, c, N)
    s = "%s r=%.6f c=%.6f N=%6d" % (tag, r, c, N)
    for w in ('RE', 'BC'):
        V = e[w]; lg = math.log(e['dr'])
        den = lg + L - tau; num = V + L
        s += "  %s=%.8f cost=%.8f margin=%+.8f bound=%.8f" % (w, V, V - m * lg, m * den - num, num / den)
    print(s, flush=True)
    return e


if __name__ == '__main__':
    N = 6000
    best = None
    for r in np.arange(0.56, 0.76, 0.02):
        for c in (2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 12.0):
            e = ev(r, c, N, False)
            cost = e['RE'] - m * math.log(e['dr'])
            if best is None or cost < best[0]:
                best = (cost, r, c)
    print('coarse best RE-cost=%.6f at r=%.3f c=%.2f  (budget %.6f)' % (best[0], best[1], best[2], bud), flush=True)
    r0, c0 = best[1], best[2]
    # refine
    for it in range(6):
        dr_, dc = 0.02 / 2 ** it, (c0 * 0.25) / 2 ** it
        cand = [(rr, cc) for rr in (r0 - dr_, r0, r0 + dr_) for cc in (max(1.05, c0 - dc), c0, c0 + dc)]
        vals = []
        for rr, cc in cand:
            e = ev(rr, cc, N, False)
            vals.append((e['RE'] - m * math.log(e['dr']), rr, cc))
        vals.sort()
        _, r0, c0 = vals[0]
    print('refined r=%.6f c=%.6f' % (r0, c0), flush=True)
    for N in (2000, 4000, 8000, 16000, 32000):
        line('opt', r0, c0, N)
