"""Optimise the archimedean template phi = x o psi for the CDT holonomy bound.

Admissibility.  Omega must be simply connected, 0 in Omega, Omega subset of the
CLOSED unit q-disc, and psi(D) (open) inside the OPEN unit disc.  Boundary
contact with |q|=1 is allowed only at a MEASURE-ZERO set of points at which the
hauptmodul stays bounded, i.e. at cusps of the relevant Gamma_0(N) other than
the cusp 0 (where x has its pole).  For X_0(p) those are the Gamma_0(p)-images
of the cusp infinity, q = zeta_p^j (j=1..p-1) being the shallowest ones.
"""
import os, sys as _sys
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in _sys.path: _sys.path.insert(0, HERE)
import numpy as np, math, sys, json
from scipy.optimize import minimize
import haupt, bcint, family
from targets import TARGETS

# admissible tangency angles (arg q / 2pi), and the "depth" parameter c of the
# cusp a/c so that the tangent disc of radius b lands on the horocycle at height
# pi(1-b)/(c^2 b).
def good_angles(key):
    if key.startswith('zeta_2') or 'X_0(2)' in key:
        return [(0.5, 2)]                                    # q=-1  (cusp inf)
    if key.startswith('zeta_3'):
        return [(1/3., 3), (2/3., 3)]                        # q=zeta_3^{1,2}
    if key.startswith('zeta_5'):
        return [(j/5., 5) for j in (1, 2, 3, 4)]
    if key.startswith('zeta_7'):
        return [(j/7., 7) for j in range(1, 7)]
    if key.startswith('L_2'):                                 # X_1(4)
        return [(0.5, 2), (0.25, 4), (0.75, 4)]               # cusp 1/2 (x=-1/16); cusp inf
    if key.startswith('3-adic'):                              # X_0(9)
        return [(1/3., 3), (2/3., 3)] + [(j/9., 9) for j in (1, 2, 4, 5, 7, 8)]
    raise ValueError(key)


def evaluate(H, specs, N, which='RE', maxq_tol=1e-9):
    psi, dr = family.compose(specs)
    if not (dr > 1e-8):
        return None
    z = np.exp(2j * np.pi * (np.arange(N) + 0.5) / N)
    q = np.array([psi(zz) for zz in z], dtype=complex)
    if np.abs(q).max() > 1.0 - maxq_tol:
        return None
    try:
        ph = H.values(q)
    except Exception:
        return None
    g = np.log(np.abs(ph))
    if not np.all(np.isfinite(g)):
        return None
    out = dict(dr=dr, RE=bcint.rearrangement(g), gmax=float(g.max()), gmin=float(g.min()),
               maxq=float(np.abs(q).max()))
    if which in ('BC', 'both'):
        out['BC'] = bcint.bc_double_fast(z, ph)
    return out


SHAPES = {
    'S1': ('SCALE', 0),
    'S2': ('SCALE', 1),
    'S3': ('SCALE', 2),
    'S4': ('SCALE', 3),
    'O1': ('OFF', 0),
    'O2': ('OFF', 1),
    'O3': ('OFF', 2),
    'T1': ('TANG', 0),
    'T2': ('TANG', 1),
    'T3': ('TANG', 2),
    'T4': ('TANG', 3),
}


def sig(x):
    return 1.0 / (1.0 + math.exp(-max(-40., min(40., x))))


def make_specs(shape, v, tang_angle=None):
    base, nb = SHAPES[shape]
    specs = []
    if base == 'SCALE':
        specs.append(('SCALE', (0.02 + 0.979 * sig(v[0]),))); i = 1
    elif base == 'OFF':
        th = v[0]; a = 0.98 * sig(v[1]); b = a + (0.999 - 2 * a) * sig(v[2]) if a < 0.4995 else a + 1e-9
        specs.append(('OFF', (2 * math.pi * th, a, b))); i = 3
    else:
        b = 0.5 + 0.4999 * sig(v[0])
        specs.append(('TANG', (2 * math.pi * tang_angle, b))); i = 1
    for j in range(nb):
        th = v[i]; c = 1.0 + math.exp(max(-30., min(12., v[i + 1])))
        specs.append(('BITE', (2 * math.pi * th, c)))
        i += 2
    return specs


def nparam(shape):
    base, nb = SHAPES[shape]
    return (1 if base == 'SCALE' else 3 if base == 'OFF' else 1) + 2 * nb


def run(T, shape, which='RE', N=1200, starts=4, seed=0, tang_angle=None):
    H, m = T['H'], T['m']
    n = nparam(shape)
    rng = np.random.default_rng(seed)
    best = (1e9, None)
    def f(v):
        sp = make_specs(shape, list(v), tang_angle)
        r = evaluate(H, sp, N, which)
        if r is None:
            return 1e6
        return r[which] - m * math.log(r['dr'])
    for s in range(starts):
        v0 = rng.normal(0, 1.2, n)
        res = minimize(f, v0, method='Nelder-Mead',
                       options=dict(maxfev=250 + 130 * n, xatol=2e-3, fatol=2e-4))
        # polish
        res = minimize(f, res.x, method='Nelder-Mead',
                       options=dict(maxfev=200 + 80 * n, xatol=5e-4, fatol=5e-5))
        if res.fun < best[0]:
            best = (res.fun, make_specs(shape, list(res.x), tang_angle))
    return best
