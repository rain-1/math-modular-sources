"""Contour optimisation for the CDT holonomy bound on p-adic zeta values."""
import os, sys as _sys
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in _sys.path: _sys.path.insert(0, HERE)
import numpy as np, math, itertools
import haupt, bcint
from targets import TARGETS

NQ = 4000


def eval_contour(H, psi, dr, N=NQ, want_bc=True):
    t, z, q, ph = bcint.sample(psi, H, N)
    g = np.log(np.abs(ph))
    RE = bcint.rearrangement(g)
    BC = bcint.bc_double_fast(z, ph) if want_bc else None
    return dict(RE=RE, BC=BC, dr=dr, gmax=g.max(), gmin=g.min(), maxq=np.abs(q).max())


def report(T, res, use='RE'):
    V = res[use]
    lg = math.log(res['dr'])
    den = lg + T['L'] - T['tau']
    num = V + T['L']
    cost = V - T['m'] * lg
    out = dict(num=num, den=den, bound=(num / den if den > 0 else float('inf')),
               margin=T['m'] * den - num, cost=cost, budget=T['budget'], logdr=lg, val=V)
    return out


# ------------------------------------------------------------------ families
def build(kind, prm):
    """returns (psi, |psi'(0)|)"""
    if kind == 'disc':
        r, = prm
        return haupt.psi_disc(r)
    if kind == 'tang':          # tangent/off-centre disc, tangency direction angle th
        b, th = prm            # a = 1-b  (tangent);  om = -exp(i th) so tangency at exp(i th)
        a = 1.0 - b
        om = -np.exp(1j * th)
        return haupt.psi_offdisc(a, b, om)
    if kind == 'offd':
        a, b, th = prm
        om = -np.exp(1j * th)
        return haupt.psi_offdisc(a, b, om)
    if kind == 'lune':          # disc radius r with a bite pointing at angle th
        r, c, th = prm
        return haupt.psi_lune(r, c, np.exp(1j * th))
    if kind == 'gob':
        r, c1, c2, th = prm
        f, d = haupt.psi_gobble(r, c1, c2, np.exp(1j * th))
        return f, d
    raise ValueError(kind)


def compose(specs):
    fs, ds = [], []
    for kind, prm in specs:
        f, d = build(kind, prm)
        fs.append(f); ds.append(d)
    def psi(z):
        w = z
        for f in reversed(fs):
            w = f(w)
        return w
    return psi, float(np.prod(ds))
