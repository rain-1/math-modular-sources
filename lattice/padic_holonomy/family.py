"""Composable conformal templates psi : D -> Omega subset closure(D), psi(0)=0.

Primitives (all send 0 -> 0, all have image inside the closed unit disc):
  SCALE(r)            z -> r z                                  |'|=r
  BITE(theta,c)       D -> D minus a lune-bite centred at e^{i theta}   |'|=(c^2-1)/(c^2+1)
                      (CDT L2chi App A.1: the two circles meet |z|=1 orthogonally)
  TANG(theta,b)       D -> the disc internally tangent to |z|=1 at e^{i theta},
                      of radius b (b in [1/2,1]);                |'|=2-1/b
  OFF(theta,a,b)      D -> D(-a e^{i(theta+pi)}, b) = off-centre disc, a+b<=1; |'|=(b^2-a^2)/b
"""
import os, sys as _sys
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in _sys.path: _sys.path.insert(0, HERE)
import numpy as np


def _h(z, c):
    c2 = c * c
    return (z * (1 + c2) - 1 - c2 + np.sqrt((1 + c2) ** 2 * (1 + z) ** 2 - 16 * c2 * z)) / (2 * (c2 - 1))


def prim(kind, prm):
    if kind == 'SCALE':
        r, = prm
        return (lambda z: r * z), abs(r)
    if kind == 'BITE':
        th, c = prm
        e = np.exp(1j * th)
        return (lambda z: -e * _h(-z / e, c)), (c * c - 1) / (c * c + 1)
    if kind == 'TANG':
        th, b = prm
        a = 1.0 - b
        e = np.exp(1j * th)
        return (lambda z: e * z * (b * b - a * a) / (b - a * z)), (b * b - a * a) / b
    if kind == 'OFF':
        th, a, b = prm
        e = np.exp(1j * th)
        return (lambda z: e * z * (b * b - a * a) / (b - a * z)), (b * b - a * a) / b
    raise ValueError(kind)


def compose(specs):
    fs, ds = [], []
    for kind, prm in specs:
        f, d = prim(kind, prm)
        fs.append(f)
        ds.append(d)

    def psi(z):
        w = z
        for f in reversed(fs):
            w = f(w)
        return w
    return psi, float(np.prod(ds))
