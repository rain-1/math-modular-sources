"""Bost-Charles double integral and the (increasing-)rearrangement integral of
CDT, for phi = x o psi on the unit circle.

  BC(phi) = int_{T^2} log|phi(z)-phi(w)| dmu dmu        (Thm 7.0.1 numerator)
  RE(phi) = int_0^1 2t g^*(t) dt, g(t)=log|phi(e^{2 pi i t})|, g^* increasing
          = int_{[0,1]^2} max(g(s),g(t)) ds dt   >=  BC(phi)   (Nazarov)
"""
import os, sys as _sys
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in _sys.path: _sys.path.insert(0, HERE)
import numpy as np


def sample(psi, haupt, N):
    t = (np.arange(N) + 0.5) / N
    z = np.exp(2j * np.pi * t)
    q = np.array([psi(zz) for zz in z], dtype=complex)
    ph = haupt.values(q)
    return t, z, q, ph


def rearrangement(g):
    N = len(g)
    gs = np.sort(g)
    w = (2 * np.arange(1, N + 1) - 1) / float(N * N)
    return float(np.dot(w, gs))


def bc_double(z, ph):
    """int int log|phi(z)-phi(w)|; diagonal handled by log|z-w| subtraction
    (int int log|z-w| = 0 exactly on the unit circle)."""
    N = len(z)
    tot = 0.0
    for i in range(N):
        d = ph[i] - ph
        dz = z[i] - z
        with np.errstate(divide='ignore', invalid='ignore'):
            v = np.log(np.abs(d)) - np.log(np.abs(dz))
        v[i] = np.log(abs((ph[i] - ph[i - 1]) / (z[i] - z[i - 1]))) if N > 1 else 0.0
        tot += float(np.sum(v))
    return tot / (N * N)


def bc_double_fast(z, ph, blk=512):
    N = len(z)
    tot = 0.0
    for i0 in range(0, N, blk):
        i1 = min(N, i0 + blk)
        d = np.abs(ph[i0:i1, None] - ph[None, :])
        dz = np.abs(z[i0:i1, None] - z[None, :])
        np.maximum(d, 1e-300, out=d)
        with np.errstate(divide='ignore', invalid='ignore'):
            v = np.log(d) - np.log(dz)
        idx = np.arange(i0, i1)
        v[np.arange(i1 - i0), idx] = 0.0
        # diagonal: log|phi'(z_i)| approximated by neighbouring difference quotient
        dd = np.abs((ph[idx] - ph[idx - 1]) / (z[idx] - z[idx - 1]))
        v[np.arange(i1 - i0), idx] = np.log(np.maximum(dd, 1e-300))
        tot += float(np.sum(v))
    return tot / (N * N)
