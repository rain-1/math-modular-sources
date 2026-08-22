"""EXACT parametrisation of all admissible archimedean templates.

The holonomy bound needs phi:(Dbar,0)->(C,0) holomorphic with phi^* f_i
meromorphic on D.  Here f_i = H, H', ..., H^{(2k)} and 1, all of which pull back
under the hauptmodul x to functions MEROMORPHIC on the whole q-disc |q|<1.
Hence phi = x o psi is admissible for ANY holomorphic psi : D -> D(0,1) with
psi(0)=0 -- psi need NOT be injective, and Omega=psi(D) need not be a Jordan
domain.

Every such psi is  psi(z) = z * G(z)  with G : D -> Dbar holomorphic, and for a
FIXED boundary modulus |G| = e^{u}, u <= 0, the choice maximising |psi'(0)| is
the OUTER function, for which

        log|psi'(0)| = log|G(0)| = mean(u)      (Jensen, equality iff G outer)

and the boundary values are  psi(e^{2 pi i t}) = exp(2 pi i t + u(t) + i*ut(t))
with ut = the conjugate function of u.  So the optimisation over ALL admissible
templates is exactly an optimisation over the single real function u <= 0.
"""
import os, sys as _sys
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in _sys.path: _sys.path.insert(0, HERE)
import numpy as np, math


def conj_fn(u):
    """u -> log G on the circle:  logG = uhat_0 + 2 sum_{k>=1} uhat_k e^{2 pi i k t}."""
    N = len(u)
    U = np.fft.fft(u) / N
    V = np.zeros(N, dtype=complex)
    V[0] = U[0]
    K = N // 2
    V[1:K] = 2 * U[1:K]
    if N % 2 == 0:
        V[K] = U[K]
    return np.fft.ifft(V) * N


def contour(u):
    """u (length N, <=0) -> boundary points q(t) of Omega and log|psi'(0)|."""
    N = len(u)
    t = np.arange(N) / N
    lg = conj_fn(u)
    q = np.exp(2j * np.pi * t + lg)
    return q, float(np.mean(u))


def rearr(g):
    N = len(g)
    gs = np.sort(g)
    w = (2 * np.arange(1, N + 1) - 1) / float(N * N)
    return float(np.dot(w, gs))


def bc(z, ph, blk=1024):
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
        dd = np.abs((ph[idx] - ph[idx - 1]) / (z[idx] - z[idx - 1]))
        v[np.arange(i1 - i0), idx] = np.log(np.maximum(dd, 1e-300))
        tot += float(np.sum(v))
    return tot / (N * N)


def evaluate(H, u, want_bc=False, cut=0.80):
    import haupt
    q, lgd = contour(u)
    if np.abs(q).max() >= 1.0:
        return None
    g = haupt.logabs_x(H, q)                  # robust everywhere in |q|<1
    if not np.all(np.isfinite(g)):
        return None
    out = dict(RE=rearr(g), logdr=lgd, gmax=float(g.max()), gmin=float(g.min()),
               maxq=float(np.abs(q).max()))
    if want_bc and out['gmax'] < 60.0:
        N = len(u)
        z = np.exp(2j * np.pi * np.arange(N) / N)
        ph = H.values_fast(q, cut)
        if np.all(np.isfinite(ph)):
            out['BC'] = bc(z, ph)
    return out
