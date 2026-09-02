# -*- coding: utf-8 -*-
"""Bost-Charles integral by an exact Jensen reduction.

For phi : (D,0) -> (C,0) analytic, continuous on the closed disc, with 0 its only zero
in D, Jensen's formula gives for every a
    int_T log|a - phi(w)| dmu(w) = log|a| + sum_{phi(w_k)=a, |w_k|<1} log(1/|w_k|),
and (a = 0)   int_T log|phi| dmu = log|phi'(0)|.  Hence
    BC(phi) = int int_{T^2} log|phi(z)-phi(w)| = log|phi'(0)| + int_T S(z) dmu(z),
    S(z) = sum_{w in D, phi(w) = phi(z)} log(1/|w|)  -- a pure multivalency term.
For phi = h o Phi with Phi : D -> Omega univalent and h the Gamma_0(2) Hauptmodul, the
interior solutions are w = Phi^{-1}(eta) for eta in the Gamma_0(2)-orbit of Phi(z)
lying inside Omega.

Numerical care: Phi^{-1} suffers catastrophic cancellation for tiny arguments (the
composed square roots collapse to u = i exactly), so for |eta| < `small` we use the
stable branch  Phi^{-1}(eta) = eta/Phi'(0) + O(eta^2)  with one Newton correction.
"""
import numpy as np, math, contour as C

def gamma_list(cdmax=14):
    out = []
    for c in range(0, cdmax+1):
        for d in range(-cdmax, cdmax+1):
            if c == 0 and d != 1: continue
            if c == 0 and d == 1: continue           # identity: skip
            if math.gcd(2*c, abs(d)) != 1: continue
            if not (c > 0 or (c == 0 and d > 0)): continue
            a, b = C._bezout(d, 2*c)
            out.append((a, -b, 2*c, d))
    return out

_GL = {}
def GL(cdmax):
    if cdmax not in _GL: _GL[cdmax] = gamma_list(cdmax)
    return _GL[cdmax]

def BC_fast(reg, logphip, N=4096, cdmax=14, rmax=None, tol=1e-6, return_parts=False,
            small=1e-4):
    if rmax is None: rmax = float(reg.R)*(1+1e-9)
    eps = 1e-9
    d0 = (complex(reg.fwd(np.array([eps]))[0]) - complex(reg.fwd(np.array([-eps]))[0]))/(2*eps)
    zs = np.exp(2j*np.pi*np.arange(N)/N)
    zeta = np.array(reg.fwd(zs))
    tau = np.log(zeta)/(2j*np.pi)
    tot = np.zeros(N); cnt = np.zeros(N); rej = [0, 0]
    for (a, b, c2, d) in GL(cdmax):
        t = (a*tau + b)/(c2*tau + d)
        ok = t.imag > 1e-14
        if not ok.any(): continue
        eta = np.where(ok, np.exp(2j*np.pi*t), 0.0)
        m = ok & (np.abs(eta) < rmax) & (np.abs(eta) > 1e-320)
        if not m.any(): continue
        w = np.zeros(N, dtype=complex)
        big = m & (np.abs(eta) >= small)
        sml = m & (np.abs(eta) < small)
        if big.any(): w[big] = reg.inv(eta[big])
        if sml.any():
            w0 = eta[sml]/d0
            w[sml] = w0 - (np.array(reg.fwd(w0)) - eta[sml])/d0
        g1 = m & (np.abs(w) < 1.0-1e-13) & (np.abs(w) > 1e-320)
        rej[0] += int((m & ~g1).sum())
        if not g1.any(): continue
        chk = np.zeros(N, dtype=complex); chk[g1] = reg.fwd(w[g1])
        good = g1 & (np.abs(chk-eta) < tol*np.maximum(1.0, np.abs(eta)))
        rej[1] += int((g1 & ~good).sum())
        tot[good] += -np.log(np.abs(w[good])); cnt[good] += 1
    shape = tot.mean()
    if return_parts: return logphip + shape, shape, cnt.mean(), rej[0]/N, rej[1]/N
    return logphip + shape
