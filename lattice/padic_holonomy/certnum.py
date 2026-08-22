"""Rigorous-style numerics for a template psi(z) = z exp(c_0 + sum_k c_k z^k).

  * admissibility:   max_t u(t) <= 0  with u(t)=c_0+sum c_k cos(2 pi k t),
                     bounded via a fine grid + the exact Lipschitz constant
                     2 pi sum k|c_k| of u.
  * |phi'(0)| = |x'(0)| |psi'(0)| = e^{c_0}  EXACTLY (x = q + O(q^2)).
  * Lipschitz constant of g(t) = log|x(psi(e^{2 pi i t}))|:
        g'(t) = Re[ 2 pi i * Lam(q) * (1 + z P'(z)) ],  Lam = q d/dq log x,
    so  |g'| <= 2 pi * max|Lam| * (1 + sum k|c_k|),  and we also compute the
    true max |g'| on a fine grid.
  * quadrature: RE_N is a left-endpoint Riemann sum of the Lipschitz function
    max(g(s),g(t)) on the torus, hence  |RE - RE_N| <= Lip(g)/N.
"""
import os, sys as _sys
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in _sys.path: _sys.path.insert(0, HERE)
import json, math, warnings
import numpy as np
warnings.filterwarnings('ignore')
import haupt, outer
from targets import TARGETS


def contour(c, N):
    t = np.arange(N) / N
    u = np.full(N, c[0])
    for k in range(1, len(c)):
        u = u + c[k] * np.cos(2 * np.pi * k * t)
    lg = outer.conj_fn(u)
    q = np.exp(2j * np.pi * t + lg)
    return t, q, u


def report(idx, c, Ns=(2**12, 2**13, 2**14, 2**15, 2**16, 2**17, 2**18, 2**19, 2**20)):
    T = TARGETS[idx]; H, m, L, tau = T['H'], T['m'], T['L'], T['tau']
    c = np.asarray(c, float)
    S1 = float(np.abs(c[1:]).sum())
    Sk = float(np.sum(np.arange(1, len(c)) * np.abs(c[1:])))
    M = 1 << 22
    tt = np.arange(M) / M
    uu = np.full(M, c[0])
    for k in range(1, len(c)):
        uu += c[k] * np.cos(2 * np.pi * k * tt)
    umax_grid = float(uu.max())
    umax_cert = umax_grid + 2 * np.pi * Sk / (2 * M)
    print("admissibility: sum|c_k| = %.6f (c0+sum|c_k| = %+.6f);  max_t u <= %.6e   [grid 2^22 + Lip]" %
          (S1, c[0] + S1, umax_cert), flush=True)
    print("|phi'(0)| = e^{c0} = %.12f   log = %.12f  (exact)" % (math.exp(c[0]), c[0]), flush=True)

    # Lipschitz constant of g
    NL = 1 << 20
    t, q, u = contour(c, NL)
    rmax = float(np.abs(q).max())
    lam = haupt.theta_logx(H, q)
    z = np.exp(2j * np.pi * t)
    zPp = np.zeros(NL, dtype=complex)
    for k in range(1, len(c)):
        zPp += k * c[k] * z ** k
    gp = np.real(2j * np.pi * lam * (1 + zPp))
    Lip_true = float(np.abs(gp).max())
    Lip_crude = 2 * np.pi * float(np.abs(lam).max()) * (1 + Sk)
    print("max|q| on contour = %.9f ;  max|Lam| = %.4f ;  sum k|c_k| = %.4f" % (rmax, np.abs(lam).max(), Sk), flush=True)
    print("Lipschitz  Lip(g) : true max|g'| = %.4f    crude product bound = %.4f  [grid 2^20]" % (Lip_true, Lip_crude), flush=True)

    prev = None
    print(" %9s %18s %14s %14s %12s" % ("N", "RE_N", "diff", "Lip/N", "bound_N"), flush=True)
    rows = []
    for N in Ns:
        t, q, u = contour(c, N)
        g = haupt.logabs_x_vec(H, q)
        RE = outer.rearr(g)
        den = c[0] + L - tau
        bd = (RE + L) / den
        d = (RE - prev) if prev is not None else float('nan')
        prev = RE
        rows.append((N, RE, bd))
        print(" %9d %18.12f %14.3e %14.3e %12.9f" % (N, RE, d, Lip_true / N, bd), flush=True)
    N, RE, bd = rows[-1]
    eps = Lip_true / N
    den = c[0] + L - tau
    print("\nCERTIFIED (modulo the grid-sampled Lip and float64 arithmetic):")
    print("  RE  in [%.9f, %.9f]" % (RE - eps, RE + eps))
    print("  bound = (RE + L)/(c0 + L - tau)  <=  (%.9f + %.9f)/%.9f = %.9f"
          % (RE + eps, L, den, (RE + eps + L) / den))
    print("  margin = m(c0+L-tau) - (RE+L)   >=  %.9f" % (m * den - (RE + eps + L)))
    print("  m = %d ;  contradiction requires bound < %d" % (m, m), flush=True)
    return dict(RE=RE, eps=eps, bound=(RE + eps + L) / den, margin=m * den - (RE + eps + L),
                Lip=Lip_true, maxq=rmax, umax=umax_cert)


if __name__ == '__main__':
    idx = int(_sys.argv[1]); jf = _sys.argv[2]
    d = json.load(open(jf))
    c = d['c'] if 'c' in d else None
    if c is None:
        import certify, freeopt
        u = freeopt.u_of(np.array(d['a']), 8192)
        c = certify.coeffs_from_u(u, int(_sys.argv[3]) if len(_sys.argv) > 3 else 24)
        mx = certify.fine_max_u(c)
        if mx > 0:
            c[0] -= mx + 1e-12
    report(idx, c)
