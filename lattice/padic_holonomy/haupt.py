"""Hauptmoduln of X_0(N) (N = p, p^2, 4) evaluated anywhere in the unit q-disc,
via Gamma_0(N)-reduction of tau, plus the archimedean contour templates of
Calegari-Dimitrov-Tang.

x_{X_0(p)}   = q prod ((1-q^{pn})/(1-q^n))^{24/(p-1)}     (p-1 | 24)
x_{X_0(p^2)} = q prod ((1-q^{p^2 n})/(1-q^n))^{24/(p^2-1)}
x_{X_1(4)}   = q prod ((1+q^n)(1+q^{2n}))^8   = (Delta(4tau)/Delta(tau))^{1/3}

All of these are Gamma_0(N)-invariant modular FUNCTIONS with integral q-expansion
q + O(q^2), so evaluation is branch-free.
"""
import numpy as np
from math import gcd, floor, ceil, pi

TWOPI = 2.0*pi


def _egcd(a, b):
    if b == 0:
        return (a, 1, 0)
    g, x, y = _egcd(b, a % b)
    return (g, y, x - (a // b) * y)


def _lattice_min(tau, N, span=4):
    """min |c*tau+d| over c in N*Z, d in Z, gcd(c,d)=1, via Gauss reduction of
    the rank-2 lattice Z*(N tau) + Z*1 followed by a small enumeration."""
    # basis vectors as (coeff_of_Ntau, coeff_of_1)
    B = [(0, 1), (1, 0)]

    def vec(b):
        return b[0] * N * tau + b[1]

    for _ in range(200):
        v0, v1 = vec(B[0]), vec(B[1])
        if abs(v0) > abs(v1):
            B[0], B[1] = B[1], B[0]
            v0, v1 = v1, v0
        if abs(v0) < 1e-300:
            break
        mu = int(round((v1.real * v0.real + v1.imag * v0.imag) / (abs(v0) ** 2)))
        if mu == 0:
            break
        B[1] = (B[1][0] - mu * B[0][0], B[1][1] - mu * B[0][1])
    best = None
    for i in range(-span, span + 1):
        for j in range(-span, span + 1):
            if i == 0 and j == 0:
                continue
            k = i * B[0][0] + j * B[1][0]
            d = i * B[0][1] + j * B[1][1]
            c = N * k
            if gcd(abs(c), abs(d)) != 1:
                continue
            v = abs(c * tau + d)
            if best is None or v < best[0]:
                best = (v, c, d)
    return best


def reduce_g0N(tau, N, K=None, iters=200):
    """Return gamma.tau for some gamma in Gamma_0(N), with Im maximised (greedily)."""
    tau = complex(tau)
    for _ in range(iters):
        tau = tau - round(tau.real)
        b = _lattice_min(tau, N)
        if b is None or b[0] >= 1.0 - 1e-14:
            return tau
        _, c, d = b
        g, u, v = _egcd(d, -c)
        assert abs(g) == 1, (c, d, g)
        if g < 0:
            u, v = -u, -v
        a, bb = u, v
        tau = (a * tau + bb) / (c * tau + d)
    return tau


class Haupt:
    """x(q) = q * prod_{n>=1} prod_j (1 - q^{m_j n})^{e_j}."""

    def __init__(self, N, prodspec, name):
        self.N = N
        self.prodspec = prodspec       # list of (m, e)
        self.name = name

    def q_series_value(self, q):
        """direct product evaluation; q complex scalar, |q| < 1."""
        aq = abs(q)
        if aq >= 1.0:
            raise ValueError("|q|>=1")
        nmax = int(40.0 / (-np.log(aq))) + 5
        val = complex(q)
        lg = 0.0 + 0.0j
        for (m, e) in self.prodspec:
            n = np.arange(1, nmax // m + 2)
            z = q ** (m * n)
            lg += e * np.sum(np.log1p(-z))
        return val * np.exp(lg)

    def value_tau(self, tau):
        t = reduce_g0N(tau, self.N)
        return self.q_series_value(np.exp(2j * pi * t))

    def value_q(self, q):
        q = complex(q)
        if abs(q) < 0.5:
            return self.q_series_value(q)
        tau = np.log(q) / (2j * pi)
        if tau.imag <= 0:
            tau = complex(tau.real, -tau.imag)
        return self.value_tau(tau)

    def values(self, qs):
        return np.array([self.value_q(q) for q in qs], dtype=complex)

    def values_fast(self, qs, cut=0.80):
        """vectorised for |q| <= cut, scalar Gamma_0(N)-reduction beyond."""
        qs = np.asarray(qs, dtype=complex)
        out = np.empty(len(qs), dtype=complex)
        a = np.abs(qs)
        lo = a <= cut
        if lo.any():
            ql = qs[lo]
            nmax = int(45.0 / (-np.log(max(a[lo].max(), 1e-12)))) + 5
            lg = np.zeros(len(ql), dtype=complex)
            for (m, e) in self.prodspec:
                n = np.arange(1, nmax // m + 2)
                z = ql[:, None] ** (m * n[None, :])
                lg += e * np.sum(np.log1p(-z), axis=1)
            out[lo] = ql * np.exp(lg)
        hi = ~lo
        if hi.any():
            out[hi] = [self.value_q(q) for q in qs[hi]]
        return out


def X0p(p):
    s = 24 // (p - 1)
    return Haupt(p, [(p, s), (1, -s)], "X_0(%d)" % p)


def X0p2(p):
    s = 24 // (p * p - 1)
    return Haupt(p * p, [(p * p, s), (1, -s)], "X_0(%d)" % (p * p))


def X14():
    # q prod (1+q^n)^8 (1+q^{2n})^8 = q prod (1-q^{2n})^8(1-q^n)^{-8}(1-q^{4n})^8(1-q^{2n})^{-8}
    #                              = q prod (1-q^n)^{-8} (1-q^{4n})^{8}
    return Haupt(4, [(4, 8), (1, -8)], "X_1(4)")


# ------------------------------------------------------------------ contours
def psi_disc(r):
    f = lambda z: r * z
    return f, r


def psi_offdisc(a, b, om=-1.0):
    """image = om * D(-a, b) ; tangent to |q|=1 at -om iff a+b=1."""
    f = lambda z: om * z * (b * b - a * a) / (b + a * z)
    return f, (b * b - a * a) / b


def h_lune(z, c):
    c2 = c * c
    num = z * (1 + c2) - 1 - c2 + np.sqrt((1 + c2) ** 2 * (1 + z) ** 2 - 16 * c2 * z)
    return num / (2 * (c2 - 1))


def psi_lune(r, c, om=1.0):
    """CDT's ps(z) = r*(-h(-z,c)) : D -> r * (D minus a bite near +1), rotated by om."""
    f = lambda z: om * (-r) * h_lune(-z, c)
    return f, r * (c * c - 1) / (c * c + 1)


def psi_gobble(r, e, f_, om=1.0):
    g = lambda z: om * (-r) * h_lune(-h_lune(z, f_), e)
    # |H'(0)| = ((e^2-1)/(e^2+1)) * ((f^2-1)/(f^2+1))
    d = r * (e * e - 1) / (e * e + 1) * (f_ * f_ - 1) / (f_ * f_ + 1)
    return g, d


# --------------------------------------------------------------------------
# Fully robust log|x| for ANY tau in H, via  log|Delta| and SL_2(Z) reduction.
#   x_N = (Delta(N tau)/Delta(tau))^{1/e},  e = 24/s   (s = the eta exponent)
#   log|x_N(tau)| = ( log|Delta(N tau)| - log|Delta(tau)| ) / e
# No branch ambiguity, no convergence problem anywhere in H.
# --------------------------------------------------------------------------
def _sl2_reduce(tau):
    """tau -> (tau_red in the standard fundamental domain, log|c tau + d|)."""
    tau = complex(tau)
    L = 0.0
    for _ in range(400):
        n = round(tau.real)
        if n:
            tau = tau - n
        a2 = abs(tau) ** 2
        if a2 >= 1.0 - 1e-15:
            break
        L += 0.5 * np.log(a2)      # |c tau + d| = |tau| for S
        tau = -1.0 / tau
    return tau, L


def log_absDelta(tau):
    t, L = _sl2_reduce(tau)
    q = np.exp(2j * pi * t)
    aq = abs(q)
    nmax = max(6, int(45.0 / max(-np.log(aq), 1e-12)) + 2)
    nmax = min(nmax, 200000)
    n = np.arange(1, nmax + 1)
    val = -TWOPI * t.imag + 24.0 * float(np.sum(np.log(np.abs(1.0 - q ** n))))
    return val - 12.0 * L


def make_logabs(N, e):
    def f(tau):
        return (log_absDelta(N * tau) - log_absDelta(tau)) / e
    return f


def logabs_x(H, qs):
    """log|x| for an array of q with 0<|q|<1, robust everywhere."""
    N = H.N
    e = 24 // [ee for (mm, ee) in H.prodspec if mm == N][0]
    out = np.empty(len(qs))
    for i, q in enumerate(qs):
        tau = np.log(complex(q)) / (2j * pi)
        if tau.imag <= 0:
            tau = complex(tau.real, -tau.imag)
        out[i] = (log_absDelta(N * tau) - log_absDelta(tau)) / e
    return out
