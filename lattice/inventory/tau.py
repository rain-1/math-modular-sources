"""tau^flat, tau^sharp, I_u^v(w): independent EXACT (Fraction) implementation of
CDT (6.0.4)-(6.0.6) and Definition 6.0.1.  Written from the formulas, not from
lattice/cdt_finder/cdt_bound.py (that file is used only as a cross-check).

Definition 6.0.1, for 0 <= max{u,1} <= v and w <= v:

 I_u^v(w) = int_{min(u,1)}^{1} max{t-w,0} dt
          + int_{max(u,1)}^{v} { sum_{h=1}^{floor((t-1)/max(1,w))} 1/h } dt
          + int_{max(u,1)}^{v} max{ t/floor((t+max(0,w-1))/max(1,w)) - w, 0 } dt .

Everything is computed exactly in Q when u,v,w are rational.
"""
from fractions import Fraction as F


def _seg(a, b, f):
    """int_a^b f, f exact on the cell; f returns a Fraction given (lo,hi)."""
    return f(a, b) if b > a else F(0)


def I_exact(u, v, w):
    u, v, w = F(u), F(v), F(w)
    mw = max(F(1), w)
    # term 1: int_{min(u,1)}^1 max(t-w,0) dt
    a, b = min(u, F(1)), F(1)
    t1 = F(0)
    if b > a:
        lo = max(a, w)
        if b > lo:
            t1 = (b - lo) ** 2 / 2
    A = max(u, F(1))
    # term 2: floor((t-1)/mw) = h  <=>  t in [1+h*mw, 1+(h+1)*mw)
    t2 = F(0); h = 0
    while True:
        c0, c1 = max(A, 1 + h * mw), min(v, 1 + (h + 1) * mw)
        if 1 + h * mw >= v:
            break
        if c1 > c0:
            H = sum(F(1, j) for j in range(1, h + 1))
            t2 += H * (c1 - c0)
        h += 1
    # term 3: floor((t+max(0,w-1))/mw) = k  <=>  t in [k*mw - s, (k+1)*mw - s)
    s = max(F(0), w - 1)
    t3 = F(0); k = 1
    while True:
        c0, c1 = max(A, k * mw - s), min(v, (k + 1) * mw - s)
        if k * mw - s >= v:
            break
        if c1 > c0:
            lo = max(c0, k * w)          # integrand t/k - w positive iff t > k w
            if c1 > lo:
                t3 += (c1 * c1 - lo * lo) / (2 * k) - w * (c1 - lo)
        k += 1
    return t1 + t2 + t3


def tau_flat(m, cols):
    """cols = [(u_j, b_j)]: sigma_m and tau^flat, exact."""
    sigma_m = sum(F(b) for (_, b) in cols)
    return sigma_m - sum(F(u) ** 2 * F(b) for (u, b) in cols) / F(m) ** 2


def tau_flat_direct(m, rows):
    """tau^flat from the definition (1/m^2) sum_i (2i-1) sigma_i, rows = list of
    sigma_i sorted nondecreasing.  Cross-check of the closed form."""
    assert len(rows) == m
    return sum(F(2 * i + 1) * F(rows[i]) for i in range(m)) / F(m) ** 2


def tau_sharp(m, e, grid=None):
    """(2/m^2) min_{xi in [0,m]} { xi*sum e_i + (max e_i) I_xi^m(xi) }.
    The objective is piecewise smooth; minimise exactly over a fine rational grid
    and then bisect.  Returns (value: float, xi: Fraction)."""
    se, me = sum(e), (max(e) if e else 0)
    if se == 0 or me == 0:
        return F(0), F(0)
    def obj(xi):
        return F(2, m * m) * (xi * se + me * I_exact(xi, m, xi))
    # coarse rational grid then local refinement (objective is convex-ish in xi;
    # a golden/ternary search is unsafe because of the floor kinks, so: grid+bisect)
    N = 2000 if grid is None else grid
    best, xb = None, None
    for i in range(N + 1):
        xi = F(m * i, N)
        v = obj(xi)
        if best is None or v < best:
            best, xb = v, xi
    step = F(m, N)
    for _ in range(40):
        step /= 2
        improved = False
        for xi in (xb - step, xb + step):
            if 0 <= xi <= m:
                v = obj(xi)
                if v < best:
                    best, xb, improved = v, xi, True
        if not improved and step < F(1, 10 ** 9):
            break
    return best, xb


def tau(m, cols, e):
    tf = tau_flat(m, cols)
    ts, xi = tau_sharp(m, e)
    return tf, ts, tf + ts, xi
