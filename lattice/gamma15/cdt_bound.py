"""Exact implementation of the Calegari-Dimitrov-Tang arithmetic holonomy bound.

Source: F. Calegari, V. Dimitrov, Y. Tang, "The linear independence of 1, zeta(2)
and L(2,chi_{-3})", arXiv:2408.15403v2.  All equation numbers below are theirs.

Theorem 6.0.2 / Theorem 7.0.1.  Data:
  * m functions f_1..f_m in Q[[x]], Q(x)-linearly independent, holonomic;
  * an m x r denominator array b, each column of the "step" shape
        0 = b_{1,j} = ... = b_{u_j,j} < b_{u_j+1,j} = ... = b_{m,j} =: b_j ,
    meaning  f_i(x) = a_{i,0} + sum_n a_{i,n} x^n / ( n^{e_i} prod_j [1..b_{i,j} n] ),
    a_{i,n} in Z;
  * an integration vector e in N^m;
  * a holomorphic phi:(D,0)->(C,0) with log|phi'(0)| > tau(b;e) and phi^* f_i
    meromorphic on D.
Then                    BC-integral(phi)
              m  <=  ---------------------------          (7.0.1 = BCbound)
                      log|phi'(0)| - tau(b;e)
with  tau(b;e) = tau_flat(b) + tau_sharp(e),
  tau_flat(b) = (1/m^2) sum_i (2i-1) sigma_i = sigma_m - (1/m^2) sum_j u_j^2 b_j   (6.0.4)
  tau_sharp(e) = (2/m^2) min_{xi in [0,m]} { xi*sum_i e_i + (max_i e_i) I_xi^m(xi) } (6.0.5)
and I_u^v(w) as in Definition 6.0.1.
"""
from fractions import Fraction
import math


# ---------------------------------------------------------------- Def 6.0.1
def I_uvw(u, v, w, nsub=200000):
    """I_u^v(w) of Definition 6.0.1.  Piecewise-exact where the floors are
    constant; the three terms are integrated exactly on each constancy cell."""
    u = float(u); v = float(v); w = float(w)
    mw = max(1.0, w)
    lo1 = min(u, 1.0)
    # term 1: int_{min(u,1)}^{1} max(t-w,0) dt
    a, b = lo1, 1.0
    t1 = 0.0
    if b > a:
        lo = max(a, w)
        if b > lo:
            t1 = 0.5*(b-lo)**2
    A = max(u, 1.0)
    # cells for term 2: floor((t-1)/mw) = h  <=>  t in [1+h*mw, 1+(h+1)*mw)
    t2 = 0.0
    h = 0
    while True:
        c0 = max(A, 1.0 + h*mw); c1 = min(v, 1.0 + (h+1)*mw)
        if c0 >= v: break
        if c1 > c0:
            H = sum(1.0/j for j in range(1, h+1))     # empty sum = 0 for h=0
            t2 += H*(c1-c0)
        h += 1
        if h > 10**6: break
    # cells for term 3: floor((t+max(0,w-1))/mw) = k <=> t in [k*mw - s, (k+1)*mw - s)
    s = max(0.0, w-1.0)
    t3 = 0.0
    k = 1
    while True:
        c0 = max(A, k*mw - s); c1 = min(v, (k+1)*mw - s)
        if k*mw - s >= v: break
        if c1 > c0:
            # integrand max(t/k - w, 0), positive for t > k*w
            lo = max(c0, k*w)
            if c1 > lo:
                t3 += (c1*c1 - lo*lo)/(2.0*k) - w*(c1-lo)
        k += 1
        if k > 10**6: break
    return t1 + t2 + t3


# ---------------------------------------------------------------- 6.0.4/6.0.5
def tau_flat(m, cols):
    """cols = list of (u_j, b_j).  Returns sigma_m and tau^flat as Fractions."""
    sigma_m = sum(Fraction(b) for (_, b) in cols)
    corr = sum(Fraction(u)**2*Fraction(b) for (u, b) in cols)/Fraction(m)**2
    return sigma_m, sigma_m - corr


def tau_sharp(m, e, ngrid=4001):
    se = sum(e); me = max(e) if e else 0
    if se == 0 and me == 0:
        return 0.0, None
    best, argbest = None, None
    for i in range(ngrid):
        xi = m*i/(ngrid-1)
        val = (2.0/m**2)*(xi*se + me*I_uvw(xi, m, xi))
        if best is None or val < best:
            best, argbest = val, xi
    return best, argbest


def tau(m, cols, e):
    sm, tf = tau_flat(m, cols)
    ts, xi = tau_sharp(m, e)
    return dict(sigma_m=sm, tau_flat=tf, tau_sharp=ts, xi=xi, tau=float(tf)+ts)


def bound(m, cols, e, log_phiprime, numerator):
    """CDT Theorem 7.0.1: m <= numerator/(log|phi'(0)| - tau).  Returns
    (bound, margin) with margin = m*(log|phi'(0)|-tau) - numerator > 0 meaning
    the bound is violated, i.e. a contradiction, i.e. the theorem follows."""
    T = tau(m, cols, e)
    den = log_phiprime - T['tau']
    return dict(tau=T, denom=den, bound=(numerator/den if den > 0 else float('inf')),
                margin=m*den - numerator, entry=den)


# ---------------------------------------------------------------- calibration
if __name__ == '__main__':
    print('=== CDT level-6 weight-3 calibration (their Section 13) ===')
    m = 14
    cols = [(1, 2), (3, 2)]                       # (13.0.2): u_1=1, u_2=3, b=2,2
    e = [0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]  # (13.0.2)
    sm, tf = tau_flat(m, cols)
    print(f'  sigma_m      = {sm}          (CDT: 4)')
    print(f'  tau^flat(b)  = {tf} = {float(tf):.6f}   (CDT: 191/49 = {191/49:.6f})')
    print(f'  I_2^14(2)    = {I_uvw(2,14,2):.6f}      (must be 21.075 exactly)')
    ts, xi = tau_sharp(m, e)
    print(f'  tau^#(e)     = {ts:.10f} at xi={xi:.4f}  (CDT: 27/80 = {27/80}, minimiser [2,13/6])')
    T = float(tf)+ts
    print(f'  tau(b;e)     = {T:.6f}                    (CDT: 16603/3920 = {16603/3920:.6f})')
    psi = Fraction(5448339453535586608000000000, 8658833407565631122430056127)
    lp = math.log(256*float(psi))
    print(f'  |psi\'(0)|    = {float(psi):.10f}         (CDT: 0.6292232680...)')
    print(f'  |phi\'(0)|    = {256*float(psi):.6f}      (= 256*|psi\'(0)|)')
    print(f'  log|phi\'(0)| = {lp:.6f}')
    for nm, num in (('Thm 7.0.1 (Bost-Charles)', 11.845),):
        r = bound(m, cols, e, lp, num)
        print(f'  {nm}: denom={r["denom"]:.6f}, m <= {r["bound"]:.4f}   (CDT: 13.9938)')
        print(f'    signed margin  m*(log|phi\'(0)|-tau) - BC = {r["margin"]:.4f}')
