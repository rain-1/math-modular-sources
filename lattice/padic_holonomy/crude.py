"""The CRUDE (sup-norm) multi-place holonomicity bound of Dimitrov
(arXiv:1912.12545 / LNT lecture notes):

    r  <=  ( sum_v log^+ S_v ) / ( sum_v log R_v - tau_K(f) ),
    S_v = sup_{|z|_v = R_v} |x_v(z)|_v ,  x_v(z) = z + O(z^2), x_v(z)=z a.e. v.

Normalised to the unit disc: phi(z) = x_inf(R_inf z), so R_inf = |phi'(0)| and
S_inf = max_{|z|=1}|phi(z)|.  tau = 2k+1 (NOT the refined tau(b)).
Contradiction  <=>   crude_cost := log^+ S_inf - m*log R_inf  <  (m-1) L - m*tau .
"""
import os, sys as _sys
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in _sys.path: _sys.path.insert(0, HERE)
import math, warnings, json
import numpy as np
warnings.filterwarnings('ignore')
from scipy.optimize import minimize
import haupt, outer, family
from targets import TARGETS


def budget_crude(m, d, L):
    return (m - 1) * L - m * d


def cost_from_g(g, logrho, m):
    S = float(g.max())
    return max(S, 0.0) - m * logrho


def eval_specs(H, specs, m, N=8000):
    psi, dr = family.compose(specs)
    z = np.exp(2j * np.pi * (np.arange(N) + 0.5) / N)
    q = np.array([psi(zz) for zz in z], dtype=complex)
    if np.abs(q).max() >= 1.0:
        return None
    g = haupt.logabs_x_vec(H, q)
    if not np.all(np.isfinite(g)):
        return None
    return cost_from_g(g, math.log(dr), m), math.log(dr), float(g.max())


if __name__ == '__main__':
    # -------- calibration: Dimitrov's zeta_2(5) contour ---------------------
    X2 = haupt.X0p(2)
    L2, m2, d2 = 12 * math.log(2), 6, 5
    print("== calibration: zeta_2(5), tau=5, m=6, L=12log2 ==")
    print("   crude budget (m-1)L - m*tau = %.6f" % budget_crude(m2, d2, L2))
    # B = {|q + 3/16| <= 5/16}: OFF(theta=0 -> centre at -a) with a=3/16,b=5/16
    r = eval_specs(X2, [('OFF', (0.0, 3 / 16., 5 / 16.))], m2, 20000)
    print("   Dimitrov's B: logR=%.6f (log(1/5)=%.6f)  log S=%.6f (Dimitrov: log 3.2316=%.6f)"
          % (r[1], math.log(0.2), r[2], math.log(3.2316)))
    print("   crude bound = (logS + L)/(logR + L - tau) = %.6f   (Dimitrov: 5.58, exact 9.4909/1.7084=%.4f)"
          % ((max(r[2], 0) + L2) / (r[1] + L2 - d2), (math.log(3.2316) + L2) / (math.log(.2) + L2 - 5)))
    print("   crude cost = %.6f  vs budget %.6f  -> %s" % (r[0], budget_crude(m2, d2, L2),
          "CONTRADICTION" if r[0] < budget_crude(m2, d2, L2) else "fails"))

    # -------- zeta_5(3) -----------------------------------------------------
    print("\n== zeta_5(3): tau=3, m=4, L=3log5 ==")
    T = TARGETS[5]; H, m, L, d = T['H'], T['m'], T['L'], T['d']
    BUD = budget_crude(m, d, L)
    print("   crude budget (m-1)L - m*tau = 3*%.7f - 4*3 = %.7f" % (L, BUD))
    print("   need  log^+ S - 4 log R  <  %.7f" % BUD)
    best = (1e9, None)
    print("\n   (a) centred discs |q|<=r:")
    for rr in (0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6):
        e = eval_specs(H, [('SCALE', (rr,))], m)
        if e: print("       r=%.2f  logR=%+.5f  logS=%+.5f  cost=%9.5f" % (rr, e[1], e[2], e[0]))
    print("\n   (b) off-centre (Mobius) discs D(-a,b), Dimitrov's family:")
    for a in (0.05, 0.1, 0.15, 0.1875, 0.25, 0.3, 0.35, 0.4, 0.45):
        for b in (a + 0.05, a + 0.1, a + 0.15, a + 0.2, a + 0.25, a + 0.3, min(0.999, a + 0.4)):
            if b <= a or a + b >= 1.0: continue
            e = eval_specs(H, [('OFF', (0.0, a, b))], m)
            if e and e[0] < best[0]: best = (e[0], ('OFF', a, b), e)
    print("       best OFF: a=%.4f b=%.4f  logR=%+.6f logS=%+.6f  cost=%.6f"
          % (best[1][1], best[1][2], best[2][1], best[2][2], best[0]))
    # refine
    def f(v):
        a = 0.98 / (1 + math.exp(-v[0])); b = a + (0.999 - a) / (1 + math.exp(-v[1]))
        e = eval_specs(H, [('OFF', (0.0, a, b))], m, 8000)
        return e[0] if e else 1e6
    r0 = minimize(f, [math.log(best[1][1] / (0.98 - best[1][1])), 0.0], method='Nelder-Mead',
                  options=dict(maxfev=800, xatol=1e-6, fatol=1e-9))
    a = 0.98 / (1 + math.exp(-r0.x[0])); b = a + (0.999 - a) / (1 + math.exp(-r0.x[1]))
    e = eval_specs(H, [('OFF', (0.0, a, b))], m, 40000)
    print("       refined : a=%.6f b=%.6f  logR=%+.6f logS=%+.6f  cost=%.6f  bound=%.6f"
          % (a, b, e[1], e[2], e[0], (max(e[2], 0) + L) / (e[1] + L - d)))
    json.dump(dict(a=a, b=b, cost=e[0], logR=e[1], logS=e[2]), open('crude_off.json', 'w'))
