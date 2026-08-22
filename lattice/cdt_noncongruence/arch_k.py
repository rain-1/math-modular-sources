"""ARCHITECTURE K: the rigorous univalent (Koebe) form of the CDT bound for a
second-order Apery row.

Hypothesis to be contradicted:  xi := lim b_n/a_n = p/q in Q.
Conditional function:           H := q B - p A  in Q[[x]],  d_n^2 H_n in Z.
  H is holomorphic at the INNER singular point x_1 = 1/lam_1 (that is exactly
  |a_n xi - b_n| ~ |lam_2|^n < |lam_1|^{-n}... i.e. radius of convergence
  |x_2| > |x_1|), and its only finite branch point is x_2 = 1/lam_2.

Domain.  Omega = C \ (the ray from x_2 away from 0).  Simply connected, contains
0 and x_1; every f_i below is single-valued and holomorphic on it.  By Koebe's
1/4 theorem the conformal radius is EXACTLY rho = 4|x_2| = 4/|lam_2| and this is
the largest possible for any simply connected Omega omitting x_2.
For a univalent phi, CDT's Bost-Charles numerator is exactly BC = log rho
(Grunsky: log((phi(z)-phi(w))/(z-w)) has constant term log phi'(0)).

Inventory (all conditional, all with rational coefficients under the hypothesis):
   f_1 = 1                                sigma = 0,  e = 0
   f_2 = H                                sigma = 2,  e = 0     (two [1..n] layers)
   f_3 = theta H = x H'                   sigma = 2,  e = 0
   f_{3+j} = theta^{-j}(H - H(0)), j>=1   sigma = 2,  e = j
The second-order inhomogeneous ODE satisfied by H gives a Q(x)-relation among
{1, H, theta H, theta^2 H}, so theta^2 H and beyond are NOT admissible; the tower
must be continued downwards by Eichler integrals, which cost e_i = j.

Optionally (only when lam_2 in Z, so that the polylogarithms have rational
integer coefficients) a PURE sub-inventory
   1, Li_1(lam_2 x), ..., Li_{p-1}(lam_2 x)   sigma = 0,  e = 0,1,...,p-1
which raises u_1 = u_2 from 1 to p.

Bound: m <= BC/(log|phi'(0)| - tau) = log rho/(log rho - tau).
Violation  <=>  m(log rho - tau) > log rho  <=>  log rho > tau * m/(m-1).
"""
import math, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'cdt_finder'))
from cdt_bound import tau_flat, tau_sharp

def inventory(J=0, p=1):
    """p pure functions (sigma=0, e=0..p-1) incl. the constant, plus the
    conditional tower H, theta H, theta^{-1}H, ..., theta^{-J}H."""
    sig = [0]*p + [2]*(2+J)
    e   = list(range(p)) + [0, 0] + list(range(1, J+1))
    m = len(sig)
    cols = [(p, 1), (p, 1)]                   # u_j = p, b_j = 1 (two [1..n] layers)
    return m, cols, e

def tau_of(m, cols, e):
    sm, tf = tau_flat(m, cols)
    ts, xi = tau_sharp(m, e)
    return float(tf), ts, float(tf)+ts

def arch_k(log_rho, J=0, p=1):
    m, cols, e = inventory(J, p)
    tf, ts, T = tau_of(m, cols, e)
    entry = log_rho - T
    return dict(m=m, u=p, J=J, tau_flat=tf, tau_sharp=ts, tau=T, entry=entry,
                bound=(log_rho/entry if entry > 0 else float('inf')),
                margin=m*entry - log_rho,
                threshold=T*m/(m-1))          # log rho needed for a contradiction

def best_arch_k(log_rho, Jmax=12, pmax=1):
    best = None
    for p in range(1, pmax+1):
        for J in range(0, Jmax+1):
            r = arch_k(log_rho, J, p)
            if best is None or r['margin'] > best['margin']:
                best = r
    return best

if __name__ == '__main__':
    print(__doc__)
    print("="*100)
    print("Thresholds:  a contradiction at inventory (p pure, J Eichler integrals) needs")
    print("             log rho = log(4/|lam_2|) > tau*m/(m-1),  i.e. score > threshold - log 4 - 2\n")
    print(f"{'p':>3} {'J':>3} {'m':>4} {'tau^flat':>9} {'tau^#':>8} {'tau':>8} "
          f"{'thresh log rho':>15} {'thresh score':>13}")
    for p in (1,):
        for J in range(0, 9):
            r = arch_k(0.0, J, p)
            print(f"{p:>3} {J:>3} {r['m']:>4} {r['tau_flat']:>9.5f} {r['tau_sharp']:>8.5f} "
                  f"{r['tau']:>8.5f} {r['threshold']:>15.5f} {r['threshold']-math.log(4)-2:>+13.5f}")
    print("\n  (p = 1: only the constant is denominator-free -- the generic case, and the")
    print("   only one available when lam_2 is not a rational integer.)")
    print("\nWith a pure polylogarithm inventory (needs lam_2 in Z):")
    print(f"{'p':>3} {'J':>3} {'m':>4} {'tau':>8} {'thresh log rho':>15} {'thresh score':>13}")
    for p in (2,3,4,7):
        for J in (0, 2, 4):
            r = arch_k(0.0, J, p)
            print(f"{p:>3} {J:>3} {r['m']:>4} {r['tau']:>8.5f} {r['threshold']:>15.5f} "
                  f"{r['threshold']-math.log(4)-2:>+13.5f}")
