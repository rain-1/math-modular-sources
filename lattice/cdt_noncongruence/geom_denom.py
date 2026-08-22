"""How CDT's framework handles a GEOMETRIC denominator lam^n.

Setting.  A modular row has a_n = lam^n [t^n] g with lam in {1,2,4}: in the
modular coordinate t the coefficients of g carry an unbounded power of 2,
c_n := [t^n] g in lam^{-n} Z.  A geometric denominator is NOT of CDT's LCM type
(6.0.9), so g(t) is not admissible as an f_i in the coordinate t.

Two ways to make it admissible; they agree exactly.

(A) RESCALE.  Put x := t/lam.  Then A(x) := g(lam x) = sum a_n x^n has a_n in Z,
    and every singular point moves x_i = t_i/lam.  A map phi_t : D -> Omega_t is
    the same thing as phi_x = phi_t/lam : D -> Omega_x, so
        log|phi_x'(0)| = log|phi_t'(0)| - log lam ,
        BC(phi_x)      = BC(phi_t)      - log lam
    (the Bost-Charles integral is log-homogeneous of degree 1 in phi).  tau is
    unchanged.  Hence
        entry_x  = entry_t - log lam                      <-- the -log lam of N1
        margin_x = margin_t - (m-1) log lam .

(B) KEEP t, DECLARE A NEGATIVE p-ADIC SLOPE.  In the coordinate t every function
    has v_p(c_n) >= -v_p(lam) n, i.e. the UNIFORM slope profile
    varsigma_i = -v_p(lam) for all i, at every p | lam.  The adelic bound
    (ADELIC_HOLONOMY.md Thm 2.2) then gives
        gamma_p = -v_p(lam) log p (1 - 1/m),   sum_p gamma_p = -(1-1/m) log lam ,
        tau_ad = tau + (1-1/m) log lam .
    Hence margin_ad = margin_t - (m-1) log lam.  IDENTICAL to (A).

So: a geometric denominator costs exactly  log lam  of ENTRY budget and
(m-1) log lam of MARGIN, and per function (m -> infinity) exactly log lam --
which is precisely the -log lam in NONCONGRUENCE_SCAN's Theorem N1,
  score = log|t_2| - k - log lam = log(1/|lam_2|) - k,   lam_2 = lam/|t_2| .
This script verifies (A) == (B) numerically.
"""
import math, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'cdt_finder'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'adelic_holonomy'))
from cdt_bound import tau_flat, tau_sharp
from adelic_bound import adelic

def check(m, cols, e, log_phi_t, BC_t, lam, p, v):
    """lam = p^v.  (A) rescale to x = t/lam ; (B) uniform slope -v at p."""
    L = math.log(lam)
    A = adelic(m, cols, e, {}, log_phi_t - L, BC_t - L)
    B = adelic(m, cols, e, {p: [-v]*m}, log_phi_t, BC_t)
    return A, B

if __name__ == '__main__':
    print(__doc__)
    print("="*94)
    print(f"{'m':>4} {'lam':>4} {'entry(A)':>10} {'entry(B)':>10} {'margin(A)':>11} {'margin(B)':>11} {'diff':>10}")
    for m in (2, 3, 5, 14, 30):
        cols = [(1,1),(1,1)]
        e = [0]*m
        for lam, p, v in ((2,2,1), (4,2,2), (8,2,3)):
            A, B = check(m, cols, e, math.log(33.9706), 11.845, lam, p, v)
            # (B)'s entry is tau_ad-based and differs by log lam / m; margins must agree
            print(f"{m:>4} {lam:>4} {A['entry']:>+10.6f} {B['entry']:>+10.6f} "
                  f"{A['margin']:>+11.6f} {B['margin']:>+11.6f} {A['margin']-B['margin']:>10.2e}")
    print("\nentry(A) - entry(B) = -log(lam)/m  (the 1/m in gamma_p's (1-1/m)); margins agree exactly.")
    print("\n--- N1 consistency: Beukers' row, modular coordinate t vs integral x ---")
    for nm, t2, lam in (("Beukers Gamma_0(6)+6", (1+math.sqrt(2))**4, 4),
                        ("level-5 Gamma_0(5)+5", 2/0.7213595499958, 2)):
        lam2 = lam/t2
        print(f"  {nm}: |t_2|={t2:.6f}, lam={lam}, |lam_2|={lam2:.8f}, "
              f"|x_2|=|t_2|/lam={t2/lam:.6f}=1/|lam_2|={1/lam2:.6f}  "
              f"score=log|t_2|-2-log(lam)={math.log(t2)-2-math.log(lam):+.4f}")
