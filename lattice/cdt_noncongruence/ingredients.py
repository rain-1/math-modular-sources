"""What each additional ingredient is worth, in score-units of threshold.

Baseline: one second-order row, one 2-term relation a + b xi = 0.  The whole
admissible inventory is {1, H, theta H} (the inhomogeneous 2nd-order ODE gives a
Q(x)-relation on {1,H,thH,th^2 H}), i.e. m = 3, and the reach is score > -0.8385.

(i)  c independent FOLD-REGULAR conditional classes (each contributing H_i and
     theta H_i, all e_i = 0):  m = 2c+1, u_j = 1, b_j = 1.
(ii) a pure polylogarithm module (p functions, sigma = 0): needs lam_2 in Z.
(iii) the normaliser descent (CDT's symmetrisation): needs s = 1/lam_2 in Q;
     ceiling 16 -> 256 but sigma_m 2 -> 4.
(iv) a p-adic slope: shifts the threshold by -m gamma/(m-1) ~ -gamma.
(v)  a Galois trace to a number field K: replaces log|x_2| by
     -log|N(lam_2)|^{1/[K:Q]} >= -log|x_2|... a strict LOSS for any decaying row.
"""
import math, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'cdt_finder'))
from cdt_bound import tau_flat, tau_sharp
from optimum import G
L4 = math.log(4)

def thresh_classes(c, sym=False):
    """c fold-regular conditional classes: m = 2c+1 functions {1} u {H_i, th H_i}."""
    m = 2*c + 1
    b = 2 if sym else 1
    sm, tf = tau_flat(m, [(1, b), (1, b)] if not sym else [(1,2),(1,2)])
    T = float(tf)
    g, r = G(m)
    if sym:
        # symmetrised: ceiling 256|s| instead of 16|s| ; no numerically computed
        # Delta available for the descended map, so quote the CEILING only.
        return None
    return (m, T, g, r, (m*T - g)/(m-1) - L4 - 2)

if __name__ == '__main__':
    print(__doc__)
    print("(i) fold-regular conditional classes, unsymmetrised Kodaira map:")
    print(f"{'classes c':>10}{'m':>5}{'tau':>9}{'G(m)':>9}{'r*':>6}{'threshold score':>17}{'gain':>8}")
    base = None
    for c in (1,2,3,4,5,6,8,10):
        m, T, g, r, th = thresh_classes(c)
        if base is None: base = th
        print(f"{c:>10}{m:>5}{T:>9.5f}{g:>9.5f}{str(r):>6}{th:>17.5f}{th-base:>+8.4f}")
    print("  asymptote (c -> oo, r -> 1): threshold -> 2 - log 16 - 2 = %.5f" % (-math.log(16)))
    print("\n(iii) the symmetrised (descent) architecture, hard ceiling, u_j = 1:")
    print("      entry log(256|s|) > tau -> sigma_m = 4, so asymptotically")
    print("      log(256|s|) > 4  <=>  normalised score > 4 - log 256 - 2 = %.5f" % (4-math.log(256)-2))
    print("      CDT's realised m = 14, u = (1,3): threshold -2.00041 (contour) / -2.49931 (ceiling)")
    print("\n(v) Galois trace: for a real-quadratic fold, log|x_2| - (-log lam2_norm)")
    for nm, l2, de in (("Beukers", 4*(17-12*2**0.5), 16), ("level-5", 44-20*5**0.5, -64),
                       ("sqrt-T", 12-8*2**0.5, 16)):
        print(f"      {nm:<10s}: log|x_2| = {-math.log(abs(l2)):+.4f}, "
              f"-log lam2_norm = {-0.5*math.log(abs(de)):+.4f}  -> loss {-math.log(abs(l2))+0.5*math.log(abs(de)):+.4f}")
