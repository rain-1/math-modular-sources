"""Cost accounting for every route to a SECOND fold-regular class on the
(Gamma_0(5)+5, t) host, in score-units of threshold.

Baseline (D): inventory {1, H, theta H}, m = 3, tau = 16/9, threshold -0.83852.
c fold-regular classes: m = 2c+1, u_j = 1, b_j = 1, all e_i = 0.
"""
import math, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'cdt_finder'))
from cdt_bound import tau_flat, tau_sharp
from optimum import G
L4 = math.log(4)

def thresh(m, cols, e=None):
    e = e or [0]*m
    T = float(tau_flat(m, cols)[1]) + tau_sharp(m, e)[0]
    g, r = G(m)
    return (m*T - g)/(m-1) - L4 - 2, T, g, r

def show(lab, m, cols, e=None, extra_cost=0.0):
    th, T, g, r = thresh(m, cols, e)
    th += extra_cost
    print(f"  {lab:<52s} m={m:2d} tau={T:7.4f} G={g:6.3f} r*={str(r):>5s} "
          f"threshold={th:+8.4f}  deficit(level-5)={-1.673382-th:+8.4f}")
    return th

if __name__ == '__main__':
    print(__doc__)
    SC = -1.673382                       # the level-5 row's elementary score
    print(f"level-5 score = {SC:+.6f}\n")
    print("BASELINE and the ideal c-class ladder (no extra cost):")
    for c in (1,2,3):
        show(f"c = {c} fold-regular classes (ideal)", 2*c+1, [(1,1),(1,1)])
    print()
    print("ROUTE (a) nebentypus M_2(Gamma_0(5),chi_5): E^{1,chi_5} has zeros in H.")
    t2, textra = 2.7725424859373686, 1.0204164998665332
    lam, lamp = 2, 4
    cost_lam = math.log(lamp/lam)
    cost_sing = math.log(t2/textra)
    print(f"    extra singular t = (-6 +- sqrt(39))/12 = {(-6+39**0.5)/12:.10f}, {(-6-39**0.5)/12:.10f}")
    print(f"    lambda' = 4 vs lambda = 2  -> common coordinate x = t/4, cost log2 = {cost_lam:.5f}")
    print(f"    |t_extra| = {textra:.7f} < |t_2| = {t2:.7f} -> contour cost log = {cost_sing:.5f}")
    print(f"    TOTAL archimedean cost = {cost_lam+cost_sing:.5f}   (gain from c=2 is +0.5852)")
    show("c = 2 via the nebentypus class", 5, [(1,1),(1,1)], extra_cost=cost_lam+cost_sing)
    print()
    print("ROUTE (b) levels 10,15,20,25: the ONLY order-2 same-singular-set element")
    print("    of M_2(Gamma_0(N)) is E_{2,5} itself (22_gauge_family.log) -- no new class.")
    print()
    print("ROUTE (c) gauge family F*P(t)^b, P = 1 - 44t - 16t^2:")
    print("    b even: g' = g*P^{b/2} in Q(t)*g  ->  A' = P(4x)^{b/2} A  ->  Q(x)-DEPENDENT.")
    print("    b = 1 : order 2, same roots, but k' = 10 (measured) -> sigma = 10 for that function:")
    show("c = 2 via F*P (one function of sigma = 10)", 4,
         [(1,1),(1,1)] + [(3,1)]*8)
    print("    b = 3 : order 4 with DOUBLED characteristic roots (ROOT_ROWS (H4)):")
    print("            the linear form converges polynomially, no fold-regular companion.")
    print()
    print("ROUTE (d) the Sym^w tower of the SAME operator L_1 (w = 2: k = 3, sigma = 3):")
    show("c = 1 + Sym^2 pair (sigma = 0,2,2,3,3)", 5, [(1,1),(1,1),(3,1)])
    print("    (and the Sym^2 conditional function needs its OWN rationality hypothesis)")
    print()
    print("ROUTE (e) the algebraic denominator-free function sqrt(P(t)) (integer")
    print("    coefficients, sigma = 0!): branches at BOTH t_1 and t_2, so it is NOT")
    print("    single-valued on the domain -- the Galois obstruction of section 10.2.")
