"""Task 2: tau and the entry test for the EMN moving-period host and for the
radial-Pade family, in the units of CDT_UNPACKED / CATALAN_MU4."""
import sys, os, math
sys.path.insert(0, '/home/ubuntu/code/math-modular-sources/lattice/cdt_finder')
from fractions import Fraction
from cdt_bound import tau_flat, tau_sharp, tau

CEIL = {
 "CDT own host P1-{0,1,oo}"          : math.log(16),
 "Catalan level 8 symmetrised"       : math.log(64),
 "mu_4 host A"                       : math.log(4),
 "EMN z-line, omit {1}, Koebe(univ)" : math.log(4),
 "EMN z-line, omit {1,2}, univ.cover": math.log(8.753758460905907),
 "EMN S/zeta-plane, C\\{+-1} at 0"    : math.log(4.376879230452953),
 "EMN zeta-line, omit {+-1/sqrt2}"   : math.log(3.094920984287841),
}

print("=== calibration: reproduce CDT's own numbers ===")
t = tau(14, [(1,2),(3,2)], [0,0,1,0,0,0,0,0,0,1,1,1,1,1])
print(f"  CDT Thm A: tau_flat={t['tau_flat']}  tau_sharp={t['tau_sharp']:.6f}  tau={t['tau']:.6f}  (expect 4.235459)")
t = tau(17, [(2,1),(4,1)], [0]*3+[1]*14)
print(f"  CDT sec.7 x-line: tau={t['tau']:.6f}  (expect 2.275626)")

print()
print("=== EMN host: the module and its denominator types ===")
print("  measured (n<=400, zero exceptions):  den(a_n) | 2*(n+1)*[1..2n+1]")
print("     => H(z) = sum a_n z^{n+1} has CDT type  n * [1..2n]  : e=1, one column with b=2")
print("     true rate log den(a_n)/n = 1.395 at n=400 (still falling; d_{2n+1} rate is 2.0)")
print("  the polylogarithm module at the single puncture z=1 is lcm-free:")
print("     Li_k(z)=sum z^n/n^k has type n^k (e=k, no lcm layer); log(1-z) type n;")
print("     products log^j(1-z) Li_k(z) carry [1..n]^j layers.")
print()
print(f"{'inventory':<52} {'m':>3} {'cols':>14} {'tau':>8} " + "  ".join(f"{k[:26]:>27}" for k in ["entry @ Koebe log4","entry @ C\\{1,2} 2.169","entry @ S-plane 1.476"]))

def row(label, m, cols, e):
    T = tau(m, cols, e)
    e1 = math.log(4) - T['tau']
    e2 = math.log(8.753758460905907) - T['tau']
    e3 = math.log(4.376879230452953) - T['tau']
    print(f"{label:<52} {m:>3} {str(cols):>14} {T['tau']:>8.4f} {e1:>27.4f} {e2:>27.4f} {e3:>27.4f}")

# inventory 1: {1, H} only
row("{1, H}", 2, [(1,2)], [0,1])
# inventory 2: {1, H, int H, int zH}  (all carry H's layer)
row("{1, H, IntH, Int zH}", 4, [(1,2)], [0,1,2,2])
# inventory 3: {1, log(1-z), Li_2, Li_3, H} - 4 lcm-free
row("{1, A, Li_2, Li_3, H}", 5, [(4,2)], [0,1,2,3,1])
row("{1,A,Li_2,Li_3,Li_4,Li_5, H}", 7, [(6,2)], [0,1,2,3,4,5,1])
row("{...10 lcm-free..., H}", 11, [(10,2)], [0,1,2,3,4,5,6,7,8,9,1])
row("{...10 lcm-free..., H, IntH, Int zH}", 13, [(10,2)], [0,1,2,3,4,5,6,7,8,9,1,2,2])
# zeta-line: index doubles so b halves
print()
print("  zeta-line (zeta^2 = z/2): the index of H doubles, so its layer becomes [1..N]: b=1")
def rowz(label, m, cols, e):
    T = tau(m, cols, e)
    print(f"{label:<52} {m:>3} {str(cols):>14} {T['tau']:>8.4f} "
          f"{math.log(3.094920984287841)-T['tau']:>27.4f} {math.log(4.376879230452953)-T['tau']:>27.4f}")
print(f"{'inventory (zeta-line)':<52} {'m':>3} {'cols':>14} {'tau':>8} {'entry @ 1.1298':>27} {'entry @ 1.4763':>27}")
rowz("{1, H}", 2, [(1,1)], [0,1])
rowz("{1, A, Li_2, Li_3, H}", 5, [(4,1)], [0,1,2,3,1])
rowz("{...10 lcm-free..., H}", 11, [(10,1)], [0,1,2,3,4,5,6,7,8,9,1])
