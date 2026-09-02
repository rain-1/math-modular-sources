"""12: the conditional function H = q*s - (p-q)*c  (i.e. s - (zeta(3)-1) c) of the
zeta(3) D1 two-variable system: its decay rate, hence its domain of convergence,
hence whether it is 'fold-regular' in the two-variable sense.

Also: the two-variable conformal ceiling via Landau's theorem.
"""
from lib2v import *
from fractions import Fraction as F
from math import comb, log
import mpmath as mp
mp.mp.dps = 200

NA = NB = 40


def build(t00, t10, t01, t11, NA=NA, NB=NB):
    t = [[F(0)]*(NB+1) for _ in range(NA+1)]
    t[0][0], t[1][0], t[0][1], t[1][1] = map(F, (t00, t10, t01, t11))
    for b in (0, 1):
        B = b*(b+1)
        for n in range(1, NA):
            mid = 2*n**3 + 3*n**2 + (3+4*B)*n + 1 + 2*B
            t[n+1][b] = (F(mid)*t[n][b] - F(n**3)*t[n-1][b])/F((n+1)**3)
    for a in range(NA+1):
        for b in range(NB-1):
            Q1 = -(2*b**3 + 9*b**2 + (15+4*a+4*a*a)*b + 9 + 6*a + 6*a*a)
            t[a][b+2] = (-F(Q1)*t[a][b+1] - F((b+1)**3)*t[a][b])/F((b+2)**3)
    return t


c = build(1, 1, 1, 5)
s = build(-1, 0, 0, 1)
xi = mp.zeta(3) - 1


def tomp(q):
    return mp.mpf(q.numerator)/mp.mpf(q.denominator)


print("=" * 78)
print("A. decay of H_{a,b} = s_{a,b} - (zeta(3)-1) c_{a,b}")
print("=" * 78)
print("   log|H|/(a+b) should tend to -log(radius) of the domain in which the")
print("   conditional function converges.  Reference values:")
print("     -log(3+2sqrt2) = %.5f  (the OUTER branch of 16xy=(1-x)^2(1-y)^2)"
      % -float(mp.log(3+2*mp.sqrt(2))))
print("     -log(3-2sqrt2) = %.5f  (the FOLD branch = radius of the row C)"
      % -float(mp.log(3-2*mp.sqrt(2))))
print()
print("   (a,b)     log|c|/(a+b)   log|H|/(a+b)")
for (a, b) in [(10, 10), (20, 20), (30, 30), (40, 40), (40, 20), (20, 40),
               (40, 13), (13, 40), (40, 8), (8, 40)]:
    H = tomp(s[a][b]) - xi*tomp(c[a][b])
    lc = log(c[a][b])/(a+b)
    lh = float(mp.log(abs(H)))/(a+b)
    print("   (%2d,%2d)     %8.5f      %8.5f" % (a, b, lc, lh))

print()
print("=" * 78)
print("B. one-variable comparison: Apery's own H_n = b_n - zeta(3) a_n")
print("=" * 78)
a3, b3 = apery3(40)
for n in [10, 20, 30, 40]:
    Hn = tomp(b3[n]) - mp.zeta(3)*mp.mpf(a3[n])
    print("   n=%2d  log|a_n|/n = %8.5f   log|H_n|/n = %8.5f"
          % (n, log(a3[n])/n, float(mp.log(abs(Hn)))/n))
print("   references: log(17+12sqrt2) = %.5f ; -log(17+12sqrt2) = %.5f"
      % (float(mp.log(17+12*mp.sqrt(2))), -float(mp.log(17+12*mp.sqrt(2)))))

print()
print("=" * 78)
print("C. the two-variable conformal ceiling (Landau)")
print("=" * 78)
print("""   The host is the pullback of the u-line along u = alpha(x) beta(y) with
   alpha(x) = x/(1-x)^2, beta(y) = y/(1-y)^2 (Koebe functions, alpha'(0)=1).
   An admissible pair (phi_1,phi_2) gives A = (alpha o phi_1)(D), B = (beta o phi_2)(D),
   both open, connected, containing 0, with
        A . B  avoiding  u = 1/16   (the branch point of h = 2F1(1/2,1/2;1;16u)).
   Since 0 is an interior point of A and of B, A contains a disc D(0,rho_A) and
   B a disc D(0,rho_B), and D(0,rho_A) . D(0,rho_B) = D(0,rho_A rho_B) must miss
   1/16, so rho_A rho_B <= 1/16.
   LANDAU's theorem (sharp constant 1/16, extremal = the modular function):
   if f is holomorphic on D, f(0)=0 and f omits w_0, then |w_0| >= |f'(0)|/16;
   equivalently the image contains D(0,|f'(0)|/16).  Hence
        rho_A >= e^{l_1}/16,   rho_B >= e^{l_2}/16,
   and therefore   e^{l_1+l_2}/256 <= 1/16, i.e.""")
print("        l_1 + l_2 <= log 16 = %.5f      (SHARP: modular maps)" % log(16))
print()
print("   Entry conditions, for the three denominator geometries:")
print("     product/split (r1,r2), r1+r2=3 : need l_1>r1 and l_2>r2, so l_1+l_2>3")
print("     min-type [1..min(a,b)]^3       : need l_1,l_2>0 and l_1+l_2>3")
print("     max-type [1..max(a,b)]^3       : need l_1>3 and l_2>3, so l_1+l_2>6")
print("   ceiling log16 = %.4f, so:" % log(16))
print("     split / min-type : deficit %.4f nats  (= the one-variable Apery deficit"
      " 3 - log 16)" % (3 - log(16)))
print("     max-type         : deficit %.4f nats" % (6 - log(16)))
print()
print("   For a rate-r row the same computation gives deficit r - log16 (min type)")
print("   and 2r - log16 (max type):")
for r in (2, 3):
    print("     r=%d :  min-type deficit %+.4f ,  max-type deficit %+.4f"
          % (r, r - log(16), 2*r - log(16)))
