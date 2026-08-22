"""Final verification of the alpha=2 row (the chi_{-3} analogue of Zudilin's Catalan row)."""
import sys, math
sys.path.insert(0,'/home/ubuntu/code/math-modular-sources/lattice/two_worlds')
from fractions import Fraction as F
from chi3_row import row
from chi3_padic import rowC, v3
import mpmath as mp
mp.mp.dps = 300
L = mp.nsum(lambda n: (mp.mpf(1)/(3*n+1)**2 - mp.mpf(1)/(3*n+2)**2), [0, mp.inf])
print("L(2,chi_-3) =", mp.nstr(L, 30))
print()
print("first terms of the row (a=2m, b=m):")
for m in range(6):
    Q,P = row(m,2*m)
    print("  m=%d  Q=%s   P=%s"%(m,Q,P))
print()
print(" m   |Q_m L - P_m|            m^3*|form|     P_m/Q_m - L      v3(Q_m)")
for m in [4,8,16,32,64]:
    Q,P = row(m,2*m)
    f = mp.mpf(Q.numerator)/Q.denominator*L - mp.mpf(P.numerator)/P.denominator
    print(" %2d   %-22s   %-12s  %-16s  %d"%(m, mp.nstr(f,8), mp.nstr(f*m**3,6),
          mp.nstr(mp.mpf(P.numerator)/P.denominator/(mp.mpf(Q.numerator)/Q.denominator)-L,6), v3(Q)))
print()
# 3-adic: xi_3 = zeta_3(2) = 2*xi_3^C
a,b = rowC(200)
xiC = b[199]/a[199]
print("3-adic limit test  v_3(P_m/Q_m - 2*xi_3^C):")
print("   m:", list(range(10,71,10)))
print("  v3:", [v3(row(m,2*m)[1]/row(m,2*m)[0] - 2*xiC) for m in range(10,71,10)])
print("  predicted 3m - 2 s_3(m) + s_3(2m)/2 - 1 =",
      [3*m - 2*sum(int(c) for c in mp.nstr(0,1)) for m in [0]][:0] or
      [3*m - 2*(lambda n: sum(int(d) for d in __import__('numpy').base_repr(n,3)))(m)
         + (lambda n: sum(int(d) for d in __import__('numpy').base_repr(n,3)))(2*m)//2 - 1
       for m in range(10,71,10)])
