import sys
from fractions import Fraction as F
sys.path.insert(0,'/home/ubuntu/code/math-modular-sources/lattice/two_worlds')
from chi3_forms import *
import mpmath as mp
mp.mp.dps = 80
L3 = mp.nsum(lambda n: (mp.mpf(1)/(3*n+1)**2 - mp.mpf(1)/(3*n+2)**2), [0, mp.inf])

def numerator(m, a):
    """N = (2t+m+1) * prod_{j=1..a}(t-j+1) * prod_{j=1..a}(t+m+j) ; sigma-invariant, deg 2a+1"""
    N = linear(2, m+1)
    for j in range(1, a+1):
        N = pmul(N, linear(1, -j+1))
    for j in range(1, a+1):
        N = pmul(N, linear(1, m+j))
    return N

def row(m, a):
    P = poles_thirds(m)
    d = linform(P, numerator(m, a))
    assert d['zeta2'] == 0 and d['pi_over_sqrt3'] == 0, (m, a, d['zeta2'], d['pi_over_sqrt3'])
    return d['L'], -d['Crat']          # Q_m , P_m  with  Q_m*L - P_m = S_m

if __name__ == "__main__":
    import sys
    amode = sys.argv[1] if len(sys.argv) > 1 else 'lin'
    for spec in [('a=0',lambda m:0), ('a=m',lambda m:m), ('a=2m',lambda m:2*m),
                 ('a=m//2',lambda m:m//2), ('a=3m//2',lambda m:(3*m)//2)]:
        name, af = spec
        print("=== %s ===" % name)
        prevQ = None
        for m in range(1, 15):
            a = af(m)
            Q, Pp = row(m, a)
            form = mp.mpf(Q.numerator)/Q.denominator*L3 - mp.mpf(Pp.numerator)/Pp.denominator
            print("  m=%2d a=%2d  v3(Q)=%3d den(Q)=%s  |Q|^(1/m)=%s  |form|^(1/m)=%s" % (
                m, a, __import__('sympy').multiplicity(3, Q.numerator) - __import__('sympy').multiplicity(3, Q.denominator) if Q else 0,
                Q.denominator, mp.nstr(mp.mpf(abs(Q.numerator))/Q.denominator, 6)+"",
                mp.nstr(abs(form)**(mp.mpf(1)/m), 8)))
