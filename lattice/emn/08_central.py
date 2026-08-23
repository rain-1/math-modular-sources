"""(g) the central higher-pole family  I_n = I_{4n,2n} = Acen_n + Bcen_n * G.
Uses the pole-lowering chain (verified in 07) to reduce to I_{2n,0}, which is
evaluated exactly by emnexact.  Measures denominators, fits the recurrence."""
import sys, time, json
sys.path.insert(0,'/home/ubuntu/code/math-modular-sources/lattice/emn')
from emnexact import linform
from fractions import Fraction as Fr
from math import comb, log, gcd

NMAX = int(sys.argv[1]) if len(sys.argv)>1 else 22

def central(n):
    """I_{4n,2n} = (a,b) exactly, via n pole-lowering steps down to I_{2n,0}."""
    if n == 0: return linform(0,0)
    a0, b0 = linform(2*n, 0)
    a, b = a0, b0
    # walk up: (2n,0) -> (2n+2,2) -> ... -> (4n,2n)
    for k in range(1, n+1):
        m, t = 2*n+2*k, 2*k
        f = Fr((m-1)**2, 4*t*(t-1)) if t>=2 else None
        add = Fr(-1, 2**(t+1)*t*(t-1)*comb(2*(m-t), m-t))
        a = add + f*a
        b = f*b
    return a, b

t0=time.time()
Acen=[]; Bcen=[]
for n in range(NMAX+1):
    a,b = central(n)
    Acen.append(a); Bcen.append(b)
    print(f"n={n:>3}  B_n=C(4n,2n)^2/256^n? {'OK' if b==Fr(comb(4*n,2*n)**2, 256**n) else 'NO '}"
          f"  den(A_n) bits={Fr(a).denominator.bit_length():>6}  num(A_n) bits={abs(Fr(a).numerator).bit_length():>6}"
          f"  t={time.time()-t0:.1f}s")
json.dump([[str(x) for x in Acen],[str(x) for x in Bcen]], open('/home/ubuntu/code/math-modular-sources/lattice/emn/central.json','w'))
