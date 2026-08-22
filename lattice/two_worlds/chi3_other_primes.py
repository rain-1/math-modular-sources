"""Do the chi_{-3} hypergeometric rows have p-adic Apery limits at p != 3?"""
import sys
from fractions import Fraction as F
sys.path.insert(0,'/home/ubuntu/code/math-modular-sources/lattice/two_worlds')
from chi3_row import row

def vp(x,p):
    if x==0: return None
    n,d=x.numerator,x.denominator; v=0
    while n%p==0: n//=p; v+=1
    while d%p==0: d//=p; v-=1
    return v

for (name,af) in [("a=2m",lambda m:2*m),("a=3m/2",lambda m:(3*m)//2)]:
    print("=== %s ===" % name)
    R=[row(m,af(m)) for m in range(0,26)]
    for p in (2,3,5,7,11,13):
        inc=[]
        for m in range(2,26):
            Q,P=R[m]; Q0,P0=R[m-1]
            if Q==0 or Q0==0: inc.append(None); continue
            inc.append(vp(P/Q-P0/Q0,p))
        print("  p=%-3d increments v_p(P_m/Q_m - P_{m-1}/Q_{m-1}) m=2..25: %s" % (p, inc))
