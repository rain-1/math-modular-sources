import sys; sys.path.insert(0,'/home/ubuntu/code/math-modular-sources/lattice/emn')
from emnexact import linform
from fractions import Fraction as Fr
from math import comb
sgn = lambda t: 1 if t%2==0 else -1
print('=== (f) b_{m,t} = (-1)^t 4^{-m} C(m,m/2) C(m,t)  [note states it for even t] ===')
bad=0; tot=0
for m in range(0,19,2):
    for t in range(-4, m+1):
        a,b = linform(m,t); tot+=1
        cl = Fr(sgn(t)*comb(m,m//2)*(comb(m,t) if 0<=t<=m else 0), 4**m)
        if b != cl: bad+=1; print('   FAIL',m,t,b,cl)
print(f'   {tot} cases (m<=18, -4<=t<=m); mismatches: {bad}')
print()
print('=== pole-lowering recurrence (note sec.8) ===')
bad=0; tot=0
for m in range(2,19,2):
    for t in range(2,m+1,2):
        a,b=linform(m,t); a2,b2=linform(m-2,t-2); tot+=1
        f = Fr((m-1)**2,4*t*(t-1))
        pa = Fr(-1, 2**(t+1)*t*(t-1)*comb(2*(m-t),m-t)) + f*a2
        pb = f*b2
        if (a,b)!=(pa,pb): bad+=1; print('   FAIL (m,t)=',m,t,(a,b),(pa,pb))
print(f'   {tot} cases; failures: {bad}')
print()
print('=== exact values (cross-check vs the 140-digit PSLQ table) ===')
for (m,t) in [(6,0),(6,2),(8,4),(10,6),(10,10)]:
    print('  I_%d,%d ='%(m,t), linform(m,t))
