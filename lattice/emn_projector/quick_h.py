from fractions import Fraction as Fr
import math
from sympy import primerange, factorint

N=400
h=[Fr(0)]*(N+1)
prev=Fr(0)
for n in range(1,N+1):
    h[n]=(Fr(1,2)+Fr((n-1)**2)*prev)/Fr(n*(2*n-1))
    prev=h[n]
print("first h:", [str(h[n]) for n in range(1,7)])
# denominators
for n in (10,20,50,100,200,400):
    d=h[n].denominator
    print(n, "log den =", math.log(d), " /n =", math.log(d)/n, "factor:", dict(list(factorint(d).items())[:8]) if d< 10**60 else "big")
