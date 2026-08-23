from fractions import Fraction as Fr
import math
from mpmath import mp, mpf, catalan, cos, sin, sec, quad, pi, log, tan, mpmathify
mp.dps=40

N=3000
h=[Fr(0)]*(N+1); prev=Fr(0)
for n in range(1,N+1):
    h[n]=(Fr(1,2)+Fr((n-1)**2)*prev)/Fr(n*(2*n-1)); prev=h[n]

# check H(1/2) = G/3 and H(3/2)=5G/3 (analytic continuation needed for 3/2) and H(1)=G
def Hval(z, M=2500):
    s=mpf(0); zz=mpf(z)
    for n in range(1,M+1):
        s+= mpf(h[n].numerator)/mpf(h[n].denominator)*zz**n
    return s
G=catalan
print("H(1/2)      =", Hval(0.5))
print("G/3         =", G/3)
print("H(1)  (slow)=", Hval(1.0,20000) if False else "skip")
# radius of convergence via |h_n|^{-1/n}
import statistics
for n in (500,1000,2000,3000):
    hn=abs(mpf(h[n].numerator)/mpf(h[n].denominator))
    print("n=",n," |h_n|^{-1/n} =", float(hn**(mpf(-1)/n)), " n^2 h_n=", float(hn*n*n))
