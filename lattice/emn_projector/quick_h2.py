from fractions import Fraction as Fr
import math
from sympy import primerange, factorint

N=600
h=[Fr(0)]*(N+1)
prev=Fr(0)
for n in range(1,N+1):
    h[n]=(Fr(1,2)+Fr((n-1)**2)*prev)/Fr(n*(2*n-1))
    prev=h[n]

n=300
d=h[n].denominator
f=factorint(d)
ps=sorted(f)
print("n=",n,"log den/n=",math.log(d)/n)
print("primes dividing den, with exponent (p<=40):", {p:f[p] for p in ps if p<=40})
print("max prime:",ps[-1], " 2n=",2*n)
# which primes in [1,2n] are ABSENT
absent=[p for p in primerange(2,2*n+1) if p not in f]
print("num primes<=2n:",len(list(primerange(2,2*n+1))),"absent:",len(absent))
print("absent sample:",absent[:20], "...", absent[-10:])
present=[p for p in ps]
print("present in (n,2n]:",len([p for p in present if n<p<=2*n]), "of", len([p for p in primerange(n+1,2*n+1)]))
print("present in (2n/3,n]:",len([p for p in present if 2*n//3<p<=n]),"of",len([p for p in primerange(2*n//3+1,n+1)]))
print("present in (n/2,2n/3]:",len([p for p in present if n//2<p<=2*n//3]),"of",len([p for p in primerange(n//2+1,2*n//3+1)]))
