from fractions import Fraction as F
from math import gcd
Nn=40
def lcmr(n):
    l=1
    for k in range(1,n+1): l=l*k//gcd(l,k)
    return l
# H_A = sum binom(2n,n) x^n  (m=1)
def binom(n,k):
    r=1
    for i in range(k): r=r*(n-i)//(i+1)
    return r
A=[F(binom(2*n,n)) for n in range(Nn+2)]
# log(1-t) = -sum_{n>=1} t^n/n  ; log(1-t)/t = -sum_{n>=0} t^n/(n+1)
LG=[F(0)]+[F(-1,n) for n in range(1,Nn+2)]
LGt=[F(-1,n+1) for n in range(Nn+2)]        # log(1-t)/t
INV=[F(1)]*(Nn+2)                            # 1/(1-t)
def mul(u,v,N=Nn+1):
    return [sum(u[i]*v[n-i] for i in range(n+1)) for n in range(N+1)]
def integ(u):   # int_0^x, coefficient shift
    r=[F(0)]*(Nn+2)
    for n in range(Nn+1): r[n+1]=u[n]/(n+1)
    return r
LGo1m=mul(LG,INV)                            # log(1-t)/(1-t)
HC=mul(A,integ(mul(A,LGt)))
HD=mul(A,integ(mul(A,LGo1m)))
print("== Task 6: m=1, independent python/Fraction re-derivation ==")
print("H_C = (1-4x)^{-1/2} int_0^x log(1-t)/(t sqrt(1-4t)) dt")
print("  first 8 coefficients a_0..a_7:")
for n in range(8): print("    a_%d = %s"%(n,HC[n]))
print("H_D = (1-4x)^{-1/2} int_0^x log(1-t)/((1-t) sqrt(1-4t)) dt")
for n in range(8): print("    a_%d = %s"%(n,HD[n]))
print()
print("  n | den(a_n) H_C | L(n)L(n/2) | L(n)^2 | den | L(n)L(n/2)?  | den | L(n)^2 ?")
for n in range(1,25):
    d=HC[n].denominator; T12=lcmr(n)*lcmr(n//2); T2=lcmr(n)**2
    print("  %2d | %-22d | %-22d | %-24d | %-5s | %-5s"%(n,d,T12,T2,T12%d==0,T2%d==0))
print()
print("  shifted conventions for H_C (is den(a_n) | L(n+s)L((n+s)/2) ?), s=0,1,2:")
for s in range(3):
    ok=[n for n in range(1,Nn+1) if (lcmr(n+s)*lcmr((n+s)//2))%HC[n].denominator!=0]
    print("    s=%d: first failures %s"%(s,ok[:8]))
print("  is den(a_n) | L(n)^2 for all n<=%d ? %s"%(Nn, all((lcmr(n)**2)%HC[n].denominator==0 for n in range(1,Nn+1))))
print("  is den(a_n) | L(n)L(n/2) for H_D for all n<=%d ? %s"%(Nn, all((lcmr(n)*lcmr(n//2))%HD[n].denominator==0 for n in range(1,Nn+1))))
