from fractions import Fraction as Fr
from math import comb, log
from emn_core import series_H
from functools import reduce
N=600
H=series_H(N)
Ht=[H[n]*Fr(2)**n for n in range(N+1)]
part=[Fr(0)]*(N+1); run=Fr(0)
for n in range(N+1):
    part[n]=run; run+=Ht[n]
B=[Fr(0)]+[part[n]/Fr(n) for n in range(1,N+1)]
def lcm(a,b):
    from math import gcd
    return a*b//gcd(a,b)
L=[1]*(N+1)
cur=1
for n in range(1,N+1):
    cur=lcm(cur,Ht[n].denominator); L[n]=cur
for n in (100,200,400,600):
    print(n,"log den(B_n)/n=",round(log(B[n].denominator)/n,4),
          " log lcm_{k<=n}den(Htil_k)/n=",round(log(L[n])/n,4))
# shape test: n * lcm(1..n)^2 * lcm(1..2n) ?
def lcmrange(m):
    r=1
    for i in range(1,m+1): r=lcm(r,i)
    return r
for n in (50,100,200):
    cand = n*lcmrange(n)**2*lcmrange(2*n)
    print(n, "B_n*cand integral?", (Fr(cand)*B[n]).denominator==1,
          " n*lcm(1..n)*lcm(1..2n):", (Fr(n*lcmrange(n)*lcmrange(2*n))*B[n]).denominator==1,
          " n*lcm(1..2n)^2:", (Fr(n*lcmrange(2*n)**2)*B[n]).denominator==1,
          " n^2*C(2n,n)*lcm(1..n)^2:", (Fr(n*n*comb(2*n,n)*lcmrange(n)**2)*B[n]).denominator==1)
