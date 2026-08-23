from fractions import Fraction as Fr
from math import comb, log
from emn_core import series_H
N=1300
H=series_H(N)
bad=[]
for n in range(1,N+1):
    x=Fr(n*n)*comb(2*n,n)*H[n]
    if x.denominator!=1: bad.append((n,x.denominator))
print("n^2 C(2n,n) h_n integral? violations:", bad[:10], "count", len(bad))
# minimal: try n^e C(2n,n)
for e in (0,1,2):
    bad=[n for n in range(1,300) if (Fr(n)**e*comb(2*n,n)*H[n]).denominator!=1]
    print("e=",e,"violations up to 300:",len(bad), bad[:6])
print("c_n:", [ (Fr(n*n)*comb(2*n,n)*H[n]) for n in range(1,11)])
for n in (100,400,1300):
    print(n, "log(n^2 C)/n =", log(n*n*comb(2*n,n))/n, " log den(h_n)/n=", log(H[n].denominator)/n)
