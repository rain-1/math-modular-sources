from fractions import Fraction as Fr
from math import comb, gcd, log
from emn_core import series_H
N=220
H=series_H(N)
Ht=[H[n]*Fr(2)**n for n in range(N+1)]
conv=[Fr(0)]*(N+1); acc=Fr(0)
for n in range(N+1):
    acc = acc + Ht[n]
    conv[n]=acc
B=[Fr(0)]+[conv[n-1]/Fr(n) for n in range(1,N+1)]
def lcm(a,b): return a*b//gcd(a,b)
def L(m):
    r=1
    for i in range(1,m+1): r=lcm(r,i)
    return r
tests = {
 'n[1..n]^2':      lambda n: n*L(n)**2,
 'n[1..2n]':       lambda n: n*L(2*n),
 'n[1..n][1..2n]': lambda n: n*L(n)*L(2*n),
 'n^2[1..n][1..2n]': lambda n: n*n*L(n)*L(2*n),
 '[1..n][1..2n]':  lambda n: L(n)*L(2*n),
 'n[1..n][1..3n/2]': lambda n: n*L(n)*L(3*n//2),
}
for name,f in tests.items():
    bad=[n for n in range(1,N+1) if (Fr(f(n))*B[n]).denominator!=1]
    print(f"{name:20s} fails at n =", bad[:6], f"({len(bad)} of {N})")
# c_n sequence
cn=[int(Fr(n*n)*comb(2*n,n)*H[n]) for n in range(1,13)]
print("c_n:", cn)
print("c_n/C(2n,n):", [cn[n-1]/comb(2*n,n) for n in range(1,9)])
print("guess 4^{n-1}? :", [4**(n-1) for n in range(1,9)])
print("c_n - 4 c_{n-1}:", [cn[n]-4*cn[n-1] for n in range(1,12)])
