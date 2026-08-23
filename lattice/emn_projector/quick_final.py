from fractions import Fraction as Fr
from math import comb, log
from emn_core import series_H
N=1600
H=series_H(N)
# (a) note's closed form a_n = (1/(2^{n+1}(n+1))) sum_k C(n,k)/(2k+1), h_{n+1}=a_n
ok=True
for n in range(0,60):
    a=Fr(1,2**(n+1)*(n+1))*sum(Fr(comb(n,k),2*k+1) for k in range(n+1))
    if a!=H[n+1]: ok=False;print("MISMATCH",n)
print("note's a_n formula matches the ODE recursion for n<60:", ok)
# (b) sigma_B at larger n
Ht=[H[n]*Fr(2)**n for n in range(N+1)]
conv=[Fr(0)]*(N+1); acc=Fr(0)
for n in range(N+1):
    acc=acc+Ht[n]; conv[n]=acc
for n in (400,800,1200,1600):
    Bn=conv[n-1]/Fr(n)
    print(f"n={n}: log den(B_n)/n = {log(Bn.denominator)/n:.4f}   log den(h_n)/n = {log(H[n].denominator)/n:.4f}")
