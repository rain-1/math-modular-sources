import importlib.util, time
from fractions import Fraction as F
from math import comb, lcm, log, gcd
import sympy as sp
spec=importlib.util.spec_from_file_location("fast","/home/ubuntu/code/math-modular-sources/lattice/sol_notes/05_emn_fast.py")
fast=importlib.util.module_from_spec(spec); spec.loader.exec_module(fast)

def L(k): 
    v=1
    for i in range(1,k+1): v=lcm(v,i)
    return v

print("E10: den(256^n a_n) | L_{12n}^2 ?   L_k = lcm(1..k)")
print()
rows=[]
NMAX=8
for n in range(0,NMAX+1):
    t0=time.time()
    r,p,g = fast.Lint(n)
    assert p==0
    A = F(256)**n * r
    den = A.denominator
    fac = sp.factorint(den)
    # smallest k with den | lcm(1..k)^2
    kk=1
    while True:
        if (L(kk)**2) % den == 0: break
        kk+=1
        if kk>4000: kk=None; break
    # smallest k with den | lcm(1..k)
    k1=1
    while True:
        if L(k1) % den == 0: break
        k1+=1
        if k1>20000: k1=None; break
    rows.append((n,A,den,fac,kk,k1,time.time()-t0))
    print(f"n={n}: 256^n a_n = {A}")
    print(f"      den = {den}")
    print(f"      factorisation = { {int(q):int(e) for q,e in fac.items()} }")
    print(f"      log(den) = {log(den):.4f}   (1/n)log(den) = {log(den)/n if n else float('nan'):.5f}")
    print(f"      smallest k with den | lcm(1..k)^2 : {kk}   (claim: 12n = {12*n})   -> divides L_{{12n}}^2 : {(L(12*n)**2)%den==0 if n else True}")
    print(f"      smallest k with den | lcm(1..k)   : {k1}")
    print(f"      [{time.time()-t0:.1f}s]")
    print()
print("summary: n, log(den(256^n a_n)), (1/n)log(den), min k for lcm(1..k)^2, min k for lcm(1..k), 12n")
for n,A,den,fac,kk,k1,tt in rows:
    print(f"  {n:>2}  {log(den):12.4f}  {(log(den)/n if n else 0):9.4f}   {str(kk):>6} {str(k1):>6}   {12*n:>4}")
