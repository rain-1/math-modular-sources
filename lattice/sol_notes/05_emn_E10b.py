import importlib.util, time
from fractions import Fraction as F
from math import comb, lcm, log
import sympy as sp
spec=importlib.util.spec_from_file_location("fast","/home/ubuntu/code/math-modular-sources/lattice/sol_notes/05_emn_fast.py")
fast=importlib.util.module_from_spec(spec); spec.loader.exec_module(fast)
def L(k):
    v=1
    for i in range(1,k+1): v=lcm(v,i)
    return v
print(" n | log den(256^n a_n) | (1/n)log den | min k: den|L_k^2 | den|L_{12n}^2 | largest prime factor | log|256^n L_n| /n")
res=[]
for n in range(1,31):
    t0=time.time()
    r,p,g = fast.Lint(n); assert p==0
    A=F(256)**n*r; B=F(256)**n*g
    den=A.denominator
    fac=sp.factorint(den); pmax=max(fac)
    ok=(L(12*n)**2)%den==0
    kk=1
    while (L(kk)**2)%den: kk+=1
    val = float(A)+float(B)*0.9159655941772190150546035149323841107741
    # careful: huge; use Fraction->mpmath
    import mpmath as mp
    mp.mp.dps=60
    v = mp.mpf(A.numerator)/mp.mpf(A.denominator) + (mp.mpf(B.numerator)/mp.mpf(B.denominator))*mp.catalan
    res.append((n,den,kk,ok,pmax,v))
    print(f"{n:>2} | {log(den):18.4f} | {log(den)/n:12.4f} | {kk:>15} | {str(ok):>13} | {pmax:>20} | {float(mp.log(abs(v))/n):8.4f}  [{time.time()-t0:.1f}s]")
print()
print("least-squares slope of log den vs n (n>=5):")
import numpy as np
ns=np.array([x[0] for x in res if x[0]>=10]); ld=np.array([log(x[1]) for x in res if x[0]>=10])
sl=np.polyfit(ns,ld,1); print("   log den ~ %.4f n + %.3f"%(sl[0],sl[1]))
print("   log(L_{12n}^2)/n -> 24 (PNT).  measured slope:", "%.4f"%sl[0])
print()
print("archimedean side: 256^n L_n ~ ?  max_Delta 256*(xy)^4(1-4g^2)^2/g^2 = 16/27 -> rate log(27/16)=%.4f"%log(27/16))
ns2=np.array([x[0] for x in res if x[0]>=10]); lv=np.array([float(__import__('mpmath').log(abs(x[5]))) for x in res if x[0]>=10])
sl2=np.polyfit(ns2,lv,1); print("   measured log|256^n L_n| ~ %.4f n + %.3f   (i.e. decay rate %.4f)"%(sl2[0],sl2[1],-sl2[0]))
print("   log 64 = %.4f ; log(27/16) = %.4f"%(log(64),log(27/16)))
print()
print("NET: (denominator cost per n) - (archimedean gain per n) = %.4f - %.4f = %.4f  >> 0"%(sl[0],-sl2[0],sl[0]+sl2[0]))
