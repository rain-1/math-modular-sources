"""(b) BARE central family I_n=I(4n,2n)=a_n+b_n G: archimedean rate (log 64) vs denominator cost."""
import numpy as np, importlib.util
from math import log, lcm
import mpmath as mp
spec=importlib.util.spec_from_file_location("fast","/home/ubuntu/code/math-modular-sources/lattice/sol_notes/05_emn_fast.py")
fast=importlib.util.module_from_spec(spec); spec.loader.exec_module(fast)
LC=[1]
def L(k):
    while len(LC)<=k: LC.append(lcm(LC[-1],len(LC)))
    return LC[k]
print("max_Delta (xy)^4/g^2 = 1/64 (E5) => |I_n| <= G/64^n, rate log 64 = %.4f"%log(64))
rows=[]
for n in range(1,41):
    r,p,g=fast.I(4*n,2*n); assert p==0
    mp.mp.dps=max(60,len(str(r.numerator))+60)
    v=mp.mpf(r.numerator)/mp.mpf(r.denominator)+(mp.mpf(g.numerator)/mp.mpf(g.denominator))*mp.catalan
    den=r.denominator
    import sympy as sp
    kk=max([int(q)**((e+1)//2) for q,e in sp.factorint(den).items()]+[1])
    rows.append((n,float(mp.log(abs(v))),log(den),kk))
    if n%5==0: print("  n=%2d  log|I_n|=%12.5f  log den(a_n)=%9.3f  (1/n)=%7.3f  min k with den|L_k^2 = %d   (4n=%d, 12n=%d)"%(n,rows[-1][1],rows[-1][2],rows[-1][2]/n,kk,4*n,12*n))
ns=np.array([r[0] for r in rows]); lv=np.array([r[1] for r in rows]); ld=np.array([r[2] for r in rows])
m=ns>=10; A=np.vstack([ns[m],np.log(ns[m]),np.ones(m.sum())]).T
s1,*_=np.linalg.lstsq(A,lv[m],rcond=None); s2,*_=np.linalg.lstsq(A,ld[m],rcond=None)
print("  fit log|I_n|     = %.5f n + %.3f log n + %.3f   (-log 64 = %.5f)"%(s1[0],s1[1],s1[2],-log(64)))
print("  fit log den(a_n) = %.5f n + %.3f log n + %.3f"%(s2[0],s2[1],s2[2]))
print("  NET bare family: den cost %.4f/n  vs archimedean gain %.4f/n  ->  %+.4f"%(s2[0],-s1[0],s2[0]+s1[0]))
