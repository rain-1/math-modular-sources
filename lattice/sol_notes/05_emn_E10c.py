"""E10 addendum: high-precision |256^n L_n| (huge cancellation between a_n and b_n G),
   denominator growth, and the net cost/gain balance."""
import importlib.util, time
from fractions import Fraction as F
from math import log, lcm
import mpmath as mp, numpy as np, sympy as sp
spec=importlib.util.spec_from_file_location("fast","/home/ubuntu/code/math-modular-sources/lattice/sol_notes/05_emn_fast.py")
fast=importlib.util.module_from_spec(spec); spec.loader.exec_module(fast)
rows=[]
print(" n |  log|256^n L_n|  | local slope | log den | (1/n)log den")
prev=None
for n in range(1,41):
    r,p,g = fast.Lint(n); assert p==0
    A=F(256)**n*r; B=F(256)**n*g
    mp.mp.dps = max(60, len(str(A.numerator))+len(str(A.denominator))+60)
    v = mp.mpf(A.numerator)/mp.mpf(A.denominator) + (mp.mpf(B.numerator)/mp.mpf(B.denominator))*mp.catalan
    lv = float(mp.log(abs(v)))
    den=A.denominator; ld=log(den)
    sl = lv-prev if prev is not None else float('nan')
    prev=lv
    rows.append((n,lv,ld))
    if n<=5 or n%5==0: print(f"{n:>2} | {lv:15.6f} | {sl:11.5f} | {ld:9.2f} | {ld/n:8.4f}")
ns=np.array([r[0] for r in rows]); lv=np.array([r[1] for r in rows]); ld=np.array([r[2] for r in rows])
for lo in (10,20,30):
    m=ns>=lo
    s1=np.polyfit(ns[m],lv[m],1); s2=np.polyfit(ns[m],ld[m],1)
    print(f"  n>={lo}: log|256^n L_n| ~ {s1[0]:.4f} n + {s1[1]:.2f}   |   log den ~ {s2[0]:.4f} n + {s2[1]:.2f}")
# fit with a log n correction: log|.| = alpha n + beta log n + c
m=ns>=10
Amat=np.vstack([ns[m], np.log(ns[m]), np.ones(m.sum())]).T
sol,*_=np.linalg.lstsq(Amat,lv[m],rcond=None)
print(f"  fit log|256^n L_n| = {sol[0]:.5f} n + {sol[1]:.3f} log n + {sol[2]:.3f}   ; log(16/27) = {log(16/27):.5f}")
sol2,*_=np.linalg.lstsq(Amat,ld[m],rcond=None)
print(f"  fit log den(256^n a_n) = {sol2[0]:.5f} n + {sol2[1]:.3f} log n + {sol2[2]:.3f}  ; log(L_12n^2)/n -> 24")
print()
print("NET per n:  denominator cost %.4f  vs archimedean gain %.4f  ->  net %+0.4f (must be < 0 for irrationality)"%(sol2[0],-sol[0],sol2[0]+sol[0]))
