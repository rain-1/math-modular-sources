"""E12: I_n := I(4n,2n) = a_n + b_n G ; r_n = -a_n/b_n ; measure v_2, v_3 of r_n - r_{n-1}.
   Uses the (verified) E3 pole-lowering recurrence twice:
     (4n,2n) -> (4n-2,2n-2) -> (4n-4,2n-4) = I_{n-1}
   step (m,t)->(m-2,t-2):  I(m,t) = const + coef*I(m-2,t-2),
     const = -1/(2^{t+1} t(t-1) C(2(m-t),m-t)),  coef = (m-1)^2/(4 t(t-1))"""
import importlib.util, time
from fractions import Fraction as F
from math import comb, log
spec=importlib.util.spec_from_file_location("fast","/home/ubuntu/code/math-modular-sources/lattice/sol_notes/05_emn_fast.py")
fast=importlib.util.module_from_spec(spec); spec.loader.exec_module(fast)

def step(m,t):
    c = -F(1, 2**(t+1)*t*(t-1)*comb(2*(m-t),m-t))
    d = F((m-1)**2, 4*t*(t-1))
    return c,d

NMAX=80
from math import comb as _C
vals={}
t0=time.time()
for n in range(1,NMAX+1):
    r,pp,g = fast.I(4*n,2*n)
    assert pp==0
    vals[n]=(r,g)
print("direct exact computation of I(4n,2n) for n=1..%d done in %.1fs"%(NMAX,time.time()-t0))
print("check b_n = C(4n,2n)^2/4^(4n) for all n:", all(vals[n][1]==F(_C(4*n,2*n)**2,4**(4*n)) for n in vals))
print()
def vp(x,p):
    if x==0: return None
    n,d = x.numerator, x.denominator
    e=0
    while n%p==0: n//=p; e+=1
    while d%p==0: d//=p; e-=1
    return e
print(" n |   v_2(r_n-r_{n-1})   v_3(r_n-r_{n-1})   v_5   v_7 |  v2/n")
rows=[]
prev=None
for n in range(1,NMAX+1):
    a,b = vals[n]
    r = -a/b
    if prev is not None:
        d = r-prev
        v2,v3,v5,v7 = vp(d,2),vp(d,3),vp(d,5),vp(d,7)
        rows.append((n,v2,v3,v5,v7))
        if n<=20 or n%10==0:
            print(f"{n:>2} | {v2:>18} {v3:>18} {v5:>5} {v7:>5} | {v2/n:6.3f}")
    prev=r
import numpy as np
ns=np.array([x[0] for x in rows]); v2=np.array([x[1] for x in rows],dtype=float); v3=np.array([x[2] for x in rows],dtype=float)
for lo in (2,10,20,40,60,80):
    msk=ns>=lo
    s2=np.polyfit(ns[msk],v2[msk],1); s3=np.polyfit(ns[msk],v3[msk],1)
    print(f"  LSQ over n>={lo}:  v_2 ~ {s2[0]:.4f} n + {s2[1]:.3f}      v_3 ~ {s3[0]:.4f} n + {s3[1]:.3f}")
print("  v_3 values: min=%d max=%d  mean=%.2f   bounded? %s"%(v3.min(),v3.max(),v3.mean(), v3.max()-v3.min()<50))
print("  full v_3 list:", [int(x) for x in v3])
print("  v_3 max over n>=40:", int(v3[ns>=40].max()), " min:", int(v3[ns>=40].min()))
print("  full v_2 list:", [int(x) for x in v2])
