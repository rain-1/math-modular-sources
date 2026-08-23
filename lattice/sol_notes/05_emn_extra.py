"""Extras: (a) refined minimax quadratic for E8; (b) archimedean rate log 64 for the BARE central
   family I(4n,2n) and its denominator cost, for the cost-vs-gain comparison in E10."""
import numpy as np, importlib.util
from numpy.polynomial import polynomial as Pn
from fractions import Fraction as F
from math import log, lcm
import mpmath as mp
spec=importlib.util.spec_from_file_location("fast","/home/ubuntu/code/math-modular-sources/lattice/sol_notes/05_emn_fast.py")
fast=importlib.util.module_from_spec(spec); spec.loader.exec_module(fast)

def sup_abs(coefs):
    ph=np.zeros(2*len(coefs)-1); ph[0::2]=coefs
    best=0.0; barg=0.0
    A=Pn.polymul([0,0,16.0],ph); dA=Pn.polyder(A)
    rts=[r.real for r in np.roots(dA[::-1]) if abs(r.imag)<1e-9] if len(dA)>1 else []
    for x in [0.0,0.5]+[r for r in rts if 0<r<0.5]:
        v=abs(Pn.polyval(x,A));  best,barg=(v,x) if v>best else (best,barg)
    N=Pn.polymul(16*np.array(Pn.polypow([1.0,-1.0],4)),ph)
    num=Pn.polysub(Pn.polymul(Pn.polyder(N),[0,1.0]),2*N)
    rts=[r.real for r in np.roots(num[::-1]) if abs(r.imag)<1e-9]
    for x in [0.5,1-1e-12]+[r for r in rts if 0.5<r<1.0]:
        v=abs(Pn.polyval(x,N)/x**2); best,barg=(v,x) if v>best else (best,barg)
    return best,barg
from scipy.optimize import minimize
print("(a) refined minimax over phi(h)=1+u h+v h^2 (Nelder-Mead, several starts):")
best=(1e9,None)
for st in [(-10,23),(-10,22.2),(-9,20),(-12,30),(-8,16)]:
    r=minimize(lambda z: sup_abs([1.0,z[0],z[1]])[0], st, method='Nelder-Mead',
               options=dict(xatol=1e-12,fatol=1e-14,maxiter=20000,maxfev=20000))
    if r.fun<best[0]: best=(r.fun,r.x)
print("    min_{u,v} 256 sup_g W(g)|1+u g^2+v g^4| = %.10f  at (u,v)=(%.8f, %.8f)"%(best[0],best[1][0],best[1][1]))
print("    q=(1,-10,23) gives %.10f ; (1-4h)^2=(1,-8,16) gives %.10f"%(sup_abs([1,-10,23])[0], sup_abs([1,-8,16])[0]))

print()
print("(b) BARE central family I_n = I(4n,2n) = a_n + b_n G")
print("    max_Delta (xy)^4/g^2 = 1/64 (E5)  =>  |I_n| <= G/64^n, rate log 64 = %.4f"%log(64))
def L(k):
    v=1
    for i in range(1,k+1): v=lcm(v,i)
    return v
rows=[]
for n in range(1,41):
    r,p,g=fast.I(4*n,2*n); assert p==0
    mp.mp.dps=max(60,len(str(r.numerator))+len(str(r.denominator))+60)
    v=mp.mpf(r.numerator)/mp.mpf(r.denominator)+(mp.mpf(g.numerator)/mp.mpf(g.denominator))*mp.catalan
    den=r.denominator
    rows.append((n,float(mp.log(abs(v))),log(den),den))
    if n%5==0:
        kk=1
        while (L(kk)**2)%den: kk+=1
        print("    n=%2d  log|I_n| = %12.5f   log den(a_n) = %10.4f   (1/n)=%8.4f   min k: den|L_k^2 = %d  (12n=%d, 4n=%d)"%(n,rows[-1][1],rows[-1][2],rows[-1][2]/n,kk,12*n,4*n))
ns=np.array([r[0] for r in rows]); lv=np.array([r[1] for r in rows]); ld=np.array([r[2] for r in rows])
m=ns>=10
A=np.vstack([ns[m],np.log(ns[m]),np.ones(m.sum())]).T
s1,*_=np.linalg.lstsq(A,lv[m],rcond=None); s2,*_=np.linalg.lstsq(A,ld[m],rcond=None)
print("    fit log|I_n|      = %.5f n + %.3f log n + %.3f   (log 64 = %.5f)"%(s1[0],s1[1],s1[2],-log(64)))
print("    fit log den(a_n)  = %.5f n + %.3f log n + %.3f"%(s2[0],s2[1],s2[2]))
print("    NET bare family: %.4f - %.4f = %+.4f"%(s2[0],-s1[0],s2[0]+s1[0]))
