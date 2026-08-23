"""Fast exact engine for I(m,t) = int_Delta (xy)^m/(1-x^2-y^2)^{t+1} dxdy.
   Values are triples (rat, pic, gc) meaning rat + pic*pi + gc*G.  Same derivation as 05_emn_core.py."""
from fractions import Fraction as F
from math import comb, lcm
import functools, sys
sys.setrecursionlimit(100000)

def add(a,b): return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
def smul(c,a): return (c*a[0], c*a[1], c*a[2])

@functools.lru_cache(maxsize=None)
def Z(k):
    if k<=0:
        n=-k
        return (sum(F(comb(n,j),2*j+1) for j in range(n+1)), F(0), F(0))
    if k==1: return (F(0), F(1,4), F(0))
    j=k-1
    p=Z(j)
    return (F(1,2*j*2**j)+F(2*j-1,2*j)*p[0], F(2*j-1,2*j)*p[1], F(0))

@functools.lru_cache(maxsize=None)
def X(m,p):
    q=m+p; assert q>=0
    r=pi_=F(0)
    for i in range(q+1):
        c=F(comb(q,i)*2**(q-i)*(-1)**i)
        z=Z(m+1-i); r+=c*z[0]; pi_+=c*z[1]
    return (r,pi_,F(0))

@functools.lru_cache(maxsize=None)
def W(k):
    if k==0: return (F(0),F(1,2),F(0))
    if k==1: return (F(1),F(0),F(0))
    p=W(k-2); return (F(k-1,k)*p[0], F(k-1,k)*p[1], F(0))

@functools.lru_cache(maxsize=None)
def A(k):
    if k==0: return (F(0),F(0),F(-1))
    return (-F(1,2*k)*sum(F((-1)**i,2*i+1) for i in range(k)), F(0), F(0))

@functools.lru_cache(maxsize=None)
def Sm(m):
    assert m%2==0
    r=m//2
    acc=smul(F(comb(2*r,r)), A(0))
    for k in range(1,r+1):
        acc=add(acc, smul(F(2*(-1)**k*comb(2*r,r-k)), A(k)))
    return smul(-F(2,4**r), acc)

@functools.lru_cache(maxsize=None)
def J(m):
    acc=(F(0),F(0),F(0))
    for i in range(m):
        acc=add(acc, smul(F((-1)**i,i+1), W(m-1-i)))
    acc=add(acc, Sm(m))     # m even so (-1)^m = +1
    return smul(F(1,2), acc)

@functools.lru_cache(maxsize=None)
def I(m,t):
    tot=(F(0),F(0),F(0))
    X0=X(m,0)
    for j in range(m+1):
        if j==t: continue
        c=F(comb(m,j)*2**(m-j)*(-1)**j, j-t)
        tot=add(tot, smul(c, add(smul(F(2)**(j-t), X0), smul(F(-1), X(m,j-t)))))
    if 0<=t<=m:
        tot=add(tot, smul(F(comb(m,t)*(-1)**t)*F(2)**(m-t), J(m)))
    return smul(F(2)**(t+1)*F(1,2)*F(1,4**m), tot)

def Lint(n):
    """L_n = int_Delta (xy)^{4n} (1-4g^2)^{2n} / g^{2n+1} dxdy,  g=1-x^2-y^2.
       (1-4g^2)^{2n} = sum_k C(2n,k)(-4)^k g^{2k}  =>  pole exponent t+1 = 2n+1-2k, i.e. t = 2n-2k."""
    tot=(F(0),F(0),F(0))
    for k in range(2*n+1):
        tot=add(tot, smul(F(comb(2*n,k)*(-4)**k), I(4*n, 2*n-2*k)))
    return tot
