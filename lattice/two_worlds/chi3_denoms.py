"""Exact denominator laws for the alpha=2 (a=2m,b=m) chi_{-3} row.
   Q_m = 27*(-256)^m*(2m)!*(1/6)_m^2*(5/6)_m^2 / ((m!)^4 (2/3)_m (4/3)_m)
   v_3(Q_m) = -3m + 2 s_3(m) - s_3(2m)/2 + 3      [proved from the closed form]
"""
import sys
from fractions import Fraction as F
sys.path.insert(0,'/home/ubuntu/code/math-modular-sources/lattice/two_worlds')
from chi3_row import row
from math import gcd

def s3(n):
    s=0
    while n: s+=n%3; n//=3
    return s
def vp(x,p):
    n,d=x.numerator,x.denominator; v=0
    while n%p==0: n//=p; v+=1
    while d%p==0: d//=p; v-=1
    return v
def Qclosed(m):
    r=F(27)
    for j in range(m):
        r *= F(-32*(2*j+1)*(6*j+1)**2*(6*j+5)**2, 9*(j+1)**3*(3*j+2)*(3*j+4))
    return r
def D(n):
    r=1
    for k in range(2,n+1): r = r*k//gcd(r,k)
    return r

print(" m  Q==closed  v3(Q) predicted  den'(P) | D_m^2 ? | D_{2m} ? | D_m^3?")
for m in range(1,26):
    Q,P=row(m,2*m)
    pred = -3*m + 2*s3(m) - s3(2*m)//2 + 3
    dP=P.denominator
    while dP%3==0: dP//=3
    print("%3d  %s  %5d %5d   %6s %6s %6s" % (m, Q==Qclosed(m), vp(Q,3), pred,
        D(m)**2 % dP == 0, D(2*m) % dP == 0, D(m)**3 % dP == 0))
