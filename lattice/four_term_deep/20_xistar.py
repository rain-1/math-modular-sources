#!/usr/bin/env python3
"""20_xistar.py -- the constant xi* = 0.7372929961855962... to high precision, from
two different rows in two different worlds, as a cross-check."""
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 320
# (i) four-term mixed row  (1/2,0;1,1,1), r=-1, (a,c,d,f,C)=(13,1,-13,-1,-1)
rp,rr = Fr(1,2),Fr(0); M,J1,J2,r,a,c,d,f,C = 1,1,1,-1,13,1,-13,-1,-1
g=C*M*M; p=g//r; s=a-r
b=int((1-rr)*r+(1-rp)*s); e=int(-rp*(2*p+r*s)-rr*r*s)
h=int(-(1+2*rp+rr)*Fr(g)); j=int(Fr(C*J1*J2))
P=lambda n:a*n*n+b*n+c; Q=lambda n:d*n*n+e*n+f; R=lambda n:g*n*n+h*n+j
lam=sorted(mp.polyroots([1,-a,d,-g],maxsteps=400,extraprec=400),key=lambda z:-abs(z))
N=int(200/mp.log10(abs(lam[0])/abs(lam[1])))+400
A=[mp.mpf(0),mp.mpf(0),mp.mpf(1)]; B=[mp.mpf(0),mp.mpf(0),mp.mpf(0),mp.mpf(1)]
for n in range(N):
    A.append((P(n)*A[n+2]-Q(n)*A[n+1]+R(n)*A[n])/mp.mpf((n+1)**2))
    if n>=1: B.append((P(n)*B[n+2]-Q(n)*B[n+1]+R(n)*B[n])/mp.mpf((n+1)**2))
x1=+(B[N+2]/A[N+2])
print("four-term  (1/2,0;1,1,1) r=-1 (13,1,-13,-1,-1)  n=%d"%N)
print("  xi =", mp.nstr(x1,200))
# (ii) five-term six-point row (0;1,1,1), (a,c,d,f,g,j,C)=(5,2,0,0,-20,-8,-16)
a5,c5,d5,f5_,g5,j5,C5 = 5,2,0,0,-20,-8,-16
b5=a5; e5=0; h5=-g5; k5=C5; l5=-C5*2; m5=C5
P5=lambda n:a5*n*n+b5*n+c5; Q5=lambda n:d5*n*n+e5*n+f5_
R5=lambda n:g5*n*n+h5*n+j5;  T5=lambda n:k5*n*n+l5*n+m5
lam5=sorted(mp.polyroots([1,-a5,d5,-g5,k5],maxsteps=500,extraprec=500),key=lambda z:-abs(z))
N5=int(200/mp.log10(abs(lam5[0])/abs(lam5[1])))+400
A=[mp.mpf(0)]*3+[mp.mpf(1)]; B=[mp.mpf(0)]*4+[mp.mpf(1)]
for n in range(N5):
    A.append((P5(n)*A[n+3]-Q5(n)*A[n+2]+R5(n)*A[n+1]-T5(n)*A[n])/mp.mpf((n+1)**2))
    if n>=1: B.append((P5(n)*B[n+3]-Q5(n)*B[n+2]+R5(n)*B[n+1]-T5(n)*B[n])/mp.mpf((n+1)**2))
x2=+(B[N5+3]/A[N5+3])
print("five-term  (0;1,1,1) (5,2,0,0,-20,-8,-16)  n=%d"%N5)
print("  xi =", mp.nstr(x2,200))
print()
print("agreeing digits:", int(-mp.log10(abs(x1-x2)/abs(x1))) if x1!=x2 else "exact")
open('out/xistar.txt','w').write("xistar_fourterm %s\nxistar_fiveterm %s\n"%(mp.nstr(x1,180),mp.nstr(x2,180)))
