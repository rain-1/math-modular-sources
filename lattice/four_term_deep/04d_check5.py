#!/usr/bin/env python3
"""Independent brute-force check of 04_fscan5.c on a small box (exact integer
arithmetic on U_n = u_n (n!)^2, no congruence shortcuts)."""
import sys, subprocess
from fractions import Fraction as Fr
from math import gcd, factorial

RN,RD,M,J1,J2 = (int(x) for x in sys.argv[1:6])
AMAX,CMAX,DMAX,FMAX,GMAX,JMAX,KMAX = (int(x) for x in sys.argv[6:13])
NEX = 16
rho = Fr(RN,RD)
def vfac(n,p):
    s=0;q=p
    while q<=n: s+=n//q; q*=p
    return s
PRIMES=[2,3,5,7,11,13]
NEED=[{p:2*vfac(n,p) for p in PRIMES} for n in range(NEX+2)]

def integral_row(a,b,c,d,e,f,g,h,j,k,l,m):
    U=[0,0,0,1]           # U_{-3..0}
    for n in range(0,NEX):
        P=a*n*n+b*n+c; Q=d*n*n+e*n+f; R=g*n*n+h*n+j; T=k*n*n+l*n+m
        nx=P*U[n+3]-n*n*Q*U[n+2]+n*n*(n-1)*(n-1)*R*U[n+1]-n*n*(n-1)**2*(n-2)**2*T*U[n]
        U.append(nx)
        if nx:
            for p,need in NEED[n+1].items():
                if need:
                    v=0; x=abs(nx)
                    while v<need and x%p==0: x//=p; v+=1
                    if v<need: return False
    return True

hits=[]
for a in range(1,AMAX+1):
  bb=Fr(a)*(1-rho)
  if bb.denominator!=1: continue
  b=int(bb)
  for d in range(-DMAX,DMAX+1):
    ee=-2*rho*Fr(d)
    if ee.denominator!=1: continue
    e=int(ee)
    for g in range(-GMAX,GMAX+1):
      if g==0: continue
      hh=-(1+3*rho)*Fr(g)
      if hh.denominator!=1: continue
      h=int(hh)
      for c in range(-CMAX,CMAX+1):
        for f in range(-FMAX,FMAX+1):
          for j in range(-JMAX,JMAX+1):
            for C in range(-KMAX,KMAX+1):
              if C==0: continue
              k=C*M*M; l=-C*M*(J1+J2); m=C*J1*J2
              if integral_row(a,b,c,d,e,f,g,h,j,k,l,m):
                  hits.append((a,c,d,f,g,j,C))
print(len(hits))
for x in sorted(hits): print(" ".join(map(str,x)))
