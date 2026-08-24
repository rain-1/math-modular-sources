#!/usr/bin/env python3
"""03b_checkmix.py -- independent exact brute force over a small mixed-class box,
to validate 03_fmix (soundness AND completeness)."""
import sys
from fractions import Fraction as Fr
def vfac(n,p):
    s=0;q=p
    while q<=n: s+=n//q; q*=p
    return s
PR=[2,3,5,7,11,13]
NEX=16
NEED=[{p:2*vfac(n,p) for p in PR} for n in range(NEX+2)]
def integral(a,b,c,d,e,f,g,h,j):
    U=[0,0,1]
    for n in range(NEX):
        nx=(a*n*n+b*n+c)*U[n+2]-n*n*(d*n*n+e*n+f)*U[n+1]+n*n*(n-1)**2*(g*n*n+h*n+j)*U[n]
        U.append(nx)
        if nx:
            for p,need in NEED[n+1].items():
                if need:
                    v=0;x=abs(nx)
                    while v<need and x%p==0: x//=p; v+=1
                    if v<need: return False
    return True
RPN,RPD,RRN,RRD,M,J1,J2 = (int(x) for x in sys.argv[1:8])
RMAX,AMAX,CMAX,DMAX,FMAX,GMAX = (int(x) for x in sys.argv[8:14])
rp=Fr(RPN,RPD); rr=Fr(RRN,RRD)
hits=[]
for r in range(-RMAX,RMAX+1):
    if r==0: continue
    for C in range(-GMAX,GMAX+1):
        if C==0: continue
        g=C*M*M
        if g % r: continue
        p=g//r
        for a in range(-AMAX,AMAX+1):
            s=a-r; d=s*r+p
            if abs(d)>DMAX: continue
            b=(1-rr)*r+(1-rp)*s
            e=-rp*(2*p+r*s)-rr*r*s
            h=-(1+2*rp+rr)*Fr(g); jj=Fr(C*J1*J2)
            if any(x.denominator!=1 for x in (b,e,h,jj)): continue
            b,e,h,jj=int(b),int(e),int(h),int(jj)
            disc=18*a*d*g-4*a**3*g+a*a*d*d-4*d**3-27*g*g
            if disc==0: continue
            for c in range(-CMAX,CMAX+1):
                for f in range(-FMAX,FMAX+1):
                    if not integral(a,b,c,d,e,f,g,h,jj): continue
                    # triviality
                    U=[Fr(0),Fr(0),Fr(1)]
                    P=lambda n:a*n*n+b*n+c; Q=lambda n:d*n*n+e*n+f; R=lambda n:g*n*n+h*n+jj
                    for n in range(9): U.append((P(n)*U[n+2]-Q(n)*U[n+1]+R(n)*U[n])/Fr((n+1)**2))
                    u=U[2:]; z=0; triv=False
                    for x in u[1:]:
                        if x==0:
                            z+=1
                            if z>=3: triv=True
                        else: z=0
                    if triv: continue
                    if a==c and b==2*a and d==f and e==2*d and h==2*g and jj==g: continue
                    hits.append((r,a,c,d,f,C))
print(len(hits))
for x in sorted(hits): print(" ".join(map(str,x)))
