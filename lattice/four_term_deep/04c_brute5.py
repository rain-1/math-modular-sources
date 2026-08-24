#!/usr/bin/env python3
"""Brute-force cross-check of 04_fscan5.c on a tiny box (exact rational arithmetic)."""
import sys, subprocess, itertools
from fractions import Fraction as Fr
sys.path.insert(0,'/home/ubuntu/code/math-modular-sources/lattice/four_term_deep')
import importlib.util
spec=importlib.util.spec_from_file_location("r5","/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/05b_report5.py")
r5=importlib.util.module_from_spec(spec); spec.loader.exec_module(r5)

RN,RD,M,J1,J2 = (int(x) for x in sys.argv[1:6])
AMAX,CMAX,DMAX,FMAX,GMAX,JMAX,KMAX = (int(x) for x in sys.argv[6:13])
NEX = 20
rho = Fr(RN,RD)
found=[]
for a in range(1,AMAX+1):
  if (Fr(a)*(1-rho)).denominator!=1: continue
  for d in range(-DMAX,DMAX+1):
    if (-2*rho*Fr(d)).denominator!=1: continue
    for g in range(-GMAX,GMAX+1):
      if g==0: continue
      if (-(1+3*rho)*Fr(g)).denominator!=1: continue
      for c in range(-CMAX,CMAX+1):
        for f in range(-FMAX,FMAX+1):
          for j in range(-JMAX,JMAX+1):
            for C in range(-KMAX,KMAX+1):
              if C==0: continue
              co=r5.coeffs(RN,RD,M,J1,J2,a,c,d,f,g,j,C)
              u=r5.seq(co,NEX)
              if not all(x.denominator==1 for x in u): continue
              # degeneracies as in the scanner
              k=co['k']
              disc_ok=True
              import mpmath as mp
              lam=mp.polyroots([1,-a,d,-g,k],maxsteps=400,extraprec=400)
              if min(abs(lam[i]-lam[q]) for i in range(4) for q in range(i+1,4))<1e-25: disc_ok=False
              if not disc_ok: continue
              zr=0; trivial=False
              for x in u[1:12]:
                  if x==0:
                      zr+=1
                      if zr>=4: trivial=True
                  else: zr=0
              if trivial: continue
              if a==c and co['b']==2*a and d==f and co['e']==2*d and g==j and co['h']==2*g and co['k']==co['m'] and co['l']==2*co['k']: continue
              found.append((a,c,d,f,g,j,C))
print("BRUTE hits:", len(found))
for x in sorted(found): print("  ",x)
