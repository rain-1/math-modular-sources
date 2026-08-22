#!/usr/bin/env python3
"""Cheap pass over ALL scan hits: characteristic roots only.
A positive score log(1/|lam_2|) - k needs |lam_2| < e^{-k} <= e^{-1} (k >= 1 for
any row with an irrational Apery limit).  So list every hit with |lam_2| < 1."""
import glob, os, math
from collections import defaultdict
per=defaultdict(int); decay=[]
for p in sorted(glob.glob('out/c_*.txt'))+['out/a0.txt']:
    if not os.path.exists(p): continue
    for line in open(p):
        f=line.split()
        if len(f)!=6: continue
        M,j1,j2,A,B,C=map(int,f)
        if (j1%M==0 and j1//M>=1) or (j2%M==0 and j2//M>=1): continue   # Casoratian degenerate
        D=C*M*M; disc=A*A-4*D
        per[(M,j1,j2)]+=1
        if disc<0:  l2=math.sqrt(abs(D))
        else:
            r1=(A+math.sqrt(disc))/2.0; r2=(A-math.sqrt(disc))/2.0
            l2=min(abs(r1),abs(r2))
        if l2<1.0: decay.append((l2,M,j1,j2,A,B,C))
decay.sort()
print("total hits:",sum(per.values()))
print("hits with |lambda_2| < 1 :",len(decay))
for d in decay[:80]:
    print("  |lam2|=%.6f  class(%d;%d,%d) A=%d B=%d C=%d"%d)
