#!/usr/bin/env python3
"""ap_analyse2.py -- refined statistics on the exact 5-adic valuations."""
import sys, math
from collections import Counter
p=sys.argv[1]; rows={}; NN=None
for line in open(p):
    line=line.strip()
    if line.startswith('#'):
        if line.startswith('# NN='): NN=int(line.split()[1][3:])
        continue
    f=line.split()
    if len(f)==5: rows[int(f[0])]=tuple(int(z) for z in f[1:])
top=int(0.8*NN)
d={n:rows[n][3]-3*n for n in range(1,top+1)}
print("=== deficit d_n = v5(h_n) - 3n, windowed extremes (n<=%d) ==="%top)
print(" window        min d_n   (at n)    max d_n   (at n)   -1.73*log(nmid)")
w=[(2,25),(26,50),(51,100),(101,200),(201,400),(401,800),(801,1600)]
for lo,hi in w:
    if lo>top: break
    hi=min(hi,top)
    sub={n:d[n] for n in range(lo,hi+1)}
    mn=min(sub.values()); mx=max(sub.values())
    amn=[n for n in sub if sub[n]==mn][0]; amx=[n for n in sub if sub[n]==mx][0]
    print(f" [{lo:4d},{hi:4d}]  {mn:6d} ({amn:5d})  {mx:6d} ({amx:5d})   {-1.73*math.log((lo*hi)**0.5):8.2f}")
print("\n#(d_n = k) histogram over 2<=n<=%d:"%top)
c=Counter(d[n] for n in range(2,top+1))
print("  ", sorted(c.items()))
print("\n=== distribution of v5(b_n), 1<=n<=%d ==="%NN)
cb=Counter(rows[n][1] for n in range(1,NN+1))
print("  v : count :", sorted(cb.items()))
print("  (a 'random' 5-adic sample would give counts ~ 1600,320,64,13,3,...)")
print("\n  v5(b_n) by n mod 4:", {r:round(sum(rows[n][1] for n in range(1,NN+1) if n%4==r)/len([n for n in range(1,NN+1) if n%4==r]),3) for r in range(4)})
print("  v5(b_n) by n mod 5:", {r:round(sum(rows[n][1] for n in range(1,NN+1) if n%5==r)/len([n for n in range(1,NN+1) if n%5==r]),3) for r in range(5)})
print("  #(v5(b_n)=0) by n mod 4:", {r:sum(1 for n in range(1,NN+1) if n%4==r and rows[n][1]==0) for r in range(4)})
print("\n  v5(a_n)==v5(b_n) for all 1<=n<=%d ? "%NN, all(rows[n][0]==rows[n][1] for n in range(1,NN+1)))
print("  max |v5(e_n)-v5(b_n)| =", max(abs(rows[n][2]-rows[n][1]) for n in range(1,NN+1)))
