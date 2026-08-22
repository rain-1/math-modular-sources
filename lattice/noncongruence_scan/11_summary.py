#!/usr/bin/env python3
"""Merge the five recurrence-class scans, dedupe, and rank."""
import glob, sys, math, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
rows=[]
for f in sorted(glob.glob('out/*_scored.txt')):
    tag=os.path.basename(f).replace('_scored.txt','').replace('eclass_','e=').replace('r3','e=2')
    for line in open(f):
        p=line.split()
        rows.append((float(p[0]),tag,int(p[1]),int(p[2]),int(p[3]),int(p[4]),
                     float(p[5]),float(p[6]),int(p[7]),int(p[8])))
rows.sort(reverse=True)
seen=set(); out=[]
for r in rows:
    key=(r[1],r[2],r[3],r[4],r[5])
    if key in seen: continue
    seen.add(key); out.append(r)
pos=[r for r in out if r[0]>0]
print("total distinct integral non-degenerate rows: %d ; positive score: %d"%(len(out),len(pos)))
print("\n--- positive score, k>=2 (the Beukers-type regime) ---")
print("score      class  (al,ga,de,ze)             lam1          lam2           k")
for r in out:
    if r[0]>0 and r[8]>=2:
        print("%+9.5f  %-5s (%d,%d,%d,%d)%s %-13.7g %-14.8g %d"%(
            r[0],r[1],r[2],r[3],r[4],r[5],' '*max(0,24-len("(%d,%d,%d,%d)"%(r[2],r[3],r[4],r[5]))),r[6],r[7],r[8]))
print("\n--- positive score, k<=1 ---   (count %d)"%len([r for r in out if r[0]>0 and r[8]<2]))
c={}
for r in out:
    if r[0]>0 and r[8]<2: c[(r[1],r[8])]=c.get((r[1],r[8]),0)+1
for k,v in sorted(c.items()): print("   class %s, k=%d : %d rows"%(k[0],k[1],v))
print("\n--- best NEGATIVE-score rows with k=2 (top 14) ---")
n=0
for r in out:
    if r[0]<=0 and r[8]==2:
        print("%+9.5f  %-5s (%d,%d,%d,%d)  lam1=%-12.7g lam2=%-13.8g"%(r[0],r[1],r[2],r[3],r[4],r[5],r[6],r[7]))
        n+=1
        if n>=14: break
