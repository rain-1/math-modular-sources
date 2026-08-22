#!/usr/bin/env python3
"""Aggregate the modular scan: rank by score = log(1/lam2) - k (k=2 by Thm R3,
measured where lam2 < 1).  Report eta exponent vectors, lambda, congruence flag."""
import json, glob, math, os, sys
os.chdir(os.path.dirname(os.path.abspath(__file__)))
d=json.load(open('eta_pairs.json'))
rows=[]
for f in glob.glob('out/mod/j*.txt'):
    for line in open(f):
        p=line.split()
        if len(p)!=10: continue
        N=int(p[0]); ti=int(p[1]); fi=int(p[2]); lam=int(p[3])
        o=int(p[4]); dg=int(p[5]); k=int(p[6])
        l1=float(p[7]); l2=float(p[8])
        v=d[str(N)]
        tr=v['t'][ti-1]['r']; tdeg=v['t'][ti-1]['deg']; fr=v['F'][fi-1]
        kk = k if k>=0 else 2
        sc = math.log(1/abs(l2))-kk if abs(l2)>1e-12 else -99
        rows.append(dict(N=N,t=tr,tdeg=tdeg,F=fr,lam=lam,order=o,deg=dg,k=k,
                         l1=l1,l2=l2,score=sc,
                         Feven=all(x%2==0 for x in fr)))
rows.sort(key=lambda r:-r['score'])
# de-duplicate identical (lam1,lam2,lam,order,deg) invariants, keep smallest level
seen={}; uniq=[]
for r in rows:
    key=(round(r['l1'],9),round(r['l2'],9),r['lam'],r['order'],r['deg'])
    if key in seen: seen[key]+=1; continue
    seen[key]=1; uniq.append(r)
for r in uniq: r['mult']=seen[(round(r['l1'],9),round(r['l2'],9),r['lam'],r['order'],r['deg'])]
json.dump(uniq, open('out/mod_table.json','w'), indent=0)
print(f"{len(rows)} rows, {len(uniq)} distinct invariant classes\n")
print("score    k lam  ord/deg  lam1          lam2          N   deg(t)  Feven  mult  t                     F")
for r in uniq[:35]:
    if r['l2']>=1: break
    print("%+8.4f %2d %3d  %d/%d   %-13.6g %-13.6g %-3d %-6d %-5s %-5d %-20s %s"%(
        r['score'],r['k'],r['lam'],r['order'],r['deg'],r['l1'],r['l2'],r['N'],r['tdeg'],
        r['Feven'],r['mult'],r['t'],r['F']))
print("\n--- all classes with lam2 < 1 (decaying linear form) ---")
for r in uniq:
    if r['l2']<1 and abs(r['l2'])>1e-12:
        print("%+8.4f k=%d lam=%d N=%-3d l1=%-12.6g l2=%-12.6g t=%-22s F=%s"%(
            r['score'],r['k'],r['lam'],r['N'],r['l1'],r['l2'],r['t'],r['F']))
