import math, json, ast
rows=[]
for ln in open('02_fricke.out').read().strip().split('\n')[1:]:
    p=ln.split('|')
    if len(p)!=5: continue
    N=int(p[0]); ordv=ast.literal_eval(p[1]); r=ast.literal_eval(p[2]); C=int(p[3]); deg=int(p[4])
    rows.append(dict(N=N,ord=ordv,r=r,C=C,deg=deg))
print("eta-quotient Fricke parameters found:", len(rows))
from collections import Counter
print("by level:", sorted(Counter(r['N'] for r in rows).items()))

def issq(n):
    s=math.isqrt(n); return s*s==n

hits=[]
for h in rows:
    C=h['C']
    if issq(C):
        s=math.isqrt(C)
        for lam2 in (1,-1):
            B=2*s+lam2
            if B<=0: continue
            lam1=B+2*s
            if abs(lam1)<=abs(lam2): continue
            hits.append(dict(**h,B=B,lam1=lam1,lam2=lam2,c=B*B-4*C,field='Q',lam2norm=1.0,rational=True))
    else:
        for cc in (1,-1):
            n=4*C+cc
            if n<0: continue
            B=math.isqrt(n)
            if B*B!=n: continue
            if B<=0: continue
            # lam2 = B - 2 sqrt(C)
            lam2=B-2*math.sqrt(C); lam1=B+2*math.sqrt(C)
            hits.append(dict(**h,B=B,lam1=lam1,lam2=lam2,c=cc,field='Q(sqrt%d)'%C,lam2norm=1.0,rational=False))
print("\nUNIT lambda_2 Fricke hosts (lambda_2^norm = 1):", len(hits))
hits.sort(key=lambda h:(h['N'],h['C'],h['B']))
seen=set()
for h in hits:
    key=(h['N'],h['C'],h['B'])
    tag='' 
    print(f"N={h['N']:4d} C={h['C']:6d} B={h['B']:4d} deg={h['deg']} lam1={h['lam1']:.6g} lam2={h['lam2']:.6g} c={h['c']:4d} field={h['field']:12s} r={h['r']}")
json.dump(hits, open('03_units.json','w'), indent=1)
