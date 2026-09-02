import collections
D=[]
for L in open('directions.txt'):
    if not L.startswith('DIR'): continue
    p=L.split()
    D.append(dict(N=int(p[1]),k=int(p[2]),u=int(p[3]),ul=int(p[4]),uo=int(p[5]),uord=int(p[6]),
                  v=int(p[7]),vl=int(p[8]),vo=int(p[9]),vord=int(p[10]),d=int(p[11]),
                  eps=int(p[12]),cls=p[13],nm=p[14]))
cnt=collections.Counter()
for r in D: cnt[(r['N'],r['k'],r['eps'])]+=1
P={}
for L in open('pdims.txt'):
    p=L.split()
    P[(int(p[1]),int(p[2]),int(p[3]))]=int(p[4])
bad=[]
for key,dim in P.items():
    if cnt.get(key,0)!=dim: bad.append((key,cnt.get(key,0),dim))
for key,c in cnt.items():
    if key not in P: bad.append((key,c,'MISSING-in-PARI'))
print("total directions:",len(D))
print("PARI nebentypus rows:",len(P))
print("mismatches:",len(bad))
for b in bad[:40]: print("  ",b)
# total dims per (N,k)
tot=collections.Counter()
for r in D: tot[(r['N'],r['k'])]+=1
ptot=collections.Counter()
for (N,k,m),dim in P.items(): ptot[(N,k)]+=dim
print("total-dim mismatches:", [(x,tot[x],ptot[x]) for x in sorted(set(tot)|set(ptot)) if tot[x]!=ptot[x]])
