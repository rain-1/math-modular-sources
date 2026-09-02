import collections, json
from fractions import Fraction as Fr
D=[]
for L in open('directions.txt'):
    if not L.startswith('DIR'): continue
    p=L.split()
    D.append(dict(N=int(p[1]),k=int(p[2]),u=int(p[3]),ul=int(p[4]),uo=int(p[5]),uord=int(p[6]),
                  v=int(p[7]),vl=int(p[8]),vo=int(p[9]),vord=int(p[10]),d=int(p[11]),
                  eps=int(p[12]),cls=p[13],nm=p[14]))
by=collections.defaultdict(list)
for r in D: by[(r['N'],r['k'])].append(r)
byeps=collections.defaultdict(list)
for r in D: byeps[(r['N'],r['k'],r['eps'])].append(r)

def rank(rows):
    rows=[list(r) for r in rows]; n=len(rows); m=len(rows[0]) if n else 0; rk=0; c=0
    for c in range(m):
        piv=None
        for i in range(rk,n):
            if rows[i][c]!=0: piv=i;break
        if piv is None: continue
        rows[rk],rows[piv]=rows[piv],rows[rk]
        pv=rows[rk][c]
        for i in range(n):
            if i!=rk and rows[i][c]!=0:
                f=Fr(rows[i][c],1)/pv
                for j in range(c,m): rows[i][j]-=f*rows[rk][j]
        rk+=1
        if rk==n: break
    return rk

OUT={}
for (N,k),L in sorted(by.items()):
    d=len(L)
    r1=[Fr(1) if r['u']==1 else Fr(0) for r in L]
    r2=[Fr(1,r['d']**k) if r['v']==1 else Fr(0) for r in L]
    rk=rank([r1,r2]); ann=d-rk
    # interesting period functionals
    per={}
    for r in L:
        if r['cls']=='I': per.setdefault(r['nm'],[]).append(r)
    surv={}
    for nm,S in per.items():
        Ip=[Fr(1,r['d']**(k-1)) if (r in S) else Fr(0) for r in L]
        surv[nm]= (rank([r1,r2,Ip])>rk)
    # elementary block projection
    elemidx=[i for i,r in enumerate(L) if r['cls']=='E']
    elemsurv=None
    if elemidx:
        # does K have nonzero projection on elem block? equivalently rank of [r1,r2] restricted to elem coords < len(elem)?
        sub=[[r1[i] for i in elemidx],[r2[i] for i in elemidx]]
        elemsurv = (len(elemidx) - rank(sub)) > 0
    OUT[(N,k)]=dict(dim=d,rk=rk,ann=ann,surv=surv,elem=elemsurv,nelem=len(elemidx))
GZ=[1,2,3,4,5,6,7,8,9,10,12,13,16,18,25]
print("== genus-zero annihilation (Gamma_1) ==")
print("N  k  dimEis rank ann  interesting periods (survives?)   elem-block-survives")
for N in GZ:
    for k in (3,4):
        if (N,k) not in OUT: 
            print(f"{N:<3}{k:<3}0      -    -")
            continue
        o=OUT[(N,k)]
        s="; ".join(f"{nm}:{'Y' if v else 'N'}" for nm,v in sorted(o['surv'].items()))
        print(f"{N:<3}{k:<3}{o['dim']:<7}{o['rk']:<5}{o['ann']:<5}{s}   | elem({o['nelem']}):{o['elem']}")
json.dump({f"{N},{k}":{'dim':o['dim'],'rk':o['rk'],'ann':o['ann'],'surv':o['surv'],'elem':o['elem'],'nelem':o['nelem']} for (N,k),o in OUT.items()}, open('ann_all.json','w'))
