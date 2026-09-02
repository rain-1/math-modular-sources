"""Score every host of this sweep with the CDT_FINDER formulas (m=14, CDT inventory)."""
import math, sys, os
sys.path.insert(0, '/home/ubuntu/code/math-modular-sources/lattice/cdt_finder')
from cdt_bound import tau_flat, tau_sharp

LOSS = math.log(0.6292232680); BC = 11.845
def tau_for(m,k,u):
    sm,tf = tau_flat(m, [(uj,2) for uj in u]); ec = round(6*m/14.0)
    ts,_ = tau_sharp(m, [1]*ec+[0]*(m-ec)); return float(tf)+ts
def score(l2n,k,m=14):
    s = 1.0/l2n
    u = [max(0,min(m,round((2*j-1)*m/14.0))) for j in range(1,k+1)]
    T = tau_for(m,k,u); ceil_ = math.log(256*s); real = ceil_+LOSS; Nbc = BC+math.log(s)
    return dict(tau=T, ceil=ceil_, entryC=ceil_-T, entryR=real-T, margin=m*(real-T)-Nbc, u=u)

def l2norm(lam1,lam2,c,rational):
    return abs(lam2) if rational else math.sqrt(abs(c))

# ---- FAMILY I : Fricke / weight-2 form / third-order rows (Theorem 3.4 shape) ----
# N, C, B, deg(u), lam1+lam2, lam1*lam2, freeint, k, period, tag
F1 = [
 (5,125,22,1, 44,  -16, 1,2,'zeta(2)/10','new to ledger'),
 (6, 81,14,2, 28, -128, 1,2,'zeta(2)/8','new to ledger'),
 (6, 72,17,1, 34,    1, 0,3,'zeta(3)/6','Apery'),
 (6, 64,20,2, 40,  144, 1,2,'L(2,chi_-3)/4','new to ledger'),
 (7, 49,13,1, 26,  -27, 1,2,'zeta(2)/7','Cooper s_7'),
 (8, 32,12,1, 24,   16, 0,3,'7 zeta(3)/32','AZ eps (12,4,16)'),
 (8, 16,24,2, 48,  512, 1,2,'G/4  (Catalan)','new to ledger'),
 (9, 27, 9,1, 18,  -27, 0,3,'L(3,chi_-3)/3','AZ zeta (9,3,-27)'),
 (10,25, 6,2, 12,  -64, 1,2,'zeta(2)/5','Cooper s_10'),
 (12, 9,10,2, 20,   64, 0,3,'7 zeta(3)/24','Domb alpha'),
 (12, 1,34,4, 68, 1152, 1,2,'5 L(2,chi_-3)/16','new to ledger'),
 (18, 1,14,4, 28,  192, 1,2,'L(2,chi_-3)/2','Cooper s_18'),
]
rows=[]
for (N,C,B,deg,s,c,fi,k,per,tag) in F1:
    sq = (int(math.isqrt(C))**2==C)
    l1 = B+2*math.sqrt(C); l2 = B-2*math.sqrt(C)
    l2n = l2norm(l1,l2,c,sq)
    fld = 'Q' if sq else 'Q(sqrt%d)'%(C if not sq else 0)
    rows.append(dict(fam='Fricke w=2', name=f'N={N} C={C} B={B}', lvl=N, deg=deg, w=2, k=k,
                     lam1=l1, lam2=l2, c=c, l2n=l2n, field=fld, period=per, tag=tag, freeint=fi,
                     sing=4, farfold=True))

# ---- FAMILY II : weight-one rows on the four-point Gamma_0(N)/Gamma_1(N) ----
r5 = math.sqrt(5)
F2 = [
 ('Gamma_0(5), pole at cusp 0','5',None,None,-22,125,'complex pair','L(3,chi_5)/2 (AZ eta)',math.sqrt(125),'Q(i sqrt..)'),
 ('Gamma_0(6), pole at cusp 1/2 = Zagier C','6',9,1,10,9,'','L(2,chi_-3)/2  [CDT]',1.0,'Q'),
 ('Gamma_0(6), pole at cusp 1/3 = Zagier A','6',8,-1,7,-8,'','zeta(2)/4',1.0,'Q'),
 ('Gamma_0(6), pole at cusp 0  = Zagier F','6',9,8,17,72,'','5 L(2,chi_-3)/8',8.0,'Q'),
 ('Gamma_0(7), pole at cusp 0','7',None,None,-13,49,'complex pair','(complex fold)',7.0,'Q(sqrt-3)'),
 ('Gamma_0(8), pole at cusp 0 or 1/2 = Zagier E','8',8,4,12,32,'','G = Catalan',4.0,'Q'),
 ('Gamma_0(8), pole at cusp 1/4','8',4,-4,0,-16,'two dominant','(degenerate)',4.0,'Q'),
 ('Gamma_0(9), pole at cusp 0 = Zagier B','9',None,None,-9,27,'complex pair','(complex fold)',math.sqrt(27),'Q(sqrt-3)'),
 ('Gamma_1(5) = Zagier D','5',(11+5*r5)/2,(11-5*r5)/2,11,-1,'','zeta(2)/5',1.0,'Q(sqrt5) unit'),
]
for (nm,lvl,l1,l2,s,c,note,per,l2n,fld) in F2:
    rows.append(dict(fam='weight-one', name=nm, lvl=lvl, deg=1, w=1, k=2, lam1=l1, lam2=l2, c=c,
                     l2n=l2n, field=fld, period=per, tag=note, freeint=0, sing=4, farfold=False))

for r in rows:
    r.update(score(r['l2n'], r['k']))
rows.sort(key=lambda r: -r['entryR'])
hdr = f"{'host':46s}{'w':>2s}{'k':>2s}{'lam2^norm':>11s}{'field':>15s}{'tau':>8s}{'ceil':>8s}{'entryC':>8s}{'entryR':>8s}{'margin':>9s}  period"
print(hdr); print('-'*len(hdr))
for r in rows:
    print(f"{r['name']:46s}{r['w']:>2d}{r['k']:>2d}{r['l2n']:>11.4f}{r['field']:>15s}"
          f"{r['tau']:>8.3f}{r['ceil']:>8.3f}{r['entryC']:>+8.3f}{r['entryR']:>+8.3f}{r['margin']:>+9.3f}  {r['period']}")
print()
print("m needed for a positive margin (entryR>0 hosts):")
for r in rows:
    if r['entryR']>0:
        # margin = m*entryR - (BC + log s); solve
        s = 1.0/r['l2n']; need = (BC+math.log(s))/r['entryR']
        print(f"  {r['name']:46s} entryR={r['entryR']:+.4f}  m >= {need:.2f}")
import json; json.dump(rows, open('11_scored.json','w'), indent=1, default=str)
