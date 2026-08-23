#!/usr/bin/env python3
"""lattice/p2_holonomic/parity.py -- Task 2/3 statistics for P2_HOLONOMIC.md.
Reads data/h2_balance200.csv, data/h6_val.csv, data/h1_hermite200.csv."""
import csv, math, sys
import numpy as np
from collections import defaultdict

D = "/home/ubuntu/code/math-modular-sources/lattice/p2_holonomic/data/"
bal = list(csv.DictReader(open(D+"h2_balance200.csv")))
val = {int(r['n']): r for r in csv.DictReader(open(D.replace('data/','data/')+"h6_val.csv"))
       } if False else None
vrows = open(D+"h6_val.csv").read().strip().split("\n")
vhdr = vrows[1].split(",")
val = {int(l.split(",")[0]): dict(zip(vhdr, l.split(","))) for l in vrows[2:]}
h1 = {(round(float(r['k']),1), int(r['n'])): r for r in csv.DictReader(open(D+"h1_hermite200.csv"))}
KS = [22.4, 23.0, 23.9]
by = {k: {int(r['n']): r for r in bal if round(float(r['k']),1)==k} for k in KS}
NS = sorted(by[22.4])

def binom_p(k, n, p=0.5):
    """two-sided exact binomial p-value"""
    from math import comb
    pr = [comb(n,i)*p**n for i in range(n+1)]
    obs = pr[k]
    return sum(x for x in pr if x <= obs*(1+1e-12))

print("="*72); print("TASK 2  --  the parity of the balance index i(n), 4 <= n <= 200")
print("="*72)
par = {k: np.array([by[k][n]['idx'] for n in NS], dtype=int) % 2 for k in KS}
for k in KS:
    p = par[k]; ev = int((p==0).sum())
    print(f"k={k}: even {ev}/{len(p)} = {ev/len(p):.4f}   binomial p = {binom_p(ev,len(p)):.3f}")
    # runs test
    runs = 1 + int((p[1:]!=p[:-1]).sum()); n1=int(p.sum()); n0=len(p)-n1
    mu = 2*n0*n1/len(p)+1; sd = math.sqrt(2*n0*n1*(2*n0*n1-len(p))/(len(p)**2*(len(p)-1)))
    print(f"      runs {runs}, expected {mu:.1f}+-{sd:.1f}, z = {(runs-mu)/sd:+.2f}")
    # longest run of odd (= outside cone)
    best=cur=0
    for x in p:
        cur = cur+1 if x==1 else 0; best=max(best,cur)
    print(f"      longest run of odd i(n) (v1 outside cone): {best}")

print("\n-- autocorrelation of the parity sequence (centred), lags 1..12 --")
for k in KS:
    x = par[k]*2.0-1
    ac = [float(np.dot(x[:-l],x[l:])/(len(x)-l)) for l in range(1,13)]
    z  = [a*math.sqrt(len(x)-l-1) for l,a in zip(range(1,13),ac)]
    print(f"k={k}: " + " ".join(f"{a:+.3f}" for a in ac))
    print(f"      z: " + " ".join(f"{v:+.2f}" for v in z))

print("\n-- dependence on n mod m (chi^2 on the 2xm table), m = 2..12 --")
from scipy import stats as st
for k in KS:
    out=[]
    for m in range(2,13):
        tab=np.zeros((2,m))
        for i,n in enumerate(NS): tab[par[k][i], n%m]+=1
        keep = tab.sum(0)>0
        chi,p_,_,_ = st.chi2_contingency(tab[:,keep]+0.0) if tab[:,keep].min()>=0 else (0,1,0,0)
        out.append(f"m={m}:p={p_:.3f}")
    print(f"k={k}: "+"  ".join(out))

print("\n-- cross-k agreement of the parity (are the three k independent?) --")
for a in range(3):
    for b in range(a+1,3):
        ag = float((par[KS[a]]==par[KS[b]]).mean())
        z = (ag-0.5)*2*math.sqrt(len(NS))
        print(f"  k={KS[a]} vs {KS[b]}: agree {ag:.4f}  z={z:+.2f}")

print("\n-- correlation with 2-adic / arithmetic invariants (point-biserial) --")
cands = {}
cands['s2(n)']       = np.array([bin(n).count('1') for n in NS],float)
cands['s2(3n)']      = np.array([bin(3*n).count('1') for n in NS],float)
cands['v2(n)']       = np.array([(n & -n).bit_length()-1 for n in NS],float)
cands['v2(XU-VY)']   = np.array([float(val[n]['v2mix']) for n in NS])
cands['v2(Xxi-Y)']   = np.array([float(val[n]['v2Zform']) for n in NS])
cands['v2(Vxi-U)']   = np.array([float(val[n]['v2Nform']) for n in NS])
cands['v2(X)']       = np.array([float(val[n]['v2X']) for n in NS])
cands['v2(D^2)']     = np.array([float(val[n]['v2D2']) for n in NS])
cands['L(n)']        = np.array([float(by[22.4][n]['len']) for n in NS])
cands['log h11/n']   = np.array([float(by[22.4][n]['logh11']) for n in NS])
cands['log h22/n']   = np.array([float(by[22.4][n]['logh22']) for n in NS])
cands['a_{i+1}']     = np.array([float(by[22.4][n]['pqnext']) for n in NS])
cands['margin']      = np.array([float(by[22.4][n]['marg']) for n in NS])
for k in KS:
    print(f"  k={k}:")
    for nm,v in cands.items():
        r,p_ = st.pointbiserialr(par[k], v)
        flag = "  <== " if p_<0.01 else ""
        print(f"     {nm:14s} r={r:+.3f}  p={p_:.3f}{flag}")

print("\n-- exact-law candidates: does parity(i(n)) equal parity of X ? --")
laws = {
 'L(n)':            np.array([int(by[22.4][n]['len']) for n in NS]),
 'n':               np.array(NS),
 'floor(k n)':      np.array([int(22.4*n) for n in NS]),
 's2(3n)':          np.array([bin(3*n).count('1') for n in NS]),
 'v2(XU-VY)':       np.array([int(val[n]['v2mix']) for n in NS]),
 'L(n)-i(n)':       np.array([int(by[22.4][n]['len'])-int(by[22.4][n]['idx']) for n in NS]),
}
for k in KS:
    for nm,v in laws.items():
        m = int((par[k]==(v%2)).sum())
        if max(m, len(NS)-m) > 0.60*len(NS):
            print(f"  k={k} parity(i) vs parity({nm}): match {m}/{len(NS)}")
    # against the other k's parities
print("  (nothing printed above 60% match means: no candidate law)")

print("\n" + "="*72); print("TASK 3 -- v_1(n) = (c_Z, c_N): size and 'entropy'")
print("="*72)
for k in KS:
    q = np.array([float(by[k][n]['logh11']) for n in NS])
    ii = np.array([int(by[k][n]['idx']) for n in NS],float)
    LL = np.array([int(by[k][n]['len']) for n in NS],float)
    r = ii/LL
    print(f"k={k}: i(n)/L(n) = {r.mean():.4f} +- {r.std():.4f}; "
          f"L(n)/n = {np.polyfit(np.array(NS,float),LL,1)[0]:.3f}; "
          f"i(n)/n = {np.polyfit(np.array(NS,float),ii,1)[0]:.3f}")
print("\nInformation content of v_1 given the rows: v_1 is the convergent at index i(n),")
print("so specifying it costs log2 L(n) bits;  L(200) ~ %d  ->  %.1f bits." %
      (int(np.mean([int(by[k][200]['len']) for k in KS])),
       math.log2(np.mean([int(by[k][200]['len']) for k in KS]))))
