#!/usr/bin/env python3
import json, os, re
from fractions import Fraction
HERE=os.path.dirname(os.path.abspath(__file__))
NAMES=['1','zeta(2)','zeta(3)','G','L(2,chi-3)','L(3,chi-3)','pi^3','pi^4']
T=json.load(open(os.path.join(HERE,'table.json')))
ident={}
for line in open(os.path.join(HERE,'lindep.out')):
    m=re.match(r'^(\d+)\s+\[(.*)\]~',line.strip())
    if not m: continue
    i=int(m.group(1)); v=[int(x) for x in m.group(2).split(',')]
    nz=[j for j in range(8) if v[j]]
    if v[8]==0 or len(nz)!=1: continue
    if max(abs(x) for x in v)>10000: continue
    j=nz[0]
    ident[i]=f"{Fraction(-v[j],v[8])}*{NAMES[j]}"
def etaq(D,e): return "".join(f"e{d}^{x} " for d,x in zip(D,e) if x)
def divisors(n): return [d for d in range(1,n+1) if n%d==0]
rows=[]
for i,z in enumerate(T):
    D=divisors(z['N'])
    rows.append((i,z,ident.get(i)))
out=[]
out.append("| # | N | t = prod eta(d)^r_d | F = prod eta(d)^s_d | w | ord p2 | (a,b,c) | roots | Lambda | budget | status | Apery limit |")
out.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
for i,z,idt in rows:
    D=divisors(z['N'])
    rr=z['roots']
    rs = "complex" if rr[0]=='cplx' else f"{rr[0]:.4f}, {rr[1]:.4f}"
    cls = "(%s,%s,%s)"%tuple(z['cls'][1:]) if z['cls'] else "-"
    out.append(f"| {i} | {z['N']} | `{etaq(D,z['r'])}` | `{etaq(D,z['s'])}` | {z['w']} | {z['degp2']} | {cls} | {rs} | {z['Lam']:.4f} | {z['budget']:+.4f} | {z['known'] or 'non-Zagier-normal'} | {idt or (z['limit'][:22] if z['limit'] else '-')} |")
open(os.path.join(HERE,'table.md'),'w').write("\n".join(out)+"\n")
print("\n".join(out))
