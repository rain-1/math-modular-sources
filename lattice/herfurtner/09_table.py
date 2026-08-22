#!/usr/bin/env python3
"""out/rows_full.json -> markdown tables for HERFURTNER_CLASSIFICATION.md"""
import json, math, sys
R=json.load(open('out/rows_full.json'))
KNOWN={
 (1,0,0,7,2,-8):'Zagier A', (1,0,0,9,3,27):'Zagier B', (1,0,0,10,3,9):'Zagier C',
 (1,0,0,11,3,-1):'Zagier D',(1,0,0,12,4,32):'Zagier E',(1,0,0,17,6,72):'Zagier F',
 (2,1,1,136,10,4):'sqrt Apery (Beukers)', (2,1,1,24,2,4):'sqrt T',
 (2,1,1,20,2,16):'sqrt Domb', (2,1,1,72,6,-108):'sqrt AZ(9,3,-27)',
 (2,1,1,88,10,500):'sqrt AZ(11,5,125)', (2,1,1,56,6,324):'sqrt AZ(7,3,81)',
 (3,1,2,26,2,-3):'sqrt Cooper s7',
 (1,0,0,16,4,0):'hypergeometric', 
}
def key(r): return (r['M'],r['j1'],r['j2'],r['A'],r['B'],r['C'])
def qstr(M,j1,j2,C):
    if (M,j1,j2)==(1,0,0): return "%d n^2"%C
    def lin(j):
        if j==0: return "n"
        return "(%dn%+d)"%(M,-j) if M!=1 else "(n%+d)"%(-j)
    return "%d %s%s"%(C,lin(j1),lin(j2))
rows=[r for r in R if not r['casdeg'] and not r['dbl']]
rows.sort(key=lambda r:(-1e9 if r['score'] is None else -r['score']))
print("| # | class $(M;j_1,j_2)$ | $(\\rho;\\delta_\\infty)$ | $P(n)$ | $Q(n)$ | $\\lambda_1$ | $\\lambda_2$ | $k$ | score | $\\mathcal I$ | Herfurtner | $u_n$ | name |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
from fractions import Fraction as F
for i,r in enumerate(rows,1):
    M,j1,j2=r['M'],r['j1'],r['j2']
    rho=F(j1+j2,2*M); dl=F(abs(j2-j1),M)
    P="%d n^2%s%+d"%(r['A'], (" %+d n"%r['be']) if r['be'] else "", r['B'])
    sc = "--" if r['score'] is None else "%+.3f"%r['score']
    l1 = "%.4f"%r['lam1'] if not r['cplx'] else "|%.4f| cplx"%r['lam1']
    l2 = "%.4f"%r['lam2'] if not r['cplx'] else "(same)"
    hf = "; ".join("#%d %s"%(h['row'],h['fib']) for h in r['herf']) or "--"
    nm = KNOWN.get(key(r),'')
    print("| %d | (%d;%d,%d) | (%s;%s) | $%s$ | $%s$ | %s | %s | %d | %s | %s | %s | %s | %s |"%(
        i,M,j1,j2,rho,dl,P,qstr(M,j1,j2,r['C']),l1,l2,r['k'],sc,r['I'],hf,
        ",".join(map(str,r['u'][:6])),nm))
print()
print("degenerate/excluded rows: %d Casoratian-degenerate or double-root, out of %d primitive integral rows"%(
    len([r for r in R if r['casdeg'] or r['dbl']]), len(R)))
