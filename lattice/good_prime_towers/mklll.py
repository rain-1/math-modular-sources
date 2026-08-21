#!/usr/bin/env python3
"""Read LAM lines from log files, emit a gp script of identification tests."""
import sys, collections
lam={}
for fn in sys.argv[1:]:
    for line in open(fn):
        t=line.split()
        if t and t[0]=="LAM":
            _,row,p,a,w,chi,v,nd,U = t[:9]
            lam[(row,int(p),int(a))]=(int(w),int(chi),int(v),int(nd),int(U))
out=[]
for (row,p,a),(w,chi,v,nd,U) in sorted(lam.items()):
    if nd < 6: continue
    out.append(f'L = {p}^({v})*({U} + O({p}^{nd}));')
    out.append(f'print("== {row} p={p} a={a} w={w} chi={chi} vLam={v} digits={nd}");')
    out.append(f'test("Q-rational      ",[1,L],{nd},{p});')
    out.append(f'test("deg<=2 algebraic",[1,L,L^2],{nd},{p});')
    out.append(f'test("deg<=3 algebraic",[1,L,L^2,L^3],{nd},{p});')
    ap = {"G":"AP4","C":"AP3","F":"AP3","B":"AP3"}.get(row)
    if ap:
        k = w+1 if row!="G" else 4   # weight: S_4 for G, S_3 for the chi_-3 rows
        out.append(f'ur = unitroot(mapget({ap},{p}),{p},{k-1},{nd});')
        out.append(f'if(ur, test("(1,unitroot,Lam) ",[1,ur,L],{nd},{p}), print("    unit root: NON-ORDINARY at p={p} (a_p=",mapget({ap},{p}),")"));')
# cross-tower pairs
for (row,p,a),(w,chi,v,nd,U) in sorted(lam.items()):
    for (row2,p2,a2),(w2,chi2,v2,nd2,U2) in sorted(lam.items()):
        if row2!=row or p2!=p or a2<=a: continue
        K=min(nd,nd2)
        if K<6: continue
        out.append(f'L1 = {p}^({v})*({U} + O({p}^{K})); L2 = {p}^({v2})*({U2} + O({p}^{K}));')
        out.append(f'print("== cross {row} p={p} a={a} vs a={a2}");')
        out.append(f'test("(1,Lam_a,Lam_a\')",[1,L1,L2],{K},{p});')
out.append("quit")
print("\n".join(out))
