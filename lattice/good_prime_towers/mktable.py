#!/usr/bin/env python3
import sys,re
lam={};agr={};udf={};bdf={}
for fn in sys.argv[1:]:
    for line in open(fn):
        t=line.split()
        if not t: continue
        if t[0]=="LAM": lam[(t[1],int(t[2]),int(t[3]))]=(int(t[4]),int(t[5]),int(t[6]),int(t[7]),t[8])
        if t[0]=="AGR": agr[(t[1],int(t[2]),int(t[3]))]=" ".join(t[4:])
        if t[0]=="UDF": udf[(t[1],int(t[2]),int(t[3]))]=" ".join(t[4:])
        if t[0]=="BDF": bdf[(t[1],int(t[2]),int(t[3]))]=" ".join(t[4:])
name={"G":"gamma (Apery zeta(3))","C":"C","D":"D (Apery zeta(2))"}
print("| row | $p$ | $a$ | $\\chi(p)$ | $v_p(\\Lambda_a)$ | digits | $\\Lambda_a\\cdot p^{-v}$ (base-$p$, l.s.d. first) | $v_p(\\rho^A_s-1)$ | $v_p(\\chi p^{w}\\rho^B_s-1)$ |")
print("|---|---|---|---|---|---|---|---|---|")
for k in sorted(lam, key=lambda k:(k[1],k[0],k[2])):
    if k[2]>3: continue
    w,chi,v,nd,U=lam[k]
    U=int(U); p=k[1]
    digs=[]
    x=U
    for i in range(min(nd,10)):
        digs.append(x%p); x//=p
    print(f"| {name.get(k[0],k[0])} | {p} | {k[2]} | {chi} | {v} | {nd} | {','.join(map(str,digs))}{'...' if nd>10 else ''} | {udf.get(k,'')} | {bdf.get(k,'')} |")
