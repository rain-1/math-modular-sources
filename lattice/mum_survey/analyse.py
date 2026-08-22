#!/usr/bin/env python3
"""Post-process survey output: apply the Euler-factor law
      r_infty = r_p * E_p(m),   E_p(m) = 1 - chi(p) p^{-m}
   to turn a p-adic Apery limit into a prediction for the archimedean one,
   and check consistency across primes and against AvSZ's tables."""
import re, sys
from fractions import Fraction as F

def kron(D,a):
    # Kronecker symbol (D/a) for fundamental discriminant D
    from math import gcd
    if gcd(D,a)!=1: return 0
    # use sympy-free implementation
    def js(a,n):
        a%=n; t=1
        while a:
            while a%2==0:
                a//=2
                if n%8 in (3,5): t=-t
            a,n=n,a
            if a%4==3 and n%4==3: t=-t
            a%=n
        return t if n==1 else 0
    if a%2==1: return js(D%a,a) if a>0 else 0
    # a even: use multiplicativity
    r=1; m=a
    while m%2==0:
        m//=2
        if D%8 in (1,7): pass
        elif D%8 in (3,5): r=-r
        else: return 0
    return r*(js(D%m,m) if m>1 else 1)

CHI = {'1':(1,lambda a:1 if a%1==0 else 0),
       'chi-4':(4,lambda a: 0 if a%2==0 else (1 if a%4==1 else -1)),
       'chi-3':(3,lambda a: 0 if a%3==0 else (1 if a%3==1 else -1))}
def chival(name,p):
    if name=='1': return 1
    if name=='chi-4': return 0 if p%2==0 else (1 if p%4==1 else -1)
    if name=='chi-3': return 0 if p%3==0 else (1 if p%3==1 else -1)
    m=re.match(r'chi(-?\d+)$',name)
    if m: return kron(int(m.group(1)),p)
    return None
ARCHNAME={('1',3):'zeta(3)',('1',2):'pi^2/6 = zeta(2)',('chi-4',2):'G',
          ('chi-3',2):'L(chi-3,2)',('chi-3',3):'L(chi-3,3)'}

def parse(fn):
    rows=[]
    for l in open(fn):
        if l.startswith('#') or '|' not in l: continue
        f=l.rstrip().split('|')
        if len(f)<10 or f[3]=='SKIP': continue
        rows.append(f)
    return rows

def padic_entries(s):
    import ast
    t = re.sub(r'(-?\d+)/(\d+)', r'"@\1/\2"', s)
    try: v = ast.literal_eval(t)
    except Exception: return []
    def num(x):
        if isinstance(x,str) and x.startswith('@'): return F(x[1:])
        return F(x)
    out=[]
    for e in v:
        out.append((int(e[0]), int(e[1]), [(h[0], num(h[1])) for h in e[2]]))
    return out

def main(fn):
    rows=parse(fn)
    print("AESZ | nn | dz | k | lam1 | lam2 | arch(measured) | p-adic data | arch(predicted from p-adics)")
    for f in rows:
        pd=padic_entries(f[9])
        if not pd: continue
        preds={}
        desc=[]
        for p,dig,hits in pd:
            for name,r in hits:
                mm=re.match(r'L_p\((\d),(.*)om\^-\d\)$',name)
                if not mm: 
                    if name=='zeta_p(3)': m_,chi=3,'1'
                    else: continue
                else:
                    m_=int(mm.group(1)); chi=mm.group(2) or '1'
                    if chi=='': chi='1'
                cv=chival(chi,p)
                if cv is None: continue
                E=F(1)-F(cv)*F(1,p**m_)
                preds.setdefault((chi,m_),[]).append((p,r*E))
                desc.append(f"p={p}(sig,{dig}d): {r}*{name}")
        if not preds: continue
        pstr=[]
        for (chi,m_),lst in preds.items():
            vals=set(v for _,v in lst)
            tag="CONSISTENT" if len(vals)==1 else "CONFLICT"
            an=ARCHNAME.get((chi,m_), f"L({chi},{m_})")
            pstr.append(f"{list(vals)[0] if len(vals)==1 else vals} * {an} [{tag}, primes {[p for p,_ in lst]}]")
        print(f"{f[0]:<10}|{f[1]:<9}|{f[2]}|{f[3]}|{float(f[4]):.4g}|{float(f[5]):.4g}| {f[7]}*{f[6]} | {'; '.join(desc)} | {'; '.join(pstr)}")

if __name__=='__main__': main(sys.argv[1] if len(sys.argv)>1 else 'survey_N240.txt')
