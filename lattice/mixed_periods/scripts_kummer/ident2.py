from decomp import *
from mpmath import pslq, mpf, sqrt, log, pi, mpc, im, re, arg, fabs
from sympy import Rational, factorint
from fractions import Fraction

def as_rat(x, maxden=10**7, tol=None):
    if tol is None: tol = mpf(10)**(-mp.dps+30)
    if abs(x) < tol: return Fraction(0)
    f = Fraction(float(x)).limit_denominator(maxden)
    if abs(mpf(f.numerator)/f.denominator - x) < tol*max(1,abs(x)): return f
    return None

def rat_ratio(x, y, maxden=10**7):
    if abs(y) < mpf(10)**-60: return None
    return as_rat(x/y, maxden)

def gather(h):
    """all algebraic numbers whose logs/args enter the closed form"""
    al=h.alphas(); k=h.k; S=[]
    for j in range(k):
        S.append(-al[j]); S.append(1-al[j])
        for l in range(k):
            if j!=l:
                d=al[j]-al[l]; S.append(d)
                S.append(1-(1-al[l])/d); S.append(1-(-al[l])/d)
                S.append((1-al[l])/d); S.append(-al[l]/d)
    return S

def logbasis(h):
    S=gather(h); P=set()
    for v in S:
        r = as_rat(abs(v)**2, 10**12)
        if r is None: raise RuntimeError(("|v|^2 not rational", v, abs(v)**2))
        for p in list(factorint(r.numerator).keys())+list(factorint(r.denominator).keys()):
            if p>1: P.add(p)
    P=sorted(P)
    return [(f"log{p}", log(p)) for p in P], P

def argbasis(h):
    S=gather(h)
    B=[("pi", pi)]
    for v in S:
        th=arg(v)
        if abs(th)<mpf(10)**-60: continue
        vec=[th]+[b[1] for b in B]
        r=pslq(vec, tol=mpf(10)**(-mp.dps+35), maxcoeff=10**6, maxsteps=100000)
        if r is None or r[0]==0:
            B.append((f"th{len(B)}[arg({mp.nstr(v,10)})]", th))
    return B

def full_basis(h, sqrtd=None):
    LB,P = logbasis(h); AB = argbasis(h)
    gens = LB+AB
    names=[];vals=[]
    for i in range(len(gens)):
        for j in range(i,len(gens)):
            names.append(gens[i][0]+"*"+gens[j][0]); vals.append(gens[i][1]*gens[j][1])
    if sqrtd:
        s=sqrt(sqrtd); base=list(zip(names,vals))
        for nm,v in base: names.append(f"sqrt{sqrtd}*"+nm); vals.append(s*v)
    return names, vals, P, [b[0] for b in AB]

def find(target, names, vals, digits, maxcoeff=10**6):
    v=[target]+list(vals)
    if any(abs(x)<mpf(10)**-80 for x in v): return None
    r=pslq(v, tol=mpf(10)**(-digits), maxcoeff=maxcoeff, maxsteps=500000)
    if r is None or r[0]==0: return None
    c0=r[0]; terms=[(Fraction(-int(c),int(c0)),nm) for nm,c in zip(names,r[1:]) if c!=0]
    resid=sum(mpf(int(c))*x for c,x in zip(r,v))
    return terms, resid
