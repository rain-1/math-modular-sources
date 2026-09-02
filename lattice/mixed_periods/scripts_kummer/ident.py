from decomp import *
from mpmath import pslq, mpf, sqrt, log, pi, mpc, im, re, nstr, fabs
from sympy import Rational, factorint, nsimplify
import itertools, math

def primeset(h):
    """primes appearing in |alpha_j|^2, |1-alpha_j|^2, |alpha_j-alpha_l|^2 (all rational)"""
    al=h.alphas(); k=h.k; S=set()
    vals=[]
    for j in range(k):
        vals.append(abs(al[j])**2); vals.append(abs(1-al[j])**2)
        for l in range(k):
            if j!=l: vals.append(abs(al[j]-al[l])**2)
    for v in vals:
        r = Rational(str(mp.nstr(v,30))).limit_denominator(10**12)
        if abs(mpf(float(r))-v) > 1e-20: raise RuntimeError(("not rational",v))
        for p in list(factorint(r.p).keys())+list(factorint(r.q).keys()):
            if p>1: S.add(p)
    return sorted(S)

def elem_basis(h, with_sqrt=None):
    P = primeset(h)
    names=[]; vals=[]
    names.append("pi^2"); vals.append(pi**2)
    for i,p in enumerate(P):
        for q in P[i:]:
            names.append(f"log{p}*log{q}"); vals.append(log(p)*log(q))
    if with_sqrt:
        s = sqrt(with_sqrt); base=list(zip(names,vals))
        for nm,v in base:
            names.append(f"sqrt{with_sqrt}*{nm}"); vals.append(s*v)
    return names, vals, P

def show_pslq(target, names, vals, tol_digits=None, maxcoeff=10**8, maxsteps=200000):
    if tol_digits is None: tol_digits = mp.dps-20
    v = [target]+list(vals)
    r = pslq(v, tol=mpf(10)**(-tol_digits), maxcoeff=maxcoeff, maxsteps=maxsteps)
    if r is None or r[0]==0: return None
    c0=r[0]
    terms=[]
    for nm,c in zip(names, r[1:]):
        if c!=0: terms.append((Rational(-int(c),int(c0)), nm))
    resid = sum(mpf(int(c))*x for c,x in zip(r,v))
    return terms, resid, r
