#!/usr/bin/env python3
"""Post-process hits.json: normalisation, characteristic roots, budget,
p-adic slopes, Apery limits (lindep via PARI)."""
import json, os, sys, math, subprocess
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from search import ptrim

HERE = os.path.dirname(os.path.abspath(__file__))
ZAG2 = {(7,2,-8),(9,3,27),(10,3,9),(11,3,-1),(12,4,32),(17,6,72)}
ZAG3 = {(7,3,81),(9,3,-27),(10,4,64),(11,5,125),(12,4,16),(17,5,1)}

def lead(p):
    p = ptrim([Fraction(x) for x in p]); return len(p)-1, p[-1]

def roots(poly):
    (d2,c2),(d1,c1),(d0,c0) = [lead(p) for p in poly]
    k = max(d2,d1,d0)
    A = float(c2) if d2==k else 0.0
    B = float(c1) if d1==k else 0.0
    C = float(c0) if d0==k else 0.0
    if A == 0: return None, k
    disc = B*B-4*A*C
    if disc < 0: return ('cplx', math.sqrt(abs(C/A))), k
    s = math.sqrt(disc)
    r1, r2 = (-B+s)/(2*A), (-B-s)/(2*A)
    return (r1, r2), k

def apery_limit(poly, prec=400, n=600):
    """b_n/a_n with a0=1,a1=?; use the recurrence itself with b0=0,b1=1."""
    import mpmath as mp
    mp.mp.dps = prec
    def pev(p, x): return sum(Fraction(c)*x**i for i,c in enumerate(p))
    p2,p1,p0 = poly
    # a: need a0=1, a1 from n=0 row if p0(0)=0
    if pev(p2,0) == 0: return None
    a1 = -pev(p1,0)/pev(p2,0)
    a = [Fraction(1), a1]; b = [Fraction(0), Fraction(1)]
    for m in range(1, n):
        d = pev(p2,m)
        if d == 0: return None
        a.append(-(pev(p1,m)*a[m]+pev(p0,m)*a[m-1])/d)
        b.append(-(pev(p1,m)*b[m]+pev(p0,m)*b[m-1])/d)
    if a[-1] == 0: return None
    r = b[-1]/a[-1]
    return mp.mpf(r.numerator)/mp.mpf(r.denominator)

CONSTS = ['1','zeta(2)','zeta(3)','Catalan','L2chi-3','L3chi-3','Pi^2','Pi^3','Pi^4']
def lindep(x, dps=200):
    gpsrc = f"""
default(realprecision,{dps});
x = {x};
L2=sumalt(n=0,(-1)^n/(3*n+1)^2)*0; \\\\ placeholder
Lchi2 = lfun(-3,2); Lchi3 = lfun(-3,3);
v = [1, zeta(2), zeta(3), Catalan, Lchi2, Lchi3, Pi^2, Pi^3, Pi^4, x];
print(lindep(v));
"""
    p = os.path.join(HERE,'_ld.gp'); open(p,'w').write(gpsrc)
    try:
        out = subprocess.run(['gp','-q',p], capture_output=True, text=True, timeout=60).stdout.strip()
    except Exception: return None
    return out

def main():
    hits = json.load(open(os.path.join(HERE,'hits.json')))
    rows = []
    for h in hits:
        poly = h['poly']
        rr, k = roots(poly)
        cls = h['cls']
        clt = None
        if cls:
            o = int(cls[0]); a,b,c = [Fraction(x) for x in cls[1:]]
            clt = (o, a, b, c)
        if rr is None: continue
        if rr[0] == 'cplx':
            Lam = rr[1]; kind='complex'
        else:
            Lam = max(abs(rr[0]), abs(rr[1])); kind='real'
        budget = math.log(Lam) - k if Lam > 0 else None
        known = None
        if clt:
            trip = tuple(int(x) if Fraction(x).denominator==1 else float(x) for x in clt[1:])
            if clt[0]==2 and trip in ZAG2: known='Zagier-2nd'
            elif clt[0]==3 and trip in ZAG3: known='AZ-3rd'
            else: known='NOT-on-list'
        rows.append(dict(N=h['N'], r=h['r'], s=h['s'], w=h['w'], flip=h['flip'],
                         a=h['a'][:8], cls=[str(x) for x in clt] if clt else None,
                         degp2=k, roots=(rr if rr[0]!='cplx' else ['cplx',rr[1]]),
                         Lam=Lam, budget=budget, known=known, poly=poly,
                         nalso=len(h.get('also',[]))))
    rows.sort(key=lambda z: -(z['budget'] if z['budget'] is not None else -99))
    json.dump(rows, open(os.path.join(HERE,'table.json'),'w'), indent=1)
    for z in rows:
        print(f"N={z['N']:<3} w={z['w']} deg={z['degp2']} cls={z['cls']} Lam={z['Lam']:.5f} budget={z['budget']:+.4f} {z['known']} a={z['a'][:6]}")

if __name__ == '__main__':
    main()
