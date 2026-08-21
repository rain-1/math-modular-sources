#!/usr/bin/env python3
"""Apery limits lim b_n/a_n for every surviving row, and lindep against a
standard constant basis (via PARI)."""
import json, os, sys, subprocess
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.dirname(os.path.abspath(__file__))
import mpmath as mp
mp.mp.dps = 90

def pev(p, n): return sum(Fraction(c)*n**i for i, c in enumerate(p))

def limit_of(poly, a0, a1, n=500):
    p2, p1, p0 = poly
    a = [mp.mpf(a0), mp.mpf(a1)]; b = [mp.mpf(0), mp.mpf(1)]
    for m in range(1, n):
        d = pev(p2, m)
        if d == 0: return None
        d = mp.mpf(d.numerator)/mp.mpf(d.denominator)
        c1 = pev(p1, m); c0 = pev(p0, m)
        c1 = mp.mpf(c1.numerator)/mp.mpf(c1.denominator)
        c0 = mp.mpf(c0.numerator)/mp.mpf(c0.denominator)
        a.append(-(c1*a[m]+c0*a[m-1])/d)
        b.append(-(c1*b[m]+c0*b[m-1])/d)
        if abs(a[-1]) > mp.mpf(10)**200:
            s = mp.mpf(10)**200
            a = [x/s for x in a]; b = [x/s for x in b]
    if a[-1] == 0: return None
    return b[-1]/a[-1]

def main():
    T = json.load(open(os.path.join(HERE, 'table.json')))
    out = []
    for z in T:
        L = None
        try: L = limit_of(z['poly'], z['a'][0], z['a'][1])
        except Exception: pass
        z['limit'] = mp.nstr(L, 60) if L is not None and mp.isfinite(L) else None
        out.append(z)
    json.dump(out, open(os.path.join(HERE,'table.json'),'w'), indent=1)
    # lindep batch
    lines = ["default(realprecision,60);",
             "c2=lfun(-3,2); c3=lfun(-3,3);",
             "B=[1,zeta(2),zeta(3),Catalan,c2,c3,Pi^3,Pi^4];"]
    idx = []
    for i, z in enumerate(out):
        if z['limit'] and abs(float(z['limit'])) > 1e-40:
            idx.append(i)
            lines.append(f'print("{i} ", lindep(concat(B,[{z["limit"]}])));')
    open(os.path.join(HERE,'_lindep.gp'),'w').write("\n".join(lines)+"\n")
    print("wrote _lindep.gp for", len(idx), "rows")

if __name__ == '__main__':
    main()
