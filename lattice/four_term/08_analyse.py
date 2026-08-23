#!/usr/bin/env python3
"""Deduplication and family detection for the four-term scan.

Reads scan hit files (RN RD M J1 J2 a c d f C per line), runs the fast report,
and sorts the survivors into
  (i)   REPEATED   : repeated characteristic root (cubic discriminant 0) -- not five points
  (ii)  DISGUISED  : an apparent singularity at some t_i or at infinity  -- four points,
                     i.e. a gauge / cusp-move image of a three-term row
  (iii) RESCALED   : u_n -> mu^n u_n image of another integral row in the box
  (iv)  CANDIDATE  : genuinely five singular points, primitive
and looks for one-parameter families among the candidates of each class.
"""
import sys, json, os, glob
from fractions import Fraction as Fr
from math import gcd
import importlib.util
spec = importlib.util.spec_from_file_location("rep", os.path.join(os.path.dirname(os.path.abspath(__file__)),"05_report.py"))
rep = importlib.util.module_from_spec(spec); spec.loader.exec_module(rep)

def is_integral(cls, row, N=60):
    try:
        co = rep.coeffs(*cls, *row)
    except Exception:
        return False
    u = rep.seq(co, N)
    return all(x.denominator == 1 for x in u)

def reduce_scale(cls, row):
    """largest mu>=2 with (a/mu, c/mu, d/mu^2, f/mu^2, C/mu^3) integral AND still an integral row."""
    a,c,d,f,C = row
    best = 1
    for mu in range(2, 40):
        if a % mu or c % mu or d % mu**2 or f % mu**2 or C % mu**3: continue
        r2 = [a//mu, c//mu, d//mu**2, f//mu**2, C//mu**3]
        if r2[4] == 0: continue
        if is_integral(cls, r2): best = mu
    return best

def canon(cls, row):
    """t -> -t : (a,c,d,f,C) -> (-a,-c,d,f,-C).  Normalise to a>0 (or C>0 if a=0)."""
    a,c,d,f,C = row
    if a < 0 or (a == 0 and C < 0):
        return [-a,-c,d,f,-C]
    return list(row)

def main(paths, out_json):
    seen = set(); recs = []
    for p in paths:
        for line in open(p):
            v = line.split()
            if len(v) != 10: continue
            v = tuple(map(int, v))
            if v in seen: continue
            seen.add(v)
            recs.append(v)
    print(f"# {len(recs)} distinct hits", file=sys.stderr)
    out = []
    for v in recs:
        cls, row = list(v[:5]), list(v[5:])
        r = rep.analyse(*cls, *row, full=False)
        if r is None:
            r = dict(cls=cls, row=row, verdict='NOT-INTEGRAL'); out.append(r); continue
        if r['disc'] == 0: r['verdict'] = 'REPEATED'
        elif any(r['apparent']) or r['inf_apparent']: r['verdict'] = 'DISGUISED'
        else:
            mu = reduce_scale(cls, row)
            r['verdict'] = f'RESCALED(mu={mu})' if mu > 1 else 'CANDIDATE'
            r['mu'] = mu
        out.append(r)
    json.dump(out, open(out_json,'w'))
    from collections import Counter
    cnt = Counter((tuple(r['cls']), r['verdict'].split('(')[0]) for r in out)
    for k in sorted(cnt): print(k, cnt[k], file=sys.stderr)
    return out

if __name__ == '__main__':
    paths = sys.argv[2:] if len(sys.argv)>2 else sorted(glob.glob('out/*.txt'))
    main(paths, sys.argv[1] if len(sys.argv)>1 else 'out/analysis.json')
