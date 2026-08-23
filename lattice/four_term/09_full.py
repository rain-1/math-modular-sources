#!/usr/bin/env python3
"""Full invariants for the CANDIDATE rows: companion b_n, sharp k, Apery limit xi,
score log(1/|lambda_2|)-k, Moebius invariants I_1 = a^3/g, I_2 = a^2/d."""
import sys, json, os
from fractions import Fraction as Fr
import importlib.util
spec = importlib.util.spec_from_file_location("rep", os.path.join(os.path.dirname(os.path.abspath(__file__)),"05_report.py"))
rep = importlib.util.module_from_spec(spec); spec.loader.exec_module(rep)
import mpmath as mp
mp.mp.dps = 160

def run(anal_json, out_json, NXI=900):
    A = json.load(open(anal_json))
    cands = [r for r in A if r.get('verdict','').startswith('CANDIDATE')]
    out = []
    for r in cands:
        cls, row = r['cls'], r['row']
        try:
            g = rep.analyse(*cls, *row, full=True, NXI=NXI)
        except Exception as ex:
            g = dict(cls=cls,row=row,ERR=str(ex))
        a,c,d,f,C = row
        M = cls[2]; gg = C*M*M
        g['I1'] = str(Fr(a**3, gg)) if gg else None
        g['I2'] = str(Fr(a**2, d)) if d else None
        out.append(g)
    json.dump(out, open(out_json,'w'), indent=0)
    return out

if __name__ == '__main__':
    out = run(sys.argv[1], sys.argv[2])
    for r in out:
        if 'ERR' in r: print(r); continue
        lam = ', '.join(f"{x[0]:.5g}{('%+.5gi'%x[1]) if abs(x[1])>1e-9 else ''}" for x in r['lam'])
        sc = f"{r['score']:+.4f}" if r['score'] is not None else '  -  '
        print(f"{str(r['cls']):<18}{str(r['row']):<26} k={r['k']} score={sc}  I1={r['I1']:<12} I2={r['I2']:<10} lam=[{lam}]")
        print(f"    u = {r['u']}   xi = {r['xi']}")
