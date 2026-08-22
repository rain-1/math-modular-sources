import os, sys as _sys
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in _sys.path: _sys.path.insert(0, HERE)
import numpy as np, math
import haupt, bcint, opt
from targets import TARGETS

N = 3000
for T in TARGETS:
    H = T['H']; m = T['m']
    print("="*100)
    print("%-26s d=%d m=%d logR=%.6f tau=%.6f budget=%.6f" % (T['key'], T['d'], m, T['L'], T['tau'], T['budget']))
    rows = []
    # family: plain discs
    for r in (0.3,0.4,0.5,0.6,0.65,0.7,0.75,0.8):
        psi, dr = opt.build('disc', (r,))
        res = opt.eval_contour(H, psi, dr, N)
        for use in ('RE','BC'):
            o = opt.report(T,res,use)
            rows.append(('disc r=%.2f'%r, use, o))
    # family: lune (bite at q ~ +r, i.e. towards the cusp 0 at q=1)
    for r in (0.5,0.6,0.7,0.8,0.9):
        for c in (1.5,2.0,2.5,4.0):
            psi, dr = opt.build('lune',(r,c,0.0))
            res = opt.eval_contour(H, psi, dr, N)
            for use in ('RE','BC'):
                rows.append(('lune r=%.2f c=%.1f'%(r,c), use, opt.report(T,res,use)))
    best = {}
    for nm,use,o in rows:
        k=(use,'margin')
        if k not in best or o['margin']>best[k][1]['margin']: best[k]=(nm,o)
        k2=(use,'bound')
        if k2 not in best or o['bound']<best[k2][1]['bound']: best[k2]=(nm,o)
    for k in sorted(best):
        nm,o = best[k]
        print("  best %-3s by %-6s : %-22s  val=%9.5f logdr=%8.5f  bound=%9.5f  margin=%9.5f  (m=%d)"%
              (k[0],k[1],nm,o['val'],o['logdr'],o['bound'],o['margin'],m))
