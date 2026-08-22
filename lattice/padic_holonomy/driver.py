"""Per-target contour optimisation driver.  usage: python3 driver.py <index> [which]"""
import os, sys as _sys
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in _sys.path: _sys.path.insert(0, HERE)
import sys, json, time, math, warnings
import numpy as np
warnings.filterwarnings('ignore')
import optimise, family
from targets import TARGETS

idx = int(sys.argv[1])
which = sys.argv[2] if len(sys.argv) > 2 else 'RE'
T = TARGETS[idx]
key = T['key']
print("### %s   d=%d m=%d logR=%.6f tau=%.6f budget=%.6f" % (key, T['d'], T['m'], T['L'], T['tau'], T['budget']), flush=True)
t0 = time.time()
out = []


def record(tag, specs):
    r = optimise.evaluate(T['H'], specs, 8000, 'both')
    if r is None:
        return
    for w in ('RE', 'BC'):
        V = r[w]; lg = math.log(r['dr'])
        den = lg + T['L'] - T['tau']
        num = V + T['L']
        rec = dict(tag=tag, which=w, specs=[[k, list(p)] for k, p in specs], val=V, logdr=lg,
                   dr=r['dr'], cost=V - T['m'] * lg, margin=T['m'] * den - num,
                   bound=(num / den if den > 1e-12 else None), gmax=r['gmax'], gmin=r['gmin'])
        out.append(rec)
        print("  %-14s %-2s val=%9.5f logdr=%+8.5f cost=%9.5f margin=%+9.5f bound=%s [%.0fs]" %
              (tag, w, V, lg, rec['cost'], rec['margin'],
               ("%9.5f" % rec['bound']) if rec['bound'] else "  n/a  ", time.time() - t0), flush=True)


for shape in ('S1', 'S2', 'S3', 'S4', 'O1', 'O2', 'O3'):
    b = optimise.run(T, shape, which, N=1200, starts=4, seed=11)
    if b[1] is not None:
        record(shape, b[1])
for ang, cdep in optimise.good_angles(key):
    for shape in ('T1', 'T2', 'T3', 'T4'):
        b = optimise.run(T, shape, which, N=1200, starts=4, seed=13, tang_angle=ang)
        if b[1] is not None:
            record('%s@%.4f' % (shape, ang), b[1])

json.dump(dict(key=key, d=T['d'], m=T['m'], L=T['L'], tau=T['tau'], budget=T['budget'], runs=out),
          open('res_%02d_%s.json' % (idx, which), 'w'), indent=1)
print("done", time.time() - t0, flush=True)
