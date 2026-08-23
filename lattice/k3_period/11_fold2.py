#!/usr/bin/env python3
"""11_fold2.py -- same as 10_fold.py but (i) real-axis continuation to the
I_6 fold at t=1, (ii) both I_3 folds, (iii) a precision self-check."""
import sys
from mpmath import mp, mpf, mpc, log, sqrt, nstr
import importlib.util
spec = importlib.util.spec_from_file_location("f0", "/home/ubuntu/code/math-modular-sources/lattice/k3_period/10_fold.py")

DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 230
mp.dps = DPS
NS  = int(sys.argv[2]) if len(sys.argv) > 2 else 2200

import types
src = open('/home/ubuntu/code/math-modular-sources/lattice/k3_period/10_fold.py').read()
src = src.split("if __name__ ==")[0]
src = src.replace("mp.dps = 230", "mp.dps = %d" % DPS)
mod = types.ModuleType("f0"); mod.__dict__['__name__'] = 'f0'
exec(compile(src, '10_fold.py', 'exec'), mod.__dict__)

start_series = mod.start_series
logcoeffs    = mod.logcoeffs

t1 = mpc('0.1', 0)
ua, ub, Av, Ad, Bv, Bd = start_series(t1, NS)
print("tail of A-series at t1:", nstr(abs(ua[-1]*t1**(len(ua)-1)), 8))

sq2 = sqrt(2)
tpm = (mpc(5) - mpc(0,1)*sq2)/27
tpp = (mpc(5) + mpc(0,1)*sq2)/27

res = {}
# --- I_3 fold at t_+ = (5 - i sqrt2)/27, path through the lower half plane
c0,c1,d0,d1 = logcoeffs(tpm, tpm - mpc(0,'0.05'), [mpc('0.1','-0.05'), tpm - mpc(0,'0.05')], t1, Av, Ad, Bv, Bd)
res['xi_I3_minus'] = d1/c1
res['c1_I3_minus'] = c1
# --- conjugate fold
c0,c1,d0,d1 = logcoeffs(tpp, tpp + mpc(0,'0.05'), [mpc('0.1','0.05'), tpp + mpc(0,'0.05')], t1, Av, Ad, Bv, Bd)
res['xi_I3_plus'] = d1/c1
res['c1_I3_plus'] = c1
# --- I_6 fold at t = 1, REAL path
path = [mpf('0.14'), mpf('0.18'), mpf('0.22'), mpf('0.3'), mpf('0.5'), mpf('0.7'), mpf('0.9')]
c0,c1,d0,d1 = logcoeffs(mpc(1), mpc('0.9'), [mpc(p) for p in path], t1, Av, Ad, Bv, Bd)
res['xi_I6'] = d1/c1
res['c1_I6'] = c1

for k in ['xi_I3_minus','xi_I3_plus','xi_I6']:
    v = res[k]
    print()
    print(k)
    print("  Re =", nstr(mp.re(v), 165))
    print("  Im =", nstr(mp.im(v), 165))
for k in ['c1_I3_minus','c1_I3_plus','c1_I6']:
    v = res[k]
    print(k, " Re=", nstr(mp.re(v), 40), " Im=", nstr(mp.im(v), 40))

with open('out/xi_dps%d.txt' % DPS, 'w') as f:
    for k in ['xi_I3_minus','xi_I3_plus','xi_I6','c1_I3_minus','c1_I3_plus','c1_I6']:
        f.write(k+" "+nstr(mp.re(res[k]), DPS-20)+" "+nstr(mp.im(res[k]), DPS-20)+"\n")
