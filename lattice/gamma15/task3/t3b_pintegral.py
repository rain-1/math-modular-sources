"""Guard: every coefficient of the 14 functions (both hosts) is p-integral at the two
primes used in T3, so the reductions mod p are legitimate."""
import sys, io, contextlib
sys.argv = ['x']
sys.path.insert(0, '/home/ubuntu/code/math-modular-sources/lattice/gamma15/task3')
import importlib.util
spec = importlib.util.spec_from_file_location('t3', 't3_indep.py')
m = importlib.util.module_from_spec(spec)
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    spec.loader.exec_module(m)
from math import gcd
ps = [m.P_SPLIT, m.P_INERT]
print("primes:", ps, " both prime:", [__import__('sympy').isprime(p) for p in ps])
for tag, (H, s) in (('CDT', m.host_CDT()), ('G15', m.host_G15())):
    G = m.symplus(H, s)
    for var in ('cdtdef', 'indep_check2'):
        bad = 0; mx = 0
        for nm, f in m.fourteen(G, var):
            for cf in f:
                for q in cf:
                    d = q.denominator
                    for p in ps:
                        if d % p == 0: bad += 1
                    mx = max(mx, d)
        print(f"  {tag}/{var}: coefficients with p | denominator: {bad};  largest denominator has {len(str(mx))} digits")
