"""14: exact basis of the joint solution space, extension to a big grid.
Reconstruction is retried with more primes until the basis verifies EXACTLY."""
import sys, json, time, pickle
from fractions import Fraction as F
from libx import *
import lib2v

name = sys.argv[1]
NBIG = int(sys.argv[2]) if len(sys.argv) > 2 else 32
N0 = int(sys.argv[3]) if len(sys.argv) > 3 else 10
recs = json.load(open("recs.json"))
r = recs[name]
reca = Rec(r["a"]["I"], r["a"]["mons"], r["a"]["coefs"], 0)
recb = Rec(r["b"]["I"], r["b"]["mons"], r["b"]["coefs"], 1)
cfun = getattr(lib2v, name)

t0 = time.time()
ok = False
for npr in (6, 8, 10, 12):
    basis, cells, free, err = grid_nullspace_exact(reca, recb, N0, N0, nprimes=npr)
    if basis is None:
        print("  nprimes=%d -> %s" % (npr, err)); continue
    Ts = [{c: v[j] for j, c in enumerate(cells)} for v in basis]
    viol = [len(check_solution(reca, recb, T, N0, N0)) for T in Ts]
    print("  nprimes=%d : dim=%d violations=%s" % (npr, len(basis), viol))
    sys.stdout.flush()
    if all(v == 0 for v in viol):
        ok = True
        break
assert ok, "could not reconstruct exactly"
freecells = [cells[j] for j in free]
print("exact nullspace at N=%d : dim=%d  free cells=%s  (%.1fs)"
      % (N0, len(basis), freecells, time.time()-t0))
for k, T in enumerate(Ts):
    print("  basis %d free values %s ; t(1,1)=%s t(2,2)=%s"
          % (k, [T[fc] for fc in freecells], T[(1, 1)], T[(2, 2)]))
sys.stdout.flush()

for k in range(len(Ts)):
    t1 = time.time()
    Ts[k], miss = propagate(reca, recb, Ts[k], NBIG, NBIG)
    bad = check_solution(reca, recb, Ts[k], NBIG, NBIG) if not miss else ["skip"]
    print("  basis %d -> N=%d : missing=%d %s ; violations=%s   (%.1fs)"
          % (k, NBIG, len(miss), miss[:6], len(bad), time.time()-t1))
    sys.stdout.flush()

cT = {(a, b): F(cfun(a, b)) for a in range(NBIG+1) for b in range(NBIG+1)}
coef = [cT[fc] for fc in freecells]
resid = max(abs(cT[(a, b)] - sum(co*T[(a, b)] for co, T in zip(coef, Ts)))
            for a in range(NBIG+1) for b in range(NBIG+1))
print("ROW c in the span:", resid == 0, " coefficients:", coef)
pickle.dump(dict(name=name, NBIG=NBIG, freecells=freecells, Ts=Ts, coef=coef),
            open("sol_%s.pkl" % name, "wb"))
print("saved sol_%s.pkl  total %.1fs" % (name, time.time()-t0))
