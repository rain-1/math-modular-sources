"""Final table: re-evaluate every stored template at high resolution and emit
the per-target results, kappa, and the numerical contour description."""
import os, sys as _sys
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in _sys.path: _sys.path.insert(0, HERE)
import json, glob, math, os, warnings, sys
import numpy as np
warnings.filterwarnings('ignore')
import haupt, outer, certify
from targets import TARGETS

BEST = {}   # idx -> (tag, coeff vector c)

for f in sorted(glob.glob('free_*.json')):
    d = json.load(open(f))
    idx = [i for i, T in enumerate(TARGETS) if T['key'] == d['key']]
    if not idx:
        continue
    idx = idx[0]
    import freeopt
    u = freeopt.u_of(np.array(d['a']), 8192)
    c = certify.coeffs_from_u(u, 60)
    mx = certify.fine_max_u(c)
    if mx > 0:
        c[0] -= mx + 1e-12
    r = certify.evaluate_c(TARGETS[idx], c, 8192, True)
    if r is None:
        continue
    cost = r['RE']['cost']
    if idx not in BEST or cost < BEST[idx][2]:
        BEST[idx] = (f, c, cost)

out = {}
for idx, (tag, c, _) in sorted(BEST.items()):
    T = TARGETS[idx]
    r = certify.evaluate_c(T, c, 16384, True)
    m, L = T['m'], T['L']
    rec = dict(key=T['key'], idx=idx, d=T['d'], m=m, L=L, tau=T['tau'], budget=T['budget'],
               src=tag, c=list(c), maxq=r['maxq'], gmax=r['gmax'], gmin=r['gmin'])
    for w in ('RE', 'BC'):
        if w not in r:
            continue
        rec[w] = r[w]
        mg = r[w]['margin']
        g = 1.0 / m
        k = None
        if mg > 0:
            s = 1 - mg / (m * L)
            if 0 <= s < 1:
                den = 1 - math.sqrt(s)
                if den > 0:
                    k = (1 - g) / den
        rec[w]['kappa'] = k
    out[T['key']] = rec
    print("%-26s m=%d budget=%8.5f | RE bound=%9.6f margin=%+9.6f kappa=%s | BC bound=%9.6f margin=%+9.6f kappa=%s | logdr=%+.6f maxq=%.4f gmax=%.3f  [%s]"
          % (T['key'], m, T['budget'], rec['RE']['bound'], rec['RE']['margin'],
             ("%.1f" % rec['RE']['kappa']) if rec['RE']['kappa'] else "  -  ",
             rec['BC']['bound'], rec['BC']['margin'],
             ("%.1f" % rec['BC']['kappa']) if rec['BC']['kappa'] else "  -  ",
             rec['RE']['logdr'], rec['maxq'], rec['gmax'], tag), flush=True)

json.dump(out, open('FINAL.json', 'w'), indent=1)

# numerical contour description
print("\n== contour profiles  rho(theta)=|psi(e^{i theta})| and log|x| ==", flush=True)
for k, rec in out.items():
    c = np.array(rec['c'])
    M = 24
    t = np.arange(M) / M
    u = np.full(M, c[0])
    for j in range(1, len(c)):
        u = u + c[j] * np.cos(2 * np.pi * j * t)
    idx = rec['idx']
    H = TARGETS[idx]['H']
    uu = certify.u_from_c(c, 4096)
    q4, _ = outer.contour(uu)
    g4 = haupt.logabs_x(H, q4)
    # sample every 4096/M
    st = 4096 // M
    print("-- %s" % k)
    print("   theta/2pi : " + " ".join("%6.3f" % v for v in t))
    print("   rho       : " + " ".join("%6.4f" % v for v in np.exp(u)))
    print("   log|x|    : " + " ".join("%6.2f" % g4[i * st] for i in range(M)), flush=True)
