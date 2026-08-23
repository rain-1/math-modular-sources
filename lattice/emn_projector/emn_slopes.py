"""p-adic slope census for the projector basis, to n = N (default 1000)."""
from fractions import Fraction as Fr
import math, sys, json
from emn_cand import build, vp

N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
NS = [N//4, N//2, N]
rows = []
for a in range(1, 13):
    for kind in ('0', 'p3'):
        F = build({(a, kind): Fr(1)}, N)
        d = math.log(F[N].denominator)/N
        sl = {}
        for p in (2,3,5,7):
            v = vp(F[N], p)
            sl[p] = (v, v/N)
        name = f"h({a}S)" if kind=='0' else f"hpair({a}S,pi/3)"
        rows.append((name, d, sl))
        print(f"{name:24s} logden/n={d:.4f}  " +
              "  ".join(f"v_{p}={sl[p][0]:>6d}({sl[p][1]:+.3f})" for p in (2,3,5,7)))
json.dump([(r[0], r[1], {str(k): v for k,v in r[2].items()}) for r in rows],
          open('/home/ubuntu/code/math-modular-sources/lattice/emn_projector/slopes.json','w'), indent=1)
