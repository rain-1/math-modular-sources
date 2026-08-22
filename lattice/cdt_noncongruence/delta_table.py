"""Delta(r) = BC(phi_r) - log|phi_r'(0)| for phi_r = lambda(r z), on a fine grid,
with a convergence check in the quadrature size N."""
import json, sys, math
from bc_multivalent import delta_of_r

rs = [round(0.20 + 0.01*i, 2) for i in range(61)]      # 0.20 .. 0.80
out = {}
for r in rs:
    d = delta_of_r(r, 256)
    out[str(r)] = d
    print(f"r={r:.2f}  Delta={d:.6f}", flush=True)
json.dump(out, open('delta_table.json','w'), indent=0)
print("\nconvergence check (N = 128, 192, 256, 384):")
for r in (0.30, 0.45, 0.60, 0.75):
    vals = [delta_of_r(r, N) for N in (128,192,256,384)]
    print("  r=%.2f  " % r + "  ".join("%.6f" % v for v in vals))
