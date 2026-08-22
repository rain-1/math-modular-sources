"""Extend Delta(r) = BC(phi_r) - log|phi_r'(0)| for the Kodaira family
phi_r = lambda(r z) to r in [0.80, 0.96], needed for G(m) at m up to ~14.
Reuses lattice/cdt_noncongruence/bc_multivalent.delta_of_r."""
import sys, json
sys.path.insert(0, '/home/ubuntu/code/math-modular-sources/lattice/cdt_noncongruence')
from bc_multivalent import delta_of_r
out = {}
for i in range(17):
    r = round(0.80 + 0.01*i, 2)
    if r > 0.96: break
    d = delta_of_r(r, 192)
    out[str(r)] = d
    print("r=%.2f  Delta=%.6f" % (r, d), flush=True)
json.dump(out, open('/home/ubuntu/code/math-modular-sources/lattice/holonomy_lindep/delta_ext.json','w'), indent=0)
