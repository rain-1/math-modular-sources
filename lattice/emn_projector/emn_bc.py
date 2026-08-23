"""Bost-Charles numerators for the EMN host P^1-{0,1/2,oo} in Z=z/2.

phi(q) = lambda(q)/2  maps (D,0) -> C-{1/2}, phi'(0) = 8, and phi^{-1}(0)={0}.
Concentric contours |q| < r.  No fold preimage constraint (the hypothesis has
deleted Z=1; lambda(q)=2 is attained but Z=1 is a regular point of the
conditional function).
"""
import sys, math
sys.path.insert(0, '/home/ubuntu/code/math-modular-sources/lattice/catalan_mu4')
from mu4_bc import BC, lam

def phi(q):
    return lam(q)/2.0

print("sanity: BC(rho z) = log rho ->", BC(lambda z: 3.0*z, 1.0, 2048), math.log(3.0))
print("\n r      log|phi'(0)|      BC        BC - log|phi'(0)|")
for r in (0.2, 0.4, 0.6, 0.75, 0.85, 0.92):
    lp = math.log(8.0*r)
    b = BC(phi, r, 2048)
    print(f" {r:<6} {lp:+.5f}      {b:+.5f}      {b-lp:+.5f}")
