"""Bost-Charles numerators on the EMN host.  BC(phi) = int int log|phi(z)-phi(w)|,
computed as int int log|(phi(z)-phi(w))/(z-w)| (the diagonal via log|phi'|)."""
import numpy as np, cmath, math, sys

def BC(phi, dphi, N=2048):
    zs = np.array([cmath.exp(2j*math.pi*j/N) for j in range(N)])
    fs = np.array([phi(z) for z in zs])
    fp = np.array([dphi(z) for z in zs])
    tot = 0.0
    for j in range(N):
        d = fs - fs[j]; dz = zs - zs[j]
        with np.errstate(divide='ignore', invalid='ignore'):
            q = np.abs(d/dz)
        q[j] = abs(fp[j])
        tot += np.sum(np.log(q))
    return tot/N**2

print("sanity: BC(rho z) must be log rho")
for rho in [0.5, 3.0, 161.081157]:
    print(f"   rho={rho:<12} BC={BC(lambda z,r=rho: r*z, lambda z,r=rho: r):.8f}  log rho={math.log(rho):.8f}")

print()
print("=== z-line host C\\{1}, base 0.  phi_R(z) = 1 - exp(-R z),  phi_R'(0)=R ===")
print("   (the host is NOT hyperbolic, so log|phi'(0)| is unbounded; the binding")
print("    quantity is BC(phi)/(log|phi'(0)| - tau), minimised over R.)")
TAU = 2.25   # module {1, H}: cols [(1,2)], e=[0,1]
print(f"{'R':>10} {'log R':>10} {'BC':>12} {'entry':>10} {'bound m <=':>12}")
best=(None,1e18)
for R in [3,5,7,9.5,10,12,15,20,30,50,100]:
    phi = lambda z,r=R: 1-cmath.exp(-r*z)
    dphi = lambda z,r=R: r*cmath.exp(-r*z)
    bc = BC(phi, dphi, N=1024)
    ent = math.log(R)-TAU
    bd = bc/ent if ent>0 else float('inf')
    if ent>0 and bd<best[1]: best=(R,bd)
    print(f"{R:>10} {math.log(R):>10.5f} {bc:>12.5f} {ent:>10.5f} {bd:>12.4f}")
print(f"   best bound over this family: m <= {best[1]:.4f} at R={best[0]}")
print(f"   the inventory has m = 2  =>  no contradiction (need m > bound).")

print()
print("=== same, with a richer inventory (tau smaller but m larger) ===")
for (label, m, TAU2) in [("{1,H}",2,2.25), ("{1,A,Li2,Li3,H}",5,2.4689), ("{1..10 lcm-free, H}",11,3.9438)]:
    best=(None,1e18)
    for R in [5,7,10,15,20,30,50,80,150,300]:
        ent = math.log(R)-TAU2
        if ent<=0: continue
        phi = lambda z,r=R: 1-cmath.exp(-r*z)
        dphi = lambda z,r=R: r*cmath.exp(-r*z)
        bc = BC(phi, dphi, N=768)
        bd = bc/ent
        if bd<best[1]: best=(R,bd)
    print(f"   {label:<24} m={m:<3} tau={TAU2:.4f}  best bound m <= {best[1]:.3f} at R={best[0]}   "
          f"{'CONTRADICTION' if m>best[1] else 'no contradiction'}")
