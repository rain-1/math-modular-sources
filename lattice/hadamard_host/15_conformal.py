"""15_conformal.py -- the conformal data of the Hadamard host, in the
normalisation of CDT_FINDER.md sec.1 (|phi'(0)| and the Bost-Charles numerator).

Domain: Omega(a,b) = C \\ ( (-inf,-a] u [b,+inf) ), 0 < a,b <= inf.
Riemann map: phi = M^{-1} o K,  K(z) = (a/b)((1+z)/(1-z))^2  (Koebe onto C\\(-inf,0]),
             M(w) = -(w+a)/(w-b),  M^{-1}(zeta) = (b*zeta - a)/(1+zeta).
   phi'(0) = 4ab/(a+b)   [exact; = 4a for b = inf, the one-slit Koebe value]
For a univalent phi, BC(phi) = int int log|phi(z)-phi(w)| = log|phi'(0)| (Grunsky),
which is the convention of CDT_NONCONGRUENCE.md sec.1.4 (architecture (K)).
Here we (i) confirm phi'(0) numerically, (ii) confirm univalence on a grid,
(iii) evaluate the rearrangement numerator int_0^1 2t (log|phi|)^* dt of CDT (6.0.15)
     on circles |z| = r < 1 (|phi| is unbounded on |z| = 1).
"""
import math, cmath
import numpy as np

def phi_factory(a, b):
    if b == float('inf'):
        return lambda z: -a*((1+z)/(1-z))**2, 4.0*a
    p = a/b
    def f(z):
        zeta = p*((1+z)/(1-z))**2
        return (b*zeta - a)/(1.0 + zeta)
    return f, 4.0*a*b/(a+b)

def check(a, b, label):
    f, rho = phi_factory(a, b)
    h = 1e-6
    d = (f(h) - f(-h))/(2*h)
    print(f"{label:<34s} a={a:<12.6g} b={b:<12.6g}  rho(formula)={rho:.10g}  rho(numeric)={abs(d):.10g}"
          f"  log rho = {math.log(rho):+.6f}")
    # univalence spot check: injectivity on a polar grid
    zs = [r*cmath.exp(2j*math.pi*t) for r in (0.3,0.6,0.85,0.95) for t in np.linspace(0,1,241)[:-1]]
    ws = [f(z) for z in zs]
    mind = min(abs(ws[i]-ws[j]) for i in range(0,len(ws),7) for j in range(i+1,len(ws),7))
    print(f"{'':34s} min |phi(z)-phi(w)| over a 4x240 grid (sampled) = {mind:.3e}  (0 would mean non-injective)")
    # rearrangement numerator on |z| = r
    for r in (0.9, 0.99, 0.999):
        N = 200000
        t = (np.arange(N)+0.5)/N
        z = r*np.exp(2j*np.pi*t)
        g = np.log(np.abs(np.array([f(zz) for zz in z])))
        gs = np.sort(g)[::-1]              # decreasing rearrangement
        tt = (np.arange(N)+0.5)/N
        val = np.sum(2*tt*gs)/N
        print(f"{'':34s} r={r:<6} rearrangement int_0^1 2t (log|phi_r|)^* dt = {val:+.6f}"
              f"   (log|phi_r'(0)| = {math.log(rho*r):+.6f})")

s1t1 =  3.4875690311312326e-07
s1t2 =  1.0971741115134826
s2t1 = -0.6488615207279722
s2t2 = -2041290.2401220216
print("Hadamard host singular points:", s1t1, s1t2, s2t1, s2t2)
print()
check(abs(s2t1), s1t2,          "W and its theta-orbit")
check(abs(s2t2), s1t2,          "COND_Z = A_Z(*)(bB_N-aA_N)")
check(abs(s2t2), float('inf'),  "DBL (doubly conditional)")
check(abs(s2t1), s1t1,          "full module (incl. A_Z(*)A_N)")
