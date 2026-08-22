"""Rigorous archimedean ceiling for the level-16 y-picture, IF phi must avoid the
extra singularity y=-1/2 (log branch of the conditional / doubly-small functions).

Ceiling = max{|phi'(0)| : phi: D -> P^1 \ {-1/2,-1,infty}, phi(0)=0}
        = |b-a| * C_01(p),  a=-1/2, b=-1, p = -a/(b-a) = -1,
C_01(p) = 2 Im(tau_p) |lambda'(tau_p)|,  lambda(tau_p)=p.
lambda(1+i) = -1 ;  lambda'(tau) = i pi lambda(1-lambda) theta_3(0|tau)^4.
"""
import mpmath as mp
mp.mp.dps = 40
def lam(tau):
    q = mp.exp(1j*mp.pi*tau)
    t2 = 2*q**mp.mpf(0.25)*mp.nsum(lambda n: q**(n*(n+1)), [0, mp.inf])
    t3 = mp.jtheta(3,0,q); t2 = mp.jtheta(2,0,q)
    return (t2/t3)**4
def dlam(tau):
    q = mp.exp(1j*mp.pi*tau); l = lam(tau)
    return 1j*mp.pi*l*(1-l)*mp.jtheta(3,0,q)**4
for tau in [mp.mpc(0,1), mp.mpc(1,1)]:
    print("tau =", tau, " lambda =", mp.chop(lam(tau)), " |lambda'| =", abs(dlam(tau)))
tau0 = mp.mpc(1,1)
C01 = 2*mp.im(tau0)*abs(dlam(tau0))
print("C_01(-1) =", C01, "  (= 2 pi^2/Gamma(3/4)^4 =", 2*mp.pi**2/mp.gamma(0.75)**4, ")")
Mp = 2.0
C = C01/Mp
print("ceiling |phi'(0)| <=", C, "   log =", mp.log(C))
print("entry (tau=4.235459) =", mp.log(C)-4.235459)
# level-8 comparison: extra point y(t_1)=-1/8, branch 4s=1
for nm,a,b in [("CDT level 6", -1/72.0, 4.0), ("level 8", -1/8.0, 1.0), ("level 16", -0.5, -1.0)]:
    print(f"  {nm:12s} extra/branch = {abs(a)/abs(b):.6f}")
