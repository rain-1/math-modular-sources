"""Task 2(2): the EMN host geometry and its conformal ceilings.
Renamed: no Phi/S/H/M."""
import mpmath as mp
mp.mp.dps = 30
from mpmath import mpf, mpc, pi, log, sqrt, exp, jtheta, nstr

# ---- lambda modular function and its derivative -------------------------
def lam_tau(tau):
    q = mp.exp(1j*pi*tau)
    return (jtheta(2,0,q)/jtheta(3,0,q))**4

def dlam_tau(tau):
    """lambda'(tau) = i pi lambda (1-lambda) theta_3(0|tau)^4."""
    q = mp.exp(1j*pi*tau)
    L = (jtheta(2,0,q)/jtheta(3,0,q))**4
    return 1j*pi*L*(1-L)*jtheta(3,0,q)**4

# sanity: numerical derivative
h = mpf(10)**-12
for t0 in [mpc(0,1), mpc(1,1)]:
    nd = (lam_tau(t0+h)-lam_tau(t0-h))/(2*h)
    print(f"  lambda({t0}) = {nstr(lam_tau(t0),12)}   lambda' formula {nstr(dlam_tau(t0),12)}  numeric {nstr(nd,12)}")

def confrad(tau0, scale=1, shift=0):
    """conformal radius at the point lambda(tau0) of C minus {0,1},
       pulled back by w -> scale*w + shift.  D -> H by tau = tau0_re + i(1+z)/(1-z)
       requires Im tau0 = 1; we use tau = tau0 + (i(1+z)/(1-z) - i) which has dtau/dz|_0 = 2i."""
    return abs(dlam_tau(tau0)*2j)*abs(scale)

print()
print("=== conformal radii (= sup |phi'(0)| over holomorphic phi: D -> X, phi(0)=basepoint) ===")
# C \ {0,1} at 1/2  (tau = i)
r_half = confrad(mpc(0,1))
print(f"  C\\{{0,1}} at 1/2      : {nstr(r_half,16)}")
# C \ {-1,1} at 0 : w = (zeta+1)/2 maps zeta -> w, {-1,1}->{0,1}, 0->1/2; dzeta/dw = 2
print(f"  C\\{{-1,1}} at 0       : {nstr(2*r_half,16)}   (log = {nstr(log(2*r_half),12)})")
print(f"     [CATALAN_MU4 sec.4 quotes 4.376879230453 for C\\{{+-i}} at 0 -- same by rotation]")
# C \ {-1/sqrt2, 1/sqrt2} at 0
print(f"  C\\{{+-1/sqrt2}} at 0  : {nstr(2*r_half/sqrt(2),16)}   (log = {nstr(log(2*r_half/sqrt(2)),12)})")
# C \ {1,2} at 0 : w = z-1, base -1 = lambda(1+i)
r_m1 = confrad(mpc(1,1))
print(f"  C\\{{0,1}} at -1       : {nstr(r_m1,16)}")
print(f"  C\\{{1,2}} at 0        : {nstr(r_m1,16)}   (log = {nstr(log(r_m1),12)})")
print()
print("=== Koebe caps (univalent phi only): |phi'(0)| <= 4*dist(0, omitted set) ===")
print(f"  z-plane, omit {{1}}          : 4        (log = {nstr(log(4),12)})")
print(f"  z-plane, omit {{1,2}}        : 4        (log = {nstr(log(4),12)})")
print(f"  zeta-plane, omit {{+-1/sqrt2}}: {nstr(4/sqrt(2),12)}  (log = {nstr(log(4/sqrt(2)),12)})")
print(f"  S-plane, omit (pi/2 + pi Z) : {nstr(4*pi/2,12)}  (log = {nstr(log(2*pi),12)})")
print()
print("=== the S-plane: sin: C\\(pi/2+piZ) -> C\\{+-1} is an unbranched covering ===")
print(f"  so the S-plane conformal radius at S=0 equals that of C\\{{+-1}} at 0 divided by |sin'(0)|=1:")
print(f"     {nstr(2*r_half,16)}   (log = {nstr(log(2*r_half),12)})")
print()
print("=== comparison table (hard ceilings, log scale) ===")
rows = [("CDT's own host P^1-{0,1,oo}, lambda template", log(mpf(16))),
        ("CDT symmetrised Y_0(2)", log(mpf(256))),
        ("Catalan level 8 unsymmetrised (16/lambda_2, lambda_2=4)", log(mpf(4))),
        ("Catalan level 8 symmetrised", log(mpf(64))),
        ("mu_4 host A (x-line), ceiling 16^{1/2}", log(mpf(4))),
        ("EMN z-line, omit {1} only, univalent Koebe", log(mpf(4))),
        ("EMN zeta-line (zeta^2=z/2), omit {+-1/sqrt2}, universal cover", log(2*r_half/sqrt(2))),
        ("EMN S-plane / C\\{+-1} at 0, universal cover", log(2*r_half)),
        ("EMN z-line, omit {1,2}, universal cover", log(r_m1))]
for n_,v in rows: print(f"  {n_:<58} {nstr(v,10)}")
