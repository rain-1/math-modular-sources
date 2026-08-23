#!/usr/bin/env python3
"""EMN note claims (b)-pure-part, (c), (d): 50-digit numerics.
Renamed to avoid builtins: Phi->phi_ser, H->hfun, S->uni."""
from mpmath import mp, mpf, catalan, asin, sqrt, cos, sin, log, tan, pi, quad, mpmathify, nstr, diff, atan
from fractions import Fraction

mp.dps = 60
GG = catalan

# ---- exact a_n as Fractions
NMAX = 4000
def amoms(N):
    out = []
    for n in range(N+1):
        s = Fraction(0)
        for k in range(n+1):
            from math import comb
            s += Fraction(comb(n,k), 2*k+1)
        out.append(s / (Fraction(2)**(n+1) * (n+1)))
    return out

# faster: use the two-term recurrence (m+1)(2m+1) a_m = m^2 a_{m-1} + 1/2
def amoms_rec(N):
    out=[Fraction(1,2)]
    for m in range(1,N+1):
        out.append((Fraction(m*m)*out[m-1] + Fraction(1,2)) / Fraction((m+1)*(2*m+1)))
    return out

A = amoms_rec(NMAX)
Achk = amoms(40)
assert A[:41]==Achk, "recurrence vs closed form mismatch"
print("a_n via two-term recurrence agrees with closed form to n=40: OK")

Amp = [mp.mpf(x.numerator)/mp.mpf(x.denominator) for x in A]

def hfun(z, N=NMAX):
    """H(z) = z*Phi(z) = sum a_n z^{n+1}, |z|<1 (slow at z=1)."""
    z = mpmathify(z); tot = mp.mpf(0); p = z
    for n in range(N+1):
        tot += Amp[n]*p; p *= z
    return tot

def hfourier(uni, M=200000):
    """script-H(S) = sum (-1)^m/(2m+1)^2 (1-cos((2m+1)S))"""
    uni = mpmathify(uni)
    return mp.nsum(lambda m: (-1)**int(m)/mp.mpf(2*int(m)+1)**2*(1-mp.cos((2*int(m)+1)*uni)),
                   [0, mp.inf], method='a')

print()
print("=== (b/3) uniformisation: script-H(S) == H(1-cos S) ===")
for uni in [mp.mpf('0.3'), mp.mpf('0.7'), mp.mpf('1.0'), mp.mpf('1.4')]:
    z = 1-cos(uni)
    lhs = hfourier(uni); rhs = hfun(z)
    print(f"  S={nstr(uni,6)}  z={nstr(z,10)}  Fourier={nstr(lhs,30)}  series={nstr(rhs,30)}  diff={nstr(lhs-rhs,5)}")

print()
print("=== (3) script-H'(S) = (1/2) log tan(pi/4 + S/2) ===")
for uni in [mp.mpf('0.3'), mp.mpf('0.9'), mp.mpf('1.3')]:
    d = diff(lambda t: hfourier(t), uni)
    tgt = log(tan(pi/4+uni/2))/2
    print(f"  S={nstr(uni,6)}  H'={nstr(d,25)}  (1/2)log tan={nstr(tgt,25)}  diff={nstr(d-tgt,5)}")

print()
print("=== (c) H(z)+H(2-z) = 2G   [equivalently script-H(S)+script-H(pi-S)=2G] ===")
for uni in [mp.mpf('0.4'), mp.mpf('1.1'), mp.mpf('1.5')]:
    lhs = hfourier(uni)+hfourier(pi-uni)
    print(f"  S={nstr(uni,6)} z={nstr(1-cos(uni),8)}: H(z)+H(2-z)={nstr(lhs,40)}   2G={nstr(2*GG,40)}  diff={nstr(lhs-2*GG,5)}")

print()
print("=== (c) H(3/2) = (5/3) G   [S = 2pi/3] ===")
h32 = hfourier(2*pi/3)
print("  H(3/2)   =", nstr(h32, 50))
print("  (5/3)G   =", nstr(mp.mpf(5)/3*GG, 50))
print("  diff     =", nstr(h32-mp.mpf(5)/3*GG, 5))

print()
print("=== (c) G = (3/2) Phi(1/2),  i.e. H(1/2)=G/3 ===")
h12 = hfun(mp.mpf(1)/2)
print("  H(1/2)         =", nstr(h12, 50))
print("  G/3            =", nstr(GG/3, 50))
print("  diff           =", nstr(h12-GG/3, 5))
print("  (3/2)Phi(1/2)  =", nstr(3*h12, 50), "  [= 3*H(1/2) since Phi(1/2)=2H(1/2)]")
print("  G              =", nstr(GG, 50))

print()
print("=== (c) R_3(z) = z(3-2z)^2  is  1-cos(3S) with z=1-cos S ===")
for uni in [mp.mpf('0.37'), mp.mpf('1.21')]:
    z = 1-cos(uni)
    print(f"  S={nstr(uni,6)}: 1-cos3S={nstr(1-cos(3*uni),30)}   z(3-2z)^2={nstr(z*(3-2*z)**2,30)}")

print()
print("=== (c) cubic distribution: H(R3(z)) + 3*sum_j H(z_j) = 10G ===")
for uni in [mp.mpf('0.25'), mp.mpf('0.6'), mp.mpf('1.05')]:
    tot = hfourier(3*uni) + 3*sum(hfourier(uni+2*pi*j/3) for j in range(3))
    zs = [1-cos(uni+2*pi*j/3) for j in range(3)]
    print(f"  S={nstr(uni,5)}: z_j={[nstr(t,8) for t in zs]}")
    print(f"     total={nstr(tot,40)}   10G={nstr(10*GG,40)}   diff={nstr(tot-10*GG,5)}")

print()
print("=== (d) PV Phi(3/2) = (10/9) G   -- direct principal-value double integral ===")
# Phi(z) = -(1/(2z)) int_0^{pi/2} log|1 - z/(cos t + sin t)^2| dt   (radial PV)
def phi_pv(zz):
    zz = mpmathify(zz)
    f = lambda t: log(abs(1 - zz/(cos(t)+sin(t))**2))
    # singularities of the log at (cos t+sin t)^2 = zz
    pts = [mp.mpf(0)]
    if 1 < zz < 2:
        a = mp.acos(sqrt(zz)/sqrt(2))
        pts += [pi/4-a, pi/4+a]
    pts.append(pi/2)
    pts = sorted(pts)
    return -(quad(f, pts))/(2*zz)

for zz, tgt, name in [(mp.mpf(1)/2, mp.mpf(2)/3*GG, "(2/3)G"),
                      (mp.mpf(1),   GG,             "G"),
                      (mp.mpf(3)/2, mp.mpf(10)/9*GG,"(10/9)G")]:
    v = phi_pv(zz)
    print(f"  PV Phi({nstr(zz,4)}) = {nstr(v,50)}")
    print(f"        target {name} = {nstr(tgt,50)}   diff={nstr(v-tgt,5)}")

print()
print("=== (6) quadratic-unit form: G = (3/2) Ti_2(2-sqrt3) + (pi/8) log(2+sqrt3) ===")
def ti2(x):
    x = mpmathify(x)
    return mp.nsum(lambda m: (-1)**int(m)*x**(2*int(m)+1)/mp.mpf(2*int(m)+1)**2, [0, mp.inf])
eps = 2-sqrt(3)
rhs = mp.mpf(3)/2*ti2(eps) + pi/8*log(2+sqrt(3))
print("  RHS =", nstr(rhs, 50))
print("  G   =", nstr(GG, 50))
print("  diff=", nstr(rhs-GG, 5))
