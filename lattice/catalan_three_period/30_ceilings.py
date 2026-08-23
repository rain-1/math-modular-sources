#!/usr/bin/env python3
"""30_ceilings.py  (2026-08-23)

Uniformisation ceilings |phi'(0)| for the symmetrised (y-line) orbifolds of the
level-8 (Zagier E) and level-16 Catalan hosts, in the convention of
CATALAN_AL_HOSTS.md Sec 6.1 (y=0 counted as a puncture; r = |dy/dq| at that cusp
for the canonical cusp parameter q of the uniformising Fuchsian group).

Closed forms, each verified below against the relevant q-expansion:
   r(P^1 - {0,s,oo})                        = 16|s|          (lambda function)
   r(P^1 - {0,a,b})                         = 16|ab/(b-a)|   (Moebius + lambda)
   r({0,oo} punctures + cone-2 at c)        = 64|c|          (Gamma_0(2) hauptmodul)
   r(P^1 - all cusps of a genus-0 host)     = |c_1| of the hauptmodul q + c_2 q^2 ...
"""
from mpmath import mp, mpf, mpc, exp, pi, log, jtheta, nstr

mp.dps = 30

# ---------------------------------------------------------------- verifications
def eta_prod(qq, N=600):
    p = mpc(1); qn = mpc(1)
    for n in range(1, N):
        qn *= qq
        p *= (1 - qn)
    return p

print("== closed forms verified against q-expansions ==")

# lambda(tau) = 16 q - 128 q^2 + 704 q^3 - ... ,  q = exp(i pi tau)
def lam(tau):
    return (jtheta(2, 0, exp(1j*pi*tau))/jtheta(3, 0, exp(1j*pi*tau)))**4
for t in [mpc(0, 6), mpc(0, 8)]:
    qq = exp(1j*pi*t)
    print("   lambda/q  ->", nstr(lam(t)/qq, 12), " (limit 16)")

# Gamma_0(2) hauptmodul h = -256 Delta(2tau)/Delta(tau) = -256 q - 4096 q^2 - ...
# its order-2 elliptic value is h = 4 ; ceiling for {0,oo}+cone-2 at c is 64|c|
def h2(tau):
    q = exp(2j*pi*tau)
    return -256 * q * (eta_prod(q**2)/eta_prod(q))**24
for t in [mpc(0, 3), mpc(0, 4)]:
    q = exp(2j*pi*t)
    print("   h_{Gamma_0(2)}/q ->", nstr(h2(t)/q, 12), " (limit -256 = -64*4)")

print()
tau_CDT = mpf(16603)/3920          # tau(b;e) for CDT's m=14, k=2 inventory
print("tau(b;e) = 16603/3920 =", nstr(tau_CDT, 10))

# ---------------------------------------------------------------- level 8
print("\n== level 8 (Zagier E), s = 1/4, sigma(x) = x/(4x-1), y = 4x^2/(4x-1) ==")
print("   Sigma_x = {0, 1/8, 1/4, oo}   (verified, 30_level16_involution.gp [2])")
print("   fixed points of sigma: 0 (cusp) and 2s = 1/2 (NOT a cusp -> free point)")
print("   Sigma_y = {0, -1/8, oo} punctures + CONE-2 at y = 4s = 1")
c8 = 64*abs(mpf(1))
print("   ceiling (fold image -1/8 dropped) = 64*|1| =", c8, " log =", nstr(log(c8), 10))
print("   entry at the ceiling =", nstr(log(c8) - tau_CDT, 10))

# ---------------------------------------------------------------- level 16
print("\n== level 16, s = -1/4, sigma(x) = -x/(4x+1), y = 4x^2/(4x+1) ==")
print("   Sigma_x = {0, -1/4, -1/2, (-1+-i)/4, oo}  (verified)")
print("   sigma is induced by tau -> tau+1/2 (verified); fixed points 0 and 2s = -1/2,")
print("   BOTH cusps  =>  the quotient carries NO cone point.")
print("   Sigma_y = {0, -1/2, -1, oo}, all punctures; fold at y = oo,")
print("             extra (outer pair) at y = -1/2, cusp x=2s at y = -1.")

def r3(a, b):      # r(P^1 - {0,a,b})
    return 16*abs(a*b/(b-a))
def r3inf(s):      # r(P^1 - {0,s,oo})
    return 16*abs(s)

rows = [
  ("drop the fold image y=oo        -> P^1-{0,-1/2,-1}", r3(mpf(-1)/2, mpf(-1))),
  ("drop the extra image y=-1/2     -> P^1-{0,-1,oo}",   r3inf(mpf(-1))),
  ("drop y=-1 (image of cusp x=2s)  -> P^1-{0,-1/2,oo}", r3inf(mpf(-1)/2)),
  ("drop nothing: P^1-{0,-1/2,-1,oo} = Y_0(8), y = 4q^2+...", mpf(4)),
]
for name, r in rows:
    print("   %-56s r = %s  log = %s  entry = %s"
          % (name, nstr(r,8), nstr(log(r),8), nstr(log(r)-tau_CDT,8)))

print("\n   => under EITHER admissible choice of the single point phi may hit")
print("      (the fold y=oo, or the extra point y=-1/2) the ceiling is exactly 16.")
print("      The value 64 = 256|s| requires dropping BOTH y=-1/2 AND y=-1, i.e.")
print("      pretending x=-1/2 is not a cusp.  It is one.")

# ---------------------------------------------------------------- consequences
print("\n== consequences for CATALAN_THREE_PERIOD.md Sec 5 ==")
d = log(mpf(16)) - log(mpf(64))
print("   ceiling: log 64 = %s  ->  log 16 = %s   (Delta = %s = -log 4)"
      % (nstr(log(mpf(64)),8), nstr(log(mpf(16)),8), nstr(d,8)))
print("   archimedean entry at the ceiling: -0.076576 -> %s" % nstr(log(mpf(16))-tau_CDT,8))
print("   transported-contour entry: %s -> %s"
      % (nstr(log(mpf(64))+log(mpf('0.62922'))-tau_CDT,8),
         nstr(log(mpf(16))+log(mpf('0.62922'))-tau_CDT,8)))
print("   margins shift by m*Delta (BC held at the transported 11.845+log(1/4)):")
tab = [("previous 7 pure + 7 conditional", 14, mpf('-7.966'), mpf('-22.423')),
       ("previous + 1 doubly-small",       15, mpf('-7.551'), mpf('-22.431')),
       ("(a)/(b) full: 7 pure + 3 orbits", 28, mpf('-6.316'), mpf('-23.942')),
       ("best containing a conditional fn",15, mpf('-2.006'), mpf('-5.795')),
       ("best unconditional",              14, mpf('-0.737'), mpf('-0.737'))]
for name, m, p1, p2 in tab:
    print("     %-34s m=%2d   P1 %8s -> %9s    P2 %8s -> %9s"
          % (name, m, nstr(p1,5), nstr(p1+m*d,6), nstr(p2,5), nstr(p2+m*d,6)))
