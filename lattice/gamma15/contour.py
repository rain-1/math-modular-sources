# -*- coding: utf-8 -*-
"""Contour machinery for the Calegari-Dimitrov-Tang arithmetic-holonomy architecture.

h(z)  = -256 z prod_{n>=1} (1+z^n)^24 = -256 Delta(2 tau)/Delta(tau),  z = q = e^{2 pi i tau}
        the Gamma_0(2) Hauptmodul h = lambda + lambda/(lambda-1) = lambda^2/(lambda-1),
        an orbifold uniformisation of P^1 minus {0, oo} with a Z/2 point at y = 4.
        (CDT, L2chi.tex, eq. (defofh).)

Slit(z,r) : (D,0) -> (D minus (-1,-r], 0),  |Slit'(0,r)| = 4r/(1+r)^2   [CDT eq. (slit)]
lune(z,c) : (D,0) -> (L(c),0),              |lune'(0,c)| = (c^2-1)/(c^2+1)  [CDT (firstgobble)]

CDT's published map:
  G(z) = -R * lune( -e(th1) Slit( e(th2) Slit( e(th3) Slit( e(th4) Slit(z,r1), r2), r3), r4), c)
  GGG(z) = G(0.995 z),  phi = h o GGG.
"""
import numpy as np
import cmath, math
from fractions import Fraction as Fr

# ---------------------------------------------------------------- h and h'
def h(z, K=400):
    z = np.asarray(z, dtype=complex)
    p = np.ones_like(z)
    zn = z.copy()
    for n in range(1, K+1):
        p *= (1.0 + zn)**24
        zn = zn*z
        if np.all(np.abs(zn) < 1e-320):
            break
    return -256.0*z*p

def hp(z, K=400):
    """h'(z) by logarithmic differentiation: h'/h = 1/z + 24 sum_n n z^{n-1}/(1+z^n)."""
    z = np.asarray(z, dtype=complex)
    s = np.zeros_like(z)
    zn = np.ones_like(z); zz = z.copy()
    for n in range(1, K+1):
        s += 24.0*n*zn/(1.0+zz)
        zn = zn*z; zz = zz*z
        if np.all(np.abs(zn) < 1e-320):
            break
    return h(z, K)*(1.0/z + s)

def _bezout(x, y):
    old_r, r = x, y
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r != 0:
        q = old_r // r
        old_r, r = r, old_r - q*r
        old_s, s = s, old_s - q*s
        old_t, t = t, old_t - q*t
    assert old_r in (1, -1)
    if old_r == -1:
        old_s, old_t = -old_s, -old_t
    return old_s, old_t

def preimages_h(target, rmax=0.999, cdmax=200, K=400):
    """All z in D with |z| < rmax and h(z) = target, from the Gamma_0(2)-orbit of the
    principal preimage.  Preimages <-> bottom rows (2c,d), gcd(2c,d)=1, up to sign.
    For a negative real target the principal preimage is the unique z in (0,1)."""
    if abs(complex(target).imag) < 1e-14 and complex(target).real < 0:
        lo, hi = 1e-300, 1.0-1e-12
        for _ in range(400):
            mid = 0.5*(lo+hi)
            if complex(h(mid, K)).real > complex(target).real: lo = mid
            else: hi = mid
        z = 0.5*(lo+hi)
    else:
        z = complex(target)/(-256.0)
    for _ in range(80):
        z = z - (complex(h(z, K)) - target)/complex(hp(z, K))
    assert abs(complex(h(z, K)) - target) < 1e-13*max(1.0, abs(target)), (z, target)
    tau0 = complex(np.log(z)/(2j*np.pi))
    out, seen = [], set()
    for c in range(0, cdmax+1):
        for d in range(-cdmax, cdmax+1):
            if (c, d) == (0, 0): continue
            if math.gcd(2*c, abs(d)) != 1: continue
            key = (c, d) if (c > 0 or (c == 0 and d > 0)) else (-c, -d)
            if key in seen: continue
            seen.add(key)
            cc, dd = key
            a, b = _bezout(dd, 2*cc)
            t = (a*tau0 - b)/(2*cc*tau0 + dd)
            zz = cmath.exp(2j*math.pi*t)
            if abs(zz) < rmax:
                out.append((zz, (2*cc, dd), abs(zz)))
    out.sort(key=lambda t: t[2])
    return out

# ---------------------------------------------------------------- slits and lunes
def Slit(z, r):
    """(D,0) -> (D minus (-1,-r], 0).  CDT eq. (slit)."""
    z = np.asarray(z, dtype=complex)
    A = (1+r)**2 - 2*(1-6*r+r*r)*z + (1+r)**2*z*z
    S = np.sqrt(A)
    num = (r+1)**2 - 2*(r-1)**2*z + (r+1)**2*z*z + (1+r)*(-1+z)*S
    return np.where(np.abs(z) < 1e-8,
                    4*r/(1+r)**2*z + 8*r*(1-r)**2/(1+r)**4*z**2
                    + 4*(1-r)**2*r*(3-14*r+3*r*r)/(1+r)**6*z**3,
                    num/(8*r*np.where(z == 0, 1.0, z)))

def Slit_inv(w, r):
    """Inverse of Slit(.,r).  Derived from CDT's own chain
       w -> v = (1+w)^2/w -> u = i sqrt((v-A)/(v-B)) -> z = (u-i)/(u+i),
       A = -(1-r)^2/r, B = 4."""
    w = np.asarray(w, dtype=complex)
    A = -(1-r)**2/r
    v = (1+w)**2/w
    u = 1j*np.sqrt((v-A)/(v-4.0))
    return (u-1j)/(u+1j)

def lune(z, c):
    """CDT eq. (firstgobble): (D,0) -> (L(c),0), L(c) = D minus a horoball at z = -1
    reaching in to -(c-1)/(c+1)."""
    z = np.asarray(z, dtype=complex)
    B = (1+c*c)**2*(1+z)**2 - 16*c*c*z
    return (z*(1+c*c) - 1 - c*c + np.sqrt(B))/(2*(c*c-1))

def lune_inv(z, c):
    z = np.asarray(z, dtype=complex)
    return z*((c*c+1) + (c*c-1)*z)/((c*c-1) + (c*c+1)*z)

def slit_rad(r): return 4*r/(1+r)**2
def lune_rad(c): return (c*c-1)/(c*c+1)
def e(th):       return cmath.exp(2j*math.pi*th)

# ---------------------------------------------------------------- CDT's published map
CDT_R  = Fr(77, 100)
CDT_c  = Fr(75, 10)
CDT_r  = [Fr(91, 100), Fr(6188, 10000), Fr(55515, 100000), Fr(772, 1000)]
CDT_th = [Fr(7977, 100000), Fr(11543, 100000), Fr(3525, 100000), Fr(-783, 10000)]
CDT_shrink = Fr(995, 1000)

def make_G(R, c, rs, ths, shrink):
    R = float(R); c = float(c)
    rs = [float(t) for t in rs]; ths = [float(t) for t in ths]; sh = float(shrink)
    def G(z):
        z = np.asarray(z, dtype=complex)*sh
        u = Slit(z, rs[0])
        u = Slit(e(ths[3])*u, rs[1])
        u = Slit(e(ths[2])*u, rs[2])
        u = Slit(e(ths[1])*u, rs[3])
        u = lune(-e(ths[0])*u, c)
        return -R*u
    return G

def G_rad(R, c, rs, shrink):
    v = Fr(shrink)*Fr(R)*(Fr(c)**2-1)/(Fr(c)**2+1)
    for r in rs:
        v *= Fr(4)*Fr(r)/(1+Fr(r))**2
    return v
