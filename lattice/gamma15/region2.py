# -*- coding: utf-8 -*-
"""Contour region  Phi = (z -> R z) o T_1 o ... o T_k,  T_j(z) = e^{i beta_j} f_j(z),
f_j = Slit(.,r) or lune(.,c).  T removes, from the CURRENT region, the image of
  a slit: the ray at angle beta+pi from radius r out to 1;
  a lune: the horoball tangent to the circle at e^{i(beta+pi)}, reaching in to (c-1)/(c+1).
Because Omega_k is contained in Omega_j for k > j and the target point p_j lies IN the
segment removed at step j, exclusion of every p_j holds BY CONSTRUCTION.
Conformal radius = R * prod (4r/(1+r)^2) * prod ((c^2-1)/(c^2+1)) -- exact."""
import numpy as np, math, cmath
from fractions import Fraction as Fr
import contour as C

class Reg2:
    def __init__(self, R, ops=None, rho=1.0):
        self.R = float(R); self.ops = list(ops or []); self.rho = float(rho)
    def _fwd_raw(self, z):
        z = np.asarray(z, dtype=complex)
        for (kind, be, par) in reversed(self.ops):
            z = np.exp(1j*be)*(C.Slit(z, par) if kind == 's' else C.lune(z, par))
        return self.R*z
    def fwd(self, z):  return self._fwd_raw(self.rho*np.asarray(z, dtype=complex))
    def _inv_raw(self, w):
        z = np.asarray(w, dtype=complex)/self.R
        for (kind, be, par) in self.ops:
            u = np.exp(-1j*be)*z
            z = C.Slit_inv(u, par) if kind == 's' else C.lune_inv(u, par)
        return z
    def inv(self, w):  return self._inv_raw(w)/self.rho
    def rad(self):
        v = self.R*self.rho
        for (kind, be, par) in self.ops:
            v *= (4*par/(1+par)**2) if kind == 's' else ((par*par-1)/(par*par+1))
        return v
    def rad_exact(self):
        v = Fr(self.R).limit_denominator(10**12)*Fr(self.rho).limit_denominator(10**12)
        for (kind, be, par) in self.ops:
            p = Fr(par).limit_denominator(10**12)
            v *= (Fr(4)*p/(1+p)**2) if kind == 's' else ((p*p-1)/(p*p+1))
        return v
    def add_lune_at(self, theta, c):
        self.ops.append(('l', theta - math.pi, float(c)))
    def add_slit_to(self, p, over=1.0):
        zeta = complex(self._inv_raw(np.array([p]))[0])
        if abs(zeta) >= 1.0: return None
        self.ops.append(('s', cmath.phase(zeta) - math.pi, abs(zeta)*over))
        return abs(zeta)
    def excluded(self, p):
        z = complex(self._inv_raw(np.array([p]))[0])/self.rho
        return abs(z) >= 1.0 - 1e-9, abs(z)
