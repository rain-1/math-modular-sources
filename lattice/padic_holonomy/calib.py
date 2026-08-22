"""Calibration against CDT ICM 6.2: the circle and the lune for zeta_2(5)."""
import os, sys as _sys
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in _sys.path: _sys.path.insert(0, HERE)
import numpy as np, math, warnings
warnings.filterwarnings('ignore')
import haupt, family, outer
X = haupt.X0p(2)
L = 12*math.log(2); tau = 175/36.; m = 6

def go(name, specs, dr_exact, offs=0.5, Ns=(4000, 8000, 16000, 32000, 64000)):
    print("== %s   |psi'(0)| = %.12f" % (name, dr_exact))
    psi, dr = family.compose(specs)
    assert abs(dr-dr_exact) < 1e-12, (dr, dr_exact)
    for N in Ns:
        z = np.exp(2j*np.pi*(np.arange(N)+offs)/N)
        q = np.array([psi(zz) for zz in z], dtype=complex)
        g = haupt.logabs_x(X, q)
        RE = outer.rearr(g)
        ph = X.values_fast(q, 0.80)
        BC = outer.bc(z, ph)
        lg = math.log(dr); den = lg+L-tau
        print("  N=%6d  RE=%.7f (bound %.6f, margin %+.6f)   BC=%.7f (bound %.6f, margin %+.6f)  max|q|=%.9f"
              % (N, RE, (RE+L)/den, m*den-(RE+L), BC, (BC+L)/den, m*den-(BC+L), np.abs(q).max()))

go("Omega_circ  |z+2/5|<=3/5   (CDT: BC=2.13322, bound 4.43206, kappa 22.0724)",
   [('TANG', (np.pi, 0.6))], 1/3.)
go("Omega_lune  (2/3)*(-h(-z,5/2))  (CDT: 3.92881, bound 4.48866, kappa 19.7439)",
   [('SCALE', (2/3.,)), ('BITE', (0.0, 2.5))], 14/29.)
