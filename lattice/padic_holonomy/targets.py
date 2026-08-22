"""Targets: (name, haupt, denominator exponent d, m = d+1, log R_p)."""
import os, sys as _sys
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in _sys.path: _sys.path.insert(0, HERE)
import math
import haupt

LOG = math.log


def budget(d, m, L):
    """max allowed value of  cost = BC - m log|psi'(0)|  for a contradiction."""
    return d * (L - m + 1.0 / m)


TARGETS = []


def add(key, H, d, L, note):
    m = d + 1
    TARGETS.append(dict(key=key, H=H, d=d, m=m, L=L,
                        tau=d * (1 - 1.0 / m ** 2), budget=budget(d, m, L), note=note))


X02, X03, X05, X07 = haupt.X0p(2), haupt.X0p(3), haupt.X0p(5), haupt.X0p(7)
X09 = haupt.X0p2(3)
X14 = haupt.X14()

add("zeta_2(3)  (2,k=1)", X02, 3, 12 * LOG(2), "Calegari Thm 3.3 (calibration)")
add("zeta_2(5)  (2,k=2)", X02, 5, 12 * LOG(2), "CDT ICM 6.2 (calibration)")
add("zeta_2(7)  (2,k=3)", X02, 7, 12 * LOG(2), "new?")
add("zeta_3(3)  (3,k=1)", X03, 3, 6 * LOG(3), "Calegari Thm 3.4 (calibration)")
add("zeta_3(5)  (3,k=2)", X03, 5, 6 * LOG(3), "new?")
add("zeta_5(3)  (5,k=1)", X05, 3, 3 * LOG(5), "PRIMARY")
add("zeta_5(5)  (5,k=2)", X05, 5, 3 * LOG(5), "new?")
add("zeta_7(3)  (7,k=1)", X07, 3, 2 * LOG(7), "new?")
add("L_2(2,chi_-4) X_1(4)", X14, 2, 8 * LOG(2), "Calegari Thm 4.2 (calibration)")
add("3-adic X_0(9) (Zagier B)", X09, 2, 3 * LOG(3), "L_3(2,chi_-3), level 9 model")
add("L_3(2,chi_-3) X_0(3)", X03, 2, 6 * LOG(3), "level-3 model of the SAME value")

if __name__ == "__main__":
    print("%-26s %2s %2s %10s %10s %10s %10s" % ("target", "d", "m", "logR", "tau", "entryfl", "budget"))
    for T in TARGETS:
        # entry floor: minimal log|psi'(0)| for the denominator to be positive
        print("%-26s %2d %2d %10.6f %10.6f %10.6f %10.6f" %
              (T['key'], T['d'], T['m'], T['L'], T['tau'], T['tau'] - T['L'], T['budget']))
