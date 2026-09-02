"""CDT score for Apery's host: the actual k=3 row and the hypothetical k=2 magnetic row."""
import sys, os, math, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'cdt_finder'))
from hosts import cdt_model, LOG_PHI_CDT, BC_CDT

print("calibration: CDT's own L(2,chi_-3) proof, k=2, p0=7, lam2=1")
r = cdt_model(2, 7, 1)
for key in ("m","k","lam2","u","tau","logphi","bc","entry","bound","margin"):
    print("   %-8s %s" % (key, r[key]))
print()
print("Apery host (lam2^norm = |c|^{1/2} = 1), p0 = 7 (CDT's own inventory):")
for k in (3, 2):
    r = cdt_model(k, 7, 1)
    print("   k=%d : tau=%.6f  entry=%.6f  BC/entry=%.4f  margin=%.4f  %s"
          % (k, r['tau'], r['entry'], r['bound'], r['margin'],
             "PASSES (margin>0)" if r['margin'] > 0 else "fails"))
print()
print("scan over the orbit size p0 (m = 2 p0), lam2 = 1:")
print("   p0 |   k=3 entry  margin  |   k=2 entry  margin")
for p0 in range(2, 21):
    a = cdt_model(3, p0, 1); b = cdt_model(2, p0, 1)
    print("   %2d | %10.5f %8.3f | %10.5f %8.3f" % (p0, a['entry'], a['margin'], p0 and b['entry'], b['margin']))
print()
print("the classical Apery-style irrationality exponent from the same data:")
import math
lam1 = 17 + 12*math.sqrt(2); lam2 = 17 - 12*math.sqrt(2)
L = math.log(lam1)
for k in (3, 2):
    num = k + L; den = L - k
    print("   k=%d : |q_n| ~ e^{%.5f n}, |q_n xi - p_n| ~ e^{%.5f n}  -> mu <= %s"
          % (k, num, -den, ("%.5f" % (1 + num/den)) if den > 0 else "no irrationality"))
