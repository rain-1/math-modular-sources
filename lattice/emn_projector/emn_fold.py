"""Numerical verification of the fold at z=2.

  scriptH(ang) = sum_m (-1)^m/(2m+1)^2 (1-cos((2m+1)ang)),  scriptH'' = (1/2)sec(ang)
  H(z) = scriptH(ang),  z = 1-cos(ang)
  H(2) = scriptH(pi) = 2G          [an honest Apery limit at the cusp z=2]
  Bfun(z) = int_0^z H(t) dt/(2-t) = int_0^ang scriptH(s) tan(s/2) ds
  Claim:  Bfun(z) + 2G log(2-z)  extends analytically across z=2.
"""
from mpmath import mp, mpf, catalan, cos, sin, tan, log, pi, quad, sec, mpmathify
mp.dps = 30
G = catalan

def scriptH(s):
    # scriptH(s) = G - (1/2)(Cl2(pi/2+s)+Cl2(pi/2-s)); use the series directly
    return quad(lambda u: quad(lambda t: 0.5/cos(t), [0, u]), [0, s]) if False else \
           mp.nsum(lambda m: (-1)**int(m)/(2*m+1)**2*(1-cos((2*m+1)*s)), [0, mp.inf])

def scriptH_fast(s):
    from mpmath import clsin
    # Cl2(x) = clsin(2, x)
    return G - (clsin(2, pi/2+s) + clsin(2, pi/2-s))/2

for s in (mp.mpf(0.5), pi/3, pi/2, 2*pi/3, pi):
    print(f"ang={float(s):.6f}  z={float(1-cos(s)):.6f}  scriptH={scriptH_fast(s)}")
print("checks:  H(1/2)=G/3 ->", scriptH_fast(pi/3), " G/3=", G/3)
print("         H(3/2)=5G/3->", scriptH_fast(2*pi/3), " 5G/3=", 5*G/3)
print("         H(1)=G     ->", scriptH_fast(pi/2), " G=", G)
print("         H(2)=2G    ->", scriptH_fast(pi), " 2G=", 2*G)

def Bfun(ang):
    return quad(lambda s: scriptH_fast(s)*tan(s/2), [0, ang])

print("\nfold test: Bfun(ang) + 2G*log(2-z),  z=1-cos(ang), ang -> pi")
prev = None
for d in (0.3, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005):
    ang = pi - d
    z = 1-cos(ang)
    val = Bfun(ang) + 2*G*log(2-z)
    print(f"  pi-ang={d:<7} z={float(z):.8f}  B+2G log(2-z) = {mp.nstr(val, 12)}"
          + (f"   diff={mp.nstr(val-prev,6)}" if prev is not None else ""))
    prev = val
print("\n(a finite smooth limit  <=>  B - 2G*A is regular at z=2, A=-log(1-z/2))")
