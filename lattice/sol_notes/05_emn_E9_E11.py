import importlib.util
from fractions import Fraction as F
from math import comb, lcm, log, gcd
import sympy as sp
spec=importlib.util.spec_from_file_location("fast","/home/ubuntu/code/math-modular-sources/lattice/sol_notes/05_emn_fast.py")
fast=importlib.util.module_from_spec(spec); spec.loader.exec_module(fast)
spec2=importlib.util.spec_from_file_location("core","/home/ubuntu/code/math-modular-sources/lattice/sol_notes/05_emn_core.py")
core=importlib.util.module_from_spec(spec2); spec2.loader.exec_module(core)

print("cross-validate fast engine vs sympy core engine:")
ok=True
for (m,t) in [(0,0),(2,1),(4,2),(6,3),(8,5),(4,-2),(6,-4),(10,10),(12,4)]:
    r,p,g = fast.I(m,t)
    rat,gc,pic = core.split(core.I(m,t))
    good = (sp.Rational(r.numerator,r.denominator)==rat and sp.Rational(g.numerator,g.denominator)==gc and p==0 and pic==0)
    ok &= good
    if not good: print("  MISMATCH",m,t,r,p,g,rat,gc,pic)
print("  all match, all pi-parts zero:", ok)

print()
print("="*72); print("E11: int_Delta (1 + 569 x^2y^2 - 2800 x^4y^4)/(1-x^2-y^2) dxdy =?= (13/2) G")
for k in [0,2,4]:
    r,p,g = fast.I(k,0)
    print(f"   Ihat({k}) = I({k},0) = {r} + ({p})*pi + ({g})*G   [log2 coordinate: 0 -- no log2 ever appears]")
tot=(F(0),F(0),F(0))
for c,k in [(1,0),(569,2),(-2800,4)]:
    tot=fast.add(tot, fast.smul(F(c), fast.I(k,0)))
print("   combination =", tot[0], "+", tot[1], "*pi +", tot[2], "*G      (13/2 =", F(13,2),")")

print()
print("="*72); print("E9: L_n = a_n + b_n G ;  256^n b_n =?= C(4n,2n) sum_j (-4)^j C(2n,j) C(4n,2n-2j)")
for n in range(0,6):
    r,p,g = fast.Lint(n)
    claim = comb(4*n,2*n)*sum((-4)**j*comb(2*n,j)*comb(4*n,2*n-2*j) for j in range(n+1))
    lhs = F(256)**n * g
    print(f"   n={n}: 256^n b_n = {lhs}   claimed = {claim}   match={lhs==claim}   (pi part {p})")
