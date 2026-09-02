"""T4: exact tau^flat, tau^sharp, tau for the MEASURED array."""
from fractions import Fraction as Fr
import sys
sys.path.insert(0, '/home/ubuntu/code/math-modular-sources/lattice/cdt_finder')
import cdt_bound
m = 14
cols = [(1, 2), (3, 2)]
e = [0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]
sm, tf = cdt_bound.tau_flat(m, cols)
print("m = 14, r = 2, b_1 = b_2 = 2, u_1 = 1, u_2 = 3")
print("e =", e, " sum e =", sum(e), " max e =", max(e))
print(f"sigma_m   = {sm}")
print(f"tau^flat  = sigma_m - (u_1^2 b_1 + u_2^2 b_2)/m^2 = 4 - (2+18)/196 = {tf} = {float(tf):.12f}")
I = cdt_bound.I_uvw(2, 14, 2)
print(f"I_2^14(2) = {I}   (843/40 = {843/40})   equal: {abs(I-843/40) < 1e-9}")
for xi in (Fr(2), Fr(13, 6), Fr(2, 1)+Fr(1, 12)):
    v = Fr(2, m*m)*(xi*sum(e)) + Fr(2, m*m)*max(e)*Fr(cdt_bound.I_uvw(float(xi), 14, float(xi))).limit_denominator(10**7)
    print(f"  xi = {xi}: 2/m^2 (xi*sum e + max e * I_xi^14(xi)) = {float(v):.12f}")
ts_exact = Fr(2, m*m)*(2*sum(e) + max(e)*Fr(843, 40))
print(f"tau^sharp = 2/196 * (2*6 + 843/40) = {ts_exact} = {float(ts_exact):.12f}   (CDT: 27/80)")
tau = tf + ts_exact
print(f"tau       = {tf} + {ts_exact} = {tau} = {float(tau):.12f}")
print(f"CDT's value 16603/3920 = {Fr(16603,3920)} ; equal: {tau == Fr(16603,3920)}")
