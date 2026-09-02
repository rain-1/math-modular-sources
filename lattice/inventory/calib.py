"""Task 1: reproduce CDT's numbers exactly from my own (6.0.4)-(6.0.6)."""
import math
from fractions import Fraction as F
from tau import I_exact, tau_flat, tau_flat_direct, tau_sharp, tau

m = 14
cols = [(1, 2), (3, 2)]                                  # u_1=1, u_2=3, b_1=b_2=2
e = [0, 0, 1] + [0] * 6 + [1] * 5                        # (13.0.2)
assert len(e) == 14 and sum(e) == 6 and max(e) == 1

# sigma_i: rows 1 has sigma=0 (the constant), rows 2,3 have sigma=2, rows 4..14 sigma=4
sig = [0] + [2, 2] + [4] * 11
tf_direct = tau_flat_direct(m, sig)
tf = tau_flat(m, cols)
print(f"tau^flat closed form  = {tf}      = {float(tf):.6f}   (CDT 191/49)")
print(f"tau^flat from (2i-1)  = {tf_direct}      = {float(tf_direct):.6f}")
assert tf == tf_direct == F(191, 49)

I = I_exact(2, 14, 2)
print(f"I_2^14(2)             = {I} = {float(I):.6f}   (CDT: 21.075)")
assert I == F(843, 40), I

ts, xi = tau_sharp(m, e)
print(f"tau^sharp             = {ts} = {float(ts):.10f} at xi = {xi} = {float(xi):.6f}   (CDT 27/80)")
assert ts == F(27, 80), ts

T = tf + ts
print(f"tau(b;e)              = {T} = {float(T):.6f}   (CDT 16603/3920)")
assert T == F(16603, 3920)

psi = F(5448339453535586608 * 10**9, 8658833407565631122430056127)
lp = math.log(256 * float(psi))
print(f"|psi'(0)|             = {float(psi):.10f}")
print(f"log|phi'(0)|          = {lp:.6f}   (|phi'(0)| = {256*float(psi):.4f})")
BC = 11.845
den = lp - float(T)
print(f"entry  = log|phi'(0)|-tau = {den:.6f}")
print(f"bound  = BC/entry         = {BC/den:.4f}   (CDT 13.9938)")
print(f"margin = m*entry - BC     = {m*den - BC:+.4f}")
print()
print(f"CDT's realised contour loss  log 0.62922 = {math.log(float(psi)):.5f}")
print(f"hard ceiling for s=1: log 256 = {math.log(256):.4f}; entry at ceiling = {math.log(256)-float(T):+.4f}")
