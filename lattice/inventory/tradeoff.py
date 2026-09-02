"""The single-layer trade-off: a functions at L=1 force e = 0,1,...,a-1."""
import math
from fasttau import tau_sharp_f
LOSS = math.log(0.6292232680)
print("Leanest inventory: 1 + a single-layer F_1..F_a (e=0..a-1) + 1 conditional generator, m=a+2")
print("| a | m | tau^flat | tau^# | tau | H1 entry@ceil | H2 entry@ceil | H2 entry@contour | H2 margin |")
print("|---|---|---|---|---|---|---|---|---|")
for a in range(0, 21):
    m = a+2
    tf = 4 - (2.0/(m*m))*(1+(1+a)**2)
    ts, _ = tau_sharp_f(m, a*(a-1)//2, max(0, a-1))
    T = tf+ts
    e1 = math.log(256)-T
    e2 = math.log(64)-T
    e2r = math.log(64)+LOSS-T
    print(f"| {a} | {m} | {tf:.4f} | {ts:.4f} | {T:.4f} | {e1:+.4f} | {e2:+.4f} | {e2r:+.4f} | "
          f"{m*e2r-(11.845+math.log(0.25)):+.3f} |")
print()
print("exact tau^sharp minimiser interval check for CDT's profile:")
import tau as TT
from fractions import Fraction as F
for xi in (F(2), F(13, 6), F(2, 1)+F(1, 12), F(9, 4)):
    v = F(2, 196)*(xi*6 + 1*TT.I_exact(xi, 14, xi))
    print(f"   xi={xi}: objective = {v} = {float(v):.10f}  ({'= 27/80' if v == F(27,80) else 'not min'})")
