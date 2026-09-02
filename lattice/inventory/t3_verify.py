"""Robust re-verification of the two saturation points with series to y^160."""
import t3_orbit as T
from t3_ode import find_relation
ONE = [1]+[0]*(T.Py-1)
Gd = [T.G]
for _ in range(6): Gd.append(T.d_(Gd[-1]))
Ints = [T.Int(T.G, i) for i in range(0, 6)]
print("\n" + "="*78)
print(f"ROBUST CHECK, series to y^{T.Py}")
for N in (3, 4, 5):
    d, nul = find_relation([ONE]+Gd[:N+1], dmax=40)
    print(f"  derivatives {{1,G..G^({N})}}: " + (f"RELATION deg<={d}" if d is not None else "INDEPENDENT"))
for I in (2, 3, 4):
    d, nul = find_relation([ONE]+Gd[:4]+Ints[:I+1], dmax=40)
    print(f"  4 derivs + {I+1} integrals    : " + (f"RELATION deg<={d}" if d is not None else "INDEPENDENT"))
d, _ = find_relation([ONE]+Gd[:4]+Ints[:3]+[T.Int(T.Int(T.G, 1), 1)], dmax=40)
print(f"  CDT's 8 + a DOUBLE integral   : " + (f"RELATION deg<={d}" if d is not None else "INDEPENDENT (e=2 cost)"))
