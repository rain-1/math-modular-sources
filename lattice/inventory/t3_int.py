"""Task 3: does the INTEGRATION tower saturate?  int G dy/y^i for i=0,1,2,...
Each has e=1.  Search for a Q(y)-relation of any degree <= 20."""
import t3_orbit as T
from t3_ode import find_relation
ONE = [1]+[0]*(T.Py-1)
Gd = [T.G]
for _ in range(3): Gd.append(T.d_(Gd[-1]))
Ints = [T.Int(T.G, i) for i in range(0, 6)]
print()
print("=" * 78)
print("INTEGRATION tower: {1,G,G',G'',G'''} + int G dy/y^i, i=0..I")
for I in range(0, 6):
    fs = [ONE]+Gd+Ints[:I+1]
    d, nul = find_relation(fs, dmax=25)
    print(f"  I={I} ({5+I+1} functions): " +
          (f"RELATION at deg <= {d}" if d is not None else "NO relation (Q(y)-INDEPENDENT)"))
