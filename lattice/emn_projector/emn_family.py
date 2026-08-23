"""The optimal-geometry projector family: all Q-combinations F = sum c_nu * g_nu
whose singular set in z lies in [1,2), i.e. all poles with theta < 1/2 cancelled.
Also: the family with all poles theta < th* cancelled, for any th*.
"""
from fractions import Fraction as Fr
import math, sys, json
from emn_span import basis_residues, nullspace

def family(A, thstar=Fr(1,2), kinds=('0','p3')):
    res = basis_residues(A)
    keys = sorted([k for k in res if k[1] in kinds])
    thetas = sorted({t for k in keys for t in res[k]})
    rows = [[res[k].get(t, Fr(0)) for k in keys] for t in thetas if t < thstar]
    ns = nullspace(rows, len(keys)) if rows else \
         [[Fr(1) if i==q else Fr(0) for i in range(len(keys))] for q in range(len(keys))]
    out = []
    for v in ns:
        rr = {t: sum(v[i]*res[keys[i]].get(t, Fr(0)) for i in range(len(keys))) for t in thetas}
        nz = sorted([t for t in thetas if rr[t] != 0])
        out.append((v, nz, rr))
    return keys, out, thetas

if __name__ == '__main__':
    A = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    keys, fam, thetas = family(A)
    print(f"A={A}: dim of the theta>=1/2 family = {len(fam)} (of {len(keys)} basis elements)")
    nonconst = [f for f in fam if f[1]]
    print(f"   nonconstant members: {len(nonconst)}")
    for v, nz, rr in nonconst[:40]:
        supp = {f"{k[0]}{k[1]}": str(v[i]) for i,k in enumerate(keys) if v[i]!=0}
        zs = [round(1-math.cos(math.pi*float(t)),5) for t in nz]
        print("   c=", supp, " poles(theta)=", [str(t) for t in nz], " z=", zs)
    print("   constant (Psi=0) members:", len(fam)-len(nonconst))
