"""The EMN projector span: monodromy residue calculus and the cancellation
linear systems.

Basis (after the Niven collapse, see EMN_PROJECTOR.md sec.2):
   g[a,'0']   = script-H(a*ang)                       (shift b=0)
   g[a,'p3']  = script-H(a*ang+pi/3)+script-H(a*ang-pi/3)   (shift b=pi/3)
Every rational-coefficient element of the span is a Q-combination of these.

Singular points: ang = pi*th with th in (0,1); z = 1-cos(pi*th).
Residue of Psi=d^2F/d(ang)^2 at ang=pi*th (in the ang variable) is rational.
"""
from fractions import Fraction as Fr
import math, itertools, json, sys

def basis_residues(A):
    """returns dict key->{theta: residue}, theta a Fraction in (0,1)."""
    out = {}
    for a in range(1, A+1):
        d = {}
        for k in range(0, a):
            th = Fr(2*k+1, 2*a)
            if 0 < th < 1: d[th] = d.get(th, Fr(0)) + Fr(a,2)*(-1)**(k+1)
        out[(a,'0')] = d
        d = {}
        for k in range(0, a+1):
            for sgn, off in ((+1, Fr(1,6)), (-1, Fr(5,6))):
                th = Fr(k,1)/Fr(a) + off/Fr(a) if False else (Fr(k)+off)/Fr(a)
                if 0 < th < 1:
                    d[th] = d.get(th, Fr(0)) + Fr(a,2)*(-1)**(k+1)
        out[(a,'p3')] = d
    return out

# ---- exact rational nullspace
def nullspace(M, ncols):
    """M: list of rows (lists of Fr).  Returns basis of nullspace as list of vectors."""
    M = [row[:] for row in M]
    nrows = len(M)
    piv = []
    r = 0
    for c in range(ncols):
        p = None
        for i in range(r, nrows):
            if M[i][c] != 0: p = i; break
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        pv = M[r][c]
        M[r] = [x/pv for x in M[r]]
        for i in range(nrows):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [x - f*y for x, y in zip(M[i], M[r])]
        piv.append(c); r += 1
        if r == nrows: break
    free = [c for c in range(ncols) if c not in piv]
    basis = []
    for fc in free:
        v = [Fr(0)]*ncols
        v[fc] = Fr(1)
        for i, c in enumerate(piv):
            v[c] = -M[i][fc]
        basis.append(v)
    return basis

def sweep(A, kinds=('0','p3'), verbose=True):
    res = basis_residues(A)
    keys = [k for k in res if k[1] in kinds]
    keys.sort()
    thetas = sorted({t for k in keys for t in res[k]})
    idx = {k:i for i,k in enumerate(keys)}
    n = len(keys)
    best = None
    results = []
    for j, thstar in enumerate(thetas):
        # cancel all poles with theta < thstar  (strictly)
        rows = []
        for t in thetas:
            if t < thstar:
                rows.append([res[k].get(t, Fr(0)) for k in keys])
        ns = nullspace(rows, n) if rows else [[Fr(1) if i==q else Fr(0) for i in range(n)] for q in range(n)]
        # among nullspace vectors, is there one with a nonzero residue somewhere?
        alive = []
        for v in ns:
            rr = {t: sum(v[i]*res[keys[i]].get(t, Fr(0)) for i in range(n)) for t in thetas}
            if any(x != 0 for x in rr.values()):
                alive.append((v, rr))
        if alive:
            # theta_min actually achieved
            tmins = []
            for v, rr in alive:
                nz = [t for t in thetas if rr[t] != 0]
                tmins.append(min(nz))
            results.append((thstar, len(ns), len(alive), max(tmins)))
            if best is None or max(tmins) > best[0]:
                # record an explicit vector realising it
                for v, rr in alive:
                    nz = [t for t in thetas if rr[t] != 0]
                    if min(nz) == max(tmins):
                        best = (min(nz), v, rr, keys); break
        else:
            results.append((thstar, len(ns), 0, None))
    return best, results, keys, thetas

if __name__ == '__main__':
    A = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    for kinds in (('0',), ('0','p3')):
        best, results, keys, thetas = sweep(A, kinds)
        if best:
            th, v, rr, keys = best
            z = 1-math.cos(math.pi*float(th))
            print(f"A={A} kinds={kinds}: best theta_min={th} -> nearest pole z={z:.6f}")
            print("   combination:", {keys[i]: str(v[i]) for i in range(len(keys)) if v[i]!=0})
            print("   surviving poles:", {str(t): str(rr[t]) for t in thetas if rr[t]!=0})
        else:
            print(f"A={A} kinds={kinds}: nothing")
