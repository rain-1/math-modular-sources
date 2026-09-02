"""21: pullback test.  A pullback C = R(x,y) h(alpha(x)beta(y)) satisfies a
FIRST-ORDER PDE  A(x) theta_x (C/R) = B(y) theta_y (C/R), i.e. a mixed
recurrence with a small shift set.  We search mod p for
   sum_{(i,j) in S} p_{ij}(a,b) c_{a+i,b+j} = 0,  S a small shift set,
   deg p_ij <= d.
Calibration: z3_D1 (a pullback) should have one; s10_D1 (not separable through
a geometric k-series) should not, at the same size."""
from lib2v import *
import sys
P = 1000003


def nullspace(rows, n, p=P):
    M = [r[:] for r in rows]
    piv, r = [], 0
    for c in range(n):
        pr = next((i for i in range(r, len(M)) if M[i][c] % p), None)
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        inv = pow(M[r][c], p-2, p)
        M[r] = [(v*inv) % p for v in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] % p:
                f = M[i][c]
                M[i] = [(M[i][j]-f*M[r][j]) % p for j in range(n)]
        piv.append(c)
        r += 1
    return n - len(piv)


def test(cf, name, S, d, A=26):
    tab = [[cf(a, b) % P for b in range(A+3)] for a in range(A+3)]
    mons = [(u, v) for u in range(d+1) for v in range(d+1) if u+v <= d]
    n = len(S)*len(mons)
    rows = []
    for a in range(1, A):
        for b in range(1, A):
            row = []
            for (i, j) in S:
                val = tab[a+i][b+j]
                for (u, v) in mons:
                    row.append((val*pow(a, u, P)*pow(b, v, P)) % P)
            rows.append(row)
    dim = nullspace(rows, n)
    print("   %-12s shifts %-28s deg<=%d : nullspace dim %d  (%d unknowns, %d eqns)"
          % (name, str(S), d, dim, n, len(rows)))
    return dim


S1 = [(0, 0), (1, 0), (0, 1), (1, 1)]
S2 = [(0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (0, 2)]
for d in (2, 3, 4):
    print("shift set S1 = %s" % S1)
    test(z3_D1, "z3_D1", S1, d)
    test(s10_D1, "s10_D1", S1, d)
    test(z3_D2, "z3_D2", S1, d)
    sys.stdout.flush()
print()
for d in (2, 3):
    print("shift set S2 = %s" % S2)
    test(z3_D1, "z3_D1", S2, d)
    test(s10_D1, "s10_D1", S2, d)
    test(z3_D2, "z3_D2", S2, d)

print()
print("low-degree (= low ODE-order) mixed relations: degree d in (a,b) <-> order d in theta")
for d in (0, 1, 2):
    for S in ([(0, 0), (1, 0), (0, 1)],
              [(0, 0), (1, 0), (0, 1), (1, 1)],
              [(0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (0, 2), (2, 1), (1, 2), (2, 2)]):
        print("  deg<=%d, |S|=%d:" % (d, len(S)), end=" ")
        for cf, nm in [(z3_D1, "z3_D1"), (s10_D1, "s10_D1"), (z3_D2, "z3_D2")]:
            tab = [[cf(a, b) % P for b in range(30)] for a in range(30)]
            mons = [(u, v) for u in range(d+1) for v in range(d+1) if u+v <= d]
            n = len(S)*len(mons)
            rows = []
            for a in range(1, 26):
                for b in range(1, 26):
                    row = []
                    for (i, j) in S:
                        val = tab[a+i][b+j]
                        for (u, v) in mons:
                            row.append((val*pow(a, u, P)*pow(b, v, P)) % P)
                    rows.append(row)
            print("%s=%d " % (nm, nullspace(rows, n)), end="")
        print()
