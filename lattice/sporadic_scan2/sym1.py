#!/usr/bin/env python3
"""For each order-3 (Sym^2) Apery-like row A_n, form the power series square root
sum c_n t^n = sqrt(sum A_n t^n), find the smallest lambda with lambda^n c_n in Z,
fit the resulting Sym^1 recurrence, and report lambda_1, lambda_2, k, score, budget."""
from fractions import Fraction
import math, json, sys
import mpmath as mp

def az_row(a, b, c, N):
    A = [Fraction(1), Fraction(b)]
    for n in range(1, N):
        A.append((Fraction((2*n+1)*(a*n*n+a*n+b))*A[n] - Fraction(c*n**3)*A[n-1])/Fraction((n+1)**3))
    return A

def sqrt_series(A, N):
    c = [Fraction(1)]
    for n in range(1, N):
        s = sum(c[k]*c[n-k] for k in range(1, n))
        c.append((A[n]-s)/2)
    return c

def smallest_scale(c, upto=40, LIM=64):
    for lam in range(1, LIM+1):
        if all((Fraction(lam)**n*c[n]).denominator == 1 for n in range(upto)):
            return lam
    return None

def fit3(a, N):
    """fit (n+1)^2 a_{n+1} = (A n^2 + A' n + B) a_n - (C n^2 + C' n + C'') a_{n-1}"""
    import itertools
    rows = []
    for n in range(1, 24):
        rows.append([Fraction(n*n*a[n]), Fraction(n*a[n]), Fraction(a[n]),
                     Fraction(-n*n*a[n-1]), Fraction(-n*a[n-1]), Fraction(-a[n-1]),
                     Fraction((n+1)**2*a[n+1])])
    # solve 6 unknowns
    M = [r[:] for r in rows]
    nu = 6
    piv = []; r0 = 0
    for col in range(nu):
        p = next((i for i in range(r0, len(M)) if M[i][col] != 0), None)
        if p is None: continue
        M[r0], M[p] = M[p], M[r0]
        pv = M[r0][col]; M[r0] = [y/pv for y in M[r0]]
        for i in range(len(M)):
            if i != r0 and M[i][col] != 0:
                f = M[i][col]; M[i] = [y-f*z for y, z in zip(M[i], M[r0])]
        piv.append(col); r0 += 1
    if r0 < nu: return None
    sol = [Fraction(0)]*nu
    for i, pc in enumerate(piv): sol[pc] = M[i][nu]
    # check
    for n in range(1, min(len(a)-1, 200)):
        lhs = Fraction((n+1)**2)*a[n+1]
        rhs = (sol[0]*n*n+sol[1]*n+sol[2])*a[n] - (sol[3]*n*n+sol[4]*n+sol[5])*a[n-1]
        if lhs != rhs: return None
    return sol

NAMES = {(17,5,1): "Apery (17,5,1)  zeta(3)/6",
         (12,4,16): "T (12,4,16)  7zeta(3)/32",
         (10,4,64): "Domb (10,4,64)  7zeta(3)/24",
         (9,3,-27): "AZ (9,3,-27)  L(3,chi-3)/3",
         (11,5,125): "AZ (11,5,125)",
         (7,3,81): "AZ (7,3,81)"}
N = 260
out = []
for (a, b, c), nm in NAMES.items():
    A = az_row(a, b, c, N)
    cc = sqrt_series(A, N)
    lam = smallest_scale(cc)
    if lam is None:
        print(f"{nm}: no integral scale <= 64"); continue
    seq = [Fraction(lam)**n*cc[n] for n in range(N)]
    ok = all(x.denominator == 1 for x in seq)
    sol = fit3(seq, N)
    if sol is None:
        print(f"{nm}: lambda={lam} integral={ok} -- no (3-term, deg 2) recurrence"); continue
    Aq, Ap, Bq, Cq, Cp, Cpp = sol
    disc = float(Aq*Aq - 4*Cq)
    if disc >= 0:
        l1 = (float(Aq)+math.sqrt(disc))/2; l2 = (float(Aq)-math.sqrt(disc))/2
    else:
        l1 = l2 = math.sqrt(float(Cq))
    if abs(l2) > abs(l1): l1, l2 = l2, l1
    # companion
    bb = [Fraction(0), Fraction(1)]
    for n in range(1, N-1):
        bb.append(((Aq*n*n+Ap*n+Bq)*bb[n] - (Cq*n*n+Cp*n+Cpp)*bb[n-1])/Fraction((n+1)**2))
    dn = [1]*N
    for n in range(1, N): dn[n] = dn[n-1]*n//math.gcd(dn[n-1], n)
    k = None
    for kk in range(0, 6):
        if all((bb[n]*dn[n]**kk).denominator == 1 for n in range(1, min(N, 201))): k = kk; break
    mp.mp.dps = 120
    n0 = N-3
    lim = None
    if seq[n0] != 0 and abs(l1) > abs(l2)*1.0000001:
        lim = mp.mpf(bb[n0].numerator)/mp.mpf(bb[n0].denominator)/(mp.mpf(seq[n0].numerator)/mp.mpf(seq[n0].denominator))
    sc = (math.log(1/abs(l2)) - k) if (k is not None and l2 != 0) else None
    bu = (math.log(abs(l1)) - k) if k is not None else None
    print(f"{nm}\n   lambda={lam}  integral={ok}  a={[str(x) for x in seq[:7]]}")
    print(f"   rec: p1=({Aq})n^2+({Ap})n+({Bq})   p0=({Cq})n^2+({Cp})n+({Cpp})")
    print(f"   lam1={l1:.5f} lam2={l2:.5f} c={float(Cq):.4g} k={k} score={sc if sc is None else round(sc,4)} budget={bu if bu is None else round(bu,4)}")
    print(f"   limit={mp.nstr(lim,50) if lim is not None else None}")
    out.append({"name": nm, "lam": lam, "a": [str(x) for x in seq[:12]],
                "rec": [str(x) for x in sol], "lam1": l1, "lam2": l2, "c": float(Cq),
                "k": k, "score": sc, "budget": bu,
                "limit": mp.nstr(lim, 100) if lim is not None else None})
json.dump(out, open("sym1.json", "w"), indent=1)
