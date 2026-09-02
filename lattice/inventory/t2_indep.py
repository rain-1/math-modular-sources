"""Task 2: Q(y)-linear independence of the single-layer candidates on the
descent orbifold, and identification of the dependencies (F_0 vs F_1 etc.)."""
from fractions import Fraction as F
from math import comb, factorial as fac
from qseries import mul, compose, w_series, to_y

Px, Py = 190, 92
s = F(1)
W = w_series(s, Px)
xmw = [a-b for a, b in zip([F(0), F(1)]+[F(0)]*(Px-2), W)]

def Li(j):
    return [F(0)]+[F(1, n**j) for n in range(1, Px)]
def symp(j):
    L = Li(j); return [a+b for a, b in zip(L, compose(L, W, Px))]
def symm(j):
    L = Li(j); return [a-b for a, b in zip(L, compose(L, W, Px))]

print("Sym+ Li_1 == 0 ?", all(v == 0 for v in symp(1)[:40]),
      "  <- the only lcm-free pure candidate vanishes identically  [PROVED: 1-w = 1/(1-x)]")

Fe = {e: [F(0)]+[F(1, n**e*comb(2*n, n)) for n in range(1, Py)] for e in range(0, 7)}
B2 = [F(0), F(0)]+[F(2*fac(n-2)*fac(n), fac(2*n)) for n in range(2, Py)]
B5 = [F(0)]+[F(fac(n-1)**2, fac(2*n-1)*(2*n-1)) for n in range(1, Py)]
L1 = to_y(mul(xmw, symm(1), Px), s, Py)       # (x-w) Sym^- Li_1     single layer, e=0
S2 = to_y(symp(2), s, Py)                     # Sym^+ Li_2           single layer, e=1
B4 = to_y(mul(xmw, symm(2), Px), s, Py)       # CDT's B_4            two layers, e=0
S3 = to_y(symp(3), s, Py)                     # Sym^+ Li_3           two layers, e=1
print("L1 = (x-w)Sym^-Li_1  first terms:", [str(v) for v in L1[1:5]])
print("S2 = Sym^+ Li_2      first terms:", [str(v) for v in S2[1:5]])
print("B4 = (x-w)Sym^-Li_2  first terms:", [str(v) for v in B4[1:5]], " (CDT: -4y+4y^2/9+31y^3/900)")

# ---- exact check of the identity F_0 = (2 F_1 + y)/(4-y)
lhs = [(4*Fe[0][n] - (Fe[0][n-1] if n >= 1 else 0)) for n in range(Py)]
rhs = [2*Fe[1][n] + (F(1) if n == 1 else F(0)) for n in range(Py)]
print("identity (4-y)F_0 = 2F_1 + y :", all(a == b for a, b in zip(lhs, rhs)))

p = (1 << 61) - 1
def modp(c): return [((v.numerator % p)*pow(v.denominator % p, p-2, p)) % p for v in c]
def rank_mod(rows, ncol):
    rows = [r[:] for r in rows]; r = 0
    for c in range(ncol):
        piv = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        inv = pow(rows[r][c], p-2, p)
        rows[r] = [x*inv % p for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                f0 = rows[i][c]; rows[i] = [(x-f0*y) % p for x, y in zip(rows[i], rows[r])]
        r += 1
        if r == len(rows): break
    return r

ONE = [F(1)]+[F(0)]*(Py-1)
pool = [('1', ONE), ('F_0', Fe[0]), ('F_1', Fe[1]), ('B_2', B2), ('L1', L1),
        ('F_2', Fe[2]), ('S2', S2), ('B_5', B5),
        ('F_3', Fe[3]), ('F_4', Fe[4]), ('F_5', Fe[5]), ('B_4', B4), ('S3', S3)]
P = {nm: modp(c) for nm, c in pool}

def test(names, D):
    L = min(Py-2, len(names)*(D+1)+14)
    rows = []
    for nm in names:
        f = (P[nm]+[0]*L)[:L]
        for j in range(D+1): rows.append([0]*j+f[:L-j])
    need = len(names)*(D+1)
    if L < need: return None, need
    return rank_mod(rows, L), need

print("\nQ(y)-independence, rank of {y^j f_i}_{j<=D} mod p=2^61-1 (series to y^%d):" % Py)
SETS = [
 ("SINGLE-LAYER e=0 pool: {1,F_0,F_1,B_2,L1}", ['1', 'F_0', 'F_1', 'B_2', 'L1']),
 ("  drop F_0:            {1,F_1,B_2,L1}",     ['1', 'F_1', 'B_2', 'L1']),
 ("  CDT's own e=0 pair:  {1,B_2}",            ['1', 'B_2']),
 ("SINGLE-LAYER e<=1:     {1,F_1,B_2,L1,F_2,S2}", ['1', 'F_1', 'B_2', 'L1', 'F_2', 'S2']),
 ("SINGLE-LAYER e<=4:     +{F_3,F_4,F_5}",     ['1', 'F_1', 'B_2', 'L1', 'F_2', 'S2', 'F_3', 'F_4', 'F_5']),
 ("all 13 candidates",                          [nm for nm, _ in pool]),
]
for lab, names in SETS:
    line = []
    for D in (0, 1, 2, 3):
        rk, need = test(names, D)
        line.append(f"D={D}: {rk}/{need}" if rk is not None else f"D={D}: n/a")
    print(f"  {lab:46s} " + "  ".join(line))
