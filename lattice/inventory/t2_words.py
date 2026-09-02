"""Task 2 (decisive): ALL iterated integrals of weight <= 4 on P^1-{0,s,infty}
holomorphic at 0 (words in {omega_0=dx/x, omega_1=dx/(x-s)} ending in 1),
symmetrised to the y-line, with their sharp denominator types and the
Q(y)-dimension of the single-layer part.

Normalisation: the pure module must be INTEGRAL, i.e. Li_j(x/s) = sum lam^n x^n/n^j.
We therefore use the coordinate X = x/s throughout (so s=1) and record that the
transported host only rescales y -> y/s, which does not change any denominator."""
from fractions import Fraction as F
from math import comb
from qseries import mul, compose, w_series, to_y
from dtypes import divides_type

Px, Py = 170, 80
s = F(1)
W = w_series(s, Px)
xmw = [a-b for a, b in zip([F(0), F(1)]+[F(0)]*(Px-2), W)]

def int0(L):      # omega_0: int L dx/x
    return [F(0)]+[L[n]/n for n in range(1, Px)][:Px-1]
def int1(L):      # omega_1: int L dx/(x-1) = -int L (1+x+x^2+...) dx
    Q = [F(0)]*Px; acc = F(0)
    for n in range(Px): acc += L[n]; Q[n] = acc      # L/(1-x)
    return [F(0)]+[-Q[n]/(n+1) for n in range(Px-1)]

words = []
for wt in range(1, 5):
    for mask in range(1 << (wt-1)):
        w = tuple((mask >> i) & 1 for i in range(wt-1)) + (1,)
        words.append(w)
Lw = {}
for w in sorted(words, key=len):
    if len(w) == 1: Lw[w] = int1([F(1)]+[F(0)]*(Px-1))
    else: Lw[w] = (int0 if w[0] == 0 else int1)(Lw[w[1:]])

def pareto(c, maxlay=4, maxe=10):
    out = []
    for L in range(0, maxlay+1):
        for e in range(0, maxe+1):
            ok, _ = divides_type(c, e, (2,)*L, 1, len(c))
            if ok: out.append((L, e)); break
    return out

p = (1 << 61) - 1
def modp(c): return [((v.numerator % p)*pow(v.denominator % p, p-2, p)) % p for v in c]
def rank_mod(rows, ncol):
    rows = [r[:] for r in rows]; r = 0
    for c in range(ncol):
        piv = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        inv = pow(rows[r][c], p-2, p); rows[r] = [x*inv % p for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                f0 = rows[i][c]; rows[i] = [(x-f0*y) % p for x, y in zip(rows[i], rows[r])]
        r += 1
        if r == len(rows): break
    return r
def qdim(series_list, Dmax=3):
    rks = []
    for D in range(Dmax+1):
        L = min(Py-2, len(series_list)*(D+1)+12)
        rows = []
        for f in series_list:
            f = (f+[0]*L)[:L]
            for j in range(D+1): rows.append([0]*j+f[:L-j])
        rks.append(rank_mod(rows, L))
    diffs = [rks[i+1]-rks[i] for i in range(len(rks)-1)]
    return rks, (diffs[-1] if diffs else None)

print("word (1=dx/(x-1), 0=dx/x; read left=outermost) | Sym^+ type | (x-w)Sym^- type")
single_p, single_m, all_ser = [], [], []
for w in sorted(words, key=lambda z: (len(z), z)):
    Ls = Lw[w]
    LsW = compose(Ls, W, Px)
    sp = to_y([a+b for a, b in zip(Ls, LsW)], s, Py)
    sm = to_y(mul(xmw, [a-b for a, b in zip(Ls, LsW)], Px), s, Py)
    pp, pm = pareto(sp), pareto(sm)
    fmt = lambda P: (", ".join(f"L{L}:e={e}" for L, e in P[:2]) or "no type <=4 layers")
    z1 = all(v == 0 for v in sp); z2 = all(v == 0 for v in sm)
    print(f"  {''.join(map(str,w)):6s} wt={len(w)}  Sym+: {'IDENTICALLY ZERO' if z1 else fmt(pp):28s}"
          f"  (x-w)Sym-: {'IDENTICALLY ZERO' if z2 else fmt(pm)}")
    if not z1:
        all_ser.append((f"S+{''.join(map(str,w))}", sp))
        if pp and pp[0][0] <= 1: single_p.append((f"S+{''.join(map(str,w))}", pp[0], sp))
    if not z2:
        all_ser.append((f"S-{''.join(map(str,w))}", sm))
        if pm and pm[0][0] <= 1: single_m.append((f"S-{''.join(map(str,w))}", pm[0], sm))

print("\nSINGLE-LAYER ([1..2n]^1 or less) symmetrised words found:")
for nm, t, _ in single_p+single_m: print(f"   {nm}: layers={t[0]}, e={t[1]}")
ONE = [F(1)]+[F(0)]*(Py-1)
Fe = {e: [F(0)]+[F(1, n**e*comb(2*n, n)) for n in range(1, Py)] for e in range(1, 6)}
lst = [modp(ONE)]+[modp(c) for _, _, c in single_p+single_m]
rks, d = qdim(lst)
print(f"   Q(y)-dim of span{{1}} + all single-layer symmetrised words = {d}   (ranks {rks})")
lst2 = [modp(ONE), modp(Fe[1]), modp(Fe[2]), modp(Fe[3])]+[modp(c) for _, _, c in single_p+single_m]
rks2, d2 = qdim(lst2)
print(f"   Q(y)-dim of span{{1,F_1,F_2,F_3}} + those = {d2}   (ranks {rks2})  "
      f"-> {'NO new directions' if d2 == 4 else 'NEW directions!'}")
rks3, d3 = qdim([modp(ONE), modp(Fe[1]), modp(Fe[2]), modp(Fe[3])])
print(f"   Q(y)-dim of span{{1,F_1,F_2,F_3}} alone = {d3}   (ranks {rks3})")
print("\nALL symmetrised words wt<=4 (both signs), Q(y)-dim:")
rks4, d4 = qdim([modp(c) for _, c in all_ser])
print(f"   {d4}   (ranks {rks4}), from {len(all_ser)} nonzero series")
