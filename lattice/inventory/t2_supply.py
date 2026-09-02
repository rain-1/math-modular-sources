"""Task 2/4 input: the PURE SUPPLY FUNCTION N(L,e) = number of Q(y)-independent
pure functions on the descent orbifold with at most L layers of [1..2n] and
exponent at most e.  Computed from all symmetrised iterated integrals of
weight <= 5 on P^1-{0,s,infty} plus the arcsine tower F_1..F_6."""
from fractions import Fraction as F
from math import comb
from qseries import mul, compose, w_series, to_y
from dtypes import divides_type

Px, Py = 168, 80
s = F(1)
W = w_series(s, Px)
xmw = [a-b for a, b in zip([F(0), F(1)]+[F(0)]*(Px-2), W)]
def int0(L): return [F(0)]+[L[n]/n for n in range(1, Px)][:Px-1]
def int1(L):
    Q = [F(0)]*Px; acc = F(0)
    for n in range(Px): acc += L[n]; Q[n] = acc
    return [F(0)]+[-Q[n]/(n+1) for n in range(Px-1)]

WMAX = 5
words = []
for wt in range(1, WMAX+1):
    for mask in range(1 << (wt-1)):
        words.append(tuple((mask >> i) & 1 for i in range(wt-1))+(1,))
Lw = {}
for w in sorted(words, key=len):
    Lw[w] = int1([F(1)]+[F(0)]*(Px-1)) if len(w) == 1 else (int0 if w[0] == 0 else int1)(Lw[w[1:]])

def minL_e(c, maxlay=5, maxe=8):
    """Pareto: for each layer count L, the minimal e (or None)."""
    out = {}
    for L in range(0, maxlay+1):
        for e in range(0, maxe+1):
            ok, _ = divides_type(c, e, (2,)*L, 1, len(c))
            if ok: out[L] = e; break
    return out

items = []   # (name, series, {L: min e})
for w in sorted(words, key=lambda z: (len(z), z)):
    Ls = Lw[w]; LsW = compose(Ls, W, Px)
    for tag, ser in (('S+', to_y([a+b for a, b in zip(Ls, LsW)], s, Py)),
                     ('S-', to_y(mul(xmw, [a-b for a, b in zip(Ls, LsW)], Px), s, Py))):
        if all(v == 0 for v in ser): continue
        items.append((tag+''.join(map(str, w)), ser, minL_e(ser)))
for j in range(1, 7):
    c = [F(0)]+[F(1, n**j*comb(2*n, n)) for n in range(1, Py)]
    items.append((f'F_{j}', c, minL_e(c)))

p = (1 << 61) - 1
def modp(c): return [((v.numerator % p)*pow(v.denominator % p, p-2, p)) % p for v in c]
def rank_mod(rows, ncol):
    rows = [r[:] for r in rows]; r = 0
    for col in range(ncol):
        piv = next((i for i in range(r, len(rows)) if rows[i][col]), None)
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        inv = pow(rows[r][col], p-2, p); rows[r] = [x*inv % p for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][col]:
                f0 = rows[i][col]; rows[i] = [(x-f0*y) % p for x, y in zip(rows[i], rows[r])]
        r += 1
        if r == len(rows): break
    return r
def qdim(lst):
    rks = []
    for D in (2, 3):
        L = min(Py-2, len(lst)*(D+1)+10)
        if L <= len(lst)*(D+1): return None
        rows = []
        for f in lst:
            f = (list(f)+[0]*L)[:L]
            for j in range(D+1): rows.append([0]*j+f[:L-j])
        rks.append(rank_mod(rows, L))
    return rks[1]-rks[0]

ONE = [F(1)]+[F(0)]*(Py-1)
print(f"{len(items)} nonzero pure candidates (symmetrised words wt<=5, both signs, + F_1..F_6)")
print("\nPURE SUPPLY N(L,e) = Q(y)-dim of span{1} u {f : f has a type with <=L layers and exponent <=e}")
print("      " + "".join(f"  e<={e}" for e in range(0, 5)))
for L in range(0, 4):
    row = []
    for e in range(0, 5):
        sel = [modp(ONE)]+[modp(c) for nm, c, pr in items
                           if any(LL <= L and ee <= e for LL, ee in pr.items())]
        d = qdim(sel)
        row.append(d)
    print(f"  L={L}: " + "".join(f"{('%4s'%v):>6s}" for v in row))
print("\n(the L=0 row must be 1 everywhere: River's theorem u_1=1)")
print("\nper-candidate Pareto (name: L->e):")
for nm, c, pr in items:
    if min(pr.keys(), default=9) <= 2:
        print(f"   {nm:9s} " + "  ".join(f"L={L}:e={e}" for L, e in sorted(pr.items()) if L <= 3))
