from fractions import Fraction as F
from math import comb, factorial as fac
from qseries import mul, compose, w_series, to_y
Px, Py = 168, 80
s = F(1); W = w_series(s, Px)
xmw = [a-b for a, b in zip([F(0), F(1)]+[F(0)]*(Px-2), W)]
def Li(j): return [F(0)]+[F(1, n**j) for n in range(1, Px)]
B4 = to_y(mul(xmw, [a-b for a, b in zip(Li(2), compose(Li(2), W, Px))], Px), s, Py)
Fe = {j: [F(0)]+[F(1, n**j*comb(2*n, n)) for n in range(1, Py)] for j in (1, 2, 3)}
B5 = [F(0)]+[F(fac(n-1)**2, fac(2*n-1)*(2*n-1)) for n in range(1, Py)]
ONE = [F(1)]+[F(0)]*(Py-1)
p = (1 << 61)-1
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
def qd(lst0, Ds=(2, 3, 4)):
    lst = [modp(c) for c in lst0]
    rks = []
    for D in Ds:
        L = min(Py-2, len(lst)*(D+1)+10)
        if L <= len(lst)*(D+1): return None, None
        rows = []
        for f in lst:
            f = (list(f)+[0]*L)[:L]
            for j in range(D+1): rows.append([0]*j+f[:L-j])
        rks.append(rank_mod(rows, L))
    return rks, rks[-1]-rks[-2]
for lab, lst in (("{1,F_1,F_2,B_4}", [ONE, Fe[1], Fe[2], B4]),
                 ("{1,F_1,F_2,B_4,B_5}", [ONE, Fe[1], Fe[2], B4, B5]),
                 ("{1,F_1,F_2,F_3,B_4,B_5}", [ONE, Fe[1], Fe[2], Fe[3], B4, B5])):
    rks, d = qd(lst)
    print(f"  {lab:28s} ranks {rks} -> Q(y)-dim {d}")
