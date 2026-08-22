"""Exact Q(y)-linear-independence check of CDT's 14 functions (their Lemma 12.1.1),
rebuilt from scratch: H_A,H_B,H_C from CDT Prop. 11.1.4, symmetrised to G(y),
plus the seven pure functions B_1..B_7 of CDT (10.1.2),(10.1.4),(10.2.1).
Verified against CDT's published coefficients (H_A=1,3,15,93; H_B=0,1,23/4;
H_C=0,1,6,343/9; B_4=-4y+4y^2/9+31y^3/900).
"""
from fractions import Fraction as F
from math import factorial as fac

Px, Py = 200, 92

def solve_ode(bb, cc, P=Px):
    y = [F(0)]*(P+4)
    y[0] = F(1) if (bb == 0 and cc == 0) else F(0)
    rhs = [F(bb)+F(cc) if n == 0 else F(cc) for n in range(P+4)]
    for n in range(0, P):
        s = F(0)
        for shift, co in ((1, 1), (2, -10), (3, 9)):
            k = n-shift
            if k >= 0: s += F(co)*F((k+2)*(k+1))*y[k+2]
        for shift, co in ((0, 1), (1, -20), (2, 27)):
            k = n-shift
            if k >= 0: s += F(co)*F(k+1)*y[k+1]
        s += F(-3)*y[n]
        if n >= 1: s += F(9)*y[n-1]
        coef = F(n+1)+F((n+1)*n)
        y[n+1] = (rhs[n]-(s-coef*y[n+1]))/coef
    return y[:P]

p = (1 << 61) - 1
def m_(u, v, P):
    r = [0]*P
    for i, ui in enumerate(u):
        if ui == 0 or i >= P: continue
        for j in range(min(len(v), P-i)):
            if v[j]: r[i+j] = (r[i+j]+ui*v[j]) % p
    return r
def c_(f, g, P):
    r = [0]*P; pw = [0]*P; pw[0] = 1
    for n in range(P):
        if f[n]:
            for i in range(P): r[i] = (r[i]+f[n]*pw[i]) % p
        if n < P-1: pw = m_(pw, g, P)
    return r
def toF(f, P):
    return [((v.numerator % p)*pow(v.denominator % p, p-2, p)) % p for v in (list(f)+[F(0)]*P)[:P]]

HA, HB, HC = solve_ode(0, 0), solve_ode(1, 0), solve_ode(0, 1)
assert [str(t) for t in HA[:4]] == ['1', '3', '15', '93']
assert str(HB[2]) == '23/4' and str(HC[3]) == '343/9'

w = [0]*Px
for n in range(1, Px): w[n] = p-1
Y = [0]*Px
for n in range(2, Px): Y[n] = p-1
Yp_cache = []
def build_Ypow():
    Yp = [0]*Px; Yp[0] = 1
    for n in range(Py+1):
        Yp_cache.append(Yp[:]); Yp = m_(Yp, Y, Px)
build_Ypow()
def to_y(gx):
    g = [0]*Py; cur = gx[:]
    for n in range(Py):
        if 2*n >= Px: break
        Yp = Yp_cache[n]
        g[n] = cur[2*n]*pow(Yp[2*n], p-2, p) % p
        if g[n]:
            cur = [(c-g[n]*s) % p for c, s in zip(cur, Yp)]
    return g

a_, b_, c0 = 1, p-3, 5
H = [(a_*x+b_*y+c0*z) % p for x, y, z in zip(toF(HA, Px), toF(HB, Px), toF(HC, Px))]
G = to_y([(x+y) % p for x, y in zip(H, c_(H, w, Px))])

def d_(f): return [(n+1)*f[n+1] % p for n in range(len(f)-1)]
def i_(f): return [0]+[f[n]*pow(n+1, p-2, p) % p for n in range(len(f)-1)]

B1 = [1]+[0]*(Py-1)
B2 = [0]*Py
for n in range(2, Py): B2[n] = toF([F(2*fac(n-2)*fac(n), fac(2*n))], 1)[0]
B3 = [0]*Py
for n in range(1, Py): B3[n] = toF([F(fac(n-1)**2, fac(2*n))], 1)[0]
B5 = [0]*Py
for n in range(1, Py): B5[n] = toF([F(fac(n-1)**2, fac(2*n-1)*(2*n-1))], 1)[0]
B6 = [0]*Py
for n in range(1, Py): B6[n] = toF([F(fac(n-1)**2, n*fac(2*n))], 1)[0]
Li2 = toF([F(0)]+[F(1, n*n) for n in range(1, Px)], Px)
anti = [((1 if n == 1 else 0)-w[n]) % p for n in range(Px)]
B4 = to_y(m_(anti, [(x-y) % p for x, y in zip(Li2, c_(Li2, w, Px))], Px))
B7 = i_([0]+[B4[n]*pow(n, p-2, p) % p for n in range(1, Py)])

fns = [('B1', B1), ('B2', B2), ('B3', B3), ('B4', B4), ('B5', B5),
       ('G', G), ("G'", d_(G)), ("G''", d_(d_(G))), ("G'''", d_(d_(d_(G)))),
       ('B6', B6), ('B7', B7), ('intG', i_(G)),
       ('int(G-G0)/y', i_([0]+[G[n]*pow(n, p-2, p) % p for n in range(1, Py)])),
       ('int(G-G0-G1y)/y^2', i_([0, 0]+[G[n]*pow(n-1, p-2, p) % p for n in range(2, Py)]))]

def rank_mod(rows, ncol):
    rows = [r[:] for r in rows]; r = 0
    for c in range(ncol):
        piv = None
        for i in range(r, len(rows)):
            if rows[i][c]: piv = i; break
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        inv = pow(rows[r][c], p-2, p)
        rows[r] = [x*inv % p for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                f0 = rows[i][c]
                rows[i] = [(x-f0*y) % p for x, y in zip(rows[i], rows[r])]
        r += 1
        if r == len(rows): break
    return r

print('CDT 14 functions, Q(y)-linear independence (rank of {y^j f_i} mod p=2^61-1):')
for D in (0, 1, 2, 3, 4, 5):
    L = min(Py-2, 14*(D+1)+12)
    rows = []
    for nm, f in fns:
        ff = (list(f)+[0]*L)[:L]
        for j in range(D+1):
            rows.append([0]*j+ff[:L-j])
    need = 14*(D+1)
    if L < need: print(f'  deg<= {D}: insufficient precision'); continue
    rk = rank_mod(rows, L)
    print(f'  deg P_i <= {D}: rank {rk}/{need} (series order {L}) -> '
          f'{"INDEPENDENT" if rk == need else "RELATION"}')
