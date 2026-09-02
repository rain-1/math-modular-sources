"""Task 3: the theta-orbit of CDT's conditional generator G = Sym^+ H.
How many Q(y)-independent members do j derivatives and i integrations give?"""
from fractions import Fraction as F
from math import comb, factorial as fac

Px, Py = 340, 160
p = (1 << 61) - 1

def solve_ode(bb, cc, P=Px):
    """CDT Prop 11.1.4: x(1-x)(1-9x)y''+(1-20x+27x^2)y'+3(3x-1)y = b + c/(1-x)."""
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

def toF(f, P): return [((v.numerator % p)*pow(v.denominator % p, p-2, p)) % p for v in (list(f)+[F(0)]*P)[:P]]
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
HA, HB, HC = solve_ode(0, 0), solve_ode(1, 0), solve_ode(0, 1)
assert [str(t) for t in HA[:4]] == ['1', '3', '15', '93'] and str(HC[3]) == '343/9'
w = [0]*Px
for n in range(1, Px): w[n] = p-1
Y = [0]*Px
for n in range(2, Px): Y[n] = p-1
Yp_cache = []
Yp = [0]*Px; Yp[0] = 1
for n in range(Py+1):
    Yp_cache.append(Yp[:]); Yp = m_(Yp, Y, Px)
def to_y(gx):
    g = [0]*Py; cur = gx[:]
    for n in range(Py):
        if 2*n >= Px: break
        g[n] = cur[2*n]*pow(Yp_cache[n][2*n], p-2, p) % p
        if g[n]: cur = [(c-g[n]*s) % p for c, s in zip(cur, Yp_cache[n])]
    return g
H = [(1*a+(p-3)*b+5*c) % p for a, b, c in zip(toF(HA, Px), toF(HB, Px), toF(HC, Px))]
G = to_y([(a+b) % p for a, b in zip(H, c_(H, w, Px))])

def d_(f): return [(n+1)*f[n+1] % p for n in range(len(f)-1)]
def Int(f, i):
    """int (f - Taylor_{<i}) dy / y^i  -> coefficients f_n/(n-i+1), n>=i (i>=1); i=0: int f dy."""
    if i == 0: return [0]+[f[n]*pow(n+1, p-2, p) % p for n in range(len(f)-1)]
    out = [0]*len(f)
    for n in range(i, len(f)):
        if n-i+1 != 0: out[n-i+1] = f[n]*pow(n-i+1, p-2, p) % p
    return out
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
def qdim(lst, Dmax=4, Ly=Py):
    rks = []
    for D in range(Dmax+1):
        L = min(Ly-2, len(lst)*(D+1)+14)
        rows = []
        for f in lst:
            f = (list(f)+[0]*L)[:L]
            for j in range(D+1): rows.append([0]*j+f[:L-j])
        rks.append(rank_mod(rows, L))
    return rks, rks[-1]-rks[-2]

ONE = [1]+[0]*(Py-1)
Gd = [G]
for _ in range(7): Gd.append(d_(Gd[-1]))
print("Q(y)-dim of the DERIVATIVE orbit of G (with 1 adjoined):  [ranks D=0..4, then dim]")
for j in range(0, 8):
    lst = [ONE]+Gd[:j+1]
    rks, d = qdim(lst)
    print(f"   {{1,G,...,G^({j})}}  ({j+1} derivs): ranks {rks}  ->  Q(y)-dim {d}"
          f"   {'(saturated)' if d < j+2 else ''}")
print()
print("Q(y)-dim of {1,G,G',G'',G'''} + integrations int G dy/y^i:")
Ints = [Int(G, i) for i in range(0, 7)]
base = [ONE]+Gd[:4]
for i in range(0, 7):
    lst = base+Ints[:i+1]
    rks, d = qdim(lst)
    print(f"   + int dy/y^0..y^{i} ({i+1} integrals): ranks {rks}  ->  Q(y)-dim {d}"
          f"   (CDT use i=0,1,2)")
