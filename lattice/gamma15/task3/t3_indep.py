"""T3: K(y)-linear independence of the fourteen functions, on
  (a) CDT's own host (s=1, K=Q -- reproduction of CDT_FINDER section 2), and
  (b) Zagier row D on Gamma_1(5) over K=Q(sqrt5), s = -phi^5,
both at a SPLIT prime (Z[phi]/p = F_p) and at an INERT prime (Z[phi]/p = F_{p^2}).
Everything goes through the same code path; the host enters only through (h_n, s)."""
from fractions import Fraction as F
from math import factorial as fac, gcd
import sys
sys.path.insert(0, '/home/ubuntu/code/math-modular-sources/lattice/gamma15/task3')
from kfield import *

NY = 100          # Y-order kept for every function
NXuse = 2*NY

# ------------------------------------------------------------------ P_n(Y)
Pc = [[0]*(NY+2) for _ in range(NXuse+1)]
Pc[0][0] = 2; Pc[1][1] = 1
for n in range(2, NXuse+1):
    for j in range(1, NY+2):
        Pc[n][j] = Pc[n-1][j-1] - Pc[n-2][j-1]

# ------------------------------------------------------------------ pure B's
def pure_series(N):
    B1 = [F(0)]*(N+1); B1[0] = F(1)
    B2 = [F(0)]*(N+1)
    for n in range(2, N+1): B2[n] = F(2*fac(n-2)*fac(n), fac(2*n))
    B3 = [F(0)]*(N+1)
    for n in range(1, N+1): B3[n] = F(fac(n-1)**2, fac(2*n))
    B5 = [F(0)]*(N+1)
    for n in range(1, N+1): B5[n] = F(fac(n-1)**2, fac(2*n-1)*(2*n-1))
    B6 = [F(0)]*(N+1)
    for n in range(1, N+1): B6[n] = F(fac(n-1)**2, n*fac(2*n))
    C = [F(fac(2*k), fac(k)**2) for k in range(N+2)]
    B4 = [F(0)]*(N+1)
    for n in range(0, N):
        B4[n+1] = 4*sum(C[k]*C[n-k]*F(1, (2*k-1)*(2*n-2*k+1)**2) for k in range(n+1))/F(16)**n
    B7 = [F(0)]*(N+1)
    for n in range(1, N+1): B7[n] = B4[n]/n
    return B1, B2, B3, B4, B5, B6, B7
B1, B2, B3, B4, B5, B6, B7 = pure_series(NY)

# ------------------------------------------------------------------ CDT host
def solve_ode(bb, cc, P):
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
    return y[:P+1]

def load(fn):
    t = open(fn).read().strip(); t = t[t.index('[')+1:t.rindex(']')]
    return [F(u.strip()) for u in t.split(',')]

def host_CDT():
    P = NXuse
    HA, HB, HC = solve_ode(0, 0, P), solve_ode(1, 0, P), solve_ode(0, 1, P)
    assert [str(t) for t in HA[:4]] == ['1', '3', '15', '93']
    assert str(HB[2]) == '23/4' and str(HC[3]) == '343/9'
    H = [K(HA[n]-3*HB[n]+5*HC[n]) for n in range(P+1)]
    return H, ONE            # s = 1

def host_G15():
    A = load('A_coeffs.txt'); BD = load('BD_coeffs.txt')
    BR3 = load('BR3_coeffs.txt'); BR4 = load('BR4_coeffs.txt')
    a = K(1); b = kadd(K(-3), PHI); c = ksub(K(5), kscal(F(2), PHI))
    H = []
    for n in range(NXuse+1):
        Bnew = kadd(K(BR3[n]), kmul(PHI5, K(BR4[n])))
        H.append(kadd(kadd(kscal(A[n], a), kscal(BD[n], b)), kmul(c, Bnew)))
    return H, S

def symplus(H, s):
    hh = []; t = ONE
    for n in range(NXuse+1): hh.append(kmul(H[n], t)); t = kmul(t, s)
    G = [ZERO]*(NY+1)
    for N in range(NY+1):
        acc = ZERO
        for n in range(N, min(2*N, NXuse)+1):
            if Pc[n][N]: acc = kadd(acc, kscal(F(Pc[n][N]), hh[n]))
        G[N] = acc
    return G

# ------------------------------------------------------------------ 14 functions (exact, in K)
def dY(f):  return [kscal(F(n+1), f[n+1]) for n in range(len(f)-1)]
def iY(f):  return [ZERO]+[kscal(F(1, n+1), f[n]) for n in range(len(f)-1)]

def fourteen(G, variant='cdtdef'):
    KB = lambda B: [K(c) for c in B]
    G1, G2, G3 = dY(G), dY(dY(G)), dY(dY(dY(G)))
    f12 = iY(G)                                        # int G dy
    if variant == 'cdtdef':
        f13 = [ZERO]+[kscal(F(1, n), G[n]) for n in range(1, len(G))]           # int (G-G0)/y dy
        f14 = [ZERO]+[kscal(F(1, n), G[n+1]) for n in range(1, len(G)-1)]       # int (G-G0-G1 y)/y^2 dy
    else:                                              # exactly indep_check2.py's coding
        f13 = iY([ZERO]+[kscal(F(1, n), G[n]) for n in range(1, len(G))])
        f14 = iY([ZERO, ZERO]+[kscal(F(1, n-1), G[n]) for n in range(2, len(G))])
    return [('B1', KB(B1)), ('B2', KB(B2)), ('B3', KB(B3)), ('B4', KB(B4)), ('B5', KB(B5)),
            ('G', G), ("G'", G1), ("G''", G2), ("G'''", G3),
            ('B6', KB(B6)), ('B7', KB(B7)), ('intG', f12), ('int(G-G0)/y', f13),
            ('int(G-G0-G1y)/y^2', f14)]

# ------------------------------------------------------------------ finite fields
class Fsplit:
    """Z[phi]/p = F_p ; p = +-1 mod 5, sqrt5 = r."""
    def __init__(self, p, r):
        self.p, self.r = p, r
        assert (r*r-5) % p == 0
        self.name = f"SPLIT  p={p} (p mod 5 = {p%5}), sqrt5 = {r}"
    def red(self, a):
        p = self.p
        u, v = a
        f = lambda q: (q.numerator % p)*pow(q.denominator % p, p-2, p) % p
        return (f(u)+self.r*f(v)) % p
    zero = 0
    def mul(self, x, y): return x*y % self.p
    def sub(self, x, y): return (x-y) % self.p
    def inv(self, x): return pow(x, self.p-2, self.p)
    def nz(self, x): return x != 0
class Finert:
    """Z[phi]/p = F_{p^2} = F_p[t]/(t^2-5) ; p = +-2 mod 5."""
    def __init__(self, p):
        self.p = p
        assert p % 5 in (2, 3)
        assert pow(5, (p-1)//2, p) == p-1, "5 must be a non-residue"
        self.name = f"INERT  p={p} (p mod 5 = {p%5}), F_p[t]/(t^2-5)"
    def red(self, a):
        p = self.p
        f = lambda q: (q.numerator % p)*pow(q.denominator % p, p-2, p) % p
        return (f(a[0]), f(a[1]))
    zero = (0, 0)
    def mul(self, x, y):
        p = self.p
        return ((x[0]*y[0]+5*x[1]*y[1]) % p, (x[0]*y[1]+x[1]*y[0]) % p)
    def sub(self, x, y): return ((x[0]-y[0]) % self.p, (x[1]-y[1]) % self.p)
    def inv(self, x):
        p = self.p
        n = pow((x[0]*x[0]-5*x[1]*x[1]) % p, p-2, p)
        return (x[0]*n % p, (-x[1])*n % p)
    def nz(self, x): return x != (0, 0)

def rank_mod(rows, ncol, FF):
    rows = [r[:] for r in rows]; r = 0
    for c in range(ncol):
        piv = None
        for i in range(r, len(rows)):
            if FF.nz(rows[i][c]): piv = i; break
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        iv = FF.inv(rows[r][c])
        rows[r] = [FF.mul(x, iv) for x in rows[r]]
        for i in range(len(rows)):
            if i != r and FF.nz(rows[i][c]):
                f0 = rows[i][c]
                rows[i] = [FF.sub(x, FF.mul(f0, y)) for x, y in zip(rows[i], rows[r])]
        r += 1
        if r == len(rows): break
    return r

def table(fns, FF, tag):
    red = [(nm, [FF.red(c) for c in f]) for nm, f in fns]
    print(f"  [{tag}] {FF.name}")
    out = []
    for D in (0, 1, 2, 3, 4, 5):
        L = min(min(len(f) for _, f in red)-2, 14*(D+1)+12)
        need = 14*(D+1)
        if L < need:
            print(f"    deg <= {D}: insufficient precision (L={L} < {need})"); continue
        rows = []
        for nm, f in red:
            ff = (list(f)+[FF.zero]*L)[:L]
            for j in range(D+1): rows.append([FF.zero]*j+ff[:L-j])
        rk = rank_mod(rows, L, FF)
        out.append((D, rk, need, L))
        print(f"    deg P_i <= {D}: rank {rk}/{need} (series order {L}) -> "
              f"{'INDEPENDENT' if rk == need else 'RELATION'}")
    return out

# ------------------------------------------------------------------ primes
def find_prime(cond, start):
    from sympy import nextprime
    p = start
    while True:
        p = int(nextprime(p))
        if cond(p): return p
def sqrt5_mod(p):
    from sympy.ntheory.residue_ntheory import sqrt_mod
    return int(sqrt_mod(5, p))

P_SPLIT = find_prime(lambda p: p % 5 in (1, 4), (1 << 61))
P_INERT = find_prime(lambda p: p % 5 in (2, 3), (1 << 61))
FS = Fsplit(P_SPLIT, sqrt5_mod(P_SPLIT))
FI = Finert(P_INERT)

print("="*78)
print("(a) CDT's own host  (s = 1, K = Q)  -- reproduction of CDT_FINDER section 2")
Hc, sc = host_CDT()
Gc = symplus(Hc, sc)
Ga = symplus([K(c) for c in solve_ode(0, 0, NXuse)], ONE)
print("  G_A(y) head =", [str(Ga[n][0]) for n in range(4)], " (CDT: 2, -27, 1014, -49536)")
for var in ('cdtdef', 'indep_check2'):
    print(f"  --- item 13,14 convention: {var} ---")
    fns = fourteen(Gc, var)
    table(fns, FS, 'F_p'); table(fns, FI, 'F_p^2')

print("="*78)
print("(b) Zagier row D on Gamma_1(5), K = Q(sqrt5), s = -phi^5")
Hg, sg = host_G15()
Gg = symplus(Hg, sg)
GA15 = symplus([K(c) for c in load('A_coeffs.txt')[:NXuse+1]], S)
print("  G_A(Y) head =", [str(GA15[n][0])+('+'+str(GA15[n][1])+'sqrt5' if GA15[n][1] else '') for n in range(5)])
for var in ('cdtdef', 'indep_check2'):
    print(f"  --- item 13,14 convention: {var} ---")
    fns = fourteen(Gg, var)
    table(fns, FS, 'F_p')
    table(fns, FI, 'F_p^2')
