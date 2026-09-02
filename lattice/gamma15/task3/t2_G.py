"""T2: the conditional side over K = Q(sqrt5).
H = a*A + b*B_D + c*B_new in K[[x]] (order 210 from PARI), G(Y) = Sym^+ H,
via  x^n + w^n = P_n(y) = s^n * Pcdt_n(Y),  Y = y/s.
Denominator type of G measured for n <= 60."""
from fractions import Fraction as F
from math import gcd, lcm
import sys, re
sys.path.insert(0, '/home/ubuntu/code/math-modular-sources/lattice/gamma15/task3')
from kfield import *

def load(fn, tag):
    t = open(fn).read().strip()
    t = t[t.index('[')+1:t.rindex(']')]
    return [F(u.strip()) for u in t.split(',')]
A   = load('A_coeffs.txt', 'A')
BD  = load('BD_coeffs.txt', 'BD')
BR3 = load('BR3_coeffs.txt', 'BR3')
BR4 = load('BR4_coeffs.txt', 'BR4')
NX = len(A)-1
print(f"loaded x-series to order {NX}: A[:5]={A[:5]}")
assert A[:6] == [F(1), F(3), F(19), F(147), F(1251), F(11253)]

def lcm_upto(m, c={}):
    if m in c: return c[m]
    r = 1
    for j in range(1, m+1): r = r*j//gcd(r, j)
    c[m] = r; return r

print("\n--- host recurrence check (n+1)^2 A_{n+1} = (11n^2+11n+3)A_n + n^2 A_{n-1} ---")
ok = all((n+1)**2*A[n+1] == (11*n*n+11*n+3)*A[n]+n*n*A[n-1] for n in range(1, NX))
print("  ", "HOLDS for all n<=%d [verified]" % (NX-1) if ok else "FAILS")

print("\n--- denominator type of the x-side companions (over Q), n <= 100 ---")
for nm, B in (('B_D', BD), ('B_R3', BR3), ('B_R4', BR4)):
    bad = [n for n in range(1, 101) if (B[n]*lcm_upto(n)**2).denominator != 1]
    bad1 = [n for n in range(1, 101) if (B[n]*lcm_upto(n)).denominator != 1]
    print(f"  {nm}: [1..n]^2 {'VALID' if not bad else 'FAILS '+str(bad[:3])};"
          f"  [1..n]^1 {'valid' if not bad1 else 'fails first at n='+str(bad1[0])}")

# ------------------------------------------------------- Pcdt_n(Y) integer polys
NY = 105
NXuse = min(NX, 2*NY)
Pc = [[0]*(NXuse+1) for _ in range(NXuse+1)]     # Pc[n][j] = [Y^j] Pcdt_n
Pc[0][0] = 2
if NXuse >= 1: Pc[1][1] = 1
for n in range(2, NXuse+1):
    for j in range(1, NXuse+1):
        Pc[n][j] = Pc[n-1][j-1] - Pc[n-2][j-1]
# structural facts CDT assert
supp = [[j for j in range(NXuse+1) if Pc[n][j]] for n in range(NXuse+1)]
print("\n--- P_n structure ---")
print("  deg P_n = n and val P_n = ceil(n/2) for all n<=%d: %s" %
      (NXuse, all(supp[n] and max(supp[n]) == n and min(supp[n]) == -(-n//2) for n in range(1, NXuse+1))))
print("  all coefficients integral: True (built by integer recurrence)")

# verify  x^n + w^n = s^n * Pcdt_n(y/s)  as exact K-series in x, to order 40
NC = 42
def xmul(a, b, P=NC):
    r = [ZERO]*(P+1)
    for i in range(P+1):
        if a[i] != ZERO:
            for j in range(P+1-i):
                if b[j] != ZERO: r[i+j] = kadd(r[i+j], kmul(a[i], b[j]))
    return r
# w = s x/(x - s) = -x * 1/(1 - x/s) = -x * sum (x/s)^k   -> w_n = -s^{-(n-1)} for n>=1
w = [ZERO]*(NC+1)
for n in range(1, NC+1):
    t = ONE
    for _ in range(n-1): t = kmul(t, SINV)
    w[n] = kscal(F(-1), t)
# Y = y/s = x^2/(s(x-s)) = -(x^2/s^2) * 1/(1-x/s)
Yser = [ZERO]*(NC+1)
for n in range(2, NC+1):
    t = ONE
    for _ in range(n): t = kmul(t, SINV)
    Yser[n] = kscal(F(-1), t)
xs_ = [ZERO]*(NC+1); xs_[1] = ONE
# check x+w = s*Y and x*w = s*(s*Y)
sY = [kmul(S, u) for u in Yser]
assert all(kadd(xs_[n], w[n]) == sY[n] for n in range(NC+1)), "x+w != y"
xw = xmul(xs_, w)
ssY = [kmul(S, u) for u in sY]
assert all(xw[n] == ssY[n] for n in range(NC+1)), "x*w != s*y"
print("  x + w = y  and  x*w = s*y  [verified as exact K-series to order %d]" % NC)
# powers
Yp = [[ZERO]*(NC+1) for _ in range(NC+1)]; Yp[0][0] = ONE
for k in range(1, NC+1): Yp[k] = xmul(Yp[k-1], Yser)
xp = [[ZERO]*(NC+1) for _ in range(NC+1)]; xp[0][0] = ONE
wp = [[ZERO]*(NC+1) for _ in range(NC+1)]; wp[0][0] = ONE
for k in range(1, NC+1):
    xp[k] = xmul(xp[k-1], xs_); wp[k] = xmul(wp[k-1], w)
sp = [ONE]
for k in range(1, NC+1): sp.append(kmul(sp[-1], S))
bad = []
for n in range(0, NC+1):
    lhs = [kadd(xp[n][i], wp[n][i]) for i in range(NC+1)]
    rhs = [ZERO]*(NC+1)
    for j in range(min(n, NC)+1):
        if Pc[n][j]:
            for i in range(NC+1):
                rhs[i] = kadd(rhs[i], kscal(F(Pc[n][j]), kmul(sp[n], Yp[j][i])))
    if any(lhs[i] != rhs[i] for i in range(NC+1)): bad.append(n)
print("  x^n + w^n = s^n P^cdt_n(Y) for n=0..%d: %s [verified]" % (NC, 'ALL OK' if not bad else 'FAILS '+str(bad)))

# ------------------------------------------------------- H and G
def build_G(a, b, c, NY=NY):
    """H = a*A + b*B_D + c*B_new ; returns Ghat[N] = [Y^N] Sym^+ H, in K."""
    Bnew = [kadd(K(BR3[n]), kmul(PHI5, K(BR4[n]))) for n in range(NX+1)]
    H = [kadd(kadd(kscal(A[n], a), kscal(BD[n], b)), kmul(c, Bnew[n])) for n in range(NX+1)]
    hh = []                       # hhat_n = h_n s^n
    t = ONE
    for n in range(NX+1):
        hh.append(kmul(H[n], t)); t = kmul(t, S)
    G = [ZERO]*(NY+1)
    for N in range(NY+1):
        acc = ZERO
        for n in range(N, min(2*N, NXuse)+1):
            if Pc[n][N]: acc = kadd(acc, kscal(F(Pc[n][N]), hh[n]))
        G[N] = acc
    return H, G

a0, b0, c0 = K(1), ksub(K(-3), kscal(F(-1), PHI)), ksub(K(5), kscal(F(2), PHI))
b0 = kadd(K(-3), PHI); c0 = ksub(K(5), kscal(F(2), PHI))
print("\ngeneric (a,b,c) = (1, -3+phi, 5-2phi) = %s %s %s" % (a0, b0, c0))
H, G = build_G(a0, b0, c0)
print("G head:", [(str(G[n][0]), str(G[n][1])) for n in range(4)])

# --- direct verification: G(Y(x)) = H(x) + H(w(x)) as x-series to order NC
Hc = H[:NC+1]
GY = [ZERO]*(NC+1)
for N in range(NC+1):
    if G[N] != ZERO:
        for i in range(NC+1): GY[i] = kadd(GY[i], kmul(G[N], Yp[N][i]))
Hw = [ZERO]*(NC+1)
for n in range(NC+1):
    if Hc[n] != ZERO:
        for i in range(NC+1): Hw[i] = kadd(Hw[i], kmul(Hc[n], wp[n][i]))
lhs = [kadd(Hc[i], Hw[i]) for i in range(NC+1)]
print("  G(Y(x)) = H(x)+H(w(x)) to order %d: %s [verified]" % (NC, all(lhs[i] == GY[i] for i in range(NC+1))))

# ------------------------------------------------------- denominator type of G
print("\n--- denominator type of G(Y), n <= 60 ---")
tests = [('generic (1,-3+phi,5-2phi)', (a0, b0, c0)), ('(1,0,0)=A only', (ONE, ZERO, ZERO)),
         ('(0,1,0)=B_D only', (ZERO, ONE, ZERO)), ('(0,0,1)=B_new only', (ZERO, ZERO, ONE)),
         ('(phi,phi^2,phi^3)', (PHI, kmul(PHI, PHI), kmul(PHI, kmul(PHI, PHI))))]
CAND = []
for e in (0, 1, 2):
    for lay in ((), (1,), (2,), (1, 1), (1, 2), (2, 2)):
        CAND.append((lay, e))
def tv(cand, n):
    lay, e = cand; v = n**e
    for b in lay: v *= lcm_upto(b*n)
    return v
def nm(cand):
    lay, e = cand
    s = ''.join('[1..n]' if b == 1 else '[1..2n]' for b in lay) or '1'
    return s + (' n' if e == 1 else ' n^2' if e == 2 else '')
NT = 60
for label, abc in tests:
    _, Gx = build_G(*abc)
    dens = [kden(Gx[n]) for n in range(1, NT+1)]
    ok = [c for c in CAND if all(tv(c, n) % dens[n-1] == 0 for n in range(1, NT+1))]
    ok.sort(key=lambda c: tv(c, NT))
    print(f"  {label:28s}: admissible types {[nm(c) for c in ok[:4]]}")
    print(f"      -> [1..2n]^2 valid: {all(lcm_upto(2*n)**2 % dens[n-1] == 0 for n in range(1,NT+1))};"
          f"  [1..n][1..2n] valid: {all(lcm_upto(n)*lcm_upto(2*n) % dens[n-1] == 0 for n in range(1,NT+1))}")
