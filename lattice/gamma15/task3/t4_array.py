"""Measured denominator array for ALL FOURTEEN functions on the Gamma_1(5) host,
then tau^flat, tau^sharp, tau via lattice/cdt_finder/cdt_bound.py."""
from fractions import Fraction as F
from math import factorial as fac, gcd
import sys
sys.path.insert(0, '/home/ubuntu/code/math-modular-sources/lattice/gamma15/task3')
sys.path.insert(0, '/home/ubuntu/code/math-modular-sources/lattice/cdt_finder')
from kfield import *
import cdt_bound

NY = 70; NXuse = 2*NY
def lcm_upto(m, c={}):
    if m in c: return c[m]
    r = 1
    for j in range(1, m+1): r = r*j//gcd(r, j)
    c[m] = r; return r
Pc = [[0]*(NY+2) for _ in range(NXuse+1)]
Pc[0][0] = 2; Pc[1][1] = 1
for n in range(2, NXuse+1):
    for j in range(1, NY+2): Pc[n][j] = Pc[n-1][j-1]-Pc[n-2][j-1]
def load(fn):
    t = open(fn).read().strip(); t = t[t.index('[')+1:t.rindex(']')]
    return [F(u.strip()) for u in t.split(',')]
A = load('A_coeffs.txt'); BD = load('BD_coeffs.txt')
BR3 = load('BR3_coeffs.txt'); BR4 = load('BR4_coeffs.txt')
a = K(1); b = kadd(K(-3), PHI); c = ksub(K(5), kscal(F(2), PHI))
H = [kadd(kadd(kscal(A[n], a), kscal(BD[n], b)),
          kmul(c, kadd(K(BR3[n]), kmul(PHI5, K(BR4[n]))))) for n in range(NXuse+1)]
hh = []; t = ONE
for n in range(NXuse+1): hh.append(kmul(H[n], t)); t = kmul(t, S)
G = [ZERO]*(NY+1)
for N in range(NY+1):
    acc = ZERO
    for n in range(N, min(2*N, NXuse)+1):
        if Pc[n][N]: acc = kadd(acc, kscal(F(Pc[n][N]), hh[n]))
    G[N] = acc
def dY(f): return [kscal(F(n+1), f[n+1]) for n in range(len(f)-1)]
def iY(f): return [ZERO]+[kscal(F(1, n+1), f[n]) for n in range(len(f)-1)]
B1 = [F(1)]+[F(0)]*NY
B2 = [F(0)]*(NY+1); B3 = [F(0)]*(NY+1); B5 = [F(0)]*(NY+1); B6 = [F(0)]*(NY+1)
for n in range(1, NY+1):
    if n >= 2: B2[n] = F(2*fac(n-2)*fac(n), fac(2*n))
    B3[n] = F(fac(n-1)**2, fac(2*n)); B5[n] = F(fac(n-1)**2, fac(2*n-1)*(2*n-1))
    B6[n] = F(fac(n-1)**2, n*fac(2*n))
Cb = [F(fac(2*k), fac(k)**2) for k in range(NY+2)]
B4 = [F(0)]*(NY+1)
for n in range(NY): B4[n+1] = 4*sum(Cb[k]*Cb[n-k]*F(1, (2*k-1)*(2*n-2*k+1)**2) for k in range(n+1))/F(16)**n
B7 = [F(0)]+[B4[n]/n for n in range(1, NY+1)]
KB = lambda B: [K(u) for u in B]
fns = [('B1', KB(B1)), ('B2', KB(B2)), ('B3', KB(B3)), ('B4', KB(B4)), ('B5', KB(B5)),
       ('G', G), ("G'", dY(G)), ("G''", dY(dY(G))), ("G'''", dY(dY(dY(G)))),
       ('B6', KB(B6)), ('B7', KB(B7)), ('intG', iY(G)),
       ('int(G-G0)/y', [ZERO]+[kscal(F(1, n), G[n]) for n in range(1, len(G))]),
       ('int(G-G0-G1y)/y^2', [ZERO]+[kscal(F(1, n), G[n+1]) for n in range(1, len(G)-1)])]

# family allowed by CDT Remark rem:overflow:  n^e prod_j [1 .. b_j n + c_j]
CAND = []
for e in (0, 1, 2, 3):
    for lay in ((), ((1,0),), ((2,0),), ((2,1),), ((2,2),), ((2,3),), ((2,4),),
                ((1,0),(1,0)), ((1,0),(2,0)), ((2,0),(2,0)), ((2,0),(2,1)),
                ((2,0),(2,2)), ((2,1),(2,1)), ((2,2),(2,2)), ((2,0),(2,3)),
                ((2,1),(2,2)), ((2,1),(2,3)), ((2,2),(2,3)), ((2,3),(2,3)),
                ((2,2),(2,4)), ((2,4),(2,4)), ((2,5),(2,5)), ((2,6),(2,6)), ((2,4),(2,5)), ((2,5),(2,6))):
        CAND.append((lay, e))
def tv(cd, n):
    lay, e = cd; v = n**e
    for bb, cc in lay: v *= lcm_upto(bb*n+cc)
    return v
def nmz(cd):
    lay, e = cd
    s = ''.join('[1..%s%s]' % ('n' if bb == 1 else '%dn' % bb, '+%d' % cc if cc else '') for bb, cc in lay) or '1'
    return s+(' n^%d' % e if e > 1 else ' n' if e == 1 else '')
NT = 50
print("=== MEASURED sharp denominator types of the fourteen (Gamma_1(5) host, n <= %d) ===" % NT)
print("   (types are for the coefficient of Y^n; Y = y/s, s a unit, so identical in y)")
rows = []
for nm, f in fns:
    dens = [kden(f[n]) for n in range(1, NT+1)]
    ok = [cd for cd in CAND if all(tv(cd, n) % dens[n-1] == 0 for n in range(1, NT+1))]
    ok.sort(key=lambda cd: (len(cd[0]), cd[1], tv(cd, NT)))
    lay22 = all(lcm_upto(2*n+6)**2 % dens[n-1] == 0 for n in range(1, NT+1))
    lay22n = all((lcm_upto(2*n+6)**2*n) % dens[n-1] == 0 for n in range(1, NT+1))
    lay2 = all(lcm_upto(2*n+6) % dens[n-1] == 0 for n in range(1, NT+1))
    lay2n = all((lcm_upto(2*n+6)*n) % dens[n-1] == 0 for n in range(1, NT+1))
    triv = all(dens[n-1] == 1 for n in range(1, NT+1))
    rows.append((nm, ok, triv, lay2, lay2n, lay22, lay22n))
    print(f"  {nm:20s} minimal admissible: {[nmz(cd) for cd in ok[:3]]}")
print()
print("  summary flags (overflow-relaxed, c=6):  trivial | [1..2n+6] | [1..2n+6]n | [1..2n+6]^2 | [1..2n+6]^2 n")
for nm, ok, triv, l2, l2n, l22, l22n in rows:
    print(f"  {nm:20s}  {str(triv):7s} {str(l2):8s} {str(l2n):9s} {str(l22):10s} {str(l22n)}")

print()
print("="*72)
print("T4: tau from cdt_bound.py")
import math
def show(tag, m, cols, e):
    sm, tf = cdt_bound.tau_flat(m, cols)
    ts, xi = cdt_bound.tau_sharp(m, e, ngrid=200001)
    print(f"  {tag}")
    print(f"    cols(u_j,b_j) = {cols},  e = {e}  (sum e = {sum(e)}, max e = {max(e)})")
    print(f"    sigma_m = {sm};  tau^flat = {tf} = {float(tf):.9f}")
    print(f"    tau^sharp = {ts:.10f}  at xi = {xi:.5f}")
    print(f"    tau = {float(tf)+ts:.9f}")
    return float(tf)+ts
mA = show("CDT array (u_1=1,u_2=3; e sum 6, max 1)", 14, [(1, 2), (3, 2)],
          [0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
print(f"    CDT's exact value 16603/3920 = {16603/3920:.9f};  27/80 = {27/80};  191/49 = {191/49:.9f}")
mB = show("ALT: B6 as [1..2n]n^2 (u_1=1,u_2=4; e sum 8, max 2)", 14, [(1, 2), (4, 2)],
          [0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1])
print(f"  --> CDT's own choice is better by {mB-mA:.6f}" if mB > mA else "  --> ALT better")
