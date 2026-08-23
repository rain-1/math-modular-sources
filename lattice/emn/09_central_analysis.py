"""(g) analysis of the central family: decay, denominators, holonomic fit."""
import sys, time
sys.path.insert(0,'/home/ubuntu/code/math-modular-sources/lattice/emn')
from emnexact import linform
from fractions import Fraction as Fr
from math import comb, log, gcd
from mpmath import mp, catalan, log as mlog, mpf, nstr
mp.dps = 400
NMAX = int(sys.argv[1]) if len(sys.argv)>1 else 80

def central(n):
    if n == 0: return linform(0,0)
    a, b = linform(2*n, 0)
    for k in range(1, n+1):
        m, t = 2*n+2*k, 2*k
        f = Fr((m-1)**2, 4*t*(t-1))
        a = Fr(-1, 2**(t+1)*t*(t-1)*comb(2*(m-t), m-t)) + f*a
        b = f*b
    return a, b

t0=time.time()
A=[]; B=[]
for n in range(NMAX+1):
    a,b = central(n); A.append(a); B.append(b)
print(f"computed n<= {NMAX} in {time.time()-t0:.1f}s")

# ---------- 1. B_n closed form ----------
bad = sum(1 for n in range(NMAX+1) if B[n] != Fr(comb(4*n,2*n)**2, 256**n))
print(f"\n[g1] B_n = C(4n,2n)^2/256^n : mismatches for n<={NMAX}: {bad}   =>  256^n B_n in Z")

# ---------- 2. decay of I_n ----------
GG = catalan
print(f"\n[g2] decay of I_n = A_n + B_n G   (claim |I_n| = 64^{{-n+o(n)}}, log 64 = {float(mlog(64)):.6f})")
print(f"{'n':>4} {'-log|I_n|/n':>14} {'log|A_n|/n':>13} {'-log|den A_n|/n':>16} {'I_n sign':>9}")
lastv=None
for n in [1,2,4,8,16,32,48,64,NMAX]:
    if n>NMAX: continue
    val = mpf(A[n].numerator)/mpf(A[n].denominator) + mpf(B[n].numerator)/mpf(B[n].denominator)*GG
    ln = -float(mlog(abs(val)))/n
    la = float(mlog(abs(mpf(A[n].numerator)/mpf(A[n].denominator))))/n
    ld = -float(mlog(mpf(A[n].denominator)))/n
    print(f"{n:>4} {ln:>14.6f} {la:>13.6f} {ld:>16.6f} {'+' if val>0 else '-':>9}")

# ---------- 3. denominators of A_n ----------
print(f"\n[g3] denominator structure of A_n")
def vp(x,p):
    v=0
    while x%p==0: x//=p; v+=1
    return v
from sympy import factorint, primerange
print(f"{'n':>4} {'log den(A_n)/n':>16} {'v_2':>7} {'largest prime':>14} {'odd part log/n':>16}")
for n in [4,8,16,24,32,48,64,NMAX]:
    if n>NMAX: continue
    d = A[n].denominator
    v2 = vp(d,2) if d%2==0 else 0
    odd = d>>v2
    fac = factorint(odd)
    lp = max(fac) if fac else 1
    print(f"{n:>4} {float(mlog(d))/n:>16.6f} {v2:>7} {lp:>14} {float(mlog(odd))/n:>16.6f}")

# lcm comparison
def dlcm(N):
    from math import lcm
    r=1
    for k in range(1,N+1): r=r*k//gcd(r,k)
    return r
print(f"\n     comparison with lcm-type denominators at n={NMAX}:")
d = A[NMAX].denominator
for lab, val in [('d_{2n}', dlcm(2*NMAX)), ('d_{4n}', dlcm(4*NMAX)),
                 ('d_{2n}^2', dlcm(2*NMAX)**2), ('d_{4n}^2', dlcm(4*NMAX)**2)]:
    print(f"     {lab:>10}: log/n = {float(mlog(val))/NMAX:>9.6f}   den(A_n) | {lab}?  {'YES' if val % d == 0 else 'no'}")
print(f"     {'den(A_n)':>10}: log/n = {float(mlog(d))/NMAX:>9.6f}")

# exact per-prime valuations, several n
print(f"\n     v_p(den A_n) profile (p, v_p) at n={NMAX}:  (also floor(log_p(4n)), floor(log_p(2n)))")
d = A[NMAX].denominator
fac = factorint(d)
import math
rows=[]
for p in sorted(fac):
    rows.append((p, fac[p], int(math.log(4*NMAX)/math.log(p)), int(math.log(2*NMAX)/math.log(p))))
print("     ", rows[:14], "..." if len(rows)>14 else "")
print(f"     number of distinct primes: {len(fac)};  max prime {max(fac)} vs 4n={4*NMAX}, 2n={2*NMAX}")
