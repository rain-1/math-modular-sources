"""Measure denominators and p-adic slopes for the EMN generator, the cubic
conditional K, the fold companion pair (A,B) and rescaled variants."""
from fractions import Fraction as Fr
import math, sys
from emn_core import (series_H, solve_from_psi, ratseries, polymul, chebT_in_z,
                      psi_basis, den_profile, padic_slope)

N = int(sys.argv[1]) if len(sys.argv) > 1 else 1200
NS = [n for n in (100,200,400,600,800,1000,1200,1600,2000,2500,3000) if n <= N]

def vp(fr, p):
    num, den = fr.numerator, fr.denominator
    if num == 0: return None
    v = 0
    while num % p == 0: num//=p; v+=1
    while den % p == 0: den//=p; v-=1
    return v

def report(name, c):
    dp = den_profile(c, NS)
    line = f"{name:34s} log den/n: " + " ".join(f"{n}:{dp[n]:.4f}" for n in NS)
    print(line)
    for p in (2,3):
        vs = [vp(c[n], p) for n in NS]
        sl = [(v/n if v is not None else None) for v,n in zip(vs,NS)]
        print(f"     v_{p}: " + " ".join(f"{n}:{v}" for n,v in zip(NS,vs)) +
              "   slope~" + (f"{sl[-1]:.4f}" if sl[-1] is not None else "-"))

H = series_H(N)
report("H(z)  (EMN generator)", H)

# --- cubic conditional K = script-H(3 ang) + 3 script-H(ang) - 4G
T3 = chebT_in_z(3)
psiK = [x+y for x,y in zip(ratseries([Fr(9,2)], T3, N), [Fr(3,2)]*(N+1))]
K = solve_from_psi(psiK, N)
report("K = H(R_3)+3H(z)  [cubic cond]", K)

# --- fold companions at z=2 in the coordinate Z = z/2 (fold Z=1, outer Z=1/2)
Ht = [H[n]*Fr(2)**n for n in range(N+1)]          # H(2Z)
report("Htil(Z)=H(2Z)", Ht)
# B(Z) = int_0^Z Htil(T) dT/(1-T)
part = [Fr(0)]*(N+1)
run = Fr(0)
for n in range(N+1):
    part[n] = run          # sum_{k<n} Htil_k
    run += Ht[n]
B = [Fr(0)]+[part[n]/Fr(n) for n in range(1,N+1)]
report("B(Z)=int H(2T)dT/(1-T)  [cond]", B)
A = [Fr(0)]+[Fr(1,n) for n in range(1,N+1)]       # -log(1-Z)
report("A(Z)=-log(1-Z)", A)
