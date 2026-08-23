"""(h) radial Pade family  Lpad_n = \\iint x^{4n} y^{4n} (1-4g^2)^{2n} / g^{2n+1},  g = 1-x^2-y^2.
   (1-4g^2)^{2n} = sum_k C(2n,k)(-4)^k g^{2k}  =>  Lpad_n = sum_k C(2n,k)(-4)^k I_{4n, 2n-2k}."""
import sys, time
sys.path.insert(0,'/home/ubuntu/code/math-modular-sources/lattice/emn')
from emnexact import linform
from fractions import Fraction as Fr
from math import comb, log
from mpmath import mp, catalan, mpf, log as mlog
from sympy import factorint
mp.dps = 300
NMAX = int(sys.argv[1]) if len(sys.argv)>1 else 16

# ---- closed form for the G-coefficient (no integration needed) ----
def Bpad_int(n):
    """256^n * (G-coefficient of Lpad_n), as an integer, from b_{m,t}=(-1)^t 4^{-m}C(m,m/2)C(m,t)."""
    tot = 0
    for k in range(0, 2*n+1):
        t = 2*n-2*k
        if t < 0: break
        tot += comb(2*n,k)*(-4)**k*comb(4*n,t)
    return comb(4*n,2*n)*tot
print("[h1] closed form:  256^n b_n^{Pade} = C(4n,2n) * sum_k C(2n,k)(-4)^k C(4n,2n-2k)   -- an INTEGER by construction")
print("     n:   256^n b_n  (integer)                                  log/n")
for n in list(range(0,9))+[16,24,32,40]:
    v = Bpad_int(n)
    print(f"     {n:>3}: {str(v)[:46]:<46} {'' if n==0 else f'{float(mlog(abs(v)))/n:>10.6f}'}   {'INTEGER OK' if isinstance(v,int) else ''}")

# ---- exact linear forms ----
t0=time.time()
AP=[]; BP=[]
for n in range(NMAX+1):
    a=Fr(0); b=Fr(0)
    for k in range(0, 2*n+1):
        t = 2*n-2*k
        c = Fr(comb(2*n,k)*(-4)**k)
        aa,bb = linform(4*n, t)
        a += c*aa; b += c*bb
    AP.append(a); BP.append(b)
    assert b == Fr(Bpad_int(n), 256**n), f"G-coefficient mismatch at n={n}"
print(f"\n[h2] exact integral computation agrees with the closed form for n<={NMAX}  ({time.time()-t0:.1f}s)")

GG = catalan
print(f"\n[h3] the bound  0 < 256^n Lpad_n <= G (16/27)^n   [log(27/16) = {float(mlog(27/16)):.7f}]")
print(f"{'n':>4} {'256^n L_n':>26} {'G(16/27)^n':>26} {'ok':>4} {'-log(256^n L_n)/n':>19}")
for n in range(NMAX+1):
    val = mpf(AP[n].numerator)/mpf(AP[n].denominator) + mpf(BP[n].numerator)/mpf(BP[n].denominator)*GG
    v = mpf(256)**n*val
    bnd = GG*(mpf(16)/27)**n
    ok = "YES" if 0 < v <= bnd else "NO"
    r = '' if n==0 else f"{-float(mlog(abs(v)))/n:>19.6f}"
    print(f"{n:>4} {mp.nstr(v,12):>26} {mp.nstr(bnd,12):>26} {ok:>4} {r}")

print(f"\n[h4] denominator cost of the rational companion a_n^{{Pade}}")
print(f"{'n':>4} {'log den(a_n)/n':>16} {'v_2(den)/n':>12} {'log(odd part)/n':>17} {'max prime':>10} {'4n':>6}")
for n in range(2,NMAX+1,2):
    d = AP[n].denominator
    v2=0; dd=d
    while dd%2==0: dd//=2; v2+=1
    fac = factorint(dd)
    print(f"{n:>4} {float(mlog(d))/n:>16.6f} {v2/n:>12.4f} {float(mlog(dd))/n:>17.6f} {(max(fac) if fac else 1):>10} {4*n:>6}")

print(f"\n[h5] the cleared Z-linear form:  Dn = lcm(den a_n, 256^n);  Dn*L_n = (Dn a_n) + (Dn b_n) G")
print(f"{'n':>4} {'log Dn/n':>12} {'log|Dn L_n|/n':>15} {'log|Dn b_n|/n':>15} {'delta = -log|l|/log|q|':>24}")
for n in range(2,NMAX+1,2):
    from math import lcm
    Dn = lcm(AP[n].denominator, 256**n)
    val = mpf(AP[n].numerator)/mpf(AP[n].denominator) + mpf(BP[n].numerator)/mpf(BP[n].denominator)*GG
    q = Fr(Dn)*BP[n]; l = mpf(Dn)*val
    ld = float(mlog(Dn))/n
    ll = float(mlog(abs(l)))/n
    lq = float(mlog(abs(mpf(q.numerator)/mpf(q.denominator))))/n
    print(f"{n:>4} {ld:>12.6f} {ll:>15.6f} {lq:>15.6f} {(-ll/lq):>24.6f}")
