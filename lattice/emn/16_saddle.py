"""(i) the homogeneous saddle-killer family
      Lsk_n(A,C) = \\iint_Delta x^{An} y^{An} (x^2-y^2)^{2Cn} / (1-x^2-y^2)^{2n+1}."""
import sys, time
sys.path.insert(0,'/home/ubuntu/code/math-modular-sources/lattice/emn')
from emngen import Tmono, linform_gen
from fractions import Fraction as Fr
from math import comb, gcd, lcm
from mpmath import mp, catalan, mpf, log as mlog, mpmathify, nstr
mp.dps = 200
GG = catalan

# ---------- the saddle constant ----------
def Msad(A,C):
    A = mpf(A); C = mpf(C)
    return mpf(4)**(-(A-1)) * (C/(A+C-2))**C * ((A-2)/(A+C-2))**(A-2)
print("[i1] the exponential saddle M(A,C) = 4^{-(A-1)}(C/(A+C-2))^C((A-2)/(A+C-2))^{A-2}")
for (A,C) in [(4,1),(10,3),(6,2),(8,2),(12,4),(14,4),(16,5)]:
    m_ = Msad(A,C)
    print(f"     (A,C)=({A},{C}):  M = {nstr(m_,12)}   1/M = {nstr(1/m_,12)}   -log M = {nstr(-mlog(m_),13)}")
m103 = Msad(10,3)
print(f"     note's value for (10,3): -log M = 18.9221280514...   computed: {nstr(-mlog(m103),16)}")
print(f"     note's margin -log M - 26 log 2 = 0.9003013569...    computed: {nstr(-mlog(m103)-26*mlog(2),16)}")
print(f"     (A,C)=(4,1): M = 1/432 ? -> 1/M = {nstr(1/Msad(4,1),12)}")

# ---------- exact linear forms ----------
def Lsk(n, A, C):
    """returns (a,b) with Lsk_n = a + b*G, exactly."""
    An = A*n; Cn2 = 2*C*n
    tot = (Fr(0),)*5
    from emngen import vadd, vmul
    for j in range(0, C*n):                       # pairs j <-> 2Cn-j
        c = comb(Cn2, j)*(-1)**j
        tot = vadd(tot, vmul(c, Tmono(An+2*(Cn2-j), An+2*j, 2*n)))
    c = comb(Cn2, C*n)*(-1)**(C*n)
    tot = vadd(tot, vmul(c, Tmono(An+Cn2, An+Cn2, 2*n)))
    return linform_gen(tot)

NM = int(sys.argv[1]) if len(sys.argv)>1 else 6
A,C = 10,3
print(f"\n[i2] exact linear forms for (A,C)=({A},{C}):  Lsk_n = ask_n + bsk_n G")
print(f"{'n':>3} {'-log|Lsk_n|/n':>15} {'-log M':>10} {'v_2(den b)':>11} {'/n':>8} "
      f"{'den(b) odd part':>16} {'log den(a)/n':>14}")
data=[]
t0=time.time()
for n in range(1, NM+1):
    a,b = Lsk(n,A,C)
    val = mpf(a.numerator)/mpf(a.denominator) + mpf(b.numerator)/mpf(b.denominator)*GG
    db = b.denominator; v2=0; dd=db
    while dd%2==0: dd//=2; v2+=1
    da = a.denominator
    print(f"{n:>3} {(-float(mlog(abs(val)))/n):>15.6f} {float(-mlog(m103)):>10.6f} {v2:>11} {v2/n:>8.3f} "
          f"{dd:>16} {float(mlog(da))/n:>14.4f}")
    data.append((n,a,b,val))
print(f"     ({time.time()-t0:.1f}s)")

print(f"\n[i3] the PRIMITIVE integer linear form:  clear denominators, then divide by gcd")
print(f"{'n':>3} {'log Dn/n':>11} {'log q_n/n':>11} {'log|l_n|/n':>12} {'gcd bits':>9} "
      f"{'primitive score':>16} {'decay -logM':>12}")
for (n,a,b,val) in data:
    Dn = lcm(a.denominator, b.denominator)
    P = Dn*a.numerator//a.denominator
    Q = Dn*b.numerator//b.denominator
    g = gcd(abs(P), abs(Q))
    Pp, Qp = P//g, Q//g
    l = mpf(Dn)/mpf(g)*val
    lq = float(mlog(abs(Qp)))/n
    ll = float(mlog(abs(l)))/n
    print(f"{n:>3} {float(mlog(Dn))/n:>11.4f} {lq:>11.4f} {ll:>12.4f} {g.bit_length():>9} "
          f"{(-ll):>16.4f} {float(-mlog(m103)):>12.4f}")
print("\n     'primitive score' = -log|primitive linear form|/n;  positive would be needed")
print("     (and then delta = 1 - F/H with F = log|l|/n, H = log q/n)")
for (n,a,b,val) in data:
    Dn = lcm(a.denominator, b.denominator)
    P = Dn*a.numerator//a.denominator; Q = Dn*b.numerator//b.denominator
    g = gcd(abs(P),abs(Q)); Qp=Q//g
    l = mpf(Dn)/mpf(g)*val
    H = float(mlog(abs(Qp)))/n; F = float(mlog(abs(l)))/n
    print(f"     n={n}:  H={H:.4f}  F={F:.4f}  delta = 1 - F/H = {1-F/H:.4f}")
