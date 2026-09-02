"""20: the closed form of the two-variable companion,
    s_{a,b} = d^(R)_{a,b} + (H^(3)_b - 1) c_{a,b},
which PROVES the max-type denominator bound (d^(R) has [1..a]^3, H^(3)_b has
[1..b]^3, c is integral), and shows it is sharp at primes in (min,max]."""
from lib2v import *
from fractions import Fraction as F
from math import comb
NA = NB = 30


def build(t00, t10, t01, t11, NA=NA, NB=NB):
    t = [[F(0)]*(NB+1) for _ in range(NA+1)]
    t[0][0], t[1][0], t[0][1], t[1][1] = map(F, (t00, t10, t01, t11))
    for b in (0, 1):
        B = b*(b+1)
        for n in range(1, NA):
            mid = 2*n**3+3*n**2+(3+4*B)*n+1+2*B
            t[n+1][b] = (F(mid)*t[n][b]-F(n**3)*t[n-1][b])/F((n+1)**3)
    for a in range(NA+1):
        for b in range(NB-1):
            Q1 = -(2*b**3+9*b**2+(15+4*a+4*a*a)*b+9+6*a+6*a*a)
            t[a][b+2] = (-F(Q1)*t[a][b+1]-F((b+1)**3)*t[a][b])/F((b+2)**3)
    return t


c = build(1, 1, 1, 5)
s = build(-1, 0, 0, 1)
dR = [[F(0)]*(NB+1) for _ in range(NA+1)]
for b in range(NB+1):
    B = b*(b+1)
    col = [F(0), F(1)]
    for n in range(1, NA):
        mid = 2*n**3+3*n**2+(3+4*B)*n+1+2*B
        col.append((F(mid)*col[n]-F(n**3)*col[n-1])/F((n+1)**3))
    for a in range(NA+1):
        dR[a][b] = col[a]


def H3(b):
    return sum(F(1, m**3) for m in range(1, b+1))


ok = all(s[a][b] == dR[a][b] + (H3(b)-1)*c[a][b]
         for a in range(NA+1) for b in range(NB+1))
print("  s_{a,b} = d^(R)_{a,b} + (H^(3)_b - 1) c_{a,b}  on 0<=a,b<=%d :" % NA, ok)
print("  den(d^(R)_{a,b}) | [1..a]^3 :",
      all(lcmrange(a)**3 % dR[a][b].denominator == 0
          for a in range(1, NA+1) for b in range(NB+1)))
print("  => den(s_{a,b}) | lcm([1..a]^3,[1..b]^3) = [1..max(a,b)]^3   (PROVED)")
print()
print("  sharpness: v_p(den s_{a,b}) for a<p<=b (should be 3):")
from sympy import primerange
bad = 0
tot = 0
for a in range(2, NA+1):
    for b in range(a+1, NB+1):
        for p in primerange(a+1, b+1):
            if p*p <= b:
                continue
            tot += 1
            d = s[a][b].denominator
            v = 0
            while d % p == 0:
                d //= p
                v += 1
            if v != 3:
                bad += 1
print("     %d of %d cells (a<p<=b, p^2>b) have v_p(den) != 3" % (bad, tot))
