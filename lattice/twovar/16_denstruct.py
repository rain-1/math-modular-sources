"""16: structure of the denominators - which primes occur, and are they inside
lcm(1..a), lcm(1..b)?  Plus the support of diagonal-vanishing solutions."""
import sys, pickle
from fractions import Fraction as F
from sympy import factorint
from libx import lcmrange
import lib2v

name = sys.argv[1]
D = pickle.load(open("sol_%s.pkl" % name, "rb"))
Ts, freecells, NBIG = D["Ts"], D["freecells"], D["NBIG"]
LCM = [lcmrange(n) for n in range(NBIG+2)]


def vp(n, p):
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


print("=" * 78)
print("CASE", name)
for k, T in enumerate(Ts):
    print("-" * 70)
    print(" basis %d  free vals %s" % (k, [T[fc] for fc in freecells]))
    nz = [(a, b) for a in range(9) for b in range(9) if T[(a, b)] != 0]
    print("   nonzero cells with a,b<=8 :", nz[:30])
    print("   values 0<=a,b<=5:")
    for a in range(6):
        print("     a=%d : %s" % (a, [str(T[(a, b)]) for b in range(6)]))
    print("   denominator anatomy (p^e ; e_max allowed by lcm(1..a),lcm(1..b)):")
    for (a, b) in [(6, 6), (10, 10), (14, 14), (20, 20), (26, 26),
                   (4, 12), (12, 4), (6, 20), (20, 6), (10, 26), (26, 10),
                   (2, 6), (2, 10), (3, 14)]:
        if a > NBIG or b > NBIG:
            continue
        den = T[(a, b)].denominator
        if den == 1:
            print("     (%2d,%2d) den = 1" % (a, b)); continue
        parts = []
        extra = []
        for p, e in sorted(factorint(den).items()):
            ea, eb = vp(LCM[a], p), vp(LCM[b], p)
            parts.append("%d^%d(a:%d,b:%d)" % (p, e, ea, eb))
            if ea == 0 and eb == 0:
                extra.append((p, e))
        print("     (%2d,%2d) den = %s   OUTSIDE lcm: %s" % (a, b, " ".join(parts), extra))
    sys.stdout.flush()
