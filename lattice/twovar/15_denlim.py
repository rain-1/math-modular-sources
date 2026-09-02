"""15: denominator profile and ray limits for each basis solution."""
import sys, pickle, json
from fractions import Fraction as F
from math import gcd
from libx import lcmrange
import lib2v
import mpmath as mp
mp.mp.dps = 80

name = sys.argv[1]
PERIODNAME = sys.argv[2] if len(sys.argv) > 2 else "zeta3"
D = pickle.load(open("sol_%s.pkl" % name, "rb"))
Ts, coef, freecells, NBIG = D["Ts"], D["coef"], D["freecells"], D["NBIG"]
cfun = getattr(lib2v, name)
cT = {(a, b): F(cfun(a, b)) for a in range(NBIG+1) for b in range(NBIG+1)}
PER = dict(zeta3=mp.zeta(3), zeta2=mp.zeta(2), G=mp.catalan)[PERIODNAME]
LCM = [lcmrange(n) for n in range(NBIG+2)]
GA = range(2, 27)


def profile(T):
    valid = {(r1, r2) for r1 in range(9) for r2 in range(9)}
    minr, maxr, wmin, wmax = 0, 0, None, None
    integral = True
    for a in GA:
        for b in GA:
            den = T[(a, b)].denominator
            if den == 1:
                continue
            integral = False
            valid = {(r1, r2) for (r1, r2) in valid
                     if (LCM[a]**r1 * LCM[b]**r2) % den == 0}
            for r in range(0, 12):
                if LCM[min(a, b)]**r % den == 0:
                    break
            else:
                r = 99
            if r > minr:
                minr, wmin = r, (a, b)
            for r in range(0, 12):
                if LCM[max(a, b)]**r % den == 0:
                    break
            else:
                r = 99
            if r > maxr:
                maxr, wmax = r, (a, b)
    mins = sorted([p for p in valid
                   if not any((q[0] <= p[0] and q[1] <= p[1] and q != p) for q in valid)],
                  key=lambda z: (z[0]+z[1], z))
    return mins, minr, maxr, wmin, wmax, integral


def show_profile(T, ind="     "):
    mins, minr, maxr, wmin, wmax, integral = profile(T)
    if integral:
        print(ind + "denominators: INTEGRAL on 2<=a,b<=26 (r1,r2)=(0,0)")
        return mins
    print(ind + "minimal (r1,r2) on 2<=a,b<=26 : %s" % mins)
    print(ind + "den | lcm(1..min(a,b))^r : r=%s (worst %s) ; den | lcm(1..max(a,b))^r : r=%s (worst %s)"
          % (minr if minr < 99 else '>11', wmin, maxr if maxr < 99 else '>11', wmax))
    return mins


def ratio(T, a, b):
    q = T[(a, b)]/cT[(a, b)]
    return mp.mpf(q.numerator)/mp.mpf(q.denominator)


def ident(x):
    if x == 0:
        return "0"
    rel = mp.pslq([mp.mpf(1), x, PER], tol=mp.mpf(10)**(-16), maxcoeff=10**7,
                  maxsteps=10000)
    if rel and rel[1] != 0:
        m0, m1, m2 = rel
        return "%s + (%s)*%s" % (F(-m0, m1), F(-m2, m1), PERIODNAME)
    # maybe rational
    rel = mp.pslq([mp.mpf(1), x], tol=mp.mpf(10)**(-16), maxcoeff=10**7)
    if rel and rel[1] != 0:
        return "%s (rational)" % F(-rel[0], rel[1])
    return "?"


def raylimits(T, ind="     "):
    diagzero = all(T[(n, n)] == 0 for n in range(NBIG+1))
    print(ind + "vanishes identically on the diagonal a=b : %s" % diagzero)
    for lam in (F(1), F(2), F(1, 2), F(3), F(1, 3), F(5), F(1, 5)):
        pts = []
        for b in range(1, NBIG+1):
            a = lam*b
            if a.denominator != 1:
                continue
            a = int(a)
            if 1 <= a <= NBIG:
                pts.append((a, b))
        if not pts:
            continue
        vals = [ratio(T, a, b) for (a, b) in pts]
        last, prev = vals[-1], (vals[-2] if len(vals) > 1 else vals[-1])
        print(ind + "lam=%-4s (a,b)=%-9s  t/c = %s   (prev pt %s)   ->  %s"
              % (lam, str(pts[-1]), mp.nstr(last, 22), mp.nstr(prev, 10), ident(last)))


print("=" * 78)
print("CASE", name, "  period", PERIODNAME, "  NBIG", NBIG)
print("   free cells", freecells, "  row c coefficients", coef)
for k, T in enumerate(Ts):
    print("-" * 70)
    print("  basis %d  free vals %s ; t(1,1)=%s t(1,2)=%s t(2,1)=%s t(2,2)=%s t(3,2)=%s"
          % (k, [T[fc] for fc in freecells], T[(1, 1)], T[(1, 2)], T[(2, 1)],
             T[(2, 2)], T[(3, 2)]))
    show_profile(T)
    raylimits(T)
    sys.stdout.flush()
