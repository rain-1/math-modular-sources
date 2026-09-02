"""17: (A) denominator type of a solution up to a BOUNDED global factor;
       (B) the 'extra' primes p>max(a,b) and the subspace killing them;
       (C) scan of small-height combinations."""
import sys, pickle, itertools
from fractions import Fraction as F
from math import gcd
from libx import lcmrange
import lib2v
import mpmath as mp
mp.mp.dps = 60

name = sys.argv[1]
PERIODNAME = sys.argv[2] if len(sys.argv) > 2 else "zeta3"
NG = int(sys.argv[3]) if len(sys.argv) > 3 else 22
Dd = pickle.load(open("sol_%s.pkl" % name, "rb"))
Ts, coef, freecells, NBIG = Dd["Ts"], Dd["coef"], Dd["freecells"], Dd["NBIG"]
DIM = len(Ts)
cfun = getattr(lib2v, name)
PER = dict(zeta3=mp.zeta(3), zeta2=mp.zeta(2), G=mp.catalan)[PERIODNAME]
LCM = [lcmrange(n) for n in range(NBIG+2)]
CELLS = [(a, b) for a in range(2, NG+1) for b in range(2, NG+1)]
BV = {ab: tuple(T[ab] for T in Ts) for ab in CELLS}
BIG = {n: LCM[n]**12 for n in range(NG+1)}
PRIMES = [q for q in range(2, 400) if all(q % r for r in range(2, int(q**.5)+1))]


def val(x, ab):
    s = F(0)
    for xi, v in zip(x, BV[ab]):
        if xi:
            s += xi*v
    return s


def extra_part(den, a, b):
    m = max(a, b)
    return den // gcd(den, BIG[m])


def has_extra(x, N=None):
    N = N or NG
    for (a, b) in CELLS:
        if a > N or b > N:
            continue
        den = val(x, (a, b)).denominator
        if den > 1 and extra_part(den, a, b) > 1:
            return (a, b, extra_part(den, a, b))
    return None


def classify(x, Ns=(10, 16, NG), rmax=4):
    pairs = [(r1, r2) for r1 in range(rmax+1) for r2 in range(rmax+1)]
    acc = {(pr, N): 1 for pr in pairs for N in Ns}
    for (a, b) in CELLS:
        den = val(x, (a, b)).denominator
        if den == 1:
            continue
        for (r1, r2) in pairs:
            X = den // gcd(den, LCM[a]**r1 * LCM[b]**r2)
            if X > 1:
                for N in Ns:
                    if a <= N and b <= N:
                        L = acc[((r1, r2), N)]
                        acc[((r1, r2), N)] = L*X//gcd(L, X)
    good = [(pr, acc[(pr, Ns[-1])]) for pr in pairs
            if acc[(pr, Ns[0])] == acc[(pr, Ns[-1])] == acc[(pr, Ns[1])]]
    gset = {p for p, _ in good}
    mins = sorted([p for p in gset
                   if not any(q[0] <= p[0] and q[1] <= p[1] and q != p for q in gset)],
                  key=lambda z: (z[0]+z[1], z))
    return mins, dict(good)


def ident(v):
    if v == 0:
        return "0"
    rel = mp.pslq([mp.mpf(1), v, PER], tol=mp.mpf(10)**(-16), maxcoeff=10**7,
                  maxsteps=20000)
    if rel and rel[1] != 0:
        return "%s + (%s)*%s" % (F(-rel[0], rel[1]), F(-rel[2], rel[1]), PERIODNAME)
    return "?"


cNN = F(cfun(NBIG, NBIG))
def limit1(x):
    q = sum(F(xi)*T[(NBIG, NBIG)] for xi, T in zip(x, Ts))/cNN
    return mp.mpf(q.numerator)/mp.mpf(q.denominator)


print("=" * 78)
print("CASE %s   dim=%d   row coefficients x_c = %s   grid 2..%d" % (name, DIM, coef, NG))
xc = [int(z) for z in coef]

print("\n(A) per basis solution")
for k in range(DIM):
    x = [1 if i == k else 0 for i in range(DIM)]
    mins, good = classify(x)
    ex = has_extra(x)
    lim = limit1(x)
    print("   B_%d : minimal (r1,r2)=%-14s extra prime? %s" % (k, str(mins), ex))
    print("         limit B_%d/c (a=b=%d) = %s = %s" % (k, NBIG, mp.nstr(lim, 20), ident(lim)))
print("   row  : minimal (r1,r2)=%s  (integral)" % str(classify(xc)[0]))

print("\n(B) extra primes: which x kill them?")
xs = []
for x in itertools.product(range(-4, 5), repeat=DIM):
    if all(v == 0 for v in x):
        continue
    g = 0
    for v in x:
        g = gcd(g, abs(v))
    xr = tuple(v//g for v in x)
    if xr[next(i for i, v in enumerate(xr) if v)] < 0:
        continue
    if xr not in xs:
        xs.append(xr)
def lcmtype(x, Ns=(10, 16, NG)):
    """cheap test: is the excess over lcm(1..a)^4 lcm(1..b)^4 BOUNDED?"""
    acc = {N: 1 for N in Ns}
    for (a, b) in CELLS:
        den = val(x, (a, b)).denominator
        if den == 1:
            continue
        X = den // gcd(den, LCM[a]**4 * LCM[b]**4)
        if X > 1:
            for N in Ns:
                if a <= N and b <= N:
                    acc[N] = acc[N]*X//gcd(acc[N], X)
    return acc[Ns[0]] == acc[Ns[1]] == acc[Ns[-1]], acc[Ns[-1]]

clean = []
for x in xs:
    ok, e = lcmtype(list(x))
    if ok:
        clean.append((x, e))
print("   %d of %d primitive small-height directions are of lcm type"
      " (excess over lcm^4 lcm^4 BOUNDED on nested grids 2..10/16/%d)"
      % (len(clean), len(xs), NG))

print("\n(C) SUMMARY over the lcm-type directions: (r1,r2) vs the period coefficient")
from collections import defaultdict
summary = defaultdict(list)
for (x, e) in clean:
    mins, good = classify(list(x))
    lim = limit1(list(x))
    rel = mp.pslq([mp.mpf(1), lim, PER], tol=mp.mpf(10)**(-16), maxcoeff=10**7,
                  maxsteps=20000) if lim != 0 else [0, 1, 0]
    if rel and rel[1] != 0:
        beta = F(-rel[2], rel[1])
        alpha = F(-rel[0], rel[1])
    else:
        beta, alpha = None, None
    summary[(tuple(mins), beta != 0 if beta is not None else "?")].append((x, alpha, beta, e))
for key in sorted(summary, key=lambda z: (str(z))):
    v = summary[key]
    print("   (r1,r2) in %-14s  period coeff nonzero: %-5s   %d directions, e.g. %s"
          % (str(list(key[0])), str(key[1]), len(v),
             ["x=%s -> %s + %s*P" % (str(a), str(b), str(c)) for (a, b, c, _) in v[:3]]))
print()
print("   BEST (r1,r2) among directions WITH a nonzero period coefficient:")
bestv = [ (sum(m), m, x, al, be) for key, v in summary.items() if key[1] is True
          for m in key[0] for (x, al, be, _) in v ]
if bestv:
    bestv.sort()
    s0 = bestv[0][0]
    n = 0
    for (s_, m, x, al, be) in bestv:
        if s_ > s0 or n >= 4:
            break
        n += 1
        print("      (r1,r2)=%s  x=%s  limit = %s + (%s)*%s" % (str(m), str(x), al, be, PERIODNAME))
    print("      (%d directions attain r1+r2=%d)" % (sum(1 for t in bestv if t[0] == s0), s0))
else:
    print("      NONE - no lcm-type direction has a nonzero period coefficient")
print()
print("   BEST (r1,r2) among ALL non-row lcm-type directions:")
allv = sorted([(sum(m), m, x, key[1]) for key, v in summary.items()
               for m in key[0] for (x, _, _, _) in v if tuple(int(z) for z in xc) != tuple(x)])
if allv:
    s0 = allv[0][0]
    n = 0
    for (s_, m, x, nz) in allv:
        if s_ > s0 or n >= 4:
            break
        n += 1
        print("      (r1,r2)=%s  x=%s  period-coeff-nonzero=%s" % (str(m), str(x), nz))
    print("      (%d directions attain r1+r2=%d)" % (sum(1 for t in allv if t[0] == s0), s0))
