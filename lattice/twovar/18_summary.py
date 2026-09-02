"""18: compact per-case summary: the row, the 'small' (limit 0) subspace, and a
companion direction, each with its denominator type + bounded excess constant."""
import sys, pickle, itertools
from fractions import Fraction as F
from math import gcd
from libx import lcmrange
import lib2v
import mpmath as mp
mp.mp.dps = 60

name, PERIODNAME = sys.argv[1], sys.argv[2]
NG = 22
Dd = pickle.load(open("sol_%s.pkl" % name, "rb"))
Ts, coef, NBIG = Dd["Ts"], Dd["coef"], Dd["NBIG"]
DIM = len(Ts)
cfun = getattr(lib2v, name)
PER = dict(zeta3=mp.zeta(3), zeta2=mp.zeta(2), G=mp.catalan)[PERIODNAME]
LCM = [lcmrange(n) for n in range(NBIG+2)]
CELLS = [(a, b) for a in range(2, NG+1) for b in range(2, NG+1)]
BV = {ab: tuple(T[ab] for T in Ts) for ab in CELLS}
cNN = F(cfun(NBIG, NBIG))


def val(x, ab):
    s = F(0)
    for xi, v in zip(x, BV[ab]):
        if xi:
            s += F(xi)*v
    return s


def classify(x, Ns=(10, 16, NG), rmax=4):
    pairs = [(r1, r2) for r1 in range(rmax+1) for r2 in range(rmax+1)]
    acc = {(pr, N): 1 for pr in pairs for N in Ns}
    for (a, b) in CELLS:
        den = val(x, (a, b)).denominator
        if den == 1:
            continue
        for pr in pairs:
            X = den // gcd(den, LCM[a]**pr[0] * LCM[b]**pr[1])
            if X > 1:
                for N in Ns:
                    if a <= N and b <= N:
                        L = acc[(pr, N)]
                        acc[(pr, N)] = L*X//gcd(L, X)
    good = {pr: acc[(pr, Ns[-1])] for pr in pairs
            if acc[(pr, Ns[0])] == acc[(pr, Ns[1])] == acc[(pr, Ns[-1])]}
    mins = sorted([p for p in good
                   if not any(q[0] <= p[0] and q[1] <= p[1] and q != p for q in good)],
                  key=lambda z: (z[0]+z[1], z))
    return mins, {m: good[m] for m in mins}


def lim(x):
    q = sum(F(xi)*T[(NBIG, NBIG)] for xi, T in zip(x, Ts))/cNN
    return mp.mpf(q.numerator)/mp.mpf(q.denominator)


def ident(v):
    if v == 0:
        return F(0), F(0)
    rel = mp.pslq([mp.mpf(1), v, PER], tol=mp.mpf(10)**(-15), maxcoeff=10**6,
                  maxsteps=20000)
    if rel and rel[1] != 0:
        return F(-rel[0], rel[1]), F(-rel[2], rel[1])
    return None, None


# beta functional
betas = []
for k in range(DIM):
    e = [1 if i == k else 0 for i in range(DIM)]
    al, be = ident(lim(e))
    betas.append((al, be))
print("=" * 78)
print("CASE %s   dim=%d   x_c=%s" % (name, DIM, [str(c) for c in coef]))
print("   basis limits: " + " ; ".join("B_%d -> %s + (%s)P" % (k, a, b)
                                       for k, (a, b) in enumerate(betas)))
xc = tuple(int(z) for z in coef)
mins, exc = classify(list(xc))
print("   ROW           x=%-12s type %-16s excess %s   limit 1" % (str(xc), str(mins), exc))

# small subspace: beta = 0 AND alpha = 0
small = []
for x in itertools.product(range(-6, 7), repeat=DIM):
    if all(v == 0 for v in x):
        continue
    g = 0
    for v in x:
        g = gcd(g, abs(v))
    xr = tuple(v//g for v in x)
    if xr[next(i for i, v in enumerate(xr) if v)] < 0 or xr in small:
        continue
    al = sum(F(xi)*betas[i][0] for i, xi in enumerate(xr))
    be = sum(F(xi)*betas[i][1] for i, xi in enumerate(xr))
    if al == 0 and be == 0:
        small.append(xr)
print("   SMALL (limit 0) directions among |x_i|<=6 :", small[:6])
for xr in small[:3]:
    mins, exc = classify(list(xr))
    print("      x=%-12s type %-16s excess %s   limit 0" % (str(xr), str(mins), exc))

# companion directions (beta != 0): find the best type
best = None
for x in itertools.product(range(-4, 5), repeat=DIM):
    if all(v == 0 for v in x):
        continue
    g = 0
    for v in x:
        g = gcd(g, abs(v))
    xr = tuple(v//g for v in x)
    be = sum(F(xi)*betas[i][1] for i, xi in enumerate(xr))
    if be == 0:
        continue
    mins, exc = classify(list(xr))
    if not mins:
        cand = (99, [], xr, exc)
    else:
        cand = (min(m[0]+m[1] for m in mins), mins, xr, exc)
    if best is None or cand[0] < best[0]:
        best = cand
print("   BEST COMPANION (period coeff != 0) over |x_i|<=4 :")
if best is None:
    print("      NO direction in the whole solution space has a nonzero period coefficient")
elif best[0] == 99:
    print("      NONE of lcm type; e.g. x=%s has UNBOUNDED excess over lcm^4 lcm^4" % str(best[2]))
else:
    al = sum(F(xi)*betas[i][0] for i, xi in enumerate(best[2]))
    be = sum(F(xi)*betas[i][1] for i, xi in enumerate(best[2]))
    print("      x=%-12s type %-16s excess %s   limit %s + (%s)P"
          % (str(best[2]), str(best[1]), best[3], al, be))
