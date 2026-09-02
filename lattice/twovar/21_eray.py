"""21: E_D1 ray a=2b and a=3b, and z3_D4 a=5b : convergent or not?"""
import pickle
from fractions import Fraction as F
import lib2v
import mpmath as mp
mp.mp.dps = 40
for name, x, lam in [("E_D1", (1, 0, 0), 2), ("E_D1", (1, 0, 0), 3),
                     ("z3_D4", (0, 1, 0), 5), ("z3_D4", (0, 1, 0), 4)]:
    Dd = pickle.load(open("sol_%s.pkl" % name, "rb"))
    Ts, NBIG = Dd["Ts"], Dd["NBIG"]
    cfun = getattr(lib2v, name)
    seq = []
    for b in range(2, NBIG//lam + 1):
        a = lam*b
        q = sum(F(xi)*T[(a, b)] for xi, T in zip(x, Ts) if xi)/F(cfun(a, b))
        seq.append(mp.mpf(q.numerator)/mp.mpf(q.denominator))
    print("%s lam=%d  b=2..%d :" % (name, lam, len(seq)+1))
    print("   ", [mp.nstr(v, 10) for v in seq[-8:]])
    if len(seq) >= 6:
        try:
            v, e = mp.richardson(lambda k: seq[int(k)-1], 1, min(len(seq), 10))
            print("    richardson ->", mp.nstr(v, 12))
        except Exception as ex:
            print("    richardson failed", ex)
        try:
            sh = mp.shanks(seq)
            print("    shanks     ->", mp.nstr(sh[-1][-1], 12))
        except Exception as ex:
            print("    shanks failed", ex)
