"""Task 2(a,b): the arcsine module F_e = sum y^n/(n^e binom(2n,n)) and its
rescaling F_e(y/s); exact denominator types to n=200."""
from fractions import Fraction as F
from math import comb
from qseries import lcm_upto
from dtypes import divides_type

N = 201
print("=== binom(2n,n) | [1..2n] ? (n <= 200) ===")
bad = [n for n in range(1, N) if lcm_upto(2*n) % comb(2*n, n)]
print("  counterexamples:", bad if bad else "none  -> ONE layer of rate 2 suffices  [verified]")
print("=== is the rate 2 sharp, i.e. binom(2n,n) does NOT divide [1..cn] for c<2 ? ===")
for c in (1, 1.5, 1.8, 1.9, 1.99):
    bad = [n for n in range(2, 60) if lcm_upto(int(c*n)) % comb(2*n, n)]
    print(f"  c={c}: first failure at n={bad[0] if bad else None}  ({'rate 2 needed' if bad else 'ok'})")

print()
print("=== F_e(y/s) = sum_{n>=1} lam^n y^n /(n^e binom(2n,n)),  lam = 1/s ===")
for lam, host in ((1, 'CDT / Apery-perfect (s=1, branch y=4)'),
                  (4, 'Catalan level 8 (s=1/4, branch y=1)')):
    print(f"-- lam_2 = {lam}: {host}")
    for e in range(0, 6):
        c = [F(0)] + [F(lam**n, n**e*comb(2*n, n)) for n in range(1, N)]
        found = None
        for ee in range(0, e+2):
            for bs in ((), (2,), (1,), (2, 2)):
                ok, _ = divides_type(c, ee, bs, 1, N)
                if ok: found = (ee, bs); break
            if found: break
        # sharpness of e: does e-1 work?
        sharp = None
        if found and found[0] > 0:
            ok, nn = divides_type(c, found[0]-1, found[1], 1, N)
            sharp = 'sharp' if not ok else f'NOT sharp (e-1 works)'
        print(f"   F_{e}: minimal type  n^{found[0]} * {'*'.join('[1..%dn]'%b for b in found[1]) or '1'}"
              f"   sigma={sum(found[1])}  {sharp or ''}")
