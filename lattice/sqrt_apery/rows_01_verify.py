"""
rows_01_verify.py -- exact verification of the Sym^1 square roots of the EIGHT
non-Apery sporadic order-3 rows (Domb, T, AZ(9,3,-27), AZ(11,5,125), AZ(7,3,81),
Cooper s7, s10, s18).

For each row:
  * build the parent A_n exactly from its order-3 recurrence to n = NMAX
  * form the formal power series sqrt(sum A_n t^n) with exact rationals
  * set a_n = lambda^n [t^n] sqrt(F); check a_n in Z
  * check a_n satisfies (n+1)^2 a_{n+1} = p1(n) a_n - p0(n) a_{n-1} exactly
  * check sharpness of lambda (lambda/p fails for each prime p | lambda)

Parent recurrences used
  AZ normalisation (a,b,c):
      (n+1)^3 A_{n+1} = (2n+1)(a n^2 + a n + b) A_n - c n^3 A_{n-1},  A_0=1, A_1=b
  Cooper (source: lattice/census/cooper.gp, itself citing
          sporadic_eisenstein_cooper_research_notes.txt lines 346-359):
      s7 : (n+1)^3 A_{n+1} = (2n+1)(13n^2+13n+4) A_n + 3 n (9n^2-1)  A_{n-1}, A_1=4
      s10: (n+1)^3 A_{n+1} = 2(2n+1)(3n^2+3n+1) A_n + 4 n (16n^2-1) A_{n-1}, A_1=2
      s18: (n+1)^3 A_{n+1} = 2(2n+1)(7n^2+7n+3) A_n - 12 n (16n^2-1) A_{n-1}, A_1=6
"""
from fractions import Fraction as Fr
import sys

NMAX = 300

def parent_AZ(a, b, c, N):
    A = [Fr(1), Fr(b)]
    for n in range(1, N):
        A.append(((2*n+1)*(a*n*n + a*n + b)*A[n] - c*n**3*A[n-1]) / Fr((n+1)**3))
    return A

def parent_cooper(name, N):
    if name == 's7':
        A = [Fr(1), Fr(4)]
        f = lambda n, u1, u0: (2*n+1)*(13*n*n+13*n+4)*u1 + 3*n*(9*n*n-1)*u0
    elif name == 's10':
        A = [Fr(1), Fr(2)]
        f = lambda n, u1, u0: 2*(2*n+1)*(3*n*n+3*n+1)*u1 + 4*n*(16*n*n-1)*u0
    elif name == 's18':
        A = [Fr(1), Fr(6)]
        f = lambda n, u1, u0: 2*(2*n+1)*(7*n*n+7*n+3)*u1 - 12*n*(16*n*n-1)*u0
    for n in range(1, N):
        A.append(f(n, A[n], A[n-1]) / Fr((n+1)**3))
    return A

def sqrt_series(A):
    """c with c^2 = A as power series, c_0 = 1 (A_0 must be 1)."""
    assert A[0] == 1
    N = len(A) - 1
    c = [Fr(1)]
    for n in range(1, N+1):
        s = sum(c[k]*c[n-k] for k in range(1, n))
        c.append((A[n] - s) / 2)
    return c

# name -> (parent kind/params, lambda, p1 coefficients (A,Ahalf,B), p0 as callable, p0 string)
ROWS = [
  ("Domb (10,4,64)",     ("AZ", 10, 4, 64),   1, (20,10,2),  lambda n: 16*(2*n-1)**2, "16(2n-1)^2"),
  ("T (12,4,16)",        ("AZ", 12, 4, 16),   1, (24,12,2),  lambda n: 4*(2*n-1)**2,  "4(2n-1)^2"),
  ("AZ(9,3,-27)",        ("AZ", 9, 3, -27),   4, (72,36,6),  lambda n: -108*(2*n-1)**2, "-108(2n-1)^2"),
  ("AZ(11,5,125)",       ("AZ", 11, 5, 125),  4, (88,44,10), lambda n: 500*(2*n-1)**2,  "500(2n-1)^2"),
  ("AZ(7,3,81)",         ("AZ", 7, 3, 81),    4, (56,28,6),  lambda n: 324*(2*n-1)**2,  "324(2n-1)^2"),
  ("Cooper s7",          ("C", "s7"),         1, (26,13,2),  lambda n: -3*(3*n-1)*(3*n-2), "-3(3n-1)(3n-2)"),
  ("Cooper s10",         ("C", "s10"),        2, (24,12,2),  lambda n: -4*(8*n-3)*(8*n-5), "-4(8n-3)(8n-5)"),
  ("Cooper s18",         ("C", "s18"),        2, (56,28,6),  lambda n: 12*(8*n-3)*(8*n-5), "12(8n-3)(8n-5)"),
]

out = {}
for name, par, lam, (A2,A1,B), p0, p0s in ROWS:
    if par[0] == "AZ":
        A = parent_AZ(par[1], par[2], par[3], NMAX)
    else:
        A = parent_cooper(par[1], NMAX)
    assert all(x.denominator == 1 for x in A), (name, "parent not integral!")
    c = sqrt_series(A)
    a = [Fr(lam)**n * c[n] for n in range(NMAX+1)]
    int_ok = all(x.denominator == 1 for x in a)
    bad = [n for n in range(NMAX+1) if a[n].denominator != 1]
    # recurrence check
    rec_ok = True; recbad = None
    for n in range(1, NMAX):
        lhs = Fr((n+1)**2) * a[n+1]
        rhs = (A2*n*n + A1*n + B)*a[n] - p0(n)*a[n-1]
        if lhs != rhs:
            rec_ok = False; recbad = n; break
    a1_ok = (a[1] == B)
    # sharpness: for each prime p | lambda, does lambda/p fail?
    sharp = []
    for p in (2,):
        if lam % p == 0 and lam > 1:
            mu = lam // p
            aa = [Fr(mu)**n * c[n] for n in range(NMAX+1)]
            firstfail = next((n for n in range(NMAX+1) if aa[n].denominator != 1), None)
            sharp.append((mu, firstfail))
    out[name] = dict(lam=lam, int_ok=int_ok, bad=bad[:5], rec_ok=rec_ok, recbad=recbad,
                     a1_ok=a1_ok, sharp=sharp, a_head=[int(x) for x in a[:6]],
                     A_head=[int(x) for x in A[:6]])
    print(f"{name:18s} lam={lam}  a_n in Z to {NMAX}: {int_ok}   recurrence exact to {NMAX}: {rec_ok}"
          + ("" if rec_ok else f" FIRST FAIL n={recbad}")
          + f"   a_1==B: {a1_ok}")
    print(f"   parent A_n = {[int(x) for x in A[:6]]}")
    print(f"   a_n        = {[int(x) for x in a[:8]]}")
    print(f"   p1 = {A2}n^2+{A1}n+{B}   p0 = {p0s}")
    for mu, ff in sharp:
        print(f"   sharpness: lambda={mu} first non-integral at n={ff}  -> lambda={lam} is "
              + ("SHARP" if ff is not None else "NOT SHARP"))
    if lam == 1:
        print("   sharpness: lambda=1 is minimal, trivially sharp")
    # denominator valuation profile of c_n (2-adic) for lam in {1,2}
    if lam in (1,2):
        vs = []
        for n in range(NMAX+1):
            d = c[n].denominator
            v = 0
            while d % 2 == 0: d//=2; v+=1
            vs.append(v)
        print(f"   v_2(den [t^n]sqrt F): max over n<=300 = {max(vs)}; "
              f"max(v - n) = {max(vs[n]-n for n in range(NMAX+1))}; "
              f"first n with v=n: {next((n for n in range(NMAX+1) if vs[n]==n), None)}")
        out[name]['v2den'] = vs
    print()
