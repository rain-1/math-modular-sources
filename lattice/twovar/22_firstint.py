"""22: the first-order (fibration) PDE of the s10_D1 and z3_D2 systems, and its
first integral u(x,y).  If u = alpha(x) beta(y) the Landau ceiling of section 3.2
applies; if u is genuinely non-separated it does not."""
from lib2v import *
from fractions import Fraction as F
import sympy as sp
import math

PRIMES = [1000003, 999983, 999979]


def nullspace_modp(rows, n, p):
    M = [r[:] for r in rows]
    piv, r = [], 0
    for c in range(n):
        pr = next((i for i in range(r, len(M)) if M[i][c] % p), None)
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        inv = pow(M[r][c], p-2, p)
        M[r] = [(v*inv) % p for v in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] % p:
                f = M[i][c]
                M[i] = [(M[i][j]-f*M[r][j]) % p for j in range(n)]
        piv.append(c)
        r += 1
    free = [c for c in range(n) if c not in piv]
    out = []
    for fc in free:
        v = [0]*n
        v[fc] = 1
        for i, pc in enumerate(piv):
            v[pc] = (-M[i][fc]) % p
        out.append(v)
    return out


def ratrec(a, m):
    bound = math.isqrt(m//2)
    r0, r1, s0, s1 = m, a % m, 0, 1
    while r1 > bound:
        q = r0//r1
        r0, r1 = r1, r0-q*r1
        s0, s1 = s1, s0-q*s1
    if s1 == 0 or math.gcd(r1, s1) != 1:
        return None
    return F(r1 if s1 > 0 else -r1, abs(s1))


S = [(i, j) for i in range(3) for j in range(3)]
MONS = [(0, 0), (1, 0), (0, 1)]     # coefficients linear in (a,b) -> order 1 in theta


def exact_rel(cf, name):
    n = len(S)*len(MONS)
    vecs = []
    for p in PRIMES:
        tab = [[cf(a, b) % p for b in range(30)] for a in range(30)]
        rows = []
        for a in range(1, 26):
            for b in range(1, 26):
                row = []
                for (i, j) in S:
                    val = tab[a+i][b+j]
                    for (u, v) in MONS:
                        row.append((val*pow(a, u, p)*pow(b, v, p)) % p)
                rows.append(row)
        ns = nullspace_modp(rows, n, p)
        if len(ns) != 1:
            print("  %s: nullspace dim %d, abort" % (name, len(ns)))
            return None
        vecs.append((p, ns[0]))
    M = 1
    for p, _ in vecs:
        M *= p
    out = []
    for k in range(n):
        x = 0
        for p, v in vecs:
            Mp = M//p
            x = (x + v[k]*Mp*pow(Mp, -1, p)) % M
        out.append(ratrec(x, M))
    return out


a, b, x, y = sp.symbols('a b x y')
tx, ty = sp.symbols('theta_x theta_y')
for cf, name in [(s10_D1, "s10_D1  sum C(a,k)^2 C(b,k)^2"),
                 (z3_D2, "z3_D2   sum C(a,k)^2 C(b+k,k)^2"),
                 (z3_D1, "z3_D1   (control: known pullback along xy/((1-x)^2(1-y)^2))")]:
    print("=" * 74)
    print(name)
    v = exact_rel(cf, name)
    if v is None:
        continue
    # relation sum_{(i,j)} (c0 + c1 a + c2 b) c_{a+i,b+j} = 0
    # generating function: c_{a+i,b+j} <-> x^{-i} y^{-j} C ; a <-> theta_x - ? ;
    # we only need A,B,E in  A theta_x C + B theta_y C + E C = 0 up to x^2y^2:
    A = sp.S(0)
    B = sp.S(0)
    E = sp.S(0)
    for idx, (i, j) in enumerate(S):
        c0, c1, c2 = v[3*idx], v[3*idx+1], v[3*idx+2]
        w = x**(2-i)*y**(2-j)
        # sum_{a,b} a c_{a+i,b+j} x^a y^b = x^{-i}y^{-j} (theta_x - i) C
        A += sp.Rational(c1.numerator, c1.denominator)*w
        B += sp.Rational(c2.numerator, c2.denominator)*w
        E += (sp.Rational(c0.numerator, c0.denominator)
              - sp.Rational(c1.numerator, c1.denominator)*i
              - sp.Rational(c2.numerator, c2.denominator)*j)*w
    A, B, E = map(sp.expand, (A, B, E))
    g = sp.gcd(sp.gcd(sp.Poly(A, x, y), sp.Poly(B, x, y)), sp.Poly(E, x, y))
    print("   A(x,y) =", sp.factor(A))
    print("   B(x,y) =", sp.factor(B))
    print("   E(x,y) =", sp.factor(E))
    print("   the foliation is  dx/(x A) = dy/(y B) ;  separable iff A/B splits:")
    print("   A/B =", sp.simplify(A/B))
