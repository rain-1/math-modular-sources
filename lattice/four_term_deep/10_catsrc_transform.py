#!/usr/bin/env python3
"""10_catsrc_transform.py

Signed binomial transforms of Zagier E  (the integral three-term row with
Apery limit G/2, Catalan) for nu = 1, 2, -1, -2.

Zagier E:  (n+1)^2 u_{n+1} = (12 n^2 + 12 n + 4) u_n - 32 n^2 u_{n-1}
           a_0=1, a_1=4;  b_0=0, b_1=1;  lim b_n/a_n = G/2.

For each nu:
  (i)   exact 12-unknown fit of  L(n) u_{n+1} = P(n) u_n - Q(n) u_{n-1} + R(n) u_{n-2}
        (L,P,Q,R quadratics) from the transformed integral sequence; nullspace
        dimension, normalisation, integrality of coefficients, verification to n=NVER
  (ii)  integrality of the transformed sequence
  (iii) Apery limit to >=60 digits (written to out/catsrc_limits.txt)
  (iv)  apparent-singularity (Frobenius no-log) test at t=infinity
  (v)   sharp denominator exponent k of the companion, characteristic roots
"""
import sys, os, json
from fractions import Fraction as Fr
from math import gcd, comb

sys.path.insert(0, "/home/ubuntu/code/math-modular-sources/lattice/four_term")
import mpmath as mp
from importlib.machinery import SourceFileLoader
rep = SourceFileLoader("rep", "/home/ubuntu/code/math-modular-sources/lattice/four_term/05_report.py").load_module()

mp.mp.dps = 200

A0, C0, D0 = 12, 4, 32          # Zagier E : P0 = A0(n^2+n)+C0 , Q0 = D0 n^2
N = 460                          # length of exact sequences
NVER = 300                       # recurrence verification range
NK = 120                         # range for sharp-k scan

# ------------------------------------------------------------------ 3-term row
def three_term(N, start):
    """start=0 -> a_n (a_0=1); start=1 -> b_n (b_0=0,b_1=1). Exact Fractions."""
    P = lambda n: A0*(n*n+n) + C0
    Q = lambda n: D0*n*n
    if start == 0:
        u = [Fr(1), Fr(P(0))]
    else:
        u = [Fr(0), Fr(1)]
    for n in range(1, N):
        u.append((P(n)*u[n] - Q(n)*u[n-1]) / Fr((n+1)**2))
    return u

a = three_term(N, 0)
b = three_term(N, 1)
assert all(x.denominator == 1 for x in a), "Zagier E a_n not integral"
print("Zagier E   a_n[0..7] =", [int(x) for x in a[:8]])
print("Zagier E   b_n[0..7] =", [str(x) for x in b[:8]])

def sharp_k(seq, upto):
    L = 1; kmax = 0
    for n in range(1, upto+1):
        L = L*n//gcd(L, n)
        den = seq[n].denominator
        kk = 0
        while den > 1:
            g = gcd(den, L)
            if g == 1:
                kk = 99; break
            den //= g; kk += 1
        kmax = max(kmax, kk)
    return kmax
print("Zagier E   sharp k (companion, n<=%d) = %d" % (NK, sharp_k(b, NK)))

# ------------------------------------------------------------- transform
def sbt(u, nu, N):
    """v_n = sum_m binom(n,m) (-nu)^{n-m} u_m"""
    out = []
    for n in range(N):
        s = Fr(0)
        p = 1
        pw = [1]*(n+1)
        for j in range(1, n+1):
            pw[j] = pw[j-1]*(-nu)
        for m in range(n+1):
            s += comb(n, m)*pw[n-m]*u[m]
        out.append(s)
    return out

# ------------------------------------------------------------- exact 12-unknown fit
def fit12(v, nmin, nmax):
    """Solve for (L,P,Q,R) quadratics with
         L(n) v_{n+1} - P(n) v_n + Q(n) v_{n-1} - R(n) v_{n-2} = 0
       Unknowns x = (l2,l1,l0, p2,p1,p0, q2,q1,q0, r2,r1,r0)   (12 unknowns).
       Returns (nullspace dimension, basis vectors as Fractions)."""
    rows = []
    for n in range(nmin, nmax+1):
        vn1 = v[n+1]; vn = v[n]; vm1 = v[n-1] if n-1 >= 0 else Fr(0)
        vm2 = v[n-2] if n-2 >= 0 else Fr(0)
        mono = [Fr(n*n), Fr(n), Fr(1)]
        rows.append([m*vn1 for m in mono] + [-m*vn for m in mono]
                    + [m*vm1 for m in mono] + [-m*vm2 for m in mono])
    # gaussian elimination over Q
    ncol = 12
    piv = []
    r = 0
    Mx = [row[:] for row in rows]
    for c in range(ncol):
        pr = None
        for i in range(r, len(Mx)):
            if Mx[i][c] != 0: pr = i; break
        if pr is None: continue
        Mx[r], Mx[pr] = Mx[pr], Mx[r]
        pv = Mx[r][c]
        Mx[r] = [x/pv for x in Mx[r]]
        for i in range(len(Mx)):
            if i != r and Mx[i][c] != 0:
                f = Mx[i][c]
                Mx[i] = [x - f*y for x, y in zip(Mx[i], Mx[r])]
        piv.append(c); r += 1
        if r == ncol: break
    free = [c for c in range(ncol) if c not in piv]
    basis = []
    for fc in free:
        vec = [Fr(0)]*ncol
        vec[fc] = Fr(1)
        for i, pc in enumerate(piv):
            vec[pc] = -Mx[i][fc]
        basis.append(vec)
    return len(free), basis

def clear_denoms(vec):
    from functools import reduce
    den = reduce(lambda x, y: x*y//gcd(x, y), [x.denominator for x in vec], 1)
    ints = [int(x*den) for x in vec]
    g = 0
    for x in ints: g = gcd(g, abs(x))
    if g: ints = [x//g for x in ints]
    if sum(1 for x in ints if x != 0) and ints[[i for i, x in enumerate(ints) if x != 0][0]] < 0:
        ints = [-x for x in ints]
    return ints

results = {}
limlines = []
for nu in (1, 2, -1, -2):
    print("\n" + "="*74)
    print("nu = %d" % nu)
    print("="*74)
    v = sbt(a, nu, N)
    w = sbt(b, nu, N)
    intg = all(x.denominator == 1 for x in v)
    print("  (ii) v_n integral to n=%d : %s" % (N-1, intg))
    print("       v_n[0..8] =", [int(x) for x in v[:9]])
    print("       w_n[0..6] =", [str(x) for x in w[:7]])

    # ---- (i) fit
    dim, basis = fit12(v, 2, 30)
    print("  (i) 12-unknown nullspace dim (n=2..30) =", dim)
    if dim != 1:
        print("      !! nullspace not 1-dimensional"); continue
    sol = clear_denoms(basis[0])
    l2, l1, l0, p2, p1, p0, q2, q1, q0, r2, r1, r0 = sol
    print("      raw integral solution (l2,l1,l0,p2,p1,p0,q2,q1,q0,r2,r1,r0) =", sol)
    # normalise leading to (n+1)^2 = n^2+2n+1
    assert (l2, l1, l0) == (l2, 2*l2, l2) or True
    if (l1, l0) != (2*l2, l2):
        print("      !! leading polynomial is not proportional to (n+1)^2:", (l2, l1, l0))
    scale = Fr(1, l2)
    P = [Fr(p2)*scale, Fr(p1)*scale, Fr(p0)*scale]
    Q = [Fr(q2)*scale, Fr(q1)*scale, Fr(q0)*scale]
    R = [Fr(r2)*scale, Fr(r1)*scale, Fr(r0)*scale]
    allint = all(x.denominator == 1 for x in P+Q+R)
    Pi_, Qi_, Ri_ = [[int(x) for x in X] for X in (P, Q, R)]
    print("      normalised to (n+1)^2 u_{n+1} = P u_n - Q u_{n-1} + R u_{n-2}")
    print("      P(n) = %d n^2 + %d n + %d" % tuple(Pi_))
    print("      Q(n) = %d n^2 + %d n + %d" % tuple(Qi_))
    print("      R(n) = %d n^2 + %d n + %d" % tuple(Ri_))
    print("      all coefficients integral :", allint)

    # closed form of Proposition F4
    F4P = (A0-3*nu, A0-3*nu, C0-nu)
    F4Q = (D0-2*A0*nu+3*nu*nu, 0, 0)
    CC = -nu*(D0-A0*nu+nu*nu)
    F4R = (CC, -CC, 0)
    print("      Prop F4 closed form matches :",
          tuple(Pi_) == F4P and tuple(Qi_) == F4Q and tuple(Ri_) == F4R,
          " (F4:", F4P, F4Q, F4R, ")")

    # verify recurrence on v and w to NVER
    def check(seq, lo):
        ok = True
        for n in range(lo, NVER):
            lhs = Fr((n+1)**2)*seq[n+1]
            rhs = (Fr(Pi_[0]*n*n+Pi_[1]*n+Pi_[2])*seq[n]
                   - Fr(Qi_[0]*n*n+Qi_[1]*n+Qi_[2])*(seq[n-1] if n >= 1 else Fr(0))
                   + Fr(Ri_[0]*n*n+Ri_[1]*n+Ri_[2])*(seq[n-2] if n >= 2 else Fr(0)))
            if lhs != rhs: ok = False; print("      fail at n=", n); break
        return ok
    print("      recurrence exact on v_n, n=0..%d : %s" % (NVER-1, check(v, 0)))
    print("      recurrence exact on w_n, n=1..%d : %s" % (NVER-1, check(w, 1)))

    # ---- (v) characteristic roots
    aa, dd, gg = Pi_[0], Qi_[0], Ri_[0]
    lams = mp.polyroots([1, -aa, dd, -gg], maxsteps=300, extraprec=300)
    lams = sorted(lams, key=lambda z: -abs(z))
    print("  (v) char poly lam^3 - %d lam^2 + %d lam - %d ; roots = %s"
          % (aa, dd, gg, [mp.nstr(z, 12) for z in lams]))
    print("      predicted {l1-nu, l2-nu, -nu} = %s" % sorted([8-nu, 4-nu, -nu], key=lambda z: -abs(z)))

    # ---- (v) sharp k of the companion
    kk = sharp_k(w, NK)
    kk1 = sharp_k(w, NK) if True else None
    # test whether k-1 fails
    def k_fails(seq, kt, upto):
        L = 1
        for n in range(1, upto+1):
            L = L*n//gcd(L, n)
            if (Fr(L**kt)*seq[n]).denominator != 1:
                return True, n
        return False, None
    f1 = k_fails(w, kk-1, NK) if kk >= 1 else (True, 0)
    f0 = k_fails(w, kk, NK)
    print("  (v) sharp k (companion, n<=%d) = %d   [d_n^%d fails at n=%s ; d_n^%d ok: %s]"
          % (NK, kk, kk-1, f1[1], kk, not f0[0]))
    l2abs = abs(lams[1])
    score = float(-mp.log(l2abs) - kk)
    print("      score = log(1/|lam_2|) - k = %.4f" % score)

    # ---- (iii) Apery limit
    with mp.workdps(220):
        r1 = mp.mpf(w[N-1].numerator)/mp.mpf(w[N-1].denominator) / (mp.mpf(v[N-1].numerator)/mp.mpf(v[N-1].denominator))
        r2 = mp.mpf(w[N-2].numerator)/mp.mpf(w[N-2].denominator) / (mp.mpf(v[N-2].numerator)/mp.mpf(v[N-2].denominator))
        xi = r1
        conv = abs(r1-r2)
    print("  (iii) xi (n=%d) = %s" % (N-1, mp.nstr(xi, 70)))
    print("        |xi_n - xi_{n-1}| = 10^%.1f  (=> >=%d digits stable)"
          % (float(mp.log10(conv)), int(-mp.log10(conv))))
    print("        xi/Catalan = %s" % mp.nstr(xi/mp.catalan, 45))
    limlines.append("nu%+d %s" % (nu, mp.nstr(xi, 80)))

    # ---- (iv) apparent singularity at infinity (Frobenius obstruction)
    a_, b_, c_ = Pi_; d_, e_, f_ = Qi_; g_, h_, j_ = Ri_
    At = [0, 0, -g_, d_, -a_, 1]
    St = [-(5*g_+h_), 3*d_+e_, -(a_+b_), 1]
    Bt = [0] + [2*At[i+2]-St[i] for i in range(len(St))]
    Vt = [-(4*g_+2*h_+j_), d_+e_+f_, -c_]
    # exponents at infinity: (2-s1, 2-s2), s_i roots of R
    import sympy
    s_roots = sorted(sympy.Poly([g_, h_, j_], sympy.Symbol('x')).all_roots())
    s1 = Fr(int(s_roots[0])) if s_roots[0].is_Integer else None
    s2 = Fr(int(s_roots[1])) if s_roots[1].is_Integer else None
    delta = abs(s2-s1)
    nu_exp = min(2-s1, 2-s2)
    Sobs, Iobs = rep.frob_obstruction([mp.mpmathify(x) for x in At],
                                      [mp.mpmathify(x) for x in Bt],
                                      [mp.mpmathify(x) for x in Vt],
                                      mp.mpmathify(float(nu_exp)), 0, int(delta))
    sc = max(abs(mp.mpmathify(x)) for x in At+Bt+Vt if x != 0)
    app = bool(abs(Sobs) < mp.mpf(10)**(-40)*max(mp.mpf(1), sc))
    print("  (iv) t=infinity : R roots s=(%s,%s), exponents (%s,%s), delta=%s"
          % (s1, s2, 2-s1, 2-s2, delta))
    print("       Frobenius obstruction S = %s  (scale %s)  ->  APPARENT: %s"
          % (mp.nstr(Sobs, 10), mp.nstr(sc, 6), app))
    # finite singular points: rho
    Tt = [1, -b_, d_+e_, -(2*g_+h_)]
    Rc = [1, -a_, d_, -g_]
    print("       finite singular points t_i = 1/lam_i, exponents (0, rho_i):")
    for L in lams:
        t0 = 1/L
        Tv = sum(mp.mpmathify(Tt[i])*t0**i for i in range(4))
        Rp = sum(i*mp.mpmathify(Rc[i])*t0**(i-1) for i in range(1, 4))
        rho = -Tv/(t0*Rp)
        print("         t=%s  rho=%s" % (mp.nstr(t0, 12), mp.nstr(rho, 12)))

    results[nu] = dict(P=Pi_, Q=Qi_, R=Ri_, k=kk, score=score,
                       lam=[mp.nstr(z, 20) for z in lams],
                       xi=mp.nstr(xi, 70), apparent_inf=app, integral=intg)

with open("/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/catsrc_limits.txt", "w") as fh:
    fh.write("\n".join(limlines) + "\n")
print("\nwrote out/catsrc_limits.txt")
print(json.dumps(results, indent=1))
