#!/usr/bin/env python3
"""
K-POLYGON: exact 2-adic Smith polygon of the E/Z-adapted Hermite evaluation map.

Implements handoff PACKET_K_RESOLUTION_RESEARCH_HANDOFF_2026-08-27.txt sections
XV (Steps 1-3) and XIX, using exact rational arithmetic (fractions.Fraction).

Deterministic, no internet.  Run:  python3 k_polygon.py
Output: verification log + polygon tables (stdout).
"""

from fractions import Fraction
import sys

# ----------------------------------------------------------------------
# 0.  basic 2-adic utilities
# ----------------------------------------------------------------------

def v2(x):
    """2-adic valuation of a nonzero rational."""
    x = Fraction(x)
    if x == 0:
        raise ValueError("v2(0) undefined")
    n, d = x.numerator, x.denominator
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    while d % 2 == 0:
        d //= 2
        v -= 1
    return v


def s2(n):
    return bin(n).count("1")


def S(D):
    return sum(s2(m) for m in range(1, D + 1))


# ----------------------------------------------------------------------
# 1.  truncated power series over Q  (lists of Fractions, index = exponent)
# ----------------------------------------------------------------------

def ser_mul(a, b, N):
    out = [Fraction(0)] * (N + 1)
    for i, ai in enumerate(a[:N + 1]):
        if ai == 0:
            continue
        for j, bj in enumerate(b[:N + 1 - i]):
            if bj:
                out[i + j] += ai * bj
    return out


def ser_inv(a, N):
    """1/a, requires a[0] != 0."""
    assert a[0] != 0
    out = [Fraction(0)] * (N + 1)
    out[0] = 1 / Fraction(a[0])
    for n in range(1, N + 1):
        acc = Fraction(0)
        for k in range(1, n + 1):
            if k < len(a) and a[k]:
                acc += a[k] * out[n - k]
        out[n] = -acc * out[0]
    return out


def ser_compose(f, w, N):
    """f(w(z)) with w[0] == 0."""
    assert w[0] == 0
    out = [Fraction(0)] * (N + 1)
    cur = [Fraction(0)] * (N + 1)
    cur[0] = Fraction(1)            # w^0
    for k in range(0, N + 1):
        if k < len(f) and f[k]:
            for i in range(N + 1):
                if cur[i]:
                    out[i] += f[k] * cur[i]
        if k < N:
            cur = ser_mul(cur, w, N)
    return out


def ser_sqrt_1_minus(c, N):
    """(1 - c*z)^(1/2) as a series in z, via the binomial series."""
    out = [Fraction(0)] * (N + 1)
    coef = Fraction(1)              # binom(1/2, k)
    for k in range(N + 1):
        out[k] = coef * Fraction(-c) ** k
        coef = coef * (Fraction(1, 2) - k) / (k + 1)
    return out


# ----------------------------------------------------------------------
# 2.  Beukers polynomials evaluated at a fixed rational x
# ----------------------------------------------------------------------

def beukers(x, N):
    """Return (q, p) lists of Fractions, q[n]=q_n(x), p[n]=p_n(x), n=0..N."""
    x = Fraction(x)
    c = 1 - x + x * x
    q = [Fraction(1), c]
    p = [Fraction(0), Fraction(1)]
    for n in range(1, N):
        A = (2 * n * (n + 1) + c)
        q.append((A * q[n] - n * n * q[n - 1]) / Fraction((n + 1) ** 2))
        p.append((A * p[n] - n * n * p[n - 1]) / Fraction((n + 1) ** 2))
    return q[:N + 1], p[:N + 1]


def binom(n, k):
    from math import comb
    return comb(n, k)


# ----------------------------------------------------------------------
# 3.  2-adic Smith exponents over Z_(2)
# ----------------------------------------------------------------------

def smith_v2(Min):
    """
    Smith exponents e_1 <= ... <= e_r over Z_(2) of a rational matrix.
    Gaussian elimination always pivoting on the entry of minimal v_2.
    Returns sorted list of pivot valuations.
    """
    M = [[Fraction(x) for x in row] for row in Min]
    n = len(M)
    m = len(M[0]) if n else 0
    rows = list(range(n))
    cols = list(range(m))
    exps = []
    while rows and cols:
        best = None
        for i in rows:
            Mi = M[i]
            for j in cols:
                e = Mi[j]
                if e:
                    v = v2(e)
                    if best is None or v < best[0]:
                        best = (v, i, j)
        if best is None:
            break                       # remaining block is zero: rank deficient
        v, pi, pj = best
        exps.append(v)
        piv = M[pi][pj]
        prow = M[pi]
        for i in rows:
            if i == pi:
                continue
            f = M[i][pj]
            if f:
                f = f / piv
                Mi = M[i]
                for j in cols:
                    if prow[j]:
                        Mi[j] -= f * prow[j]
        rows.remove(pi)
        cols.remove(pj)
    return sorted(exps)


def matmul(A, B):
    n = len(A); k = len(B); m = len(B[0])
    C = [[Fraction(0)] * m for _ in range(n)]
    for i in range(n):
        Ai = A[i]; Ci = C[i]
        for t in range(k):
            a = Ai[t]
            if a:
                Bt = B[t]
                for j in range(m):
                    if Bt[j]:
                        Ci[j] += a * Bt[j]
    return C


def det_frac(A):
    """Exact determinant by fraction-free-ish Gaussian elimination."""
    n = len(A)
    M = [[Fraction(x) for x in row] for row in A]
    det = Fraction(1)
    for c in range(n):
        piv = None
        for r in range(c, n):
            if M[r][c]:
                piv = r
                break
        if piv is None:
            return Fraction(0)
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            det = -det
        det *= M[c][c]
        inv = 1 / M[c][c]
        for r in range(c + 1, n):
            if M[r][c]:
                f = M[r][c] * inv
                for j in range(c, n):
                    if M[c][j]:
                        M[r][j] -= f * M[c][j]
    return det


# ----------------------------------------------------------------------
# 4.  reporting helpers
# ----------------------------------------------------------------------

LOG = []


def say(s=""):
    LOG.append(s)
    print(s)
    sys.stdout.flush()


PASS = []


def check(name, ok, detail=""):
    PASS.append((name, ok))
    say(("  PASS  " if ok else "  **FAIL**  ") + name + ("   " + detail if detail else ""))
    return ok


# ======================================================================
#                       MAIN
# ======================================================================

def main():
    NMAX = 80          # coefficient horizon for beta
    DMAX = 12          # largest D for the full four-formulation tables
    DEXT = 34          # largest D for closed-form / extended scans (Parts VI-VII)

    say("=" * 78)
    say("PART I.  VERIFICATION OF EVERY DEFINITION USED")
    say("=" * 78)

    # ---- 4.1 Beukers half-integer valuations -------------------------
    say()
    say("[1] v_2(q_n(1/2-m)) = -4n + 2 s_2(n),  n<=24, m=0..12")
    q_half, p_half = beukers(Fraction(1, 2), NMAX)
    ok = True
    bad = []
    for m in range(0, 13):
        qm, pm = beukers(Fraction(1, 2) - m, 24)
        for n in range(0, 25):
            if v2(qm[n]) != -4 * n + 2 * s2(n):
                ok = False
                bad.append((m, n, v2(qm[n]), -4 * n + 2 * s2(n)))
    check("A1/A5 denominator valuation (325 cases)", ok,
          "" if ok else str(bad[:5]))

    # also A3.1: v_2(p_n(1/2)) = -4n+2s_2(n)+2 for n>=1
    ok = all(v2(p_half[n]) == -4 * n + 2 * s2(n) + 2 for n in range(1, 25))
    check("A17 numerator valuation v_2(p_n(1/2)) = -4n+2s_2(n)+2 (n=1..24)", ok)

    # ---- 4.2 sqrt(1-16z) Y_0(16z) = P(64z) ---------------------------
    say()
    say("[2] sqrt(1-16z) * Y_0(16z) = P(64z) = sum binom(2n,n)^2 z^n, to z^24")
    N = 24
    Y0_at16z = [q_half[n] * Fraction(16) ** n for n in range(N + 1)]
    Y1_at16z = [p_half[n] * Fraction(16) ** n for n in range(N + 1)]
    sq = ser_sqrt_1_minus(16, N)
    lhs = ser_mul(sq, Y0_at16z, N)
    Pz = [Fraction(binom(2 * n, n)) ** 2 for n in range(N + 1)]
    ok = all(lhs[n] == Pz[n] for n in range(N + 1))
    check("K2.8 fixed-fibre bridge (exact, 25 coefficients)", ok,
          "" if ok else "first mismatch n=%d" % next(n for n in range(N + 1) if lhs[n] != Pz[n]))

    # ---- 4.3 beta ----------------------------------------------------
    say()
    say("[3] beta(z) = R0(64z)/P(64z) = (1/8) Y_1(16z)/Y_0(16z),  b_n for n<=%d" % NMAX)
    Y0f = [q_half[n] * Fraction(16) ** n for n in range(NMAX + 1)]
    Y1f = [p_half[n] * Fraction(16) ** n for n in range(NMAX + 1)]
    beta = ser_mul(Y1f, ser_inv(Y0f, NMAX), NMAX)
    beta = [c / 8 for c in beta]
    b = beta                                   # b[0] should be 0
    check("b_0 = 0", b[0] == 0)
    check("b_1 = 2  (matches R0(Y)=Y/32+O(Y^2), i.e. 64/32)", b[1] == 2,
          "b_1 = %s" % b[1])

    # cross-check beta against the explicit r_n of Section E1
    Sn = [Fraction(0)] * (NMAX + 1)
    acc = Fraction(0)
    for n in range(1, NMAX + 1):
        m = n - 1
        acc += Fraction(4) ** m / (Fraction((2 * m + 1) ** 2) * binom(2 * m, m))
        Sn[n] = acc
    r = [Fraction(0)] * (NMAX + 1)
    for n in range(1, NMAX + 1):
        r[n] = Fraction(1, 2) * Fraction(binom(2 * n, n)) ** 2 * Sn[n]
    beta2 = ser_mul(r, ser_inv(Pz + [Fraction(binom(2 * n, n)) ** 2
                                     for n in range(N + 1, NMAX + 1)], NMAX), NMAX)
    ok = all(b[n] == beta2[n] for n in range(NMAX + 1))
    check("beta from (1/8)Y_1/Y_0  ==  beta from Section E1 r_n/p_n (n<=%d)" % NMAX, ok)

    # ---- 4.4 parity theorem ------------------------------------------
    say()
    say("[4] parity theorem: v_2(b_n) >= 1, with equality iff n is a power of 2")
    ok = True
    bad = []
    for n in range(1, 61):
        vv = v2(b[n])
        want_unit = (s2(n) == 1)
        if vv < 1 or ((vv == 1) != want_unit):
            ok = False
            bad.append((n, vv))
    check("E2 ratio parity theorem (n=1..60)", ok, "" if ok else str(bad[:6]))
    # stronger: v_2(b_n) == 2 s_2(n) - 1 ?
    strong = all(v2(b[n]) == 2 * s2(n) - 1 for n in range(1, 61))
    say("      [extra] v_2(b_n) = 2 s_2(n) - 1 for all n=1..60 : %s" % strong)

    # ---- 4.5 M_D determinant -----------------------------------------
    say()
    say("[5] v_2(det M_D) = D  for D = 1..%d,  (M_D)_{ij} = b_{D+i-j}" % (DEXT))

    def bb(i):
        return b[i] if i >= 1 else Fraction(0)

    def M_block(D):
        return [[bb(D + i - j) for j in range(D)] for i in range(D)]

    ok = True
    detrow = []
    for D in range(1, DEXT + 1):
        dv = v2(det_frac(M_block(D)))
        detrow.append((D, dv))
        if dv != D:
            ok = False
    check("E3 binary Hankel determinant (D=1..%d)" % (DEXT), ok, str(detrow))

    # T_D block, for completeness of Psi_D
    def T_block(D):
        return [[bb(i - j) for j in range(D)] for i in range(D)]

    # ---- 4.6 Delta_n --------------------------------------------------
    say()
    say("[6] Delta_n = sigma_n - rho_n,   sigma_m = 2 P_m/Q_m,  rho_m = 4 B_{2m}/A_{2m}")
    say("    Q_m = q_m(1/2-m),  P_m = Q_m sigma^alt_m + (-1)^m p_m(1/2-m)/8   (A6)")
    say("    A(z) = (1-8z) Y_0(16z(1-4z)),  16 B(z) = (1-8z) Y_1(16z(1-4z))   (A13/A14)")

    NA = 2 * DEXT
    w = [Fraction(0)] * (NA + 1)
    w[1] = Fraction(16)
    if NA >= 2:
        w[2] = Fraction(-64)
    Y0t = [q_half[n] for n in range(NA + 1)]
    Y1t = [p_half[n] for n in range(NA + 1)]
    g = [Fraction(0)] * (NA + 1)
    g[0] = Fraction(1)
    g[1] = Fraction(-8)
    A_ser = ser_mul(g, ser_compose(Y0t, w, NA), NA)
    B_ser = [c / 16 for c in ser_mul(g, ser_compose(Y1t, w, NA), NA)]

    # A15/A16 sanity
    okA = all(v2(A_ser[n]) == 2 * s2(n) for n in range(0, NA + 1))
    okB = all(v2(B_ser[n]) == 2 * s2(n) - 2 for n in range(1, NA + 1))
    check("A15  v_2(A_n) = 2 s_2(n)  (n=0..%d)" % NA, okA)
    check("A16  v_2(B_n) = 2 s_2(n)-2 (n=1..%d)" % NA, okB)

    Delta = {}
    sigalt = Fraction(0)
    for m in range(1, DEXT + 1):
        qm, pm = beukers(Fraction(1, 2) - m, m)
        Qm = qm[m]
        sigalt = sum(Fraction((-1) ** j, (2 * j + 1) ** 2) for j in range(m))
        Pm = Qm * sigalt + Fraction((-1) ** m, 8) * pm[m]
        sigma = 2 * Pm / Qm
        rho = 4 * B_ser[2 * m] / A_ser[2 * m]
        Delta[m] = sigma - rho

    ok = True
    row = []
    for m in range(1, DEXT + 1):
        vv = v2(Delta[m])
        row.append((m, vv, 8 * m - 4 * s2(m)))
        if vv != 8 * m - 4 * s2(m):
            ok = False
    check("A23  v_2(Delta_n) = 8n - 4 s_2(n)  (n=1..%d)" % (DEXT), ok,
          str([(a, bq) for a, bq, c in row]))

    # A24 cross-check
    ok24 = True
    for m in range(1, DEXT + 1):
        qm, pm = beukers(Fraction(1, 2) - m, m)
        Qm = qm[m]
        sigalt = sum(Fraction((-1) ** j, (2 * j + 1) ** 2) for j in range(m))
        Pm = Qm * sigalt + Fraction((-1) ** m, 8) * pm[m]
        val = Pm * A_ser[2 * m] - 2 * B_ser[2 * m] * Qm
        if v2(val) != 4 * m - 1:
            ok24 = False
    check("A24  v_2(P_m A_{2m} - 2 B_{2m} Q_m) = 4m-1", ok24)

    # ---- 4.7 K4.6 -----------------------------------------------------
    say()
    say("[7] K4.6:  Delta_j = u_j h_j^2,  h_j = 4^j (j!)^2,  u_j a 2-adic unit")
    import math as _m
    h = {j: Fraction(4) ** j * Fraction(_m.factorial(j)) ** 2 for j in range(0, DEXT + 1)}
    ok = all(v2(Delta[j] / h[j] ** 2) == 0 for j in range(1, DEXT + 1))
    check("K4.6 unit quotient u_j = Delta_j / h_j^2 (j=1..%d)" % (DEXT), ok)

    # ---- 4.8 Omega_D --------------------------------------------------
    say()
    say("[8] Omega_D = sum_{n<=D} v_2(Delta_n) = 4D^2+4D-4S_D,  and Omega_D = 2 A_D")
    Omega = {}
    ok = True
    okA2 = True
    for D in range(1, DEXT + 1):
        Omega[D] = sum(v2(Delta[n]) for n in range(1, D + 1))
        if Omega[D] != 4 * D * D + 4 * D - 4 * S(D):
            ok = False
        AD = sum(v2(h[j]) for j in range(1, D + 1))
        if Omega[D] != 2 * AD:
            okA2 = False
    check("A25  Omega_D = 4D^2+4D-4S_D (D=1..%d)" % (DEXT), ok,
          str([(D, Omega[D]) for D in range(1, min(DMAX, 10) + 1)]))
    check("K4.9  Omega_D = 2 A_D,  A_D = v_2(det H_D) = 2D(D+1)-2S_D", okA2)

    # ---- 4.9 L_D, C_D^(0) unimodular ---------------------------------
    say()
    say("[9] K4.2/K4.4: L_D and C_D^(0) lie in GL_D(Z_(2))")

    qmov = {}
    for i in range(0, DEXT + 1):
        qmov[i] = beukers(Fraction(1, 2) - i, DEXT)[0]

    def L_mat(D):
        return [[(h[j] * qmov[i][j] if j <= i else Fraction(0))
                 for j in range(1, D + 1)] for i in range(1, D + 1)]

    def C0_mat(D):
        return [[(h[j] * q_half[j] if i == j else Fraction(0))
                 for j in range(1, D + 1)] for i in range(1, D + 1)]

    okL = True
    okC = True
    for D in range(1, DMAX + 1):
        LD = L_mat(D)
        C0 = C0_mat(D)
        if v2(det_frac(LD)) != 0:
            okL = False
        if v2(det_frac(C0)) != 0:
            okC = False
        for rr in LD:
            for e in rr:
                if e and v2(e) < 0:
                    okL = False
    check("L_D lower-triangular, 2-integral, v_2(det L_D)=0 (D=1..%d)" % DMAX, okL)
    check("C_D^(0) diagonal of 2-adic units, v_2(det)=0 (D=1..%d)" % DMAX, okC)

    # ---- 4.10 Smith implementation validation -------------------------
    say()
    say("[10] Smith implementation validation")
    Ddiag = [[Fraction(0)] * 3 for _ in range(3)]
    Ddiag[0][0] = Fraction(1); Ddiag[1][1] = Fraction(2); Ddiag[2][2] = Fraction(8)
    # deterministic unimodular matrices over Z_(2) (odd determinant).
    # Built as (unitriangular) x (unitriangular) x diag(odd units) so that
    # membership in GL_3(Z_(2)) is manifest.
    def _F(m):
        return [[Fraction(x) for x in row] for row in m]

    Uu = _F([[1, 3, 5], [0, 1, 7], [0, 0, 1]])
    Ul = _F([[1, 0, 0], [2, 1, 0], [6, 4, 1]])
    Uo = _F([[3, 0, 0], [0, 5, 0], [0, 0, 7]])
    U = matmul(matmul(Uu, Ul), Uo)
    Vu = _F([[1, 0, 0], [5, 1, 0], [7, 9, 1]])
    Vl = _F([[1, 11, 2], [0, 1, 13], [0, 0, 1]])
    Vo = _F([[9, 0, 0], [0, 1, 0], [0, 0, 3]])
    V = matmul(matmul(Vu, Vl), Vo)
    du, dv = det_frac(U), det_frac(V)
    say("      det U = %s (v2=%d), det V = %s (v2=%d)" % (du, v2(du), dv, v2(dv)))
    check("test fixtures U,V really lie in GL_3(Z_(2))", v2(du) == 0 and v2(dv) == 0)
    got = smith_v2(matmul(matmul(U, Ddiag), V))
    check("diag(1,2,8) conjugated by GL_3(Z_(2)) -> [0,1,3]", got == [0, 1, 3], str(got))
    # a rational (negative-valuation) test
    Ddiag2 = [[Fraction(1, 4), 0, 0], [0, Fraction(3, 1), 0], [0, 0, Fraction(16, 5)]]
    Ddiag2 = [[Fraction(x) for x in row] for row in Ddiag2]
    got2 = smith_v2(matmul(matmul(U, Ddiag2), V))
    check("diag(1/4,3,16/5) conjugated -> [-2,0,4]", got2 == [-2, 0, 4], str(got2))
    # rectangular / permutation invariance
    got3 = smith_v2(matmul(matmul(V, Ddiag), U))
    check("order-swapped conjugation -> [0,1,3]", got3 == [0, 1, 3], str(got3))

    say()
    nf = sum(1 for _, o in PASS if not o)
    say(">>> %d checks run, %d failures" % (len(PASS), nf))
    if nf:
        say(">>> WARNING: failures above; downstream numbers are suspect.")

    # ==================================================================
    say()
    say("=" * 78)
    say("PART II.  THE POLYGONS")
    say("=" * 78)

    def H_mat(D):
        return [[(h[j] if i == j else Fraction(0)) for j in range(1, D + 1)]
                for i in range(1, D + 1)]

    def C_adapted(D):
        """C_D = L_D C_D^(0) H_D^{-2}  (E/Z-adapted quotient lattice basis)."""
        LC = matmul(L_mat(D), C0_mat(D))
        return [[LC[i][j] / (h[j + 1] ** 2) for j in range(D)] for i in range(D)]

    results = {}

    # ---------------- Formulation 1 -----------------------------------
    say()
    say("-" * 78)
    say("FORMULATION 1 (source-adapted, quotient only):  N_D = M_D * C_D")
    say("  M_D  : (b_{D+i-j})_{0<=i,j<D}       target functionals [z^{D+i}]")
    say("  C_D  : L_D C_D^(0) H_D^{-2}          E/Z-adapted source lattice on")
    say("         the quotient basis beta, z beta, ..., z^{D-1} beta")
    say("  coordinate charge: everything is in z = Y/64; no 16^k or 64^k")
    say("         diagonal is applied here.  (The Y-coordinate 64-scaling of")
    say("         K5.3 is deliberately NOT charged in this formulation.)")
    say("-" * 78)
    f1 = {}
    for D in range(1, DEXT + 1):
        N_D = matmul(M_block(D), C_adapted(D))
        e = smith_v2(N_D)
        f1[D] = e
    results["F1"] = f1
    report_table("Formulation 1", f1, Omega, DMAX)

    # ---------------- Formulation 2 -----------------------------------
    say()
    say("-" * 78)
    say("FORMULATION 2 (target-adapted, reversed-Pade moving functionals):")
    say("  phi_m(F) = sum_{k=0}^m q_k(1/2-m) a_k     for F = sum a_k t^k")
    say("  psi_m    = Delta_m^{-1} phi_m,   m = 1..D")
    say("  source   = quotient columns z^j beta, j = 0..D-1")
    say("  CONVENTION (see report):  a_k = 16^k c_k where c_k = [z^k] of the")
    say("         column.  Entry(m,j) = Delta_m^{-1} sum_{k<=m} q_k(1/2-m) 16^k b_{k-j}")
    say("-" * 78)

    def F2_mat(D, sgn=+1):
        out = []
        for m in range(1, D + 1):
            qm = qmov[m]
            row = []
            for j in range(0, D):
                acc = Fraction(0)
                for k in range(0, m + 1):
                    if k - j >= 1:
                        acc += qm[k] * Fraction(16) ** (sgn * k) * b[k - j]
                row.append(acc / Delta[m])
            out.append(row)
        return out

    f2 = {}
    for D in range(1, DEXT + 1):
        f2[D] = smith_v2(F2_mat(D, +1))
    results["F2"] = f2
    report_table("Formulation 2 (a_k = 16^k c_k)", f2, Omega, DMAX)

    say()
    say("  [variant] the opposite coordinate reading a_k = 16^{-k} c_k")
    say("            (i.e. literally substituting t = 16z).  Reported for")
    say("            transparency; see report Section 'coordinate trap'.")
    f2b = {}
    for D in range(1, DEXT + 1):
        f2b[D] = smith_v2(F2_mat(D, -1))
    results["F2neg"] = f2b
    report_table("Formulation 2' (a_k = 16^{-k} c_k)", f2b, Omega, DMAX)

    # ---------------- Formulation 3 (both-adapted) --------------------
    say()
    say("-" * 78)
    say("FORMULATION 3 (bonus: source AND target adapted):  psi_D-moving o C_D")
    say("  = F2_mat(D) * C_D")
    say("-" * 78)
    f3 = {}
    for D in range(1, DEXT + 1):
        f3[D] = smith_v2(matmul(F2_mat(D, +1), C_adapted(D)))
    results["F3"] = f3
    report_table("Formulation 3", f3, Omega, DMAX)

    # ---------------- structure hunting -------------------------------
    say()
    say("=" * 78)
    say("PART III.  STRUCTURE SEARCH")
    say("=" * 78)

    for name in ["F1", "F2", "F2neg", "F3"]:
        e = results[name]
        say()
        say("### %s" % name)
        # D -> 2D multiset comparison
        for D in range(1, DMAX // 2 + 1):
            if 2 * D in e:
                a = e[D]; c = e[2 * D]
                doubled = sorted(2 * x for x in a)
                say("  D=%2d -> 2D=%2d :  sorted(e_2D) = %s" % (D, 2 * D, c))
                say("                   2*e_D        = %s" % doubled)
                # E_2D(2k) - 2 E_D(k)
                diffs = []
                for k in range(0, D + 1):
                    E2 = sum(c[:2 * k])
                    E1 = sum(a[:k])
                    diffs.append(E2 - 2 * E1)
                say("       E_2D(2k)-2E_D(k), k=0..D: %s" % diffs)
        # candidate formulas against the sorted list of the largest D
        D = DMAX
        ee = e[D]
        say("  candidate fits at D=%d (exponents sorted ascending, index i=1..D):" % D)
        cands = {
            "-8i+4s2(i)": [-8 * i + 4 * s2(i) for i in range(1, D + 1)],
            "-8i+6s2(i)+1": [-8 * i + 6 * s2(i) + 1 for i in range(1, D + 1)],
            "-8i+2s2(i)": [-8 * i + 2 * s2(i) for i in range(1, D + 1)],
            "-4i": [-4 * i for i in range(1, D + 1)],
            "-8i": [-8 * i for i in range(1, D + 1)],
        }
        for label, vals in cands.items():
            say("     %-16s sorted = %s   match=%s"
                % (label, sorted(vals), sorted(vals) == ee))
        # consecutive gaps
        gaps = [ee[i + 1] - ee[i] for i in range(len(ee) - 1)]
        say("  gaps of sorted exponents at D=%d: %s" % (D, gaps))

    # ---------------- limit shape -------------------------------------
    say()
    say("=" * 78)
    say("PART IV.  NORMALIZED POLYGON  E_D(k)/D^2  vs  x = k/D")
    say("=" * 78)
    for name in ["F1", "F2", "F2neg", "F3"]:
        say()
        say("### %s,  D = %d" % (name, DMAX))
        ee = results[name][DMAX]
        D = DMAX
        say("   k    x=k/D     E_D(k)      E_D(k)/D^2")
        acc = 0
        for k in range(0, D + 1):
            if k:
                acc += ee[k - 1]
            say("  %3d   %.4f   %8d    %+0.6f" % (k, k / D, acc, acc / (D * D)))

    # ---------------- endpoint table ----------------------------------
    say()
    say("=" * 78)
    say("PART V.  ENDPOINT CHECKS")
    say("=" * 78)
    say("  D  Omega_D  D-Omega_D   sum(F1)   sum(F2)   sum(F2')   sum(F3)  S_D")
    for D in range(1, DMAX + 1):
        say("  %2d %8d %10d %9d %9d %10d %9d %4d"
            % (D, Omega[D], D - Omega[D], sum(f1[D]), sum(f2[D]),
               sum(f2b[D]), sum(f3[D]), S(D)))
    say()
    say("  predicted F2 endpoint  -4D^2-3D+6S_D :  %s"
        % [(-4 * D * D - 3 * D + 6 * S(D)) for D in range(1, DMAX + 1)])
    say("  F1 endpoint - F2 endpoint            :  %s"
        % [sum(f1[D]) - sum(f2[D]) for D in range(1, DMAX + 1)])
    say("  2*S_D                                :  %s"
        % [2 * S(D) for D in range(1, DMAX + 1)])

    # ==================================================================
    say()
    say("=" * 78)
    say("PART VI.  RIGIDITY: SMITH POLYGON vs DIAGONAL PROFILE")
    say("=" * 78)

    def closed(fn, D):
        return sorted(fn(j) for j in range(1, D + 1))

    say()
    say("[F1] closed form  e_j = 1 - v_2(Delta_j) = 1 - 8j + 4 s_2(j) ?")
    okf = all(f1[D] == closed(lambda j: 1 - 8 * j + 4 * s2(j), D)
              for D in range(1, DEXT + 1))
    check("F1 exponents == sorted{1-8j+4s_2(j)}  (D=1..%d)" % DEXT, okf)

    say()
    say("[F2'] closed form  e_j = -16j + 6 s_2(j) + 1 ?")
    okf = all(f2b[D] == closed(lambda j: -16 * j + 6 * s2(j) + 1, D)
              for D in range(1, DEXT + 1))
    check("F2' exponents == sorted{-16j+6s_2(j)+1}  (D=1..%d)" % DEXT, okf)

    say()
    say("[F3] closed form  e_j = -16j + 10 s_2(j) + 1 ?")
    okf = all(f3[D] == closed(lambda j: -16 * j + 10 * s2(j) + 1, D)
              for D in range(1, DEXT + 1))
    check("F3 exponents == sorted{-16j+10s_2(j)+1}  (D=1..%d)" % DEXT, okf)

    say()
    say("[F2] diagonal profile of the (triangular) F2 matrix is")
    say("     d_m = -v_2(Delta_m) + v_2(q_m(1/2-m)) + 4m + v_2(b_1)")
    say("         = -8m + 6 s_2(m) + 1.")
    say("     Smith polygon vs sorted diagonal profile:")
    say("     D | k | E_smith(k) | E_diag(k) | defect = E_diag - E_smith")
    maxdef = {}
    for D in range(1, DMAX + 1):
        dg = closed(lambda j: -8 * j + 6 * s2(j) + 1, D)
        es, ed = 0, 0
        defs = []
        for k in range(1, D + 1):
            es += f2[D][k - 1]
            ed += dg[k - 1]
            defs.append(ed - es)
        maxdef[D] = max(defs)
        say("     %2d | defects k=1..D: %s   (max %d)" % (D, defs, max(defs)))
    say()
    say("     F2 exponents == its diagonal profile?  %s"
        % [f2[D] == closed(lambda j: -8 * j + 6 * s2(j) + 1, D)
           for D in range(1, DMAX + 1)])
    say("     max interior defect by D: %s"
        % [maxdef[D] for D in range(1, DMAX + 1)])

    say()
    say("  LIMIT SHAPE.  If e_j = -c j + O(log j) then, sorting descending in j,")
    say("  E_D(k)/D^2 -> -(c/2)(2x - x^2) with x = k/D.  Compare:")
    for name, c in [("F1", 8), ("F2", 8), ("F2neg", 16), ("F3", 16)]:
        ee = results[name][DMAX]
        D = DMAX
        acc = 0
        rows = []
        for k in range(1, D + 1):
            acc += ee[k - 1]
            x = k / D
            pred = -(c / 2) * (2 * x - x * x)
            rows.append("%.3f/%.3f" % (acc / (D * D), pred))
        say("   %-6s obs/pred: %s" % (name, "  ".join(rows)))

    # ==================================================================
    say()
    say("=" * 78)
    say("PART VII.  EXTENDED F2 SCAN,  D = 1..%d" % DEXT)
    say("=" * 78)
    say("  Only Formulation 2 shows off-diagonal mixing (Part VI), so it is")
    say("  tabulated here over the full extended range.")
    say("  Columns: max interior defect (E_diag - E_smith), largest exponent,")
    say("  total positive Smith mass, total negative Smith mass.")
    say()
    say("   D | maxdef | top e | pos mass | neg mass | endpoint | D-Omega_D")
    f2ext = f2
    for D in range(1, DEXT + 1):
        e = f2ext[D]
        dg = sorted(-8 * j + 6 * s2(j) + 1 for j in range(1, D + 1))
        es = ed = 0
        defs = []
        for k in range(D):
            es += e[k]
            ed += dg[k]
            defs.append(ed - es)
        say("  %2d | %6d | %5d | %8d | %8d | %8d | %9d"
            % (D, max(defs), e[-1], sum(x for x in e if x > 0),
               sum(x for x in e if x < 0), sum(e), D - Omega[D]))
    results["F2ext"] = f2ext

    say()
    say("  Observed step positions.")
    steps_def = [D for D in range(2, DEXT + 1)
                 if max_defect(f2ext, D) != max_defect(f2ext, D - 1)]
    steps_top = [D for D in range(2, DEXT + 1)
                 if f2ext[D][-1] != f2ext[D - 1][-1]]
    say("    max-defect steps at D = %s" % steps_def)
    say("    top-exponent steps at D = %s" % steps_top)
    say("    2^k-1 up to %d      : %s"
        % (DEXT, [2 ** k - 1 for k in range(2, 8) if 2 ** k - 1 <= DEXT]))
    say("    3*2^(k-1)-1 up to %d: %s"
        % (DEXT, [3 * 2 ** (k - 1) - 1 for k in range(2, 8)
                  if 3 * 2 ** (k - 1) - 1 <= DEXT]))
    say()
    say("    conjecture  top e_{D,D} = 2*floor(log2(D+1)) - 5   (D>=3):")
    import math as _mm
    say("      obs : %s" % [f2ext[D][-1] for D in range(3, DEXT + 1)])
    say("      pred: %s" % [2 * int(_mm.log2(D + 1)) - 5 for D in range(3, DEXT + 1)])
    say("      match: %s"
        % all(f2ext[D][-1] == 2 * int(_mm.log2(D + 1)) - 5
              for D in range(3, DEXT + 1)))

    return results


def max_defect(f2ext, D):
    e = f2ext[D]
    dg = sorted(-8 * j + 6 * s2(j) + 1 for j in range(1, D + 1))
    es = ed = 0
    m = 0
    for k in range(D):
        es += e[k]
        ed += dg[k]
        m = max(m, ed - es)
    return m


def report_table(name, e, Omega, DMAX):
    say()
    say("  %s -- Smith exponents" % name)
    say("   D | exponents (ascending)")
    for D in range(1, DMAX + 1):
        say("  %2d | %s" % (D, e[D]))
    say()
    say("   D | neg mass | pos mass | zeros | endpoint sum | D-Omega_D")
    for D in range(1, DMAX + 1):
        ee = e[D]
        neg = sum(x for x in ee if x < 0)
        pos = sum(x for x in ee if x > 0)
        z = sum(1 for x in ee if x == 0)
        say("  %2d | %8d | %8d | %5d | %12d | %9d"
            % (D, neg, pos, z, sum(ee), D - Omega[D]))
    say()
    say("   D | partial sums E_D(k), k=1..D")
    for D in range(1, DMAX + 1):
        ee = e[D]
        ps = []
        acc = 0
        for x in ee:
            acc += x
            ps.append(acc)
        say("  %2d | %s" % (D, ps))


if __name__ == "__main__":
    res = main()
