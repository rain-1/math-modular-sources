#!/usr/bin/env python3
"""
catalan_rows_090_fast.py

Experimental generator for the 0.902526602856971... Catalan construction.

Mathematical normalization
--------------------------
Zudilin row:          index m = 5n
modular E row:        index k = 8n
common LCM square:    S_n = D_(10n)^2
safe 2-adic modulus:  T_n = 2^(40n - 3 - 2 floor(log2(5n)))
total division scale: M_n = S_n T_n

Rows:
    X1 = S A_(8n)
    Y1 = 2 S B_(8n)

    X2 = 2^e S Q_(5n)
    Y2 = 2^e S P_(5n)

where
    e = min(6m, 4m+3+floor(log2(2m-1))), m=5n.

The coefficient lattice is exactly the double-congruence lattice

    a1*c1 + a2*c2 == 0 mod T
    Y1*c1 + Y2*c2 == 0 mod S*T

with
    a1=A_(8n), a2=2^e Q_(5n).

For c in this lattice,
    q=(c1 X1+c2 X2)/(S T),
    p=(c1 Y1+c2 Y2)/(S T)
are integers.

The proof only asserts existence via two successive minima.  For experiment,
this script constructs the exact lattice, thins it to covolume ~ S*T, performs
a weighted 2D Gauss reduction using the theorem's coefficient anisotropy, and
then tests a small collection of reduced vectors.  Thus the CSV gives genuine
integer pairs from the proof lattice, but the finite-n "best" choice is an
experimental selector, not part of the theorem.

Outputs
-------
    zudilin_5n_rows.csv
    modular_E_8n_rows.csv
    combined_090_rows.csv
    cross_090_diagnostics.csv
    zudilin_denominators_090.csv
"""

import argparse
import csv
import math
import sys
from fractions import Fraction
from math import gcd, floor, log2
from pathlib import Path

import mpmath as mp

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(1_000_000)


# ---------------------------------------------------------------------
# basic helpers
# ---------------------------------------------------------------------
def to_mpf(x):
    if isinstance(x, Fraction):
        return mp.mpf(x.numerator) / mp.mpf(x.denominator)
    return mp.mpf(x)


def digits_correct(approx, ref, cap=100000):
    diff = abs(to_mpf(approx) - ref)
    if diff == 0:
        return float(cap)
    return max(0.0, min(float(cap), float(-mp.log10(diff))))


def linear_form_digits(q, p, ref):
    ell = abs(to_mpf(q) * ref - to_mpf(p))
    if ell == 0:
        return float("inf")
    return float(-mp.log10(ell))


def v2_int(n):
    n = abs(int(n))
    if n == 0:
        return None
    return (n & -n).bit_length() - 1


def v2_fraction(x):
    x = Fraction(x)
    if x == 0:
        return None
    return v2_int(x.numerator) - v2_int(x.denominator)


def s2(n):
    return int(n).bit_count()


def lcm_table(M):
    D = [1] * (M + 1)
    for m in range(1, M + 1):
        D[m] = D[m - 1] * m // gcd(D[m - 1], m)
    return D


def e_bound(m):
    return min(6 * m, 4 * m + 3 + floor(log2(2 * m - 1)))


def ceil_div(a, b):
    return -(-a // b)


# ---------------------------------------------------------------------
# Zudilin recurrence
# ---------------------------------------------------------------------
def phi_z(m):
    return 20 * m * m - 8 * m + 1


def zudilin_recurrence(M):
    Q = [Fraction(0) for _ in range(M + 1)]
    P = [Fraction(0) for _ in range(M + 1)]
    Q[0] = Fraction(1)
    P[0] = Fraction(0)
    if M >= 1:
        Q[1] = Fraction(7, 4)
        P[1] = Fraction(13, 8)

    for m in range(1, M):
        L = (2*m+1)**2 * (2*m+2)**2 * phi_z(m)
        C = (
            3520*m**6 + 5632*m**5 + 2064*m**4
            - 384*m**3 - 156*m**2 + 16*m + 7
        )
        R = (2*m-1)**2 * (2*m)**2 * phi_z(m+1)
        Q[m+1] = (C*Q[m] + R*Q[m-1]) / L
        P[m+1] = (C*P[m] + R*P[m-1]) / L

    return Q, P


# ---------------------------------------------------------------------
# modular E recurrence
# ---------------------------------------------------------------------
def modular_E_recurrence(N):
    A = [Fraction(0) for _ in range(N + 1)]
    B = [Fraction(0) for _ in range(N + 1)]

    A[0], B[0] = Fraction(1), Fraction(0)
    if N >= 1:
        A[1], B[1] = Fraction(4), Fraction(1)

    for n in range(1, N):
        c = 12*n*(n+1) + 4
        den = (n+1)**2
        A[n+1] = (c*A[n] - 32*n*n*A[n-1]) / den
        B[n+1] = (c*B[n] - 32*n*n*B[n-1]) / den

    return A, B


# ---------------------------------------------------------------------
# exact congruence lattices
# ---------------------------------------------------------------------
def ext_gcd(a, b):
    """Return positive gcd g and x,y with x*a+y*b=g."""
    old_r, r = int(a), int(b)
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q*r
        old_s, s = s, old_s - q*s
        old_t, t = t, old_t - q*t
    if old_r < 0:
        old_r, old_s, old_t = -old_r, -old_s, -old_t
    return old_r, old_s, old_t


def congruence_lattice_basis(a, b, modulus):
    """
    Basis for {(x,y) in Z^2 : a*x+b*y == 0 mod modulus}.

    Returns basis vectors and the exact index.
    """
    g1, p, q = ext_gcd(a, b)
    if g1 == 0:
        return (1, 0), (0, 1), 1

    aprime, bprime = a // g1, b // g1
    g = gcd(g1, modulus)
    idx = modulus // g

    w1 = (idx*p, idx*q)
    w2 = (bprime, -aprime)
    return w1, w2, idx


def lincomb(v1, z1, v2, z2):
    return (z1*v1[0] + z2*v2[0], z1*v1[1] + z2*v2[1])


def det2(v1, v2):
    return v1[0]*v2[1] - v1[1]*v2[0]


def double_congruence_basis(a1, a2, y1, y2, S, T):
    """
    Exact basis for
        a1*c1+a2*c2 == 0 mod T,
        y1*c1+y2*c2 == 0 mod S*T.
    """
    w1, w2, idx1 = congruence_lattice_basis(a1, a2, T)

    b1 = y1*w1[0] + y2*w1[1]
    b2 = y1*w2[0] + y2*w2[1]

    z1, z2, idx2 = congruence_lattice_basis(b1, b2, S*T)

    k1 = lincomb(w1, z1[0], w2, z1[1])
    k2 = lincomb(w1, z2[0], w2, z2[1])

    idx = abs(det2(k1, k2))
    if idx != idx1 * idx2:
        raise ArithmeticError(
            f"double-congruence index mismatch: det={idx}, product={idx1*idx2}"
        )
    return k1, k2, idx


# ---------------------------------------------------------------------
# weighted reduction / experimental selection
# ---------------------------------------------------------------------
def weighted_gauss_reduce(v1, v2, weight_bits):
    """
    Gauss reduction for norm
        (2^weight_bits * c1)^2 + c2^2.

    At the 5:8 optimum, coefficient c1 is the modular-E coefficient and must
    be exponentially smaller than c2.  The theorem predicts
        B/A = exp((E1-E2)n),
    so weighting c1 by this factor converts the anisotropic Minkowski
    rectangle into an approximately round one.
    """
    W2 = 1 << (2 * max(0, int(weight_bits)))

    def dot(x, y):
        return W2*x[0]*y[0] + x[1]*y[1]

    v1, v2 = list(v1), list(v2)
    for _ in range(10000):
        if dot(v2, v2) < dot(v1, v1):
            v1, v2 = v2, v1
        n1 = dot(v1, v1)
        if n1 == 0:
            return tuple(v1), tuple(v2)
        mu = Fraction(dot(v1, v2), n1)
        m = round(mu)
        if m == 0:
            return tuple(v1), tuple(v2)
        v2[0] -= m*v1[0]
        v2[1] -= m*v1[1]
    raise RuntimeError("weighted Gauss reduction did not terminate")


def output_pair(c, X1, Y1, X2, Y2, M):
    c1, c2 = c
    nq = c1*X1 + c2*X2
    np_ = c1*Y1 + c2*Y2
    if nq % M or np_ % M:
        raise ArithmeticError("candidate does not give integral divided outputs")
    return nq // M, np_ // M


def add_linear_form_reductions(v1, v2, X1, Y1, X2, Y2, M, G):
    """
    Return lattice coefficient vectors:
      basis vectors plus small nearest-linear-form reductions in both directions.
    """
    coeffs = [tuple(v1), tuple(v2)]

    def ell(v):
        q, p = output_pair(v, X1, Y1, X2, Y2, M)
        return to_mpf(q)*G - to_mpf(p)

    l1, l2 = ell(v1), ell(v2)

    if l1 != 0:
        m = int(mp.nint(l2/l1))
        for j in range(-3, 4):
            coeffs.append((v2[0]-(m+j)*v1[0], v2[1]-(m+j)*v1[1]))

    if l2 != 0:
        m = int(mp.nint(l1/l2))
        for j in range(-3, 4):
            coeffs.append((v1[0]-(m+j)*v2[0], v1[1]-(m+j)*v2[1]))

    # Remove duplicate vectors.
    out = []
    seen = set()
    for c in coeffs:
        if c != (0, 0) and c not in seen:
            seen.add(c)
            out.append(c)
    return out


def candidate_stats(c, X1, Y1, X2, Y2, M, G):
    q, p = output_pair(c, X1, Y1, X2, Y2, M)
    if q == 0:
        return None

    ratio = to_mpf(p) / to_mpf(q)
    digits = digits_correct(ratio, G)
    log10q = float(mp.log10(abs(to_mpf(q)))) if abs(q) > 1 else 0.0
    delta = digits/log10q if log10q > 0 else float("-inf")
    lfd = linear_form_digits(q, p, G)

    return {
        "c1": c[0],
        "c2": c[1],
        "q": q,
        "p": p,
        "digits": digits,
        "log10q": log10q,
        "delta": delta,
        "linear_form_digits": lfd,
    }


# ---------------------------------------------------------------------
# main
# ---------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-max", type=int, default=100)
    parser.add_argument("--out-dir", type=Path, default=Path("."))
    parser.add_argument(
        "--dps",
        type=int,
        default=None,
        help="mpmath decimal precision; default 16*n_max+400",
    )
    args = parser.parse_args()

    N = args.n_max
    outdir = args.out_dir
    outdir.mkdir(parents=True, exist_ok=True)

    mp.mp.dps = args.dps or (16*N + 400)
    G = +mp.catalan

    print("Generating exact recurrences...")
    Q, P = zudilin_recurrence(5*N + 2)
    A, B = modular_E_recurrence(8*N + 2)
    D = lcm_table(10*N + 2)

    # Asymptotic constants used only for weighted experimental reduction.
    golden = (1 + mp.sqrt(5))/2
    E1 = 20 + 16*mp.log(2)
    E2 = 20 + 20*mp.log(2) - 25*mp.log(golden)
    weight_bits_per_n = float((E1-E2)/mp.log(2))

    theory_delta = float(
        50*mp.log(golden)
        / (10 - 2*mp.log(2) + mp.mpf(75)/2*mp.log(golden))
    )

    z_rows = []
    e_rows = []
    c_rows = []
    d_rows = []
    den_rows = []

    for m in range(1, 5*N + 1):
        den_rows.append([
            m,
            v2_int(Q[m].denominator),
            e_bound(m),
            Q[m].denominator.bit_length(),
            Q[m].numerator.bit_length(),
        ])

    print(f"theory delta = {theory_delta:.15f}")
    print(f"weighted reduction anisotropy = {weight_bits_per_n:.8f} bits/n")

    for n in range(1, N+1):
        m = 5*n
        k = 8*n
        e = e_bound(m)
        S = D[10*n]**2

        # modular E row
        X1f = S*A[k]
        Y1f = 2*S*B[k]
        if X1f.denominator != 1 or Y1f.denominator != 1:
            raise ArithmeticError(f"E-row integrality failed at n={n}")
        X1, Y1 = int(X1f), int(Y1f)

        # Zudilin row
        X2f = (1 << e)*S*Q[m]
        Y2f = (1 << e)*S*P[m]
        if X2f.denominator != 1 or Y2f.denominator != 1:
            raise ArithmeticError(f"Zudilin-row integrality failed at n={n}")
        X2, Y2 = int(X2f), int(Y2f)

        a1 = int(A[k])
        a2f = (1 << e)*Q[m]
        if a2f.denominator != 1:
            raise ArithmeticError(f"a2 integrality failed at n={n}")
        a2 = int(a2f)

        # Cross determinant before and after integer scaling.
        Delta = P[m]*A[k] - 2*Q[m]*B[k]
        v2_delta = v2_fraction(Delta)
        exact_digit_bound = 20*n - 1 - 2*abs(s2(5*n)-s2(n))
        safe_delta_bound = 20*n - 3 - 2*floor(log2(5*n))

        h = a1*Y2 - a2*Y1
        v2_h = v2_int(h)

        T_exp = 40*n - 3 - 2*floor(log2(5*n))
        T = 1 << T_exp
        M = S*T

        if h % T:
            raise ArithmeticError(
                f"safe cross modulus failed at n={n}: v2(h)={v2_h}, Texp={T_exp}"
            )

        # Exact double-congruence lattice.
        k1, k2, indexK = double_congruence_basis(a1, a2, Y1, Y2, S, T)
        if indexK > M:
            raise ArithmeticError(
                f"index theorem failed at n={n}: indexK > S*T"
            )

        # Thinning.  Try both ways of scaling a basis vector; both are valid.
        rthin = ceil_div(M, indexK)
        thinnings = [
            ((rthin*k1[0], rthin*k1[1]), k2),
            (k1, (rthin*k2[0], rthin*k2[1])),
        ]
        indexL = indexK*rthin
        if not (M <= indexL < 2*M):
            raise ArithmeticError(f"thinning bound failed at n={n}")

        weight_bits = round(weight_bits_per_n*n)

        all_stats = []
        for t1, t2 in thinnings:
            b1, b2 = weighted_gauss_reduce(t1, t2, weight_bits)
            for c in add_linear_form_reductions(
                b1, b2, X1, Y1, X2, Y2, M, G
            ):
                st = candidate_stats(c, X1, Y1, X2, Y2, M, G)
                if st is not None and math.isfinite(st["delta"]):
                    all_stats.append(st)

        if not all_stats:
            raise ArithmeticError(f"no candidate at n={n}")

        # Experimental finite-n selector.
        # To avoid tiny-height artifacts, first retain candidates with at least
        # half the predicted height scale; among those maximize empirical delta.
        phi = (1 + math.sqrt(5))/2
        H_theory = (
            10 - 2*math.log(2) + 37.5*math.log(phi)
        )
        predicted_log10q = H_theory*n/math.log(10)
        eligible = [
            s for s in all_stats if s["log10q"] >= 0.5*predicted_log10q
        ]
        if not eligible:
            eligible = all_stats
        best = max(eligible, key=lambda s: (s["delta"], s["log10q"]))

        # Raw-row diagnostics.
        rz = to_mpf(P[m]) / to_mpf(Q[m])
        re = 2*to_mpf(B[k]) / to_mpf(A[k])

        z_rows.append([
            n, m, digits_correct(rz, G),
            X2.bit_length(), Y2.bit_length(), str(X2), str(Y2),
        ])
        e_rows.append([
            n, k, digits_correct(re, G),
            X1.bit_length(), Y1.bit_length(), str(X1), str(Y1),
        ])

        c_rows.append([
            n,
            best["digits"],
            best["log10q"],
            best["delta"],
            best["linear_form_digits"],
            best["q"].bit_length(),
            best["c1"].bit_length() if best["c1"] else 0,
            best["c2"].bit_length() if best["c2"] else 0,
            indexK.bit_length(),
            indexL.bit_length(),
            T_exp,
            theory_delta,
            str(best["c1"]),
            str(best["c2"]),
            str(best["q"]),
            str(best["p"]),
        ])

        d_rows.append([
            n,
            v2_delta,
            exact_digit_bound,
            safe_delta_bound,
            v2_h,
            T_exp,
            math.log2(abs(h))/n,
            math.log2(S)/n,
            math.log2(T)/n,
            math.log2(M)/n,
            math.log2(indexK)/n,
            math.log2(indexL)/n,
            math.log2(abs(best["c1"]))/n if best["c1"] else float("-inf"),
            math.log2(abs(best["c2"]))/n if best["c2"] else float("-inf"),
        ])

        if n <= 5 or n % 10 == 0:
            print(
                f"n={n:3d}  v2Delta={v2_delta:5d}  v2h={v2_h:5d}  "
                f"T={T_exp:5d}  digits={best['digits']:8.2f}  "
                f"delta={best['delta']:.8f}"
            )

    def write(name, header, rows):
        path = outdir / name
        with path.open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(header)
            w.writerows(rows)
        print(f"wrote {path}")

    write(
        "zudilin_5n_rows.csv",
        ["n","m","digits_correct","X2_bitlen","Y2_bitlen","X2","Y2"],
        z_rows,
    )
    write(
        "modular_E_8n_rows.csv",
        ["n","k","digits_correct","X1_bitlen","Y1_bitlen","X1","Y1"],
        e_rows,
    )
    write(
        "combined_090_rows.csv",
        [
            "n","best_digits_correct","log10_abs_q","empirical_delta",
            "minus_log10_abs_linear_form","q_bitlen","c1_bitlen","c2_bitlen",
            "indexK_bitlen","indexL_bitlen","T_exp","theory_delta",
            "c1","c2","q_n","p_n",
        ],
        c_rows,
    )
    write(
        "cross_090_diagnostics.csv",
        [
            "n","v2_Delta","exact_digit_bound","safe_Delta_bound",
            "v2_h","T_exp","log2_abs_h_per_n","log2_S_per_n",
            "log2_T_per_n","log2_M_per_n","log2_indexK_per_n",
            "log2_indexL_per_n","log2_abs_c1_per_n","log2_abs_c2_per_n",
        ],
        d_rows,
    )
    write(
        "zudilin_denominators_090.csv",
        ["m","v2_denom_Qm","e_m_bound","denom_bitlen","numerator_bitlen"],
        den_rows,
    )

    print("Done.")


if __name__ == "__main__":
    main()
