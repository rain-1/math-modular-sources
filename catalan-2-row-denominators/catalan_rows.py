"""
Two Apery-like rows for Catalan's constant G = beta(2), their LCM-square
combination, and a study of the p=2 denominator arithmetic.

Uses mpmath.hyper for exact, ultra-fast evaluation of J_n.
Includes instant disk-syncing checkpointing, true zero-recomputation resume, 
and full sequence integer values saved directly in the CSV files.
"""

import csv
import os
import sys
from fractions import Fraction
from math import gcd, log2, floor, factorial
import mpmath as mp

# Expand Python's string conversion limit for huge integers (Python 3.11+)
if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(1_000_000)


def to_mpf(x):
    if isinstance(x, Fraction):
        return mp.mpf(x.numerator) / mp.mpf(x.denominator)
    return mp.mpf(x)


def digits_correct(approx, ref, cap=10000):
    a = to_mpf(approx) if not isinstance(approx, mp.mpf) else approx
    diff = abs(a - ref)
    if diff == 0:
        return cap
    d = float(mp.log10(1 / diff))
    return max(0.0, min(cap, d))


def lcm_table(M):
    D = [1] * (M + 1)
    for m in range(1, M + 1):
        D[m] = D[m - 1] * m // gcd(D[m - 1], m)
    return D


def val2(n):
    if n == 0: return None
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def is_power_of_two(n):
    return n > 0 and (n & (n - 1)) == 0


# =====================================================================
# PART 1 -- Zudilin's row
# =====================================================================
def phi(t):
    return 20 * t * t - 8 * t + 1


def zudilin_recurrence(M):
    Q = [Fraction(1), Fraction(7, 4)] + [None] * (M - 1)
    P = [Fraction(0), Fraction(13, 8)] + [None] * (M - 1)
    for m in range(1, M):
        a = (2 * m + 1) ** 2 * (2 * m + 2) ** 2 * phi(m)
        b = (3520 * m**6 + 5632 * m**5 + 2064 * m**4
             - 384 * m**3 - 156 * m**2 + 16 * m + 7)
        c = (2 * m - 1) ** 2 * (2 * m) ** 2 * phi(m + 1)
        Q[m + 1] = (b * Q[m] + c * Q[m - 1]) / a
        P[m + 1] = (b * P[m] + c * P[m - 1]) / a
    return Q, P


def e_bound(m):
    return min(6 * m, 4 * m + 3 + floor(log2(2 * m - 1)))


def zudilin_row_at_3n(Q, P, D, n):
    m = 3 * n
    e = e_bound(m)
    Dval = D[6 * n - 1]
    Xf = (Fraction(2) ** e) * Dval * Dval * Q[m]
    Yf = (Fraction(2) ** e) * Dval * Dval * P[m]
    return int(Xf), int(Yf)


# =====================================================================
# PART 2 -- Nesterenko's (4,7) row
# =====================================================================
def A2_coeff(n, j):
    num = factorial(8 * n + 2 * j) * factorial(j) * factorial(6 * n)
    den = (factorial(4 * n) * factorial(4 * n + j)
           * factorial(3 * n - j) ** 2 * factorial(2 * j) ** 2)
    f = Fraction(num, den)
    e = -14 * n + 2 * j + 1
    if e >= 0:
        f *= (1 << e)
    else:
        f /= (1 << (-e))
    return f


def B_n(n):
    return sum(A2_coeff(n, j) for j in range(0, 3 * n + 1))


def _t0_series(n):
    return Fraction(4 ** (7 * n + 2) * factorial(8 * n) * factorial(6 * n)
                     * factorial(7 * n + 1) ** 2,
                     factorial(14 * n + 2) ** 2)


def J_n_hyper(n, dps):
    """Computes J_n closed-form via 3F2 hypergeometric series acceleration."""
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        t0f = _t0_series(n)
        t0 = mp.mpf(t0f.numerator) / mp.mpf(t0f.denominator)
        
        a1 = mp.mpf(4 * n + 1)
        a2 = mp.mpf(4 * n + 1)
        a3 = mp.mpf(4 * n) + mp.mpf('0.5')
        
        b1 = mp.mpf(7 * n) + mp.mpf('1.5')
        b2 = mp.mpf(7 * n) + mp.mpf('1.5')
        
        # 3F2 hypergeometric series at z=1
        h = mp.hyper([a1, a2, a3], [b1, b2], 1)
        val = t0 * h
    finally:
        mp.mp.dps = old
    return val


def nesterenko_row(n, D, dps_guard=50):
    Dval = D[6 * n]
    b = B_n(n)
    Vf = Fraction(4) ** (7 * n + 1) * Dval * Dval * b
    V = int(Vf)

    # Calculate exact working precision needed based on digit length of V_n
    v_digits = int(V.bit_length() * 0.30103)
    dps = v_digits + dps_guard

    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        g_ref = mp.catalan
        J = J_n_hyper(n, dps=dps)
        lam = mp.mpf(4) ** (7 * n) * mp.mpf(Dval) ** 2 * J
        Uf = mp.mpf(V) * g_ref - lam
        U = int(mp.nint(Uf))
        residual = abs(Uf - U)
    finally:
        mp.mp.dps = old

    if residual > mp.mpf('1e-5'):
        print(f"    WARNING n={n}: residual {mp.nstr(residual, 3)} high; retrying with +100 dps...")
        dps_ret = dps + 100
        mp.mp.dps = dps_ret
        try:
            g_ref = mp.catalan
            J = J_n_hyper(n, dps=dps_ret)
            lam = mp.mpf(4) ** (7 * n) * mp.mpf(Dval) ** 2 * J
            Uf = mp.mpf(V) * g_ref - lam
            U = int(mp.nint(Uf))
        finally:
            mp.mp.dps = old

    return V, U


# =====================================================================
# PART 3 -- combining the two rows
# =====================================================================
def ext_gcd(a, b):
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r != 0:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    return old_r, old_s, old_t


def congruence_lattice_basis(a, b, m):
    g1, p, q = ext_gcd(a, b)
    aprime, bprime = a // g1, b // g1
    g = gcd(g1, m)
    m2 = m // g
    w1 = (m2 * p, m2 * q)
    w2 = (bprime, -aprime)
    return w1, w2, m2


def gauss_reduce(v1, v2):
    v1, v2 = list(v1), list(v2)
    def dot(x, y):
        return x[0] * y[0] + x[1] * y[1]
    while True:
        if dot(v2, v2) < dot(v1, v1):
            v1, v2 = v2, v1
        n1 = dot(v1, v1)
        if n1 == 0:
            return tuple(v1), tuple(v2)
        mu = Fraction(dot(v1, v2), n1)
        mu_round = round(mu)
        if mu_round == 0:
            return tuple(v1), tuple(v2)
        v2 = [v2[0] - mu_round * v1[0], v2[1] - mu_round * v1[1]]


def combine_rows(X, Y, V, U, S, g_ref):
    w1, w2, d = congruence_lattice_basis(Y, U, S)
    b1, b2 = gauss_reduce(w1, w2)
    candidates = []
    for c in (b1, b2):
        c1, c2 = c
        num_q, num_p = c1 * X + c2 * V, c1 * Y + c2 * U
        q, p = num_q // S, num_p // S
        candidates.append((q, p))

    def lin_form(q, p):
        return to_mpf(q) * g_ref - to_mpf(p)

    (q1, p1), (q2, p2) = candidates
    l1, l2 = lin_form(q1, p1), lin_form(q2, p2)
    if l1 != 0:
        m_n = int(mp.nint(l2 / l1))
        qprime, pprime = q2 - m_n * q1, p2 - m_n * p1
        lprime = lin_form(qprime, pprime)
        if qprime != 0 and lprime != 0:
            candidates.append((qprime, pprime))
    return candidates, d


# =====================================================================
# CSV CHECKPOINTING HELPERS
# =====================================================================
def load_csv_checkpoint(filename):
    rows = {}
    if os.path.exists(filename):
        with open(filename, "r") as f:
            reader = csv.reader(f)
            header = next(reader, None)
            for row in reader:
                if row:
                    n = int(row[0])
                    rows[n] = row
    return rows


# =====================================================================
# MAIN REPORT
# =====================================================================
if __name__ == "__main__":
    import time

    NEST_N_MAX = 98
    
    ZUDILIN_TERMS = 3 * NEST_N_MAX + 10
    MAX_LCM_INDEX = 6 * NEST_N_MAX + 10
    
    DENOM_TABLE_M = [1, 2, 3, 5, 8, 10, 15, 20, 30, 40, 60, 80, 100]

    print("=" * 78)
    print("PART 1 -- Zudilin's row: exact rational recurrence")
    print("=" * 78)
    t0 = time.time()
    Q, P = zudilin_recurrence(ZUDILIN_TERMS)
    print(f"Computed Q_m, P_m exactly for m=0..{ZUDILIN_TERMS} in {time.time()-t0:.1f}s.")
    
    D = lcm_table(MAX_LCM_INDEX)

    # Save Part 1 Denominators CSV
    denom_data = []
    for m in DENOM_TABLE_M:
        if m > ZUDILIN_TERMS: continue
        denom = Q[m].denominator
        denom_data.append([m, is_power_of_two(denom), val2(denom), e_bound(m)])

    with open("zudilin_denominators.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["m", "is_power_of_two", "v2_denom", "e_m_bound"])
        writer.writerows(denom_data)

    # Reference value for Catalan's constant at high precision for output check
    mp.mp.dps = int(NEST_N_MAX * 10) + 100
    G_REF = +mp.catalan

    # Save Part 1 Zudilin Rows CSV
    zudilin_data = []
    for n in range(1, NEST_N_MAX + 1):
        if 3 * n > ZUDILIN_TERMS: continue
        X, Y = zudilin_row_at_3n(Q, P, D, n)
        ratio = to_mpf(Y) / to_mpf(X)
        zudilin_data.append([n, digits_correct(ratio, G_REF), X.bit_length(), str(X), str(Y)])

    with open("zudilin_rows.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["n", "digits_correct", "X_n_bit_length", "X_n", "Y_n"])
        writer.writerows(zudilin_data)

    print("\n" + "=" * 78)
    print("PART 2 -- Nesterenko's (4,7) row (Fast Hypergeometric Evaluation)")
    print("=" * 78)
    print(f"{'n':>4} {'digits correct':>16} {'V_n bit-length':>16} {'time (s)':>10} {'status':>10}")

    nest_file = "nesterenko_rows.csv"
    existing_nest = load_csv_checkpoint(nest_file)
    
    if not os.path.exists(nest_file) or len(existing_nest) == 0:
        with open(nest_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["n", "digits_correct", "V_n_bit_length", "time_s", "V_n", "U_n"])

    nest_rows = {}
    with open(nest_file, "a", newline="") as f_out:
        writer = csv.writer(f_out)
        
        for n in range(1, NEST_N_MAX + 1):
            if n in existing_nest and len(existing_nest[n]) >= 6:
                row = existing_nest[n]
                dc, v_bits, dt = float(row[1]), int(row[2]), float(row[3])
                V, U = int(row[4]), int(row[5])
                nest_rows[n] = (V, U)
                print(f"{n:4d} {dc:16.2f} {v_bits:16d} {dt:10.1f} {'[resumed]':>10}")
            else:
                t0 = time.time()
                V, U = nesterenko_row(n, D)
                dt = time.time() - t0
                nest_rows[n] = (V, U)
                ratio = to_mpf(U) / to_mpf(V)
                dc = digits_correct(ratio, G_REF)
                v_bits = V.bit_length()
                
                writer.writerow([n, dc, v_bits, round(dt, 2), str(V), str(U)])
                f_out.flush()
                os.fsync(f_out.fileno())
                print(f"{n:4d} {dc:16.2f} {v_bits:16d} {dt:10.2f} {'[computed]':>10}")

    print("\n" + "=" * 78)
    print("PART 3 -- combined row via LCM-square lattice")
    print("=" * 78)
    print(f"{'n':>4} {'best digits correct':>20} {'best |q_n| bit-length':>22} {'status':>10}")

    comb_file = "combined_rows.csv"
    existing_comb = load_csv_checkpoint(comb_file)

    if not os.path.exists(comb_file) or len(existing_comb) == 0:
        with open(comb_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["n", "best_digits_correct", "best_q_n_bit_length", "q_n", "p_n"])

    with open(comb_file, "a", newline="") as f_out:
        writer = csv.writer(f_out)
        
        for n in range(1, NEST_N_MAX + 1):
            if n in existing_comb and len(existing_comb[n]) >= 5:
                row = existing_comb[n]
                best_dg, best_bits = float(row[1]), int(row[2])
                print(f"{n:4d} {best_dg:20.2f} {best_bits:22d} {'[resumed]':>10}")
            else:
                X, Y = zudilin_row_at_3n(Q, P, D, n)
                V, U = nest_rows[n]
                S = D[6 * n] ** 2
                cands, d_n = combine_rows(X, Y, V, U, S, G_REF)
                best = None
                for q, p in cands:
                    if q == 0: continue
                    ratio = to_mpf(p) / to_mpf(q)
                    dg = digits_correct(ratio, G_REF)
                    if best is None or dg > best[0]:
                        best = (dg, q, p)
                best_bits = best[1].bit_length()
                best_q, best_p = best[1], best[2]
                
                writer.writerow([n, best[0], best_bits, str(best_q), str(best_p)])
                f_out.flush()
                os.fsync(f_out.fileno())
                print(f"{n:4d} {best[0]:20.2f} {best_bits:22d} {'[computed]':>10}")

    print("\nDone. All computations saved/resumed seamlessly.")