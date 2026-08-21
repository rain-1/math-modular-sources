#!/usr/bin/env python3
"""
plot_catalan_arithmetic_structure.py

Arithmetic diagnostics for the Catalan CSV files emitted by
catalan_rows_fast.py.

Expected files in --csv-dir:
    zudilin_rows.csv
    nesterenko_rows.csv
    combined_rows.csv
    zudilin_denominators.csv

The script reconstructs the scaled cross determinant

    H_n = X_n U_n - Y_n V_n,

and displays:

1. v_2(H_n)/n, with a tail-fitted asymptotic slope;
2. total versus 2-adic versus odd-part logarithmic mass of H_n;
3. LCM-square mass, reconstructed congruence-lattice index, and the common-gcd saving;
4. the exact Zudilin denominator v_2 and the proved e_m clearing bound.

It also writes a derived CSV with the determinant/lattice diagnostics.

This script deliberately does not hard-code a theorem slope.  It lets the
finite data reveal the slope.  Add comparisons later only after matching the
normalization used by the theorem.
"""

import argparse
import csv
import math
import sys
from math import gcd
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(1_000_000)


def read_csv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def v2_integer(n):
    n = abs(int(n))
    if n == 0:
        return None
    return (n & -n).bit_length() - 1


def lcm_table(M):
    D = [1] * (M + 1)
    for m in range(1, M + 1):
        D[m] = D[m - 1] * m // gcd(D[m - 1], m)
    return D


def log2_integer(n):
    n = abs(int(n))
    if n == 0:
        return float("-inf")
    # math.log2 handles arbitrary Python ints without requiring decimal conversion.
    return math.log2(n)


def tail_slope(x, y, tail_fraction):
    if len(x) < 4:
        return None
    start = max(0, int((1.0 - tail_fraction) * len(x)))
    xx = np.asarray(x[start:], dtype=float)
    yy = np.asarray(y[start:], dtype=float)
    slope, intercept = np.polyfit(xx, yy, 1)
    return float(slope), float(intercept), xx


def load_by_n(path):
    return {int(r["n"]): r for r in read_csv(path)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv-dir", type=Path, default=Path("."))
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("catalan_arithmetic_structure.png"),
    )
    parser.add_argument(
        "--derived-csv",
        type=Path,
        default=Path("catalan_arithmetic_diagnostics.csv"),
    )
    parser.add_argument("--tail-fraction", type=float, default=0.5)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    required = [
        "zudilin_rows.csv",
        "nesterenko_rows.csv",
        "combined_rows.csv",
        "zudilin_denominators.csv",
    ]
    for name in required:
        if not (args.csv_dir / name).exists():
            raise SystemExit(f"missing required file: {args.csv_dir / name}")

    z = load_by_n(args.csv_dir / "zudilin_rows.csv")
    ne = load_by_n(args.csv_dir / "nesterenko_rows.csv")
    co = load_by_n(args.csv_dir / "combined_rows.csv")
    zd = read_csv(args.csv_dir / "zudilin_denominators.csv")

    common_n = sorted(set(z) & set(ne) & set(co))
    if not common_n:
        raise SystemExit("No common n values across row CSV files.")

    D = lcm_table(6 * max(common_n))

    derived = []
    for n in common_n:
        X = int(z[n]["X_n"])
        Y = int(z[n]["Y_n"])
        V = int(ne[n]["V_n"])
        U = int(ne[n]["U_n"])

        H = X * U - Y * V
        if H == 0:
            continue

        vh = v2_integer(H)
        odd = abs(H) >> vh

        S = D[6 * n] ** 2

        # This is the exact index returned by congruence_lattice_basis():
        #
        #   g1 = gcd(Y,U),
        #   d  = S / gcd(g1,S)
        #      = S / gcd(Y,U,S).
        #
        # Reconstruct it from the raw row data rather than depending on a
        # particular combined_rows.csv schema.
        common_gcd = gcd(gcd(abs(Y), abs(U)), S)
        lattice_index_reconstructed = S // common_gcd

        # Older/newer generator revisions may or may not write the index.
        # If a recognizable field is present, use it only as an audit check.
        csv_index = None
        csv_index_field = None
        for candidate in (
            "lattice_index_d_n",
            "lattice_index",
            "d_n",
            "index_d_n",
        ):
            if candidate in co[n] and str(co[n][candidate]).strip() != "":
                csv_index = int(co[n][candidate])
                csv_index_field = candidate
                break

        derived.append(
            {
                "n": n,
                "v2_cross": vh,
                "log2_abs_cross": log2_integer(H),
                "log2_odd_part": log2_integer(odd),
                "log2_S": log2_integer(S),
                "log2_lattice_index": log2_integer(lattice_index_reconstructed),
                "log2_common_gcd_gain": log2_integer(common_gcd),
                "csv_index_field": csv_index_field or "",
                "csv_index_matches": (
                    "" if csv_index is None
                    else int(csv_index == lattice_index_reconstructed)
                ),
            }
        )

    if not derived:
        raise SystemExit("All reconstructed cross determinants vanished.")

    args.derived_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(derived[0].keys())
    with args.derived_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(derived)

    n = np.array([r["n"] for r in derived], dtype=float)
    v2h = np.array([r["v2_cross"] for r in derived], dtype=float)
    log2h = np.array([r["log2_abs_cross"] for r in derived], dtype=float)
    log2odd = np.array([r["log2_odd_part"] for r in derived], dtype=float)
    log2S = np.array([r["log2_S"] for r in derived], dtype=float)
    log2idx = np.array(
        [r["log2_lattice_index"] for r in derived], dtype=float
    )
    log2g = np.array([r["log2_common_gcd_gain"] for r in derived], dtype=float)

    fit = tail_slope(n, v2h, args.tail_fraction)

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    ax_v2, ax_mass, ax_lattice, ax_zden = axes.flat

    ax_v2.plot(n, v2h / n, marker=".", linewidth=1.3, label=r"$v_2(H_n)/n$")
    if fit:
        slope, intercept, xx = fit
        ax_v2.axhline(
            slope, linestyle="--", linewidth=1.2,
            label=f"tail slope of v2(H): {slope:.6g}"
        )
        print(f"tail fit: v2(H_n) ~ {slope:.10f} n + {intercept:.6f}")
    ax_v2.set_title("2-adic reconciliation slope")
    ax_v2.set_xlabel("n")
    ax_v2.set_ylabel(r"$v_2(H_n)/n$")
    ax_v2.legend()

    ax_mass.plot(n, log2h / n, marker=".", linewidth=1.3,
                 label=r"$\log_2|H_n|/n$")
    ax_mass.plot(n, v2h / n, marker=".", linewidth=1.3,
                 label=r"$v_2(H_n)/n$")
    ax_mass.plot(n, log2odd / n, marker=".", linewidth=1.3,
                 label=r"$\log_2|H_n^{\rm odd}|/n$")
    ax_mass.set_title("Cross determinant: total, dyadic, and odd mass")
    ax_mass.set_xlabel("n")
    ax_mass.set_ylabel("bits per n")
    ax_mass.legend()

    ax_lattice.plot(n, log2S / n, marker=".", linewidth=1.3,
                    label=r"$\log_2 D_{6n}^2/n$")
    ax_lattice.plot(n, log2idx / n, marker=".", linewidth=1.3,
                    label="congruence-lattice index / n")
    ax_lattice.plot(n, log2g / n, marker=".", linewidth=1.3,
                    label="common-gcd saving / n")
    ax_lattice.set_title("LCM-square entropy and lattice index")
    ax_lattice.set_xlabel("n")
    ax_lattice.set_ylabel("bits per n")
    ax_lattice.legend()

    m = np.array([int(r["m"]) for r in zd], dtype=float)
    v2den = np.array([int(r["v2_denom"]) for r in zd], dtype=float)
    ebound = np.array([int(r["e_m_bound"]) for r in zd], dtype=float)

    ax_zden.plot(m, v2den, linewidth=1.2, label=r"$v_2(\mathrm{den}\,Q_m)$")
    ax_zden.plot(m, ebound, linewidth=1.2, label=r"$e_m$ clearing bound")
    ax_zden.plot(m, ebound - v2den, linewidth=1.2,
                 label="clearing margin")
    ax_zden.set_title("Zudilin denominator arithmetic")
    ax_zden.set_xlabel("m")
    ax_zden.set_ylabel("power of 2")
    ax_zden.legend()

    audited = [
        r for r in derived if r["csv_index_matches"] != ""
    ]
    mismatches = sum(
        1 for r in audited if r["csv_index_matches"] == 0
    )
    audit_note = (
        f"CSV index audit: {len(audited)} rows, {mismatches} mismatches"
        if audited
        else "lattice index reconstructed exactly from X,Y,U,V data"
    )
    fig.suptitle(
        "Catalan arithmetic reconciliation\n"
        rf"$H_n=X_nU_n-Y_nV_n$; " + audit_note,
        fontsize=14,
    )
    fig.tight_layout()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=180, bbox_inches="tight")

    print(f"saved {args.out}")
    print(f"wrote {args.derived_csv}")
    if audited:
        field_names = sorted(
            {r["csv_index_field"] for r in audited if r["csv_index_field"]}
        )
        print(
            "optional CSV index audit using field(s) "
            + ", ".join(field_names)
            + f": {mismatches} mismatches in {len(audited)} rows"
        )
    else:
        print(
            "combined_rows.csv has no lattice-index column; "
            "this is fine.  The index was reconstructed exactly as "
            "S/gcd(Y,U,S)."
        )

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
