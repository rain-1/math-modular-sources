#!/usr/bin/env python3
"""
plot_catalan_090.py

Six-panel visualization of the 5:8 Catalan construction proving the
0.902526602856971... scalar worthiness bound.

Reads:
    zudilin_5n_rows.csv
    modular_E_8n_rows.csv
    combined_090_rows.csv
    cross_090_diagnostics.csv

Panels
------
1. raw and combined approximation gain;
2. empirical worthiness, theorem value, and irrationality barrier;
3. raw and integer-scaled 2-adic cross slopes;
4. LCM / dyadic / total arithmetic entropy;
5. coefficient-lattice anisotropy versus theoretical Minkowski scales;
6. the scalar optimization curve delta(r), showing its unique maximum at
   r=k/i=8/5.

The theoretical lines are derived from the proved asymptotic formulas; the
solid data are the exact finite CSV outputs.
"""

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_csv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def col(rows, name):
    return np.array([float(r[name]) for r in rows], dtype=float)


def tail_fit(x, y, frac=0.5):
    start = max(0, int((1-frac)*len(x)))
    xx, yy = x[start:], y[start:]
    if len(xx) < 2:
        return None
    a, b = np.polyfit(xx, yy, 1)
    return float(a), float(b), xx


def delta_curve(r):
    phi = (1 + math.sqrt(5))/2
    H = (
        np.maximum(2.0, r)
        + (2.0+r)*math.log(2)
        - 0.5*np.minimum(8.0, 5.0*r)*math.log(2)
        + 7.5*math.log(phi)
    )
    return 10*math.log(phi)/H


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv-dir", type=Path, default=Path("."))
    ap.add_argument("--out", type=Path, default=Path("catalan_090_geometry.png"))
    ap.add_argument("--tail-fraction", type=float, default=0.5)
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()

    z = read_csv(args.csv_dir/"zudilin_5n_rows.csv")
    e = read_csv(args.csv_dir/"modular_E_8n_rows.csv")
    c = read_csv(args.csv_dir/"combined_090_rows.csv")
    d = read_csv(args.csv_dir/"cross_090_diagnostics.csv")

    n = col(c, "n")
    nz = col(z, "n")
    ne = col(e, "n")
    nd = col(d, "n")

    dz = col(z, "digits_correct")
    de = col(e, "digits_correct")
    dc = col(c, "best_digits_correct")
    empirical_delta = col(c, "empirical_delta")
    log10q = col(c, "log10_abs_q")
    lf_digits = col(c, "minus_log10_abs_linear_form")
    theory_delta = float(c[-1]["theory_delta"])

    v2Delta = col(d, "v2_Delta")
    exact_bound = col(d, "exact_digit_bound")
    safe_bound = col(d, "safe_Delta_bound")
    v2h = col(d, "v2_h")
    Texp = col(d, "T_exp")

    log2S = col(d, "log2_S_per_n")
    log2T = col(d, "log2_T_per_n")
    log2M = col(d, "log2_M_per_n")
    log2idxK = col(d, "log2_indexK_per_n")
    log2idxL = col(d, "log2_indexL_per_n")

    c1bits = col(d, "log2_abs_c1_per_n")
    c2bits = col(d, "log2_abs_c2_per_n")

    phi = (1 + math.sqrt(5))/2
    H = 10 - 2*math.log(2) + 37.5*math.log(phi)
    F = 10 - 2*math.log(2) - 12.5*math.log(phi)
    digits_slope_theory = (H-F)/math.log(10)
    height_slope_theory = H/math.log(10)

    kappa = 20 + 40*math.log(2)
    x = 10 + 22*math.log(2) - 12.5*math.log(phi)
    c1_bits_theory = x/math.log(2)
    c2_bits_theory = (kappa-x)/math.log(2)

    fig, axes = plt.subplots(3, 2, figsize=(14, 13))
    ax_gain, ax_delta, ax_v2, ax_entropy, ax_coeff, ax_opt = axes.flat

    # 1. approximation slopes
    ax_gain.plot(nz, dz, marker=".", linewidth=1.2, label="Zudilin 5n")
    ax_gain.plot(ne, de, marker=".", linewidth=1.2, label="modular E 8n")
    ax_gain.plot(n, dc, marker=".", linewidth=1.3, label="combined")
    ax_gain.plot(
        n, digits_slope_theory*n, "--", linewidth=1.2,
        label=f"theory combined slope {digits_slope_theory:.5f} digits/n"
    )
    ax_gain.set_title("Real approximation gain")
    ax_gain.set_xlabel("n")
    ax_gain.set_ylabel(r"$-\log_{10}|G-p/q|$")
    ax_gain.legend()

    # 2. empirical delta
    ax_delta.plot(n, empirical_delta, marker=".", linewidth=1.3,
                  label="experimental lattice selector")
    ax_delta.axhline(theory_delta, linestyle="--", linewidth=1.4,
                     label=rf"theorem $\delta={theory_delta:.9f}$")
    ax_delta.axhline(1.0, linestyle=":", linewidth=1.4,
                     label="irrationality barrier")
    ax_delta.set_title("Worthiness and the scalar ceiling")
    ax_delta.set_xlabel("n")
    ax_delta.set_ylabel(r"$-\log|G-p/q|/\log|q|$")
    ax_delta.legend()

    # 3. 2-adic cross slopes
    ax_v2.plot(nd, v2Delta/nd, marker=".", linewidth=1.2,
               label=r"$v_2(\Delta_{5n,8n})/n$")
    ax_v2.plot(nd, exact_bound/nd, linewidth=1.1,
               label="digit-sum lower bound / n")
    ax_v2.axhline(20.0, linestyle=":", linewidth=1.2,
                  label="raw asymptotic slope 20")
    ax_v2.plot(nd, v2h/nd, marker=".", linewidth=1.2,
               label=r"$v_2(h_n)/n$")
    ax_v2.plot(nd, Texp/nd, "--", linewidth=1.2,
               label=r"safe $T_n$ exponent / n")
    ax_v2.axhline(40.0, linestyle=":", linewidth=1.2,
                  label="integer cross slope 40")
    ax_v2.set_title("The 5:8 dyadic ridge")
    ax_v2.set_xlabel("n")
    ax_v2.set_ylabel("powers of 2 per n")
    ax_v2.legend(fontsize=8)

    # 4. entropy budget
    ax_entropy.plot(nd, log2S, marker=".", linewidth=1.2,
                    label=r"$\log_2 S_n/n$")
    ax_entropy.plot(nd, log2T, marker=".", linewidth=1.2,
                    label=r"$\log_2 T_n/n$")
    ax_entropy.plot(nd, log2M, marker=".", linewidth=1.2,
                    label=r"$\log_2(S_nT_n)/n$")
    ax_entropy.plot(nd, log2idxK, linewidth=1.2,
                    label="raw double-congruence index / n")
    ax_entropy.plot(nd, log2idxL, linewidth=1.2,
                    label="thinned covolume / n")
    ax_entropy.axhline(20/math.log(2), linestyle=":", linewidth=1.0)
    ax_entropy.axhline(40, linestyle=":", linewidth=1.0)
    ax_entropy.set_title("Arithmetic entropy budget")
    ax_entropy.set_xlabel("n")
    ax_entropy.set_ylabel("bits per n")
    ax_entropy.legend(fontsize=8)

    # 5. coefficient anisotropy
    ax_coeff.plot(nd, c1bits, marker=".", linewidth=1.2,
                  label=r"$\log_2|c_1|/n$ (modular coefficient)")
    ax_coeff.plot(nd, c2bits, marker=".", linewidth=1.2,
                  label=r"$\log_2|c_2|/n$ (Zudilin coefficient)")
    ax_coeff.axhline(c1_bits_theory, linestyle="--", linewidth=1.2,
                     label=f"predicted c1 scale {c1_bits_theory:.3f}")
    ax_coeff.axhline(c2_bits_theory, linestyle="--", linewidth=1.2,
                     label=f"predicted c2 scale {c2_bits_theory:.3f}")
    ax_coeff.set_title("Minkowski coefficient anisotropy")
    ax_coeff.set_xlabel("n")
    ax_coeff.set_ylabel("bits per n")
    ax_coeff.legend(fontsize=8)

    # 6. exact scalar optimization
    r = np.linspace(0.45, 3.0, 1200)
    dd = delta_curve(r)
    ax_opt.plot(r, dd, linewidth=1.5)
    ax_opt.axvline(8/5, linestyle="--", linewidth=1.2,
                   label=r"unique optimum $k/i=8/5$")
    ax_opt.axhline(theory_delta, linestyle="--", linewidth=1.2,
                   label=f"maximum {theory_delta:.9f}")
    ax_opt.axhline(1.0, linestyle=":", linewidth=1.2,
                   label="barrier 1")
    ax_opt.scatter([8/5], [theory_delta], zorder=5)
    ax_opt.set_title("Scalar architecture: exact ratio optimization")
    ax_opt.set_xlabel(r"$r=k/i$")
    ax_opt.set_ylabel(r"$\delta(r)$")
    ax_opt.legend()

    # Console tail fits
    print("Tail-fit summary")
    print("="*78)
    for name, yy in [
        ("combined digits/n", dc),
        ("combined log10|q|/n", log10q),
        ("raw v2 Delta slope", v2Delta),
        ("integer v2 h slope", v2h),
    ]:
        fit = tail_fit(n if len(yy)==len(n) else nd, yy, args.tail_fraction)
        if fit:
            a, b, _ = fit
            print(f"{name:28s}: {a:.10f}")

    fitd = tail_fit(n, dc, args.tail_fraction)
    fith = tail_fit(n, log10q, args.tail_fraction)
    if fitd and fith:
        print(f"tail slope ratio (worthiness): {fitd[0]/fith[0]:.12f}")
    print(f"theorem worthiness:              {theory_delta:.12f}")
    print(f"theory digits slope:             {digits_slope_theory:.10f}")
    print(f"theory height slope:             {height_slope_theory:.10f}")

    fig.suptitle(
        "Catalan 0.902526 construction: Zudilin 5n + modular E 8n\n"
        "solid/markers = exact finite data; dashed/dotted = theorem asymptotics",
        fontsize=15,
    )
    fig.tight_layout()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=180, bbox_inches="tight")
    print(f"saved {args.out}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
