#!/usr/bin/env python3
"""
plot_catalan_approximation_geometry.py

Visualize the approximation geometry of the Catalan CSV files emitted by
catalan_rows_fast.py.

Expected files in --csv-dir:
    zudilin_rows.csv
    nesterenko_rows.csv
    combined_rows.csv

The four panels are chosen to match the Diophantine mathematics:

1. digits correct versus n;
2. denominator height log_10 |q_n|;
3. empirical worthiness
       delta_n = -log|G-p_n/q_n| / log|q_n|;
4. linear-form scale
       -log_10 |q_n G - p_n|
       = digits_correct - log_10 |q_n|.

Panel 3 includes the barrier delta=1.  A theorem line may optionally be added
with --theory-delta, e.g.

    --theory-delta 0.85791445247

Do not add such a line unless the CSV normalization actually realizes that
theorem; the uploaded generator itself combines the rows using S=D_(6n)^2.
"""

import argparse
import csv
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(1_000_000)


def log10_abs_decimal_string(s):
    """Accurate log10(|integer represented by s|) without converting all digits."""
    s = str(s).strip()
    if s.startswith("-"):
        s = s[1:]
    s = s.lstrip("0")
    if not s:
        return float("-inf")
    k = min(16, len(s))
    lead = int(s[:k])
    return math.log10(lead) + (len(s) - k)


def read_csv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def load_series(csv_dir):
    specs = [
        (
            "Zudilin 3n",
            csv_dir / "zudilin_rows.csv",
            "digits_correct",
            "X_n",
        ),
        (
            "Nesterenko (4,7)",
            csv_dir / "nesterenko_rows.csv",
            "digits_correct",
            "V_n",
        ),
        (
            "Combined",
            csv_dir / "combined_rows.csv",
            "best_digits_correct",
            "q_n",
        ),
    ]

    series = []
    for name, path, digit_col, q_col in specs:
        if not path.exists():
            print(f"warning: missing {path}; skipping {name}", file=sys.stderr)
            continue

        rows = read_csv(path)
        n = np.array([int(r["n"]) for r in rows], dtype=float)
        digits = np.array([float(r[digit_col]) for r in rows], dtype=float)
        log10q = np.array(
            [log10_abs_decimal_string(r[q_col]) for r in rows], dtype=float
        )

        mask = np.isfinite(log10q) & (log10q > 0)
        n = n[mask]
        digits = digits[mask]
        log10q = log10q[mask]

        delta = digits / log10q
        minus_log10_linear_form = digits - log10q

        series.append(
            {
                "name": name,
                "n": n,
                "digits": digits,
                "log10q": log10q,
                "delta": delta,
                "lin": minus_log10_linear_form,
            }
        )
    return series


def tail_fit(x, y, tail_fraction):
    if len(x) < 4:
        return None
    start = max(0, int(math.floor((1.0 - tail_fraction) * len(x))))
    xx = x[start:]
    yy = y[start:]
    if len(xx) < 2:
        return None
    slope, intercept = np.polyfit(xx, yy, 1)
    return float(slope), float(intercept), xx


def add_fit(ax, x, y, tail_fraction, label_prefix="tail fit"):
    fit = tail_fit(x, y, tail_fraction)
    if fit is None:
        return None
    slope, intercept, xx = fit
    ax.plot(xx, slope * xx + intercept, "--", linewidth=1.2,
            label=f"{label_prefix}: slope={slope:.5g}")
    return slope


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--csv-dir",
        type=Path,
        default=Path("."),
        help="directory containing the three CSV files",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("catalan_approximation_geometry.png"),
        help="output image",
    )
    parser.add_argument(
        "--tail-fraction",
        type=float,
        default=0.5,
        help="fraction of the largest-n data used for asymptotic linear fits",
    )
    parser.add_argument(
        "--theory-delta",
        type=float,
        default=None,
        help="optional horizontal theorem/reference line in the delta panel",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="show the figure interactively after saving",
    )
    args = parser.parse_args()

    series = load_series(args.csv_dir)
    if not series:
        raise SystemExit("No input CSV files were found.")

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    ax_digits, ax_height, ax_delta, ax_linear = axes.flat

    print("Tail-fit summary")
    print("=" * 72)

    for s in series:
        name = s["name"]
        n = s["n"]

        ax_digits.plot(n, s["digits"], marker=".", linewidth=1.4, label=name)
        fit_digits = tail_fit(n, s["digits"], args.tail_fraction)
        if fit_digits:
            sd, id_, xx = fit_digits
            ax_digits.plot(xx, sd * xx + id_, "--", linewidth=1.0)

        ax_height.plot(n, s["log10q"], marker=".", linewidth=1.4, label=name)
        fit_height = tail_fit(n, s["log10q"], args.tail_fraction)
        if fit_height:
            sh, ih, xx = fit_height
            ax_height.plot(xx, sh * xx + ih, "--", linewidth=1.0)
        else:
            sh = None

        ax_delta.plot(n, s["delta"], marker=".", linewidth=1.4, label=name)
        ax_linear.plot(n, s["lin"], marker=".", linewidth=1.4, label=name)

        if fit_digits and fit_height and abs(sh) > 0:
            asymptotic_delta = sd / sh
            print(
                f"{name:20s}  "
                f"digits/n ~ {sd: .8f}   "
                f"log10|q|/n ~ {sh: .8f}   "
                f"ratio ~ {asymptotic_delta: .10f}"
            )

    ax_digits.set_title("Approximation gain")
    ax_digits.set_xlabel("n")
    ax_digits.set_ylabel(r"$-\log_{10}|G-p_n/q_n|$  (digits correct)")
    ax_digits.legend()

    ax_height.set_title("Arithmetic height")
    ax_height.set_xlabel("n")
    ax_height.set_ylabel(r"$\log_{10}|q_n|$")
    ax_height.legend()

    ax_delta.set_title("Empirical worthiness")
    ax_delta.set_xlabel("n")
    ax_delta.set_ylabel(
        r"$\delta_n=-\log|G-p_n/q_n|\,/\,\log|q_n|$"
    )
    ax_delta.axhline(1.0, linestyle=":", linewidth=1.3, label="irrationality barrier")
    if args.theory_delta is not None:
        ax_delta.axhline(
            args.theory_delta,
            linestyle="--",
            linewidth=1.3,
            label=f"reference {args.theory_delta:.8f}",
        )
    ax_delta.legend()

    ax_linear.set_title("Size of the integer linear form")
    ax_linear.set_xlabel("n")
    ax_linear.set_ylabel(r"$-\log_{10}|q_nG-p_n|$")
    ax_linear.axhline(0.0, linestyle=":", linewidth=1.3)
    ax_linear.legend()

    fig.suptitle(
        "Catalan two-row approximation geometry\n"
        "solid = exact CSV data, dashed = tail linear fits",
        fontsize=14,
    )
    fig.tight_layout()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=180, bbox_inches="tight")
    print(f"\nsaved {args.out}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
