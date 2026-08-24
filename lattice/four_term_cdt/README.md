# `lattice/four_term_cdt` — the CDT entry and margin for the four-term Catalan hosts

Report: `consolidation/FOUR_TERM_CDT.md`.

The six mixed-exponent five-singular-point four-term rows of `FOUR_TERM_DEEP.md` §6.4
(class $(-\tfrac12,0;1,0,0)$) carry Catalan- and $L(2,\chi_{-3})$-class periods on hosts
with **irrational conjugate singular points** — evasion (E1) of `CATALAN_OBSTRUCTION.md`
§3 — so the width law, which is proved only for covers of $\mathbf P^1\setminus\{0,1,\infty\}$,
does not apply and the conformal size has to be computed honestly.

| file | what it does |
|---|---|
| `01_geometry.py` | exact geometry: characteristic cubic, the five singular points, local exponents, Jordan blocks, which point is the fold, the post-hypothesis set, the monotonicity ceiling → `out/geometry.json` |
| `02_periods.sh` | re-verifies the full period matrices with `four_term_deep/11_foldmix.py` (Frobenius + Taylor continuation, two precisions) → `out/02_periods.log` |
| `03_foldreg.py` | fold-regularity by **radius enlargement**: $\log\lvert b_n-\xi a_n\rvert/n$ against $-\log\lvert t_2^{\rm post}\rvert$, with a perturbed-$\xi$ control → `out/foldreg.json` |
| `04_denom.py` | the denominator array $\mathbf b$ (sharp $k$, best $(k,e)$), $p$-adic slopes of $a_n$, and the **pure algebraic functions** $\sqrt{1-\sigma t+\pi t^2}$, $\sqrt{1-rt}$ → `out/denom.json` |
| `05_uniformise.py` | **numerical Fuchsian uniformisation** of a 4-point orbifold on $\mathbf P^1$ (cusp at 0, cusp at ∞, two more points of type cusp / order-2): accessory parameter by the "one common horizontal line" condition on the monodromy, cusp size $=e^{-2\pi h}$.  `python3 05_uniformise.py test` runs the five exact validations ($\lambda$; the $z\mapsto z^2$ cover; $X_0(6)$; $\Gamma_0(2)$; and CDT's $256s$) |
| `06_sizes.py` | the post-hypothesis conformal sizes of the six rows → `out/sizes.json`, `out/06_sizes.log` |
| `07_padic.py` | $R_p$ for Theorem A of `ADELIC_HOLONOMY.md` §2.6: slopes of $a_n$, $b_n$, a conditional combination, and the polylogarithm module → `out/padic.json` |
| `08_entry.py` | entry and margin in `CDT_FINDER.md` conventions, with the benchmarks → `out/entry.json` |
| `09_second.py`, `09b_second_slow.py` | is there a **second** conditional direction?  fold constants of the solutions of $Ly=t^k$, $k=1,2,3$, against $1,\zeta(2),G$ / $L(2,\chi_{-3})$ by PSLQ → `out/second.json`, `out/second_slow.json` |
| `10_symm.py` | all involutions of $\mathbf P^1$ fixing 0 and permuting the singular set, before and after the hypothesis; why symmetrisation is unavailable |

Headline: entry at the hard uniformisation ceiling is $-1.914$ (row 1, $\xi=\tfrac14G$),
$-0.703$ (rows 2, 3), $-4.081$ (row 4), $+0.008$ (rows 5, 6), against level 8's
$-0.077$.  **No host carrying $G$ beats level 8**; the two that do carry
$L(2,\chi_{-3})$, where the target is already CDT's theorem.
