# `lattice/thmB_exact/` — the exact rational prefactor of Theorem B

Scripts for `consolidation/THEOREM_B_EXACT.md` and
`paper/sections/02_sources_thmB_exact.tex`.  PARI/GP; run with a large stack,
e.g. `gp -q -s 12000000000 01_limits.gp`.

| script | what it does |
|---|---|
| `common.gp` | the twelve rows; Frobenius -> nome -> exact integral `q`-expansions of `t` and `F` from the recurrence alone; closed-form Eisenstein coefficients `c(m)`; Mellin polynomial `P(s)`; `L(Phi,s)`; endpoint criterion |
| `01_limits.gp` | direct exact `B_n/A_n` at `n=2600` versus `L(Phi,w+1)` for the nine real-fold rows (no modular input on the left) |
| `02_cuspvalue.gp` | `Theta(iy) -> L(Phi,w+1)` as `y -> 0`, all twelve rows; the endpoint criterion table |
| `03_fricke_fold.gp` | the six third-order rows: Fricke eigenvalues of `Phi`, `F`, invariance of `t`, the fold at `i/sqrt(N)`, and the fold-connection value |
| `04_complex_folds.gp` | the complex folds of `delta`, `eta` (interior, CM) and of `B` (cusp `1/6`), to 50-90 digits, with the closed forms |

`build(r,a,b,c,NT)` returns `[t(q), F(q)]` as exact power series with integer
coefficients, obtained from the recurrence by
`dlog q/dlog t = 1/K`, `K = P_2 y_0^2` (r=2) resp. `K = y_0 sqrt(P_3)` (r=3),
series reversion, and `F = sqrt(K_q/P(t))` resp. `K_q/sqrt(P(t))` — no eta
quotient or modular input is used.
