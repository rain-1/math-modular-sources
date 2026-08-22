# lattice/padic_irrationality

Scripts for `consolidation/PADIC_IRRATIONALITY_CENSUS.md`.

| script | what it does | log |
|---|---|---|
| `calegari_rows.gp` | rebuilds Calegari (math/0408214) §3 rows on $X_0(p)$, $p=2,3,5,7,13$ and §4 row on $X_1(4)$ from their modular definitions; reproduces his printed $a_n,b_n$; measures $\sigma_p$, $\log\lambda_1$, $k$ | `calegari_rows.log` |
| `census_scores.gp` | six Zagier + six AZ sporadics to $n=500/420$: $k$, $\kappa_p$, measured $\sigma_p$ at $n=400$, archimedean rate, $S_p$, $\theta_p$ | `census_scores.log` |
| `cooper_scores.gp` | Cooper $s_7,s_{10},s_{18}$ to $n=420$ | `cooper_scores.log` |
| `zudilin_check.gp` | Zudilin's Catalan row to $m=420$; the $\kappa_p>0$ worked example | `zudilin_check.log` |
| `scores.gp` | assembles the whole table sorted by $S_p$; cusp-move placement scan | `scores.log` |

Run with `gp -q <absolute path>`.  Note (GP gotchas hit while writing these):
outside a `{ }` block gp ends a statement at end of line, and a function body
written `f(x) = a; b;` swallows `b` into the closure.
