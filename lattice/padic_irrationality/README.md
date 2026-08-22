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

## Follow-up: the shifted-normalisation scan (section 9 of the write-up)

| script | what it does | log |
|---|---|---|
| `20_scan.c` → `20_scan` | C/OpenMP sieve of $(n+1)^2u_{n+1}=(a_qn^2+b_qn+c_q)u_n-c(n-r)^2u_{n-1}$ over $\sim1.5\cdot10^{11}$ triples; emits rows whose $A$-solution survives $\ge12$ exact divisibility steps | `20_scan.out`, `20_scan.err` |
| `20_verify.gp` | exact PARI re-derivation to $n=60$ then $n=300$; $k$, $\kappa_p$, measured $\sigma_p$, $S_p$, $\theta_p$ | `20_verify.log` |
| `20_analyse.gp` | diagnostics for the $S_p>0$ rows: true archimedean rate, $\max v_p(A_n)$, nonvanishing, rationality | `20_analyse.log` |
| `20_ident.gp` | affine $p$-adic `lindep` identification of $\xi_2$ against Kubota–Leopoldt values | `20_ident.log` |
| `20_hunt.gp`, `20_hunt2.gp` | deep identification attempts for the cells `20_ident.gp` could not match | `20_hunt.log`, `20_hunt2.log` |
| `20_family.gp` | the family $(n+1)^2u_{n+1}=4(2n+1)u_n-cn^2u_{n-1}$: which $c$ are integral, and the deep check of $c=-4^m$ | `20_family.log` |

Build the C sieve with
`gcc -O2 -fopenmp -o /home/ubuntu/code/math-modular-sources/lattice/padic_irrationality/20_scan /home/ubuntu/code/math-modular-sources/lattice/padic_irrationality/20_scan.c -lm`.
