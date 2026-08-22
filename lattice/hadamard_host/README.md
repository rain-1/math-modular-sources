# `lattice/hadamard_host/` — the Hadamard host of (Zudilin, Nesterenko)

Computation for `consolidation/HADAMARD_HOST.md` (task set in `CDT_UNPACKED.md` §6).
PARI/GP 2.15 + python3.  All paths absolute inside the scripts; intermediate data and
logs in `out/`.  Total runtime ≈ 2 minutes on the dev box (no `snake` needed).

## Order of execution

| script | what it does | output |
|---|---|---|
| `00_csv2gp.py` | converts the published exact integer rows (`catalan-2-row-denominators/*.csv`) to GP vectors | `out/csvrows.gp` |
| `01_rows.gp` | Zudilin $Q_m,P_m$ to $m=1200$; Nesterenko $B_n$ from the closed form to $n=400$; $C_n$ from the published $U_n$; cross-checks every published integer row | `out/rows_raw.gp` |
| `02_bridge.gp` | reproduces $v_2(X_nU_n-V_nY_n)\ge24n-O(\log n)$ and the identity $h_n=2^{e(3n)+14n}D_{6n}^4w_n$ | `out/minor.gp`, `out/02_bridge.log` |
| `lib_fit.gp` | recurrence-fitting library (mod-$p$ sweep, binary search on the degree, CRT reconstruction, exact verification, singular set) — adapted from `lattice/two_prime_holonomy/02_fit_recur.gp` | — |
| `03_rowrecur.gp` | minimal recurrences of the two rows; extends $C_n$ from $98$ to $400$ | `out/rows_full.gp` |
| `04_hadamard.gp` | builds $A_Z\odot A_N$, $A_Z\odot B_N$, $A_N\odot B_Z$, $B_Z\odot B_N$, $W$ to $n=400$; growth rates | `out/had.gp`, `out/04_hadamard.log` |
| `05_fitprod.gp` | order search: the joint module has recurrence order exactly $4$, minimal degree $96$ | `out/05_fitprod.log` |
| `06_operator.gp` | exact order-$4$ operator over $\mathbf Q$; leading polynomial, irreducible, roots $=$ the four products (Hadamard's theorem) | `out/op4.gp`, `out/06_operator.log` |
| `07_denoms.gp` | exact denominator types, $v_p$ profiles for $p\le47$, $2$-adic slopes | `out/cond.gp`, `out/07_denoms.log` |
| `08_odeorder.gp`, `08b_odeorder.gp` | order–degree curve (differential order $\le32$) | `out/08b.log` |
| `09_indep.gp` | $\mathbf Q(x)$-independence rank checks to degree $5$ mod $2^{61}-1$ | `out/09_indep.log` |
| `10_folds.gp` | which singular points each function is regular at (exact $G$) | `out/10_folds.log` |
| `11_inventory.gp` | the full candidate inventory with growth / slope / $k$ / $b$ | `out/11_inventory.log` |
| `13_baselines.gp` | the same measurements on each source row alone (the comparison baseline) | `out/13_baselines.log` |
| `14_pairs.gp` | the `POSITIVITY_PROGRAM.md` adjacent pair $(j_0,j_0+1)$, $n\le24$ | `out/14_pairs.log` |
| `15_conformal.py` | conformal radii $4ab/(a+b)$ of the two-slit domains, univalence check, rearrangement numerators | `out/15_conformal.log` |
| `12_adelic.py` | the adelic arithmetic holonomy bound (reuses `lattice/adelic_holonomy/adelic_bound.py`): entry, margin, deficit for every architecture and inventory | `out/12_adelic.log` |

## Headline

Entry fails on every configuration.  Best: architecture (K), $\{1,W,\theta W\}$,
$\log|\varphi'(0)|=+0.4891$, $\tau^\flat=10.6667$, $\gamma_2=-0.9242$,
**entry $=-11.10$**, margin $-33.79$, deficit $-16.90$ per function
(modular Catalan host: $-0.613$).  At the hard Landau ceiling the deficit is still
$-9.25$.  See `consolidation/HADAMARD_HOST.md`.

## Gotchas

* `lattice/positivity/rows_pos.gp` cannot be `read()` from another GP script: its own
  `default(parisizemax,...)` aborts the read.  `14_pairs.gp` therefore carries a copy of
  `mom(m,j)`.
* In GP, `t << k` with negative `k` on a `t_FRAC` truncates; use `t * 2^k`.
* Multi-line `for(...)`/vector literals outside `{ }` are parsed line by line and silently
  break.
