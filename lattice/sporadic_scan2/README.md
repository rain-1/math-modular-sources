# `lattice/sporadic_scan2` — scan 2 for Apéry-like rows

Report: `consolidation/SPORADIC_SCAN2.md`. Everything here is regenerable; see §9 of
the report for the exact command sequence.

## Pipeline

| file | what it does |
|---|---|
| `qser.py` | integer $q$-series, eta quotients, Ligozat cusp orders |
| `gen_t.py` | enumerates eta-quotient parameters $t=q+O(q^2)$ by **divisor** (cusp-order vector), $\deg t\le4$, all $N\le120$ → `t_list.json` (3979 parameters, 51 levels) |
| `gen_spaces.gp` | integral lattice of weight-$w$ forms with rational $q$-expansion in one nebentypus Galois orbit (Galois traces + `matkerint` saturation); `saturate_fast` saturates at 210 coefficients and transports to `PREC` |
| `run_spaces.gp`, `w_c*.gp`, `e_s*.gp` | drivers → `spaces/`, `spaces600/` |
| `prep.gp`, `run_prep.gp`, `p_c*.gp`, `e_p*.gp` | splits each lattice into $F_0$ (constant term 1) + LLL-reduced $a_0=0$ sublattice → `prep/`, `prep600/` |
| `scan2.py` | the **bilinear solve**: nullspace of the linearised (recurrence $\otimes$ $F$) system mod $p=99999989$, rank-one extraction by a matrix pencil + Cantor–Zassenhaus, rational reconstruction of $m$. No search box. → `hits*.jsonl` |
| `verify.py` | exact recomputation: $a_n$, minimal $(r,D)$ recurrence over $\mathbf Q$ verified to $n=200$, $\lambda_i$, $c$, both companions (Eichler $B=F D^{-(w+1)}(F\,Dt)$ and the census $b_0{=}0,b_1{=}1$ solution), sharp $k$, slopes $\sigma_p$, Apéry limit to hundreds of digits → `rows*.jsonl` |
| `dedup.py` | de-duplicates by $a_n$ (up to $(-1)^n$), flags the census rows, ranks → `table.json`, `limits.txt` |
| `ident.gp`, `identlib.gp` | `lindep` against constants and against $L(f,s)$ for every newform of weight $w+2$ at every level dividing $M$ or $2M,3M,4M$ → `ident.out` |
| `build_tables.py`, `fmt_period.py` | → `table.md` (all 85 rows) and `top_table.md` (top by budget, $k\ge1$); identifications are re-keyed by the canonical $a_n$ via `table_prev.json`, so they survive a re-run of `dedup.py` |
| **`sym1.py`** | **Part II**: power-series square root of each order-three sporadic row, minimal integralising scale $\lambda$, refitted recurrence, $k$, score, budget, limit → `sym1.json` |
| **`family_scan.py`** | **Part II**: independent Zagier-style integrality sweep of the normalisation $(n+1)^2a_{n+1}=(An^2+\tfrac A2n+B)a_n-C(2n-1)^2a_{n-1}$ → `family_hits.json` |
| `sym1_limits.py`, `ident_sym1.gp` | high-precision limits and identification for the $\operatorname{Sym}^1$ rows |
| `drive_verify.sh`, `drive_verify_local.sh` | shard `verify.py` over cores |

## Raw results kept here

`t_list.json` (parameters), `prep/` (the integral form lattices, cosets and LLL bases —
the one non-regenerable-in-seconds intermediate), `hitsA_*` ($\deg t\le2$, precision 205),
`hitsB_*` ($\deg t\in\{3,4\}$), `hitsC_*` ($\deg t\le2$, precision 400, larger sublattice),
`rowsA/B/C_*` (verified rows), `table.json`/`table.md` (85 distinct rows, ranked),
`top_table.md`, `ident.out`, `sym1.json`, `sym1_cooper.json`, `family_hits.json`,
`sym1_limits.txt`, `ident_sym1.out`.

`spaces/` and `spaces600/` are deleted after `prep/` is built (regenerate with the gp
drivers if needed).
