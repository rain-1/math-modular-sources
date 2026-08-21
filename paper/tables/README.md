# paper/tables — regeneration notes

Two LaTeX table fragments, built from exact PARI/GP computation:

- `census.tex` (\label{tab:census}) — Table 1, one row per known lattice sequence.
- `budgets.tex` (\label{tab:budgets}) — Table 2, worked two-row alignment/budget cases.

Both are standalone `table`/`tabular` environments (plain `tabular`, no `booktabs` — nothing
else in `paper/sections/*.tex` uses `booktabs`) meant to be pulled in with `\input{tables/census}`
/ `\input{tables/budgets}` from a section file, exactly as `paper/sections/03_archimedean.tex`
already does for `census`. Compiled successfully as a standalone fragment against
`amsmath,amssymb` (see below); braces/environments are balanced.

## How to regenerate

All generator scripts live in `lattice/paper_tables/`. Every `.gp` file ends in `\q` and is run as

```
timeout 110 gp -q lattice/paper_tables/<script>.gp
```

Scripts and what they produce (all numbers appearing in `census.tex`/`budgets.tex` were
transcribed by hand from these runs; nothing was typed in without a corresponding gp
computation, except the items explicitly marked "cited" below):

- `gen_census.gp` — the main driver. Computes, from exact rational companion recurrences:
  - Zagier A–F (order-2, `lattice/zagier_padic.gp`-style rows, N=300),
  - the six third-order AZ-family rows AZ(7,3,81), AZ(9,3,−27), Domb, η=(11,5,125), T=(12,4,16),
    Apéry (`lattice/third_order.gp`-style rows, N=400),
  - Cooper s7/s10/s18 with the **corrected** initial conditions A₁=4,2,6 respectively (from
    `consolidation/SLOPE_CENSUS.md`; the naive A₁=b IC is wrong for s10, s18),
  - the weight-3 cusp-form row L(f,2) (recurrence `(n+1)^2 a_{n+1} = (20n^2+10n+2)a_n -
    16(2n-1)^2 a_{n-1}`, a₀=1,a₁=2, from `consolidation/SPORADIC_SEARCH.md`'s ADDENDUM).
  For each row it prints characteristic roots, the archimedean limit (+ `lindep` against
  {1, ζ(2), L(2,χ₋₃)} where applicable), the denominator exponent `k` with a sharpness check
  (`d_n^k b_n ∈ Z` for n≤200, and `k−1` fails), and measured `p`-adic slopes
  `v_p(b_n/a_n − b_{n−1}/a_{n−1})` at n=300 (or 400) for p=2,3,5,7. Output saved verbatim in
  `census_raw.txt`.
- `gen_zeta7.gp` — copy of `lattice/census/zeta7_level24.gp`, rerun verbatim (confirms the
  33-digit agreement with (1463/13824)ζ(7), integrality of `d_n^7 B_n`, and "no slope at
  p=2,3,5,7").
- `cooper_align.gp` — copy of `lattice/census/cooper_align.gp`, rerun verbatim: cross-determinant
  slopes for (s10,D), (s7,A), (s18,C), (s18,B), (s18,F).
- `helpers/zagier_roots.gp` — sanity check that fixed a bug (see Limitations): correct
  characteristic equation for the **order-2** Zagier family is `x²−ax+c=0` (not `x²−2ax+c=0`,
  which is only right for the order-3 "third-order" family where the recurrence's leading
  coefficient is `2a`, not `a`).
- `helpers/coop_k.gp`, `helpers/coop_denom.gp` — verify Cooper s7/s10/s18 have `A_n ∈ ℤ`
  (κ=0, confirmed to n=300 under the corrected ICs — contradicting an earlier, wrong-IC-based
  claim in `SLOPE_CENSUS.md` that s10's denominators grow) and `k=2` sharp for `B_n` (matches
  the "free integration" claim in `paper/sections/03_archimedean.tex`).
- `helpers/cross_bcf.gp` — cross-determinant slopes for (B,C), (B,F), (C,F) at p=3, with the
  correct archimedean-matching scalars (r=1,1 for B,C; r=5,4 for B,F and C,F).
- `helpers/cross_dombT.gp` — cross-determinant slope for (Domb,T) at p=2: confirms
  `v_2(3·a^T_n b^D_n − 4·a^D_n b^T_n) = 4n−3` exactly for n=50..400.
- `helpers/zud_slope.gp` — reruns Zudilin's exact `Q_m,P_m` recurrence
  (`lattice/catalan_audit/rows.gp`) and independently confirms κ₂=4 and measured σ₂=8.
- `helpers/check_lf2.gp` — diagnostic run isolating the L(f,2) slope discrepancy (see
  Limitations).

To regenerate from scratch:
```
cd lattice/paper_tables
timeout 110 gp -q gen_census.gp | tee census_raw.txt
timeout 110 gp -q gen_zeta7.gp
timeout 110 gp -q cooper_align.gp
timeout 110 gp -q helpers/zagier_roots.gp
timeout 110 gp -q helpers/coop_k.gp
timeout 110 gp -q helpers/coop_denom.gp
timeout 110 gp -q helpers/cross_bcf.gp
timeout 110 gp -q helpers/cross_dombT.gp
timeout 110 gp -q helpers/zud_slope.gp
```
then transcribe the printed numbers into `../census.tex` / `../budgets.tex` (this step is
manual; there is no automated gp→LaTeX formatter in this pass — see Limitations).

## Normalisation

Companion recurrence normalisation `b_0=0, b_1=1` is used throughout **except**:
- the Zudilin Catalan row, whose own natural normalisation is `Q_0=1,Q_1=7/4,P_0=0,P_1=13/8`
  (Zudilin's own integral-representation constants, not (0,1));
- the Nesterenko (4,7) row, built from a Padé/partial-fraction determinant construction, not a
  simple companion recurrence, so `(b_0,b_1)` is not meaningful for it;
- the ζ(7) level-24 row, whose `B_n` is obtained as a $q$-series composition
  `B(z)=A(z)Ψ(q(z))`, not a two-term recurrence with a chosen `(b_0,b_1)`.

## What was NOT independently (re)computed — read before citing this table

- **Nesterenko (4,7) row**: `λ1,λ2,k,σ_p,κ_p` are **cited from
  `consolidation/CATALAN_AUDIT.md`**, not independently decomposed into a single
  `(λ1,λ2,k)` triple here. The audit gives exact rates (`log B_n/n→7.6507`,
  `log|4B_nG−C_n|/n→−7.3109`, `log|V_nG−U_n|/n→E_2=14.3931`) for its own two-integer-sequence
  construction, but disentangling these into the census's `(λ1,λ2,k,σ_p)` columns would require
  re-deriving/rerunning the nested partial-fraction construction in `lattice/catalan_audit/nest_core.gp`
  (heavy: solves a growing linear system per `n`), which was judged out of budget for this task.
  Reported as "cited, not decomposed" in the table footnotes.
- **AZ(9,3,−27) period = L(3,χ₋₃)/3** and **η=(11,5,125) period = ½L(3,χ₅)**: both cited from
  `consolidation/SLOPE_CENSUS.md` §2 (which in turn cites archive source text / the book), not
  independently confirmed by `lindep` against an explicit L-value basis in this pass (would need
  `lfuncreate` for χ₋₃ weight 3 / χ₅ weight 3 and high precision; not attempted here).
- **ζ(7) level-24 row**: `k`, `c`, score, budget are left as "?" — the order-7 recurrence's
  polynomial coefficients are not available in closed form in the source notes (per
  `consolidation/SLOPE_CENSUS.md` §3's own caveat), so `k` cannot be checked by the same
  `d_n^k b_n ∈ ℤ` method used for every other row. `λ1,λ2` are the paper's own analytically
  stated values, independently confirmed only via the empirical decay rate (33-digit match to
  target, decay exponent 0.7059 vs. predicted 1/√2=0.7071).
- **Budget pairs (Zudilin,Zagier E), (Zudilin,Nesterenko)**: the cross-determinant slopes
  (`≈24.13n`, `≥(24−k)n−O(log n)`) and quality values 0.9025 / "Lean 1−ε" are **cited** from
  `consolidation/ARCHIVE_CATALOGUE.md`, `consolidation/ZETA3_TWO_LATTICE.md`, and
  `consolidation/CATALAN_AUDIT.md` respectively — reproducing the full 5:8-sampled construction
  (Zagier E at index 8n against Zudilin at index 3n) or the Lean-certified 2-adic package was out
  of scope for the time budget of this task.
- **(Nesterenko, Zagier E) pair**: no certified construction for this specific pairing was found
  in the consolidation notes (`CATALAN_AUDIT.md` line 91 explicitly states the modular row at 8n
  was aligned against *Zudilin*, not Nesterenko). Reported as an open item ("?"), not fabricated.

## Known discrepancies (flagged, not silently resolved)

- **L(f,2) row, σ₂.** `consolidation/SPORADIC_SEARCH.md` states a predicted σ₂=6 (from `v₂(c)`,
  c=64). My own fresh computation (`helpers/check_lf2.gp`, `gen_census.gp`) measures σ₂≈2
  (raw valuation 590 at n=300, slope per-step ≈1.9–2.1 across n=50..300), with `a_n` confirmed
  integral (κ₂=0) throughout. The discrepancy is real: this row's recurrence coefficient of
  `a_{n−1}` is `−16(2n−1)²`, which is only asymptotically `−cn²` (it has lower-order terms), so
  the exact Casoratian identity underlying the naive `σ_p=v_p(c)` slope law
  (`THEORY_NOTES_03_lattices.md` §1, eq. 1) does not literally apply to this recurrence — the
  formula there was derived for the exact palindromic form `−cn²u_{n−1}`. This is reported as
  measured in the table (with the doc's predicted 6 in parentheses) rather than silently
  overwritten. The **budget** value (+0.7726, using only λ1 and k, not σ₂) is unaffected and
  matches `SPORADIC_SEARCH.md` exactly.
- **Zagier third-order family characteristic roots.** An early version of `gen_census.gp` used
  the same `x²−2ax+c=0` characteristic-equation formula for both the order-2 Zagier family and
  the order-3 AZ family; this is only correct for the order-3 family (whose recurrence's leading
  coefficient is `2a` from the `(2n+1)(an²+an+b)` factor). Fixed via `helpers/zagier_roots.gp`
  before the numbers in `census.tex` were transcribed; all six Zagier roots now match
  `consolidation/THEORY_NOTES_03_lattices.md` §1's table exactly (cross-checked digit-for-digit
  on λ1, λ2, score, budget for all six rows).
- **Cooper `A_n` integrality.** `consolidation/SLOPE_CENSUS.md` states "`denom(A_n^{(s10)})`
  grows (1,4,18,144,600,…)" — but that statement was made while scanning wrong initial
  conditions. Under the corrected IC (A₁=2), `A_n` for s10 is confirmed integral to n=300
  (`helpers/coop_denom.gp`), consistent with `s_10(n)=Σ_k C(n,k)^4` being manifestly integral.
  κ=0 for all three Cooper rows is used in the table on this basis.

## Score/budget convention for complex-root rows

For rows with complex-conjugate characteristic roots (Zagier B; AZ(7,3,81); η=(11,5,125)), there
is no real archimedean limit, so "score" is reported as n/a per the task instructions. "Budget"
is still well-defined and is computed as `log|λ| − k` using the common modulus
`|λ1|=|λ2|=√|c|` of the conjugate pair.
