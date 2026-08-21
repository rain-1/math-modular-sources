# Consolidation plan (started 2026-08-21)

## 1. What exists (inventory)

| Strand | Canonical artefacts | Status |
|---|---|---|
| **S1. Hecke–Eisenstein "three-layer" book** | `book/` V8 (+V11 patch) | Concrete claims verified (see `verify/`); framework over-formalised, no code |
| **S2. Companion inversion / Eisenstein sources for the 15 sporadics** | ledger C01–C05, E01–E05; `packages/phase 0 .../eisenstein/EISENSTEIN_SOURCE_THEOREM.tex` + Lean project (mis-filed under `catalan/`) | Proved core; closest to publishable |
| **S3. Domb cusp-form apparatus** | ledger E06–E14; `phase 0/domb/*` | Gaps E07/E09/E10 claimed closed in phase-0 PDFs — needs audit |
| **S4. Harmonic jets / Catalan companion / second-row Lucas** | ledger H01–H21; `S7_LUCAS_COMPLETION.tex` | Proved core + Lean for H13/H17 |
| **S5. Catalan two-row lattice ("multiple-lattice")** | ledger D01–D04; `phase 0/catalan/CATALAN_TWO_ROW_FULL_LCM_PAPER.tex` (0.65887) + Lean (mis-filed under `eisenstein/`); `catalan-2-row-denominators/5-8 theorem/` (0.90253, 2-adic) | 0.65887 proof-shaped; 0.90253 has an "idea of proof" gap in the key limit lemma |
| **S6. Brown–Zudilin ζ(5) Frobenius / cellular jets** | ledger F01–F07, B01; `ZETA5_RESEARCH_ALL_IMPORTANT_DISCOVERIES`, `BZ_KAPPA_MINIMALIZATION_AUDIT.tex` | Mixed; κ closed forms numerical only |
| **S7. p-adic Apéry limits / Λ-algebra** | ledger P01–P15 | Proved core + conditional (IND_p) |
| **S8. ζ(5)/ζ(7) modular searches** (X0+(169), level 12/24/30/60) | `zeta5_*`, `level30_*`, chat archive | Many rigorous no-go certificates; no positive result |

Detailed per-file inventory: `ARCHIVE_CATALOGUE.md` (chat archive, papers, catalan dirs) and the ledger dump (`LEDGER_DUMP.md`).

## 2. Diagnosis

1. **The mathematics is real.** Every concrete claim I re-derived (ζ(7) source 209/1728, ζ(5) level-16 projector, β(4) inner/outer vectors, β(4) Lambert identity to 40 digits, v₂(Δ_{5n,8n}) lemma at n=100) checks out.
2. **Provenance is the weak point.** The chat archive stores message text only; tool code/outputs are gone. "Certificate PASSED" statements are unverifiable as they stand.
3. **The book's framework outruns its theorems.** V8/V11 wrap ~12 genuinely new concrete results in Layer-A/B/C, adelic, free-depth, syntomic language much of which is definitional or open. Merging V11 into V8 would deepen this, not fix it.
4. **The ledger is the best map but is stale**: E01/E02, E07/E09/E10, H18, D03 all have later "closure" documents.
5. Duplicates / mis-filed items (Lean tarballs swapped; `unknown packages/` redundant).

## 3. Plan

### Phase A — Reproducible verification layer (`verify/`)  ← started
One script per concrete claim, runnable with `python3`/`gp`, printing PASS/FAIL. Order of attack = ledger Action Queue + the book's concrete tables:
- A1 Dirichlet-polynomial/period-annihilation facts (done: `verify/period_annihilation.gp`)
- A2 β(4) Lambert identity (done: `verify/beta4_lambert.py`)
- A3 Eisenstein-source identities for the 12 sporadics to a Sturm bound (closes E02 honestly)
- A4 Domb cusp source identity + d_n³B*_n (E07/E09) and L(f6,3)/2 numerics (E10)
- A5 Catalan 5:8 valuation lemma for all n≤200 and the 0.65887 / 0.90253 quality computations
- A6 β(4), L(4,χ₋₃), ζ(7) level-24/36 Apéry pairs: recompute sequences, limits, denominators
- A7 CDT free-depth numbers (14,3,1), τ_F=191/49 re-derived from the CDT paper

### Phase B — Update the ledger
Re-grade each row from the verification outputs; mark closures from phase-0 papers; retire D03→0.65887.

### Phase C — New book ("Version 1 of the real thing")
Not a V8+V11 merge. Structure by *theorems*, with the three-layer picture as an organising narrative, not as load-bearing formalism:
1. Period annihilation as finite Hecke linear algebra (barycentric theorem, Gaussian-binomial projectors, Fricke count) — fully proved, elementary.
2. The Eisenstein Source Theorem for the sporadic families + companion inversion (S2).
3. Predicted systems: ζ(3) prime completions, ζ(5) level 12, ζ(7) parent, β(4) level 24, L(4,χ₋₃) — each with a script.
4. Benchmarks: CDT L(2,χ₋₃) and the Tate controls.
5. Catalan: the two-row lattice method (S5) — state 0.65887 as proved, 0.90253 as conditional with the gap named.
6. Obstructions / no-go census (ζ(5) X0+(169), Γ₁(7), level-60 ζ(7)).
7. Programme part (clearly labelled): adelic Layer C, free-depth, p-adic functoriality, BZ bridge.
Paper bundles A–I from the ledger remain the publication units; the book is the connective tissue.

### Phase D — Expository application (plan.txt)
Run Zudilin's three Catalan papers through the classification (which Eisenstein source/orientation is his recurrence; what the permutation group is geometrically).

## 4. Housekeeping
- Fix swapped Lean tarball directories; delete `unknown packages/`; gitignore the 56M archive or commit it deliberately.
- Recover "the lost exact asymptotic formula for the linear-form size" (tasks.txt) by grepping the archive.
