# Genuine contributions ledger — full programme, non-Catalan focus (2026-08-29, rev. 2)

Scope: everything of genuine value **outside** the closed Catalan/PTPAH strand. Sources:
`SURVEY_NONCATALAN_NOTES.md` (strand notes, with the 27-item conflict catalogue in its
§7) and `SURVEY_OUTPUTS.md` (papers, book, artefacts). The Catalan-specific salvage
(rigidity theorem, obstruction package, 2-adic cluster) is inventoried separately in
`catalan/PTPAH_K_final_resolution_2026-08-27.md` §4 and rev. 1 of this file (git history).

Status codes: **[T]** proved, **[V]** verified computationally, **[C]** conjecture with
evidence, **[M]** method/artefact, **[D]** drafted paper.

---

## Tier 1 — shareable essentially as-is

### 1.1 Classification paper **[T][D — ready]**
*Integral second-order Apéry-like recurrences and rational elliptic surfaces* (21 pp,
`paper/classification/main.tex`). Answers Zagier's integrality question: twist parity
forces $12 \mid \deg\mathcal J$, so no Kodaira-dimension hypothesis is needed and the
surface is one of Beauville's six; cusp placements reproduce Zagier's seven rows; the
exhaustive scan ($1.27\times10^6$ hits) adds two new elliptic-surface rows. Self-contained,
inline bibliography, compiles. **Submit.**

### 1.2 AESZ 207 correction **[V — 1139 digits]**
The published Almkvist–van Straten–Zudilin Apéry limit for AESZ 207 is **wrong from the
4th significant figure** (their value was un-converged — which explains why no PSLQ
identification was ever found). Corrected value to 1139 digits with a validated pipeline
(`AESZ207.md`). Cheapest high-value item in the repo: a short correction note serves the
CY-operator community immediately.

### 1.3 K3 row note **[T][D — needs novelty sweep only]**
`paper/k3_row/` (3–4 pp): the first Apéry-like sequence whose Apéry constant is a
critical cusp-form $L$-value — $L(g,2)$ for the CM newform 32.3.d.a. Filed as
Catalan-adjacent but stands alone; highest novelty-per-page in the tree.

### 1.4 Improved p-adic irrationality measures **[T]**
$\mu(\zeta_3(3))\colon 22.3 \to 10.0$; $\mu(\zeta_2(5))\colon 19.744 \to 19.39$;
$\mu(\zeta_2(3))\colon 7.18 \to 6.3$. Independent of the ζ₅(3) priority situation —
these improve *known* results regardless of who proved what first. Short paper or a
section of 2.1.

## Tier 2 — near-submittable / paper-shaped

### 2.1 Main paper **[D — blockers listed]**
*Modular Apéry systems: sources, slopes, and the two-row design rule* (51 pp,
`paper/main.tex` = `paper_draft8.pdf`). Theorem A (all twelve modular sporadic families
have holomorphic-Eisenstein source, zero cuspidal projection); Theorem B\* (the Apéry
limit **is** $L(\Phi,w+1)$, prefactor from the Mellin polynomial, no case-by-case
constants; endpoint condition separates the nine real-fold from three complex-fold
rows); Prop C (slope law $\sigma_p = v_p(c) + 2\kappa_p$ with product formula); Theorem
E (closed-form two-row quality + design rule); Theorem F (below); Conjecture D proved
for Zagier B, C, F at $p=3$ — note the strand survey found **row B's gap is closed in a
later update** that the paper may not yet reflect. Blockers: three `\todo`s, five cited
companion papers absent from the tree, and the row-B descent status to reconcile.

### 2.2 Theorem F + Dwork-crystal reformulation **[T + the deepest idea in the corpus]**
Euler-factor criterion: a $p$-adic Apéry limit exists **iff** $\mathcal E_p(s) =
1-\psi(p)p^{-s}$ divides the Mellin polynomial, with Kubota–Leopoldt value
$\xi_p = -Q(w+1)\kappa_p$; verified to $\ge p^{2991}$; transfers verbatim to rank 4.
The reformulation (`CRYSTAL_THEOREM_F.md`): the criterion is Frobenius-slope separation
on an extension of a MUM Dwork crystal — so it applies to **non-modular** rows (e.g.
Zudilin's Catalan row has a Dwork crystal with tower eigenvalue $\chi_{-4}(p)p^{-2}$).
This deserves to be the conceptual centerpiece of 2.1 or its own paper. (Conflict to
resolve first: `CRYSTAL_THEOREM_F.md` vs `SOURCES_S18_ZUDILIN.md` §5 on whether
Zudilin's row has a geometric model.)

## Tier 3 — methods, technique, artefacts

- **ζ₅(3) technique residue [T, overlap unquantified].** The theorem itself is out but
  was independently proved by two others (one covering further $p$-adic zeta values).
  Candidate residual novelty: the archimedean template optimization over holomorphic
  self-maps of the $q$-disc, and the §2.6 multi-place bound with weight $1-1/m$
  (`ADELIC_HOLONOMY.md`). **Action: compare against the two independent proofs before
  claiming either as new.**
- **Falsification methodology [M].** `CATALAN_AUDIT.md`'s control-experiment discipline
  (the `bestappr` control; "at $F_n<0$ Minkowski *is* Dirichlet") — would save others
  months; publishable as a short methods note despite the Catalan provenance.
- **`verify/` certified-numerics suite [M].** Re-run clean during the survey (11/11,
  error $4.6\times10^{-41}$). Reusable.
- **Censuses [V].** Slope census, $p$-adic irrationality census, noncongruence scans,
  Herfurtner window, sporadic scans — database-grade material backing 2.1.

## Tier 4 — reference material and stranded items

- **The book** (`book/`: *Hecke–Eisenstein Extensions and Modular Apéry Systems*, 75 pp
  v8 + unmerged 15 pp v11 supplement; compiles). Honest about its own status; framework
  outruns theorems; some certificates absent, one no-go certificate float64-only.
  Working-draft shareable with caveats — treat as the prehistory/companion of 1.1 + 2.1.
- **Stranded ζ(5) linear-form asymptotic.** `CLAUDE_FINDINGS.md` holds an exact
  asymptotic for the linear-form size, in no paper, blocked on one undetermined CM
  constant $K$. This appears to be the "awesome formula that's gone missing" from
  `chats/tasks.txt`. Recover it.

## Hygiene — fix before sharing anything (from the conflict catalogue, §7 of the survey)

1. **$\zeta_p$ notation collision**: `MUM_SURVEY.md` uses $\zeta_p(k) =
   L_p(k,\omega^{1-k})$; `EULER_CRITERION.md` uses $L_p(s,\mathbf 1)$ — they disagree at
   exactly $p=5$, $k=3$. The released ζ₅(3) paper uses the standard convention (safe),
   but `MUM_SURVEY.md`'s "$-40\zeta_p(3)$" is a different number. Unify before any
   further release.
2. `RIGIDITY_PROOF.md` "Theorem R" is withdrawn by its own header — cite
   `CONJ_D_PROOF.md` + `ACF_ONE_SURFACE.md` instead.
3. `CDT_FINDER.md` §4/§6's scoring of $X_1(5)\,\mathrm{Sym}^2$ is vacuous ($\lambda_1$
   double) and it misidentifies Beukers' row.
4. The $\delta>1$ "contradiction" in `ZETA3_TWO_LATTICE.md` §16 is resolved by
   `CATALAN_AUDIT.md` — annotate it so no one rediscovers the scare.
5. 23 further minor conflicts catalogued in `SURVEY_NONCATALAN_NOTES.md` §7.

---

## Recommended order of operations

1. **Submit 1.1** (classification paper) and post **1.2** (AESZ 207 correction note) —
   both are days of work, not months.
2. **Novelty-sweep and release 1.3** (K3 row note).
3. **Unblock 2.1**: reconcile row B, clear the three `\todo`s, decide the companion-paper
   question; fold in 1.4 and make 2.2 the centerpiece.
4. **Overlap-check the ζ₅(3) technique pieces** against the two independent proofs;
   whatever survives becomes a methods note.
5. Run the hygiene list before any of the above goes out.
6. Formalization candidates (in Lean-tractability order): Theorem F's Euler-factor
   criterion (algebraic, finite checks), the slope law Prop C, then the classification
   paper's parity/degree argument; the analytic material (B\*, holonomy bounds) is not
   currently formalizable at reasonable cost.
