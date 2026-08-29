# Survey of polished outputs (non-Catalan, non-$\zeta_5(3)$)

*Compiled 2026-08-29. Scope: everything in the repo that is paper-shaped or artefact-shaped,
excluding the released $\zeta_5(3)$ flagship (state confirmed only) and excluding the closed
Catalan strands (recorded as one-liners for completeness).*

**Headline finding.** The repository's most finished non-Catalan output is not the book. It is
the pair `paper/main.tex` (51 pp, "Modular Apéry systems") and `paper/classification/main.tex`
(21 pp, "Integral second-order Apéry-like recurrences and rational elliptic surfaces"), both of
which compile cleanly with zero undefined references, carry per-claim evidence tags, and are
supported by a verification suite I re-ran successfully. The book (`book/`) is the older,
superseded framework layer that these two papers grew out of.

---

## 0. Orientation: root files

| File | Size | Content |
|---|---|---|
| `README.md` | 0 bytes | **Empty.** The repo has no top-level orientation document. |
| `CLAUDE_FINDINGS.md` | 55 lines | A single dated findings entry (2026-08-21): the recovered full leading asymptotics of the $\zeta(5)$ linear form from the level-12 CM-isogeny construction. Result: with $W=V-C_5U$, $C_5=\tfrac{11}{144}\zeta(5)$, Fricke invariance makes $W$ regular at $x_+=7-4\sqrt3$, giving $d_n/c_n-C_5\sim\tfrac{13}{96K}(-1)^n(7-4\sqrt3)^n n^{1/2}(\log n)^4$ with $K\approx0.6853718485$. **Open piece:** no closed form for $K$ (identified as a CM period at discriminant $-48$). Not incorporated into any paper. |
| `plan.txt` | 23 lines | Early (2026-08-19) consolidation sketch. Largely executed: the "geometric modular source theory" became the main paper, "multiple-lattice" became §5 of it. The expository idea (run Zudilin's three Catalan papers through the classification) is **not done**. |

The real orientation documents are `consolidation/CONSOLIDATION_PLAN.md` (2026-08-21 inventory +
diagnosis of eight strands S1–S8) and `consolidation/CONTRIBUTIONS_LEDGER.md` (2026-08-29,
five-tier grading of what survives). The latter is Catalan/$\zeta_5(3)$-focused and **does not
mention the classification paper or the main modular-Apéry paper at all** — a gap this survey
fills.

---

## 1. `book/` — the Hecke–Eisenstein three-layer monograph

**Title.** *Hecke–Eisenstein Extensions and Modular Apéry Systems: A Three-Layer Classification
of Dirichlet $L$-Value Realizations*, Version 8, dated 18 August 2026. No author line.

**Artefacts.** `HECKE_EISENSTEIN_APERY_THEORY_V8.pdf` (**75 pp**) + LaTeX in `v8/` (21 files,
3970 lines); `HECKE_EISENSTEIN_APERY_V11_INTEGRATION_SUPPLEMENT.pdf` (**15 pp**) + `v11/`
(833 lines); two source zips. Only 4 files of `book/` are git-tracked — the PDFs and zips are,
the loose `.tex` trees appear untracked.

**Both compile cleanly** — verified by re-running `pdflatex` twice on copies: v8 exits 0 at 75
pages, v11 exits 0 at 15 pages, both byte-plausible against the shipped PDFs. Only cosmetic
hyperref Unicode warnings.

### Scope and table of contents (v8, condensed)

1. Introduction — 2. The three-layer object — 3. Parity dichotomy / opposite-parity Eisenstein
extensions — 4. Hecke projectors, interpolation, Fricke folding — 5. Franke critical kernels and
Hecke–Eisenstein purification — 6. Worked example: the level-100 Franke–$\eta$ boundary test —
7. The Tate side: central projectors and even-zeta controls — 8. Affine extension geometry and
boundary purification — 9. Companion descent, exactness, vector-valued realizations — 10. Adelic
Layer C: arithmetic holonomy and Frobenius filtration — 11. Applications portfolio A–P —
12–16. Standalone application chapters ($L(2,\chi_{-3})$ CDT benchmark; $\beta(4)$ at level 24;
$L(3,\chi_5)$ elliptic folds; Shvets Belyi pullbacks) — 17. The parity checkerboard —
18. Research programme and publication architecture — 19. Appendix — 20. References.

### What it claims

The organising claim is that a modular Apéry phenomenon is a three-layer object:

- **Layer A** — a pure local system $\widetilde{\mathbb V}_{r,\chi}=\Sym^{r-1}\widetilde{\mathbb H}\otimes\mathbb L_\chi$
  on an orientation-resolved modular cover, together with its finite direct images under a
  rational quotient coordinate.
- **Layer B** — a mixed extension $0\to\mathbb V_{r,\chi}\to\mathbb E_\nu\to\mathbf 1\to0$ whose
  distinguished regulator is the opposite-parity Dirichlet value $L(r,\chi)$.
- **Layer C** — an *adelic* Apéry realization: marked de Rham section, companion descent,
  coefficient coordinate, arithmetic lattice, continuation leaf, and finite-place
  Frobenius/free-depth data.

Concrete theorem-level content: finite interpolation theorems for critical-period annihilation;
explicit prime-power projectors; Fricke/Weyl rank formulas; a cohomological criterion for scalar
vs. vector companion descent; a pullback-defect formula; a quantitative Frobenius-filtered
refinement of CDT's holonomy penalty (if the $j$th LCM layer is free on a subspace of dimension
$u_j$ inside an $m$-dimensional holonomy space, the penalty drops by exactly $b_ju_j^2/m^2$);
the level-100 Franke–$\eta$ boundary computation (orientation-resolved critical kernel is the
unique line $(1,-10,16)$, full boundary kernel zero); the $L(3,\chi_5)$ conjugate-fold geometry
(branch constants $\Lambda/2\pm2i\pi^3/(125\sqrt5)$); a $\beta(4)$ level-24 stress test with a
degree-four modular direct image of Picard–Fuchs rank 16; and a Shvets Belyi/symmetric-cube
external test giving a rational depth-four realization with target $7L(4,\chi_{-3})/32$.

**v11** is a *supplement*, not a merge: it claims (i) an exact boundary-completed replacement of
Franke's critical kernel, (ii) a reconstruction of the CDT $L(2,\chi_{-3})$ proof inside Layer C
(free-depth flag $14\to3\to1$, flat penalty $\tau_F=191/49$, new horizontal-integration invariant
$e$ with $\tau_{\rm hor}=27/80$, total $\tau_{\rm ad}=16603/3920$), and (iii) a monotone
entry-obstruction principle giving a computer-assisted level-20 no-go for $L(3,\chi_5)$.

### State

35 numbered results (14 theorems, 20 propositions, 1 corollary), 23 `\begin{proof}` blocks; the
remainder are argued in unlabelled prose. Zero `TODO`s. The 4-page appendix is derivation notes
and a dictionary, **not** a proof repository. The book's own §1 status section is unusually
honest, triaging every claim as *Literature* / *Programme* (own certificates, "require a normal
publication audit before novelty claims") / *Open*, and explicitly stating: **"The paper does not
claim irrationality of Catalan's constant, $\zeta(5)$, $\zeta(7)$, $L(3,\chi_5)$, or any other
individual special value whose irrationality is currently open."** No dependence on $\zeta_5(3)$
(that value never appears). Catalan appears only as Application F, a *structural* result
(intrinsic rank-two orientation monodromy; irreducible dihedral orientation of order 8, hence no
scalar normal function), with the irrationality question explicitly deferred.

### Judgment — **NEEDS WORK / partly superseded**

The repo's own diagnosis (`CONSOLIDATION_PLAN.md` §2.3) is the right one and I endorse it: *"The
book's framework outruns its theorems. V8/V11 wrap ~12 genuinely new concrete results in Layer-A/B/C,
adelic, free-depth, syntomic language much of which is definitional or open. Merging V11 into V8
would deepen this, not fix it."* Specific blockers:

1. **v11 is unmerged and contradicts v8's open list** (it converts v8's open boundary-Franke
   problem into a theorem). Shipping both together is confusing.
2. **Several headline theorems are certificate-level or conditional**: the $\beta(4)$
   degree-four elimination and rank-16 cyclicity; $L(3,\chi_5)$ irreducibility; the
   Frobenius–holonomy corollary (a formal consequence only *once* an all-$n$ free-depth flag is
   proved).
3. **Programme claims cite companion papers/certificates not in the tree** — an outside reader
   cannot verify the $q^{70}$ elimination or the denominator-lattice computations from the text.
4. **v11's $L(3,\chi_5)$ no-go certificate is float64**, not interval arithmetic; the text itself
   flags that an Arb/MPFI rerun is wanted.

Shareable *as a working draft to a close colleague* with the "v8 draft, v11 patch pending merge,
several items certificate-level" caveat — which its own status page already says. Not shareable
as a submission. The planned "Phase C new book" of `CONSOLIDATION_PLAN.md` was, in effect,
executed instead as the two papers in §2 below; the book is now best read as their prehistory.

---

## 2. `paper/` — the working paper tree

`papers/` (28 M) is **not** ours: it is a literature archive (Beukers 1987, Zagier's Apéry-like
survey, Calegari math/0408214, Dimitrov's LNT $\zeta(2),\zeta(5)$ notes, Bouw–Möller accessory
parameters, Movasati–Reiter Heun, and the CDT arXiv sources 2408.15403 / 2510.04156 with their
`.tex`). No own drafts. **Reference material, not an output.**

Every PDF in `consolidation/` is a byte-identical copy of the corresponding `paper/*/main.pdf`
(md5-verified for `paper_draft8.pdf` ↔ `paper/main.pdf` and `classification_paper.pdf` ↔
`paper/classification/main.pdf`). The `paper/` tree is the live source; `consolidation/` is the
snapshot shelf.

### 2.1 `paper/main.tex` — **Modular Apéry systems: sources, slopes, and the two-row design rule**

*= `consolidation/paper_draft8.pdf`. 51 pp, amsart, no author line (provenance paragraph names
River as first author with GPT-5/Claude assistance). Compiles clean: `Output written on main.pdf
(51 pages)`, **zero** undefined references, **zero** undefined citations.*

**Structure.** §1 Introduction — §2 Modular Apéry systems and their sources (+ Theorem B exact)
— §3 Square roots of the sporadic rows and Beukers' Theorem 3 — §4 Classification of second-order
rows: Herfurtner's list — §5 Archimedean geometry and the single-row score (with census table) —
§6 Roots of $\Sym^w$ rows and the price of a free integration — §7 $p$-adic slopes and
extension-class rigidity (incl. the Euler-factor criterion, Conjecture D at $p=3$, good primes,
the weight drop) — §8 Two-row constructions: master formula and design rule — §9 Past the
threshold — §10 Problems — Appendix A verification scripts.

**Main results.**

- **Theorem A (Source Theorem)** `[PROVED]`. For each of the twelve modularly identifiable
  sporadic families $\mathbf A$–$\mathbf F$, $\alpha,\gamma,\delta,\varepsilon,\zeta,\eta$, the
  companion source $\Phi=F\,\theta_q t$ is exactly the holomorphic Eisenstein series of weight
  $r+1$ at the family's level; in particular every cuspidal projection is zero.
- **Theorem B / B\* (Apéry limit $=$ critical value)** `[PROVED]`. $\xi_\infty=\lim B_n/A_n=L(\Phi,w+1)$,
  and for an Eisenstein source $\Phi=P(V)E^{\psi,\varphi}_{w+2}$ satisfying the endpoint condition
  $\delta(\varphi)P(w+2)=0$ this is $P(w+1)L(\psi,w+1)L(\varphi,0)$ — so the census multiples
  $\tfrac16,\tfrac7{32},\tfrac58,\dots$ are **read off the Mellin polynomial**, with no
  case-by-case prefactors. Gives $\zeta(3)$ for $\gamma,\alpha,\varepsilon$; $\zeta(2)$ for
  $\mathbf A,\mathbf D$; $L(2,\chi_{-3})$, $L(3,\chi_{-3})$ for $\mathbf C,\mathbf F,\zeta$;
  Catalan $G$ for $\mathbf E$. Cuspidal sources give cusp-form $L$-values.
- **Beukers' Theorem 3, recovered** `[PROVED]`. Apéry's own square root is a *non-congruence*
  row; its period $\xi=L(\Psi,2)=0.10018744922933940616\ldots$, $\Psi=\tfrac14(\eta_1\eta_6)^{9/2}(\eta_2\eta_3)^{-3/2}$,
  is irrational with $\mu(\xi)\le50.66$. The paper is scrupulous that this **is** Beukers 1987
  rediscovered; what is new is the framing as the symmetric-square root, the explicit recurrence
  and Casoratian, the measure, and a formalization blueprint.
- **Proposition C (slope law)** `[PROVED]`. For a three-term row with constant term
  $c=\lambda_1\lambda_2$ the Casoratian is $-c^{n-1}/n^{w+1}$, so
  $\sigma_p=v_p(c)+2\kappa_p$; when $\sigma_p>0$, $b_n/a_n\to\xi_p$ in $\mathbb Q_p$ with
  $v_p(\xi_p-b_n/a_n)=\sigma_pn+O(\log n)$. **Product formula:** for integral rows
  $\sum_p\sigma_p\log p=\log|\lambda_1\lambda_2|$ — archimedean imperfection is redistributed as
  $p$-adic convergence.
- **Theorem F (Euler-factor criterion)** `[PROVED]` under four hypotheses. A $p$-adic Apéry limit
  exists iff the Euler factor $1-\psi(p)p^{-s}$ divides the source polynomial, and then
  $\xi_p=-Q(w+1)\kappa_p$ with $\kappa_p$ the Kubota–Leopoldt constant term. Fourteen census
  identifications hold to $\sim3000$ $p$-adic digits.
- **Conjecture D (rigidity)** — **now a theorem in the covered cases**, by Calegari's
  overconvergence method: Zagier's rows C and F at $p=3$ ($\xi_3=\tfrac12\zeta_3(2),\tfrac58\zeta_3(2)$),
  the Domb/$\varepsilon$ pair at $p=2$, the level-16 $\zeta(5)$ system. Row **B** is $[OPEN]$:
  the value $\tfrac12\zeta_3(2)$ is numerical to $3^{328}$, but the descent of the overconvergent
  function to $\mathbb P^1_{t_B}$ (B's Hauptmodul has degree 6 on genus-one $X_0(36)$) is not
  written. Notable: B has complex-conjugate characteristic roots so has **no** archimedean limit,
  yet its 3-adic limit coincides with C's — the $p$-adic limit survives where the archimedean one
  does not.
- **Good primes.** At good primes the extension splits along towers, tower limits are
  $\Gamma_p$-assemblies carrying no $L$-value, and no rigidity holds — a trichotomy read off the
  tower Frobenius eigenvalue $(\psi(p)p^{-w},1,0)$ across 144 cells.
- **Weight drop.** Cooper's rows are weight-one rows in disguise: their free integration
  $(n+1)\mid A_n$ *is* their source, $\Xi=\theta_q^{-1}\Phi$, and $\xi_\infty=L(\Xi,2)$.
- **Theorem E (two-row master formula)** `[PROVED modulo the selection theorem]`. Closed form for
  the quality $\delta$ of a correlated-lattice two-row construction, reproducing $0.857914$,
  $0.9025$, $0.6589$, the Lean-verified $1-\varepsilon$, and the new Domb$\times\varepsilon$ value
  $0.9010$ for $\zeta(3)$. **Design rule:** $\delta>1\iff\log\Lambda_{\rm dec}>k\max(r,1)+r\log\rho_2^{\rm eng}$;
  the binding constraint is a *decaying* partner with $p$-power denominators, which only
  hypergeometric rows supply — hence "cross-world pairs".
- **§9 (Catalan, past the threshold).** Nothing valuation-theoretic can exclude the
  rationality-kernel vector; the residue is Lemma P2, itself of irrationality strength.

**Evidence discipline.** 74 `[PROVED]`, 48 `[VERIFIED: range]`, 3 `[EXCLUDED]`, 11 `[OPEN]`, 3
`[LIT]` markers, plus a documented convention that these are load-bearing and never silently
upgraded. Seven `problem` environments. **Three `\todo` markers survive** and are honest debts:
(i) §2 the weight-three CM cusp-limit identification is 60-digit-verified only, with no written
source identity; (ii) §2 the $\beta(4)$ package's degree-four elimination / rank-16 cyclicity are
exact-computational certificates whose modular-function upgrade is unwritten; (iii) §8 the
$0.857914$ construction is reconstructed from the master formula because its original source is a
research conversation, not a paper.

**Dependencies.** No dependence on $\zeta_5(3)$ machinery (`\zeta_5` appears only in
`04_padic_euler_draft.tex` as a passing $p$-adic-zeta reference). Depends on external companion
project papers cited as `\cite{ProjSources,ProjCatalanLCM,ProjCatalan58,LeanCatalan,CatalanLean}`
— **these are not in this repo tree**, and that is the paper's chief external liability.
`sections/04_padic_rigidity_draft.tex` is orphaned (superseded by `04_padic_conjD_proof.tex`).

**Judgment — SHAREABLE NOW as a draft, near-submittable.** The claim envelope is explicitly and
correctly bounded ("No new irrationality theorem"), the evidence grading is exemplary, the
computations are reproducible (I re-ran the appendix's two `verify/` scripts: 11/11 + 1/1 PASS),
and the LaTeX is clean. To submit: resolve the three `\todo`s (or downgrade the claims they
guard), supply or inline the cited project papers, and either finish or explicitly park row B.

### 2.2 `paper/classification/main.tex` — **Integral second-order Apéry-like recurrences and rational elliptic surfaces**

*= `consolidation/classification_paper.pdf`. 21 pp. Author line: "Claude (Fable), with River".
Self-contained, inline bibliography. Compiles: 21 pages, 4 warnings, no errors.*

The most self-contained and most conventionally publishable item in the repo. It answers Zagier's
question — for which $(P,Q)$ does $(n+1)^2u_{n+1}=P(n)u_n-Q(n)u_{n-1}$ have all $u_n\in\mathbb Z$?
— in two tiers.

- **Theorem A (conditional classification)** `[PROVED]`, given Beauville's and Herfurtner's
  classifications. For the Zagier shape $(an^2+an+b,\;cn^2)$: if $\ker L$ is, up to rational gauge
  and quadratic twist, $R^1\pi_*\mathbb Q$ of an elliptic surface with non-constant $\mathcal J$,
  then (a) all four exponent differences vanish, so all four singular fibres are semistable;
  (b) $\deg\mathcal J\le12$ by Riemann–Hurwitz **and** $12\mid\deg\mathcal J$ by a twist-parity
  lemma, so $\deg\mathcal J=12$ exactly and $\mathcal J$ is Belyi — hence the untwisted member is
  a *rational* elliptic surface, one of Beauville's six. **No Kodaira-dimension hypothesis is
  needed**, because a quadratic twist raises $\chi$ (two branch points give a K3) without changing
  $\mathcal J$, the local system, or the row. (c) The row is the period expansion at a cusp in the
  integral coordinate; (d) the cross-ratio $\mathcal I=a^2/c$ lies in the explicit finite set
  $\{0,3,\tfrac92,-\tfrac{49}8,\tfrac{100}9,\tfrac{289}{72},-121\}$. The complete list of cusp
  placements is Zagier's seven rows $\mathbf A$–$\mathbf G$, with $\mathbf A,\mathbf C,\mathbf F$
  the three placements of $I_6I_3I_2I_1$ and $\mathbf E,\mathbf G$ the two of $I_1I_1I_2I_8$.
  The analogue for the root-row shape $(An^2+\tfrac A2n+B,\;C(2n-1)^2)$ gives $\deg\mathcal J=6$,
  three surviving Herfurtner configurations, exactly two integral rows.
  Beukers' row (the square root of Apéry's $\zeta(3)$ recurrence) is **not** among them: it lives
  on $X_0(6)/W_6$, over which the universal curve does not descend.
- **Theorem B (unconditional necessary conditions)** `[PROVED]`, no geometric hypothesis:
  normalisation ($\rho_1=\rho_2$ when $a^2-4d\notin\mathbb Q^{\times2}$); Kodaira admissibility
  (every local exponent difference in $\{0,\tfrac12,\tfrac13,\tfrac23\}\bmod1$); **integrality
  rigidity** ($\lambda_1\in\mathbb Q\iff\lambda_1,\lambda_2\in\mathbb Z\iff a^2-4d$ square, whence
  $\mathrm{score}\le-k\le-1$ on every rational cusp-move orbit — so a positive-score row has
  irrational characteristic roots); and the Riemann–Hurwitz bound
  $\deg\mathcal J\le6c+4e_3+3e_2-12\le12$.
- **Theorem C (the scans)** `[VERIFIED: boxes as stated]`. Exhaustive integrality scan of all 26
  Kodaira-admissible normalisation classes, $1\,270\,065$ hits, integrality prime-tested to $n=30$
  and re-verified exactly to $n=300$: 28 primitive non-degenerate rows in the rigid classes, of
  which exactly ten have an elliptic-surface local system — Zagier's six, the two root rows, and
  **two new rows** on $I_1I_7\,II\,II$ (monodromy $\Gamma_0(7)$) and on $I_3\,III\,III\,III$ (the
  latter with Chowla–Selberg period $\operatorname{Im}\xi=-\Gamma(1/4)^4/(2^{7/2}3^{11/4}\pi)$).
  The only rows anywhere in the scan with positive score and $k\ge2$ are Apéry's $\zeta(2)$ row
  $(+0.40606)$ and Beukers' row $(+0.13920)$. An independent modular scan of $3\,036\,854$
  eta-quotient pairs on $\Gamma_0(N)$, $N\le60$, gives exactly two classes with $|\lambda_2|<1$.

Tags: 25 `[PROVED]`, 24 `[VERIFIED]`, 8 `[NUMERICAL]`, 9 `[OPEN]`. No $\zeta_5(3)$ or Catalan
dependence. The unconditional classification is explicitly identified as Zagier's conjecture
(Bombieri–Dwork in rank two with four singular points), with the reduction made precise.

**Judgment — SHAREABLE NOW; the single best candidate for external release.** Self-contained,
correctly attributed (Zagier §7 + Beauville; new content is the mechanism, the twist-parity
lemma, the cusp-move theorem, the completeness, and the two new rows), inline bibliography, clean
build, honest tags. Two things a referee would want: the scans' scripts packaged reproducibly
(they live scattered in `lattice/`), and a decision on whether the 8 `[NUMERICAL]` items are
load-bearing.

### 2.3 `paper/four_term/main.tex` — **Apéry-like recurrences beyond the rigid case: four-term rows, K3 surfaces, and new periods**

*= `consolidation/four_term_paper.pdf`. 19 pp, "Claude (Fable), with River". Compiles: 19 pages,
21 warnings, no errors.* Companion/sequel to §2.2, moving from four singular points to five (two
accessory parameters). **Only partly Catalan-adjacent** — most of it is general classification.

Theorems A (dictionary and normalisation classes, equal- and mixed-exponent), B (the census),
C (no arithmetic gain), plus structure theorems D1 (a decaying row needs an irreducible
characteristic cubic), D2 (a real decayer needs a totally real one), D3 (mixed normalisation
forms), D5 (five-term dictionary), and Conjecture D6. Riemann–Hurwitz bound
$\deg\mathcal J\le6c+4e_3+3e_2-12\le18$. Four findings: **near-total rigidity persists** (of 2178
integral rows, 98.2% have an apparent singularity and are gauge/cusp-move images of three-term
rows); **a finite sporadic list appears** (30 primitive genuine five-point rows, 28 with $k=2$);
**K3 geometry enters exactly here** (one row is the PF system of an elliptic K3 with five fibres
$I_4I_6I_3I_3IV^*$, impossible with four points); **new periods appear** (mixed-exponent rows with
Apéry limits $\tfrac14G$, $\tfrac12G-\tfrac3{16}\zeta(2)$, $\tfrac38\zeta(2)-\tfrac12G$,
$L(2,\chi_{-3})$-combinations, gauge partners at $\Gamma(\tfrac14)^4/(64\pi)$, and a six-point
five-term row with limit $\tfrac{\sqrt2}3L(g,2)$). Negative headline: **no new
irrationality-relevant arithmetic** — no genuine five-point row has positive score, no
unequal-exponent class can contain a decayer, and the 2-adic content of the new Catalan rows is
the ledger's own. Tags: 27 `[PROVED]`, 3 `[NUMERICAL]`, **15 `[OPEN]`**.

**Judgment — NEEDS WORK (light).** Paper-shaped and compiling, but the highest `[OPEN]` density of
any item here, and the census is box-complete rather than complete. Best released after §2.2 as
its explicit sequel.

### 2.4 `paper/k3_row/main.tex` — **An Apéry-like sequence attached to a singular K3 surface**

*= `consolidation/k3_row_note.pdf`. **3 pp**, "Claude (Fable), with River", 23 August 2026 —
draft. Compiles clean, zero warnings.* Filed as Catalan-adjacent but **it is not Catalan**; it is
the sharpest self-contained novelty claim in the repo.

The sequence $1,4,16,64,250,928,3136,8704,11866,-79400,\dots$ from
$(n+1)^2u_{n+1}=(11n^2+11n+4)u_n-(37n^2+3)u_{n-1}+3(3n-1)(3n-2)u_{n-2}$ is the period sequence of
an elliptic K3 with five singular fibres $I_4,I_6,I_3,I_3,IV^*$ and Picard number 20. Its
transcendental motive is that of the weight-3 CM newform `32.3.d.a` (CM by $\mathbb Q(\sqrt{-2})$),
and its Apéry constants are critical $L$-values of that form: fold constants
$\tfrac{\sqrt2\mp i}3L(g,2)$ at the conjugate singularities $\tfrac{5\mp i\sqrt2}{27}$, and
$\tfrac{2\sqrt2}3L(g,2)$ at $t=1$, with
$L(g,2)=\tfrac\pi{32}\Gamma(\tfrac18)\Gamma(\tfrac38)/\Gamma(\tfrac58)\Gamma(\tfrac78)$ by
Chowla–Selberg. **Claimed novelty:** the first Apéry-like sequence whose Apéry constant is a
critical value of a cusp form — with the reason three-term rows can never produce one being
precisely the Riemann–Hurwitz/twist-parity bound of §2.2.

**Judgment — SHAREABLE NOW as a short note, after a novelty sweep.** Three pages, one crisp
claim, complete argument chain, and it is the natural "teaser" for the classification paper. It
contains only one `\begin{proof}` environment across three sections, so a referee would want the
period computation written out. The novelty claim ("to our knowledge the first") is the one thing
that must be audited against the Apéry-like literature before release.

### 2.5 Catalan-adjacent (one line each, per instruction)

- **`paper/catalan_1eps/`** = `consolidation/catalan_1eps_proof.pdf`, **4 pp**: *The near-critical
  two-row construction for Catalan's constant: a self-contained proof of the worthiness
  $1-\varepsilon$ theorem* — for every $0\le k<k^*=22.3513\ldots$, combinations of Zudilin's and
  Nesterenko's approximations over $D_{6n}^2 2^{\lfloor kn\rfloor}$ give worthiness
  $\delta=D/(D+F(k))>1-\varepsilon$ with $F(k)=\tfrac{\log2}2(k^*-k)$; matches the Lean project
  `lean-catalan-worthiness`. Complete and compiling; explicitly says what it does *not* prove.
- **`paper/catalan_cdt/`** = `consolidation/catalan_cdt_obstruction.pdf`, **7 pp**, draft v2:
  *Catalan's constant and the arithmetic holonomy method: an obstruction* — every genus-zero
  modular host of the $\chi_{-4}$ Eisenstein class has $|t_2|\le\tfrac14$, so the CDT entry
  condition fails by $\ge0.077$ nats/index even at the uniformisation ceiling; the 2-adic repair
  is blocked by Calegari's theorem; AL-quotient hosts put the fold at an orbifold point; Hadamard
  two-host constructions admit no pure module; EMN's realisation clears entry only
  self-defeatingly. The ledger wants a v3 folding in the PTPAH closure.
- **`consolidation/zeta5_3_draft.pdf`** (6 pp, 23 Aug) is **superseded** by
  `paper/zeta5_3/main.pdf` (**7 pp**, 27 Aug, "corrected draft for audit", compiles with zero
  errors and zero warnings, last touched by commit `e2aa7ec` "resolve error in 4.3"). Released
  flagship; content not re-surveyed here.

### 2.6 `consolidation/paper_draft0..8.pdf` — the draft series

All nine are the **same evolving paper** as §2.1, confirmed by page-1 extraction: drafts 0 and 4
open with the identical title *MODULAR APÉRY SYSTEMS: SOURCES, SLOPES, AND THE TWO-ROW DESIGN
RULE* and identical §1 opening; draft 8 adds the abstract. Growth 26 → 26 → 28 → 31 → 37 → 39 →
44 → 45 → **51** pages over 21–22 August. `paper_draft8.pdf` is md5-identical to
`paper/main.pdf`. **Drafts 0–7: SUPERSEDED — safe to delete or archive.**

---

## 3. Artefact directories

| Directory | Size | Contents | Verdict |
|---|---|---|---|
| **`verify/`** | 12 K | Two self-contained checkers with run instructions in their docstrings: `period_annihilation.gp` (PARI/GP, exact rational arithmetic, ~11 period-annihilation facts — Beukers $\zeta(3)$ level 6, level-8/16 projectors, the level-12 $\zeta(7)$ parent hitting $\tfrac{209}{1728}$, $\beta(4)$ inner/outer vectors, the $\chi_{-3}$ purified vector $(1,-10,16)$) and `beta4_lambert.py` (mpmath, `mp.dps=40`, $\beta(4)$ Lambert CM identity). **I re-ran both: 11/11 PASS and 1/1 PASS (error $4.6\times10^{-41}$).** | **REUSABLE ARTEFACT** — small, exact, self-verifying, independently confirmed working today. |
| **`lattice/`** | 2.5 G | The main compute programme: 815 `.gp`, 255 `.py`, 809 `.log` across 60 subdirs (36 with their own README), exact integer/rational PARI/GP throughout. Covers the $p$-adic census, Euler criterion ($\ge p^{2991}$ tables), Zagier/Domb alignments, four-term scans, the P2$'$ work, CDT finder, non-congruence scans. 3380 files git-tracked. Size is dominated by one file: `p2_scale/data/rows_scale.txt` (2.4 G), correctly gitignored and documented as regenerable in 474 s via `build_rows.gp`. Logs contain genuine `FAIL`/`FIRST FAILURE` lines that are themselves findings, not breakage. | **REUSABLE (with friction)** — real, documented, re-runnable exact arithmetic; but 60 semi-independent experiment folders and no single unified test entry point, so reproducing a *specific* claim means finding the right script and prepend chain. This is the biggest under-packaged asset in the repo. |
| **`certificates/`** | 1.7 M | 32 `.tex` documents pulled from ChatGPT sessions, filed by strand (`eisenstein/`, `catalan/`, `harmonic/`, `padic/`, `misc/`), with `INDEX.md` recording per-document "Ran? latexmk / Result: PASS, N pp." and `WANTED.md` a to-do list with a liaison workflow. `documents/` is a byte-identical flat mirror of the strand files. **No PDFs or compile logs are stored** — the PASS markers are recorded claims, not local evidence. | **PARTIALLY REUSABLE** — valuable dated provenance index; nothing here is independently re-checkable as it stands. |
| **`tools/`** | 108 K | Two unrelated groups: `build_lineage_index.py` (stdlib-only, reconstructs branch lineage across the ChatGPT archive into `lineage_index.{json,md}`) with `SNAKE.md` (remote PARI/GP compute box: SSH, rsync patterns for offloading long `lattice/` jobs) and `summarization_spec.md`; plus `chatgpt-archiver-extension/`, a full unpacked Chrome extension for archiving conversations. | **MIXED / infrastructure** — useful tooling, not a mathematical artefact. |
| **`catalan-2-row-denominators/`** | 3.0 M | Python/mpmath research scratch on the two Catalan rows: `catalan_rows.py` (checkpointed high-precision row generation), plot scripts + PNGs, and a genuine data set (`zudilin_rows.csv`, `nesterenko_rows.csv`, `combined_rows.csv`, `catalan_arithmetic_diagnostics.csv`). `5-8 theorem/` holds a self-contained LaTeX proof **with a compiled PDF** — real build evidence, unlike `certificates/`. `zudilin note.txt` is a candid strategy memo flagging that $\delta\to1$ is not an irrationality proof. No top-level README. | **PARTIALLY REUSABLE** — the CSVs are a usable database with regenerating scripts present, but the directory is informally organised active scratch. |
| **`catalan_rows/`** | 4 K | **Empty and untracked.** Superseded by `catalan-2-row-denominators/`. | **SCRATCH — delete.** |
| **`packages/`** | 3.0 M, untracked | Informal named batches ("phase 0 math campaign", "phase 1 zeta math package", "unknown packages") of raw research notes (23 `.txt`), 6 `.tex` drafts, `Odd_Zeta_Research_Master_Ledger.xlsx`, and archives. `certificates/WANTED.md` flags a mis-filed Lean tarball here. Source of the `CLAUDE_FINDINGS.md` $\zeta(5)$ asymptotics recovery. | **SCRATCH / raw source material** — mine it, don't ship it. |
| **`papers/`** | 28 M | External literature only (see §2 preamble). | **Reference material.** |

---

## 4. Ranked list — most shareable outputs (excluding $\zeta_5(3)$ and Catalan)

1. **`paper/classification/` — *Integral second-order Apéry-like recurrences and rational
   elliptic surfaces* (21 pp).** Ready now. A clean answer to a question of Zagier's, with a real
   new mechanism (twist parity forcing $12\mid\deg\mathcal J$, so no Kodaira-dimension hypothesis
   is needed), a completeness result, two new rows, and a $1.27\times10^6$-hit exhaustive scan.
   Self-contained with inline bibliography, honest evidence tags, compiles clean. Ship first.

2. **`paper/main.tex` — *Modular Apéry systems: sources, slopes, and the two-row design rule*
   (51 pp).** The programme's synthesis, and the paper the book was trying to be. Theorem B\*
   (Apéry limit *is* the critical value, prefactor read off the Mellin polynomial), the slope law
   and its product formula, the Euler-factor criterion, and Conjecture D proved in named cases are
   each independently interesting. Blockers before submission: three surviving `\todo` debts,
   citations to five companion project papers that are not in this tree, and the unfinished
   descent for Zagier's row B.

3. **`paper/k3_row/` — *An Apéry-like sequence attached to a singular K3 surface* (3 pp).** The
   highest novelty-per-page item here: the first Apéry-like sequence whose Apéry constant is a
   critical cusp-form value ($L(g,2)$ for the CM newform `32.3.d.a`), with the classification
   paper's twist-parity bound explaining why three-term rows never do this. Needs one literature
   sweep on the novelty claim and the period computation written out, then it's a letter.

*Runners-up:* `verify/` (small but genuinely reusable and confirmed working — it should be grown
into the reproducibility backbone for items 1–3); `paper/four_term/` (19 pp, paper-shaped, but
15 `[OPEN]`s and best released as the explicit sequel to item 1); the `lattice/` suite (a large
real asset that needs a single documented entry point per published claim before it can be cited).

*Not on the list, and why:* the **book** (framework outruns its theorems; v11 unmerged;
certificate-level headline claims) and the **`CLAUDE_FINDINGS.md` $\zeta(5)$ asymptotic** (a
genuinely nice exact result — $d_n/c_n-C_5\sim\tfrac{13}{96K}(-1)^n(7-4\sqrt3)^nn^{1/2}(\log n)^4$
— that is stranded in a findings log, has one undetermined constant $K$, and appears in no paper;
worth either closing $K$ via the discriminant $-48$ CM period or folding into item 2).
