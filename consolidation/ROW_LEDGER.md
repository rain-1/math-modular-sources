# Row ledger — census of Apéry-like rows scored for irrationality potential

*Computed 2026-08-21. Scripts: `lattice/census/master12.gp` (this session, exact PARI,
verified for the 12 order-2/order-3 sporadics), reusing prior verified work in
`consolidation/SLOPE_CENSUS.md`, `consolidation/THEORY_NOTES_03_lattices.md`,
`consolidation/ZETA3_TWO_LATTICE.md`, `CLAUDE_FINDINGS.md`. Every number below is
tagged `[proved]`, `[verified to n=N]`, or `[numerical]`. Nothing here is guessed;
where a limit could not be identified this is stated as "unidentified", not filled in.*

## Definitions used

Order-2 normalisation: $(n+1)^2u_{n+1}=(an^2+an+b)u_n-cn^2u_{n-1}$, $a_0=1,a_1=b,b_0=0,b_1=1$.
Order-3: $(n+1)^3u_{n+1}=(2n+1)(an^2+an+b)u_n-cn^3u_{n-1}$, same initial data.
$\lambda_1,\lambda_2$: roots of $x^2-ax+c$ (order 2) or $x^2-2ax+c$ (order 3) — the
$n\to\infty$ char. equation of the recurrence's leading coefficients (exact, [proved]
trivially from the recurrence shape). $c=\lambda_1\lambda_2$; $\sigma_p=v_p(c)$.
$k$ = smallest exponent with $d_n^k b_n\in\mathbf Z$ ($d_n=\mathrm{lcm}(1,\dots,n)$),
checked by exhaustive integrality scan.
**score** $=\log(1/|\lambda_2|)-k$ (row alone proves irrationality of its period iff $>0$).
**budget** $=\log\lambda_1-k$ (best conceivable score after full $p$-adic harvesting).
**headroom** $=$ budget $-$ score $=\sum_p\sigma_p\log p$ (order-2/3 two-root product case only;
flagged separately where the row has more than 2 characteristic roots).

## Table 1 — the twelve order-2/order-3 sporadics (this session, exact PARI, `master12.gp`)

All six order-2 and six order-3 rows below were run to $n=300$ (order 2) / $n=220$
(order 3), primes $p\in\{2,3,5,7,11,13\}$. In every case the *only* prime with
positive measured slope is the one predicted by $v_p(c)$ — this reproduces
`THEORY_NOTES_03` §1 exactly and extends it to $p=11,13$ (all zero, as expected since
$11,13\nmid c$ for every row here). $k=2$ (order 2) / $k=3$ (order 3) for all twelve,
confirmed integral to the stated $n$ — no exceptions, no free denominator savings anywhere.

| row | $(a,b,c)$ | $\lambda_1$ | $\lambda_2$ | limit (basis, `lindep`) | $\sigma_2,\sigma_3,\sigma_5,\sigma_7$ | $k$ | score | budget | headroom |
|---|---|---|---|---|---|---|---|---|---|
| Zagier A | $(7,2,-8)$ | 8 | $-1$ | $\zeta(2)/4$ [verified n=280/300, 2 indep. truncations agree] | $3,0,0,0$ | 2 [verified n≤300] | $-2.000$ | $+0.079$ | $+2.079$ |
| Zagier B | $(9,3,27)$ | $|\lambda|=\sqrt{27}$ (complex) | — | **no real archimedean limit** (complex roots) [proved: discriminant $<0$] | $0,3,0,0$ | 2 [verified] | n/a (complex) | $-0.352$ | n/a |
| Zagier C | $(10,3,9)$ | 9 | 1 | $L(2,\chi_{-3})/2$ [verified n=280/300] | $0,2,0,0$ | 2 [verified] | $-2.000$ | $+0.197$ | $+2.197$ |
| Zagier D | $(11,3,-1)$ | $11.0902$ | $-0.0902$ | $\zeta(2)/5$ [verified n=280/300] | $0,0,0,0$ | 2 [verified] | $+0.4061$ | $+0.4061$ | $\approx0$ |
| Zagier E | $(12,4,32)$ | 8 | 4 | $G/2$ (Catalan) [verified n=280/300] | $5,0,0,0$ | 2 [verified] | $-3.386$ | $+0.079$ | $+3.466$ |
| Zagier F | $(17,6,72)$ | 9 | 8 | $\tfrac58L(2,\chi_{-3})$ [numerical, 16 digits agree at n=300; `lindep` search up to n=900/350 digits found no low-height relation because decay is only $(8/9)^n$ — confirmed instead by direct comparison] | $3,2,0,0$ | 2 [verified] | $-4.079$ | $+0.197$ | $+4.277$ |
| AZ $\delta$ | $(7,3,81)$ | $|\lambda|=9$ (complex) | — | **no real limit** (complex, on Beauville-IV surface w/ Domb/T/Apéry per archive; 3-adic limit not aligned w/ Domb/T/Apéry at any of 22 rational scalars tested [SLOPE_CENSUS]) | $0,4,0,0$ | 3 [verified] | n/a | $-0.803$ | n/a |
| AZ $(9,3,-27)$ | $(9,3,-27)$ | $19.392$ | $-1.392$ | $L(3,\chi_{-3})/3$ [verified n=200/220] | $0,3,0,0$ | 3 [verified] | $-3.331$ | $-0.035$ | $+3.296$ |
| Domb | $(10,4,64)$ | 16 | 4 | $\tfrac7{24}\zeta(3)$ [verified n=200/220] | $6,0,0,0$ | 3 [verified] | $-4.386$ | $-0.227$ | $+4.159$ |
| AZ $\eta$ | $(11,5,125)$ | $|\lambda|=\sqrt{125}$ (complex) | — | **no real limit** (complex); $=\tfrac12L(3,\chi_5)$ per book `08c_chi5_application.tex` [source-cited, not re-derived] | $0,0,3,0$ | 3 [verified] | n/a | $-0.586$ | n/a |
| T | $(12,4,16)$ | $23.314$ | $0.6863$ | $\tfrac7{32}\zeta(3)$ [verified n=200/220] | $4,0,0,0$ | 3 [verified] | $-2.624$ | $+0.149$ | $+2.773$ |
| Apéry | $(17,5,1)$ | $33.971$ | $0.02944$ | $\zeta(3)/6$ [verified n=200/220] | $0,0,0,0$ | 3 [verified] | $+0.5255$ | $+0.5255$ | $\approx0$ |

## Table 2 — other rows (from prior verified work in this repo; not re-derived this session unless noted)

| row | order | period | $\lambda_1$ | $\lambda_2$(-ish) | $k$ | score | budget | headroom | notes |
|---|---|---|---|---|---|---|---|---|---|
| Cooper $s_{10}=\sum\binom nk^4$ | 2 | $\zeta(2)/5$ | 16 | $-4$ | **undetermined** ($A_n\notin\mathbf Z$, denominator grows; the usual $d_n^kb_n$ test doesn't directly apply) [SLOPE_CENSUS] | n/a | n/a | n/a | measured $p$-adic slope **0** at $p=2,3,5,7$ [verified n≤500] — arithmetically inert despite sharing $\zeta(2)/5$ with Zagier D |
| Cooper $s_7=\sum\binom nk\binom{2k}k\binom{2n-2k}{n-k}/\text{(Cooper's def.)}$ | 2 | $\zeta(2)/7$ | 27 | $-1$ | undetermined | n/a | n/a | n/a | slope 0 at $p=2,3,5,7$ [verified n≤500]; $|\lambda_2|=1$ so even archimedeanly it is at best borderline |
| Cooper $s_{18}$ | 2 | $\tfrac12L(2,\chi_{-3})$ | 16 | 12 | undetermined (but row happens to stay integral through $n\le29$) | n/a | n/a | n/a | **only** Cooper row with positive slope: $\sigma_3=1$ [verified n≤500, matches $v_3(192)=1$]; aligned with Zagier B, C at $p=3$ [verified] |
| ζ(5), level 12 (CM-isogeny, Sym$^4$ of $\lambda^2-14\lambda+1$) | 5 | $\tfrac{11}{144}\zeta(5)$ | $\alpha^4\approx37634.0$, $\alpha=7+4\sqrt3$ | second-dominant root $\alpha^2\approx193.99$ | $k=5$ [as stated by task/archive source; not independently re-derived here] | $-10.27$ (using 2nd-dominant-root decay law) | $+5.536$ | $+15.80$ **nominal, but see caveat** | $c=\prod_{i=1}^5\lambda_i=\alpha^4\cdot\alpha^2\cdot1\cdot\alpha^{-2}\cdot\alpha^{-4}=1$ exactly [proved from stated char. poly], so **true $\sum_p\sigma_p\log p=0$ — all slopes are ZERO**, exactly as the task brief states. The naive headroom $=\log\alpha^2=5.27$ above is *not* $p$-adically harvestable; it comes from the three intermediate roots $(1,\alpha^{-2},\alpha^{-4})$ of this order-5 system, not from a prime. This is flagged as an open structural point: the two-root "budget=headroom+score" identity used elsewhere is a genuine order-2 fact and does **not** transfer to $\mathrm{Sym}^{>2}$ rows without modification. Slopes at $p=2,3,5,7$: **not measured** (order-14 recurrence not transcribed/run) [caveat honestly stated in SLOPE_CENSUS]. |
| ζ(7), level 24 (Sym$^6$, eta-quotient $r(\tau)$) | 7 | $\tfrac{1463}{13824}\zeta(7)$ | $\lambda_1=4+2\sqrt2\approx6.828$ | $\lambda_2=\lambda_1/\sqrt2\approx4.828$ (2nd-dominant, boxed rate $1/\sqrt2$ in source paper) | 7 [claimed sharp by source paper; **verified numerically** $d_n^7B_n\in\mathbf Z$ for all $n\le218$] | $-8.575$ | $-5.079$ | $+3.496$ **but see caveat** | Measured slopes at $p=2,3,5,7$: **all zero/no growth**, checked to $n=218$ [verified] — this row is "Apéry-perfect": nothing to harvest at these primes; the naive headroom above is again an order$>2$ artifact from intermediate roots, not real $p$-adic budget. Limit convergence $B_n/A_n\to\tfrac{1463}{13824}\zeta(7)$ verified to 33 digits at $n=218$, decay rate matches $1/\sqrt2$ to 4 sig figs. |
| Zudilin Catalan row ($Q_m,P_m$) | 2, non-Zagier normalisation | $G$ | — (not recomputed this session) | $\lambda_2=\varphi^{-5}<1$ [source: `CATALAN_0902526_SHORT_PROOF.tex`] | — | — | — | — | Not independently re-derived here; used only via the already-published pair with Zagier-type modular row (see Table 3). |

## Table 3 — pairs with shared period and positive slope at a common prime; two-lattice quality $\delta$

Only two configurations in the whole census satisfy the pairing precondition
(shared period $\Theta$, both rows with $\sigma_p>0$ at a common $p$):

### (i) Domb + T, period $\zeta(3)$, $p=2$ — **already computed exactly** in `ZETA3_TWO_LATTICE.md`

Balanced sampling $r=\sigma_{\rm dec}/\sigma_{\rm eng}=6/4=3/2$ (Domb engine at $2n$, T
decayer at $3n$), $k=3$:
$$\delta = \frac{12\log\Lambda_T-24\log2}{9+9\log\Lambda_T-20\log2}=0.9009531686558563\ldots \quad[\text{proved algebra, conditional on }\varepsilon=0\text{ (§3 "Open" of that file)}]$$
$\delta<1$: **not** an irrationality proof, quality is $\approx10\%$ short. Design-rule
check (§5.1 shopping-list table, $k=3,\sigma_{\rm eng}=6,\sigma_{\rm dec}=4$ cell): need
$\Lambda_{\rm dec}>50.6$; T has $\Lambda_T=23.31$ — **fails** the fast filter by roughly a
factor of 2, consistent with $\delta<1$.

### (ii) B, C, F (Zagier) + Cooper $s_{18}$, period $L(2,\chi_{-3})/2$-family, $p=3$

All four share $p=3$ slope: $\sigma_3(B)=3$, $\sigma_3(C)=2$, $\sigma_3(F)=2$, $\sigma_3(s_{18})=1$
[verified, this session for B/C/F, `SLOPE_CENSUS.md` for $s_{18}$]. **But the two-lattice
architecture requires one row to play the "decayer" role with $|\lambda_2|<1$, and
*none* of these four does**: $C$ has $\lambda_2=1$ exactly (only polynomial decay, not
exponential), $F$ has $\lambda_2=8$ (growing), $B$ and $s_{18}$ have $|\lambda_2|\ge\sqrt{27}$
and $12$ respectively (growing). **The design-rule fast filter cannot even be
evaluated** ($\log\Lambda_{\rm dec}$ has no candidate row) — this family is **ruled out**
for the Domb/T-style construction entirely, independent of any $\delta$ computation.
No $\delta>0.9$ pair exists here; none was computed since the precondition (a genuine
decayer) fails.

No other pair in the census shares both a period and a common positive-slope prime:
Cooper $s_7$/$s_{10}$ align with nothing (slope 0 everywhere tested); AZ $\delta$ aligns
with nothing at $p=3$ (22 rational scalars tried, all failed [SLOPE_CENSUS]); AZ
$(9,3,-27)$'s period $L(3,\chi_{-3})/3$ has no second row in this census; AZ $\eta$'s
period $\tfrac12L(3,\chi_5)$ has no second row; the level-12/level-24 higher rows have
zero measured slope at every tested prime so cannot pair with anything by this
criterion.

## What this rules out / what survives

**Ruled out (this census):**
- No single row in Tables 1–2 has score $>0$ **and** a genuine exponential-decay
  linear form combined with nonzero $p$-adic slope — the only rows with score $>0$
  (Zagier D at $+0.406$, Apéry at $+0.526$) are exactly the two "Apéry-perfect" rows
  with $c=\pm1$ and **zero** slope: their positive score is pure archimedean decay,
  already known (D gives nothing new — $\zeta(2)$ is trivially rational-adjacent
  known-irrational via other means; Apéry **is** the classical theorem). Nothing new
  is proved by any single row here.
- The only two candidate *pairs* satisfying the shared-period/shared-slope
  precondition both fail: Domb+T gives $\delta=0.9010<1$ (proved algebra, conditional
  on one open $\varepsilon=0$ claim); the $L(2,\chi_{-3})/2$ family (B,C,F,$s_{18}$) has
  no viable decayer at all, so the construction cannot even be attempted.
- Cooper's $s_7,s_{10}$ and AZ's $\delta$ carry **no p-adic slope at any tested small
  prime** despite sharing periods/geometry with rows that do — confirming the
  "extension-class rigidity" law of `THEORY_NOTES_03` (alignment needs the *same*
  extension class, not just the same period).
- The level-12 $\zeta(5)$ and level-24 $\zeta(7)$ rows are (as far as tested, $p\le7$)
  "Apéry-perfect": all decay is archimedean, nothing $p$-adically harvestable — same
  situation as Apéry's classical row, not the Catalan/Domb-type split-budget rows.

**Survives / open:**
- Apéry's own row remains the only unconditional single-row proof in the whole
  census (of course — that is the 1979 theorem, not a new finding).
- The Domb+T $\zeta(3)$ construction is a complete, honest calibration: it is proved
  conditional on one 2-adic constant-equality claim ($\varepsilon=0$) that is itself
  believed provable via Eichler-integral/Eisenstein-series methods but **not done**;
  even if closed, $\delta=0.901<1$ so it would not yield a new proof of Apéry's
  theorem, only validate the machine.
- The design rule's own "cheapest lever" observation ($k{=}2,\sigma{=}3{:}3$ engines
  with slow-growing $\rho_2^{\rm eng}=2$ need only $\Lambda_{\rm dec}>14.8$) is not
  instantiated by any row found in this census — no row here has both slope 3 at a
  shared prime *and* a genuine $\lambda_2<1$ partner at that prime. This is the
  concrete open target the ledger points to for future search.
- $\zeta(5)$/$\zeta(7)$ higher-Sym rows: slopes at $p\le7$ genuinely not measured for
  the full order-14/order-7 systems (only the Sym$^4$/Sym$^6$ eigenvalue structure and
  a truncated $q$-series check of the ratio were done) — an honest gap, not a claimed
  negative result, for primes/structure beyond what was checked.

## Scripts
- `lattice/census/master12.gp` — this session: exact integer models for all twelve
  order-2/order-3 sporadics, denominator-exponent scan, slopes at $p\in\{2,3,5,7,11,13\}$,
  `lindep` limit identification with cross-truncation stability check (basis pruned to
  remove the redundant $\pi^2\equiv6\zeta(2)$ relation that was corrupting `lindep` in
  an earlier draft of this script — flagged and fixed in-session).
- `lattice/census/cooper.gp`, `cooper_align.gp`, `az_families.gp` — prior session,
  Cooper rows and the complex-root AZ family (reused, not rerun).
- `lattice/zeta3_lattice/rows.gp`, `consolidation/ZETA3_TWO_LATTICE.md` — prior
  session, the Domb+T exact $\delta$ computation (reused, not rerun).

---

## CORRECTION (audit by the parent agent, same session)

Two entries above are wrong and one definition needs generalising. The source is
`packages/phase 1 zeta math package/modular_apery_cm_isogeny_research_notes.txt`,
which states the characteristic polynomials **explicitly**; they were not read.

### C1. The order-$m$ product formula

The identity `budget = score + headroom` used above is an **order-2 fact**. In general,
for a scalar recurrence with characteristic roots $\lambda_1,\dots,\lambda_m$,

$$\sum_p\sigma_p\log p=\log\Bigl|\prod_{i=1}^m\lambda_i\Bigr|\qquad\text{(not }\log|\lambda_1\lambda_2|\text{)},$$

so the correct definitions for any order are

$$\textbf{score}=\log\tfrac1{|\lambda_2|}-k,\qquad
\textbf{harvestable budget}=\textbf{score}+\log\Bigl|\prod_i\lambda_i\Bigr| .$$

For $m=2$ this reduces to $\log|\lambda_1|-k$ as used elsewhere. **Using
$\log\lambda_1-k$ for a higher-order row overstates the budget.**

### C2. ζ(5), level 12 (CM-isogeny) — row corrected

The ledger models this as $\mathrm{Sym}^4$ of $\lambda^2-14\lambda+1$ with dominant root
$\alpha^4\approx37634$ and budget $+5.536$. The source (line 511) gives the characteristic
polynomial of the homogeneous recurrence as

$$(\lambda+1)^4(\lambda^2-14\lambda+1)^5,\qquad\text{roots }-1,\ 7-4\sqrt3,\ 7+4\sqrt3,$$

"which are exactly the reciprocals of the cusp/branch singular values", and for the joint
$c_n,d_n$ recurrence $(\lambda+1)^5(\lambda^2-14\lambda+1)^5=(\lambda^3-13\lambda^2-13\lambda+1)^5$.

$$\lambda_1=7+4\sqrt3=13.9282,\quad\lambda_2=-1,\quad\textstyle\prod_i\lambda_i=-1 .$$

- $\log|\lambda_1|=2.63392$, **not** $\log(37634)=10.535$. The dominant root is $\alpha$, not $\alpha^4$.
- This is independently confirmed by `CLAUDE_FINDINGS.md`: $c_n\sim K(7+4\sqrt3)^n n^{-3/2}$, and the
  three stated asymptotics there are mutually consistent only with $\lambda_1=7+4\sqrt3$, $|\lambda_2|=1$.
- $\prod_i\lambda_i=-1$ **proves** all $\sigma_p=0$ — no measurement needed, and the ledger's
  "not measured" caveat can be discharged.

$$\textbf{score}=0-5=-5,\qquad\textbf{harvestable budget}=-5+0=\boxed{-5}.$$

Not $+5.536$. This row is dead for the two-lattice architecture, and dead for a reason that
is *proved*, not measured: it is Fricke-symmetric ($\prod\lambda_i=-1$), the exact ζ(5)
analogue of Apéry's ζ(3) row, and there is no $p$-adic slope anywhere to harvest.
(The parent agent's earlier figure of $-2.366$ used $\log\lambda_1-k$ and is also wrong,
for the reason in C1; $-5$ is correct.)

### C3. ζ(7), level 24 — the slope claim contradicts the product formula

The ledger states "measured slopes at $p=2,3,5,7$: **all zero/no growth**, checked to
$n=218$ [verified]". The source (line 1233) gives

$$\bigl[(\lambda-1)(\lambda^2-8\lambda+8)(\lambda^2+4\lambda-4)\bigr]^7,$$

so per block $\prod_i\lambda_i=1\cdot8\cdot(-4)=-32$, hence

$$\sum_p\sigma_p\log p=\log32=5\log2\ \Longrightarrow\ \boxed{\sigma_2=5\ne0}.$$

Roots: $1$, $4\pm2\sqrt2=(6.8284,\,1.1716)$, $-2\pm2\sqrt2=(0.8284,\,-4.8284)$; so
$\lambda_1=4+2\sqrt2$, $|\lambda_2|=4.8284$.
$$\textbf{score}=\log\tfrac1{4.8284}-7=-8.5745,\qquad
\textbf{harvestable budget}=-8.5745+5\log2=-5.1088 .$$
The **conclusion** (hopeless, budget $\ll0$) is unchanged, but the slope measurement is
wrong and should not be relied on: a nonzero $\sigma_2$ is forced by the characteristic
polynomial. The likely cause is that the order-2 Casoratian slope test
($v_p$ of successive increments of $b_n/a_n$) was applied to an order-5/order-7 system,
where the relevant Wronskian is the full $m\times m$ one.

### C4. What survives

The census's two headline conclusions are **not** affected: Domb+T remains the only pair in
the entire census meeting the pairing precondition, and the $L(2,\chi_{-3})$ family
(Zagier B/C/F + Cooper $s_{18}$) is still structurally ruled out for lack of any decayer
($|\lambda_2|\ge1$ throughout). Both corrections make the ζ(5)/ζ(7) rows *worse*, not better.
