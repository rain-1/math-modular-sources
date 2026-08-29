# Survey: the non-Catalan contributions

*Inventory of everything of genuine value in `consolidation/` other than (a) the two closed
Catalan strands and (b) the finished $\zeta_5(3)$ theorem. Compiled 2026-08-29 by Fable
(Opus 5) from the `.md` working notes, cross-checked against `LEDGER_DUMP.md`,
`SOL_NOTES_DIGEST.md`, `ARCHIVE_CATALOGUE.md`, the `paper/` drafts and the git history.*

**Scope.** Excluded from deep reading (inventoried elsewhere): `CATALAN_*`, `FOUR_TERM_*`,
`HADAMARD_HOST.md`, `EMN_*`, `PADIC_HOLONOMY_CENSUS.md` §11 and `ADELIC_HOLONOMY.md` §2.6.
One-line notes on those appear in §7. Everything else in the directory is covered.

---

## 0. How to read this

**Status vocabulary** (the notes' own tags, which are load-bearing and are *not* silently
upgraded anywhere in the corpus — a genuinely good practice):
**[proved]** exact identity or written-out argument; **[verified: R]** exact rational or
$p$-adic computation over a stated range; **[measured]** floating-point rate, not converged
to a proof; **[conj]**; **[open]**; **[failed]**.

**Three global caveats.**

1. **Provenance split.** The material has two layers. The *older* layer (`LEDGER_DUMP.md`,
   rows C/E/H/D/F/P/B) is ChatGPT-generated; `ARCHIVE_CATALOGUE.md` §0 warns that the chat
   archive stores message text only, so every "certificate PASSED" there is self-reported and
   unverifiable as it stands. The *newer* layer (everything dated 2026-08-21 onward, written
   by Fable/Opus) has reproducible scripts under `lattice/` — 87 directories, exact PARI/GP or
   Python — and is the layer whose numbers I would stand behind. Where the two disagree, the
   newer layer wins and usually says so explicitly.
2. **Nothing here is an irrationality theorem about a classical constant.** Almost every
   note carries an explicit "No irrationality claim is made anywhere" line. The genuine
   Diophantine output is (i) the $\zeta_5(3)$ theorem (excluded here), (ii) *improved
   irrationality measures* for already-known $p$-adic values (§1.1), and (iii) a large set of
   quantified obstructions (§6).
3. **Novelty is nowhere audited.** Every "new" below is the file's own claim. Literature
   sweeps were done for Beukers 1987, Calegari 2005, Beukers 2008, CDT, AvSZ, Zagier,
   Herfurtner, Beauville; they were *not* done systematically, and `[CDT24]` in particular was
   unavailable to the project.

---

## 1. $p$-adic irrationality results

### 1.1 Improved irrationality measures for known $p$-adic $L$-values **[the best separable result in the corpus]**

**Claim.** Running the CDT holonomy recipe with the archimedean template optimised over
*all* holomorphic self-maps $\psi:(\mathbb D,0)\to(\mathbb D,0)$ — not just a chosen domain —
improves the published irrationality measures of several $p$-adic values:

| value | previous best | here | source of previous |
|---|---|---|---|
| $\zeta_3(3)$ | $22.28$ | $\boldsymbol{10.0}$ | Calegari 2005, $\theta=1.04699$ |
| $\zeta_2(5)$ | $19.7439$ | $\boldsymbol{19.39}$ | CDT ICM (their own printed value) |
| $\zeta_2(3)$ | $7.18$ | $\boldsymbol{6.3}$ | Calegari 2005, $\theta=1.16188$ |
| $L_3(2,\chi_{-3})=\zeta_3(2)$ on $X_0(3)$ | — | $\boldsymbol{6.0}$ | no previous value |
| $L_2(2,\chi_{-4})=\zeta_2(2)$ | $7.18$ | $7.2$ (*worse*) | Calegari's elementary argument is already efficient here |

with the closed form
$$\kappa=\frac{1-\gamma}{1-\sqrt{1-\text{margin}/(mL)}},\qquad \gamma=\tfrac1m,\quad
\text{margin}=m(\log\rho+L-\tau(\mathbf b))-(\mathrm{BC}+L).$$

**Status.** Computationally verified, and *calibrated exactly*: the formula reproduces CDT's
own printed $22.0724$ and $19.7439$ to $10^{-3}$ (their rounding), and Calegari's five
$\theta$ values to all printed digits from independently rebuilt rows. Margins are numerical
(converged, not interval-certified) except in the $\zeta_5(3)$ cell.

**Novelty.** The $\theta$-improvements are new *as numbers*; the mechanism (Bost–Charles vs
rearrangement numerator, outer-function templates) is a refinement of published CDT. Two
subsidiary observations are claimed new and are cheap wins: (i) the admissible search space
is larger than "pick a domain", since $H,H',\dots$ pull back meromorphically to the whole
$q$-disc so $\psi$ need not be injective; (ii) the identification of which integral each of
CDT's two printed numbers is — **circle $\to$ Bost–Charles, lune $\to$ rearrangement**.

**Dependencies.** Uses the multi-place bound of `ADELIC_HOLONOMY.md` §2.6 (excluded here,
but note: it is proved only with weight $(1-1/m)^2$ in the *rearrangement*-numerator form and
$1-1/m$ in the Bost–Charles form; the two must not be conflated). Independent of the Catalan
strand.

*Files:* `PADIC_HOLONOMY_CENSUS.md` §§1–8 (esp. §7 the measure table, §3 calibration,
§4 template optimisation). *Judgment:* **HIGH** — this is publishable on its own as a short
note ("optimised templates improve the CDT/Calegari $p$-adic measures"), and it is *entirely
independent of the $\zeta_5(3)$ novelty question*, so it survives whatever the audit says.

### 1.2 A new infinite family of 2-adic constants with $\theta_2\nearrow2$

**Claim.** The exhaustive scan of Calegari's *shifted, non-self-dual* normalisation
($1.5\times10^{11}$ triples, 2132 cells) returns 769 cells with $S_2>0$, all at $p=2$, and
every one with an identified $\xi_2$ recomputes $\zeta_2(2)$ — Calegari again. What is new is
the family
$$(n+1)^2u_{n+1}=4(2n+1)u_n+4^m n^2u_{n-1}\qquad(m\ge2),$$
integral for exactly $c=-4^m$ and no other $|c|\le30000$, with
$$S_2^{(m)}=m\log2-2,\qquad \theta_2^{(m)}=\frac{2m\log2}{2+m\log2}\ \nearrow\ 2 .$$
$m=3$ computes $-\tfrac12\zeta_2(2)$; **$m\ge4$ compute constants that are not rational, not
algebraic of degree $\le4$ (height $10^{14}$), and not affine in any $L_2(s,\chi)$ tested**
($\chi$ even quadratic, conductor $\le104$, $s\le5$, height $10^{14}$), at 2-adic precision
$2^{3800}$.

**Status.** $\sigma_2=2m$ is **[proved]** outright from the exact Casoratian
$W_n=c^n/(n+1)^2$, given only $v_2(A_n)=O(\log n)$. $k=2$ and $v_2(A_n)=O(\log n)$ are
**[measured]** to $n=400$–$600$ — the same two hypotheses the whole census runs on. So the
irrationality statement is *conditional on two measured hypotheses*.

**The file is scrupulously honest about the limitation** (§9.6): only *quadratic* characters
were tested, so "unidentified" means "not a rational-affine combination of the values
tested", **not "new"**. A weight-3 Eisenstein source with non-quadratic $\psi$ (Zagier's row
$\mathbf D$ already has quartic characters mod 5) would land outside the tested basis. The
file also notes $S_2\to\infty$ is not alarming: $\theta_2$ is the scale-free quantity and it
converges to the Roth/Ridout threshold 2, which is completely ordinary.

**Structure.** Rescaling turns the family into one operator
$L_\varepsilon=\theta^2-\varepsilon t(2\theta+1)+16t^2(\theta+1)^2$ with $\varepsilon=2^{4-m}\to0$
2-adically — a Heun equation with singular points $0,\pm\tfrac18,\infty$, the same shape as
the Zagier rows. The guess that the $\xi_2^{(m)}$ are values of a 2-adic analytic family
along $\varepsilon$ (Coleman weight-space families) *is labelled a guess*.

*Files:* `PADIC_IRRATIONALITY_CENSUS.md` §9, esp. §9.5–§9.6. *Judgment:* **MEDIUM-HIGH** —
worth a conjecture note, but the identification gap must be closed first (test non-quadratic
$\psi$). If $\xi_2^{(m)}$ turns out to be a recognisable family of $L$-values this becomes
high-value; if it stays unidentified the theorem has no content independent of the recurrence.

### 1.3 The Calegari criterion, calibrated and censused

**Claim.** Calegari's criterion in the project's normalisation is
$$S_p=\bigl(v_p(c)+\kappa_p\bigr)\log p-k-\log\lambda_1,\qquad S_p>0\Rightarrow\xi_p\notin\mathbb Q,$$
and $S_p>0\iff\theta_p>1$ in Calegari's own normalisation **[proved]**. Over the whole
census: **no cell has $S_p>0$ except Calegari's own three.** Census maximum
$S_3(\mathbf B)=-0.3521$. Best near-miss: AZ $\eta$ at $p=5$, $\xi_5=\tfrac12\zeta_5(3)$,
$\theta_5=0.8917942081$ — *numerically identical* to Calegari's own failed $X_0(5)$ cell,
because $\eta$'s level-lowering twin **is** a $\Gamma_0(5)$ row.

**The obstruction is exact and structural** (§7): for an order-2 row
$S_p\le\tfrac12v_p(c)\log p+\kappa_p\log p-k$, with equality iff the row is *Fricke-balanced*
($|\lambda_1|=|\lambda_2|$, $|c|=p^{v_p(c)}$). Calegari's rows are balanced with
$v_p(c)=\tfrac{12}{p-1}$ resp. $\tfrac{24}{p^2-1}$; the census's balanced rows have
$v_p(c)=3$ or $4$ — simply too small. **Neither lever helps:** the cusp-move orbit leaves
$v_p(c)\log p-\log\lambda_1$ invariant on every placement keeping the slope prime, and the
measured $\sigma_p$ *equals* the overconvergence-radius bound in every case including
Calegari's — there is no slack.

*Files:* `PADIC_IRRATIONALITY_CENSUS.md` §§1–8. *Judgment:* **MEDIUM-HIGH** as a
census-database + a clean no-go; the "exact structural obstruction" is the reusable part.

### 1.4 The $X_0(9)$ vs $X_0(3)$ host correction

$\xi_3^{\mathbf B}=\tfrac12L_3(2,\chi_{-3})$. On the $X_0(9)$ host the holonomy bound
**fails** (margin $-0.390$); the same value on $X_0(3)$ — $\chi_{-3}$ has conductor 3 — has
$R_3=3^6$ instead of $3^3$, budget $1.258\to7.850$, and the contradiction is easy (margin
$+4.127$). The *value* is already known irrational (Beukers 2008 Cor 7.3); what was wrong
was the census's identification of the host. A clean methodological lesson: **the right host
is the one at the conductor, not the one the row happens to be written on.**

*Files:* `PADIC_HOLONOMY_CENSUS.md` §5.4 (explicitly a correction to
`PADIC_IRRATIONALITY_CENSUS.md` §6.1). *Judgment:* **MEDIUM** — a method lemma / cautionary
remark, one paragraph in a paper.

### 1.5 $\xi_2^{\rm Zud}=\zeta_2(2)$, **proved**

**Claim.** For Zudilin's Catalan row $(Q_m,P_m)$ (arXiv:math/0201024 Thm 1),
$$\lim_{m\to\infty}P_m/Q_m=\zeta_2(2)\ \text{in}\ \mathbb Q_2,\qquad
v_2\bigl(\zeta_2(2)-P_m/Q_m\bigr)=8m-1-4s_2(m)\ \text{exactly}.$$

**Status. [proved]**, via a new exact identity relating Zudilin's row to Beukers' Padé pair:
with $x_m=\tfrac12-m$, $G_m=\sum_{j\le m}(-1)^{j-1}/(2j-1)^2$,
$$Q_m=q_m(x_m),\qquad P_m=\tfrac{(-1)^m}{8}p_m(x_m)+G_mq_m(x_m).$$
The proof is explicitly **not** the argument sketched in the 5:8 paper — that argument is
invalid as written — and the repair is an exact identity, not an estimate.

**Consequence.** The "common 2-adic Catalan period" lemma of the 5:8 paper becomes a theorem,
so the $0.9025$ construction is unconditional. (It is still $<1$, hence still not an
irrationality proof — see §6.1.)

*Files:* `ZUDILIN_2ADIC.md`; supersedes `SOURCES_S18_ZUDILIN.md` §6.3–§6.5.
*Judgment:* **HIGH** — a short, self-contained, checkable arithmetic theorem with an exact
valuation formula. Excellent Lean target. Catalan-adjacent but the statement is about
Zudilin's row and $\zeta_2(2)$, and it is not in `CONTRIBUTIONS_LEDGER.md`.

---

## 2. Holonomy method extensions (reusable method lemmas)

### 2.1 Linear independence costs exactly what irrationality costs

**Claim.** A hypothesised relation $a_0+\sum a_i\xi_i=0$ among the $r-1$ companion periods of
**one** row costs exactly **one** conditional generator $H$, whose $\theta$-orbit has size
$R$ = the holonomic rank. So $m=R+1$: *the linear-independence statement is bought at the
price of the irrationality statement.* **[proved]**

### 2.2 Fold-regularity is a codimension-$\mu$ condition

**Claim.** Fold-regularity of $H$ is codimension $\mu$ = multiplicity of the dominant
characteristic root; one rational relation supplies one condition. Hence **the whole
architecture applies only to rows with a simple $\lambda_1$.** **[proved]**, verified three
ways (AESZ 207 with $\mu=1$ folds sharply; $\mathrm{Sym}^3(\text{Zagier }\mathbf E)$ with
$\mu=3$ and $X_1(5)\,\mathrm{Sym}^2$ with $\mu=2$ do not).

**This disqualifies every $\mathrm{Sym}^w$ multi-period row in the corpus** — $\beta(4)@24$,
$\zeta(7)$ level 24, $X_1(5)\,\mathrm{Sym}^2$ — and 82 of the 227 rank-4 AESZ rows of
$z$-degree $\ge3$. It also **voids `CDT_FINDER.md` §4/§6's scoring of $X_1(5)\,\mathrm{Sym}^2$**:
that host has $\lambda_1=\varphi^5$ *double*, its companion's linear form does not decay
($|b_n-\xi a_n|^{1/n}\to11.09$, not $0.0902$), so its quoted margin $-24.4$ and its tempting
deficit $+0.34$ are both vacuous. The file calls this "the single most consequential
correction here" — and it is a **status conflict**: `CDT_FINDER.md` §6 as written is wrong.

*Files:* `HOLONOMY_LINDEP.md` §§1–2. *Judgment:* **HIGH** — a genuine structural constraint
on who can use the CDT method at all. Method-lemma paper or a section of one.

### 2.3 The holonomy irrationality measure is finite where the classical one is vacuous

**Claim.** The classical measure is $1+\log Q/\mathrm{score}$ and blows up as
$\mathrm{score}\to0^+$; the holonomy measure is $\approx2EmL/\Delta_0$ with
$\Delta_0=(m-1)(\mathrm{score}-\theta)$, $\theta=-0.8385$, and **stays finite and bounded on
the whole window $-0.8385<\mathrm{score}\le0$ where the classical argument proves nothing at
all.** **[proved]**

**Numbers.** $\mu_{\rm eff}\le\mathbf{18.16}$ for Beukers' row against the classical
$50.654$ — better by $2.79\times$; $\mu_{\rm eff}\le\mathbf{10.56}$ for $\zeta(2)$ against
Apéry's $11.851$ (better by $1.12$) but against Rhin–Viola's $5.095$ (**worse** by $2.07$);
control $\zeta(3)$: $12.11$ vs Apéry's $13.418$ vs Rhin–Viola's $5.514$. Calibrated by
reproducing CDT's own $24781$ exactly.

*Files:* `HOLONOMY_LINDEP.md` §6, §6.3. *Judgment:* **HIGH** — the structural law is the
interesting part; the Beukers-row number is a concrete improvement. Pairs naturally with §1.1.

### 2.4 A second prime is never worth anything (two proved cancellation identities)

**(a) The pure-module multi-prime identity.** For the polylogarithm module $\mathrm{Li}_j(x/s)$,
$\sum_p\varsigma_p\log p=\log|\lambda_2|$ **exactly**, so
$$\underbrace{\log\tfrac{256}{|\lambda_2|}}_{\text{archimedean}}+\underbrace{\sum_p\varsigma_p\log p}_{\text{all primes}}=\log256$$
on *every* host. The two halves enter the entry condition at weights $1$ and
$(u-1)^2/m^2=0.1837$, so $\lambda_2=\pm1$ is optimal and a second prime in $\lambda_2$ is
never worth anything. **[proved]** + **[verified]**

**(b) Two-prime decayers are inadmissible and cost entry.** Adding a well-poised decayer's
generating function contributes a *negative* $\gamma_2+\gamma_3$ (between $-0.25$ and
$-5.69$), because two-prime-ness in the lattice sense is carried by $p$-power *denominators*
($\kappa_p>0$), which the holonomy bound reads as negative overconvergence slopes.
**[proved]** + **[verified]**

**(c) The underlying distinction.** $\sigma_p>0$ is a Casoratian invariant of a **pair**
$(A,B)$, realised only in a direction with $p$-adically irrational coefficients, whereas the
holonomy bound sees only individual $\mathbb Q$-rational functions. Concretely: all fourteen of
CDT's functions have $p$-adic slope $0$ at $p=2$ and $p=3$ (verified to $y$-order 830), so
$\gamma_2=\gamma_3=0$ and the adelic margin on CDT's own host is exactly CDT's $+0.0053$,
unchanged to the last digit. The host *does* carry a 3-adic resource $\sigma_3=v_3(9)=2$,
but only in the direction $B-\xi_3A$ whose coefficients are 3-adically irrational; a scan of
all 865 rational $(a,b,c)$ of height $\le6$ finds no direction with positive slope. Were it
realisable it would be worth $\gamma_3=+0.807$ and $+11.3$ of margin — an order of magnitude
more than everything else. **It is not realisable.**

*Files:* `TWO_PRIME_HOLONOMY.md` §§1–6, `MULTI_PRIME_LATTICE.md`. *Judgment:* **HIGH** —
this is the sharpest general statement in the corpus about why the "use more primes" instinct
fails, and it is proved. Method-lemma paper.

### 2.5 Geometric denominators cost exactly $\log\lambda$

**Claim.** On non-congruence hosts, geometric denominators cost exactly $\log\lambda$
(adelic slope $-v_p(\lambda)$); the Kodaira multivalent architecture is worth $+0.839$ on the
score. Beukers' row is **re-proved by the CDT route** with margin $+0.98$.

*Files:* `CDT_NONCONGRUENCE.md`. *Judgment:* **MEDIUM-HIGH** — a clean exchange rate, and a
nice independent re-proof of a 1987 theorem by a 2024 method.

### 2.6 Proposition G and the cusp ladder

**Claim.** (i) Galois-conjugate singular points $\Rightarrow$ everything is conditional
(Proposition G). (ii) The number of conditional generators is $c\le\#\text{cusps}-1$; CDT's
own host is $c=3$; on $\Gamma_0(5)+5$ there is only one cusp so $c=1$ and no second
fold-regular class exists. **[proved]**

*Files:* `CDT_NONCONGRUENCE.md` §10. *Judgment:* **MEDIUM** — a genuinely useful cheap
screening test.

### 2.7 Adelic Padé: the right two-place formulation, and Lemma AP

**Claim.** The correct two-place problem is *not* "impose a 2-adic congruence and an
archimedean condition". It is: find $e$ with
$$\text{(a) } v_2(A_e)\ge kn,\qquad
\text{(c) } \varepsilon(e)\le24n-\psi(n),\ \psi\to\infty,\qquad
\text{(b) } |A_eG-B_e|\le e^{Fn},\ F<0,$$
where $\varepsilon(e)$ is the cancellation excess. **(a) and (c) are both $G$-free explicit
congruence conditions**, and they are compatible precisely because $k_*=22.3513<24$ — so the
window is real and the sieve is **not** circular. That is a genuinely useful reformulation
given §6.1's circularity diagnosis.

**Lemma AP [proved modulo the Lean inputs].** The nonvanishing implication, in a sharper and
$M_n$-free form. It needs only $G_2\notin\mathbb Q$ (Calegari/Beukers) and margin $\to\infty$.
**Correction to the record:** "$v_2(q_n)=o(n)$" is *not* the right sufficient condition; the
right one is $v_2(G_2-p_n/q_n)\to\infty$, and $(24-k)n$ is exactly its lower bound.

**Exact closed-form 2-adic structure of the moment family.** With $m=3n$, $M=m+j$:
$$v_2(A_{n,j})=3-2(m+M)+s_2(m)+s_2(M),\quad v_2(B_{n,j})=v_2(A_{n,j})-1,$$
$$v_2\bigl(G_2-B_{n,j}/A_{n,j}\bigr)=4(m+M)-1-2s_2(m)-2s_2(M)=24n+4j-1-\cdots,$$
verified with **zero exceptions** for all $(n,j)$, $n\le5$ (all $3n+1$ moments each).

*Files:* `ADELIC_PADE.md`. *Judgment:* **MEDIUM-HIGH** — the formulation and Lemma AP are the
reusable parts; the exact valuation formulas are clean and Lean-friendly. Catalan-flavoured
but the *shape* of the argument transfers to any two-place construction.

### 2.8 `CDT_UNPACKED.md`: the bound rederived in one page

Not a new theorem, but a clean, correct, one-page derivation of the arithmetic holonomy bound
from Siegel's lemma + Liouville + Cauchy, with the two conditions (*entry*
$\log|\varphi'(0)|>\tau$; *margin* $m>$ bound) isolated, an explanation of where the
rearrangement $\tau(\mathbf b)=\frac1{m^2}\sum_i(2i-1)\sigma_i$ comes from, and the crisp
statement that **the fold is the only place the hypothesis enters**. It also contains a
sharp observation flagged as new: *the $p$-adic avatar of the Catalan hypothesis is false*
(under $G=a/b$, $h_n\approx(b\xi_2-a)a_n$ and $\xi_2=\tfrac12\zeta_2(2)$ is irrational by
Calegari 2005, so $H$ is 2-adically as large as $A$) — Calegari's theorem is exactly what
kills the 2-adic route to the archimedean Catalan.

*Judgment:* **MEDIUM-HIGH** as an expository/method note; genuinely the best short account of
the method in the repo and the natural §2 of any paper using it.

---

## 3. Apéry-limit structure theory

This is the spine of `paper/main.tex` ("Modular Apéry systems: sources, slopes, and the
two-row design rule") and, in my judgement, the strongest coherent body of work in the
non-Catalan material.

### 3.1 Theorem B\*: the Apéry limit **is** the critical value, exactly

**Claim.** For a modular Apéry row with Eisenstein source $\Phi=P(V)E^{\psi,\varphi}_{w+2}$,
under the **endpoint condition** $\delta(\varphi)P(w+2)=0$ (equivalently: $\Phi$ vanishes at
the cusp $0$; equivalently $L(\Phi,s)$ is regular at $s=w+2$),
$$\xi_\infty=\lim_n\frac{B_n}{A_n}=L(\Phi,w+1)=P(w+1)L(\psi,w+1)L(\varphi,0),$$
so $r_\infty=P(w+1)L(\varphi,0)$, $=-\tfrac12P(w+1)$ when $\varphi=\mathbf 1$.
**There is no further rational factor.**

**Status. [proved].** The endpoint condition holds for **exactly the nine real-fold rows**
and fails for **exactly the three complex-fold rows** $\mathbf B,\delta,\eta$ — a clean
dichotomy over the twelve. Two fold geometries are identified (cusp for the (R2) rows, an
Atkin–Lehner fixed point for the (R3) rows); complex-fold constants are identified with
$\mathrm{Re}=L(\Phi,w+1)$.

**Supersession to note:** `THEOREM_B_EXACT.md` itself flags that an earlier statement
asserting a real fold for all fifteen rows is **false for the (R2) six**, and that
`ver:complex` is superseded by its own §5.

*Files:* `THEOREM_B_EXACT.md`; `paper/sections/02_sources_thmB_exact.tex`.
*Judgment:* **HIGH** — theorem paper. This is what makes the sporadic Apéry limits *explained*
rather than *observed*.

### 3.2 Theorem F: the Euler-factor criterion for $p$-adic Apéry limits

**Claim.** With $\mathcal E_p(s)=1-\psi(p)p^{-s}$ the $p$-Euler factor of the first $L$-factor:

> The row has a $p$-adic Apéry limit **exactly when** $\mathcal E_p(s)$ divides the Mellin
> polynomial $P(s)=\sum_dc_dd^{-s}$ in $\mathbb Q[p^{-s}]$. With $Q=P/\mathcal E_p$,
> $$\xi_p=-Q(w+1)\kappa_p,\qquad
> \kappa_p=\begin{cases}\tfrac12L_p(w+1,\psi\omega^{-w})&\varphi=\mathbf1\\ 0&\varphi\ne\mathbf1\end{cases}$$
> For a **cuspidal** source with $U_p$-eigenvalue $a_p$ the Euler factor is $1-a_pp^{-s}$,
> $\kappa_p=0$, hence $\xi_p=0$ whenever it exists.

**Corollary (= Conjecture D).** If $\xi_\infty=r_\infty\Lambda$ with $\Lambda=L(\psi,w+1)$
then $\xi_p=r_p\Lambda_p$ with $r_p=r_\infty/\mathcal E_p(w+1)$ — **the correction factor
depends only on $(\psi,p,w)$, not on the row**. Hence any two rows sharing $E$ have
$r_p/r'_p=r_\infty/r'_\infty$.

**Status.** Criterion and value formula **[verified exactly]** on the whole census, every row
with a positive slope at every slope prime, to $\ge p^{2991}$; the slope set
$\{p:\sigma_p>0\}=\{p:p\mid c\}$ is reproduced for all twelve Eisenstein families at
$p\in\{2,3,5\}$ with no exception. Derivation = Calegari's overconvergence argument with one
new input ($V_p$ preserves overconvergence). **Hypothesis (b) (descent to $\mathbb P^1_t$) is
recorded Open in the four rows where $Q$ still contains $V_p$**; hypothesis (b) is proved for
11/14 rows via a Mellin-shift eigenvalue lemma and (d) is weakened to a liminf.

**Reach beyond rank 2.** Theorem F transfers verbatim to rank 4: across the 297 AESZ order-4
operators with a second solution, the law
$\xi_p=\frac{r_\infty}{\mathcal E_p(m)}L_p(m,\chi\omega^{1-m})$ holds in **every** case where
both places are computable (15 operators), and $\xi_p=0$ exactly when $\chi\omega^{1-m}$ is
odd — which kills all $\pi^2$ limits and all $L(\chi_{-3},3)$ limits.

**Withdrawal to note:** the Theorem-F entry for $s_{18}$ in `EULER_CRITERION.md` §4.1 **must
be withdrawn** (per `SOURCES_S18_ZUDILIN.md` §6); the repair is the weight-drop mechanism (§3.5).

*Files:* `EULER_CRITERION.md`, `THEOREM_F_HYPOTHESES.md`, `CRYSTAL_THEOREM_F.md`,
`MUM_SURVEY.md` §4; `paper/sections/04_padic_euler_draft.tex`.
*Judgment:* **HIGH — the single best theorem in the non-Catalan material.** Clean statement,
exact criterion, verified to 3000 digits, transfers across ranks, and it *explains* Conjecture D
rather than merely checking it. Theorem paper.

### 3.3 Conjecture D proved for Zagier's rows B, C, F at $p=3$

**Claim.** $\xi_3^{\mathbf B}=\xi_3^{\mathbf C}=\tfrac12\zeta_3(2)$,
$\xi_3^{\mathbf F}=\tfrac58\zeta_3(2)$, with $\zeta_3=L_3(\cdot,\mathbf1)$ — the rational
factors $\tfrac12,\tfrac12,\tfrac58$ being *exactly* the archimedean factors of
$\Lambda=L(2,\chi_{-3})$.

**Status. [proved]** for C and F modulo (i) two standard citations (Coleman's Eisenstein
family and its overconvergence; Katz weight rigidity) and (ii) input **H1**
($v_3(a_n)=O(\log n)$, verified $n\le400$, bounded by 13). **Row B's gap is now closed too**
(update in `CONJ_D_PROOF.md` §10, from `THEOREM_F_HYPOTHESES.md` §3): $\mathbf B$ has a
second modular presentation at level 9 with $t=\eta_9^3/\eta_1^3$, $F=\eta_1^3/\eta_3$, source
$E^{\chi_{-3},\mathbf1}_3$ primitive, $\deg t=1$ so hypothesis (b) is vacuous. Numerically
$\xi_3^{\mathbf B}=\tfrac12\zeta_3(2)$ to $3^{328}$.

**Also proved:** H3 (the 3-adic period identity), previously listed OPEN — it is the statement
that the overconvergent form $\Theta_{\mathbf F}-\xi_3^{\mathbf F}$ is a $W_4$
Atkin–Lehner eigenvector on $X_0(12)$.

**Status conflict, important.** `RIGIDITY_PROOF.md` presents "Theorem R" as a proof; its own
review header **withdraws** it — Theorem R is a *reduction*, not a proof, and the
complex-analytic argument for H3 is withdrawn (`PADIC_PERIOD.md` §191: the archimedean
analogue is *not* a vanishing statement, $\Theta_{\mathbf F}|_{-1}\gamma-\Theta_{\mathbf F}\ne0$).
`ACF_ONE_SURFACE.md` then shows the degree-2 cover of `RIGIDITY_PROOF.md` §3 is an **artefact
of Zagier's level-12 normalisation** and re-proves $\xi_3^{\mathbf F}=\tfrac54\xi_3^{\mathbf C}$
from H1 plus the C-side Calegari input alone, deleting the H3 dependency entirely.
**Cite `CONJ_D_PROOF.md` + `ACF_ONE_SURFACE.md`, never `RIGIDITY_PROOF.md` §3–§4.**

*Judgment:* **HIGH** — theorem paper, jointly with §3.2.

### 3.4 A, C, F are one system in three coordinates

**Claim (Theorem 1).** Every second-order Apéry system admits two *cusp moves* — a Möbius
change of Hauptmodul plus an explicit **linear** gauge — and
$$t_{\mathbf A}=\frac{t_{\mathbf C}}{1-t_{\mathbf C}},\quad
t_{\mathbf F'}=\frac{t_{\mathbf C}}{1-9t_{\mathbf C}},\quad
F_{\mathbf A}=(1-t_{\mathbf C})F_{\mathbf C},\quad
F_{\mathbf F'}=(1-9t_{\mathbf C})F_{\mathbf C}.$$
So $\mathbf A,\mathbf C,\mathbf F$ are literally **one** Picard–Fuchs system on $X_0(6)$
written in the three Hauptmoduls that put the $I_2$, $I_3$, $I_6$ fibre at $t=\infty$.
**[proved symbolically]**

**The relating map is *not* Atkin–Lehner:** $\{1,W_2,W_3,W_6\}$ acts simply transitively on
the four cusps of $X_0(6)$, so no nontrivial AL element fixes the MUM cusp. The cusp move is
a change of Hauptmodul, not an automorphism.

**And it cannot reach A** — an exact obstruction, not a gap. The source is not
Möbius-covariant: $\Phi\mapsto\Phi/(1-\lambda t)$ moves $\Phi$ *between the two weight-3
Eisenstein families*, $\Phi_{\mathbf C}=(1-8V_2)\mathcal S$, $\Phi_{\mathbf F'}=(1+V_2)\mathcal S$
but $\Phi_{\mathbf A}=(1-V_2)\mathcal T$ with $\mathcal S=E_{3,\chi_{-3},\mathbf1}$ and
$\mathcal T=E_{3,\mathbf1,\chi_{-3}}$. So $\mathbf A$ carries $\zeta(2)$, its alignment prime
is 2 not 3, and $\xi_2^{\mathbf A}=0$ (verified to $2^{2392}$).

**Numerics.** $v_3(4\xi_3^{\mathbf F}-5\xi_3^{\mathbf C})=1586$ at $N=800$ (the row's own
Cauchy precision is 1581–1582) — twice the 786 of `RIGIDITY_PROOF.md`.

*Files:* `ACF_ONE_SURFACE.md`. *Judgment:* **HIGH** — beautiful, fully proved, and it is what
lets §3.3 drop its weakest hypothesis.

### 3.5 The weight drop: Cooper's free integration is a source

**Claim.** Cooper's rows are $w=1$ rows whose operative source is not $\Phi$ but
$$\Xi:=\theta_q^{-1}\Phi=\int_0^tF\,dt=\sum_{n\ge0}\frac{A_n}{n+1}t^{n+1}$$
(Bogner integrality), with $\xi=L(\Xi,2)$; $s_{18}$ has the source data of Zagier $\mathbf C$.
This resolves $s_{18}$'s previously open verdict and *repairs* an apparent counterexample to
the slope law.

**Self-correction to note (V1).** The originally proposed hypothesis — that the weight drop
*is* the square-root phenomenon — **is false in its specific form and true only in its
slogan**: the slogan ("one integration fewer; the operative source is weight three, not the
meromorphic weight-four $\Phi$") is now a theorem, but the identification of that weight-three
object with the square-root row's source $\Psi_{\rm root}=g^3u$ is **wrong**. The file says so
explicitly and repairs `SOURCES_S18_ZUDILIN.md` §0.2/§3.4 and `EULER_CRITERION.md` §4.1 note 8.

*Files:* `WEIGHT_DROP.md`; `paper/sections/04_padic_weightdrop_draft.tex`. Explicitly
supersedes `SOURCES_S18_ZUDILIN.md` §3.4.
*Judgment:* **MEDIUM-HIGH** — a genuine structural correction; a section, not a paper.

### 3.6 Roots of $\mathrm{Sym}^w$ rows: Theorems R1–R4

* **R1 (integrality, sharp).** For every $G\in1+x\mathbb Z[[x]]$ and $w\ge2$,
  $\lambda_w^n[x^n]G^{1/w}\in\mathbb Z$ with $\lambda_w=w\cdot\mathrm{rad}(w)=\prod_{p\mid w}p^{v_p(w)+1}$,
  and $\lambda_w$ is **minimal** ($\lambda$ works for all $G$ iff $\lambda_w\mid\lambda$).
  Graded form $\lambda(G)=\prod_{p\mid w}p^{\max(0,v_p(w)+1-e_p)}$ with
  $e_p=\min_nv_p([x^n]G)$.
* **R2 (descent).** The $w$-th root of a $\mathrm{Sym}^w$ solution solves the second-order
  equation, for all $w$.
* **R3 (free integration).** The companion costs $d_n^2$ whatever $w$ is — i.e. **every**
  third-order sporadic row has an integral square root with one integration fewer ($k=2$
  instead of $k=3$).
* **R4 (score bookkeeping).**

**Census.** Apéry's square root is the only positive-score root ($+0.1392$); $\beta(4)@24$ is
$\mathrm{Sym}^3$ of Zagier $\mathbf E$; obstructions at $w\ge4$ (interior zero of $F$;
fractional exponents).

**Attribution, correctly handled.** The irrationality of the square-root period is
**Beukers 1987, Theorem 3** — *not* new; the project initially thought it was and corrected
itself (commit `d31aef0`). What *is* the project's own: the explicit recurrence
$(136n^2+68n+10\mid4(2n-1)^2)$, the integrality theorems R1–R4, the period identification
$$\xi=L(\Psi,2),\qquad \Psi=\tfrac14(\eta_1\eta_6)^{9/2}(\eta_2\eta_3)^{-3/2}$$
(a **non-congruence** weight-3 eta quotient), the measure $\mu\le50.66$ (improved to
$\le18.16$ by §2.3), and a Lean formalization blueprint whose single external input is nome
integrality.

*Files:* `ROOT_ROWS.md`, `SQRT_APERY.md`, `SQRT_APERY_FORMAL.md`, `NOME_INTEGRALITY.md`;
`paper/sections/02_root_rows_draft.tex`, `02_sqrt_rows.tex`.
*Judgment:* **HIGH** for R1–R3 + the period identification; **MEDIUM** for the Lean blueprint
(good target, single external input: $\sqrt{1-34t+t^2}\in\mathbb Z[[t]]$, which is elementary,
and Theorem 4.2 $\iff$ Gauss congruences for $1/(F\sigma)$ via Dwork's lemma, with
Krattenthaler–Rivoal as the only external $p$-adic input).

### 3.7 The cusp move as a group action: four theorems and a rigidity obstruction

* **Theorem 1 (general move).** For **every** second-order row (not just Zagier's
  normalisation), $s=t/(1-\lambda t)$, $z=(1-\lambda t)^\alpha y$ gives again a three-term
  row of degree $\le2$ **iff** $(\alpha-1+r_1)(\alpha-1+r_2)=0$, i.e.
  $\alpha\in\{1-r_1,1-r_2\}$, the two local exponents at $t=\infty$; and
  $L^\sharp z=(1-\lambda t)^{1-\alpha}Ly$. **[proved symbolically]**
* **Theorem 2 (the companion).** The naive guess $B^\sharp=(1-\lambda t)^\alpha B$ is
  **false**; correctly $B^\sharp=(1-\lambda t)^\alpha\widehat B$ with
  $L\widehat B=t(1-\lambda t)^{-\alpha}$ — *the inhomogeneity is replaced by the new
  Hauptmodul*, i.e. $B^\sharp=F^\sharp D^{-2}(\Phi/(1-\lambda t))$. **[proved]**, verified in
  all 22 available cases to $n=40$.
* **Theorem 3 (orbits).** Exactly three placements, up to four gauges: orbits of size $\le12$.
  Measured sizes $2,3,8,12$; all 55 listed members (49 distinct rows) integral to $n=200$
  after rescaling $c\in\{1,2,4,9\}$, all with $k\in\{1,2\}$, sharp.
* **Theorem 4 (rigidity — the fatal obstruction). [proved]** For an integral row with
  $a=\mathrm{lc}\,P$, $d=\mathrm{lc}\,Q$, $a^2\ne4d\ne0$:
  $$a^2-4d\ \text{a square}\ \Longrightarrow\ |\lambda_2|\ge1\ \text{at every placement}
  \ \Longrightarrow\ \mathrm{score}\le-k\le-1 .$$
  **No rational cusp-move orbit can contain a positive-score row.** Conversely a
  positive-score row necessarily has irrational characteristic roots, so its only equally-good
  placement is its **Galois conjugate** over $\mathbb Q(\sqrt{a^2-4d})$.

**The programme's headline hypothesis is half true.** *False archimedeanly*: across the 8
rational orbits the only rationally proportional pair is $(\mathbf C,\mathbf F')$ with ratio
$-\tfrac54$ (80 digits); every other pair has no rational relation of height $\le10^{12}$.
Strikingly, Beukers' Galois-conjugate placement approximates a **different** real number
($0.1002150533\ldots$ vs $0.1001874492\ldots$). *True $p$-adically, non-trivially*: inside the
$\sqrt{\text{Domb}}$ orbit, placements 2 and 3 satisfy $\xi_p^{(2)}=-2\xi_p^{(3)}$
**simultaneously at $p=2$ and $p=3$** (to 792 and 386 digits) while their archimedean periods
are $\mathbb Q$-unrelated.

**Two genuinely new identified rows** in the $\sqrt{\text{Domb}}$ orbit:
$$(20,10,64):\ \xi=\tfrac1{16}(15L(2,\chi_{-3})-6\zeta(2)),\qquad
(20,4,64):\ \xi=\tfrac1{32}(6\zeta(2)+15L(2,\chi_{-3})),$$
both integral to $n=200$, $k=2$ sharp, $\lambda_1=16$, $\lambda_2=4$. These are the **first
rows in the project whose single Apéry limit is a non-trivial $\mathbb Q$-combination of
$\zeta(2)$ and $L(2,\chi_{-3})$** — exactly the shape a Nesterenko argument for
$\dim\langle1,\zeta(2),L(2,\chi_{-3})\rangle$ wants. Their parent carries the **cuspidal**
value $L(f_{12},2)$, so the cusp move pushed a weight-3 source out of the cuspidal part into
the Eisenstein part.

*Files:* `CUSP_MOVE_PROGRAM.md`. *Judgment:* **HIGH** — Theorems 1–4 are proved, the census is
reusable, and Theorem 4 is a real structural obstruction that saves future effort.

### 3.8 Good-prime tower trichotomy

**Claim.** A 144-cell census gives a trichotomy for the tower Frobenius: the extension
**splits** at good primes and is **unipotent** at slope primes.

**Serious caveat, self-flagged.** `GOOD_PRIME_TOWERS.md` says the originally proposed
*mechanism* "is wrong in every particular — there is no unit root, no [...]" and that the
claimed structure (a specific fibre statement) "is false". So the trichotomy survives as a
*measured classification*; the explanation attached to it did not. `AESZ207.md` §5 then shows
the trichotomy **cannot be used** for AESZ 207: the ratios $t_{s+1}/t_s$ do not converge and
their valuations are not even constant (at $p=3$ they alternate $-3,-4$, mean $-\tfrac72$,
non-integral), so $\chi(3),\chi(5),\chi(7)$ cannot be read off.

*Judgment:* **LOW-MEDIUM** — keep as a census observation with the mechanism removed; do not
build on it.

### 3.9 Theorem F in Dwork-crystal language — the conceptually deepest item here

**Claim.** *Theorem F is not about Eisenstein series.* The existence criterion is **not** the
Euler factor but a **Frobenius-slope separation** $\sigma_p>0$ on an extension of a MUM
crystal by $\mathbb Q_p(0)$, with $\sigma_p=v_p(c)$ for every (R2)/(R3) row **[proved]** (exact
Casoratian) **+ [verified]** on 84 cells. The Eisenstein/Euler-factor formulation is the
modular *computation* of that separation — and it computes the **value**, which crystal theory
does not.

The trichotomy is likewise a statement about the *extension* crystal, with eigenvalue
$\varepsilon(p)p^{-(w+1)}$ where $w+1$ is the **Hodge depth (denominator exponent), not the
modular weight** — verified on the 144 old cells plus **20 new non-modular cells**.

**New verified consequences** (all *new*, per the file):
* **Zudilin's Catalan row has a Dwork crystal.** $R_m:=16^mQ_m\in\mathbb Z$ satisfies the full
  Dwork tower $R_{mp^s}\equiv R_{mp^{s-1}}\pmod{p^s}$ — 2288 tests at $p=2,3,5,7$, **0
  failures**. *This contradicts the reading of `SOURCES_S18_ZUDILIN.md` §5 that the row has no
  geometric model* (a status conflict; the crystal file wins).
* Its good-prime tower eigenvalue is $\chi_{-4}(p)p^{-2}$ **exactly** ($p=3,5,7,11$) — the
  character of $G=L(2,\chi_{-4})$, **read off a non-modular row**. That is a striking way to
  recover the character.
* Brown–Zudilin's $\zeta(5)$ row is also Dwork (1716 tests, 0 failures), eigenvalue $p^{-5}$
  exactly at $p=2,3,5,7$ (and $p^{-3}$ on the $\zeta(3)$ direction), trivial character —
  refining the "no tower limit" of `ZETA5_TWO_ROW.md` §4.
* **The cuspidal row is *not* a Dwork row** (531 failures, first at $p=2,s=2,m=1$) yet still
  has $\xi_2=0$: **the two hypotheses are independent.**
* The crystal argument reproduces the exact Zudilin tail $v_2=8m-1-4s_2(m)$ **[proved]** from
  the Casoratian + the denominator law (re-verified to $m=200$) — but the **value**
  $\xi_2^{\rm Zud}=\zeta_2(2)$ stays **[open]** in crystal terms, because the value is the
  syntomic regulator. §1.5's Padé proof remains the only proof of the value.

**Supporting exposition.** `DWORK_CRYSTALS_PRIMER.md` is a genuinely good primer (Cartier
operator, unit root, Dwork's lemma, the Beukers–Vlasenko Frobenius matrix, and where this
project sits inside it), written from six papers read in full with theorem numbers, plus
citation corrections for `CRYSTAL_THEOREM_F.md` and a verification that the Apéry polytope is
reflexive. Not a research contribution, but the best on-ramp in the repo.

*Files:* `CRYSTAL_THEOREM_F.md`, `DWORK_CRYSTALS_PRIMER.md`, `THEORY_NOTES_05_syntomic.md`,
`THEORY_NOTES_06_deeper_object.md`. *Judgment:* **HIGH** — this is the right conceptual
framing of §3.2 and it makes the theorem apply to non-modular rows. It should be §2 of the
Theorem F paper, or a companion note. Note honestly that the *value* half is modular-only.

### 3.10 The Casoratian sign theorem and the Hankel-moment criterion

* **Positivity is a Casoratian phenomenon, free for the whole corpus. [proved]** For any row
  in Zagier normalisation the exact Casoratian $a_nb_{n-1}-a_{n-1}b_n=-c^{n-1}/n^{w+1}$
  telescopes to
  $$a_N\xi-b_N=a_N\sum_{n>N}\frac{c^{n-1}}{n^{w+1}a_na_{n-1}},$$
  whose summands have constant sign $\mathrm{sign}(c)^{n-1}$. So $c>0,a_n>0\Rightarrow
  a_N\xi-b_N>0$ for **every** $N$, unconditionally, **with no integral anywhere**. Verified
  in exact arithmetic for all twelve sporadic rows, $n\le110$. $\mathrm{sign}(c)$ classifies
  all 12.
* **The right criterion for "carries a Beukers-type positive integral" is Hankel
  positivity.** A Beukers/Zudilin/Nesterenko kernel is $r_n=\int u^n\,d\mu(u)$, $\mu\ge0$; so
  the row admits one iff $(a_n\xi-b_n)$ is a Hamburger moment sequence — decidable. All nine
  corpus rows with a real archimedean limit pass ($\det H_d>0$, $d\le10$), with the measure on
  $[0,\infty)$ exactly when $c>0$ and on $(-\infty,0]$ exactly when $c<0$. **[verified]**

**A useful negative:** "products of Zudilin and Nesterenko kernels" is *empty* — the product
of two half-integer Catalan kernels has integer exponents throughout, so it is a $\zeta(2)$-world
Beukers integral and produces a linear form in $1,\zeta(2)$, not in $1,G$.

*Files:* `POSITIVITY_PROGRAM.md` §§2, 6. *Judgment:* **MEDIUM-HIGH** — the Casoratian theorem
is a one-line proof of something people usually establish by exhibiting integrals; the Hankel
criterion turns "does a Beukers integral exist?" into a decidable question. Method lemma.

### 3.11 Theorem E: the two-row master formula and the design rule

**Claim.** The construction quality of a two-row lattice construction has a closed form
$\delta=(H-F)/H$ in $(\lambda_1,\lambda_2,k,\kappa_p,\sigma_p)$ and the sampling ratio, which
**reproduces every number the project has produced**: $0.857914$, $0.9025$, $0.9010$,
$0.6589$, and the Lean $1-\varepsilon$ at $k^*$. The **design rule** is
$$\delta>1\iff\log\Lambda_{\rm dec}>k\max(r,1)+r\log\rho_2^{\rm eng}$$
(corrected for $\kappa$; a further correction $w+\kappa_p>0$ comes from
`ONE_CLASS_TWO_WORLDS.md`). **[proved modulo the papers' selection theorem]** — and note
§6.1: that selection theorem is exactly where the method fails.

**The $\zeta(3)$ calibration.** Running the construction verbatim on the aligned pair
Domb $(10,4,64)$ + $\mathbf T$ $(12,4,16)$ gives
$$\delta=\frac{12\log(12+8\sqrt2)-24\log2}{9+9\log(12+8\sqrt2)-20\log2}=0.9009531686558563\ldots<1,$$
so it is **not** a second proof of Apéry's theorem — it is a clean calibration showing the
method lands $\approx10\%$ short on $\zeta(3)$ exactly as it does on Catalan ($0.90253$). The
file explicitly says the near-coincidence is an accident of the numbers and is useful only as a
transcription check. Good discipline.

**The alignment-prime census** (`MULTI_PRIME_LATTICE.md` §1, **[verified]**). *Every* decayer
previously known to the project has **exactly one** alignment prime: Zudilin's row $\{2\}$
($\kappa_2=4$, $\sigma_2=8$); Nesterenko's $(4,7)$ row $\{2\}$ ($\kappa_2=14$, $\sigma_2=28$);
the $\chi_{-3}$ conductor-3 hypergeometric family $\{3\}$ only; $\mathbf T$, Domb,
Brown–Zudilin, $\eta$, AESZ 184 integral with one prime or none. The unique row among the
twelve sporadic families with **two** slope primes is $\mathbf F$ $(17,6,72)$, $c=72=2^3\cdot3^2$
— **and it is integral, hence not a decayer**. A new conductor-6 two-prime decayer was built
($\sigma_2=8$, $\sigma_3=3$, both $p$-adic limits exactly as Theorem F predicts, 352/951
digits) and the two-prime lattice reads $\delta:0.51\to0.65\to1.05$ — **but $F<0$ there and a
rational surrogate reproduces it verbatim**, so it carries no information; the honest $F>0$
best is $0.94$.

*Files:* `ZETA3_TWO_LATTICE.md`, `MULTI_PRIME_LATTICE.md`, `ONE_CLASS_TWO_WORLDS.md`,
`paper/sections/05_two_row.tex`. *Judgment:* **MEDIUM-HIGH** — the master formula is a real
unification and the alignment census is strong evidence for §5's conjecture; but every
$\delta>1$ reading in the corpus is neutralised by the rational-surrogate control, and the
files say so up front, repeatedly and correctly.

---

## 4. Period / motive censuses and classifications

### 4.1 Classification of second-order integral Apéry-like rows **[the flagship census]**

Already written up as `paper/classification/main.tex` (21 pp, `classification_paper.pdf`).

**Theorem A (conditional classification, proved).** If the Picard–Fuchs operator of an
integral row $(n+1)^2u_{n+1}=P(n)u_n-Q(n)u_{n-1}$, $\deg P,\deg Q\le2$, is that of an
elliptic surface over $\mathbb P^1_t$, then for the Zagier shape $(an^2+an+b,\,cn^2)$ all four
singular fibres are semistable, $\deg\mathcal J=12$ **exactly** — by Riemann–Hurwitz together
with a **twist-parity constraint $12\mid\deg\mathcal J$** — and the configuration is, up to
quadratic twist, one of Beauville's six. No Kodaira-dimension hypothesis is needed. The
complete list of cusp placements is Zagier's seven rows $\mathbf A$–$\mathbf G$, with
$\mathbf A,\mathbf C,\mathbf F$ the three placements of $I_6I_3I_2I_1$ and $\mathbf E,\mathbf G$
the two of $I_1I_1I_2I_8$.

**Theorem B (unconditional necessary conditions).** Kodaira admissibility of exponent
differences; **integrality rigidity** (rational characteristic roots are integers, so the
irrationality score is $\le-k$ on every rational cusp-move orbit — this is §3.7 Theorem 4);
Riemann–Hurwitz. These exclude the square roots of Cooper's $s_{10},s_{18}$ and an infinite
$\arctan$-Padé family.

**The honest framing.** "The unconditional classification is exactly Zagier's conjecture that
any integral solution has a modular parametrisation, i.e. Bombieri–Dwork in rank two with four
singular points; we make the reduction precise and give the evidence." That is the right
level of claim, and the paper makes it.

**Supporting structure theorems** (`HERFURTNER_CLASSIFICATION.md`):
* **H1.** $\rho_1=\rho_2=\rho$ iff $b=(1-\rho)a$, $e=-2\rho d$. Hence **"Zagier's
  normalisation is not a choice — it is the statement that all four fibres are semistable"**,
  and the root-row shape $(An^2+\tfrac A2n+B,\,C(2n-1)^2)$ is the statement that
  $t_1,t_2$ carry $III/III^*$ and $\infty$ an $I_m^*$. A genuinely illuminating result.
* **H2 (Kodaira admissibility).** Every local exponent difference lies in
  $\{0\}\cup\{\tfrac12\}\cup\{\tfrac13,\tfrac23\}\bmod1$. Excludes three census objects
  outright — precisely three of the rows whose periods the identification batteries had failed
  to recognise.
* **H3 (cross-ratio invariant, finiteness).** $\mathcal I=(1+z)^2/z=a^2/d$ with
  $z=t_1/t_2=\lambda_2/\lambda_1$. Since $a,d\in\mathbb Z$, $\mathcal I\in\mathbb Q$, and every
  Herfurtner configuration with irrational or non-real $\mathcal I$ is **excluded by
  integrality alone**. Consequence: the rational values available in the all-cusp class are
  $\mathcal I\in\{0,3,\tfrac92,-\tfrac{49}8,\tfrac{100}9,\tfrac{289}{72},-121\}$ and Zagier's
  $\mathbf B,\mathbf E,\mathbf A,\mathbf C,\mathbf F,\mathbf D$ realise the last six.
  **So "the six" of Beauville–Zagier is, on the recurrence side, the count of rational
  cross-ratios in Herfurtner's four-cusp block.** This is a lovely reframing.
* **H4.** Of the nine $\mathrm{Sym}^1$ square roots, exactly **two** are Picard–Fuchs systems
  of a rational elliptic surface over the row's own $\mathbb P^1_t$ ($\sqrt\eta$ on $I_1I_5\,III\,III$,
  $\sqrt{\mathrm{AZ}(9,3,-27)}$ on $I_3I_3\,III\,III$). Apéry's own square root, $T$, Domb's,
  $\mathrm{AZ}(7,3,81)$'s and $s_7$'s match **no** Herfurtner configuration: their projective
  monodromy is not conjugate into $\mathrm{PSL}_2(\mathbb Z)$, they live on Atkin–Lehner
  quotients $X_0(N)/W$. **This is the precise sense in which Beukers' row is
  "non-congruence", and it is not the CM/non-CM dichotomy.**
* **H5 (the scan).** 26 Kodaira-admissible classes, $|a|\le3000$, $|c|\le150$, $|d|\le12000$,
  $1.27\times10^6$ integrality hits, tested prime-by-prime to $n=30$ and re-verified exactly
  to $n=300$. Recovers every census row and produces **two genuinely new Herfurtner rows**:
  $(117n^2+78n+21\mid441(3n-1)^2)$ on $I_1I_7\,II\,II$ (monodromy $\Gamma_0(7)$), and
  $(72n^2+36n+6\mid108(4n-1)(4n-3))$ on $I_3\,III\,III\,III$. Both $k=2$, complex-conjugate
  roots, no archimedean limit. The second has
  $$\mathrm{Im}\,\xi=-\frac{\Gamma(1/4)^4}{2^{7/2}3^{11/4}\pi}=-\frac{\varpi^2}{\sqrt2\,3^{11/4}}
  \qquad\textbf{[verified to }10^{-249}\textbf{]},$$
  the **lemniscatic Chowla–Selberg period** of discriminant $-4$ — and it is *the first
  period in the project's census that is a Chowla–Selberg period rather than a critical
  $L$-value*. Everything else about it resists a 340-target sweep.
* **H6 (the only positive scores).** Over the whole scan the only Casoratian-non-degenerate
  rows with $\log(1/|\lambda_2|)>k\ge2$ are **Apéry's $\zeta(2)$ row** $(11,3,-1)$ at
  $+0.4061$ and **Beukers' square root of Apéry's row** $(136,10,4)$ at $+0.1392$. Every
  other positive score has $k=1$ and is a classical Legendre–Padé row for a logarithm. The
  test is a condition on $(a,d)$ alone, so **this pass is complete over all 1 270 065 hits.**
* **Proposition K1 (`HERFURTNER_K3_WINDOW.md`) — the $K3$ worry, closed. [proved]** For *any*
  elliptic surface (rational, $K3$, or higher $\chi$) whose Picard–Fuchs local system is that
  of a second-order row with four singular points, Riemann–Hurwitz forces
  $$\deg\mathcal J\ \le\ 6c+4e_3+3e_2-12\ \le\ 12,$$
  with $c,e_2,e_3$ counting the singular points whose projective local monodromy is unipotent,
  of order 2, of order 3; equality only in the all-cusp class. **So the $\deg\le12$ window of
  H4 was complete all along** — a $K3$ row is not excluded geometrically, but any such row
  would already have been found. Re-running at $\deg\le24$ agrees. This closes open item 4 of
  `HERFURTNER_CLASSIFICATION.md` §8 and is exactly the bound quoted in the K3 note (§4.2) to
  show three-term recurrences can never give a cusp-form Apéry constant.
* **Structural dichotomy.** Rigid Herfurtner configurations support finitely many integral
  rows; one-parameter families (exactly those containing $I_m^*$ or $I_0^*$) support infinitely
  many, and those are the classical Legendre–Padé rows. **Zagier's finiteness statement is a
  statement about the rigid block.**

*Judgment:* **HIGH — already paper-shaped and, with §3.2, the thing I would send out first
after $\zeta_5(3)$.** The attribution is handled correctly (Zagier 2009 §7 + Beauville for the
correspondence; new = the mechanism, the twist-parity lemma, the cross-ratio invariant, the
cusp-move enumeration, and completeness).

### 4.2 The K3 row: the first cusp-form Apéry constant

**Claim.** The four-term row
$(n+1)^2u_{n+1}=(11n^2+11n+4)u_n-(37n^2+3)u_{n-1}+3(3n-1)(3n-2)u_{n-2}$,
$u_n=1,4,16,64,250,928,3136,8704,11866,-79400,\dots$, is the period sequence of an elliptic
$K3$ with fibres $I_4,I_6,I_3,I_3,IV^*$ and $\rho=20$. Its transcendental motive is that of
the weight-3 CM newform **32.3.d.a** (CM by $\mathbb Q(\sqrt{-2})$), and its Apéry constants
are critical $L$-values: fold constants at $\tfrac{5\mp i\sqrt2}{27}$ are
$\tfrac{\sqrt2\mp i}3L(g,2)$ and the constant at $t=1$ is $\tfrac{2\sqrt2}3L(g,2)$, with
$L(g,2)=\tfrac\pi{32}\Gamma(\tfrac18)\Gamma(\tfrac38)/\Gamma(\tfrac58)\Gamma(\tfrac78)$ by
Chowla–Selberg (230 digits).

**"To our knowledge this is the first Apéry-like sequence whose Apéry constant is a critical
value of a cusp form; three-term integral recurrences can never produce one, by a
Riemann–Hurwitz/twist-parity bound that confines them to rational elliptic surfaces."**

**Status.** Written up (`paper/k3_row/main.tex`, `k3_row_note.pdf`, 4 pp). Full data: explicit
Weierstrass model; $\rho=20$ (singular $K3$), $\mathrm{MW}$ rank 0, $\mathrm{MW}_{\rm tors}=\mathbb Z/3$
(section over $\mathbb Q(\sqrt3)$); transcendental lattice $T\cong\mathrm{diag}(6,12)$,
$\det T=72$; newform $\texttt{32.3.d.a}=\texttt{8.3.d.a}\otimes\chi_{-4}$, weight 3, level 32,
nebentypus $\chi_{-8}$; $a_p$ match for 76 primes; $L(g,2)$ to 230 digits.
**Open: integrality of the row.** Note $\xi_2=0$ with a log rate, so **Proposition C is
Eisenstein-only**.

**Overlap:** this is also covered by the four-term paper (excluded from this survey); the
standalone note is the additional artefact.

*Judgment:* **HIGH** — a clean, striking, already-drafted result with a proved impossibility
statement attached to it. Send it.

### 4.3 The rank-4 MUM / AESZ survey

* **The quintic has no Apéry limit at all** (its recurrence has order 1) — this is *forced*,
  and its $\zeta(3)$ lives one level up, in the connection coefficients / $\widehat\Gamma$-class:
  $V=\frac{\sqrt5}{4\pi^2}(\hat y_3-10\zeta(2)\hat y_1-40\zeta(3)\hat y_0)$ with
  $-40=\chi/H^3$, verified to 116 digits. So "the quintic is the wrong test case".
* **The $p$-adic $\widehat\Gamma$-class.** $[\rho^3]\log\frac{\Gamma_p(1+5\rho)}{\Gamma_p(1+\rho)^5}=-40\zeta_p(3)$
  for every $p$ including the slope prime $p=5$, via the identity
  $(\log\Gamma_p(1+x))^{(k)}|_0=-(k-1)!\zeta_p(k)$ for odd $k$, $0$ for even $k$ — the exact
  mirror of $\log\Gamma(1+x)$ with the even part deleted, verified at seven primes and three
  weights. **So the rational number $\chi/H^3$ is the same at every place.** Nice.
  **NOTATION WARNING — see §7.**
* **Theorem F transfers verbatim to rank 4** (see §3.2).
* **The $p$-adic side is strictly stronger than the archimedean one.** 25 operators whose
  archimedean limit the numerics cannot reach — **five of them because it does not exist** —
  have $p$-adic limits identified to 185–2441 digits. Five resulting archimedean predictions
  were confirmed against AvSZ's independent PSLQ table, and **two new period classes appear
  that are not in AvSZ's list at all: $L(\chi_8,2)$ and $L(\chi_5,3)$.**

*Files:* `MUM_SURVEY.md`. *Judgment:* **HIGH** — "compute the $p$-adic limit when the
archimedean one is unreachable or nonexistent" is a genuinely useful new tool, and the two new
period classes are a concrete deliverable.

### 4.4 AESZ 207: an erratum against a published table

**Claim.** Almkvist–van Straten–Zudilin (*Apéry limits…*, Fields Inst. Commun. **54** (2008))
§4.8 print, for #207,
$-0.00050462505145900474057831709244307528529622730007723$. **This is wrong.** The true limit
is
$$\xi_\infty(207)=-0.00050455459344136708862542545797071516117215446057914345028\ldots$$
**[verified to 1139 digits]** ($N=5000$). They disagree from the 4th significant figure
($7\times10^{-8}$). The reason is visible: AvSZ write "these limits converge very fast", but
for #207 $|\lambda_2/\lambda_1|=0.5947411$, i.e. **0.2257 digits per index — the slowest in
their table**. Their fifty digits are un-converged. *This alone explains why no PSLQ
identification was ever found: the target number was wrong.*

**The pipeline is validated**: it reproduces AvSZ's own identified entries (AESZ 16, 28, 29,
42, 58, 182, 185) to $10^{-347}$, and their #17, #34 to 50 and 31 digits.

**Further.** Geometry located: #207 is the reflection at $z=\infty$ of AESZ 99 (degree-13
$5\times5$-Pfaffian in $\mathbb P^5$); conifolds over $\mathbb Q(\sqrt{17})$; Riemann scheme
re-derived independently. **The two places see different singularities**: the archimedean
limit is governed by the $\mathbb Q(\sqrt{17})$ conifold while $\sigma_2=12=v_2(53248)$
comes from the *rational* singularity. $\xi_2$ known to **5958** 2-adic digits (the survey's
2449 was an artefact). Exclusions are sharp: neither constant is a rational multiple or
2-term rational combination of any of $\sim350$ archimedean / 188 2-adic Dirichlet or
Kubota–Leopoldt targets, nor in the $\mathbb Q$- or $\mathbb Q(\sqrt{17})$-span of the
$\widehat\Gamma$-class basis at the dominant conifold, nor algebraic of degree $\le8$ (arch.)
or $\le6$ (2-adic); and $\xi_2\notin\mathbb Q(\sqrt{17})$.

*Judgment:* **HIGH** — an erratum against a published table, with a validated pipeline and a
sharp exclusion set, is exactly the kind of thing that is cheap to write, easy to check, and
useful to the AESZ community. Short note. **This is the most immediately shareable item in
the whole non-Catalan corpus.**

### 4.5 Sporadic scans

* **Scan 1** (`SPORADIC_SEARCH.md`): all 14 genus-zero levels, **520 005** pairs $(t,F)$ with
  both eta quotients. Reproduces all six AZ third-order rows and five of six Zagier rows with
  their limits (the acceptance test). **No new sporadic row in Zagier normalisation**; every
  extra hit is a two-term degeneration.
* **Scan 2** (`SPORADIC_SCAN2.md`): $5.5\times10^5$ exact bilinear solves, 3979 parameters,
  51 levels $N\le104$, weights 1–3, the *whole* integral form space (not an eta-quotient
  ansatz), fixing scan 1's two named gaps. **No new $\mathrm{Sym}^2$ or $\mathrm{Sym}^3$ row.**
  What it produces: a family of weight-one rows over Atkin–Lehner quotients, one genuinely new
  and good — over the $\zeta(3)$ parameter $T$, the weight-one CM theta series of
  $\mathbb Z[\sqrt{-2}]$ gives $a_n=1,2,18,236,3634,\dots$ with $\lambda_1=12+8\sqrt2$,
  $\lambda_2=12-8\sqrt2<1$, $k=2$, budget $+1.149$, period exactly $L(f,2)$ for the weight-3
  level-8 CM newform $f=\eta_1^2\eta_2\eta_4\eta_8^2$.
* **Part II (the main result):** *every* sporadic order-three row (six AZ + three Cooper) has
  an integral $\mathrm{Sym}^1$ square root with a three-term quadratic-coefficient recurrence
  and $k=2$ instead of $k=3$ — a free integration. (= §3.6 R3.)
* **A clean dichotomy.** Three of the nine $\mathrm{Sym}^1$ rows have a weight-one form that
  is an integral CM theta series of class number one, with period $L(f,2)$ for the
  corresponding weight-3 CM newform, **exactly and with no rational factor**: disc $-7$
  ($s_7$) $\to L(\eta_1^3\eta_7^3,2)$; disc $-8$ ($T$) $\to L(\eta_1^2\eta_2\eta_4\eta_8^2,2)$;
  disc $-3$ (Domb) $\to L(\eta_2^3\eta_6^3,2)$. The other four resist a 1792-element `lindep`.

*Judgment:* **MEDIUM-HIGH** as a census-database + negative result; the CM-theta dichotomy is
the publishable nugget. Note `SPORADIC_SCAN2.md` §252 flags one of its own normalisation
artefacts.

### 4.6 Non-congruence scan: Theorems N1–N3

* **N2 (finiteness). [proved]** (Gauss–Bonnet) **+ computed exhaustively.** Four special
  points $\Rightarrow$ covolume index $\le12$; there are exactly **28 four-special-point
  genus-zero subgroups of $\mathrm{PSL}_2(\mathbb Z)$, of which 12 are non-congruence**, all
  enumerated. **Check:** the index-12 four-*cusp* row returns exactly six groups, widths
  $[1,1,5,5],[1,2,3,6],[1,1,2,8],[1,1,1,9],[2,2,4,4],[3,3,3,3]$, **all congruence** —
  Beauville's six families / Zagier's six rows, **recovered from scratch**. Hence
  ***no* non-congruence group has four cusps: every non-congruence host has an elliptic point.**
* **N3 (unequal exponents kill the score). [proved]** If the exponent differences at
  $t_1,t_2$ differ then $t_1,t_2\in\mathbb Q$, hence (for the monic integral rows of the
  census) $\lambda_{1,2}$ are **rational integers**; a positive score needs $|\lambda_2|<e^{-k}\le1$
  and the only integer of absolute value $<1$ is $0$, which degenerates the recurrence to
  first order. A genuinely elegant argument. **Consequence:** all four index-7 non-congruence
  groups are excluded outright, and index 9/10 reduce to a class already scanned — so
  **a non-congruence host *group* contributes nothing beyond the existing scan; the
  interesting non-congruence-ness is not in the group.**
* **Beukers' Theorem 3 row is the unique positive-score non-congruence row.** A second
  Beukers-type row on $\Gamma_0(5)+5$ is exhibited; the $\arctan$-Padé family identifies the
  previously unrecognised §4.4 candidate.

*Judgment:* **MEDIUM-HIGH** — a small, complete, checkable classification with two proved
theorems; N3 is a nice self-contained lemma.

### 4.7 The row ledger

A census database of Apéry-like rows scored for irrationality potential, with explicit
definitions (score $=\log(1/|\lambda_2|)-k$; budget; headroom $=\sum_p\sigma_p\log p$), every
number tagged, and "unidentified" left unfilled rather than guessed. **It carries its own
audit section correcting itself** (see §7). Genuinely reusable.

*Judgment:* **MEDIUM-HIGH** as a database artefact; publish as a table/appendix, not a paper.

### 4.8 The number-field template

* **Theorem NF-1.** For a row over $\mathcal O_K$, $[K:\mathbb Q]=g$, with sharp lcm exponent
  $k$: the criterion is $\frac1g\sum_{v\mid\infty}d_v\log\rho_v>k+\log D$, where $\rho_v$ is
  the radius at $v$ of the *linear-form* generating function. **The whole content of a
  particular theorem is which places admit the fold.**
* **Beukers' Theorem 4 folds at one place only, and that is forced.** At the second place the
  fold constant of the conjugate construction is the image of the *formula* under
  $\sqrt5\mapsto-\sqrt5$, not the Galois conjugate of the number assumed to lie in $K$ — they
  agree only if one already knows what one is trying to prove. Numerically
  $\sqrt{493.9654530919\times0.8351896967}=20.3114464521\ldots>e^3=20.0855369231\ldots$:
  **Theorem 4 is true by 1.12%, and the missing percent is exactly the
  Galois-conjugate-versus-analytic-conjugate gap.** Had the conjugate fold been available the
  margin would be 21%.
* **Theorem NF-0, a complete no-go for complex-fold rows.** A row with $|\lambda_1|=|\lambda_2|$
  has the fold cancel one of two equal-modulus singularities, so $\rho_v=1/|\lambda|$ and the
  linear form does not decay at all. **No Diophantine statement of Beukers' shape is available
  for $L(\chi_5,3)$, $L(\chi_{-3},2)$-via-$\mathbf B$, $\zeta(3)$-via-$\delta$, or either new
  Herfurtner row.**
* **Three corrections to `CDT_FINDER.md`** (see §7), including: the field norm is the wrong
  object, because $\lambda_2\notin K$ — Beukers' row's characteristic roots generate the
  **cyclic quartic field $\mathbb Q(\zeta_{20})^+$**.

*Judgment:* **HIGH** — turning Beukers' Theorem 4 into a reusable template with an explicit
Mode I/II taxonomy, plus the observation that it is true by 1%, is a genuine contribution and
a nice piece of mathematical archaeology.

---

## 5. Conjectures with evidence

| conjecture | statement | evidence | file |
|---|---|---|---|
| **Alignment-prime conjecture** | Two-lattice gain lives exactly at primes dividing the level/conductor of the extension: for $p\mid N$ the modular row is an overconvergent $p$-adic modular form (Calegari's mechanism), the hypergeometric row converges by Lucas/Kummer laws, different rates $\Rightarrow$ large $v_p$ of the cross determinant. **Prediction:** $L(3,\chi_5)$ aligns at $p=5$; $\beta(4)$ at $p=2$; $\zeta(2k+1)$ at **no** prime — so for odd zeta values the two-lattice method is *structurally* unavailable. | Catalan aligns at 2; $L(2,\chi_{-3})$ at 3; $\zeta(5)$ (conductor 1) 2+3 fusion died by a defect theorem; all $p=3$ Catalan slopes were artefacts | `THEORY_NOTES_02.md` §3 |
| **Zagier's conjecture, made precise** | The unconditional classification of §4.1 is exactly Zagier's conjecture that any integral solution has a modular parametrisation = Bombieri–Dwork in rank 2 with four singular points | the reduction is proved; $1.27\times10^6$ + $3.0\times10^6$ scan hits | `paper/classification` |
| **Scalar barrier** | No single scalar pipeline beats $\zeta(2),\zeta(3)$ | E14 bounded scan ($N\le25$); H6's complete $(a,d)$ pass | `paper` §3.4, `HERFURTNER_CLASSIFICATION.md` H6 |
| **P2$'$** | see §6.3 | 39 988 exact instances to $n=10^4$ | `P2_SCALE.md` |
| **F2 binary mixing** | max Smith defect steps by $+2$ exactly at $D=2^k-1$ and $3\cdot2^{k-1}-1$; unique positive exponent $2\lfloor\log_2(D+1)\rfloor-5$ | to $D=34$ | ledger Tier 5 |
| **$\xi_2^{(m)}$ family** | see §1.2 | $2^{3800}$ | `PADIC_IRRATIONALITY_CENSUS.md` §9 |
| **AESZ 207 identity** | $\xi_\infty,\xi_2$ are periods of the $\mathbb Q(\sqrt{17})$ conifold outside all tested catalogues | 1139 / 5958 digits of exclusions | `AESZ207.md` §4 |
| **Rate–Purity Conservation (F07)** | $r=1+\lfloor m/2\rfloor$ links recurrence order to multiplicity of exponent 0 | observed families only | `LEDGER_DUMP.md` F07 |
| **Master $\chi$-twisted valuation law (H21)** | $v_p(p^wB_nA_q-\chi(p)B_qA_n)\ge w$ for all fifteen families | verified broadly, sharp in tested cells; theorem-level only for 4 families single-digit (H17) | `LEDGER_DUMP.md` H21 |
| **BZ $\kappa$ closed forms (F02/F03)** | $\kappa_2,\dots,\kappa_5$ are explicit zeta combinations; $\lambda_5=(514/87)\zeta(5)$ | 434-digit PSLQ; the *relation* $\lambda_5=\kappa_5-\kappa_2\kappa_3$ is **proved** | `LEDGER_DUMP.md` F01–F03 |
| **Brown–Zudilin bridge (B01)** | compact two-variable harmonic formula for the BZ top row $P_n$ | heavily verified, unproved; the flagship open problem of the old ledger | `LEDGER_DUMP.md` B01 |
| **Unidentified shared constant** | $\xi^*=0.7372929961855962$ shared across the six-point rows | a shared constant usually means an unrecognised period | ledger Tier 5 |

---

## 6. Dead ends, with the quantified obstructions (these are contributions)

### 6.1 The two-row lattice method cannot prove irrationality — and *why*, exactly

`CATALAN_AUDIT.md` is, independently of Catalan, the best piece of methodology in the corpus.
The rows, CSVs, 2-adic divisor, congruence lattice and integrality are **all correct**; the
inference is not. Four separable, reusable findings:

1. **The control experiment.** Let $G^*=\mathrm{bestappr}(G,10^{320})$ — a *rational* number.
   The lattice never uses $G$, so it produces the same $c_n$. At $n=98$,
   $|q_nG-p_n|=e^{-78.051}$ and $|q_nG^*-p_n|=e^{-78.050}$, non-zero. **Identical numerics
   from a number that is rational.** So "$|\ell_n|$ is exponentially small and non-zero"
   proves nothing at any finite $n$. *This falsification device should be run on every future
   lattice construction before any claim.*
2. **The pairs are worse than trivial.** A CF convergent with $\log|q|=812$ gives
   $e^{-813}$; the bridge's pair with $\log|q|=1388$ gives only $e^{-78}$.
3. **The precise broken step.** At $F_n<0$ the proof can deliver "$\ell_n\ne0$" *or*
   "$|\ell_n|\le e^{F_nn}$", never both provably; proving $c_1\Lambda_1+c_2\Lambda_2\ne0$
   requires knowing $\Lambda_2/\Lambda_1\notin\mathbb Q$. **Circular.**
4. **Modulus double-use, resolved.** $[\mathbb Z^2:\mathcal K_n]\approx M_n$ *is* genuinely
   true — so the box area equals the covolume, and **Minkowski's first theorem in that box is
   Dirichlet's theorem.** Any construction pushed to $\sigma_n>E_1+E_2$ lands there and loses
   all discriminating power.

**This resolves the apparent $\delta>1$ contradiction in `ZETA3_TWO_LATTICE.md` §16.**

*Judgment:* **HIGH** — write it up as a short "how two-row lattice constructions fail"
methodological note. It would save other people months.

### 6.2 The multislope programme: premise false, cleanly

**The $p$-adic slope is a property of the *row*, not of the individual companion**: within a
row either every companion has a limit at $p$ with the *same* $\sigma_p$, or none does (AESZ
207's three Cauchy increments at $n=1600$ are $19156,19144,19129$). **No row in the corpus has
two companions with limits at two different primes.** The file also self-corrects several of
its own sub-claims (§§443, 819, 824–825) and marks one script superseded. A cleanly killed
programme.

*Judgment:* **MEDIUM** — publish as one paragraph of "things that don't work".

### 6.3 Lemma P2 / P2$'$ — the exact remaining content past the threshold

**What P2$'$ asserts** (as reduced): the cone-minimum/true-minimum ratio
$\rho$ satisfies $\frac1n\log\rho\to0$, equivalently there are infinitely many even balance
indices; equivalently the selected linear form does not vanish.

**Evidence at scale** (`P2_SCALE.md`): **39 988 exact instances to $n=10\,000$**, four values
of $k$. Cone-ratio slope $-2.2\times10^{-6}\pm1.4\times10^{-6}$ — zero, and the bound is
improved **790×**. Run-length spectrum **exactly fair-coin** (longest run 13, matching
$\log_2N$). The quotient *at* the balance index is **not** Gauss–Kuzmin (biased) but is
**independent of parity** — which closes the one remaining P2$'$ failure route. Exact
valuation law $v_2(q\xi_2-p)=(24-k)n$ exactly. Earlier statistics on $13\,332$ pair lattices
give $\frac1m\log\rho=-0.00076\pm0.00082$ and median $\rho=\sqrt2$ to four figures at every
$m$ from 8 to 44 (and $\ell^1/\ell^2\in[1,\sqrt2]$ explains the range).

**Corrections between the three files, explicitly:** "no bad run longer than 5" was a
small-$N$ artefact (`P2_SCALE.md` §39, §522 — it is 13 at $n\le10^4$ and grows like
$\log_2N$); the "$24.02$" was a finite-$n$ artefact (§460); the $i/L$ law constant was
corrected to $0.5215$; `P2_HOLONOMIC.md` corrects a "2-adic residue structure" claim as false
in both halves and the framing "$S_n$ = the 2-adic modulus" as wrong ($S_n=D_{6n}^2$ is the
*denominator*); `P2_STRUCTURE.md` §520 notes `cone80.gp` scanned the wrong index and corrects
51 of 231 cone minima.

*Judgment:* **MEDIUM-HIGH** — the statistics are unusually careful (fair-coin run-lengths,
parity-independence) and the reduction of an irrationality question to a parity statement
about one continued-fraction index is genuinely interesting. Conjecture note.

### 6.4 Other quantified obstructions worth recording

* **$\zeta(7)$ level 60**: the archive's most promising margin (filtered entropy $19/6$ vs
  $\log R=3.4639$) dies on a Koszul depth-1 lift obstruction. Audited: **no Apéry system
  exists** (no AL-eigen purified source on the genus-zero quotient); the archive arithmetic is
  correct but the analytic radius is conditional on seven unproved cancellations; level-12
  parent score $-8.79$.
* **$\zeta(5)$**: level-16 system built, $\xi_2=\tfrac7{32}\zeta_2(5)$ to 372 digits
  (Theorem F in rank 5); Conjecture D for $\zeta(5)$ at $p=2$; the Brown–Zudilin row has **no
  2-adic resource**; a level-12 counterexample to *sufficiency* of the Euler criterion.
  The 2+3 "fusion" died by a **fusion-defect theorem**: $\mathrm{Sym}^{a-1}\otimes\mathrm{Sym}^{b-1}$
  never hits $\mathrm{Sym}^4$.
* **$L(3,\chi_5)$**: the ETA programme is closed (determinant-16 descent collapses, monodromy
  order 768; $V_3$ route obstructed, "$0=2$"). **Correction:** AESZ 184 is *not* a new period
  class — $A^{184}_n=\binom{2n}nA^\eta_n$, $B^{184}_n=\tfrac12\binom{2n}nB^\eta_n$ exactly, the
  $\eta$ row Hadamard-multiplied by $(1-4z)^{-1/2}$; it has no archimedean Apéry limit but a
  complex fold $\tfrac14L(\chi_5,3)+i\tfrac\pi{20}L(\chi_5,2)$.
* **One class, two worlds:** a hypergeometric $L(2,\chi_{-3})$ row (third-integer poles) is
  **proved**; 3-adic alignment with Zagier $\mathbf C$ exact to $3^{199}$ (Theorem F's
  prediction), $\kappa_3=3$, best $\delta=0.9191$; design-rule correction $w+\kappa_p>0$.
* **$\Gamma_1(7)$**: no Cooper freeness; Plücker minors nonzero $\Rightarrow$ scalar wedge
  descent obstructed; $S_3$ rational descent capped at radius $0.148$.
* **Conductor-12 two-prime Catalan row** (new, `TWO_PRIME_HOLONOMY.md` §4): exists, is proved,
  confirms Theorem F at both primes ($\xi_2=\zeta_2(2)$, $\xi_3=\tfrac9{10}L_3(2,\chi_{12})$,
  to $v_2=463$, $v_3=176$) — **and has nothing to bridge against** ($\sigma_3=0$ on Zagier
  $\mathbf E$ and on Zudilin's row), and forcing it into the holonomy bound makes the margin
  $-7.97\to-67.4$.
* **Cusp moves buy no measure gain**: the merged two-sequence bound is
  $\mu\le1+\min_i\sigma_i/\delta_i$ — the better of the two placements, with no
  convex-combination gain at leading order.

---

## 7. Status conflicts, supersessions and corrections

Everything below is a place where the corpus contradicts itself. Most are self-flagged; a few
are not.

| # | conflict | resolution |
|---|---|---|
| 1 | **$\zeta_p$ notation collision.** `MUM_SURVEY.md` §3 uses $\zeta_p(k)=L_p(k,\omega^{1-k})$ (the $\Gamma$-class convention); `EULER_CRITERION.md` §313 uses $\zeta_p(s)=L_p(s,\mathbf 1)$ (Washington). **These disagree at $p=5$, $k=3$** — exactly the flagship's value. | The file itself flags it. **The $\zeta_5(3)$ paper uses the standard Washington/Calegari convention** ("interpolating $(1-p^{-s})\zeta(s)$ at negative integers of the right parity"), so the flagship is safe. But `MUM_SURVEY.md`'s "$-40\zeta_p(3)$" is $L_5(3,\omega^{-2})$, a *different number*. **Do not conflate; fix the notation before publishing either.** |
| 2 | `RIGIDITY_PROOF.md` Theorem R presented as a proof | **Withdrawn by its own header**: it is a *reduction*. The complex-analytic proof of H3 is withdrawn (`PADIC_PERIOD.md`). Superseded by `CONJ_D_PROOF.md` (H3 proved) and `ACF_ONE_SURFACE.md` (H3 eliminated). |
| 3 | `RIGIDITY_PROOF.md` §3's degree-2 cover | An **artefact** of Zagier's level-12 normalisation of $\mathbf F$ (`ACF_ONE_SURFACE.md`). |
| 4 | `CDT_FINDER.md` §4/§6 scoring of $X_1(5)\,\mathrm{Sym}^2$ | **Vacuous** — $\lambda_1$ is double, the companion's linear form does not decay. `HOLONOMY_LINDEP.md` §2.3, "the single most consequential correction here". |
| 5 | `CDT_FINDER.md` §3 normalised number-field budget $\lambda_2^{\rm norm}=|N(\lambda_2)|^{1/g}$ | Correct only in **Mode II**; Beukers' Thm 4 is precisely a case where Mode II's hypothesis is unavailable. And **the field norm is the wrong object**: $\lambda_2\notin K$ (`NUMBER_FIELD_PROGRAM.md`). |
| 6 | `CDT_FINDER.md` §6 "$X_1(5)\ \mathrm{Sym}^2$" identified as Beukers' row | **It is not.** Beukers' row is over $\mathbb Z[\varphi]$, Apéry (R3) shape, $c=1$, roots of $\lambda^2-(248+110\sqrt5)\lambda+1$. |
| 7 | `EULER_CRITERION.md` §4.1 entry for $s_{18}$ | **Must be withdrawn**; repaired by `WEIGHT_DROP.md`. |
| 8 | `SOURCES_S18_ZUDILIN.md` §3.4 and §6.3–§6.5 | **Superseded** by `WEIGHT_DROP.md` and `ZUDILIN_2ADIC.md` respectively; $\xi_2^{\rm Zud}=\zeta_2(2)$ is now **proved**, not open. |
| 9 | `ADELIC_HOLONOMY.md` §4.2 tables; "the entry test flipped on the level-8 Catalan host" | **SUPERSEDED / withdrawn** (§7): the conditional orbit's $y$-slope is $-2$, not $0$; the claimed integrality floor is not available because the descent is not integral. |
| 10 | `SLOPE_CENSUS.md` "$\mathbf F$ vs $s_{18}$ irregular, no clean $5/4$ law" | **False negative of the test statistic**; the clean test gives $v_3(4\xi^{\mathbf F}-5\xi^{s_{18}})=589$ at $N=600$ (`ACF_ONE_SURFACE.md` §7). |
| 11 | `SLOPE_CENSUS.md` "$s_{10}$'s $A_n$ is not integral" | **Artefact of the initial conditions** (`SQRT_APERY.md` §337). |
| 12 | `ROW_LEDGER.md`: $\zeta(5)$ level 12 budget $+5.536$; $\zeta(7)$ level 24 "all slopes zero" | **Both wrong**, self-corrected in its own audit section: the $\zeta(5)$ dominant root is $\alpha$ not $\alpha^4$, so score = budget = $-5$ (and $\prod\lambda_i=-1$ *proves* all $\sigma_p=0$); the $\zeta(7)$ row has $\sigma_2=5\ne0$ forced by $\prod\lambda_i=-32$. Also: `budget = score + headroom` is an **order-2 fact**; in general $\sum_p\sigma_p\log p=\log|\prod_i\lambda_i|$. Conclusions unchanged (both get *worse*). |
| 13 | `MULTISLOPE_PROGRAM.md` §§443, 819, 824–825 | Several self-corrections; premise false; `row1_arch.gp` superseded (its `lindep` output "must not be quoted"). |
| 14 | `GOOD_PRIME_TOWERS.md` | The proposed mechanism "is wrong in every particular"; a fibre statement "is false". Trichotomy survives as a measurement only. |
| 15 | `THEOREM_B_EXACT.md` | Its own §5 supersedes `ver:complex`; the claim of a real fold for all fifteen rows is **false for the (R2) six**. |
| 16 | `P2_SCALE.md` vs `P2_STRUCTURE.md`/`P2_HOLONOMIC.md` | "no run $>5$" and "$24.02$" are finite-$n$ artefacts; $i/L$ constant corrected to $0.5215$; 51 of 231 cone minima corrected; two framings in `P2_HOLONOMIC.md` called wrong. |
| 17 | `PADIC_IRRATIONALITY_CENSUS.md` §6.1 host for $L_3(2,\chi_{-3})$ | Corrected by `PADIC_HOLONOMY_CENSUS.md` §5.4 ($X_0(3)$, not $X_0(9)$). |
| 18 | `MUM_SURVEY.md` §5.5 "AESZ 184 is a new period class" | **Not new** — Hadamard of the $\eta$ row with $(1-4z)^{-1/2}$ (`L3CHI5_TWO_WORLDS.md`). The "drift $5\times10^{-2}$" is non-convergence, not slow convergence. |
| 19 | `MUM_SURVEY.md` §5.6 "$\xi_2$ to 2449 digits" | Artefact of a $0.85\sigma N$ estimate; the true figure is 5958 (`AESZ207.md`). |
| 20 | `ZETA3_TWO_LATTICE.md` §16 apparent $\delta>1$ | **Resolved** by `CATALAN_AUDIT.md` (§6.1 above): at $F_n<0$ the method is Dirichlet's theorem and reproduces itself verbatim for a rational surrogate. |
| 21 | Ledger `D04` (Catalan $0.588274$) | **SUPERSEDED / DO NOT CITE**; use `D03` $=0.581983473693357\ldots$ |
| 22 | `L3CHI5_TWO_WORLDS.md` §295 | Flags an earlier statement as "wrong as stated" — AESZ 184 carries the project's own $\eta$ row's class. |
| 23 | `ONE_CLASS_TWO_WORLDS.md` §344, §402 | `chi3_asym.py` is wrong by a constant (sign change in $Q_m=9\sum_jA_j$); kept as a record, marked "**wrong for $\Lambda$**". |
| 24 | `SOL_NOTES_DIGEST.md` N3 | "$V\simeq\Lambda^2_0W$" is **false as stated**; needs the quadratic twist $\chi=(1-z)^{1/2}$. Everything downstream survives *except possibly* $\deg_{\rm par}F^{2,0}=1/6$, hence $\lambda=1$, hence $\delta_1=0$ and its Theorem B. |
| 25 | `SPORADIC_SCAN2.md` §252 | A normalisation artefact of the scan, not a fact about the curve. |
| 26 | `SOURCES_S18_ZUDILIN.md` §5 read as "Zudilin's row has no geometric model" | **Contradicted** by `CRYSTAL_THEOREM_F.md`: $R_m=16^mQ_m\in\mathbb Z$ satisfies the full Dwork tower (2288 tests, 0 failures), so the row *does* have a Dwork crystal. The crystal file is later and computational; it wins. |
| 27 | `HERFURTNER_CLASSIFICATION.md` §8 item 4 (the "$K3$ window might have been missed" worry) | **Closed** by Proposition K1: $\deg\mathcal J\le6c+4e_3+3e_2-12\le12$ always, so the $\deg\le12$ window was complete. |

**Excluded-file one-liners** (for completeness): `EMN_VERIFICATION.md` — verification of the
Eskandari–Murty–Nemoto Catalan-motive note (12 of 13 claims verified exactly, 1 corrected: §13
is unconditional; the EMN moving-period host is the first Catalan architecture whose CDT
entry is not vacuous, $+0.096$, but self-defeating). `EMN_PROJECTOR.md` — the adelic-projector
search terminates and returns nothing new (Niven collapse of the shift lattice); its by-product
is a genuine $G$-fold at $z=2$ with $\mathrm{Hf}(2)=2G$, priced at entry $-0.85$.
`HADAMARD_HOST.md` — the Hadamard host of (Zudilin, Nesterenko); the one reusable lemma is that
**the Hadamard product does not double $k$** (type stays $[1..6n]^2$), but entry fails by $-11.10$
because the four product points form a single Galois orbit so no pure module exists.

---

## 8. Top ten, ranked by (proved × novel × shareable)

1. **Theorem F, the Euler-factor criterion for $p$-adic Apéry limits**, together with
   **Conjecture D proved for Zagier's B, C, F at $p=3$** ($\xi_3=\tfrac12,\tfrac12,\tfrac58$
   times $\zeta_3(2)$) **and its Dwork-crystal reformulation** (the criterion is really
   Frobenius-slope separation on an extension of a MUM crystal, so it applies to non-modular
   rows — Zudilin's row has a Dwork crystal with tower eigenvalue $\chi_{-4}(p)p^{-2}$).
   Exact criterion, closed-form value, verified to $\ge p^{2991}$, transfers verbatim to
   rank 4, and *explains* the alignment phenomenon instead of observing it. Hypothesis (b)
   open in 4 rows, and the *value* half is modular-only — state both.
   **Theorem paper.** (§3.2, §3.3, §3.9)
2. **The classification of second-order integral Apéry-like rows** — Theorem A (conditional,
   proved, via the twist-parity lemma forcing $\deg\mathcal J=12$), Theorem B (unconditional
   necessary conditions), the cross-ratio invariant reframing "Beauville's six" as the count
   of rational cross-ratios, the $1.27\times10^6$-hit scan, two new elliptic-surface rows, and
   Corollary H6 (only Apéry's $\zeta(2)$ row and Beukers' square root have positive score, a
   *complete* pass). **Already drafted at 21 pp.** (§4.1)
3. **AESZ 207: an erratum against the published AvSZ table** — their printed limit is wrong
   from the 4th significant figure; correct value to 1139 digits, $\xi_2$ to 5958, geometry
   located, sharp exclusion set, validated pipeline. Cheapest high-value thing here. **Short
   note; send this week.** (§4.4)
4. **Improved irrationality measures for known $p$-adic $L$-values** by optimising the
   template over all holomorphic self-maps of the disc: $\zeta_3(3)$ $22.3\to10.0$,
   $\zeta_2(5)$ $19.744\to19.39$, $\zeta_2(3)$ $7.18\to6.3$. Calibrated exactly against CDT's
   and Calegari's printed numbers, and **independent of the $\zeta_5(3)$ novelty audit.**
   (§1.1)
5. **Theorem B\*: the Apéry limit is exactly $L(\Phi,w+1)$**, with the endpoint condition
   $\delta(\varphi)P(w+2)=0$ separating the nine real-fold from the three complex-fold rows,
   and no further rational factor. The archimedean half of the same picture as (1).
   **Theorem paper, with (1).** (§3.1)
6. **$\xi_2^{\rm Zud}=\zeta_2(2)$ proved**, with the exact tail $v_2=8m-1-4s_2(m)$, via a new
   exact identity to Beukers' Padé pair. Self-contained, short, checkable, ideal Lean target,
   and it makes a previously conditional construction unconditional. (§1.5)
7. **The K3 row**: the first Apéry-like sequence whose Apéry constant is a critical value of a
   cusp form (32.3.d.a), with Chowla–Selberg evaluation *and* a proved impossibility statement
   (three-term integral recurrences can never produce one). **Already drafted at 4 pp.** (§4.2)
8. **The holonomy method-lemma cluster**: linear independence costs exactly one generator;
   fold-regularity is codimension $\mu$ so only simple-$\lambda_1$ rows qualify (this voids a
   prior scoring); the holonomy measure stays finite on the whole window
   $-0.8385<\mathrm{score}\le0$ where the classical argument is vacuous, giving
   $\mu\le18.16$ for Beukers' row against $50.65$; and the two proved cancellation identities
   showing a second prime is never worth anything. **Method-lemma paper.** (§2.1–2.4)
9. **The cusp-move theorems and the $\mathbf A/\mathbf C/\mathbf F$ one-surface theorem** —
   Theorem 1 (the move, iff condition), Theorem 2 (the companion transforms via the *new
   Hauptmodul*, not a gauge), Theorem 4 (rigidity: no rational orbit contains a positive-score
   row), the 49-row orbit census, the simultaneous 2- and 3-adic relations
   $\xi_p^{(2)}=-2\xi_p^{(3)}$, and two new rows whose limit is a genuine
   $\mathbb Q$-combination of $\zeta(2)$ and $L(2,\chi_{-3})$. (§3.4, §3.7)
10. **How two-row lattice constructions fail** — the $G^*=\mathrm{bestappr}$ control
    experiment, the "at $F_n<0$ Minkowski *is* Dirichlet" identification, and the precise
    circularity (one can prove $\ell_n\ne0$ *or* bound $|\ell_n|$, never both). Negative, fully
    quantified, immediately useful, and it resolves the corpus's one apparent
    $\delta>1$ contradiction. **Short methodological note.** (§6.1)

**Just outside:** the number-field template + "Beukers' Theorem 4 is true by 1.12%" (§4.8);
the rank-4 MUM survey's two new period classes $L(\chi_8,2)$, $L(\chi_5,3)$ and the
"$p$-adic limit where the archimedean one does not exist" tool (§4.3); the Casoratian sign
theorem and the Hankel-moment decidability criterion (§3.9); Theorems R1–R3 on roots of
$\mathrm{Sym}^w$ rows (§3.6); and P2$'$ at $n=10^4$ (§6.3).

---

## 9. Recommended order of operations

1. **Send AESZ 207** (§4.4). It is done, it is short, it corrects a published table, and it
   costs a day.
2. **Fix the $\zeta_p$ notation collision** (§7 item 1) before anything $p$-adic goes out.
3. **The structure paper**: Theorem B\* + Theorem F + Conjecture D at $p=3$ + the
   $\mathbf A/\mathbf C/\mathbf F$ one-surface theorem (§§3.1–3.4). This is the coherent
   mathematical story and it is largely written in `paper/sections/`.
4. **Release the classification paper** (§4.1) and the K3 note (§4.2) — both drafted.
5. **The method-lemma note** (§2) and **the "how lattices fail" note** (§6.1) — short,
   negative, and the most immediately useful to other people working on this.
6. **Lean**: $\xi_2^{\rm Zud}=\zeta_2(2)$ (§1.5) and Theorem R1 (§3.6) are the two cleanest
   targets outside the Catalan 2-adic cluster.
