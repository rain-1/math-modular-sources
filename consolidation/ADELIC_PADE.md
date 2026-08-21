# Adelic Padé (Problem 7.5): formulation, mechanism, and where it stops

**No claim of irrationality is made anywhere below.** Every number is a finite-$n$
measurement; the `control.gp` argument of `CATALAN_AUDIT.md` §4(a) applies verbatim.

---

## Verdict

1. **Formulated, and the formulation is clean.** The correct two-place problem is *not*
   "impose a 2-adic congruence and an archimedean condition". It is: find $e$ with
   $$\text{(a) } v_2(A_e)\ \ge\ kn\ \ (\text{integrality}),\qquad
     \text{(c) } \varepsilon(e)\ \le\ 24n-\psi(n),\ \psi\to\infty\ \ (\text{separation}),\qquad
     \text{(b) } |A_eG-B_e|\le e^{Fn},\ F<0 ,$$
   where $\varepsilon(e)$ is the *cancellation excess* in the numerator, defined below.
   **(a) and (c) are both $G$-free explicit congruence conditions**, and they are compatible
   precisely because $k_*=22.3513<24$. So the window is real and the sieve is not circular.
2. **Proved (modulo the Lean inputs): the nonvanishing implication, in a sharper and
   $M_n$-free form.** See §3, Lemma AP. It needs only $G_2\notin\mathbb Q$ (Calegari/Beukers)
   and the margin $\to\infty$; the "$v_2(q_n)=o(n)$" of Sol's route 1 is *not* the right
   condition — the right one is $v_2(G_2-p_n/q_n)\to\infty$, and Sol's $(24-k)n$ is exactly
   the lower bound for it. Verified: the identity
   $$\boxed{\ v_2\!\left(G_2-\tfrac{p_n}{q_n}\right)\;=\;\bigl(24n-1-4s_2(3n)\bigr)-\varepsilon(c)\ }$$
   holds **exactly**, at every $n\in\{4,6,8,10,12,14,16,20\}$, for the LLL-shortest vector.
3. **Computed: the entire 2-adic structure of the moment family, in closed form** (§2).
   With $m=3n$, $M=m+j$,
   $$v_2(A_{n,j})=3-2(m+M)+s_2(m)+s_2(M),\quad v_2(B_{n,j})=v_2(A_{n,j})-1,\quad
     v_2(A_{n,j}G_2-B_{n,j})=5-v_2(A_{n,j}),$$
   hence $v_2\bigl(G_2-B_{n,j}/A_{n,j}\bigr)=4(m+M)-1-2s_2(m)-2s_2(M)=\mathbf{24n+4j}-1-\dots$.
   Verified with **zero exceptions** for all $(n,j)$, $n\le5$ (all $3n+1$ moments each).
   Specialisations: $j=0$ gives $24n-1-4s_2(3n)$ (= Lean `val_ratio_cross_N`), $j=n$ gives
   $28n-1-2s_2(n)-2s_2(3n)$ (= Lean `nestJform_v2` + `nest_four_B_v2`). **New:** the
   2-adic approximation quality of the $j$-th moment *grows linearly in $j$*, so interior and
   high moments are 2-adically strictly better than the Zudilin/Nesterenko pair, and the
   admissible window widens from $k<24$ to $k<24+4j_0/n$.
4. **Computed: the "24n 2-adic bridge" is an ultrametric triviality.** For every pair,
   $$v_2(h_{j_0,j_1})=\min\bigl(v_2(a_{j_1})+v_2(a_{j_0}G_2-b_{j_0}),\ v_2(a_{j_0})+v_2(a_{j_1}G_2-b_{j_1})\bigr)
   =2v_2(d)+5-\bigl|v_2(A_{j_0})-v_2(A_{j_1})\bigr| ,$$
   with **equality at every one of the $\binom{3n+1}{2}$ pairs, $n\le4$** (0 mismatches).
   The mixed minor carries *no arithmetic information beyond the denominators of $A$*.
   This is the cleanest form of the Absorption theorem (§6.2 of the paper) yet found.
5. **Open, and unchanged: (b).** The route does *not* die of circularity, and it does *not*
   need positivity — but the reduction it achieves is the same one the cone achieves: it makes
   nonvanishing free and leaves an archimedean-size lemma of irrationality strength (§4).
   Its version, **Lemma P2′**, is formally *weaker* than P2 (no cone, only removal of an
   exponentially thin sublattice) — that is the one genuine gain.
6. **Refuted: the Thue–Siegel determinant idea, structurally.** The moment family has
   $\mathbb Q$-rank **2** — verified, $\operatorname{rank}=2$ and an
   *identically-vanishing* relation space of dimension $3n-1$ for $n\le4$. There is no third
   independent period, so there is no $3\times3$ approximant determinant: every "system of
   three approximants" is degenerate and its determinant is $0$. §5.
7. **Route status: alive but not advanced.** Nothing here is closer to a theorem than
   `CATALAN_POSITIVITY.md` was. What is new is (3), (4), (6) and the exact form of the
   nonvanishing criterion.

---

## 1. Mechanism: what Padé problem, at which place

**Archimedean.** Both rows are the Beukers $\zeta(2)$ double integral at half-integer
parameters. Zudilin (`papers/0201024v3.pdf` Thm 3, eq. (6)):
$u_nG-v_n=\frac{(-1)^n}{4}\iint x^{n-1/2}(1-x)^ny^n(1-y)^{n-1/2}(1-xy)^{-n-1}$, and the
Nesterenko $(4,7)$ row is the same integrand at index $3n$ times $w^n$, $w=xy/(1-xy)$
(Mechanism A). The Hermite–Padé data is the *very-well-poised* rational function
$$R_n(t)=n!(2t+n+1)\frac{(t-n+1)_n\,(t+n+1)_n}{\bigl((t+\tfrac12)_{n+1}\bigr)^{3}},$$
a triple pole at each of the $n+1$ half-integers $t=-k-\tfrac12$, simple zeros at
$t=0,\dots,n-1$ and $t=-n-1,\dots,-2n$, $R_n(t)=O(t^{-2})$. The *contact conditions are the
$2n$ zeros*; the Padé/Rivoal reading is the one in `0210423v1.pdf` §1.

**2-adic — and this is the load-bearing correction.** There is **no independent 2-adic
interpolation condition.** From the Lean development (`NesterenkoHypergeometric.lean`,
`NesterenkoWeightedBridge.lean`, `Zudilin.lean`, `Limits.lean`) the mechanism is:

* $\mathcal G_2$ is *defined* as the 2-adic limit of the Zudilin ratios
  (`Limits.lean:62`, `exists_GZ2`), not as a Beukers/Calegari $L$-value; the identification
  with $\xi_{\rm Cat}/8$ and with the level-8 modular constant is separate
  (`Beukers.lean:116`, `CrossIdentity.lean:220`).
* The Zudilin rate is a **telescoping** argument: the Wronskian
  $8(k+1)^2(2k+1)^2W_k=(-1)^k\varphi(k+1)$ (`Zudilin.lean:297`) plus the exact
  $v_2(Q_m)=-4m+2s_2(m)$ (`:206`) give the increment valuation $8m-5-4s_2(m)-4v_2(m)$
  (`val_ratio_diff`, `:375`), strictly increasing, hence convergence *and* the exact tail
  $v_2(\mathcal G_2-P_m/Q_m)=8m-1-4s_2(m)$ (`val_GZ2_sub`, `Limits.lean:71`). Reproduced
  exactly here for $m\le20$ (`g2.gp`).
* The Nesterenko rate comes from the **beta prefactor** $4^{7n+2}(8n)!(6n)!((7n+1)!)^2/((14n+2)!)^2$
  whose valuation is $14n+2-s_2(n)-s_2(3n)$ by Legendre/Kummer (`nestPrefactor_v2`), times a
  ${}_3F_2$ that is a 2-adic **unit** ($\|{\rm term}_k\|_2\le2^{-k}$, `norm_nestHyperTerm_le`,
  because each half-odd Pochhammer factor has norm exactly $2$ and $(4n+1)_k^2/k!$ is an
  integer). The factor $(14n+2)!^2$ in the denominator exists only because the numerator
  $\prod_{r<4n}(t-r)$ of $R_n$ forces the series to start at $k=4n$.

**Consequence.** *The 2-adic slope is the archimedean order of vanishing, read through
Legendre's formula.* There is no free 2-adic contact order to tune: raising the 2-adic rate
requires raising the archimedean vanishing order, which raises the archimedean size in
lockstep. This is exactly why §2 gives quality $24n+4j$ — one number, one parameter — and why
the $F$-improvement of interior pairs in `CATALAN_EXPLICIT.md` §5 came from denominators, not
from $v_2(h)$. **A genuinely two-place Padé problem does not exist in this family.**

*(Caveat on sources: the 5:8 paper's "Rivoal's moving Padé specialisation / Beukers'
$\mathcal T(x)=2R(x)-R(x/2)$, $x_m=\tfrac12-m$" paragraph is headed "Idea of proof" and is a
plausibility narrative: $R(x)$ is never defined there, the Padé identity is never stated, and
the $8m$ rate is not derived from the moving point — it is asserted from
$v_2(Q_m)=-4m+2s_2(m)$, which is stronger than anything in either Zudilin paper. Only the Lean
route above is load-bearing.)*

## 2. Closed forms for the moment family ( `tj.gp`, `lawcheck.gp` )

$M_{n,j}=\iint K_Z^{(3n)}w^j=A_{n,j}G-B_{n,j}$, $0\le j\le 3n$. Put $m=3n$, $M=m+j$.
Verified **exactly, no exceptions, for every $(n,j)$ with $n\le5$** ($\sum(3n+1)=65$ moments),
$\mathcal G_2$ known mod $2^{351}$:

| quantity | closed form |
|---|---|
| $v_2(A_{n,j})$ | $3-2(m+M)+s_2(m)+s_2(M)$ |
| $v_2(B_{n,j})$ | $v_2(A_{n,j})-1$  (this proves Lemma E1 of `CATALAN_EXPLICIT.md` §6) |
| $\tau_j:=v_2(A_{n,j}\mathcal G_2-B_{n,j})$ | $5-v_2(A_{n,j})$ |
| 2-adic quality $v_2(\mathcal G_2-B_{n,j}/A_{n,j})$ | $4(m+M)-1-2s_2(m)-2s_2(M)=24n+4j-1-2s_2(3n)-2s_2(3n+j)$ |

Derivation of the first line: the general moment is the Nesterenko integral with
$(m,M)=(3n,4n)$ replaced by $(m,m+j)$, so its prefactor is
$4^{m+M+2}(2M)!(2m)!((m+M+1)!)^2/((2m+2M+2)!)^2$; Legendre gives
$v_2=2(m+M)+2-s_2(m)-s_2(M)$ for the form and $5-\text{that}$ for $A$. This is a **proof
sketch, not a proof** — the Lean file proves the case $j=n$ only.

The relation $\tau_j=5-v_2(A_{n,j})$ is the sign flip that drives everything: $A_{n,j}$ is
2-adically huge ($v_2\approx-13n$) exactly to the extent that the form is 2-adically tiny.

## 3. The nonvanishing lemma, in its $M_n$-free form

> **Lemma AP.** Let $e_n\in\mathbb Z^{3n+1}$ with $A_{e_n}:=\sum_je_{n,j}A_{n,j}\ne0$, and set
> $B_{e_n}=\sum_je_{n,j}B_{n,j}$. If
> $$\mu(e_n):=v_2\bigl(A_{e_n}\mathcal G_2-B_{e_n}\bigr)-v_2\bigl(A_{e_n}\bigr)
>   \;=\;v_2\!\left(\mathcal G_2-\frac{B_{e_n}}{A_{e_n}}\right)\ \longrightarrow\ +\infty,$$
> then $A_{e_n}G-B_{e_n}\ne0$ for all large $n$.
>
> *Proof.* Otherwise $G=B_{e_n}/A_{e_n}\in\mathbb Q$ for infinitely many $n$; $G$ is a *fixed*
> number, so that ratio equals one fixed $a/b$ along the subsequence, and
> $v_2(\mathcal G_2-a/b)=\mu(e_n)\to\infty$ forces $\mathcal G_2=a/b\in\mathbb Q$,
> contradicting Calegari/Beukers. $\square$

Three things to note.

* **Sol's route 1 is this, with the wrong sufficient condition.** "$v_2(q_n)=o(n)$" is neither
  necessary nor quite sufficient as stated; the right invariant is $\mu$, which does not
  mention $M_n$, $S_n$, $T_n$ or $k$ at all. Sol's $\mu\ge(24-k)n-O(\log n)$ is the correct
  *lower bound* for $\mu$ in the two-row case (see the identity below).
* **The criterion is exactly kernel-exclusion, and is not vacuous.** If $G=a/b$ then the
  rationality-kernel vector has $A_eG-B_e=0$, so $A_e\mathcal G_2-B_e=A_e(\mathcal G_2-a/b)$ and
  $\mu=v_2(\mathcal G_2-a/b)=O(1)$: the kernel vector is precisely characterised by a *bounded*
  margin. So $\{\mu\ge C\}$ is a $G$-free congruence sublattice containing the kernel vector.
* **Cancellation excess.** In the two-row lattice write $\varepsilon(c)=v_2(c_ZX_n+c_NV_n)
  -\min\bigl(v_2(c_ZX_n),v_2(c_NV_n)\bigr)$. Since $X_n$ has 2-adic quality $q_Z=24n-1-4s_2(3n)$
  and $V_n$ has $q_N=28n-1-2s_2(n)-2s_2(3n)$, ultrametricity gives
  $\mu(c)\ge\min(q_Z,q_N)-\varepsilon(c)=q_Z-\varepsilon(c)$. Measured (`eps.gp`, $k=22.4$, LLL
  vector; $\mathcal G_2$ mod $2^{1107}$):

  | $n$ | $v_2(M_n)$ | $\varepsilon$ | $\varepsilon/n$ | $q_Z$ | $\mu$ | $q_Z-\varepsilon$ | $v_2(q_n)$ |
  |---|---|---|---|---|---|---|---|
  | 4 | 97 | 79 | 19.75 | 87 | **8** | 8 | 1 |
  | 8 | 189 | 167 | 20.88 | 183 | **16** | 16 | 1 |
  | 12 | 280 | 255 | 21.25 | 279 | **24** | 24 | 1 |
  | 16 | 370 | 345 | 21.56 | 375 | **30** | 30 | 1 |
  | 20 | 460 | 432 | 21.60 | 463 | **31** | 31 | 1 |

  The inequality is an **equality at every $n$**, and $\mu\approx(24-k)n$ as predicted.
  A 2600-point ball of short vectors (`sieve.gp`, $n\le16$) contains **no** vector with
  $\mu\ge q_Z$: the minimum $\mu$ over the ball is $-2,2,5,3,14,12,18$ at $n=4,\dots,16$,
  i.e. already $\gtrsim1.1n$ and rising. **Every short vector of $\mathcal K_n$ currently
  satisfies the sieve.** (Which proves nothing: for $G^*=\mathrm{bestappr}(G,10^{320})$ the
  kernel vector has length $\asymp b\,e^{2Fn}$ and is simply not in the ball yet.)

**The two-place problem, stated.** Find $c\in\mathcal K_n$ with
(a) $v_2(c_ZX_n+c_NV_n)\ge v_2(M_n)$ [integrality — the "2-adic contact" condition, order
$\lfloor kn\rfloor+v_2(S_n)$]; (c) $\varepsilon(c)\le q_Z-\psi(n)$ [separation]; (b)
$|c_Z\Lambda_Z+c_N\Lambda_N|\le e^{Fn}$. (a) and (c) are congruences on $c$; (b) is the only
condition that sees $G$. **One congruence suffices** for (a) because $T_n\mid a^Z_nU_n-a^N_nY_n$
(Lean `dvd_reduced_cross_N`, $2^{24n}$ divides the reduced cross determinant), so the $p$-side
integrality follows from the $q$-side.

## 4. Where it stops (this is the honest part)

$\{c:\varepsilon(c)\ge q_Z-C\}$ is a sublattice $\mathcal B_n\subset\mathcal K_n$ of index
$\approx2^{(24-k)n}=e^{1.14n}$ at $k=k_*$. Because $\mathcal B_n$ is a *group*, if the first
minimum $\lambda_1\in\mathcal B_n$ then some vector of length $\le\lambda_2$ lies outside it.
Hence:

> **Lemma P2′ (adelic).** For infinitely many $n$,
> $\min\{|c_Z\Lambda_Z+c_N\Lambda_N| : c\in\mathcal K_n\setminus\mathcal B_n\}\le
> e^{o(n)}\sqrt{\operatorname{covol}}$.

P2 $\Rightarrow$ P2′ (a positive-cone vector with $\mu\to\infty$ is outside $\mathcal B_n$),
and P2′ is formally weaker: it removes a thin sublattice instead of half the plane.
**That is the entire gain of the adelic formulation.** It is still of irrationality strength:
if $G=a/b$, the kernel vector is in $\mathcal B_n$ and, once $n>\log b/|F|$, is the first
minimum with $\lambda_1\asymp b\,e^{2Fn}$; then by Minkowski's second theorem the sieved
minimum is $\asymp\operatorname{covol}/\lambda_1\asymp1/b$, which does not tend to $0$.
So proving P2′ proves $\mathrm{den}(G)\ge e^{cn}$ — legitimate, not circular, but not generic:
any proof must bound $\lambda_1$ **from below**, i.e. must supply an irrationality-measure
statement for $G$. Confirming the task's expectation: **the classical Thue–Siegel step reduces
here to "the second minimum is small", and that is available only when the first minimum is not
too small.**

A product-formula ("truly adelic") substitute does not help: for $G=a/b$ and
$N(e)=aA_e-bB_e$, $\prod_v|N|_v=1$ gives a lower bound on $|N|_\infty$ that involves $b$, which
is unbounded. This is exactly why $\mathcal G_2\notin\mathbb Q$, and not the product formula, is
the input in Lemma AP.

## 5. Three or more approximants: structurally impossible ( rank check )

For a determinant argument one needs $r+1$ forms in $r+1$ independent quantities. Here every
$M_{n,j}$ lies in $\mathbb Q+\mathbb QG$: verified,
$\operatorname{rank}_\mathbb Q\begin{pmatrix}A_{n,0}&\cdots&A_{n,3n}\\B_{n,0}&\cdots&B_{n,3n}\end{pmatrix}=2$
for $n=1,2,3,4$, with an **identically-vanishing relation space of dimension $3n-1$**. Hence:

* Any structured $(\ge3)\times(\ge3)$ linear system built on $\{M_{n,j}\}$ is singular; its
  determinant is $0$, not a small nonzero number. The "Hermite–Padé matrix" of the family is
  rank 2 and its kernel is mostly *trivial* relations $A_e=B_e=0$, which are useless
  (they are also automatically inside $\mathcal B_n$, $\varepsilon=\infty$ — a consistency
  check on the sieve).
* The only surviving determinant is the $2\times2$ minor $h_{j_0,j_1}$, and by §4 of the
  Verdict its valuation is a pure denominator statement. A nonzero $2\times2$ determinant says
  only that the two forms are not *both* zero — precisely `CATALAN_AUDIT.md` §4(c).
* Gaining a third row would require a *third period* independent of $1,G$ over $\mathbb Q$
  (e.g. $\pi$, $\log2$, or an $L$-value), i.e. a different construction entirely.
  $\mathcal G_2$ cannot serve: it lives in $\mathbb Q_2$, not in a $\mathbb Q$-determinant.

## 6. Analytic estimate for the explicit $e$ (task 4)

No closed-form $e$ exists to estimate, and §5 explains why: the conditions (a),(c) are
congruences on a rank-2 image, so they do not determine $e$ — they determine a coset structure,
and $e$ remains a congruence-selected object. Characterisation of the actual selected vectors
(from `CATALAN_EXPLICIT.md` §4 plus `eps.gp` here): $n^{-1}\log|c_Z|\to14.0$,
$n^{-1}\log|c_N|\to12.9$, $v_2(c_Z)\in\{0,\dots,4\}$, $v_2(c_N)\ge5$, $v_2(q_n)\in\{1,2\}$,
$\varepsilon/n\uparrow$ toward $k$ from below ($19.75\to21.60$ over $n=4\to20$ at $k=22.4$),
and $\mu=q_Z-\varepsilon\approx(24-k)n$. Nothing matched a recognisable closed form.
The one *provable* scale statement is the balance law $c_N/c_Z\approx\Lambda_Z/\Lambda_N$
(`CATALAN_EXPLICIT.md` §4); the congruence class is P2/P2′ and is not prescribed.

Against the $F<0$ requirement: at $k=22.4$, $F=-0.0169$, and the measured
$n^{-1}\log|{\rm form}|$ for the LLL vector runs $-1.55,-0.89,-0.80,-0.54,-0.59$ at
$n=4,8,12,16,20$ — negative but drifting up to $F$ from below, as $\kappa_n\uparrow\sigma$.
Unchanged from `CATALAN_POSITIVITY.md` §3.2.

## 7. What to do next, if anything

* **Cheap and provable.** Prove the §2 closed forms for all $(m,j)$ (generalise
  `nestPrefactor_v2` from $(3n,4n)$ to $(m,m+j)$). This upgrades Lemma E1 and gives the whole
  family's 2-adic data unconditionally, and it widens the admissible $k$-window to
  $k<24+4j_0/n$ for a pair with lower index $j_0$ — combine with the interior-pair $F$ gain of
  `CATALAN_EXPLICIT.md` §5, which was $\approx2$ per $n$ and came from denominators.
* **Do not pursue** three-row determinants in this family (§5), a two-place Padé problem with
  independently tunable contact orders (§1), or any further valuation sieve as a route past
  P2′ (§4).
* The only live question is still an *archimedean* one: an explicit short vector, or a lower
  bound on $\lambda_1(\mathcal K_n)$.

Scripts: `lattice/adelic_pade/{g2.gp, tj.gp, lawcheck.gp, sieve.gp, eps.gp}` (see the README
there for invocation); longest single run $\approx110$ s.
