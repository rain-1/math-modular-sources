# Higher rank, several companions: multi-slope Dwork crystals and the two-period question

*Claude (Fable), 2026-08-22.  Scripts: `lattice/multislope/`.
Sources read: `consolidation/CRYSTAL_THEOREM_F.md`, `DWORK_CRYSTALS_PRIMER.md`,
`AESZ207.md`, `GOOD_PRIME_TOWERS.md`, `ZETA5_TWO_ROW.md`, `ZETA7_LEVEL60.md`,
`SLOPE_CENSUS.md`, `ROOT_ROWS.md`, `paper/sections/05_two_row.tex` §5.
Tags: **[proved]**, **[verified over range R]**, **[numerical, D digits]**,
**[open]**, **[failed]**.*

---

## 0. Verdict

1. **The programme's premise is false, and cleanly so.**  The $p$-adic slope is
   a property of the **row**, not of the individual companion: within a row
   either every companion has a limit at $p$ with the *same* $\sigma_p$, or none
   does (AESZ 207's three Cauchy increments at $n=1600$ are $19156,19144,19129$).
   **No row in the corpus has two companions with limits at two different
   primes.**  §2.5(C).
2. **What replaces it is better: AESZ 207 is a genuine *multi-period* row.**  Its
   three companion limits are $\mathbf Q$-linearly independent together with $1$
   **at both places** — no relation of height below $\approx10^{200}$
   archimedeanly (1118 digits) or $\approx10^{1400}$ $2$-adically (19129 digits),
   and none over $\mathbf Q(\sqrt{17})$.  The rank-four extension really carries
   three periods.  §2.2.
3. **AESZ 207 is a Dwork row, and satisfies a supercongruence.**
   $A_{mp^s}\equiv A_{mp^{s-1}}\pmod{p^{2s}}$ ($p$ odd) and $\pmod{2^{3s}}$,
   $5985+5985$ tests, $0$ failures.  Since the bare Dwork half is what a
   Laurent-polynomial model gives, this is evidence that AESZ 207 *has* such a
   model — hypothesis (F1) of Theorem F$'$.  §4.2.
4. **Three structural statements of the rank-two theory do not survive.**
   $\sigma_p=\mu_{\det}-2\mu_1$ and $\sigma_p=v_p(c)$ fail (three
   counterexamples); the tower exponent $w+1$ is **not** the denominator
   exponent $k$ ($\mathrm{Sym}^3\mathbf E$: $k=4$ but $w+1=2$ at $p=2$ and $4$ at
   $p=3$); and my own Newton-polygon existence criterion is **refuted** by the
   $\zeta(5)$ level-16 row.  Only Proposition 3.1 carries over.  **At rank
   $\ge3$, $\sigma_p$ is measured, not derived.**  §2.0.
5. **No supercongruence for the companions at slope primes.**  With
   $E(m,s)=v_p(x_{mp^s}-\xi_pa_{mp^s})-\sigma_pmp^s$, $E$ is bounded on every
   row tested and *exactly constant* for AESZ 207 — Conjecture MS-1 gives the
   law $v_2(x^{(j)}_n-\xi_2^{(j)}A_n)=12n+c_j$ for even $n$, $c_1=-4$,
   $c_3=-28$, with a textbook slope-separation interference at $v_2(n)=4$ for
   $j=2$.  The modulus is enormous but grows **linearly**, not with an extra
   $p^s$.  §4.3.
6. **The design question is answered negatively twice.**  Zagier $\mathbf F$
   ($c=72=2^33^2$) is the *unique* row in the whole second- and third-order
   census with two slope primes, and at order two the multi-prime rate is
   already counted by $\log|c|$ — it can only be *lost*, never gained, and
   $\mathbf F$'s $\Lambda/\lambda=9/8$ caps any construction at
   $\delta\le0.111$.  At rank $\ge3$ the product formula that made the design
   rule work is gone, and all four higher-rank rows have deeply negative
   harvestable budgets ($-5.4$ to $-8.6$).  §3.
7. **Two new tower conjectures**, and the first instance in the project where a
   *later* companion carries information the first does not: for
   $\mathrm{Sym}^3\mathbf E$ at $p=3$, $X^{(3)},X^{(4)},X^{(5)}$ have a tower
   limit and $X^{(1)},X^{(2)}$ do not.  §4.4.
8. Corrections forced on `SLOPE_CENSUS.md` §5, `CRYSTAL_THEOREM_F.md` §2.2/§3,
   `ZETA5_TWO_ROW.md` §0.2 and `AESZ207.md` §5 are collected in §5.

---

## 1. The framework: one row, several extension classes

### 1.1 Companions are extensions, and there are $r-1$ of them

Let $L=\sum_{i=0}^{r}z^iP_i(\theta)$, $\theta=z\,d/dz$, be MUM of rank $m$ at
$z=0$ (so $P_0(\theta)=\theta^m$), with coefficient recurrence
$$\sum_{i=0}^{r}P_i(n-i)\,u_{n-i}=0 .$$
The recurrence has order $r$ **in $n$**, which is the $z$-degree of $L$ and has
nothing to do with the rank $m$ of $L$ (AESZ 184: $m=4$, $r=2$; AESZ 207:
$m=4$, $r=4$; Apéry $\zeta(3)$: $m=3$, $r=2$).  Set

* $A$: $a_0=1$, $a_n=0$ for $n<0$, recurrence for $n\ge1$;
* $X^{(j)}$ for $1\le j\le r-1$: $x^{(j)}_i=\delta_{ij}$ for $i\le j$, recurrence
  for $n>j$.

**Lemma 1.1.** **[proved]** $L\bigl(\sum_nx^{(j)}_nz^n\bigr)=j^{m}z^{j}$.

*Proof.* The $z^n$-coefficient of $L(\sum x_nz^n)$ is $\sum_iP_i(n-i)x_{n-i}$,
which vanishes for $n>j$ by construction and for $n<j$ because all $x_i=0$
there; at $n=j$ only the $i=0$ term survives, giving $P_0(j)x_j=j^m$. $\square$

So $X^{(j)}$ is the period of an extension
$$0\longrightarrow\mathcal M\longrightarrow\widetilde{\mathcal E}_j
\longrightarrow\mathbf Q(0)\longrightarrow0,$$
$\mathcal M$ the Picard–Fuchs (iso)crystal, with class represented by $z^j$.
**A row of recurrence order $r$ therefore carries $r-1$ a priori different
extension classes, hence up to $r-1$ different periods at every place.**  This
is the precise sense in which "the extension has two periods" at rank 4: the
usual $B=X^{(1)}$ is only the first of them.

### 1.2 Slopes

Following `CRYSTAL_THEOREM_F.md` §3, for each companion put
$$\mu_1=\lim_n\frac{v_p(a_n)}n,\qquad
\sigma_p\bigl(X^{(j)}\bigr)=\lim_n\frac{v_p\bigl(x^{(j)}_{n+1}/a_{n+1}-x^{(j)}_n/a_n\bigr)}n .$$
Proposition 3.1 there gives $\sigma_p(X^{(j)})=\mu_{\det}^{(j)}-2\mu_1$ where
$\mu^{(j)}_{\det}$ is the slope of the $2\times2$ Casoratian
$a_{n+1}x^{(j)}_n-a_nx^{(j)}_{n+1}$, and $\xi_p^{(j)}:=\lim_nx^{(j)}_n/a_n$
exists in $\mathbf Q_p$ iff $\sigma_p(X^{(j)})>0$.  At order $r=2$ the
Casoratian is the full Wronskian and $\sigma_p=v_p(c)$; **at order $r\ge3$ the
$2\times2$ Casoratian of the pair $(A,X^{(j)})$ is not the full Wronskian and
$\sigma_p$ must be measured, one companion at a time.**  This is the first
technical point at which rank $\ge3$ genuinely differs.

Independently, the **Frobenius slopes** at $p$ are the Newton polygon of the
characteristic polynomial $\chi(x)=\sum_i(\mathrm{lc}\,P_i)\,x^{r-i}$, i.e. the
multiset $\{v_p(\lambda_i)\}$; the product formula
$\sum_p\sigma_p\log p=\log\bigl|\prod_i\lambda_i\bigr|$ holds only in Zagier
normalisation (`L3CHI5_TWO_WORLDS.md` §1.3 gives AESZ 184 as a counterexample,
failing by exactly $4\log2$).

### 1.3 The period rank, and the multi-slope question

> **Definition.**  The **period rank** of a row at a place $w$ (archimedean or
> $p$-adic) is
> $$\operatorname{prk}_w:=\dim_{\mathbf Q}\operatorname{span}_{\mathbf Q}
> \bigl\{1,\ \xi_w^{(1)},\dots,\xi_w^{(r-1)}\bigr\}-1 .$$

$\operatorname{prk}_w\le r-1$ always; $\operatorname{prk}_w=0$ means every
companion limit is rational, $\operatorname{prk}_w=1$ means all companions
compute *the same* period up to $\mathbf Q$-multiples (so the higher companions
are redundant), and $\operatorname{prk}_w\ge2$ is the interesting case: the row
genuinely carries several periods at one place.

There are two reasons to expect collapse.  (i) The classes $[z^j]$ live in a
$\mathbf Q$-vector space $\operatorname{Ext}^1$ computed from the de Rham
cohomology of $L$; if that space is one-dimensional the $\xi^{(j)}$ are
$\mathbf Q$-proportional.  (ii) Even when the classes are independent, the
syntomic regulator is one $\mathbf Q_p$-linear functional; the images can still
degenerate.  §2 answers this experimentally.

**The design question (§3).**  Two-row lattice constructions
(`paper/sections/05_two_row.tex` Theorem E) harvest a *hidden divisor* rate
$$G=\sum_p\min\bigl(\sigma^{\rm eng}_p\alpha,\ \sigma^{\rm dec}_p\gamma\bigr)\log p,$$
and the sum over $p$ is legitimate (the Smith-normal-form step needs only
$T_n\mid h_n$).  So a row with slopes at **two different primes** contributes
$G$ from both.  The multi-slope hope is therefore concrete:

> *Can one row supply a decayer/engine pair from **two of its own companions**,
> with the two $p$-adic resources at two different primes adding in $G$?*

The master formula
$$F=\tfrac12\Bigl[\underbrace{k\max(\alpha,\gamma)+\eta+\alpha\log\rho_2^{\rm eng}}_{\rm COST}
-\underbrace{\bigl(\gamma\log\tfrac1\lambda+G\bigr)}_{\rm RESOURCE}\Bigr],\quad
H=F+\gamma\log\frac{\Lambda_{\rm dec}}\lambda,\quad \delta=1-\frac FH$$
and the order-$m$ score of `05_two_row.tex` eq. (orderm)
$$\text{score}=\log\frac1{|\lambda_2|}-k,\qquad
\text{harvestable budget}=\text{score}+\log\Bigl|\prod_i\lambda_i\Bigr|$$
are the two numbers to compute for every candidate pair.

---

## 2. The rows, measured

*(scripts `lattice/multislope/row1_*.gp`, `row2_*.gp`, `row3_*.gp`, `row4_*.gp`,
`mp_primesupport.gp`, `mp_twoprime.gp`, logs alongside.)*

### 2.0 What governs existence at higher rank — and what does not

The four rows were chosen to test whether the rank-two theory of
`CRYSTAL_THEOREM_F.md` §3 transports.  It half does, and the half that fails
fails hard.  Summarised first, evidence in §2.1–§2.4.

**(a) Proposition 3.1 carries over verbatim.**  $v_p(x_{n+1}/a_{n+1}-x_n/a_n)=n\sigma_p+O(\log n)$,
convergence iff $\sigma_p>0$, and then $v_p(\xi_p-x_n/a_n)=n\sigma_p+O(\log n)$.
Its proof uses only the $2\times2$ Casoratian of the pair $(a,x)$ and is
rank-independent.  Confirmed in every cell that has a slope: AESZ 207's
increments $573,1170,2367,4765,9560,14362,19156$ at $n=50,\dots,1600$ are $12n$
with residual $\pm27$ over the whole range.  **[verified]**

**(b) $\sigma_p=\mu_{\det}-2\mu_1$ does NOT carry over**, because at rank $r>2$
the $2\times2$ Casoratian of $(a,x)$ is a **minor** of the system, not its
determinant, and minors have no product formula.  With
$\mu_{\det}=v_p(\prod\text{roots})$ and $\mu_1=\lim v_p(a_n)/n$ ($=0$ in all four
rows, since $v_p(a_n)=O(\log n)$ throughout):

| row | $v_2(\prod\lambda_i)$ | rank-normalised $/r$ | measured $\sigma_2$ |
|---|---|---|---|
| $\mathrm{Sym}^3\mathbf E$ | $15$ | $2.5$ | $\mathbf 0$ |
| AESZ 207 | $48$ | $12$ | $\mathbf{12}$ |
| $\zeta(7)$ L24 | $35$ | $1$ | $\mathbf 0$ |

Neither column is the rule.  Correspondingly **Corollary 3.2
($\sigma_p=v_p(c)$) does not extend**: its naive higher-rank analogue
$\sigma_p=v_p(\prod\lambda_i)$ predicts $15,48,35$ against measured $0,12,0$ —
**three independent counterexamples**.  (None of these rows is in $(R2)$/$(R3)$
normalisation, which is exactly the hypothesis Cor. 3.2 carries.)

**(c) A Newton-polygon criterion was conjectured and then refuted.**  On three
rows the pattern "the $p$-adic limit exists iff the lowest Newton slope of
$\chi$ at $p$ is simple" holds:

| row | $\chi$ | slopes at $p=2$ (mult.) | lowest simple? | limit? |
|---|---|---|---|---|
| $\mathrm{Sym}^3\mathbf E$ | $(x-4)^3(x-8)^3$ | $2^{(3)},3^{(3)}$ | no | **no** |
| AESZ 207 | $(x-53248)^2(x^2{+}89344x{-}2^{24})$ | $8^{(1)},12^{(2)},16^{(1)}$ | yes | **yes**, $\sigma_2=12$ |
| $\zeta(7)$ L24 | $q_5^{\,7}$, $q_5=(x-1)(x^2{-}8x{+}8)(x^2{+}4x{-}4)$ | $0^{(7)},1^{(14)},\tfrac32^{(14)}$ | no | **no** |
| **$\zeta(5)$ level 16** | genuine roots $-2,-4,2\pm4\sqrt2$ | lowest $1$, mult. $\ge3$ | no | **YES**, $\sigma_2=1$ |

> The $\zeta(5)$ row is an outright **counterexample**, so the criterion is
> **[refuted]** and must not be quoted as a rule.  What the four rows *do*
> support is weaker and more specific:
>
> **Observation 2.1 (perfect-power degeneracy).**  When $\chi$ is a perfect
> power — $\bigl((x-8)(x-4)\bigr)^3$ for $\mathrm{Sym}^3\mathbf E$, $q_5^{\,7}$
> for $\zeta(7)$ level 24 — the row has **no $p$-adic limit at any prime**, and
> simultaneously its archimedean convergence degrades from geometric to
> **logarithmic**.  Both were measured directly.  Rows whose $\chi$ is not a
> perfect power (AESZ 207, the $\zeta(5)$ row) both have slopes.

**(d) The best available predictor is the geometric reading, and it is fitted,
not derived.**  $\sigma_p=-\max\{v_p(s):s\ \text{singular},\ s\ne0,\ |s|_p\ne1,\
s\ \text{not the $p$-adically distinguished one}\}$ reproduces both rows that
have slopes — AESZ 207 ($v_2(z)=-8,-12,-16$, giving $12$) and the $\zeta(5)$ row
($v_2(t)=-2,-1,-1,-1$, giving $1$) — but *which* singular point is excluded is
chosen differently in the two cases ($-8$ there, $-2$ here), so this is a fit,
not a rule.  It also overshoots on the degenerate rows.

> **Net.  At rank $\ge3$, $\sigma_p$ is measured, not derived.**  The only
> proved input that survives is Proposition 3.1.  This should be recorded in
> `CRYSTAL_THEOREM_F.md` §3, whose Cor. 3.2 is currently stated without the
> rank-two hypothesis being emphasised.

Two further facts recurred and are worth isolating.

* **The companion need not satisfy the $A$-row's own recurrence.**  For the
  $\zeta(5)$ row the residual $\sum_iQ_i(n)B_{n-i}$ is nonzero at **every**
  $n=16,\dots,399$ (384 of 384) **[exact]**; the minimal joint operator is
  $(r,D)=(18,11)$.  Same for $\zeta(7)$ level 24 (minimal $A$-operator order 22,
  joint operator order 47).  So "the $r-1$ companions of the $A$-recurrence" is
  not automatically the right object at higher rank — Lemma 1.1's picture
  applies to the *joint* operator.
* **Arithmetic normalisation can fail for the higher companions.**  For the
  $\zeta(5)$ row's joint operator, $X^{(1)},X^{(2)},X^{(3)}$ admit **no
  $k\le14$** with $\operatorname{lcm}(1..n)^kx_n\in\mathbf Z$ ($n\le150$)
  **[verified]**, because forward iteration divides by $Q_0(n)$, a degree-11
  polynomial with $\sim$220-digit coefficients whose primes are essentially
  random.  They are genuine solutions with genuine $2$-adic limits, but they are
  **not extension periods on the same footing as $B$**.  For AESZ 207 and for
  $\zeta(7)$ level 24 the normalisation *does* survive ($k=4$ resp. $7$, sharp,
  for every companion).

### 2.1 Row 1 — $\beta(4)$ at level 24 $=\mathrm{Sym}^3$ of Zagier $\mathbf E$

$A_n=[t^n]F(t)^3$ with $F=\sum e_nt^n$ the Zagier $\mathbf E$ row
$(12,4,32)$: $A_n=1,12,108,880,6876,52752,401744,\dots$

> **Correction (important).**  The characteristic roots are $\mathbf 4$ and
> $\mathbf 8$, **not** $\{512,256,128,64\}$.  $A_n$ is the *Cauchy* cube of $F$,
> not the $\mathrm{Sym}^3$ of the local monodromy: convolution does not move the
> singularities, it only raises the exponents there.  Concretely
> $$\chi_A(x)=(x-4)^3(x-8)^3,\qquad A_n\asymp 8^n\ \text{(not }512^n),$$
> and $|\lambda_2/\lambda_1|=\tfrac12$, matching `ROOT_ROWS.md` §5.4's
> "geometric, $\lambda_2/\lambda_1=1/2$".  **[verified]**

* **Minimal recurrence: shift order $6$, polynomial degree $4$** — $P_0=n^4$,
  $P_6=2^{15}(n-3)^4$, coefficients in `lattice/multislope/sc_rows.gp` (`R4cf`);
  verified on $n=6..148$ with $0$ failures and reproducing the series cube
  exactly to $n=120$.  Minimality proved by scan: **no** recurrence of order
  $\le5$ and degree $\le14$ exists ($52$ (order,degree) cells, $200$ exact
  terms, $0$ hits).  Since $\deg P_i=4$ the underlying **ODE has rank $4$** with
  $z$-degree $6$ — so there are **five** companions $X^{(1)},\dots,X^{(5)}$.
* $k=4$ sharp for every companion **[verified, $N=1500$]**.
* **Newton polygon at $p=2$: slopes $[2,3],[3,3]$** — lowest slope $2$ with
  multiplicity $\mathbf 3$; $\chi$ is a perfect cube, so Observation 2.1
  applies, and indeed:
* **$\sigma_p=0$ at $p=2,3,5,7,13$ for all five companions** (increments bounded
  and negative at $n=700..800$).  **$\beta(4)$ at level 24 has no $p$-adic Apéry
  limit at any small prime.**  **[verified]**
* Product formula: $\prod\lambda_i=4^38^3=2^{15}$, so
  $\log|\prod\lambda_i|=15\log2=10.397$, while
  $\sum_p\sigma_p\log p=0$.  **The $\Lambda$-ceiling of `THEORY_NOTES_03` §1 fails
  here by the whole $10.397$** — the $2\times2$ Casoratian of $(A,X^{(j)})$ is
  not the $6\times6$ Wronskian, and only the former controls $\sigma_p$.  (Same
  phenomenon as AESZ 184, `L3CHI5_TWO_WORLDS.md` §1.3, but total rather than
  partial.)
* Score and budget: $\log(1/|\lambda_2|)-k=\log\tfrac14-4=-5.386$; harvestable
  budget $=-5.386+0=-5.386$.  Not a candidate for anything.
* Dwork tower and good-prime towers: see §4.2, §4.4 — the row passes Dwork at
  $p=2$ only, and it is where the higher companions first behave *better* than
  the first (Conjecture MS-3).

### 2.2 Row 2 — AESZ 207 (rank 4, recurrence order 4, three companions)

**Confirmed against `AESZ207.md`.**  $A_1=-2320$, $A_2=57601296$,
$A_3=-2373661139200$; $B_2=-29761$, $B_3=\tfrac{32129592496}{27}$; $A_n$ integral
for all $n\le1600$; $\sigma_2=12$; and
$\xi_\infty(B)=-0.000504554593441367088625425457970715161172154460579143450288\ldots$
reproduced exactly.  **[verified / exact]**

* $\chi=(x-53248)^2(x^2+89344x-2^{24})$, roots $-89531.38920672014678\ldots$,
  $+187.38920672014678\ldots$, $53248$ (double), product $-2^{48}13^2$
  **[exact]**.
* **$k=4$ sharp for all three companions** ($k=3$ fails) **[verified $n\le200$]** —
  new for $C$ and $D$.
* **$\sigma_2=12$ for all three** (slopes $11.985/11.99/11.99$ over
  $n\in[1200,1600]$), and a scan over **every prime $p<60$** ($n\le600$) shows
  $p=2$ is the **only** prime with a slope, for every companion.  $v_2(A_n)=12$
  at every sampled $n$, so $\mu_1=0$.  **[verified $n\le1600$]**
* Newton polygon at $2$: $\{8,16,12,12\}$ (using $\sqrt{17}\in\mathbf Q_2$,
  $17\equiv1\bmod8$); at $13$: $\{1,1,0,0\}$; at $17$: flat.  Bad primes $2,13$
  (from $\prod=-2^{48}13^2$) and $17$ as the conifold field
  ($\operatorname{disc}(x^2+89344x-2^{24})=2^{16}5^217^3$).  **[exact]**
* $2$-adic limits, Cauchy precision $19156/19144/19129$ digits **[verified]**:
  $$v_2(\xi_2(B))=v_2(\xi_2(C))=-4,\qquad v_2(\xi_2(D))=-8,$$
  $16\xi_2(B)\equiv360776799631$, $16\xi_2(C)\equiv460323768881$,
  $256\xi_2(D)\equiv449896537247\ \ (\mathrm{mod}\ 2^{40})$; the $B$ line
  reproduces `AESZ207.md` §3 exactly.
* Archimedean limits at $N=5000$ (exact rationals, $1400$ digits working
  precision; certified by $|r_N-r_{N-100}|\approx10^{-1118}$, $10^{-1121}$,
  $10^{-1125}$) **[floating, 1118/1121/1125 digits]**:
  $$\xi_\infty(C)=+1.501998871618802820919660703567693867874\ldots,\qquad
  \xi_\infty(D)=-3.026195781429694336873872900571897312700\ldots$$

> **Result 2.2 (the headline of Program B).  AESZ 207 is a genuine
> *multi-period* row.**  The three companion limits are $\mathbf Q$-linearly
> independent together with $1$, **at both places**:
> * archimedean: `lindep` on $[1,\xi_B,\xi_C,\xi_D]$ at $1100$ digits returns
>   $220$–$224$-digit coefficients with residual $3.83\cdot10^{-671}$; a genuine
>   relation would leave $\approx10^{-1118}$, and $-671=-1118+2\cdot224$ is
>   exactly the spurious signature.  Every $3$-term subset gives $\approx290$-digit
>   coefficients at residual $10^{-596}$; every $2$-term subset $\approx545$ digits.
>   **No $\mathbf Q$-relation of height below $\approx10^{200}$.**
> * $2$-adically at $19000$ digits: $4$-term `lindep` returns $1431$-digit
>   coefficients against a spurious level of $10^{1430}$; pairs $\approx2860$
>   (level $10^{2860}$); triples $\approx1907$.  With $\sqrt{17}$ adjoined,
>   likewise.  **No $\mathbf Q$- and no $\mathbf Q(\sqrt{17})$-relation of height
>   below $\approx10^{1400}$.**
> * $\xi_\infty(C)/\xi_\infty(B)=-2976.88077988600456148\ldots$ and
>   $\xi_\infty(D)/\xi_\infty(B)=+5997.75687461135022676\ldots$ are neither
>   rational nor quadratic (`algdep` degree 2 returns $298$-digit heights at
>   $1100$ digits), hence not in $\mathbf Q(\sqrt{17})$.
>
> So the rank-four extension really carries **three independent periods at each
> place** — the answer to §1.3's period-rank question is
> $\operatorname{prk}_\infty=\operatorname{prk}_2=3$, the maximum.

This was **not** the expected answer.  The double root $53248$ spans a
two-dimensional $\mathbf Q$-subspace of the solution space whose members decay
relative to $A$ at both places, which would have forced two simultaneous
$\mathbf Q$-relations.  The computation refutes it, and the reason is instructive:
"$u_n=O(53248^n)$" is a *growth* condition, and the subspace of
$\mathbf Q$-rational solutions satisfying it is generically $0$, not
two-dimensional — the coefficient-of-$\lambda_1$ functionals are transcendental
Stokes constants whose kernels meet $\mathbf Q^4$ only in $0$.

Extending `AESZ207.md` §4.2 to the higher companions: $\xi_2(C)$ and $\xi_2(D)$
are not rational multiples of $\zeta_2(m)$ for $m=2,\dots,7$, not in
$\mathbf Q+\mathbf Q\zeta_2(5)$, not in
$\mathbf Q\zeta_2(3)+\mathbf Q\zeta_2(5)$, not rational, and `algdep` of degree 4
returns heights $\ge10^{58}$ — all at the spurious level.  No single-constant
identification for $\xi_\infty(C),\xi_\infty(D)$ against
$\{\zeta(3),\zeta(5),\zeta(7),\pi^2..\pi^6,\pi^2\zeta(3),G,\log2,\log^22,\log^32,
\log13,\log17,\sqrt{13},\sqrt{17},\zeta(3)^2,\log^42\}$ either.  (The
Kubota–Leopoldt code is Washington Thm 5.11 with $F=4$ and passed the exact
interpolation test $L_2(1-n,\mathbf1)=-(1-2^{n-1})B_n/n$ for $n=2,4,\dots,12$
before use.)

### 2.3 Row 3 — the $\zeta(5)$ host at level 16

**New exact structure.**  The *minimal* recurrence for $A_n$ alone is
**order 16, degree 11** (the repo had only the over-parametrised $(17,10)$
member); reconstructed exactly over $\mathbf Q$ by kernel mod a 2200-digit prime
plus rational reconstruction, and verified by exact substitution on all $A_n$,
$n=16,\dots,399$ **[exact]** (`row3_rec.txt`).  Its characteristic polynomial is
$(x+2)^5(x+4)^5(x^2-4x-28)^2q(x)$ with $q$ of enormous height coming from
apparent singularities; the genuine roots are $-2,-4,2\pm4\sqrt2$, confirming
`ZETA5_TWO_ROW.md` §1.2a.  $B_n$ does **not** satisfy this recurrence (residual
nonzero at every $n=16,\dots,399$); the minimal **joint** operator is
$(r,D)=(18,11)$, reconstructed exactly and verified on both $A$ and $B$ for
$n=18,\dots,399$ **[exact]** (`row3_rec18.txt`).

* $B$ has $k=5$ sharp **[verified $n\le150$]**, as recorded.
* $\sigma_2=1$ for $B$ and for the three joint-operator companions
  $X^{(1)},X^{(2)},X^{(3)}$ (measured $0.991,0.987,0.988,0.986$ at $N=2600$,
  working precision $2^{30000}$); no slope at $p=3,5,7$ for any of them.
  $v_2(\xi_2)=-8,-10,-8,-13$.  **[measured]**
* **Pipeline validation, and a sharpening of the repo's headline result:**
  `lindep([\xi_2(B),\zeta_2(5)])=[32,-7]`, i.e.
  $$\xi_2(B)=\tfrac7{32}\zeta_2(5)\quad\text{verified to }\mathbf{2576}\ 2\text{-adic digits}$$
  (`ZETA5_TWO_ROW.md` §0.2 had $372$).  **[verified]**
* The higher companions do **not** reduce to $\zeta_2(5)$: at $\approx1900$
  $2$-adic digits none is a rational multiple of $\zeta_2(m)$, $m=2..7$, none is
  in $\mathbf Q+\mathbf Q\zeta_2(5)$ or
  $\mathbf Q\zeta_2(3)+\mathbf Q\zeta_2(5)$, and
  $\{1,\xi_2(B),\xi_2(X^{(1)}),\xi_2(X^{(2)}),\xi_2(X^{(3)})\}$ has **no
  $\mathbf Q$-relation** (5-term `lindep` returns 117-digit coefficients against
  a spurious level $10^{114}$).  **[verified]**
* **Caveat that must travel with this row:** the $X^{(j)}$ admit no $k\le14$
  (§2.0), so they are solutions with $2$-adic limits but not arithmetically
  normalised extension periods.  Building a $d_n^k$-integral second companion
  for this host **[failed]** — the natural extension periods beyond $B$ are not
  produced by seeding the fitted operator, and no route to them was found.

Candidate (ii), the level-12 $\mathrm{Sym}^4$ system, needs no computation: its
roots $\alpha^{\pm4},\alpha^{\pm2},1$ have product $1$, so every root is a
$p$-adic unit at every $p$, the Newton polygon is flat, and no $p$-adic Apéry
limit is possible.  **[structural]**



### 2.4 Row 4 — $\zeta(7)$ at level 24 ($\mathrm{Sym}^6$)

**[verified]**  Built from `lattice/census/zeta7_level24.gp`'s construction.

* Minimal recurrence: shift order **22**, polynomial degree **23** (kernel
  dimension 1, verified exactly for $n=22,\dots,1597$, no failures); the
  $\theta$-form operator has $\theta$-order $7$ (the $\mathrm{Sym}^6$ rank) and
  $z$-degree $35$; the joint $(A,B)$ operator has order $47$, degree $8$.  All
  three reconstructed exactly over $\mathbf Q$ by multi-prime CRT plus rational
  reconstruction.
* $\chi_A=(x-1)^6(x^2-8x+8)^4(x^2+4x-4)^4$,
  $\chi_B=(x-1)^7(x^2-8x+8)^7(x^2+4x-4)^7$; both are powers of the base quintic
  $q_5=(x-1)(x^2-8x+8)(x^2+4x-4)$, roots $\{1,\,4\pm2\sqrt2,\,-2\pm2\sqrt2\}$.
  Hence
  $$\lambda_1=4+2\sqrt2=6.82842\ldots,\qquad\lambda_2=-2-2\sqrt2=-4.82842\ldots,
  \qquad\Bigl|\frac{\lambda_2}{\lambda_1}\Bigr|=\frac1{\sqrt2}\ \textbf{exactly},$$
  matching the measured archimedean decay $0.70847$ and confirming
  `SLOPE_CENSUS.md` §5's observed $1/\sqrt2$ — but **not** its explanation.
* $d_n^kB_n\in\mathbf Z$ with $k=7$ **sharp, verified to $n\le1597$** ($k=6$
  fails).  $B_n/A_n$ agrees with $\tfrac{1463}{13824}\zeta(7)$ to **105 digits**.
* Newton polygons `[slope, multiplicity]`:

| | $p=2$ | $p=3$ | $p=5$ | $p=7$ |
|---|---|---|---|---|
| $q_5$ | $[0,1],[1,2],[\tfrac32,2]$ | $[0,5]$ | $[0,5]$ | $[0,5]$ |
| $\chi_A$ (order 22) | $[0,6],[1,8],[\tfrac32,8]$ | $[0,22]$ | $[0,22]$ | $[0,22]$ |
| $\chi_B$ (order 35) | $[0,\mathbf 7],[1,14],[\tfrac32,14]$ | $[0,35]$ | $[0,35]$ | $[0,35]$ |

  At $p=3,5,7$ the polygon is **flat**: every characteristic root is a $p$-adic
  unit, there is no $p$-adically distinguished mode, and no limit is possible.
  At $p=2$ the polygon is not flat, but the lowest slope $0$ has multiplicity
  $6$ resp. $7$, and $\chi_B=q_5^{\,7}$ is a perfect power — Observation 2.1
  again, and it is sharper than the census's stated reason.
* The higher companions $X^{(1)},X^{(2)}$ of the order-35 operator also have
  $k=7$ sharp and still no slope at $p=2,3,5,7$.  (The order-47 joint companions
  need $k=41$; Row 3's companions have no $k\le14$ at all.)

> **Correction to `SLOPE_CENSUS.md` §5.**  That section says the row has no
> $p$-adic slope because "the product of characteristic roots appears to be
> coprime to $2,3,5,7$", and conjectures the ladder
> $\beta^{\pm6},\dots,\beta^{\pm2},1$ with $\beta=2^{1/4}$.  Both are wrong:
> $\prod(\text{roots})=-32$ per quintet ($2^{20}$ for $\chi_A$, $-2^{35}$ for
> $\chi_B$), so $v_2=5$ per quintet, **not** $0$; and the root set
> $\{1,4\pm2\sqrt2,-2\pm2\sqrt2\}$ is not inversion-stable, so it is not a
> $\mathrm{Sym}^6$ ladder — the only ratio ladder present is two steps of
> $\sqrt2$, which is what produces $|\lambda_2/\lambda_1|=1/\sqrt2$.  The
> conclusion (no slope at $2,3,5,7$) survives; the reason is the perfect-power
> degeneracy of Observation 2.1.

### 2.5 Summary table, and the three answers

| row | rank / rec. order | char. roots | bad primes | $k$ (sharp) | slope prime(s) | $\sigma_p$ | $\operatorname{prk}$ |
|---|---|---|---|---|---|---|---|
| $\beta(4)@24=\mathrm{Sym}^3\mathbf E$ | $4$ / $6$ (5 companions) | $8^{(3)},4^{(3)}$ | $2$ | $4$ | **none** | $0$ at $2,3,5,7,13$ | — |
| AESZ 207 | $4$ / $4$ (3 companions) | $-89531.4,\,187.4,\,53248^{(2)}$ | $2,13,17$ | $4$ | $\mathbf 2$ | $12$, all three | $\mathbf 3$ at both places |
| $\zeta(5)$ level 16 | — / $16$; joint $18$ | $-2,-4,2\pm4\sqrt2$ | $2$ | $5$ ($B$); none $\le14$ ($X^{(j)}$) | $\mathbf 2$ | $1$, all four | $\ge2$ ($2$-adically) |
| $\zeta(7)$ level 24 | $7$ / $22$; $\theta$-form $35$; joint $47$ | $1,4\pm2\sqrt2,-2\pm2\sqrt2$ (mult. 7) | $2$ | $7$ | **none** | $0$ at $2,3,5,7$ | — |

**(A) Do the higher companions carry new $p$-adic information?**  **Yes — this
is the strongest positive result of the programme.**  AESZ 207's three companion
limits admit no $\mathbf Q$-relation with $1$ of height below $\approx10^{1400}$
$2$-adically and $\approx10^{200}$ archimedeanly, and no
$\mathbf Q(\sqrt{17})$-relation either; the $\zeta(5)$ host's higher companions
are independent of $\xi_2(B)=\tfrac7{32}\zeta_2(5)$ and are not
Kubota–Leopoldt values.  The rank-four extension has **three** periods, not one.

**(B) Which primes give limits for which companion?**

* $\mathrm{Sym}^3\mathbf E$: none, any companion, $p=2,3,5,7$ **[verified $n\le800$]**.
* AESZ 207: $p=2$ only, and it gives a limit for **all three**, all with
  $\sigma_2=12$; nothing at any prime $p<60$ **[verified]**.
* $\zeta(5)$ level 16: $p=2$ only, for $B$ and all three $X^{(j)}$, all
  $\sigma_2=1$ **[measured to $N=2600$]**.
* $\zeta(7)$ level 24: none, either companion, $p=2,3,5,7$ **[verified $n\le260$]**.

**(C) Is there a row where two different companions have limits at two different
primes?**  **No — and the reason is systematic.**

> **Result 2.3.  The slope belongs to the row, not to the companion.**  Within a
> row, either every companion has a limit at $p$ with the *same* $\sigma_p$, or
> none does.  In AESZ 207 the three increments at $n=1600$ are
> $19156,19144,19129$ — the same slope to within $27$ parts in $19000$.  The
> multi-*slope* phenomenon this programme was designed to find **does not occur
> in any of the four rows.**

What AESZ 207 exhibits instead is the adjacent phenomenon, worth naming:

> **multi-period at a single slope** — several $\mathbf Q$-independent $p$-adic
> numbers arriving at the same prime at the same rate.

A genuine multi-slope row would need an operator whose Newton polygon at two
*different* primes each has an isolated distinguished slope, with the companion
basis adapted differently at the two primes.  None of these four is such, and
§3 explains why one should not expect the pay-off even if one were found.

---


---

## 3. The design-rule numbers, and what a second slope prime actually buys

*(script `lattice/multislope/mp_twoprime.gp`, log alongside.)*

### 3.1 At order two the multi-prime resource is already counted — and can only be lost

Measured at $N=600$ over the six Zagier rows **[verified]** (the measured
$\sigma_p$ agree with $v_p(c)$ to $O(\log N/N)$, as Corollary 3.2 of
`CRYSTAL_THEOREM_F.md` predicts):

| row | $c$ | $v_p(c)$, $p=2,3,5,7$ | measured $\sigma_p$ |
|---|---|---|---|
| $\mathbf A$ | $-8$ | $3,0,0,0$ | $2.960,\,0,\,0,\,0$ |
| $\mathbf B$ | $27$ | $0,3,0,0$ | $0,\,2.975,\,0,\,0$ |
| $\mathbf C$ | $9$ | $0,2,0,0$ | $0,\,1.975,\,0,\,0$ |
| $\mathbf D$ | $-1$ | $0,0,0,0$ | $0,0,0,0$ |
| $\mathbf E$ | $32$ | $5,0,0,0$ | $4.957,\,0,\,0,\,0$ |
| **$\mathbf F$** | $72=2^33^2$ | $\mathbf{3,2,0,0}$ | $\mathbf{2.960,\,1.975,\,0,\,0}$ |

> **Zagier $\mathbf F$ is the unique row in the second-order census with two
> slope primes.**  Cooper's $s_{18}$ has $c=192=2^63$ but is not in $(R2)$/$(R3)$
> normalisation and measures $\sigma_2=0$, $\sigma_3=1$
> (`SLOPE_CENSUS.md` §1, `CRYSTAL_THEOREM_F.md` §3); every other row has $c$ a
> prime power (or a unit).

**What that buys: nothing, and here is why.**  For an order-two row the product
formula $\Lambda\lambda=c$ gives
$\gamma\log\tfrac1\lambda+G=\gamma\log\Lambda_{\rm dec}$ *whatever the prime
factorisation of $c$*, provided the partner row has the same slope profile.  So
the design rule
$$\delta>1\iff\log\Lambda_{\rm dec}>k\max(r,1)+r\log\rho_2^{\rm eng}$$
is blind to how $c$ splits.  Taking $\mathbf F$ as both decayer and engine (two
samplings of one row, $r=\sigma_{\rm dec}/\sigma_{\rm eng}=1$):
$$\log9=2.1972\ \not>\ 2\cdot1+\log8=4.0794 ,$$
a deficit of $1.882$ per index.  Equivalently, $\mathbf F$'s budget
$\log\Lambda-k=\log9-2=+0.1972$ is positive but $\Lambda/\lambda=9/8$, so
$H-F=\gamma\log(\Lambda/\lambda)=0.1178$ is tiny and
$$\delta=1-\frac FH\le0.1112\quad(\text{at }\alpha=\gamma=1),$$
with $\delta=0.0732,0.0883,0.1112,0.0567,0.0380$ at $\alpha/\gamma=0.5,0.75,1,1.5,2$
**[computed]**.

> **Proposition 3.1 (the multi-prime rate is a loss mechanism, not a gain).**
> For order-two rows the total harvestable $p$-adic resource is
> $\sum_p\sigma_p\log p=\log|c|$, independent of the factorisation of $c$.  The
> multi-prime form $G=\sum_p\min(\sigma_p^{\rm eng}\alpha,\sigma_p^{\rm dec}\gamma)\log p$
> can only *reduce* it, and does so exactly when the two rows' slope supports
> differ.  In the corpus every pair sharing a period shares at most **one** slope
> prime: $\{\mathbf B,\mathbf C,\mathbf F,s_{18}\}$ all carry
> $L(2,\chi_{-3})$ and all align at $p=3$ only; $\mathbf F$'s extra $p=2$
> resource ($3\log2=2.079$ per index) is unharvestable because no other
> $L(2,\chi_{-3})$ row has a $2$-adic slope.

### 3.2 At higher rank the product formula itself fails

At order two the design rule is blind to the factorisation of $c$ because
$\Lambda\lambda=c$ turns the $p$-adic resource into the numerator growth
automatically.  At rank $\ge3$ **that identity is gone**, and with it the whole
mechanism:

| row | $\log\bigl|\prod\lambda_i\bigr|$ | $\sum_p\sigma_p\log p$ (measured) | gap |
|---|---|---|---|
| $\mathrm{Sym}^3\mathbf E$ | $15\log2=10.397$ | $0$ | $-10.397$ |
| AESZ 207 | $48\log2=33.271$ | $12\log2=8.318$ | $-24.953$ |
| $\zeta(7)$ level 24 | $35\log2=24.260$ | $0$ | $-24.260$ |
| $\zeta(5)$ level 16 | — (apparent singularities) | $\log2=0.693$ | — |

(The AESZ 184 anomaly of `L3CHI5_TWO_WORLDS.md` §1.3 — the $\Lambda$-ceiling
failing by $4\log2$ outside Zagier normalisation — is the mild case; here it
fails by everything.)  Consequently the correct row score, `05_two_row.tex`
eq. (orderm), must be used with a **measured** $\sum_p\sigma_p\log p$:

| row | $\log(1/|\lambda_2|)$ | $k$ | score | $+\sum_p\sigma_p\log p$ | **harvestable budget** |
|---|---|---|---|---|---|
| $\mathrm{Sym}^3\mathbf E$ | $-1.386$ | $4$ | $-5.386$ | $+0$ | $\mathbf{-5.386}$ |
| AESZ 207 | $-10.883$ | $4$ | $-14.883$ | $+8.318$ | $\mathbf{-6.565}$ |
| $\zeta(7)$ level 24 | $-1.575$ | $7$ | $-8.575$ | $+0$ | $\mathbf{-8.575}$ |
| $\zeta(5)$ level 16 | $-1.386$ | $5$ | $-6.386$ | $+0.693$ | $\mathbf{-5.693}$ |

Every one is deeply negative, and none is a decayer ($|\lambda_2|>1$ in all four),
so none can even play the decayer role in a two-row construction.  Using
$\log\Lambda-k$ instead would report $+7.4$ for AESZ 207 — the exact error
`05_two_row.tex` correction (i) warns against, here worth $22$ per index.

> **Result 3.2.  There is no best pair.**  The multi-slope design question is
> answered in the negative twice over: (i) no row in the corpus has two
> harvestable slope primes shared with a partner carrying the same period
> (§3.1), and (ii) the four higher-rank rows have no $p$-adic resource worth
> harvesting at all — three have $\sigma_p=0$ everywhere, and AESZ 207's
> $\sigma_2=12$ is dwarfed by its $k=4$ and its $|\lambda_2|=53248$.  The
> quantity that would have to improve is $\log(1/|\lambda_2|)$, i.e. the
> *archimedean* separation, and at rank $4$ it is worse than at rank $2$, not
> better.

**What the higher rank does buy** is not Diophantine but structural: three
independent periods at one prime (§2.5(A)), and companions that can be strictly
better behaved than the first one (Conjecture MS-3).  Those are the results to
keep.

---

## 4. Supercongruences for the companions

### 4.1 What is being tested, and why the normalisation matters

Three statements are separated here, and they are logically independent.

* **(D-pure) Dwork/Gauss tower for the $A$-row.**
  $a_{mp^s}\equiv a_{mp^{s-1}}\pmod{p^s}$.  A theorem for every row with a
  Laurent-polynomial model, at every $p$ including $p=2$ (Mellit–Vlasenko;
  `CRYSTAL_THEOREM_F.md` §1.3, §2.1).  For AESZ 207 and the $\mathrm{Sym}^3$ row
  no Laurent model is on record, so it is a genuine test.
* **(S) the slope law**, $v_p(x_n-\xi_pa_n)=n\sigma_p+O(\log n)$ — this is
  Proposition 3.1 of `CRYSTAL_THEOREM_F.md` and is *proved*; the question is only
  the size of the $O(\log n)$ term.
* **(SC) a genuine supercongruence**: an *extra* $p$-power along the tower
  $n=mp^s$, i.e.
  $$E(m,s):=v_p\bigl(x_{mp^s}-\xi_p\,a_{mp^s}\bigr)-\sigma_p\,m\,p^{s}$$
  unbounded and growing in $s$.  A bounded $E$ means "no supercongruence, just
  the slope law"; $E\sim cs$ means the Dwork tower is *sharpening* the generic
  bound by a factor $p^{cs}$.

The precision budget is the trap: $\xi_p$ is only known to the precision of the
largest $n$ used to compute it, so $E(m,s)$ is meaningful only while
$\sigma_pmp^s+E\ll\sigma_pN_{\max}$.  Every table below states its budget.

### 4.2 T1 — the Dwork tower for the $A$-rows

`sc_t1.gp`, `sc_t1b.gp`, `sc_t1c.gp`; $N=3000$, $p=2,3,5,7,13$, $5985$ tests per
row.  Rows: R1 Apéry $\zeta(3)$ (control), R2 Zagier $\mathbf C$, R3 AZ $\eta$,
R4 $\mathrm{Sym}^3$ Zagier $\mathbf E$, R5 AESZ 207.

| row | failures of $a_{mp^s}\equiv a_{mp^{s-1}}\ (p^s)$ | sharpened statement verified |
|---|---|---|
| R1 (control) | $0/5985$ | $p^{3s}$ at $p=5,7,13$; $p^{3s-1}$ at $p=2,3$ — **reproduces Beukers–Coster**, validating the pipeline |
| R2 | $0/5985$ | $p^{2s}$, all five primes |
| R3 | $0/5985$ | $p^{3s}$ at $p=3,5,7,13$; $2^{3s-2}$ |
| **R4** $\mathrm{Sym}^3\mathbf E$ | **$2043/5985$** — $0$ at $p=2$, $889/560/361/233$ at $p=3,5,7,13$ | $2^{2s}$ at $p=2$, **sharp** ($2^{2s+1}$ fails at $m=1,s=2$) |
| **R5** AESZ 207 | $0/5985$ | $p^{2s}$ at odd $p$; $\mathbf{2^{3s}}$ at $p=2$ |

Two results here.

> **Result 4.1 (new).  AESZ 207 is a Dwork row, with a supercongruence.**
> $$A_{mp^s}\equiv A_{mp^{s-1}}\pmod{p^{2s}}\ (p\ \text{odd}),\qquad
> A_{m2^s}\equiv A_{m2^{s-1}}\pmod{2^{3s}} .$$
> **[verified: $N=3000$, $p=2,3,5,7,13$, $5985+5985$ tests, $0$ failures;
> sharp — $p^{2s+1}$ fails at odd $p$, and $v_2(A_{m2^s}-A_{m2^{s-1}})=3s+6$
> exactly for $s\ge3$.]**  `AESZ207.md` never tested this.  The bare Dwork half
> is what a Laurent-polynomial model would give (Mellit–Vlasenko), and no such
> model is known for AESZ 207 — **this is evidence that one exists**, which
> in turn is hypothesis (F1) of Theorem F$'$ (`CRYSTAL_THEOREM_F.md` §4).

> **Result 4.2.  A Cauchy cube inherits the Dwork tower at $p=2$ only.**
> $\mathrm{Sym}^3$ Zagier $\mathbf E$ passes at $p=2$ ($2993$ tests, $0$
> failures, and in fact mod $2^{2s}$, sharp) and **fails at every odd prime,
> already at the mod-$p$ (Lucas) level**: first failure $p=3$, $m=1$, $s=1$,
> $a_3-a_1=880-12=868$, $v_3=0$.  **[verified]**  So "$A_n=[t^n]F(t)^3$ with $F$
> a Dwork row" does not inherit Dwork's congruences — the Cauchy power is not a
> constant-term construction.

### 4.3 T2 — no supercongruence for the companions at slope primes

Slope primes in the row set (re-measured, `sc_diag_run.gp`): R2 at $p=3$
($\sigma_3=2$), R3 at $p=5$ ($\sigma_5=3$), R5 at $p=2$ ($\sigma_2=12$ for
**all three** companions, secant slopes $12.000/12.002/12.000$ over
$n=500..1000$).  **R4 has $\sigma_p=0$ at every $p\in\{2,3,5,7,13\}$ for all five
companions** (increments bounded and negative at $n=700..800$): it has no
$p$-adic Apéry limit anywhere, so T2 is empty for it.

Precision budget: Cauchy precision $3986$ ($R2$, $p=3$), $5974$ ($R3$, $p=5$),
$71951$ ($R5$, $p=2$, $N=6000$); every $E$-value is read at $n\le N/2$, so the
margin is at least half the precision — no quantity used comes within an order
of magnitude of the cut.

**The answer is negative and clean.**  With
$E(m,s)=v_p\bigl(x_{mp^s}-\xi_pa_{mp^s}\bigr)-\sigma_p\,mp^s$:

* R2 at $p=3$: $E\in\{-4,\dots,-1\}$, $|E|\le8$ over $182$ cells ($m\le300$,
  $s\ge2$, $N=4000$);
* R3 at $p=5$: $E\in\{-5,\dots,-1\}$, $|E|\le5$ over $35$ cells;
* R5 at $p=2$: $E$ **exactly constant**.

$E$ is bounded, **not** growing in $s$: *there is no extra $p^s$, and the
hoped-for supercongruence does not exist.*  For R2/R3 the mechanism is explicit
and was verified as an identity,
$$a_nx_{n+1}-a_{n+1}x_n=\frac{9^n}{(n+1)^2}\ (\text{R2}),\qquad
=\frac{125^n}{(n+1)^3}\ (\text{R3})$$
**[verified, 4000 checks each, 0 failures]**, whence
$\xi-t_n=\sum_{j\ge n}c^j/\bigl((j+1)^ka_ja_{j+1}\bigr)$ and $E$ is a bounded
expression in $v_p(a_j)$; the first-term-dominant prediction
$E(m,s)=-k\,v_p(mp^s+1)-v_p(a_{mp^s+1})$ holds in $490/630$ (R2) and $395/399$
(R3) cells, the exceptions being exactly where two tail terms have equal
valuation.

> **Conjecture MS-1 (exact $2$-adic linear-form law for AESZ 207).**
> With $\sigma_2=12$ and $\xi_2^{(j)}=\lim_nx^{(j)}_n/A_n$,
> $$v_2\bigl(x^{(j)}_n-\xi_2^{(j)}A_n\bigr)=12n+c_j\ \textbf{exactly, for every even }n,
> \qquad c_1=-4,\quad c_3=-28,$$
> while for $j=2$,
> $v_2\bigl(x^{(2)}_n-\xi_2^{(2)}A_n\bigr)=12n+\max\bigl(-12+v_2(n),\,-8\bigr)$
> for $v_2(n)\ne4$, with interference exactly at $v_2(n)=4$.
> **[verified: $N=6000$, all $n\le3000$; $1501$ even-$n$ cells per companion,
> $0$ exceptions; plus $744$ tower cells ($m\le250$, $s\ge2$) with $0$
> deviations; precision budget $71951$ against a maximum quantity $\approx36000$.]**

The $j=2$ shape is a textbook **slope-separation signature**: $X^{(2)}$ is a sum
of two Frobenius components, of valuations $12n-12+v_2(n)$ and $12n-8$, which
cross precisely at $v_2(n)=4$ — the only stratum where the value is not
constant.  Stratified by $v_2(n)$ (all $n\le3000$):

| $v_2(n)$ | $\#n$ | $X^{(1)}$ | $X^{(2)}$ | $X^{(3)}$ |
|---|---|---|---|---|
| $0$ | 1499 | $[-3,7]$ | $\mathbf{-12}$ | $[-27,-16]$ |
| $1$ | 750 | $\mathbf{-4}$ | $\mathbf{-11}$ | $\mathbf{-28}$ |
| $2$ | 375 | $\mathbf{-4}$ | $\mathbf{-10}$ | $\mathbf{-28}$ |
| $3$ | 188 | $\mathbf{-4}$ | $\mathbf{-9}$ | $\mathbf{-28}$ |
| $4$ | 94 | $\mathbf{-4}$ | $[-7,2]$ | $\mathbf{-28}$ |
| $5..11$ | 47,23,12,6,3,1,1 | $\mathbf{-4}$ | $\mathbf{-8}$ | $\mathbf{-28}$ |

(bold = a single value with **zero** exceptions).  The Dwork-type congruence for
the companion, T2(iii), follows: e.g. R5 $X^{(1)}$ has
$v_2(t_{m2^{s+1}}-t_{m2^s})=12\cdot m2^s-8$ for $m=1$ and $36\cdot2^s-12$ for
$m=3$, i.e. $\sigma_pmp^s+c_j-v_p(a_{mp^s})$ throughout.

### 4.4 T3 — good-prime towers for the higher companions

Control (R1 $\zeta(3)$, $N=8192$, $p=2,3,5,7,13$, $a=1,2,3$, 15 cells):
$v_p(t_{s+1}/t_s)\equiv-3$ in every cell, agreement depths grow by exactly $3$
per level, limit exactly $+1\cdot p^{-3}$, i.e. $\varepsilon=+1$, $w+1=3=k$
**[verified]** — the pipeline reproduces the trichotomy of
`GOOD_PRIME_TOWERS.md` §6.

**AESZ 207** ($N=6000$–$8192$):

| $p$ | $X^{(1)}$ | $X^{(2)}$ | $X^{(3)}$ | verdict |
|---|---|---|---|---|
| $3$ | $-3,-4,-3,-4,\dots$ | $-2,-5,-1,-6,-2,-5,-2$ | $1,-8,-2,-5,-2,-5,-1$ | all three have mean $-\tfrac72$: **no limit** |
| $5$ | $\to-3$, units agree only mod $5$ | $\to-3$, same | alternates $-4,-2$ | **no limit** |
| $7$ | $\to-4$ | $\to-4$ | $\to-4$; units agree only to $7^2$–$7^3$ | **no limit** |
| $11$ | valuation not constant | same | same | **no limit** |
| $\mathbf{13}$ | $-2$ | $-2$ | $-2$ | **limit $\Lambda=13^{-2}$**, $\varepsilon=+1$ |

> **Result 4.3.**  `AESZ207.md` §5's negative at $p=3$ (alternating $-3,-4$,
> mean $-\tfrac72$, non-integral) is **not a property of the first companion**:
> all three companions have consecutive increments summing to $-7$, so the
> half-integral $3$-adic slope belongs to the whole extension.  **[verified]**
> At $p=13$ a tower limit *does* exist — $\varepsilon(13)=+1$, $w+1=2$ (not
> $k=4$), uniform over all bases ($\min_av_{13}(13^2t_{13n}/t_n-1)=0,2,5$ at
> $s=0,1,2$, over $630+48+3$ cells) — but $13$ is a **bad** prime for the
> operator ($53248=2^{12}\!\cdot\!13$, $2^{28}\!\cdot\!13z^3$,
> $2^{44}\!\cdot\!13^2z^4$), so it does not rescue the character-determination
> programme.  **[numerical, $s\le3$]**

**$\mathrm{Sym}^3$ Zagier $\mathbf E$** ($N=8192$, $15625$; five companions):

> **Conjecture MS-2 ($2$-adic tower of $\mathrm{Sym}^3\mathbf E$).**  For every
> companion $X^{(j)}$, $j=1..5$, and every $n\ge1$,
> $$\frac{x^{(j)}_{2n}}{A_{2n}}\equiv-\tfrac14\,\frac{x^{(j)}_n}{A_n}
> \pmod{2^{\varepsilon(s)}},\qquad\varepsilon(s)\to\infty,$$
> so $\Lambda^{(j)}_a=\lim_s(-4)^st_{a2^s}$ exists with
> $\varepsilon(2)=-1$, $w+1=2$.
> **[verified: $N=8192$, all bases $a=1..4096$, levels $s=0..6$; $8128$ cells per
> companion $\times5$ companions $=40640$ cells;
> $\min_av_2(-4t_{2n}/t_n-1)=1,2,6,7,9,11,13$, monotone.]**

> **Conjecture MS-3 ($3$-adic tower of the *higher* companions only).**  For
> $j=3,4,5$, $\;x^{(j)}_{3n}/A_{3n}\equiv-3^{-4}x^{(j)}_n/A_n\pmod{3^{2s-1}}$,
> $n=a3^s$: $\varepsilon(3)=-1$, $w+1=4$.  For $j=1,2$ no limit exists
> ($v_3(t_{s+1}/t_s)$ alternates $-4,0$).
> **[verified: $N=15625$, all bases, $s=1..6$;
> $\min_av_3(-3^4t_{3n}/t_n-1)=1,3,5,7,9,11$; $4091$ cells per companion.]**

Two things in MS-2/MS-3 are worth flagging.

1. **$w+1\ne k$.**  `CRYSTAL_THEOREM_F.md` §2.2 identifies the tower exponent
   $w+1$ with the denominator exponent $k$ (the "Hodge depth").  Here $k=4$ but
   the $2$-adic tower exponent is $2$, and the $3$-adic one (for $j\ge3$) is $4$.
   **The identification $w+1=k$ fails on this row**, and the two primes even
   disagree with each other.  This is a genuine counterexample to the §2.2
   reading and should be recorded there.
2. **The higher companions can be strictly better behaved than the first.**  At
   $p=3$, $X^{(1)},X^{(2)}$ have no tower limit while $X^{(3)},X^{(4)},X^{(5)}$
   do.  This is the first instance in the project where a *later* companion
   carries information the first one does not — the one genuinely new
   phenomenon that higher rank supplies.
3. Suggestive but unsettled: limits appear at $p=3,7,11$ ($\chi_{-4}(p)=-1$) and
   not at $p=5,13,17$ ($\chi_{-4}(p)=+1$), which would match
   $\beta(4)=L(4,\chi_{-4})$.  At $p=5$ this is **[verified]** to $s=6$
   ($N=15625$, agreement depth stuck at $\le2$, $u\equiv16\bmod25$); at
   $p\ge13$ only $2$–$3$ levels are reachable, so **[inconclusive]**.



---

## 5. Status ladder and open problems

| statement | status |
|---|---|
| Lemma 1.1: $L(y^{(j)})=j^mz^j$, so a recurrence of order $r$ carries $r-1$ extension classes | **[proved]** |
| Proposition 3.1 of `CRYSTAL_THEOREM_F.md` holds at rank 4 | **[verified]**, AESZ 207 residual $\pm27$ over $n\le1600$ |
| $\sigma_p=\mu_{\det}-2\mu_1$ and $\sigma_p=v_p(c)$ extend to rank $\ge3$ | **[refuted]**, three counterexamples (§2.0(b)) |
| "lowest Newton slope simple $\iff$ limit exists" | **[refuted]** by the $\zeta(5)$ level-16 row (§2.0(c)) |
| Perfect-power $\chi$ $\Rightarrow$ no $p$-adic limit at any $p$, and logarithmic archimedean convergence | **[verified]** on two rows |
| AESZ 207 is a Dwork row, with $A_{mp^s}\equiv A_{mp^{s-1}}$ mod $p^{2s}$ ($p$ odd), mod $2^{3s}$ | **[verified]**, $5985{+}5985$ tests, 0 failures |
| A Cauchy cube inherits Dwork at $p=2$ only | **[verified]**, $\mathrm{Sym}^3\mathbf E$ fails mod $3$ already |
| AESZ 207 has three $\mathbf Q$-independent periods at each place | **[verified]** to heights $10^{200}$ (arch.), $10^{1400}$ ($2$-adic) |
| No supercongruence for companions at slope primes: $E(m,s)$ bounded | **[verified]** on R2, R3, R5 |
| Conjecture MS-1 (exact $2$-adic law for AESZ 207) | **[conjecture, verified over $n\le3000$]** |
| Conjecture MS-2, MS-3 (towers of $\mathrm{Sym}^3\mathbf E$) | **[conjecture, 40640 resp. 4091 cells]** |
| $w+1=k$ (tower exponent $=$ denominator exponent, `CRYSTAL_THEOREM_F.md` §2.2) | **[refuted]** on $\mathrm{Sym}^3\mathbf E$: $k=4$, $w+1=2$ at $p=2$ and $4$ at $p=3$ |
| The slope belongs to the row, not the companion | **[verified]** on all four rows |
| A genuine multi-slope row exists in the corpus | **[negative]** |

**Corrections this note forces on earlier documents.**

1. `SLOPE_CENSUS.md` §5: the $\zeta(7)$ level-24 root set and the claim
   "$\prod\lambda_i$ coprime to $2$" are wrong; the correct roots are
   $\{1,4\pm2\sqrt2,-2\pm2\sqrt2\}$ each with multiplicity $7$, and
   $v_2(\prod)=35$.  The conclusion (no slope) survives; the mechanism is the
   perfect-power degeneracy.  (§2.4)
2. `CRYSTAL_THEOREM_F.md` §3: Corollary 3.2 must carry its $(R2)$/$(R3)$
   hypothesis explicitly; the naive higher-rank analogue is false on three rows.
   §2.2's identification of the tower exponent with $k$ is false on
   $\mathrm{Sym}^3\mathbf E$.
3. `ZETA5_TWO_ROW.md` §0.2: $\xi_2=\tfrac7{32}\zeta_2(5)$ is now verified to
   $2576$ $2$-adic digits (was $372$), and the minimal recurrence for $A_n$ is
   order 16 degree 11 (the repo carried the over-parametrised $(17,10)$).
4. `AESZ207.md` §5: the non-integral $3$-adic tower slope $-\tfrac72$ is a
   property of the whole extension — all three companions have consecutive
   increments summing to $-7$ — not of the first companion; and a tower limit
   *does* exist at the bad prime $p=13$ ($\varepsilon=+1$, $w+1=2$).

**Open problems, in priority order.**

1. **A formula for $\sigma_p$ at rank $\ge3$.**  Prop. 3.1 gives existence and
   rate in terms of the $2\times2$ Casoratian; what is missing is an expression
   for that minor's $p$-adic growth in terms of the operator.  The two rows with
   slopes both match a "second distinguished singularity" reading, but with
   different choices of which singularity is distinguished (§2.0(d)).
2. **Identify AESZ 207's three periods.**  With $1118$ archimedean and $19129$
   $2$-adic digits and every standard basis excluded, the next step is the
   connection problem at the conifold $z_-$ (`AESZ207.md` §6.2(1)), now with
   three constants to express rather than one.
3. **Prove Conjecture MS-1**, the exact law $v_2(x^{(j)}_n-\xi_2^{(j)}A_n)=12n+c_j$.
   The $j=2$ interference at $v_2(n)=4$ is a clean slope-separation signature and
   should be provable from the Frobenius decomposition.
4. **Explain Conjecture MS-3**: why do $X^{(3)},X^{(4)},X^{(5)}$ of
   $\mathrm{Sym}^3\mathbf E$ have a $3$-adic tower limit when $X^{(1)},X^{(2)}$
   do not?  This is the only place in the project where a later companion
   carries information the first does not.
5. **Settle the $\chi_{-4}$ pattern** for $\mathrm{Sym}^3\mathbf E$: limits at
   $p=3,7,11$ ($\chi_{-4}(p)=-1$) and none at $p=5$ ($\chi_{-4}(p)=+1$,
   verified to $s=6$); $p\ge13$ needs a larger $N$.
6. **Row 1's archimedean limits.**  Convergence is $O(1/\log n)$, so
   $n=3\cdot10^5\to10^6$ buys $0.1$ decimal digit; only $2$–$3$ digits are
   available and no identification is possible without a different method
   (connection coefficients at $t=1/8$).



---

## 6. Scripts

All in `lattice/multislope/`, each with a `.log` alongside.

**Rank-4 / higher-rank census.**

| file | what it does |
|---|---|
| `row1_fit.gp`, `row1_fit2.gp` | $\mathrm{Sym}^3\mathbf E$: kernel-dimension scan mod $2^{61}-1$ over $(r,D)$; exact reconstruction of the minimal $(6,4)$ recurrence over $\mathbf Q$; verification $n=6..420$; char poly.  Writes `row1_rec.txt` |
| `row1_main.gp` | builds $A$ and all five companions exactly; sharp $k$; $v_p(A_n)$; $\sigma_p$ at $p=2,3,5,7$; archimedean ratios.  Writes `row1_data.txt` |
| `row1_asym.gp`, `row1_logfit.gp`, `row1_logfit2.gp` | float runs to $N=10^6$; establishes $A_n\sim C\,8^n(\log n)^2/n$; $1/\log n$ extrapolation of the five limits |
| `row1_arch.gp` | **superseded** — a $1/n$ Richardson fit, invalid here; its `lindep` output must not be quoted |
| `row1_v2.gp` | $v_2(A_n)$ boundedness; Newton polygons |
| `row2_padic.gp` | AESZ 207: operator, char poly, Newton polygons at $2,13,17$, $\sqrt{17}\in\mathbf Q_2$; builds $A,B,C,D$ to $N=1600$; sharp $k$; $\sigma_p$; $v_2(\xi_2)$.  Writes `row2_data.txt` |
| `row2_arch.gp`, `row2_resid.gp` | $N=5000$ exact archimedean run (1118–1125 digits); `lindep` on all subsets of $[1,\xi_B,\xi_C,\xi_D]$; the explicit residual $3.83\cdot10^{-671}$ proving the negative.  Writes `row2_arch.txt` |
| `row2_lindep.gp`, `row2_kl.gp` | $2$-adic `lindep` at 19000 digits, with $\sqrt{17}$; Kubota–Leopoldt $\zeta_2$ (Washington Thm 5.11, $F=4$) with its exact interpolation validation |
| `row2_primes.gp` | slope scan over every prime $p<60$, all three companions |
| `row2_archid.gp` | archimedean identification sweep for $\xi_\infty(C),\xi_\infty(D)$ |
| `row3_fit.gp`, `row3_joint.gp`, `row3_comp.gp` | $\zeta(5)$ level 16: exact minimal $(16,11)$ recurrence for $A$; the check that $B$ does not satisfy it; the minimal joint $(18,11)$ operator.  Write `row3_rec.txt`, `row3_rec18.txt` |
| `row3_slopes.gp`, `row3_padic.gp` | the companions $X^{(1..3)}$, the failed $k\le14$ test, $\sigma_2=1$, and $\xi_2(B)=\tfrac7{32}\zeta_2(5)$ to 2576 digits.  Writes `row3_data.txt` |
| `row4_build.gp`, `row4_fit*.gp`, `row4_exA/exB/exC.gp`, `row4_char.gp`, `row4_comp.gp`, `row4_C2.gp`, `row4_bchk.gp`, `row4_resid.gp` | $\zeta(7)$ level 24: data, minimal $(22,23)$ / $\theta$-form $(35,7)$ / joint $(47,8)$ operators, char polys and Newton polygons, companions and their $k$ and slopes |

**Design rule and prime support.**

| file | what it does |
|---|---|
| `mp_primesupport.gp` | slope-prime support of every second- and third-order census row; establishes Zagier $\mathbf F$ as the unique row with two slope primes |
| `mp_twoprime.gp` | the two-sampling design numbers for Zagier $\mathbf F$ ($F$, $H$, $\delta$ at five ratios) |

**Dwork congruences and supercongruences.**

| file | what it does |
|---|---|
| `sc_rows.gp` | row library (R1–R5), including the fitted $\mathrm{Sym}^3\mathbf E$ recurrence `R4cf` |
| `sc_kexp.gp` | sharp $k$ for every companion of every row |
| `sc_t1.gp`, `sc_t1b.gp`, `sc_t1c.gp` | T1: Dwork tower and its sharpenings, $N=3000$, $p=2,3,5,7,13$ |
| `sc_diag_run.gp`, `sc_slope4.gp` | slope diagnostics; the $\sigma_p=0$ verdict for $\mathrm{Sym}^3\mathbf E$ |
| `sc_t2*.gp` | T2: $E(m,s)$ tables, the Casoratian identities, the stratification by $v_2(n)$ |
| `sc_t3*.gp` | T3: good-prime towers, uniform-over-all-bases tests, the $p=13$ positive for AESZ 207 and the $p=2,3$ positives for $\mathrm{Sym}^3\mathbf E$ |

