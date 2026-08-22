# The cusp move as a group action: orbits, companions, and what several
# realisations of one class at different cusps do and do not buy

*Claude (Opus 5), 2026-08-22. Scripts and logs: `lattice/cusp_move/`.
Executes the programme "several realisations of one extension class at different
cusps: same $p$-adic limits, different $\lambda_2$". Builds on
`ACF_ONE_SURFACE.md` (Theorem 1, the $\mathbf A/\mathbf C/\mathbf F$ move),
`HERFURTNER_CLASSIFICATION.md` §3 and Theorem H5 (the exponent normal form),
`ROOT_ROWS.md` (Theorems R1–R4), `SQRT_APERY.md`, `NONCONGRUENCE_SCAN.md` §5
(the $\Gamma_0(5)+5$ row), `paper/sections/03_archimedean.tex` (score, measures)
and `paper/sections/05_two_row.tex` (Theorem E, the design rule).*

---

## 0. Verdict first

Four theorems, one census of $49$ rows, and one sharp negative answer.

* **Theorem 1 (the general cusp move).** *Proved symbolically.* For **every**
  second-order Apéry-like row — not just the Zagier normalisation
  $Q(n)=dn^2$ — the move
  $$s=\frac{t}{1-\lambda t},\qquad z(s)=(1-\lambda t)^{\alpha}y(t)$$
  produces again a three-term row of degree $\le 2$ **if and only if**
  $$\boxed{\ (\alpha-1+r_1)(\alpha-1+r_2)=0,\qquad\text{i.e.}\qquad
  \alpha\in\{1-r_1,\,1-r_2\},\ }$$
  where $r_1,r_2$ are the roots of $Q$, so that $1-r_i$ are exactly the two local
  exponents of $L$ at $t=\infty$. The transformed row is given in closed form
  (§2), and the operators are related by the clean identity
  $L^{\sharp}z=(1-\lambda t)^{1-\alpha}\,Ly$. `01_general_move.py`.

* **Theorem 2 (the companion, and how $\xi$ moves).** *Proved, and verified
  exactly in all $22$ available (row, $\lambda$, $\alpha$) cases to $n=40$.*
  The naive guess $B^{\sharp}=(1-\lambda t)^{\alpha}B$ is **false**. The correct
  statement is
  $$\boxed{\ B^{\sharp}=(1-\lambda t)^{\alpha}\widehat B,\qquad
  L\widehat B=t\,(1-\lambda t)^{-\alpha}=\Bigl(\tfrac{\partial t}{\partial t^{\sharp}}\Bigr)^{\!-1}\!\!\cdot\text{(gauge)},\ }$$
  i.e. **the inhomogeneity $t$ is replaced by the new Hauptmodul** rather than
  the companion being gauged. In modular terms this is
  $B^{\sharp}=F^{\sharp}D^{-2}\bigl(\Phi/(1-\lambda t)\bigr)$, which is exactly
  the mechanism behind the $\Xi$-correction of `ACF_ONE_SURFACE.md` (5.1);
  the $-\tfrac54$ there is the price of that correction, not a gauge factor.
  `00_selftest.gp`, `11_companion_general.gp`.

* **Theorem 3 (orbit structure).** The three finite non-MUM singular points of a
  row are $v\in\{0,\lambda,\mu\}$ in the coordinate $v=1/t$, and the cusp move is
  translation on that line. Hence there are exactly **three placements**
  (which of the three sits at $t=\infty$) and, on each, up to **four gauges**,
  giving an orbit of at most $12$ rows. Measured orbit sizes in the corpus:
  $2,3,8,12$. All $55$ listed orbit members ($49$ distinct rows: the
  $\mathbf A$, $\mathbf C$ and $\mathbf F$ orbits coincide) are integral to
  $n=200$ after an explicit rescaling $c\in\{1,2,4,9\}$, and all have $k\in\{1,2\}$, sharp.
  `02_orbits.gp`, `03_census.gp`.

* **Theorem 4 (rigidity — the programme's fatal obstruction).** *Proved.*
  Let a row be integral with $a=\operatorname{lc}P$, $d=\operatorname{lc}Q\in\mathbf Z$,
  $a^2\ne4d$, $d\ne0$. Then $\lambda,\mu$ are algebraic integers, so
  $$\lambda\in\mathbf Q\iff \lambda,\mu\in\mathbf Z\iff a^2-4d\ \text{is a square},$$
  and in that case each of the three gaps $|\lambda|,|\mu|,|\lambda-\mu|$ is a
  **positive integer**. Since $|\lambda_2|$ at a placement is the smaller of the
  two gaps at that vertex,
  $$\boxed{\ a^2-4d\ \text{a square}\ \Longrightarrow\ |\lambda_2|\ge1\ \text{at every placement}
  \ \Longrightarrow\ \operatorname{score}\le-k\le-1 .\ }$$
  **No rational cusp-move orbit can contain a positive-score row.** Conversely a
  positive-score row ($|\lambda_2|<e^{-k}<1$) *necessarily* has irrational
  characteristic roots, so its only other equally-good placement is its **Galois
  conjugate**, defined over the real quadratic field
  $K=\mathbf Q(\sqrt{a^2-4d})$ and not over $\mathbf Q$. Verified on the whole
  corpus (`05_placements.gp`): the best score in a rational orbit is exactly
  $-2$, attained by $\mathbf A$, $\mathbf C$, $\sqrt{s_7}$ and their partners.

* **The programme's headline hypothesis is half true and half false, and the
  split is the interesting part.**
  * *False archimedeanly.* Different placements almost never carry the same
    period. Across the $8$ rational orbits ($49$ distinct rows) the **only** rationally
    proportional pair is $(\mathbf C,\mathbf F')$ with ratio $-\tfrac54$
    (verified to $80$ decimal digits); every other pair has no rational relation
    of height $\le10^{12}$ (§5). Even for the two positive-score rows, the
    Galois-conjugate placement approximates a **different** real number
    ($0.3278170944\ldots$ against $\zeta(2)/5$; $0.1002150533\ldots$ against
    Beukers' $0.1001874492\ldots$), §7. `04_periods.gp`, `09_conj_xi.gp`.
  * *True $p$-adically, and non-trivially so.* Inside the $\sqrt{\text{Domb}}$
    orbit, placements $2$ and $3$ satisfy $\xi_p^{(2)}=-2\,\xi_p^{(3)}$
    **simultaneously at $p=2$ and $p=3$** (to $792$ and $386$ digits), and
    $\xi_p^{(4)}=\xi_p^{(6)}$ at both primes ($789$, $390$ digits), while their
    archimedean periods are $\mathbf Q$-unrelated. This is the phrase
    "same $p$-adic limits, different $\lambda_2$" realised exactly.
    `06_padic.gp`.

* **Two genuinely new identified rows.** The $\sqrt{\text{Domb}}$ orbit contains
  $$(20,10,64):\ (n{+}1)^2u_{n+1}=(20n^2{+}26n{+}10)u_n-64n^2u_{n-1},\qquad
    \xi=\tfrac1{16}\bigl(15L(2,\chi_{-3})-6\zeta(2)\bigr),$$
  $$(20,4,64):\ (n{+}1)^2u_{n+1}=(20n^2{+}14n{+}4)u_n-64n^2u_{n-1},\qquad
    \xi=\tfrac1{32}\bigl(6\zeta(2)+15L(2,\chi_{-3})\bigr),$$
  both integral to $n=200$, $k=2$ sharp, $\lambda_1=16$, $\lambda_2=4$,
  agreement with the target to the full truncation precision
  ($2.8\cdot10^{-124}$ resp. $2.7\cdot10^{-118}$ at $n=200$, predicted
  truncation $3.9\cdot10^{-121}$). Their parent $\sqrt{\text{Domb}}$ carries the
  **cuspidal** value $L(f_{12},2)$: the cusp move has pushed the weight-three
  source out of the cuspidal part into the Eisenstein part, extending the
  $\mathcal S$-plane/$\mathcal T$-plane observation of `ACF_ONE_SURFACE.md` §4
  from two Eisenstein families to (cusp form $\to$ Eisenstein). These are the
  **first rows in the project whose single Apéry limit is a non-trivial
  $\mathbf Q$-combination of $\zeta(2)$ and $L(2,\chi_{-3})$** — exactly the
  shape a Nesterenko argument for $\dim\langle1,\zeta(2),L(2,\chi_{-3})\rangle$
  wants. `10_domb_ident.gp`.

* **Irrationality measures (Task 2): no improvement, quantified.** The merged
  two-sequence bound is $\mu\le1+\min_i\sigma_i/\delta_i$ — the better of the two
  placements, with no convex-combination gain at leading order. Numbers:
  Apéry's $\zeta(2)$ row gives $\mu\le11.85078219105$ and its conjugate
  placement $11.87072446756$ (worse); Beukers gives $50.65365917218$ and its
  conjugate $50.64743124451$ (better by $0.00622793$) — but the conjugate rows
  have coefficients in $\mathcal O_K$, the norm rate is
  $\log|\lambda_2\bar\lambda_2|+2k=4$ resp. $6.7717$, i.e. $>0$, so no descent to
  $\mathbf Q$ is possible, and in any case they approximate a different number.
  For the one same-period pair $(\mathbf C,\mathbf F')$ both $\delta<0$
  ($-2$ and $-4.0794$), so the bound is vacuous. §7. `07_measures.gp`.

* **Simultaneous approximation and Nesterenko (Task 3): a deficit of
  $10.394$ nats per step.** From the pair $(\mathbf A,\mathbf C)$ on the common
  surface $I_6I_3I_2I_1$ one gets integer linear forms in $1,\zeta(2),L(2,\chi_{-3})$
  with $\sigma=k+\log\lambda_1=4.19722$ and $\log|L_n|/n=+2$, so
  $$\tau=\delta/\sigma=-0.4765053580,\qquad
  \dim_{\mathbf Q}\langle1,\zeta(2),L(2,\chi_{-3})\rangle\ \ge\ 1+\tau=0.5234946420,$$
  vacuous. Deficit to $\dim\ge2$: $2$ nats/step. Deficit to $\dim\ge3$ (the CDT
  theorem): $\mathbf{10.394449}$ nats/step. Even the *budget* version
  (all $p$-adic resource harvested, $|\lambda_2|\to1/\lambda_1$) reaches only
  $\tau=0.018927$, $\dim\ge1.018927$. §8. `07_measures.gp`, `10_domb_ident.gp`.

---

## 1. Notation

A **row** is
$$(n+1)^2u_{n+1}=P(n)\,u_n-Q(n)\,u_{n-1},\qquad u_0=1,\ u_{-1}=0,\qquad
\deg P,\deg Q\le2,$$
so $u_1=P(0)$. Write $P(n)=p_2n^2+p_1n+p_0$, $Q(n)=q_2n^2+q_1n+q_0$,
$a=p_2$, $d=q_2$. The Picard–Fuchs operator is
$$L=\theta^2-t\,P(\theta)+t^2Q(\theta+1),\qquad \theta=t\frac{d}{dt},$$
equivalently $Ly=R_2t^2y''+(R_2+R_1)ty'+R_0y$ with
$R_2=1-p_2t+q_2t^2$, and $A=\sum a_nt^n$ solves $LA=0$, the companion
$B=\sum b_nt^n$ ($b_0=0$, $b_1=1$) solves $LB=t$.

**Characteristic roots** $\lambda,\mu$: the roots of $x^2-ax+d$, so
$R_2=(1-\lambda t)(1-\mu t)$ and the finite singular points are
$t_1=1/\lambda$, $t_2=1/\mu$. **$Q$-roots** $r_1,r_2$: $Q(n)=d(n-r_1)(n-r_2)$.
By `HERFURTNER_CLASSIFICATION.md` Theorem H5(i) the exponents of $L$ are
$(0,0)$ at $t=0$, $(0,\rho_1)$ and $(0,\rho_2)$ at $t_1,t_2$, and
$(1-r_1,1-r_2)$ at $t=\infty$, with $\rho_1+\rho_2=r_1+r_2$ (Fuchs).

$\xi=\lim b_n/a_n$ (archimedean), $\xi_p$ the same limit in $\mathbf Q_p$,
$d_n=\operatorname{lcm}(1,\dots,n)$, $k=\min\{j:d_n^jb_n\in\mathbf Z\}$,
$\operatorname{score}=\log(1/|\lambda_2|)-k$ with
$|\lambda_1|\ge|\lambda_2|$ (`paper/sections/03_archimedean.tex`).

Rows are considered modulo the **scaling** $u_n\mapsto c^nu_n$, i.e.
$(P,Q)\mapsto(cP,c^2Q)$; the tuple
$\bigl(a^2/d,\;p_1/a,\;p_0/a,\;q_1/d,\;q_0/d\bigr)$ is a complete invariant and
is used as the orbit key throughout.

---

## 2. Theorem 1: the cusp move for an arbitrary row

> **Theorem 1.** Let $y$ solve $Ly=0$ for a row $(P,Q)$ with characteristic roots
> $\lambda,\mu$ and $Q$-roots $r_1,r_2$. Put
> $$s=\frac{t}{1-\lambda t}\ \Bigl(\text{so } t=\frac{s}{1+\lambda s}\Bigr),\qquad
> z(s)=(1-\lambda t)^{\alpha}y(t).$$
> Then $z$ satisfies a second-order equation whose $s$-coefficients are
> polynomials of degree $\le2$ — i.e. $z$ is again the analytic solution of a
> three-term Apéry-like row — **iff**
> $$(\alpha-1+r_1)(\alpha-1+r_2)=0 .$$
> With $\alpha=1-r_1$ the new row is
> $$\boxed{\begin{aligned}
> P^{\sharp}(n)&=(\mu-2\lambda)\,n^2+\bigl(p_1-3\lambda+2\lambda r_1\bigr)n
>   +\bigl(p_0-\lambda+\lambda r_1\bigr),\\[1mm]
> Q^{\sharp}(n)&=\lambda(\lambda-\mu)\,n^2
>   +\lambda\bigl(\lambda+\mu-p_1-2\lambda r_1+\mu r_1-\mu r_2\bigr)n\\
> &\qquad+\lambda\bigl(\lambda r_1^2-\lambda r_1+\mu r_1r_2-\mu r_1+p_1r_1\bigr),
> \end{aligned}}$$
> its characteristic roots are $-\lambda$ and $\mu-\lambda$
> ($R_2^{\sharp}=(1+\lambda s)(1+(\lambda-\mu)s)$), its $Q$-roots are
> $r_1^{\sharp}=1-\alpha=r_1$ and $r_2^{\sharp}=1-\alpha-\rho_1$, and the two
> operators are related by
> $$\boxed{\ L^{\sharp}z=(1+\lambda s)^{\alpha-1}\,Ly=(1-\lambda t)^{1-\alpha}\,Ly.\ }$$

*Proof.* Direct computation in `01_general_move.py` (sympy): substituting
$t=s/(1+\lambda s)$ and $y=(1+\lambda s)^{\alpha}z$ into $Ly$ and dividing by
$\rho=(1+\lambda s)^{\alpha-1}$ gives $R_2^{\sharp},R_1^{\sharp}$ polynomial of
degree $2$ unconditionally, and $R_0^{\sharp}$ polynomial with remainder
$$\frac{\alpha^2\lambda\mu-2\alpha\lambda\mu-\alpha q_1+\lambda\mu+q_0+q_1}{\lambda^2}
=\frac{d}{\lambda^2}\,(\alpha-1+r_1)(\alpha-1+r_2)$$
after substituting $q_1=-d(r_1+r_2)$, $q_0=dr_1r_2$. Setting this to zero is the
stated condition; the boxed coefficients are then read off. $\square$

**Reading.** The condition $\alpha\in\{1-r_1,1-r_2\}$ is exactly "gauge away the
exponent of $L$ at $t=\infty$": the move sends $t=\infty$ to the finite point
$s=-1/\lambda$, where an Apéry-type row needs local exponents $(0,\rho^{\sharp})$;
the gauge $(1-\lambda t)^{\alpha}$ shifts the exponents there by $-\alpha$, and
the two admissible choices are the two exponents $1-r_i$. This is the exponent
bookkeeping of Theorem H5 turned into a group action.

**Specialisations.**

| class | $r_1,r_2$ | $\alpha$ | move |
|---|---|---|---|
| Zagier, $Q=dn^2$ | $0,0$ | $1$ | $(a,b,d)\mapsto(\mu-2\lambda,\;b-\lambda,\;\lambda^2-\lambda\mu)$ |
| root rows, $Q=c(2n-1)^2$ | $\tfrac12,\tfrac12$ | $\tfrac12$ | $P^{\sharp}=(\mu-2\lambda)n^2+(p_1-2\lambda)n+(p_0-\tfrac{\lambda}2)$ |
| $\sqrt{s_7}$, $Q=-3(3n-1)(3n-2)$ | $\tfrac23,\tfrac13$ | $\tfrac13$ or $\tfrac23$ | two gauges |
| $\sqrt{s_{10}},\sqrt{s_{18}}$ | $\tfrac58,\tfrac38$ | $\tfrac38$ or $\tfrac58$ | two gauges |
| $\Gamma_0(5)+5$, Herfurtner #45 | $\tfrac34,\tfrac14$ | $\tfrac14$ or $\tfrac34$ | two gauges |
| Herfurtner #30, $Q=441(3n-1)^2$ | $\tfrac13,\tfrac13$ | $\tfrac23$ | one gauge |

The first line **is** Theorem 1 of `ACF_ONE_SURFACE.md`, recovered as the
$r_1=r_2=0$ case; the remaining lines are new — they are exactly the rows that
`HERFURTNER_CLASSIFICATION.md` §6.2 flagged as invisible to Zagier's
normalisation.

### 2.1 The coefficient form

Since $A^{\sharp}(s)=(1+\lambda s)^{-\alpha}A\bigl(s/(1+\lambda s)\bigr)$,
$$\boxed{\ a_n^{\sharp}=\sum_{m=0}^{n}\binom{n+\alpha-1}{\,n-m\,}(-\lambda)^{\,n-m}\,a_m\ }$$
which for $\alpha=1$ is the **signed binomial transform**
$a_n^{\sharp}=\sum_m\binom nm(-\lambda)^{n-m}a_m$. Consequences:

* for $\alpha=1$ and $\lambda\in\mathbf Z$, integrality of $a^{\sharp}$ is
  *automatic* — the Zagier-class orbit needs no rescaling at all (confirmed:
  scale $c=1$ throughout the $\mathbf A/\mathbf C/\mathbf F$ and $\mathbf E$
  orbits);
* for $\alpha\notin\mathbf Z$ the binomials carry denominators supported on the
  primes dividing $\operatorname{den}(\alpha)$, and the rescaling of Theorem R1
  of `ROOT_ROWS.md` reappears: measured $c\in\{1,2,4,9\}$ (§4).

The closed formula of Theorem 1 was cross-checked against this transform for
every corpus row, both roots and both gauges, exactly to $n=40$
(`02_orbits.gp`, `11_companion_general.gp`: 22/22 **VERIFIED**).

---

## 3. Theorem 2: the companion, and the exact transformation of $\xi$

> **Theorem 2.** In the situation of Theorem 1, let $\widehat B$ be the solution of
> $$L\widehat B=t\,(1-\lambda t)^{-\alpha}
>  =\sum_{k\ge1}\binom{\alpha+k-2}{k-1}\lambda^{k-1}t^{k},\qquad \widehat b_0=0 .$$
> Then the companion of the moved row is
> $$\boxed{\ B^{\sharp}=(1-\lambda t)^{\alpha}\,\widehat B,\qquad\text{i.e.}\qquad
> b_n^{\sharp}=\sum_{m=0}^n\binom{n+\alpha-1}{n-m}(-\lambda)^{n-m}\,\widehat b_m .\ }$$

*Proof.* By Theorem 1, $L^{\sharp}\bigl[(1-\lambda t)^{\alpha}X\bigr]
=(1-\lambda t)^{1-\alpha}LX$ for any $X$. Taking $X=\widehat B$ gives
$L^{\sharp}B^{\sharp}=(1-\lambda t)^{1-\alpha}\cdot t(1-\lambda t)^{-\alpha}
=t/(1-\lambda t)=s=t^{\sharp}$, which is the defining inhomogeneity of the moved
row's companion; and $B^{\sharp}=s+O(s^2)$, matching $b^{\sharp}_0=0$,
$b^{\sharp}_1=1$. Uniqueness of the analytic solution with $b_0=0$ finishes. $\square$

**Verification.** `11_companion_general.gp` recomputes $b_n^{\sharp}$ from the
formula and compares with the recurrence companion of the moved row: **OK in all
22 cases to $n=40$**, and the naive $B^{\sharp}=(1-\lambda t)^{\alpha}B$ **fails
in all 22**, as predicted.

**Why $\xi$ is not preserved.** Write $\widetilde X=\widehat B-B$, so
$L\widetilde X=t^{\sharp}-t=\lambda t^2/(1-\lambda t)$. Both $A$ and $B$ are
gauged the same way, so
$$B^{\sharp}-\xi^{\sharp}A^{\sharp}=(1-\lambda t)^{\alpha}\bigl(\widehat B-\xi^{\sharp}A\bigr),$$
and $\xi^{\sharp}$ is fixed by regularity of $\widehat B-\xi^{\sharp}A$ at the
*new* dominant singularity. Two things change at once: the dominant singularity
(the move sends one of $t_1,t_2,\infty$ to $\infty$) **and** the second solution
($\widehat B\ne B$). If only the first changed, moving the *subdominant* root
would preserve $\xi$ exactly; the measured failure ($\mathbf C\to\mathbf A$ takes
$\tfrac12L(2,\chi_{-3})$ to $\tfrac14\zeta(2)$) is precisely the size of
$\widetilde X$.

**Modular reading.** With $\Phi=F\,Dt$ (weight 3), $\Theta=D^{-2}\Phi$, $B=F\Theta$,
$$B^{\sharp}=F^{\sharp}\,D^{-2}\!\left(\frac{\Phi}{1-\lambda t}\right),\qquad
B^{\sharp}-(1-\lambda t)^{\alpha}B=F^{\sharp}D^{-2}\bigl[\Phi^{\sharp}-\Phi\bigr].$$
The correction is the Eichler integral of the difference of the two weight-three
sources. In `ACF_ONE_SURFACE.md` §5 this correction is the series $\Xi$ and the
relation $\Theta_{\mathbf F'}=-\tfrac54\Theta_{\mathbf C}+\tfrac94\Theta_{\Psi_0}$;
Theorem 2 says that identity is not a coincidence of level 6 but the general
shape of the companion transform. Multiplication by the modular unit
$(1-\lambda t)^{-1}$ is not Hecke-equivariant, which is why the source can jump
between Hecke eigenspaces — between the two Eisenstein families at level 6
(`ACF_ONE_SURFACE.md` §4) and, as §5.3 below shows, between the cuspidal and the
Eisenstein part at level 12.

---

## 4. Theorem 3: orbits, and the census

In the coordinate $v=1/t$ the four singular points are $v=\infty$ (MUM),
$v=0$ (the point at $t=\infty$), $v=\lambda$ and $v=\mu$; the cusp move by
$\lambda$ is the translation $v\mapsto v-\lambda$. Hence:

> **Theorem 3.** The cusp-move orbit of a row consists of the **three placements**
> obtained by translating one of $\{0,\lambda,\mu\}$ to $v=0$, each carrying at
> most **four gauges** (the exponent at the image of $t=\infty$ can be normalised
> in two ways at each of the two finite singular points). The orbit therefore has
> at most $12$ elements modulo the scaling $u_n\mapsto c^nu_n$, and the placement
> at $v=v_j$ has $|\lambda_2|=\min_{i\ne j}|v_i-v_j|$.

The three placements of the Zagier normalisation are the three fibres of
`HERFURTNER_CLASSIFICATION.md` §3.1 that can be sent to $\infty$; Theorem 3 makes
"which fibre at $\infty$" into a $\mathbf Z$-translation and adds the gauge
direction, which §3.1 does not see.

### 4.1 The census (rational-root rows, `03_census.gp`)

Eight orbit computations, $55$ listed members ($49$ distinct rows --- the
$\mathbf A$, $\mathbf C$, $\mathbf F$ orbits coincide), **all integral to $n=200$**, all $k$ sharp
(smallest $j$ with $d_n^jb_n\in\mathbf Z$, checked $n\le200$).
$c$ is the minimal rescaling $u_n\mapsto c^nu_n$ needed for integrality.

| row | $c$ | $(a,b,d)$ | $P(n)$ | $Q(n)$ | $\lambda_1$ | $\lambda_2$ | $k$ | score |
|---|---|---|---|---|---|---|---|---|
| A.1 $=\mathbf A$ | 1 | $(7,2,-8)$ | $7n^2+7n+2$ | $-8n^2$ | $8$ | $-1$ | 2 | $-2.00000$ |
| A.2 $=\mathbf F'$ | 1 | $(-17,-6,72)$ | $-17n^2-17n-6$ | $72n^2$ | $-9$ | $-8$ | 2 | $-4.07944$ |
| A.3 $=\mathbf C$ | 1 | $(10,3,9)$ | $10n^2+10n+3$ | $9n^2$ | $9$ | $1$ | 2 | $-2.00000$ |
| E.1 $=\mathbf E$ | 1 | $(12,4,32)$ | $12n^2+12n+4$ | $32n^2$ | $8$ | $4$ | 2 | $-3.38629$ |
| E.2 | 1 | $(0,0,-16)$ | $0$ | $-16n^2$ | $4$ | $-4$ | 2 | $-3.38629$ |
| D.1 $=\sqrt{\text{Domb}}$ | 1 | $(20,2,64)$ | $20n^2+10n+2$ | $64n^2-64n+16$ | $16$ | $4$ | 2 | $-3.38629$ |
| D.2 | 1 | $(-28,-6,192)$ | $-28n^2-22n-6$ | $192n^2-96n$ | $-16$ | $-12$ | 2 | $-4.48491$ |
| D.3 | 1 | $(8,0,-48)$ | $8n^2+2n$ | $-48n^2+24n$ | $12$ | $-4$ | 2 | $-3.38629$ |
| D.4 | 1 | $(8,6,-48)$ | $8n^2+14n+6$ | $-48n^2-24n$ | $12$ | $-4$ | 2 | $-3.38629$ |
| **D.5** | 1 | $(20,10,64)$ | $20n^2+26n+10$ | $64n^2$ | $16$ | $4$ | 2 | $-3.38629$ |
| D.6 | 1 | $(-28,-12,192)$ | $-28n^2-34n-12$ | $192n^2+96n$ | $-16$ | $-12$ | 2 | $-4.48491$ |
| **D.7** | 1 | $(20,4,64)$ | $20n^2+14n+4$ | $64n^2$ | $16$ | $4$ | 2 | $-3.38629$ |
| D.8 | 1 | $(20,12,64)$ | $20n^2+30n+12$ | $64n^2+64n+16$ | $16$ | $4$ | 2 | $-3.38629$ |
| $s_7$.1 $=\sqrt{s_7}$ | 1 | $(26,2,-27)$ | $26n^2+13n+2$ | $-27n^2+27n-6$ | $27$ | $-1$ | 2 | $-2.00000$ |
| $s_7$.6 | 4 | $(104,62,-432)$ | $104n^2+160n+62$ | $-432n^2+12$ | $108$ | $-4$ | **1** | $-2.38629$ |
| $s_7$.10 | 4 | $(104,6,-432)$ | $104n^2+48n+6$ | $-432n^2+12$ | $108$ | $-4$ | 2 | $-3.38629$ |
| $s_7$.12 | 1 | $(26,15,-27)$ | $26n^2+39n+15$ | $-27n^2-27n-6$ | $27$ | $-1$ | 2 | $-2.00000$ |
| $s_{18}$.3 | 1 | $(-40,-14,256)$ | $-40n^2-44n-14$ | $256n^2-64n-12$ | $-32$ | $-8$ | **1** | $-3.07944$ |

(Full table, 55 listed members: `out/03_census.log`; the remaining rows of the
$\sqrt{s_7},\sqrt{s_{10}},\sqrt{s_{18}}$ orbits all have $k=2$ and score
between $-5.47$ and $-4.08$.) Orbit sizes: $\mathbf A/\mathbf C/\mathbf F$: $3$;
$\mathbf E$: $2$; $\sqrt{\text{Domb}}$: $8$; $\sqrt{s_7},\sqrt{s_{10}},\sqrt{s_{18}}$: $12$ each.

**Two members have $k=1$** — a *second* free integration on top of the one of
`ROOT_ROWS.md` Theorem R3. It is not free: both need a rescaling ($c=4$ for
$\sqrt{s_7}.6$), and $\log c=1.386$ buys at most one integer drop of $k$, so the
score changes by $(k_{\rm old}-k_{\rm new})-\log c=1-1.386=-0.386$. The
bookkeeping is exactly Theorem R1/R4 of `ROOT_ROWS.md`.

### 4.2 Theorem 4 and the gap picture (`05_placements.gp`)

For each corpus row the three gaps $|\lambda|,|\mu|,|\lambda-\mu|$ and the best
score over the orbit ($k=2$ throughout):

| row | $a^2-4d$ | roots | gaps | best $|\lambda_2|$ | best score | rational orbit? |
|---|---|---|---|---|---|---|
| Zagier $\mathbf A$ | $81=9^2$ | $8,-1$ | $8,1,9$ | $1$ | $-2.0000$ | **yes**, size 3 |
| Zagier $\mathbf C$ | $64=8^2$ | $9,1$ | $9,1,8$ | $1$ | $-2.0000$ | **yes**, size 3 |
| Zagier $\mathbf F$ | $1=1^2$ | $9,8$ | $9,8,1$ | $1$ | $-2.0000$ | **yes**, size 3 |
| Zagier $\mathbf E$ | $16=4^2$ | $8,4$ | $8,4,4$ | $4$ | $-3.3863$ | **yes**, size 2 |
| Zagier $\mathbf B$ | $-27$ | complex | — | — | none | no arch. limit |
| **Zagier $\mathbf D$** (Apéry $\zeta(2)$) | $125$ | $\tfrac{11\pm5\sqrt5}2$ | $11.0902,\,0.09017,\,11.1803$ | $0.09017$ | $\mathbf{+0.40606}$ | **no**, $K=\mathbf Q(\sqrt5)$ |
| **Beukers $=\sqrt{\text{Apéry}}$** | $18432$ | $4(17\pm12\sqrt2)$ | $135.882,\,0.11775,\,135.765$ | $0.11775$ | $\mathbf{+0.13920}$ | **no**, $K=\mathbf Q(\sqrt2)$ |
| $\sqrt T$ | $512$ | $12\pm8\sqrt2$ | $23.314,\,0.68629,\,22.627$ | $0.68629$ | $-1.6235$ | no, $K=\mathbf Q(\sqrt2)$ |
| $\sqrt{\text{Domb}}$ | $144=12^2$ | $16,4$ | $16,4,12$ | $4$ | $-3.3863$ | **yes**, size 8 |
| $\sqrt{\mathrm{AZ}(9,3,-27)}$ | $6912$ | $36\pm24\sqrt3$ | $77.569,\,5.5692,\,83.138$ | $5.5692$ | $-3.7173$ | no, $K=\mathbf Q(\sqrt3)$ |
| $\sqrt{\mathrm{AZ}(11,5,125)}$ | $-256$ | complex | — | — | none | no arch. limit |
| $\sqrt{\mathrm{AZ}(7,3,81)}$ | $-2048$ | complex | — | — | none | no arch. limit |
| $\sqrt{s_7}$ | $784=28^2$ | $27,-1$ | $27,1,28$ | $1$ | $-2.0000$ | **yes**, size 12 |
| $\sqrt{s_{10}}$ | $1600=40^2$ | $32,-8$ | $32,8,40$ | $8$ | $-4.0794$ | **yes**, size 12 |
| $\sqrt{s_{18}}$ | $64=8^2$ | $32,24$ | $32,24,8$ | $8$ | $-4.0794$ | **yes**, size 12 |
| $\Gamma_0(5)+5$ non-congruence | $8000$ | $44\pm20\sqrt5$ | $88.721,\,0.72136,\,89.443$ | $0.72136$ | $-1.6734$ | no, $K=\mathbf Q(\sqrt5)$ |
| Herfurtner #30 $I_1I_7II\,II$ | $-2187$ | complex | — | — | none | no arch. limit |
| Herfurtner #45 $I_3III\,III\,III$ | $-1728$ | complex | — | — | none | no arch. limit |

Every rational-root row has all three gaps integral and $\ge1$, as Theorem 4
demands; the best score in a rational orbit is exactly $-2$. **The two
positive-score rows of the whole project both sit in the "no rational orbit"
column, and this is forced, not accidental.**

### 4.3 The placements over $K$ (`08_quadratic.gp`)

For the five real-quadratic rows both non-trivial placements were computed
exactly over $K$, the sequences checked to be **algebraic-integer sequences with
$c=1$** (traces and norms in $\mathbf Z$, $n\le60$), and the scores recomputed:

| row | placement at $v=\lambda$ | placement at $v=\mu$ |
|---|---|---|
| Zagier $\mathbf D$ | $\lambda_1=-11.1803$, $\lambda_2=-11.0902$, score $-4.40606$ | $\lambda_1=11.1803=5\sqrt5$, $\lambda_2=0.090170$, score $\mathbf{+0.40606}$ |
| Beukers | $\lambda_1=-135.882$, $\lambda_2=-135.765$, score $-6.91092$ | $\lambda_1=135.765$, $\lambda_2=-0.117749$, score $\mathbf{+0.13920}$ |
| $\sqrt T$ | score $-5.11916$ | score $-1.62355$ |
| $\sqrt{\mathrm{AZ}(9,3,-27)}$ | score $-6.35117$ | score $-3.71725$ |
| $\Gamma_0(5)+5$ | score $-6.48550$ (both gauges) | score $-1.67338$ (both gauges) |

So the "second realisation with the same quality" **does exist** for both
positive-score rows — with the *same* score to the last digit — but only over
$K$, and §7 shows it approximates a different number.

---

## 5. The archimedean periods across an orbit

`04_periods.gp` computes $\xi$ for every orbit member to $\ge45$ decimals
(the truncation index $n$ is chosen so that $|\lambda_2/\lambda_1|^n<10^{-80}$;
$n$ up to $1564$), then tests (i) rational proportionality inside the orbit at
height $\le10^{12}$ and (ii) `lindep` against $\zeta(2),L(2,\chi_{-3}),G,\zeta(3)$.

### 5.1 The $\mathbf A/\mathbf C/\mathbf F$ orbit

| member | $(a,b,d)$ | $\xi$ | identification |
|---|---|---|---|
| $\mathbf A$ | $(7,2,-8)$ | $0.411233516712056609118103791662$ | $\zeta(2)/4$ (`lindep` $[-4,1,0,0,0]$) |
| $\mathbf C$ | $(10,3,9)$ | $0.390651206448243148433593714812$ | $\tfrac12L(2,\chi_{-3})$ ($[-2,0,1,0,0]$) |
| $\mathbf F'$ | $(-17,-6,72)$ | $-0.488314008060303935541992143515$ | $-\tfrac58L(2,\chi_{-3})$ ($[8,0,5,0,0]$) |

Ratios: $\xi^{\mathbf F'}/\xi^{\mathbf C}=-\tfrac54$ with error
$9.3\cdot10^{-81}$ (**exact**); $\xi^{\mathbf A}/\xi^{\mathbf C}$ has no rational
approximation of height $\le10^{12}$ better than $7.7\cdot10^{-24}$, i.e. the
random level. **Two of the three placements share a period, the third does not** —
the recurrence-level form of `ACF_ONE_SURFACE.md` Theorem 3(a).

### 5.2 The other orbits

Across the $\sqrt{\text{Domb}}$ ($8$), $\sqrt{s_7}$, $\sqrt{s_{10}}$,
$\sqrt{s_{18}}$ ($12$ each) orbits — $44$ further rows, $226$ pairs — there is
**no** rational proportionality: every best rational of height $\le10^{12}$
misses by $10^{-28}$ to $10^{-23}$, exactly the accuracy a random real would
give. So: **within a cusp-move orbit the archimedean period is generically not
preserved, even up to $\mathbf Q^\times$.**

### 5.3 Two new identified rows in the $\sqrt{\text{Domb}}$ orbit

`lindep` returned exact relations for two members
(`04_periods.log`, confirmed at $200$ digits in `10_domb_ident.gp`):

$$\boxed{\ \sqrt{\text{Domb}}.5=(20,10,64):\quad
\xi=\frac{15\,L(2,\chi_{-3})-6\,\zeta(2)}{16}=0.1156207370223709896358\ldots\ }$$
$$\boxed{\ \sqrt{\text{Domb}}.7=(20,4,64):\quad
\xi=\frac{6\,\zeta(2)+15\,L(2,\chi_{-3})}{32}=0.6746606435792704084951\ldots\ }$$

*Status:* **numerical**, at the full truncation precision — agreement
$2.8\cdot10^{-124}$ resp. $2.7\cdot10^{-118}$ at $n=200$ against a predicted
truncation error $(\lambda_2/\lambda_1)^{200}=3.9\cdot10^{-121}$. Both rows are
integral to $n=200$ with $k=2$ sharp.

This is structurally new. The parent $\sqrt{\text{Domb}}$ has
$\xi=L(f_{12},2)$, the critical value of the weight-three **newform** of level
$12$ (`ROOT_ROWS.md` §5.1). Its cusp-move relatives $.5$ and $.7$ have
Eisenstein periods. Since Theorem 2 identifies the change of source as
multiplication by the modular unit $(1-\lambda t)^{-1}$ followed by
$D^{-2}$, and multiplication by a unit is not Hecke-equivariant, the cusp move
can and does mix the cuspidal and Eisenstein parts of $M_3(\Gamma_0(12),\chi_{-3})$.
`ACF_ONE_SURFACE.md` §4 saw the phenomenon inside the Eisenstein space (the
$\mathcal S$-plane versus the $\mathcal T$-plane); this is the same mechanism one
step further.

*Consequence used in §8:* rows $.5$ and $.7$ give **single** rows whose linear
form involves $1,\zeta(2)$ **and** $L(2,\chi_{-3})$ simultaneously, with
coefficient vectors $(-6,15)$ and $(6,15)$ in the $(\zeta(2),L)$-plane —
linearly independent. This is the shape a Nesterenko argument needs, and to our
knowledge the project had no such row before.

---

## 6. The $p$-adic side: where the programme is right

`06_padic.gp`, $N=400$, measuring the Cauchy slope $v_p(\xi_N-\xi_{N-1})$ and
then testing $\xi_p^{(i)}/\xi_p^{(j)}$ against rationals of height $\le10^6$
with $O(p^{40})$ precision. A relation is called *real* only when the residual
valuation is at the row's own Cauchy precision (hundreds of digits), not at the
$\sim40$ digits a spurious `bestappr` produces.

| orbit | relation | $v_p$ of the residual | remark |
|---|---|---|---|
| $\mathbf A/\mathbf C/\mathbf F$ | $\xi_3^{\mathbf F'}=-\tfrac54\,\xi_3^{\mathbf C}$ | $786$ at $p=3$ | reproduces `ACF_ONE_SURFACE.md` Thm 2 |
| $\mathbf A/\mathbf C/\mathbf F$ | $\xi_2^{\mathbf A}=0$ | $v_2\ge1191$, precision $1170$ | reproduces `ACF_ONE_SURFACE.md` Thm 3(c) |
| $\sqrt{\text{Domb}}$ | $\xi_p^{(2)}=-2\,\xi_p^{(3)}$ | $792$ at $p=2$ **and** $386$ at $p=3$ | **new** |
| $\sqrt{\text{Domb}}$ | $\xi_p^{(4)}=\xi_p^{(6)}$ | $789$ at $p=2$ **and** $390$ at $p=3$ | **new** |
| $\sqrt{\text{Domb}}$ | $\xi_2^{(5)}=2\,\xi_2^{(7)}$ | $2387$ at $p=2$ | **new** |
| $\sqrt{\text{Domb}}$ | $\xi_2^{(1)}=0$ | $v_2\ge795$, precision $784$ | **new**: the parent's $2$-adic limit vanishes |
| $\sqrt{s_7},\sqrt{s_{10}},\sqrt{s_{18}}$ | none | residuals $38$–$43$ only | no relation of height $\le10^6$ |

**The archimedean control (`12_domb_pairs.gp`).** None of the three
$p$-adically related pairs is archimedean-rationally related: with $\xi$ computed
to $\ge300$ digits, the best rational of height $\le10^{15}$ misses by
$9.3\cdot10^{-31}$, $5.5\cdot10^{-31}$, $1.9\cdot10^{-30}$ respectively --- the
random level $\approx$ height$^{-2}$. For the pair $(5,7)$ this is provable from
§5.3: $\xi^{(5)}/\xi^{(7)}=2(15L-6\zeta(2))/(6\zeta(2)+15L)$ is irrational unless
$\zeta(2)/L(2,\chi_{-3})\in\mathbf Q$.

**Reading.**

1. The pairs $(2,3)$ and $(4,6)$ of the $\sqrt{\text{Domb}}$ orbit are the exact
   realisation of the programme's slogan: their **archimedean** periods are
   $\mathbf Q$-independent (§5.2) while their **$p$-adic** limits agree up to
   $-2$ resp. $1$ at *two different primes at once*, and their $\lambda_2$'s
   differ ($-12$ versus $-4$; $-12$ versus $-4$). *Same $p$-adic limits,
   different $\lambda_2$* — but with different archimedean periods, which is not
   what the programme hoped for.
2. The pair $(5,7)$ shows the same at $p=2$ with a doubled slope: the Cauchy
   precision of $.5,.7$ is $\approx6N$ ($2368$ at $N=400$) against $\approx2N$
   for the parent, so the gauge changes the $p$-adic slope even though $v_2(d)=6$
   throughout. Their archimedean periods, by §5.3, are the two *different*
   Eisenstein combinations $\tfrac1{16}(15L-6\zeta(2))$ and
   $\tfrac1{32}(6\zeta(2)+15L)$ — related $2$-adically by the factor $2$ but not
   over $\mathbf R$.
3. $\xi_2^{\sqrt{\text{Domb}}}=0$ is the $L(f_{12},2)$ analogue of
   $\xi_2^{\mathbf A}=0$; whether it has the same Euler-factor explanation
   (`EULER_CRITERION.md` Theorem F: $P(s)$ equal to the Euler factor plus an odd
   co-divisor character) was **not** checked here and is an open item.

---

## 7. Task 2: irrationality measures — no gain, quantified

For a single row with denominator exponent $k$, the linear form
$L_n=d_n^k(a_n\xi-b_n)$ has integer coefficients of size $e^{\sigma n}$ and
absolute value $e^{-\delta n}$ with
$$\sigma=k+\log\lambda_1,\qquad \delta=-k-\log|\lambda_2|=\operatorname{score},$$
whence $\mu(\xi)\le1+\sigma/\delta$ whenever $\delta>0$.
(The formula quoted in the task statement, $1+\log\lambda_1/(-\log\lambda_2-k)$,
omits the $+k$ in the numerator, which the $d_n^k$ contributes.)

**The two-sequence lemma.** Given two families with parameters
$(\sigma_i,\delta_i)$ and integer coefficients, merge them into one sequence
ordered by denominator. At scale $Q$ the best available form has
$|L|\le Q^{-\max_i(\delta_i/\sigma_i)}$, so
$$\mu(\xi)\ \le\ 1+\min_i\frac{\sigma_i}{\delta_i},$$
i.e. **exactly the better of the two, with no convex-combination gain at leading
order.** (A gain would require combining the two forms arithmetically — the
correlated-lattice construction of `paper/sections/05_two_row.tex` §5.1 — which
is a different mechanism; see below.)

### 7.1 The numbers (`07_measures.gp`)

| construction | $\lambda_1$ | $\lambda_2$ | $k$ | $\sigma$ | $\delta$ | $\mu\le$ |
|---|---|---|---|---|---|---|
| Apéry $\zeta(2)$ (Zagier $\mathbf D$), $v=0$ | $11.09017$ | $-0.0901699$ | 2 | $4.406059$ | $0.4060591$ | $\mathbf{11.85078219105}$ |
| its conjugate placement, $v=\mu$ (over $\mathbf Q(\sqrt5)$) | $11.18034$ | $0.0901699$ | 2 | $4.414157$ | $0.4060591$ | $11.87072446756$ |
| Beukers $\sqrt{\text{Apéry}}$, $v=0$ | $135.88225$ | $0.1177490$ | 2 | $6.911789$ | $0.1391999870$ | $\mathbf{50.65365917218}$ |
| its conjugate placement (over $\mathbf Q(\sqrt2)$) | $135.76450$ | $-0.1177490$ | 2 | $6.910922$ | $0.1391999870$ | $50.64743124451$ |

Merged bound: $11.85078219105$ for $\zeta(2)$ (gain $0$) and $50.64743124451$
for Beukers (gain $0.00622793$). The $\zeta(2)$ value reproduces Apéry's
classical $11.85078\ldots$ exactly, which calibrates the formula.

### 7.2 Why even the $0.0062$ is not collectable — three independent reasons

1. **The conjugate placement approximates a different number** (`09_conj_xi.gp`,
   $n=260$, Cauchy tail $<10^{-231}$):
   $$\xi^{\sharp}_{\mathbf D}=0.32781709442032838316549\ldots\ \ne\ \zeta(2)/5=0.32898681336964528729448\ldots,$$
   $$\xi^{\sharp}_{\text{Beukers}}=0.10021505332616528071910\ldots\ \ne\ 0.10018744922933940616776\ldots,$$
   with `lindep` on $[\xi^{\sharp},\xi,1]$ returning $20$-digit junk in both
   cases (and on $[\xi^{\sharp},\zeta(2),1]$ for $\mathbf D$). So the two
   placements are not two sequences for one target at all.
2. **Coefficients live in $\mathcal O_K$, and the norm rate is positive.** For
   an $\mathcal O_K$-form to prove a rational value irrational one needs
   $|N(\ell_n)|\to0$. The two conjugate placements have
   $|\lambda_2^{(v=\mu)}\cdot\lambda_2^{(v=\lambda)}|=|\mu|\cdot|\lambda|=|d|=1$
   for Zagier $\mathbf D$ and $|\mu(\lambda-\mu)|=15.98614$ for Beukers, and each
   place costs $d_n^{2}$, so the norm rate is
   $$\log|d|+2k=4\quad(\text{Apéry }\zeta(2)),\qquad
     \log 15.98614+4=6.771722\quad(\text{Beukers}),$$
   both $>0$: the norm grows, no descent to $\mathbf Q$.
3. **Theorem 4**: no *rational* second placement exists in either case, by
   arithmetic necessity.

### 7.3 The one same-period pair, $(\mathbf C,\mathbf F')$

$\xi^{\mathbf F'}=-\tfrac54\xi^{\mathbf C}$ exactly, so these two rows **are**
two sequences of approximations to $L(2,\chi_{-3})$ with different $\lambda_2$
($1$ and $8$). But
$$\sigma_{\mathbf C}=\sigma_{\mathbf F'}=2+\log9=4.197225,\qquad
\delta_{\mathbf C}=-2,\quad \delta_{\mathbf F'}=-4.079442,$$
both negative: the merged bound is vacuous. Deficits: $2$ nats/step to make
$\mathbf C$ decay, $4.079$ for $\mathbf F'$. And the project's own two-row
machinery does not start either: the design rule
(`paper/sections/05_two_row.tex` Prop. "design rule") needs a **decayer** with
$|\lambda|<1$, and by Theorem 4 **no member of any rational orbit in the corpus
has $|\lambda_2|<1$** — the minimum over all $49$ rows is exactly $1$. The only
rows in the corpus with $|\lambda_2|<1$ are Apéry $\zeta(2)$ ($0.09017$),
Beukers ($0.11775$), $\sqrt T$ ($0.68629$) and $\Gamma_0(5)+5$ ($0.72136$), and
all four have irrational characteristic roots.

---

## 8. Task 3: simultaneous approximation and the Nesterenko exponent

Calegari–Dimitrov–Tang prove the linear independence of $1,\zeta(2),L(2,\chi_{-3})$
over $\mathbf Q$, and their host is precisely Zagier $\mathbf C$'s curve
$Y_0(6)$ with $t_1=\tfrac19$, $t_2=1$ (`CDT_FINDER.md` §2–§3). The question here
is how far a classical Nesterenko-criterion argument built from the cusp-move
orbit gets.

**Criterion used.** If integer linear forms $L_n=q_{0,n}+q_{1,n}\theta_1+\dots+q_{m,n}\theta_m$
satisfy $\max_j|q_{j,n}|=e^{\sigma n(1+o(1))}$ and $|L_n|=e^{-\delta n(1+o(1))}$
with matching upper and lower bounds, then
$\dim_{\mathbf Q}\langle1,\theta_1,\dots,\theta_m\rangle\ge1+\tau$, $\tau=\delta/\sigma$.

### 8.1 The pair $(\mathbf A,\mathbf C)$

Row $\mathbf A$ gives $L^{\mathbf A}_n=d_n^2a_n^{\mathbf A}\zeta(2)-4d_n^2b_n^{\mathbf A}\in\mathbf Z[1,\zeta(2)]$,
row $\mathbf C$ gives $L^{\mathbf C}_n=d_n^2a_n^{\mathbf C}L(2,\chi_{-3})-2d_n^2b_n^{\mathbf C}$;
both are forms in $1,\zeta(2),L(2,\chi_{-3})$, and their coefficient vectors in
the $(\zeta(2),L)$-plane are independent.

$$\sigma_{\mathbf A}=2+\log8=4.079442,\quad \sigma_{\mathbf C}=2+\log9=4.197225,\quad
\log|L^{\mathbf A}_n|/n=\log|L^{\mathbf C}_n|/n=+2 .$$
$$\boxed{\ \sigma=4.197225,\qquad \delta=-2,\qquad \tau=-0.4765053580,\qquad
\dim\ \ge\ 1+\tau=0.5234946420\ }$$
— vacuous. **Deficits.**

| target | requirement | deficit (nats per unit $n$) |
|---|---|---|
| $\dim\ge2$ (one of them irrational) | $\delta>0$ | $\mathbf{2.000000}$ |
| $\dim\ge3$ (CDT) | $\delta\ge2\sigma=8.394449$ | $\mathbf{10.394449}$ |

**The best conceivable version.** Replace $|\lambda_2|$ by the Fricke-perfect
$1/\lambda_1$, i.e. use the *budget* of `paper/tables/census`:
$\text{budget}_{\mathbf A}=\log8-2=0.0794415$,
$\text{budget}_{\mathbf C}=\log9-2=0.1972246$. Then
$\tau=0.0189272$ and $\dim\ge1.0189272$: enough for **one** irrationality (which
is known), still short of $\dim\ge3$ by $1.98107$ in $\tau$, i.e. $8.315008$
nats/step.

### 8.2 The single-row version from §5.3

Rows $\sqrt{\text{Domb}}.5$ and $.7$ each give a form involving all three
numbers at once, with independent coefficient directions $(-6,15)$ and $(6,15)$:
$$16\,d_n^2\bigl(a_n\xi-b_n\bigr)=d_n^2a_n\bigl(15L-6\zeta(2)\bigr)-16\,d_n^2b_n,\qquad
32\,d_n^2\bigl(a_n\xi-b_n\bigr)=d_n^2a_n\bigl(6\zeta(2)+15L\bigr)-32\,d_n^2b_n .$$
Here $\sigma=2+\log16=4.772589$, $\delta=-2-\log4=-3.386294$,
$$\tau=-0.7095299,\qquad \dim\ge0.2904701,$$
deficit $3.386294$ to $\dim\ge2$ and $12.931472$ to $\dim\ge3$ — structurally the
right object, numerically worse than $(\mathbf A,\mathbf C)$ because
$\lambda_1=16$ costs more than it gains.

### 8.3 Simultaneous approximation with a common denominator

With $q_n=d_n^2a_n^{\mathbf A}a_n^{\mathbf C}$,
$p_n^{(1)}=4d_n^2b_n^{\mathbf A}a_n^{\mathbf C}$,
$p_n^{(2)}=2d_n^2b_n^{\mathbf C}a_n^{\mathbf A}$ (all integers),
$$\frac{\log q_n}{n}\to2+\log72=6.276666,\qquad
\frac{\log|q_n\zeta(2)-p^{(1)}_n|}{n}\to 2+\log9=4.197225,$$
$$\frac{\log|q_nL(2,\chi_{-3})-p^{(2)}_n|}{n}\to 2+\log8=4.079442 .$$
Simultaneous-approximation exponent $-0.6687029$; the dual Nesterenko criterion
in dimension $2$ needs $>\tfrac12$. Again the forms grow.

### 8.4 Honest comparison with CDT

CDT's margin on this very row is $+0.0053$ (their $m\le13.9938<14$), reproduced
independently in `CDT_FINDER.md` §2. Their input is not a linear form of
Nesterenko type at all but the Bost–Charles arithmetic-holonomy bound applied to
a $14$-dimensional function module on $Y_0(6)$. **The classical Padé/Nesterenko
route from the cusp-move orbit misses by $10.394$ nats per step; there is no
sense in which the orbit "gets anywhere near".** What the orbit *does* supply
that CDT's setup does not is §5.3: explicit integral rows whose Apéry limit is a
$\mathbf Q$-combination of the two periods, which is a natural input for a
future two-row/correlated-lattice attempt (`paper/sections/05_two_row.tex`
Theorem E) rather than for a scalar criterion.

---

## 9. What this does **not** prove

* **Theorem 1 is proved; Theorems 3 and 4 are proved; the census is verified,
  not certified.** Integrality is checked to $n=200$, sharpness of $k$ to
  $n=200$, algebraic integrality over $\mathcal O_K$ to $n=60$. No induction is
  written out.
* **The two new identifications of §5.3 are numerical.** They hold to the full
  truncation precision available ($\sim120$ digits), which is overwhelming, but a
  proof would require identifying the weight-three source of $(20,10,64)$ and
  $(20,4,64)$ as an explicit Eisenstein combination in
  $M_3(\Gamma_0(12),\chi_{-3})$ and running Theorem B$^*$. That eta-quotient
  identification was **not** done here and is the obvious next step.
* **All $p$-adic statements of §6 are numerical**, at $N=400$; the residual
  valuations quoted are at the rows' own Cauchy precision, which is the strongest
  numerical status available, but no $p$-adic transfer theorem in the style of
  `ACF_ONE_SURFACE.md` Theorem 2 was proved for the $\sqrt{\text{Domb}}$
  relations. In particular the $\mathbf Q(\sqrt{\phantom{x}})$-free statement
  "$\xi_p^{(2)}=-2\xi_p^{(3)}$ at $p=2$ and $p=3$" is a measurement.
* **Theorem 4 assumes $a,d\in\mathbf Z$ after rescaling.** It says nothing about
  rows whose leading coefficients are not integral in any normalisation, nor
  about higher-order ($w\ge2$) systems, where the "three placements" picture has
  more than three finite singular points and the argument does not apply.
* **The failure of the two-sequence lemma to improve $\mu$ is a leading-order
  statement.** It does not exclude a genuinely arithmetic combination of two
  placements — the correlated congruence lattice of
  `paper/sections/05_two_row.tex` §5.1 — but §7.3 shows that construction has no
  decayer to work with inside any rational orbit, so it cannot even be set up
  here.
* **$\xi_2^{\sqrt{\text{Domb}}}=0$ was measured, not explained.** Whether
  `EULER_CRITERION.md` Theorem F applies (as it does to $\xi_2^{\mathbf A}=0$)
  is open.
* **Nothing here is an irrationality theorem.** Every score computed in this
  document is negative except the two that were already known, and for those the
  cusp move produces no new rational data.

---

## 10. Reproduction

```
cd lattice/cusp_move && ./run_all.sh          # writes out/*.log
```

| script | what | runtime |
|---|---|---|
| `lib.gp` | rows, sequences, the transform $a\mapsto a^{\sharp}$, recurrence fitting, $k$, scaling | — |
| `rows.gp` | the corpus: Zagier A–F, the six $\sqrt{\ }$-rows, Cooper $\sqrt{s_7},\sqrt{s_{10}},\sqrt{s_{18}}$, $\Gamma_0(5)+5$, Herfurtner #30 and #45 | — |
| `00_selftest.gp` | Theorem 1 and Theorem 2 on $\mathbf A,\mathbf C,\mathbf F$ | seconds |
| `01_general_move.py` | **symbolic proof of Theorem 1** (sympy); the $\alpha$-condition | ~1 min |
| `02_orbits.gp` | closed formula vs. binomial transform, all rows, over $\mathbf Q$ and over $K$ | seconds |
| `03_census.gp` | orbit closure, integrality to $n=200$, $k$ sharp, score, $\xi$, $p$-adic | ~10 min |
| `04_periods.gp` | $\xi$ to $\ge45$ digits, rational relations, `lindep` identification | ~15 min |
| `05_placements.gp` | the gap picture and Theorem 4 | seconds |
| `06_padic.gp` | $p$-adic Apéry limits and pairwise relations, $N=400$ | ~10 min |
| `07_measures.gp` | irrationality measures, two-sequence bound, Nesterenko exponents | seconds |
| `08_quadratic.gp` | placements over $K$, algebraic integrality, corrected scores | ~1 min |
| `09_conj_xi.gp` | the conjugate placements' Apéry limits | ~1 min |
| `10_domb_ident.gp` | the two new identifications at $200$ digits | ~1 min |
| `11_companion_general.gp` | Theorem 2 in all 22 (row, $\lambda$, $\alpha$) cases | seconds |
| `12_domb_pairs.gp` | archimedean ratios of the $p$-adically related $\sqrt{\text{Domb}}$ pairs | ~2 min |

PARI/GP 2.15.4, exact rational (and $\mathbf Q(y)/(y^2-ay+d)$) arithmetic
throughout; `realprecision` up to $400$ where floating point is used.
