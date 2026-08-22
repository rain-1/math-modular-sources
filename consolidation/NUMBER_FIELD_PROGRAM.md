# Number-field periods and mixed statements: the Beukers-1987-Theorem-4 mechanism, made into a template

*Claude (Fable), 2026-08-22.  Scripts: `lattice/number_field/`.
Sources read: `papers/beukers1987.txt` (Theorems 1–5, in full),
`consolidation/L3CHI5_TWO_WORLDS.md`, `AESZ207.md`, `MUM_SURVEY.md`,
`CDT_FINDER.md`, `HERFURTNER_CLASSIFICATION.md`, `THEOREM_B_EXACT.md`,
`ROOT_ROWS.md`, `SLOPE_CENSUS.md`, `paper/sections/05_two_row.tex`.
Tags: **[proved]**, **[verified over range R]**, **[numerical, D digits]**, **[open]**.*

---

## 0. Verdict

1. **The mechanism of Beukers' Theorem 4 is a norm inequality, and it can be
   written down once and for all.**  For an Apéry-type row whose coefficients lie
   in the ring of integers of a number field $K$ of degree $g$, with sharp
   lcm-denominator exponent $k$, the criterion is
   $$\boxed{\ \frac1g\sum_{v\mid\infty}d_v\log\rho_v\;>\;k+\log D\ }$$
   where $\rho_v$ is the radius of convergence, at the archimedean place $v$, of
   the *linear-form* generating function, $d_v\in\{1,2\}$ is the local degree and
   $D$ is the exponential denominator rate ($D=1$ throughout this corpus).
   The whole content of a particular theorem is *which places admit the fold*:
   at a folded place $\rho_v=1/|\lambda_2^{(v)}|$, at an unfolded place
   $\rho_v=1/|\lambda_1^{(v)}|$.  §1, Theorem NF-1.
2. **Beukers' Theorem 4 folds at one place only, and that is forced.**  At the
   second place the fold constant of the conjugate construction is the image of
   the *formula* under $\sqrt5\mapsto-\sqrt5$, not the Galois conjugate of the
   number one has assumed to lie in $K$; the two agree only if one already knows
   what one is trying to prove.  So the second place contributes the crude
   radius.  Numerically **[verified, exact]**
   $$\sqrt{\rho_{v_1}\rho_{v_2}}=\sqrt{493.9654530919\ldots\times0.8351896967\ldots}
     =20.311446452183879\ldots\;>\;e^3=20.085536923187667\ldots,$$
   a margin of $1.1247\%$, i.e. Mode-I margin $+0.01118459177540946495$ per place
   per index ($0.02236918$ over the two).  Had the conjugate fold been available the figure would be
   $\sqrt{591.44\ldots}=24.32$, a $21\%$ margin.  **Theorem 4 is true by one
   percent, and the missing percent is exactly the Galois-conjugate-versus-
   analytic-conjugate gap.**  §1.4.
3. **A complete no-go for the complex-fold rows.**  A row over $\mathbf Q$ (or a
   place of $K$) with $|\lambda_1|=|\lambda_2|$ has two singularities of equal
   modulus; the fold cancels one of them and the radius is unchanged.  Hence
   $\rho_v=1/|\lambda|$ and the linear form does not decay at all.  The fold
   constant $\xi\in\mathbf C$ of Theorem B\* is still well defined and still a
   period, but it is approached only at rate $O(1/n)$, by the
   $\mathbf Q(\lambda)$-rationals $(b_{n+1}-\bar\lambda b_n)/(a_{n+1}-\bar\lambda a_n)$.
   **No Diophantine statement of Beukers' shape is available for
   $L(\chi_5,3)$, $L(\chi_{-3},2)$-via-$\mathbf B$, $\zeta(3)$-via-$\delta$, or
   for either of the two new Herfurtner rows.**  §1.5, Theorem NF-0.
4. **Two corrections to earlier documents.**  (i) `CDT_FINDER.md` §6 analyses
   "$X_1(5)\ \mathrm{Sym}^2$" — the Cauchy square of Zagier $\mathbf D$ over
   $\mathbf Q$, characteristic roots $\varphi^5,\varphi^5,-\varphi^{-5},-\varphi^{-5}$ —
   but *that is not Beukers' row*.  Beukers' row is defined over
   $\mathbf Z[\varphi]$, is of Apéry $(R3)$ shape with $c=1$, and has
   characteristic roots $\lambda^2-(248+110\sqrt5)\lambda+1$.  (ii) The
   "normalised" number-field budget $\lambda_2^{\rm norm}=|N(\lambda_2)|^{1/g}$
   of `CDT_FINDER.md` §3 is the *both-places-folded* accounting; it is the
   correct one only in **Mode II** of §1.3, and Beukers' Theorem 4 is
   precisely an instance where Mode II's hypothesis is unavailable.  (iii) The
   field norm is the wrong object anyway: $\lambda_2\notin K$ here — the
   characteristic roots of Beukers' row generate the **cyclic quartic field
   $\mathbf Q(\zeta_{20})^+$** — so what enters is the *place-by-place* product
   $\prod_v|\lambda_2^{(v)}|^{d_v}$, which coincides with a field norm only when
   $\lambda_2\in K$.  §1.4, §2.
5. **The two new Herfurtner rows are uniformised exactly, and one of their fold
   constants is identified — as a CM period, not an $L$-value.**  Row #30
   ($\Gamma_0(7)$) satisfies
   $A(t)(1-117t+3969t^2)^{1/3}=\theta_{-7}$ with $t=-\tfrac19(\eta(7\tau)/\eta(\tau))^4$;
   row #45 satisfies $j=(24t-1)^3/(8t^3)$.  Their folds sit at the classical CM
   points $j=0$ (disc $-3$) and $j=1728$ (disc $-4$), and
   $$\operatorname{Im}\xi(\#45)=-\frac{\Gamma(1/4)^4}{2^{7/2}3^{11/4}\pi}
   =-\frac{\varpi^2}{\sqrt2\,3^{11/4}}\qquad\textbf{[verified to }10^{-249}\textbf{]},$$
   the lemniscatic Chowla–Selberg period.  Everything else resists a $340$-target
   sweep, and the reason is structural: the sources are **not** Eisenstein — but
   the untwisted combination is, and
   $t\,\theta_{-7}^3/(1-117t+3969t^2)=-\tfrac19\eta(\tau)^3\eta(7\tau)^3$, the
   newform 7.3.b.a **[verified to $q^{80}$]**.  §3.
6. **Beukers' row is unique, and the integral locus has only seven seeds.**  An
   exhaustive scan of the two Apéry shapes over all real quadratic fields
   $\mathbf Q(\sqrt D)$, $2\le D\le97$ — $\approx4.26\cdot10^{14}$ parameter
   tuples, $8017$ integral classes re-verified exactly to $n=200$, zero false
   positives — finds **exactly one $\mathrm{GL}_1$-orbit with a positive Mode-I
   margin**, and it is Beukers' Theorem-4 row, reporting
   $+0.011184591775409465$ along the whole orbit.  For the Zagier shape the count
   is **zero**.  The runner-up is $-0.682$ (Beukers rescaled by a non-unit) and
   everything else $\le-3.0$.  Structurally, *everything* integral in the box is a
   $K$-rescaling of one of **seven seeds**: the six rational AZZ $\zeta(3)$-type
   rows, or Zagier's six, plus Beukers' single genuinely quadratic row.  §4.2a.
7. **The corpus is scored in §2**; §4 is the shopping list — what a row would
   have to look like for the template to give something new.

---

## 1. The template

Throughout: $d_n=\operatorname{lcm}(1,\dots,n)$, so $d_n=e^{n(1+o(1))}$ by the
prime number theorem.  $K$ is a number field of degree $g=[K:\mathbf Q]$ with
archimedean places $v\mid\infty$ and local degrees $d_v$ ($\sum_vd_v=g$);
$|\cdot|_v$ is the *absolute value* (not the normalised one), so
$|N_{K/\mathbf Q}(x)|=\prod_v|x|_v^{d_v}$.

### 1.1 Hypotheses

Let $(a_n)_{n\ge0}$, $(b_n)_{n\ge0}$ be sequences in $K$ satisfying one common
linear recurrence over $K$ with polynomial coefficients, and let $\xi$ be a
complex number.  Write $\ell_n:=b_n-\xi a_n$.

* **(NF-int) Integrality.**  $D^na_n\in\mathcal O_K$ and $D^nd_n^k\,b_n\in\mathcal O_K$
  for fixed $D\in\mathbf Z_{\ge1}$ and $k\in\mathbf Z_{\ge1}$ ($k$ *sharp*: the
  statement fails for $k-1$).
* **(NF-rad) Place-by-place radii.**  For each $v\mid\infty$ there is $\rho_v>0$ with
  $$\limsup_n|\ell_n|_v^{1/n}\le\rho_v^{-1}.$$
* **(NF-ne) Non-vanishing.**  $\ell_n\ne0$ for infinitely many $n$.

In every instance below (NF-rad) is supplied by singularity analysis of the
Picard–Fuchs operator: the generating function $\sum_n\ell_nt^n$ is
$y_B(t)-\xi\,y_0(t)$ and $\rho_v$ is its radius of convergence in the embedding
$v$, so $\rho_v$ is the modulus of the nearest singularity of $y_B-\xi y_0$ at
$v$.  Since the finite singularities of the operator are $t_i=1/\lambda_i$
($\lambda_i$ the characteristic roots),

$$\rho_v=\begin{cases}
1/|\lambda_2^{(v)}|,&\text{if the fold cancels the inner singularity at }v,\\[2pt]
1/|\lambda_1^{(v)}|,&\text{otherwise,}
\end{cases}
\qquad |\lambda_1^{(v)}|\ge|\lambda_2^{(v)}| .$$

("The fold" is Beukers' §1 mechanism: an involution of the $t$-line fixing the
inner singularity, under which $E\cdot(f-\xi)$ is invariant, so that the
would-be branch point is not one.  Archimedeanly it is the Fricke fold of
`THEOREM_B_EXACT.md` Theorem II for third-order rows and the cusp mechanism of
Theorem I for second-order rows.)

### 1.2 Theorem NF-1 (one target, norm criterion)

> **Theorem NF-1.** **[proved]**  Assume (NF-int), (NF-rad), (NF-ne) and
> $$\frac1g\sum_{v\mid\infty}d_v\log\rho_v\;>\;k+\log D .$$
> Then $\xi\notin K$.

*Proof.*  Suppose $\xi\in K$ and let $\Delta\in\mathbf Z_{\ge1}$ with
$\Delta\xi\in\mathcal O_K$.  Put $L_n:=\Delta\,D^n d_n^k\,\ell_n$.  By (NF-int)
$L_n\in\mathcal O_K$, and by (NF-ne) $L_n\ne0$ for infinitely many $n$, whence
$|N_{K/\mathbf Q}(L_n)|\ge1$ for those $n$.  On the other hand, for every
$\varepsilon>0$ and all large $n$,
$$|N(L_n)|=\prod_v|L_n|_v^{d_v}
\le|N(\Delta)|\cdot\bigl(D^ne^{kn(1+\varepsilon)}\bigr)^{g}
\prod_v\rho_v^{-nd_v(1-\varepsilon)}
=|N(\Delta)|\Bigl[\frac{D\,e^{k}}{\bigl(\prod_v\rho_v^{d_v}\bigr)^{1/g}}\Bigr]^{gn(1+O(\varepsilon))},$$
which tends to $0$ under the hypothesis.  Contradiction. $\square$

**Remarks.**
* $g=1$, one real place, $D=1$: the criterion is $\log(1/|\lambda_2|)>k$ — the
  project's **score** $\log(1/|\lambda_2|)-k>0$ of `ROOT_ROWS.md` Theorem R4 and
  `HERFURTNER_CLASSIFICATION.md` Corollary H6.  So NF-1 *is* the score, read
  adelically.
* $\Delta$ enters only as a constant, not exponentially; this is why the
  criterion is insensitive to the (unknown) denominator of the hypothetical
  $\xi$.  $D$ *does* enter exponentially.
* The theorem gives $\xi\notin K$, not linear independence.  See §1.6.

### 1.3 Which places fold — the three modes

The fold at $v$ cancels the inner singularity of
$y^{(v)}_B-\xi^{(v)}y^{(v)}_0$, where $\xi^{(v)}$ is the **analytic** fold
constant of the $v$-embedded construction (the Apéry limit one would measure
from the $v$-embedded recurrence).  What the *hypothesis* supplies at $v$ is
$\sigma_v(\beta)$ for the assumed $\beta\in K$.  So

$$\boxed{\ v\ \text{folds}\iff \xi^{(v)}=\sigma_v(\beta).\ }$$

Since $\xi^{(v_0)}=\beta$ is the assumption at the distinguished place, one may
choose *which* places to assume this at, and each choice is a different theorem.

> **Mode I (Beukers, Theorem 4).**  Assume only $\beta:=\xi^{(v_0)}\in K$.  Then
> $v_0$ folds and no other place does:
> $$\rho_{v_0}=1/|\lambda_2^{(v_0)}|,\qquad \rho_v=1/|\lambda_1^{(v)}|\ (v\ne v_0).$$
> Conclusion of NF-1: $\xi^{(v_0)}\notin K$.  **Strongest conclusion, weakest
> criterion.**

> **Mode II (pair / Minkowski version).** **[proved]**  Assume the whole vector
> $\bigl(\xi^{(v)}\bigr)_{v\mid\infty}\in\prod_vK_v\cong\mathbf R^{g}$ is the
> Minkowski image of a single $\beta\in K$, i.e. $\xi^{(v)}=\sigma_v(\beta)$ for
> every $v$.  Then **every** place folds, $\rho_v=1/|\lambda_2^{(v)}|$, and
> NF-1's criterion becomes
> $$\lambda_2^{\rm norm}:=\bigl|N_{K/\mathbf Q}(\lambda_2)\bigr|^{1/g}<D^{-1}e^{-k}.$$
> Conclusion: **no such $\beta$ exists** — the point $(\xi^{(v)})_v$ is not in
> the image of $K$ under the Minkowski embedding.  **Weakest conclusion,
> strongest criterion.**  This is exactly the "normalised" accounting of
> `CDT_FINDER.md` §3, and it is legitimate *only* in this mode.

> **Mode III (partial).**  Assume $\xi^{(v)}=\sigma_v(\beta)$ for $v$ in a subset
> $T\ni v_0$ of the places; fold on $T$, crude bound off $T$.  Mode I is
> $T=\{v_0\}$ and Mode II is $T=$ all.

**What Mode II says for Beukers' pair.**  There
$\xi^{(v_1)}=-\tfrac{16}5\zeta(3)+4\sqrt5\,L(3,\chi_5)$ and
$\xi^{(v_2)}=-\tfrac{16}5\zeta(3)-4\sqrt5\,L(3,\chi_5)$, so writing
$\beta=x+y\sqrt5$ the Mode-II hypothesis reads
$x=-\tfrac{16}5\zeta(3)\in\mathbf Q$ **and** $y=4L(3,\chi_5)\in\mathbf Q$.  Its
negation is therefore only *"$\zeta(3)$ and $L(3,\chi_5)$ are not both
rational"* — true but weak (Apéry).  So on this row Mode II buys a fat margin
for a cheap conclusion and Mode I buys a thin margin for the real theorem.
**The two modes are genuinely different statements and both are correct; the
mistake to avoid is quoting Mode II's margin for Mode I's conclusion.**

### 1.3a What is actually approximated

It is worth saying explicitly what object the row approximates, because the two
cases of the task ("the pair $(\xi,\xi^\sigma)$" and "$(\operatorname{Re}\xi,
\operatorname{Im}\xi)$") are the *same* object seen at a real and at a complex
place.

Let $\iota:K\hookrightarrow K\otimes\mathbf R=\prod_{v\mid\infty}K_v\cong\mathbf R^{r_1}\times\mathbf C^{r_2}$
be the Minkowski embedding.  The approximants are $\iota(b_n/a_n)$ and the
target is the vector of analytic fold constants $\bigl(\xi^{(v)}\bigr)_v$.

* **$K$ real quadratic** ($r_1=2$): the target is the pair
  $(\xi^{(v_1)},\xi^{(v_2)})\in\mathbf R^2$ — Beukers' $(\xi,\xi^\sigma)$ — and
  the quality of the simultaneous approximation is governed by
  $\prod_v\rho_v$, i.e. by the *norm* of the linear form.  Mode II is the
  statement that this point misses $\iota(K)$.
* **row over $\mathbf Q$, complex characteristic roots** ($r_2=1$ for
  $\mathbf Q(\lambda)$): the target is the single complex number
  $\xi=\operatorname{Re}\xi+i\operatorname{Im}\xi$, and the approximants
  $\xi_n=(b_{n+1}-\bar\lambda b_n)/(a_{n+1}-\bar\lambda a_n)$ lie in the
  imaginary quadratic field $\mathbf Q(\lambda)$.  Formally the same picture —
  but the local degree is $d_v=2$ and $|\lambda_1^{(v)}|=|\lambda_2^{(v)}|$, so
  the criterion degenerates to $\rho^2>e^{2k}$, i.e. $\rho>e^k$ with
  $\rho=1/|\lambda|$: *the complex place buys exactly nothing*, and the rate is
  only $O(1/n)$ (Theorem NF-0).  The number field is present, and useless.

### 1.4 Beukers' Theorem 4 as the model instance, and why Mode II is unavailable there

Beukers works on $\Gamma_1(5)$ with the two conjugate hauptmoduls $t(\tau)$,
$s(\tau)$ whose $q$-coefficients are conjugate algebraic integers of
$K=\mathbf Q(\sqrt5)$.  The row he records is
$$(n+1)^3d_{n+1}=\bigl\{(124+55\sqrt5)\,n(n+1)+34+15\sqrt5\bigr\}(2n+1)\,d_n-n^3d_{n-1},$$
$d_0=1$, $d_1=34+15\sqrt5$, $d_2=7111+3180\sqrt5$, $d_3=2040334+912465\sqrt5$ —
an $(R3)$ Apéry-shaped row over $\mathbf Z[\varphi]$ with
$$\chi(\lambda)=\lambda^2-(248+110\sqrt5)\lambda+1,\qquad c=\lambda_1\lambda_2=1
\ \text{(Apéry-perfect, a unit)} .$$
At the two real places **[verified, 50 digits]**

| place | $\lambda_1^{(v)}$ | $\lambda_2^{(v)}$ | inner sing. | fold? | $\rho_v$ |
|---|---|---|---|---|---|
| $v_1$ ($\sqrt5\mapsto+2.2360\ldots$) | $493.9654530919042$ | $0.0020244330726788$ | $0.0020244\ldots$ | **yes** | $493.9654530919042$ |
| $v_2$ ($\sqrt5\mapsto-2.2360\ldots$) | $1.1973327782587006$ | $0.8351896967644328$ | $0.8351897\ldots$ | **no** | $0.8351896967644328$ |

(the analytic fold point is the *inner* singularity at **both** places: at $v_1$
it is $t(i/\sqrt5)=0.0020244\ldots$, at $v_2$ it is $s(i/\sqrt5)=0.83519\ldots$,
each the smaller of the two branch values there — so the failure at $v_2$ is not
a geometric one).  Hence

$$\tfrac12\log(\rho_{v_1}\rho_{v_2})=\log 20.311446452\ldots=3.0111\ldots>3
=\log e^3,\qquad\text{margin }+0.011184591775409465\ \text{per place per index}.$$

Beukers' own bookkeeping is the same statement: he multiplies his (5) and (6)
and weighs the product against $[1,\dots,n]^{6}=[1,\dots,n]^{2k}$, "the LCM cost
appearing once per place", and quotes $20.30$ against $e^3$.

**Why $v_2$ does not fold.**  The conjugate construction at $v_2$ is Beukers'
second triple $(F,E,f)$, obtained from the first by $\sqrt5\mapsto-\sqrt5$ *in
the defining formulas*; its fold constant is
$$A^{\rm an}_{v_2}=-\tfrac25\bigl(8\zeta(3)+5\sqrt5\,L(3,\chi_5)\bigr),$$
whereas what the hypothesis "$A\in K$" supplies at $v_2$ is $\sigma(A)$, the
Galois conjugate of an *assumed* element of $K$.  There is no reason for
$A^{\rm an}_{v_2}=\sigma(A)$ — an equality of that kind would already be a
statement about $\zeta(3)$ and $L(3,\chi_5)$.  Beukers is therefore obliged to
take $\theta\in\mathbf C$ arbitrary in his (6) and use the crude radius; and
this is **not** a slip:

> **Cost of the missing fold.**  $\log(1.19733/0.83519)=0.36089$ per index at
> $v_2$, i.e. the effective radius drops from $\sqrt{591.44}=24.3195$ to
> $\sqrt{412.55}=20.3114$.  With $e^3=20.0855$ sitting between the two, Theorem 4
> is true because of the crude bound and would be *comfortably* true with the
> conjugate fold.  Any attempt to strengthen Theorem 4 (CDT amplification,
> `CDT_FINDER.md` §6 (C2)) must either supply the conjugate fold or work with
> $20.31$.

**Correction to `CDT_FINDER.md` §6.**  That section computes with the Cauchy
square of Zagier $\mathbf D$ ($A_n=[t^n]F(t)^2=1,6,47,408,3745,\dots$, roots
$\varphi^5,\varphi^5,-\varphi^{-5},-\varphi^{-5}$, $\lambda_2^{\rm norm}=1$).
That is a row over $\mathbf Q$ and it is **not** Beukers' row.  Beukers' row has
$\lambda_2^{\rm norm}=|N(\lambda_2)|^{1/2}=|0.0020244\times0.83519|^{1/2}
=0.04111916395174473$, so the Mode-II (both-folds) budget is
$\log(1/0.04111916395174473)-3=+0.191280989954746221$ and the Mode-I (one-fold) budget is
$\tfrac12\log(412.5549)-3=+0.011184591775409465$.  Both are positive — which is why
Theorem 4 exists — whereas the $\mathrm{Sym}^2$ object of `CDT_FINDER.md`
scores $\log\varphi^{5}-3=-0.5936$ and could never have supported it.  The
CDT-transported margins of `CDT_FINDER.md` §6 for the "priority host" should be
recomputed on the correct row.

### 1.5 Theorem NF-0 (complex folds give nothing)

> **Theorem NF-0.** **[proved]**  Let a row over $K$ have, at an archimedean
> place $v$, two singularities of equal modulus (equivalently
> $|\lambda_1^{(v)}|=|\lambda_2^{(v)}|$; over $\mathbf R$ this means the
> characteristic roots are a complex-conjugate pair).  Then:
> 1. the fold at $t_c$ cancels one of $t_c,\bar t_c$ and leaves the other, so
>    $\rho_v=1/|\lambda^{(v)}|$ **unchanged**: the linear form
>    $b_n-\xi a_n$ does not decay;
> 2. consequently the contribution of $v$ to the criterion of NF-1 is
>    $-\log|\lambda^{(v)}|$, which is negative for every row in this corpus, and
>    no Diophantine conclusion follows;
> 3. the fold constant $\xi\in\mathbf C$ nevertheless exists and is the
>    connection constant of Theorem B\* ($\operatorname{Re}\xi=L(\Phi,w+1)$);
>    it is approximated by
>    $$\xi_n=\frac{b_{n+1}-\bar\lambda\,b_n}{a_{n+1}-\bar\lambda\,a_n}\in\mathbf Q(\lambda),
>    \qquad \xi_n-\xi=O(1/n)\ \text{times an oscillating factor},$$
>    i.e. at **polynomial** rate; and what is approximated is the *pair*
>    $(\operatorname{Re}\xi,\operatorname{Im}\xi)$ simultaneously.

*Proof.*  (1) is the definition of the radius of convergence: the folded
function $E(f-\xi)$ is regular at $t_c$ but retains the branch point $\bar t_c$
with $|\bar t_c|=|t_c|$.  (2) is NF-1 read off.  (3): write
$a_n=\alpha s^+_n+\bar\alpha s^-_n$, $b_n=\beta s^+_n+\bar\beta s^-_n$ with
$s^{\pm}_n\sim\lambda_{\pm}^nn^\theta(1+c_1/n+\cdots)$, $\lambda_\pm=\lambda,\bar\lambda$.
Then $s^-_{n+1}-\bar\lambda s^-_n=\bar\lambda^{n+1}\bigl[(n+1)^\theta-n^\theta\bigr]
=O(|\lambda|^nn^{\theta-1})$ while
$s^+_{n+1}-\bar\lambda s^+_n\sim\lambda^{n+1}(1-\bar\lambda/\lambda)n^\theta$;
so $\xi_n\to\beta/\alpha=\xi$ with relative error $O(1/n)$. $\square$

**Consequences for the corpus.**  Rows $\mathbf B$, $\delta$, $\eta$
(`THEOREM_B_EXACT.md` §5), AESZ 184, $\sqrt{\mathrm{AZ}(7,3,81)}$ and the two
new Herfurtner rows are all in this class.  This is the precise form of the
obstruction diagnosed in `L3CHI5_TWO_WORLDS.md` §5.2 ("for the class
$L(\chi_5,3)$ at $p=5$ every known realisation is a complex fold"): the
obstruction is *archimedean and structural*, not a failure of search.  A
degenerate real analogue is worth recording: a real-quadratic row with
$|\lambda_2^{(v_1)}|=|\lambda_1^{(v_2)}|$ has
$\bigl(\prod\rho_v\bigr)^{1/2}=1/|\lambda_2^{(v_1)}|^{1/2}\cdot|\lambda_1^{(v_2)}|^{-1/2}$
and the number field buys nothing; the complex case is the limit of this.

### 1.6 What the mechanism cannot do: linear independence

One might hope to convert "$\xi\notin K$" into "$1,\xi_1,\xi_2$ linearly
independent over $\mathbf Q$" by taking traces.  It fails, and the failure is
quantitative.  Suppose (best case) the folds work at both places, so
$$|\ell_n|_{v_1}\lesssim\rho_{v_1}^{-n},\qquad|\ell_n|_{v_2}\lesssim\rho_{v_2}^{-n},
\qquad \rho_{v_1}\gg1\gg\rho_{v_2}^{-1}.$$
$\ell_n$ is a $K$-linear form in $(1,\xi_1,\xi_2)$; its two $\mathbf Q$-shadows
$\operatorname{Tr}(\ell_n)$ and $\operatorname{Tr}(\sqrt5\,\ell_n)/5$ are
$\mathbf Q$-linear forms in $(1,\xi_1,\xi_2)$, but
$$|\operatorname{Tr}\ell_n|\le|\ell_n|_{v_1}+|\ell_n|_{v_2}
\lesssim\max_v\rho_v^{-n}=\rho_{v_2}^{-n},$$
i.e. the *worst* place dominates.  For Beukers' row that is
$1.19733^{-n}$ (both folds) or $0.83519^{-n}$ (as actually available), against a
denominator cost $e^{3n}=20.0855^n$: the trace forms are useless by a factor
$e^{3}/1.19733=16.78$ per index.

> **Proposition NF-3.** **[proved]**  The norm route yields "$\xi\notin K$" and
> nothing more; a linear-independence statement over $\mathbf Q$ or over $K$
> requires two $\mathbf Q(t)$-independent small forms, i.e. a genuinely rank-$\ge2$
> lattice.  This is exactly why `CDT_FINDER.md` §7 (C2) has to go through the
> Calegari–Dimitrov–Tang holonomy bound and not through Beukers' argument.

### 1.7 Summary of the template

| ingredient | symbol | where it comes from |
|---|---|---|
| field, places | $K$, $v\mid\infty$, $d_v$ | the coefficient ring of the row |
| integrality | $D^na_n\in\mathcal O_K$, $D^nd_n^kb_n\in\mathcal O_K$ | measured, sharp |
| radii | $\rho_v$ | singularity analysis + which places fold |
| criterion | $\frac1g\sum_vd_v\log\rho_v>k+\log D$ | Theorem NF-1 |
| conclusion | Mode I: $\xi^{(v_0)}\notin K$; Mode II: $(\xi^{(v)})_v\notin\iota(K)$ | §1.3 |
| no-go | $|\lambda_1^{(v)}|=|\lambda_2^{(v)}|$ at some $v$ | Theorem NF-0 |

---

## 2. The corpus, scored

*(computed in `lattice/number_field/`: `beukers_thm4*.gp`, `beukers_field.gp`,
`corpus.gp`, `corpus2.gp`, `aesz207*.gp`, `rootfields.gp`, `phi_id.gp`, each with
a `.log`; the exact algebra over $K$ is done in `Mod(y,y^2-y-1)$, no floats.)*

### 2.1 Beukers' row, re-derived

**[verified, exact]** The printed values $d_1=34+15\sqrt5$, $d_2=7111+3180\sqrt5$,
$d_3=2040334+912465\sqrt5$ all come out of the recurrence.  Stronger than
Beukers states: the $\sqrt5$-coefficient of $d_n$ is even for every $n\le300$, so
$d_n\in\mathbf Z[\sqrt5]$, not merely $\mathbf Z[\varphi]$.

**Characteristic polynomial, derived not assumed.**  Dividing by $n^3$ and letting
$n\to\infty$ gives $\lambda^2-(248+110\sqrt5)\lambda+1$, and

$$\boxed{\ 248+110\sqrt5=4\varphi^{10}+2,\qquad
\lambda_1=\bigl(\varphi^5+\sqrt{\varphi^{10}+1}\bigr)^2=493.965453091904187844\ldots\ }$$

**[exact identity]** — Beukers' branch value $(\lambda_2-\sqrt{1+\lambda_2^2})^2$ with
his $\lambda_2=-\tfrac{11}2-\tfrac52\sqrt5=-\varphi^5$.  So the whole row is built
out of Zagier $\mathbf D$'s $\varphi^5$.

**The field.**  $\lambda\notin K$: the discriminant $(248+110\sqrt5)^2-4=122000+54560\sqrt5$
is not a square in $K$ (`nfroots` empty) **[exact]**.  The minimal polynomial of
$\lambda$ over $\mathbf Q$ is $x^4-496x^3+1006x^2-496x+1$, irreducible, **cyclic
quartic**, totally real, of discriminant $2000=2^4 5^3$, with
$\mathrm{polredabs}=x^4-5x^2+5$: $\lambda$ generates $\mathbf Q(\zeta_{20})^+$, and
$N_{\mathbf Q(\zeta_{20})^+/\mathbf Q}(\lambda)=1$ ($\lambda$ is a unit) **[exact]**.
Consequently the shorthand $|N_{K/\mathbf Q}(\lambda_2)|^{1/g}$ is **not defined**
here; the correct object is the place-by-place product

$$\lambda_2^{\rm norm}:=\Bigl(\prod_{v\mid\infty}\bigl|\lambda_2^{(v)}\bigr|^{d_v}\Bigr)^{1/g}
=\bigl(0.00202443307\ldots\times0.83518969676\ldots\bigr)^{1/2}=0.0411191639517447\ldots$$

**[floating, 60 digits]**, and its reciprocal is exactly the $24.3195606110461522$ of
Mode II.  In particular Beukers' row is **not** Apéry-perfect in the
$|N(t_2)|=1$ sense — unlike the $\mathrm{Sym}^2$ object of `CDT_FINDER.md` §6,
where $|N|=1$ holds and the budget is negative.

**Denominators.**  Companion $c_n$ ($c_0=0,c_1=1$), $c_2=\tfrac{423}4+\tfrac{375}8\sqrt5$;
$\operatorname{lcm}(1..n)^3c_n\in\mathbf Z[\varphi]$ for every $n\le200$ and
$\operatorname{lcm}(1..n)^2c_2\notin\mathbf Z[\varphi]$: **$k=3$, sharp**
**[verified $n\le200$]**.  Beukers' integrality claim and its sharpness confirmed.

**The archimedean limit.**  At $v_1$, $c_{300}/d_{300}=0.0148059474735170089225665\ldots$
(stable to 200 digits from $n=200$ on), and `lindep` at 90 digits, reconfirmed at
200, returns $[4,0,-8,5]$:

$$\boxed{\ \lim_{n}\frac{c_n}{d_n}\Big|_{v_1}=2\zeta(3)-\tfrac54\sqrt5\,L(3,\chi_5)
=\tfrac14\bigl(8\zeta(3)-5\sqrt5\,L(3,\chi_5)\bigr)\ }$$

**[verified, 200 digits]** with **no rational constant term** — so the limit is
exactly $-\tfrac58$ times Beukers' $A$, and Theorem 4's statement
"$8\zeta(3)-5\sqrt5L(3,\chi_5)\notin\mathbf Q(\sqrt5)$" is the statement
$\lim c_n/d_n\notin K$.

### 2.2 The scoring table

Every row below **except Beukers' has rational coefficients**, hence exactly one
archimedean place: for those, the criterion of Theorem NF-1 collapses to the
project's single-place score $\log(1/|\lambda_2|)-k$, and $D=1$ throughout
($u_n\in\mathbf Z$, checked to $n\le150$–$260$).  $k$ is **measured** in every
row (companion $b_0=0,b_1=1$; for the rank-4 rows the initial vector
$(0,1,0,0)$), verified for all $n\le150$ together with the first $n$ at which
$\operatorname{lcm}^{k-1}b_n$ fails.

| # | row | char. poly | $\lvert\lambda_1\rvert$ | $\lvert\lambda_2\rvert$ | root field | $k$ | $\rho$ | **margin** | verdict |
|---|---|---|---|---|---|---|---|---|---|
| 1 | Zagier $\mathbf B$ | $x^2-9x+27$ | $5.19615$ | $=\lvert\lambda_1\rvert$ | $\mathbf Q(\sqrt{-3})$, 1 cplx | 2 | $1/5.19615$ | $-3.6479$ | FAIL (NF-0) |
| 2 | Zagier $\mathbf D$ | $x^2-11x-1$ | $11.09017$ | $0.09017$ | $\mathbf Q(\sqrt5)$ | 2 | $11.09017$ | $\mathbf{+0.40606}$ | **PASS** |
| 3 | Zagier $\mathbf E$ | $(x-8)(x-4)$ | $8$ | $4$ | $\mathbf Q$ (splits) | 2 | $1/4$ | $-3.38629$ | FAIL |
| 4 | **Herfurtner #30** | $x^2-117x+3969$ | $63$ | $=\lvert\lambda_1\rvert$ | $\mathbf Q(\sqrt{-3})$, 1 cplx | 2 | $1/63$ | $-6.1431$ | FAIL (NF-0) |
| 5 | **Herfurtner #45** | $x^2-72x+1728$ | $41.5692$ | $=\lvert\lambda_1\rvert$ | $\mathbf Q(\sqrt{-3})$, 1 cplx | 2 | $1/41.5692$ | $-5.7274$ | FAIL (NF-0) |
| 6 | $\sqrt{\mathrm{AZ}(7,3,81)}$ | $x^2-56x+1296$ | $36$ | $=\lvert\lambda_1\rvert$ | $\mathbf Q(\sqrt{-2})$, 1 cplx | 2 | $1/36$ | $-5.5835$ | FAIL (NF-0) |
| 7 | Beukers $\sqrt{\text{Apéry}}$ | $x^2-136x+16$ | $135.88225$ | $0.117749$ | $\mathbf Q(\sqrt2)$ | 2 | $8.49264$ | $\mathbf{+0.13920}$ | **PASS** |
| 8 | Apéry $\zeta(3)$ | $x^2-34x+1$ | $33.97056$ | $0.0294373$ | $\mathbf Q(\sqrt2)$ | 3 | $33.97056$ | $\mathbf{+0.52549}$ | **PASS** |
| 9 | AZ $\eta=(11,5,125)$ | $x^2-22x+125$ | $11.18034$ | $=\lvert\lambda_1\rvert$ | $\mathbf Q(i)$, 1 cplx | 3 | $1/11.18034$ | $-5.4142$ | FAIL (NF-0) |
| 10 | AZ $\delta=(7,3,81)$ | $x^2-14x+81$ | $9$ | $=\lvert\lambda_1\rvert$ | $\mathbf Q(\sqrt{-2})$, 1 cplx | 3 | $1/9$ | $-5.1972$ | FAIL (NF-0) |
| 11 | AZ$(9,3,-27)$ | $x^2-18x-27$ | $19.39230$ | $1.39230$ | $\mathbf Q(\sqrt3)$ | 3 | $0.718234$ | $-3.3310$ | FAIL |
| 12 | $\mathrm{Sym}^2$ Zagier $\mathbf D$ | $(x^2-11x-1)^2$ | $11.09017$ | $0.09017$ | $\mathbf Q(\sqrt5)$ | 3 | $11.09017$ | $-0.59394$ | FAIL |
| 13 | AESZ 184 | $x^2-88x+2000$ | $44.72136$ | $=\lvert\lambda_1\rvert$ | $\mathbf Q(i)$, 1 cplx | 3 | $1/\sqrt{2000}$ | $-6.8005$ | FAIL (NF-0) |
| 14 | AESZ 207 | $(x-53248)^2(x^2{+}89344x{-}2^{24})$ | $89531.389$ | $53248$ | $\mathbf Q$, $\mathbf Q(\sqrt{17})$ | 4 | $1/53248$ | $-14.8827$ | FAIL |
| — | **Beukers Thm 4** | $x^2-(248{+}110\sqrt5)x+1$ | see §2.1 | — | $\mathbf Q(\zeta_{20})^+$ | 3 | two places | $\mathbf{+0.0111846}$ | **PASS** |

First failure of $\operatorname{lcm}^{k-1}b_n$: $n=2$ (rows 1,2,8,9,10,11,13), $n=3$ (row 3),
$n=5$ (rows 5,7,12,14), $n=11$ (row 6), $n=13$ (row 4).  All $k$'s agree with
`HERFURTNER_CLASSIFICATION.md` §6.1 and `CDT_FINDER.md` §6; **$k=4$ for AESZ 207
is measured here for the first time**.

**Reading.**  Only four rows in the entire corpus pass, and three of them are the
already-known rational ones (Apéry $\zeta(2)$, Apéry $\zeta(3)$, Beukers'
$\sqrt{\text{Apéry}}$ = his own Theorem 3).  **The unique genuinely number-field
instance is Beukers' Theorem 4 row, and it passes only because of the second
real place.**

### 2.3 Two convention points that bite

**(a) $\rho$ for a repeated root (row 12).**  With a repeated $\lambda_2$ the two
readings of "$\lambda_2$" differ.  The correct one is the *second-nearest
singular point* (the fold cancels a singular point, whatever its multiplicity),
and this was checked numerically: the genuine small form
$q_n=[t^n](\alpha F-G)^2$ with $\alpha=\zeta(2)/5$ exactly has
$\log|q_n/q_{n-1}|\to-\log\varphi^5=-2.406059$ (measured $-2.41422,-2.41025,
-2.40889,-2.40820$ at $n=100,200,300,400$), so $\rho=\varphi^5$ and the margin is
$-0.5939$, not $-5.4061$.  **[measured]**

> **Caveat, and it matters.**  A *two-term* form $\beta A_n-B_n$ inside the
> $\mathrm{Sym}^2$ recurrence does **not** reach that rate: its measured growth is
> $+\log\varphi^5$, i.e. only the dominant root.  Reaching $\varphi^{-5n}$ needs
> the genuine quadratic $\alpha^2F^2-2\alpha FG+G^2$, whose coefficients are
> *dependent*.  So the $-0.5939$ is the optimistic figure and no linear form
> realises it.

**(b) AESZ 207's second conjugate root is positive.**  $-2^7(349-85\sqrt{17})=+187.389\ldots$
(since $85\sqrt{17}=350.463>349$), and the product of the conjugate pair is
$-2^{24}$.  $|N(2^7(349+85\sqrt{17}))|^{1/2}=2^{12}=4096$ **[exact]**.

### 2.4 AESZ 207: the $\mathbf Q(\sqrt{17})$ conjugate pair cannot be separated

This is the one row where a genuine two-place structure is visible without the
coefficients leaving $\mathbf Q$, so it deserved a test.  **[verified]**

1. The subdominant solution of growth $187.389^n$ **exists**: backward iteration
   ($n=520\to0$ at 2600 digits) converges to one direction from three independent
   seeds, agreeing to $\approx1250$ digits, with normalised initial vector
   $(1,\;11.96532007154343270608,\;713.25388423109355112515,\;64609.37448842018931453)$;
   forward iteration gives $\log|u_n/u_{n-1}|\to\log187.38920672=5.233188$.
2. Its coefficients are **not** in $\mathbf Q(\sqrt{17})$ nor in $\mathbf Q$:
   `lindep([1,√17,u_j/u_0])` at 200 digits returns heights $10^{65}$–$10^{70}$;
   over $\mathbf Q$, $10^{95}$.  The subdominant line is a **connection quantity**,
   transcendental in exactly the way $\zeta(3)$ is in Apéry's construction.
3. Every $\mathbf Q$-rational solution grows at the dominant rate
   $\log89531.389=11.402345$ (all four coordinate vectors, indistinguishably).
4. A two-term Apéry form $\alpha A_n-B_n$ ($\alpha=-5.359831688\cdot10^{-6}$) caps
   out at $\log53248=10.882716$ exactly, confirming $\rho=1/53248$.

> **Conclusion.**  The asymptotic solution basis at $n=\infty$ is defined over
> $\mathbf Q(\sqrt{17})$ and Galois swaps the two conifold lines, so no
> Galois-equivariant construction can isolate one of them; what survives is the
> norm, $|N|^{1/2}=2^{12}=4096$, and $\log4096=12\log2$ **is exactly** the
> measured $2$-adic resource $\sigma_2\log2$, $\sigma_2=12$ (`AESZ207.md` §0.3).
> But the obstruction to *using* the small line is stronger than Galois: the
> connection map from $n=0$ data to it is transcendental, so it is unreachable by
> an algebraic choice of coefficients over **any** number field.  The
> counterfactual two-place margin $\log4096-4=+4.3178$ is therefore not a
> Diophantine quantity; it is the product formula in disguise.

---

## 3. The two new Herfurtner rows: uniformisation, and the complex folds

### 3.1 Both rows uniformised exactly (new)

*(script `lattice/number_field/nf_herfurtner_uniformisation.gp`, log alongside.)*

`HERFURTNER_CLASSIFICATION.md` §6.2 produced the two rows and identified their
fibre configurations and monodromy groups; it did not produce the uniformising
maps.  Both are now explicit.

> **Row #30** ($I_1I_7\,II\,II$, monodromy $\Gamma_0(7)$), 
> $(n+1)^2u_{n+1}=(117n^2+78n+21)u_n-441(3n-1)^2u_{n-1}$, $u_1=21$:
> $$\boxed{\ t=-\frac19\Bigl(\frac{\eta(7\tau)}{\eta(\tau)}\Bigr)^4,\qquad
> A(t)\,\bigl(1-117t+3969t^2\bigr)^{1/3}=\theta_{-7}(\tau),\ }$$
> where $\theta_{-7}=1+2\sum_{n\ge1}\bigl(\sum_{d\mid n}\chi_{-7}(d)\bigr)q^n$ is
> the weight-one theta series of $\mathbf Q(\sqrt{-7})$ (the unique form in
> $M_1(\Gamma_0(7),\chi_{-7})$).  Equivalently
> $$q\frac{dt}{dq}\Big/t=\theta_{-7}^{\,2}\in M_2(\Gamma_0(7)).$$
> **[verified to $q^{90}$, exact rational/series arithmetic]**

> **Row #45** ($I_3\,III\,III\,III$, the index-three level-three congruence
> group), $(n+1)^2u_{n+1}=(72n^2+36n+6)u_n-108(4n-1)(4n-3)u_{n-1}$, $u_1=6$:
> $$\boxed{\ j=\frac{(24t-1)^3}{8t^3}\quad\Longleftrightarrow\quad
> t=\frac{-1}{2\bigl(j^{1/3}-12\bigr)},\qquad j^{1/3}=\frac{E_4}{\eta^8},\ }$$
> and
> $$3\,q\frac{dt}{dq}\Big/t=A(t)^2\bigl(1-72t+1728t^2\bigr)^{1/2},$$
> a *meromorphic* weight-two form (pole at the elliptic point $\tau=i$, which is
> the $III$ fibre placed at $t=\infty$); $t$ is a series in $w=q^{1/3}$,
> $t=-\tfrac12w-6w^2-72w^3-740w^4-\cdots$.
> **[verified to $w^{45}$]**

The two behave differently, and the difference is exactly the exponent:

* at #30 the gauge exponent is $1/3$ and the two $II$ fibres are elliptic points
  of order $3$, so $R(t)^{1/3}$ is single-valued on $\mathbf H$ and
  $A\cdot R^{1/3}$ is a genuine holomorphic modular form;
* at #45 the gauge exponent is $1/4$ while the three $III$ fibres are elliptic
  points of order $2$, so $R(t)^{1/4}$ is **not** single-valued: the row's "$F$"
  is a square root of a weight-two object, i.e. the quadratic twist allowed by
  `HERFURTNER_CLASSIFICATION.md` Theorem H2, and there is no weight-one form.

**Consequence for Theorem B\*.**  Neither row is a modular Apéry row in the
Eisenstein-source sense: the Eichler ratio $\Theta=B(t)/A(t)$ has
$\Phi:=D^2\Theta$ with unbounded $3$-power denominators
($\Phi=-\tfrac19q-\tfrac4{27}q^2-\tfrac{17}{81}q^3+\tfrac{364}{729}q^4-\cdots$
for #30) and is **not** in $M_3(\Gamma_0(7),\chi_{-7})$ (a three-dimensional
space: the two Eisenstein orientations plus $\eta(\tau)^3\eta(7\tau)^3$)
**[verified: the $q^1,q^2,q^3$ fit leaves a nonzero residual from $q^4$ on]**.
So `THEOREM_B_EXACT.md`'s formula $\operatorname{Re}\xi=L(\Phi,w+1)$ has **no
source to apply to here**, and the fold constants are genuinely unknown targets
rather than predicted ones.

**Where the folds sit.**  **[exact]**  For #45, with $t_c=(3-i\sqrt3)/144$ one
computes $(24t_c-1)^3=-i\sqrt3/9$ and $8t_c^3=-i\sqrt3/15552$, so
$$j(\tau_*)=1728 :$$
the fold is the CM point $\tau_*\sim i$, **discriminant $-4$** (as it must be —
the $III$ fibres are the order-two elliptic points, all three lying over
$j=1728$).  For #30 the $II$ fibres are the order-three elliptic points of
$\Gamma_0(7)$, lying over $j=0$: the fold is $\tau_*\sim\rho=e^{2\pi i/3}$,
**discriminant $-3$**.  So both complex folds are at the classical CM points, and
a Chowla–Selberg period ($\Gamma(1/4)$ resp. $\Gamma(1/3)$) is a legitimate
candidate in the identification, alongside $L$-values.

### 3.2 The fold constants, to 245 digits

*(script `lattice/number_field/cfold.gp` with drivers `run_all.gp`, `run_v23.gp`,
`run_t.gp`; identification sweeps `hunt3.gp`–`hunt9.gp`; logs
`log_validation.txt`, `log_v23.txt`, `log_targets.txt`, `log_hunt*.txt`,
`xi_targets.txt`.)*

**Method.**  Frobenius local solutions $u,v$ (and an analytic particular
solution $p$ for $L[y_B]=t$) at $t_c$, matched against the exact $y_0,y_B$ at
$t_m=t_c(1-0.4)$, then $\xi=C_-/A_-$; the branch of $s^\rho$ cancels.  Two
safeguards: residuals of $L[u],L[v],L[p]-t$ at $t_m$ (all $\le10^{-250}$) and an
independent re-solve of the connection problem at a second point $t_c(1-0.28)$.

**Validation.**  Zagier $\mathbf B$, $\mathbf E$, $\mathbf D$ all reproduce their
known constants to $\ge10^{-232}$:
$\xi_{\mathbf B}=\overline{\tfrac12L(\chi_{-3},2)+i\tfrac{2\pi^2}{27\sqrt3}}$
(agreement $2\cdot10^{-232}$), $\xi_{\mathbf E}-G/2=1.8\cdot10^{-249}$ with
$\operatorname{Im}\xi_{\mathbf E}=0$ exactly, $\xi_{\mathbf D}-\zeta(2)/5=1.9\cdot10^{-249}$.

> **Correction to the method.**  All three controls have $\rho=0$, not a generic
> exponent: $T(t)\equiv R(t)$ whenever $b=a$ and $e=0$, so the second local
> solution is **logarithmic**, $v=u\log s+w$.  This is exactly §1 of
> `THEOREM_B_EXACT.md` ("$r=2\Rightarrow\rho^2=0$").  Both branches are
> implemented.  The two new rows are genuinely different: $\rho=\tfrac13$ and
> $\tfrac12$, **verified symbolically in $\mathbf Q(\sqrt{-3})$**
> ($t_cR'(t_c)=-(39y+27)/98$, $T(t_c)=(13y+9)/98$, $y=\sqrt{-3}$, for #30).

**The constants.**  **[numerical, $\ge245$ digits; first 60 quoted, full values
in `lattice/number_field/xi_targets.txt`]**

$$\text{#30}:\quad \xi=
0.067867219542878261666469989645637217931570579747852538116276\ldots
-\,0.069078199633222255146707979074349191907178580841245953874001\ldots\,i$$
$$\text{#45}:\quad \xi=
0.164208734159057458471422597020694706583012299953628740569431\ldots
-\,0.236965975165911201071114329980666567666966531886373451808002\ldots\,i$$

with the Eichler value at the fold, $\Theta^*=C_+/A_+=y_B/y_0$ at $t_c$:

$$\text{#30}:\ \Theta^*=0.0302889702955968991649\ldots-0.0323795536843940643341\ldots i,$$
$$\text{#45}:\ \Theta^*=0.0318417514250226942900\ldots-0.0488294605540965137333\ldots i.$$

My independent window-averaged estimator (§5, `est*.gp`, $N=10^5..2\cdot10^5$)
agrees to $\approx10^{-10}$ on both, and the agent's own estimator run
($n=1200..2400$) gives $3.2\cdot10^{-7}$ and $4.3\cdot10^{-7}$; in both cases the
estimator lands on $\xi$ itself (the fold with $\operatorname{Im}t_c<0$), not on
$\bar\xi$.

### 3.3 One exact identification — and it is a CM period, not an $L$-value

> **Result.**  $$\boxed{\ \operatorname{Im}\xi(\#45)
> =-\frac{\Gamma(1/4)^4}{2^{7/2}\,3^{11/4}\,\pi}
> =-\frac{\varpi^2}{\sqrt2\;3^{11/4}}\ }$$
> where $\varpi=\Gamma(1/4)^2/(2\sqrt{2\pi})=2.62205755429211981\ldots$ is the
> lemniscate constant, the Chowla–Selberg period of discriminant $-4$.
> Equivalently $4\cdot3^{11}(\operatorname{Im}\xi)^4=\varpi^8$: `algdep` returns
> the minimal polynomial $2902376448\,x^4-1$ for
> $\operatorname{Im}\xi/(\Gamma(1/4)^4\pi^{-1})$, with
> $2902376448=2^{14}3^{11}$.
> **[VERIFIED to $1.0\cdot10^{-249}$ — every one of the $\approx245$ correct
> digits; coefficient height $2\cdot10^9$, degree 4.]**
> Since $|\lambda|^2=1728$ this reads
> $\operatorname{Im}\xi=-\tfrac2{\sqrt3}\,\Gamma(1/4)^4/(\pi|\lambda|^{3/2})$.
> With the conjugate fold the sign is $+$.

This is exactly what §3.1's CM computation predicted: the fold of #45 is
$\tau_*\sim i$, $j(\tau_*)=1728$, discriminant $-4$ — and the constant is the
disc-$-4$ period, not a Dirichlet $L$-value.  **It is the first constant in this
project's census that is a Chowla–Selberg period rather than a critical
$L$-value**, and it is the direct consequence of the source $\Phi$ not being
modular (§3.1).

### 3.4 The negatives, and why they are structural

**The Theorem-B\* pattern fails outright.**  $\operatorname{Im}\xi(\#30)$ is
*not* $q\,\pi^2/\sqrt7$ and $\operatorname{Im}\xi(\#45)$ is not $q\,\pi^2/\sqrt3$
for any rational $q$ (the near miss
$\operatorname{Im}\xi(\#30)\cdot54\sqrt7/\pi^2=0.99996327\ldots$ is definitively
$\ne1$ at 245 digits); $\operatorname{Re}\xi$ is not a rational multiple of
$L(\chi_{-7},2)$, $L(\chi_{-3},2)$ or $\zeta(2)$.  Searched, all at 130–200
digits working precision against 245-digit targets, on
$\operatorname{Re}\xi,\operatorname{Im}\xi,\operatorname{Re}\Theta^*,
\operatorname{Im}\Theta^*$ for both rows:

1. **Newform/Dirichlet basis**, 2- and 3-term `lindep`, with $1$ adjoined:
   $\{\pi^2,\pi^3,\pi^4,\pi^2/\sqrt7,\pi^3/\sqrt7,\pi^2/\sqrt3,\pi^3/\sqrt3,
   L(\chi_{-7},2),L(\chi_{-3},2),\zeta(3)\}$ together with $L(f,1),L(f,2)$ for
   $f=$ 7.3.b.a $=\eta(\tau)^3\eta(7\tau)^3$ and $f=$ 12.3.c.a, and
   $\sqrt7,\sqrt3,\sqrt{21},\pi$ multiples.  **No hit** (heights $<10^6$,
   residual $<10^{-85}$).
2. **Chowla–Selberg monomial scans**: $\Gamma(1/3)^a\pi^b3^{c/2}7^{d/2}$ (308
   monomials), $\Gamma(1/4)^a\pi^b2^{c/2}3^{d/2}$ (396), and
   $\Omega^j\pi^k\times\{1,\sqrt2,\sqrt3,\sqrt6,\sqrt7,\sqrt{21},3^{1/3},3^{2/3},
   7^{1/3},7^{2/3},2^{1/4},3^{1/4},6^{1/4},L\text{-values}\}$ (1188 / 1404).
   **No hit besides the one above.**
3. **`algdep` with algebraic coefficients** (strictly stronger than `bestappr`):
   target $/(\Gamma(1/3)^a\pi^b M)$ for $a\in[-18,18]$, $b\in[-8,8]$;
   $\Gamma(1/4)^a\pi^b$, $a\in[-20,20]$; and $A_7^{a/2}\pi^b$ with
   $A_7=\prod\Gamma(j/7)^{\chi_{-7}(j)}$ the disc-$-7$ Chowla–Selberg period;
   $M$ ranging over ten $L$-values and logs; degrees $\le18$, heights $<10^7$.
   **The only hit in the entire sweep is $\operatorname{Im}\xi(\#45)$.**
4. **3-term $\mathbf Q$-`lindep` over 560 monomials**
   $\Gamma(1/4)^a\pi^b2^{c/4}3^{d/4}$ (#45) and 441 monomials
   $\Gamma(1/3)^a\pi^b3^{c/3}7^{d/3}$ (#30) — $156\,000$ resp. $97\,000$ pairs.
   **No hit** (height $<10^7$, residual $<10^{-100}$).
5. **Ratios** $\operatorname{Re}\xi/\operatorname{Im}\xi$,
   $\operatorname{Re}\Theta^*/\operatorname{Im}\Theta^*$, $\xi/\Theta^*$,
   $\operatorname{Re}\xi/\operatorname{Re}\Theta^*$: not algebraic of degree
   $\le26$ with height $<10^8$ (the apparent hits at degree 24–26 sit exactly at
   the LLL noise floor).

> **A correction worth recording.**  There is **no** newform in
> $S_3(\Gamma_0(9),\chi_{-3})$ — PARI gives dimension $0$, and $\eta(3\tau)^6$ is
> not a modular form there ($q$-order $3/4$).  The genuine weight-3 CM newform
> for $\mathbf Q(\sqrt{-3})$ is **12.3.c.a**, and that is what was used.  Also
> $L(f_7,1)/L(f_7,2)=\sqrt7/(2\pi)$ exactly ($\varepsilon=+1$), so those two are
> not independent targets.

**Why the $L$-value basis cannot work — and one new exact identity.**  The
integrating factor of the inhomogeneous equation is $\mu=t\,R(t)^\kappa$ with
$\kappa=(2d+e)/(2d)$, and the residual $-b+a\kappa$ vanishes identically for both
rows: $\kappa=\tfrac23$ for #30, $\tfrac12$ for #45 — the same twist as the
$R^{-1/3}$, $R^{-1/2}$ gauge in $y_0$ (§3.1).  Hence for #30, with
$\theta=\theta_{-7}$ and $y_0=\theta R^{-1/3}$,
$$\Phi=t\,\theta^3R(t)^{-2/3},$$
which is holomorphic on $\mathbf H$ ($R$ has a *triple* zero in the local
$\tau$-coordinate at an order-3 elliptic point) but acquires $t^{-1/3}$ at the
cusp $0$ — so $\Phi\notin M_3(\Gamma_0(7),\chi_{-7})$, confirming §3.1's fit
failure from the other side.  Its natural home is the cyclic **cubic** cover of
$X_0(7)$ branched at $t_c,\bar t_c,\infty$, of **genus one with $j=0$** (CM by
$\mathbf Q(\sqrt{-3})$).  But the *un*-twisted combination is a genuine cusp
form, and $S_3(\Gamma_0(7),\chi_{-7})$ is one-dimensional:

> $$\boxed{\ \frac{t\,\theta_{-7}(\tau)^3}{1-117t+3969t^2}=-\tfrac19\,\eta(\tau)^3\eta(7\tau)^3,
> \qquad t=-\tfrac19\Bigl(\frac{\eta(7\tau)}{\eta(\tau)}\Bigr)^4\ }$$
> — the newform **7.3.b.a**.  **[verified independently to $q^{80}$, exact
> series arithmetic]**

Consequently the *bilinear* combination $C_-A_+-A_-C_+$ **is** a period of
7.3.b.a, from the cusp $\infty$ to the elliptic point:
$$\int_0^{t_c}y_0(x)R(x)^{-2/3}dx=2\pi i\int_{i\infty}^{\tau_*}\Psi\,d\tau
=\frac{t_c}3\bigl(3969(t_c-\bar t_c)\bigr)^{2/3}\bigl(C_-A_+-A_-C_+\bigr).$$
What resists identification is $\xi=C_-/A_-$, which needs $A_+$ and $A_-$
*separately*: $A_+=W(\tau_*)$ is a CM value, but $A_-\propto W'(\tau_*)$ is a
**Maass–Shimura derivative** of a weight-one form, so
$A_-\in\overline{\mathbf Q}\,\Omega^3+\overline{\mathbf Q}\,\Omega/(4\pi y_*)$ —
a quasi-period term no monomial basis in $\Gamma$-values reaches.

> **This is very likely why #45's imaginary part came out clean and #30's did
> not**: #45's fold is $\tau_*=i$, so $y_*=\operatorname{Im}\tau_*=1$ and the
> quasi-period correction is rational; #30's fold is an order-3 point of
> $\Gamma_0(7)$ with $y_*$ irrational.  The next basis to try is
> $\Omega_{-3}^k/(\pi\operatorname{Im}\tau_*)^m$ for #30, and for
> $\operatorname{Re}\xi(\#45)$ the set $\{\varpi^2,\pi/\varpi^2,\pi\}$ **with the
> $1/(4\pi)$ Maass correction adjoined as a separate basis vector**.



---

## 4. Shopping list: what a new instance would have to look like

### 4.1 The three levers

Theorem NF-1 has exactly three inputs, and only three:

$$\frac1g\sum_v d_v\log\rho_v\;>\;k+\log D .$$

1. **$k$** — the sharp Eichler denominator exponent, $=w+1$ unless a free
   integration drops it (`ROOT_ROWS.md` Theorem R3).  $k=2$ costs $e^2=7.39$ per
   place, $k=3$ costs $20.09$, $k=4$ costs $54.60$.  This is the dominant term
   and the reason nothing above weight one survives over $\mathbf Q$.
2. **$\rho_{v_0}=1/|\lambda_2^{(v_0)}|$ at the folded place** — the project's
   score.  Over $\mathbf Q$ the complete answer is known: by
   `HERFURTNER_CLASSIFICATION.md` Corollary H6 the *only* second-order rows with
   $\log(1/|\lambda_2|)>k\ge2$ are Zagier $\mathbf D$ and Beukers'
   $\sqrt{\text{Apéry}}$, over the whole Kodaira-admissible world.
3. **the other places** — the number-field lever, worth
   $\frac1g\sum_{v\ne v_0}d_v\log\rho_v$.  For Beukers' row this is
   $\tfrac12\log0.83519=-0.0900$: the second place *costs* $0.09$ and is worth
   having only because the first place gains $\tfrac12\log493.965=3.101$ instead
   of the $\log$ of a rational row's radius.  **The number field is not a
   free lunch; it halves the good place's contribution and adds a bad one.**
   It pays only when the folded place is spectacularly good, which happens
   exactly when $\lambda_2^{(v_0)}$ is a very small unit.

### 4.2 Where Beukers' row comes from, and how to look for more

The mechanism is visible in Beukers §3: on $\Gamma_1(5)$ the hauptmodul $y$ takes
the values $-\tfrac{11}2\pm\tfrac52\sqrt5$ at the cusps $0$ and $\tfrac12$, so
the four cusps fall into a Galois orbit of size two over $K=\mathbf Q(\sqrt5)$;
the two coordinates $t,s$ built from them are conjugate, and so are their
$q$-coefficients.  The source $F$ carries $\pm5\sqrt5F_4(\chi,\tau)$ and is
therefore genuinely over $K$.  Hence the recipe:

> **Recipe.**  A genus-zero modular curve (or Atkin–Lehner quotient) whose
> cusps — or whose singular fibres — form a Galois orbit of size $g\ge2$ over a
> real field $K$, together with an Eisenstein/cuspidal source defined over $K$
> and *not* over $\mathbf Q$.

Rows produced this way are invisible to every earlier scan in the project, all
of which work over $\mathbf Z$; §4.2a closes that gap.  The two families to
sweep are the $(R3)$ shape
$$(n+1)^3u_{n+1}=(2n+1)(An^2+An+B)u_n-Cn^3u_{n-1},\qquad A,B,C\in\mathcal O_K,$$
(Beukers: $K=\mathbf Q(\sqrt5)$, $A=124+55\sqrt5$, $B=34+15\sqrt5$, $C=1$) and
its $(R2)$ analogue.  Since Beukers' $C=1$ and his $\lambda$ is a unit of
$\mathbf Q(\zeta_{20})^+$, the place to look first is $|N(C)|=1$.

### 4.2a The scan: Beukers' row is unique

*(scripts `lattice/number_field/rq_*`: `rq_scan.c` (C + OpenMP, exact modular
$\mathcal O_K$ arithmetic), `rq_cgen*.gp`, `rq_driver*.sh`, with every survivor
re-verified exactly in PARI by `rq_verify.gp` / `rq_beukers.gp` / `rq_detail.gp`.)*

Both shapes were scanned over all $60$ squarefree $D$ with $2\le D\le97$, writing
$\alpha=a_0+a_1\omega$, $\omega=\tfrac{1+\sqrt D}2$ or $\sqrt D$.  The integrality
test is division-free — $v_n=(n!)^ku_n$ satisfies
$v_{n+1}=(2n{+}1)(An^2{+}An{+}B)v_n-Cn^6v_{n-1}$ (R3) resp.
$v_{n+1}=(An^2{+}An{+}B)v_n-Cn^4v_{n-1}$ (R2), and $u_n\in\mathcal O_K$ iff
$(n!)^k\mid v_n$ — sieved modularly to $n=34$, survivors re-checked exactly to
$n=200$.  **Zero false positives at the exact re-check in every pass.**

| pass | $C$ restricted to | box | tuples tested | integral |
|---|---|---|---|---|
| R3 unit | $\pm\varepsilon^m$, $|m|\le6$; $|a|\le400$, $|b|\le60$ | | $1.4654\cdot10^{13}$ | $1288$ |
| R2 unit | same | | $1.4654\cdot10^{13}$ | $1804$ |
| R3 small norm | $1\le|N(C)|\le16$, $|c|\le20$; $|a|\le200$, $|b|\le60$ | | $6.0458\cdot10^{12}$ | $1972$ |
| R2 small norm | same | | $6.0458\cdot10^{12}$ | $2620$ |
| R3 rational $C$ (validation) | $C\in\mathbf Z$, $|C|\le130$; $|a|\le60$, $|b|\le40$ | | $1.4985\cdot10^{12}$ | $5022$ |
| R2 rational $C$ | same | | $1.4985\cdot10^{12}$ | $7988$ |
| R3 norm-$V$ | $C=\pm1$; $|a|\le3000$, $|b|\le600$, margin prefilter | | $1.5205\cdot10^{13}$ | $\mathbf 4$ (Beukers only) |
| R2 norm-$V$ | $C=\pm1$; $|a|\le2000$, $|b|\le200$ | | $5.5779\cdot10^{12}$ | $0$ |
| R3/R2 norm-$W$ | $C=\pm\varepsilon$; $|a|\le500/400$ | | $4.0172\cdot10^{12}$ | $0$ |
| R3/R2 norm-$T$ | $|N(C)|\le16$ mod $\varepsilon^2$; $|a|\le500/400$ | | $2.3257\cdot10^{12}$ | $4/0$ |

Completed passes: $7.152\cdot10^{13}$ tuples; with two further wide runs aborted
after producing nothing new, **grand total $\approx4.26\cdot10^{14}$ tuples**.
After deduplication: **R3 — $3020$ integral classes, $750$ non-degenerate, $\mathbf 5$
with Mode-I $>0$, $105$ with Mode-II $>0$; R2 — $4997$ integral classes, $1199$
non-degenerate, $\mathbf 0$ with Mode-I $>0$, $107$ with Mode-II $>0$.**  Every
non-degenerate row has sharp $k=3$ (R3) / $k=2$ (R2), with $k-1$ failing at $n=2$.

> **Result (uniqueness).** **[verified over the box above]**  Exactly one
> $\mathrm{GL}_1$-orbit in the entire scan has a positive Mode-I margin, and it
> is **Beukers' Theorem-4 row**.  The five surviving classes are
> $\pm\varepsilon^m\cdot(\text{Beukers})$, $\varepsilon=\tfrac{1+\sqrt5}2$, in the
> basis $\alpha=a_0+a_1\omega$:
>
> | $D$ | $A$ | $B$ | $C$ | $m$ |
> |---|---|---|---|---|
> | $5$ | $-179+110\omega$ | $-49+30\omega$ | $1$ | $0$ (Beukers, Galois- and sign-flipped) |
> | $5$ | $41+69\omega$ | $11+19\omega$ | $2-\omega=\varepsilon^{-2}$ | $-1$ |
> | $5$ | $-28-41\omega$ | $-8-11\omega$ | $5-3\omega=\varepsilon^{-4}$ | $-2$ |
> | $5$ | $-13-28\omega$ | $-3-8\omega$ | $13-8\omega=\varepsilon^{-6}$ | $-3$ |
> | $5$ | $-110-179\omega$ | $-30-49\omega$ | $1+\omega=\varepsilon^{2}$ | $+1$ |
>
> each reporting Mode-I $=+0.011184591775409465$ and Mode-II
> $=+0.191280989954746221$ — the margins are *exactly* equal along the orbit,
> because both are invariant under the unit rescaling $u_n\mapsto\mu^nu_n$.
> **For the R2 (Zagier) shape the count is zero.**

**How isolated is it.**  The runner-up in the whole R3 scan is $-0.682$ —
Beukers' row rescaled by the *non-unit* $\mu=2\varepsilon^{-1}$, which costs
exactly $\tfrac12\log|N(2)|=0.693$ — and everything else is $\le-3.0$.  So there
is a gap of $0.69$ below Beukers and then a further gap of $2.3$.  The reason is
visible in the template: for unit $C$, Mode-I $>0$ forces
$$\frac{|\lambda_1^{(v_{\rm other})}|}{|\lambda_1^{(v_0)}|}<e^{-6}\approx\frac1{403},$$
i.e. $A$ must be extraordinarily skew between the two places.  Beukers sits at
ratio $0.004115$ against the threshold $2e^{-6}=0.004958$ — **it clears by 20%,
and nothing else in a $4\cdot10^{14}$-tuple box clears at all.**

**Structure of the integral locus.**  Everything integral in the box is a
$K$-scaling $(A,B,C)\mapsto(\mu A,\mu B,\mu^2C)$ of exactly **seven seeds**: the
six rational Almkvist–Zagier–Zudilin $\zeta(3)$-type rows (R3) or Zagier's six
(R2), plus Beukers' one genuinely quadratic row.  **There are no new sporadic
quadratic solutions at all.**  Corollaries: every positive-Mode-I row has $C$ a
unit; every positive-Mode-II row has $|N(C)|\in\{1,4\}$; Mode-I $\le$ Mode-II
always, the two differing by $\tfrac12\log(|\lambda_1|/|\lambda_2|)$ at the
unfolded place, and both dropping by $\tfrac12\log|N(\mu)|$ under a non-unit
rescaling.  Beukers' row has $|N(A)|=251$ and $|N(B)|=31$, both prime.

The Mode-II list is likewise degenerate: only **four** distinct values for R3 —
$0.525494348$ (90 rows: Apéry's $\zeta(3)$ row and its unit twists inside $K$),
$\mathbf{0.191280990}$ (5 rows: **the Beukers orbit, the only genuinely quadratic
entries**), $0.178920758$ (10 rows: Apéry scaled by $|N(\mu)|=2$) — and three for
R2, headed by $0.406059125$ (97 rows: Zagier $\mathbf D$ and its twists).

**Validation.**  The rational-$C$ pass independently re-derives the complete
classical lists and nothing else: R3 $(17,5,1)$, $(11,5,125)$, $(7,3,81)$,
$(9,3,-27)$, $(10,4,64)$, $(12,4,16)$; R2 $(7,2,-8)$, $(9,3,27)$, $(10,3,9)$,
$(11,3,-1)$, $(12,4,32)$, $(17,6,72)$.  Nothing extra, nothing missing.

**Honest limits of the scan.**

1. **$|N(C)|>16$ with $C$ irrational and non-unit is not covered** — in
   particular there is no quadratic analogue of the $C=125$, $81$, $27$ rows.
   *This is the biggest gap.*
2. **$B$ ceiling.**  Mode-I positivity depends only on $(A,C)$, so a
   positive-margin $(A,C)$ whose only integral $B$ has height $>600$ would be
   missed (the unrestricted unit pass has $|b_i|\le60$).
3. The $\varepsilon$- and small-norm classes reach only $|a|\le500$, because
   fields with huge fundamental units ($D=94$: $\varepsilon\approx4.3\cdot10^6$)
   make the prefilter vacuous.  The class $C=\pm1$ where Beukers lives is
   covered to $|a|\le3000$.
4. The margin prefilter keeps $(A,C)$ with
   $-\tfrac12(\log|\lambda_2^{(v_0)}|+\log|\lambda_1^{(v_{\rm other})}|)\ge k-0.35$,
   so it is complete for sharp $k\ge2.65$.  For unit $C$ this is provably
   lossless ($k\le2$ would force $4\mid C$ in $\mathcal O_K$, impossible for a
   unit); for the $|N(C)|\le16$ pass a row with sharp $k=2$ and a margin in a
   narrow window could in principle escape.
5. Only $D$ squarefree $\le97$, and only the shapes (R2), (R3).

Rescaling detection (`ROOT_ROWS.md` Theorem R1) needed no separate
implementation: the scan ranges over all of $\mathcal O_K^3\cap\text{box}$, so any
$K$-rescaling landing back in $\mathcal O_K^3$ is itself scanned — which is
exactly how the five Beukers entries and the $\mu=2\varepsilon^{-1}$ runner-up
appear.

### 4.3 What would be genuinely new

* **Mode I with a new target.**  Any $K$-row whose fold constant is a mixed
  combination $\alpha\zeta(2m+1)+\beta\sqrt D\,L(2m+1,\chi)$ not already covered.
  The pay-off is a statement of Beukers' exact shape — but §4.2a shows there is
  no second example in either Apéry family over any real quadratic field of
  small height, so a new instance must change the *shape*, not the field.  The
  one concrete gap left by the scan is **$C$ irrational and non-unit with
  $|N(C)|>16$** — there is no quadratic analogue of the $C=125,81,27$ rows in
  the corpus, and that is the first place to look.
* **Mode II with a rational target.**  A $K$-row whose *fold constants agree at
  all places* — see §1.3 — would give $\xi\notin\mathbf Q$ under the much weaker
  criterion $\lambda_2^{\rm norm}<e^{-k}$.  No such row is known; the obstruction
  is that a $K$-rational hauptmodul with $t_c\notin\mathbf Q$ forces the fold
  constants apart.  **Deciding whether such a row can exist is the sharpest open
  question this note raises.**
* **Number-field CDT.**  `CDT_FINDER.md` §7 (C2) needs "a clean
  $m\le\mathrm{BC}/(\overline{\log|\varphi'(0)|}-\tau)$ over $K$".  §1.3 here
  says which normalisation such a theorem must use: the *unfolded* radius at
  every place where the target's Galois conjugate is not the analytic fold
  constant.  Transporting CDT's contour with $\lambda_2^{\rm norm}$ (Mode II) at
  a Mode-I conclusion is the error to avoid, and it is worth a factor
  $e^{0.18}$ per index — enough to flip a verdict.

---

## 5. Scripts

All in `lattice/number_field/`, most with a `.log` alongside.  Run with
`gp -q <file>`; every script begins `default(parisizemax,...)`.

**Beukers' Theorem 4 and the accounting.**

| file | what it does |
|---|---|
| `beukers_thm4.gp`, `beukers_thm4b.gp`, `beukers_thm4c.gp` | the row over $\mathbf Z[\varphi]$ (exact, `Mod(y,y^2-y-1)`): re-derives $d_1,d_2,d_3$, integrality to $n\le300$, the characteristic polynomial, both embeddings to 50 digits, the two radii, the Mode-I and Mode-II margins, the companion $c_n$ and sharp $k=3$, and the 200-digit `lindep` giving $\tfrac14(8\zeta(3)-5\sqrt5L(3,\chi_5))$ |
| `beukers_field.gp`, `rootfields.gp` | the field generated by $\lambda$: minimal polynomial $x^4-496x^3+1006x^2-496x+1$, cyclic quartic, `polredabs` $=x^4-5x^2+5$, i.e. $\mathbf Q(\zeta_{20})^+$; $\lambda$ a unit; $\lambda\notin\mathbf Q(\sqrt5)$ |
| `corpus.gp`, `corpus2.gp` | the fourteen rational rows of §2.2: characteristic polynomials, roots, root fields, measured sharp $k$ (with the first failing $n$), radii, scores, margins |
| `twoterm.gp`, `twoterm2.gp` | the $\mathrm{Sym}^2$ convention check of §2.3(a): the genuine quadratic form decays at $\varphi^{-5n}$, a two-term form does not |
| `aesz207.gp`, `aesz207_growth.gp`, `aesz207_growth2.gp` | AESZ 207: exact factorisation of the characteristic polynomial, $k=4$, the backward-iteration construction of the subdominant $187.389^n$ solution, and the `lindep` tests showing its coordinates lie in neither $\mathbf Q$ nor $\mathbf Q(\sqrt{17})$ |

**The complex folds.**

| file | what it does |
|---|---|
| `cfold.gp` | the connection machinery: Frobenius solutions $u,v$ (generic and logarithmic branches) and an analytic particular solution at $t_c$, matched at $t_m$, giving $\xi=C_-/A_-$ and $\Theta^*=C_+/A_+$; residual and second-point safeguards |
| `run_v1.gp`, `run_v23.gp` | the three controls (Zagier $\mathbf B$, $\mathbf E$, $\mathbf D$) |
| `run_all.gp`, `run_t.gp`, `run_t2.gp` | the two targets to $\ge245$ digits; writes `xi_targets.txt` |
| `lvals.gp`, `ident_setup.gp`, `ident.gp`, `ident2.gp` | the $L$-value and newform target library (7.3.b.a, 12.3.c.a) and the first identification sweeps |
| `hunt1.gp`–`hunt9.gp` | the Chowla–Selberg monomial scans, the `algdep`-with-algebraic-coefficient sweeps, the 3-term $\mathbf Q$-`lindep` over 560/441 monomials, and the ratio tests |
| `phi_id.gp` | the check that $\Phi=D^2(B/A)\notin M_3(\Gamma_0(7),\chi_{-7})$ |
| `nf_herfurtner_uniformisation.gp` | the two exact uniformisations of §3.1, verified to $q^{90}$ resp. $w^{45}$, and the cusp-form identity $t\theta_{-7}^3/R=-\tfrac19\eta^3\eta_7^3$ to $q^{80}$ |

**The real-quadratic scan.**

| file | what it does |
|---|---|
| `rq_scan.c` | the exhaustive $(R2)$/$(R3)$ integrality scan over $\mathcal O_K$, compiled to `rq_scan` |
| `rq_cgen.gp`, `rq_cgen2.gp` | generation of the fields and of the admissible $C$ (units and small norms) — `rq_fields.txt` |
| `rq_driver.sh`, `rq_driver2.sh`, `rq_final.sh`, `rq_aggregate.sh` | the scan passes; outputs `rq_out_*.txt`, logs `rq_log_*.txt` |
| `rq_beukers.gp`, `rq_verify.gp`, `rq_detail.gp`, `rq_run_verify.sh` | exact re-verification in PARI of the survivors and of Beukers' own row — `rq_ver_*.txt` |

