# Linear independence and irrationality measures from the arithmetic holonomy bound

*Claude (Opus 5), 2026-08-22. Scripts and logs: `lattice/holonomy_lindep/`.
Calibrated on and reusing `lattice/cdt_finder/` (`CDT_FINDER.md`),
`lattice/cdt_noncongruence/` (`CDT_NONCONGRUENCE.md`), `lattice/multislope/`
(`MULTISLOPE_PROGRAM.md`), `lattice/mum_survey/` (`MUM_SURVEY.md`).
Literature: CDT, arXiv:2408.15403v2 (`papers/cdt/cdt2/L2chi.tex`) and the ICM
survey `papers/cdt/cdt/icm/ICM.tex` (Theorem `hol_bound`, Proposition
`theomainreprise`, Remark `withintegrationsremark`, Theorem `theotrue`, §6.1).
Tags: **[proved]**, **[verified exact over range R]**, **[measured]**,
**[estimated]**, **[open]**.*

**No irrationality or linear-independence theorem is claimed anywhere in this note.**

---

## 0. Verdict

| claim | verdict |
|---|---|
| A hypothesised relation $a_0+\sum_i a_i\xi_i=0$ among the $r-1$ companion periods of **one** row costs exactly **one** conditional generator $H$, whose $\theta$-orbit has size $=R$, the holonomic rank. So $m=R+1$ and the *linear-independence* statement is bought at the price of the *irrationality* statement | **[proved]** §1 |
| **Fold-regularity of $H$ is a codimension-$\mu$ condition, $\mu=$ multiplicity of the dominant characteristic root.** One rational relation supplies one condition. Hence the whole architecture applies **only to rows with a simple $\lambda_1$** | **[proved]** §2 $+$ **[verified]** three ways: AESZ 207 ($\mu=1$) folds sharply, $\mathrm{Sym}^3(\text{Zagier }\mathbf E)$ ($\mu=3$) and $X_1(5)\,\mathrm{Sym}^2$ ($\mu=2$) do not |
| This **disqualifies every $\mathrm{Sym}^w$ multi-period row in the corpus** — $\beta(4)@24$, $\zeta(7)$ level 24, $X_1(5)\,\mathrm{Sym}^2$ — and 82 of the 227 rank-4 AESZ rows of $z$-degree $\ge3$ | **[verified]** §2, §4 |
| **It also voids `CDT_FINDER.md` §4/§6's scoring of $X_1(5)\,\mathrm{Sym}^2$.** That host has $\lambda_1=\varphi^5$ *double*; its census companion's linear form does **not** decay ($\lvert b_n-\xi a_n\rvert^{1/n}\to 11.09$, not $0.0902$), and $b_n/a_n-\xi\asymp n^{-3/2}$. Its quoted margin $-24.4$ is vacuous, and so is the tempting deficit $+0.34$ it would otherwise have had in architecture (D) | **[measured]** §2.3 — *the single most consequential correction here* |
| **No multi-period row in the corpus has a positive margin, and none can be rescued within `CDT_FINDER.md` §4's inventory ceiling.** The condition $\log C>\tfrac34k$ ($C$ = the hyperbolic/Kodaira ceiling of the post-fold domain), necessary under that cap, fails for every one of them; best headroom $-1.61$ | **[computed]** §5 |
| Best realistic deficits: **AESZ 388** $-3.41$, **AESZ 243** $-3.54$, $\zeta(5)$ level 16 $-4.63$, **AESZ 207** $-12.86$. The only positive multi-period configuration in the corpus is **CDT's own** (Zagier $\mathbf C$ on $\Gamma_0(6)$, architecture (S), $+0.0004$) | **[computed]** §5 |
| **AESZ 207 is the richest object and the worst host**: $\mathrm{prk}_\infty=\mathrm{prk}_2=3$ (maximal), $k=4$ sharp, fold-regular with $h_n=c\cdot53248^n/n^4$ — and $\lambda_2=53248$, so $\mathrm{score}=-14.88$ | **[verified]** §3 |
| Over $\mathbf Q(\sqrt{17})$ the margin is **identical** to the margin over $\mathbf Q$: the row is $\mathbf Q$-rational, so the two archimedean places carry the same $\varphi$ and the `CDT_FINDER.md` §3 per-degree normalisation is a no-op. Linear independence over the larger field is free in this accounting | **[proved]** §3.4 |
| **Task B.** The holonomy route gives $\mu_{\rm eff}\le\mathbf{18.16}$ for Beukers' row against the classical $50.654$ — **better by a factor $2.79$** — and $\mu_{\rm eff}\le\mathbf{10.56}$ for $\zeta(2)$ against Apéry's $11.851$ (better by $1.12$) but against Rhin–Viola's $5.095$ (**worse** by $2.07$). Control: $\zeta(3)$, $12.11$ vs Apéry's $13.418$ vs Rhin–Viola's $5.514$ | **[computed]** §6, calibrated by reproducing CDT's own $24781$ exactly |
| The structural law behind that: the classical measure is $1+\log Q/\mathrm{score}$ and blows up as $\mathrm{score}\to0^+$; the holonomy measure is $\approx 2Em L/\Delta_0$ with $\Delta_0=(m-1)(\mathrm{score}-\theta)$, $\theta=-0.8385$, and stays **finite and bounded on the whole window** $-0.8385<\mathrm{score}\le0$ where the classical argument proves nothing at all | **[proved]** §6.3 |

One sentence: *the holonomy bound makes linear independence exactly as cheap as
irrationality — one conditional generator either way — but it can only see rows
with a simple dominant root and denominator exponent $k\le3$, and every genuine
multi-period row in the corpus has $k\ge3$ together with a $\lambda_2$ far too
large; what the bound is actually worth on the rows it does reach is a better
irrationality measure than Apéry's, by a factor that grows without bound as the
elementary score approaches zero.*

---

## 1. The architecture for a multi-period row

### 1.1 One relation, one generator

Let $L=\sum_{i=0}^{r}z^iP_i(\theta)$ be MUM of holonomic rank $R=\deg_\theta P_0$ at
$z=0$, $P_0=\theta^R$, with coefficient recurrence $\sum_iP_i(n-i)u_{n-i}=0$ of
shift order $r$. Following `MULTISLOPE_PROGRAM.md` §1.1 put $A$ ($a_0=1$) and the
$r-1$ companions $X^{(j)}$ ($x^{(j)}_i=\delta_{ij}$, $i\le j$), so that
$L\bigl(\sum_nx^{(j)}_nz^n\bigr)=j^Rz^j$: each $X^{(j)}$ is the period of an
extension of $\mathbf Q(0)$ by the Picard–Fuchs crystal, with archimedean Apéry
limit $\xi_j=\lim x^{(j)}_n/a_n$.

**Hypothesis to be contradicted.** $a_0+\sum_{j}a_j\xi_j=0$ with
$(a_0,\dots,a_{r-1})\in\mathbf Q^{r}\setminus 0$ (or in $\mathcal O_K^r$).

> **Lemma 1.1 [proved].** Set $H:=\sum_ja_jX^{(j)}+a_0A$. Then
> $h_n=\sum_ja_j\bigl(x^{(j)}_n-\xi_ja_n\bigr)$; $H\in\mathbf Q[\![z]\!]$; $H$ has
> denominator type $[1,\dots,n]^{k}$ with $k=\max_jk_j$; and
> $L(H)=\sum_ja_jj^Rz^j$ is a **polynomial**, so
> $$\{1,\ H,\ \theta H,\ \dots,\ \theta^{R-1}H\}$$
> spans the whole admissible inventory generated by $H$: $\theta^RH$ is
> $\mathbf Q(z)$-dependent on it. Hence $m=R+1$.

*The point.* $m$ does not depend on $r$: **hypothesising a relation among $r$
numbers costs exactly what hypothesising the rationality of one costs.** This is
the sense in which the holonomy bound gets linear independence for free — it is
also why CDT obtain $1,\zeta(2),L(2,\chi_{-3})$ from a single conditional
generator $\mathcal G=\mathrm{Sym}^+(aH_A+bH_B+cH_C)$ (`CDT_UNPACKED.md` §3).

### 1.2 Denominators

All the $\theta^jH$ have the same type as $H$ ($nh_n$ has the same denominator),
and the constant $1$ is denominator-free. So in CDT's array (`CDT_FINDER.md` §1)
there are $r=k$ columns of rate $b_j=1$ with $u_j=1$, no integrations, and

$$\boxed{\ \tau=\tau^\flat=k\Bigl(1-\frac1{m^2}\Bigr),\qquad \tau^\sharp=0 .\ }$$

(For $k=2,m=3$ this is $16/9$, CDT_NONCONGRUENCE's second-order value.)

### 1.3 The domain and the templates

Let $x_1=1/\lambda_1$ be the fold and $\Sigma'$ the set of finite singular points
of $L$ other than $0$ and $x_1$. Under the hypothesis every $\theta^jH$ is
single-valued holomorphic on $\Omega=\mathbf C\setminus\Sigma'$.

* **(K) Koebe slit.** If $\Sigma'$ lies on one ray from $0$, $\Omega$ contains the
  cut plane and $\rho_\Omega=4\min|s|$ exactly; $\mathrm{BC}=\log|\varphi'(0)|$
  (Grunsky). Rigorous end to end.
* **(D) Kodaira multivalent** $\varphi_r(z)=x_2\lambda(rz)$. **Available only when
  $\#\Sigma'=1$**: $\lambda$ omits the single value $1$ on $\mathbf D$, so
  $\varphi_r$ omits $x_2$ and nothing else. With two or more remaining
  singularities $\varphi_r$ *hits* the others and the pullback is not meromorphic.
* **(C) ceiling.** If $\#\Sigma'\ge2$, $\Omega$ is hyperbolic and
  $|\varphi'(0)|\le C:=r_\Omega(0)$ for **every** admissible $\varphi$
  (Schwarz–Pick); we compute $C$ exactly through the modular $\lambda$-function
  (`lattice/cdt_finder/conformal.py`). If $\#\Sigma'=1$, $\Omega$ is not
  hyperbolic and we use Kodaira's ceiling $C=16|x_2|$ instead.
* **(S) symmetrised.** Needs $s=1/\lambda_2\in\mathbf Q$, $\lambda_2\in\mathbf Z$
  *and* $\Sigma'=\{x_2\}$ (otherwise the involution $w=sx/(x-s)$ throws the extra
  points onto new points $\varphi$ must also avoid).

**Readout**, in `CDT_NONCONGRUENCE.md` units:
$$\text{entry}=\log|\varphi'(0)|-\tau,\qquad
\text{margin}=m\cdot\text{entry}-\mathrm{BC},\qquad
\boxed{\text{deficit}=\frac{\text{margin}}{m-1}}$$
(nats per function; equals score $-$ threshold in (K)/(D) and margin$/13$ in (S)).

---

## 2. Fold-regularity is codimension $\mu$ — the gate nearly everything fails

> **Proposition 2.1 [proved].** Let $\lambda_1$ be the dominant characteristic
> root, of multiplicity $\mu$. Near $x_1=1/\lambda_1$ the local solution space of
> $L$ carries a Jordan block of size $\mu$: solutions
> $y_1,\ y_1\log,\ \dots,\ y_1\log^{\mu-1}$ (up to holomorphic corrections).
> A solution has radius $>|x_1|$ iff **all $\mu$** of its local coordinates on
> that block vanish. A hypothesised $\mathbf Q$-linear relation among the
> $\xi_j$ makes exactly the **leading** one vanish. Hence:
> $$\text{one relation}\ \Longrightarrow\ H\ \text{fold-regular}\iff\mu=1 .$$

For $\mu\ge2$ one gets only $x_n^{(j)}/a_n-\xi_j\to0$ *polynomially*, with
$|h_n|^{1/n}\to|\lambda_1|$ — the linear form does not decay at all, and there is
nothing for the holonomy bound (or for Beukers' criterion) to act on.

### 2.1 AESZ 207: $\mu=1$, folds sharply **[verified, $N=6000$]**

$\chi=(x-53248)^2(x^2+89344x-2^{24})$, $\lambda_1=-2^7(349+85\sqrt{17})=-89531.3892$
**simple**. Measured $|x^{(j)}_n-\xi_jA_n|^{1/n}$ at $n=500,1000,2000,3000,5000$:

| $n$ | $B$ | $C$ | $D$ |
|---|---|---|---|
| $500$ | $49429.55$ | $48587.65$ | $47818.48$ |
| $2000$ | $52121.88$ | $51898.50$ | $51691.93$ |
| $5000$ | $52755.97$ | $52665.42$ | $52581.48$ |

against $|A_n|^{1/n}\to89531.39$ ($89165.5$ at $n=5000$). The approach to
$\lambda_2=53248$ from below is *exactly* a power of $n$: $n^4h_n/53248^n$ is
constant to $5$–$6$ digits over $n=500\ldots5000$, $\to-4.34112\cdot10^{-6}$ ($B$),
$+8.07514\cdot10^{-10}$ ($C$), $+2.77484\cdot10^{-13}$ ($D$). The double root
$53248$ produces the $n^{-4}$ and no logarithm. Generic rational combinations
$(1,-3,5)$, $(2,3,-5)$, $(1,1,1)$ show no further cancellation.
**The radius enlarges from $1/89531.39$ to $1/53248$, as required.**

### 2.2 $\beta(4)@24=\mathrm{Sym}^3(\text{Zagier }\mathbf E)$: $\mu=3$, does **not** fold **[verified]**

$\chi=(x-8)^3(x-4)^3$. Measured at $n=10^6$: $|A_n|^{1/n}=7.99993$ and
$|x^{(j)}_n-\xi_jA_n|^{1/n}=7.99992,\,7.99992,\,7.99991,\,7.99990,\,7.99989$ for
$j=1..5$ — i.e. $\to8$, **not** $4$. The linear form is only logarithmically
damped: $\log(n)\,h^{(j)}_n/A_n\to\kappa_j$ with
$\kappa_1=\kappa_2=-4.1249268649$, $\kappa_3=-1.8510099534$,
$\kappa_4=-0.5795032908$, $\kappa_5=-0.1481087332$, so
$h^{(j)}_n\sim\kappa_j\,(24/\pi^3)\,8^n\log n/n$. Folding here needs the $\log^3$,
$\log^2$ and $\log$ coefficients all to vanish (codimension 3); an LLL search over
the six rational sequences finds the shortest integer vector at height
$10^{25.78}$ against a spurious level $10^{27.5}$ — **no rational fold**, not even
a partial one (2 conditions: $10^{17.00}$ vs $10^{18.3}$).

*By-product* **[measured, 55 digits]**: the row is nevertheless genuinely
multi-period, $\mathrm{prk}_\infty=3$, with exactly two $\mathbf Q$-relations
$$2\xi_1-3\xi_2=0,\qquad
16875+281250\,\xi_1-2120000\,\xi_3+4218750\,\xi_4-3538944\,\xi_5=0$$
(residuals $-2.9\cdot10^{-69}$, $5.9\cdot10^{-64}$). The Cauchy criterion gives
only $2$–$4$ digits here ($\mu=3$ ⇒ $1/\log n$ convergence); the digits come from
an exact fit to $u_n=(8^n/n)\sum_kn^{-k}(c_{k,2}L^2+c_{k,1}L+c_{k,0})$, $L=\log n$,
validated against the exact prediction $c^{(A)}_{0,2}=24/\pi^3$ to $68$ digits.

### 2.3 $X_1(5)\,\mathrm{Sym}^2$: $\mu=2$, does **not** fold — a correction to `CDT_FINDER.md`

`lattice/holonomy_lindep/x15a.gp`, `x15b.gp`. The Cauchy square of the Zagier
$\mathbf D$ row, $A_n=1,6,47,408,3745,35598,\ldots$, has minimal recurrence of
shift order $4$ and degree $3$
$$n^3A_n+(-22n^3{+}33n^2{-}23n{+}6)A_{n-1}+(119n^3{-}357n^2{+}366n{-}128)A_{n-2}
+(22n^3{-}99n^2{+}155n{-}84)A_{n-3}+(n-2)^3A_{n-4}=0$$
**[verified exact]**, characteristic polynomial $(x^2-11x-1)^2$, so
$\lambda_1=\varphi^5=11.0902$ and $\lambda_2=-\varphi^{-5}$ are **both double**
(as `CDT_FINDER.md` §6 records). Holonomic rank $R=3$, three companions, $k=3$
sharp for all three **[verified, $n\le300$]**. Measured at $n=100,\dots,600$:

| | $n=100$ | $200$ | $300$ | $400$ | $500$ | $600$ |
|---|---|---|---|---|---|---|
| $\lvert x^{(1)}_n-\xi_1A_n\rvert^{1/n}$ | $10.330$ | $10.643$ | $10.762$ | $10.824$ | $10.861$ | $10.881$ |
| $\lvert A_n\rvert^{1/n}$ | $10.659$ | $10.840$ | $10.910$ | $10.948$ | $10.972$ | $10.988$ |
| $\lvert x^{(1)}_n/A_n-\xi_1\rvert$ | $4.35\cdot10^{-2}$ | $2.57\cdot10^{-2}$ | $1.66\cdot10^{-2}$ | $1.06\cdot10^{-2}$ | $6.22\cdot10^{-3}$ | $2.80\cdot10^{-3}$ |

The linear form **grows like $\lambda_1^n$**; the ratio decays like $n^{-3/2}$; the
Cauchy criterion certifies $\approx3$ digits of $\xi_1$.

> **Correction.** `CDT_FINDER.md` §4 lists "$X_1(5)\,\mathrm{Sym}^2$,
> $\lambda_2^{\rm norm}=1$, entryC $-0.435$, margin $-24.4$" and §6/§7 (C2) builds
> a whole candidate programme on it. Those numbers presuppose an Apéry row with a
> geometrically decaying linear form. **This host has none**: its dominant root is
> double and folding is codimension 2. The correct statement is that the CDT
> architecture does not apply to the $\mathrm{Sym}^2$ row at all, and candidate
> (C2) has to be re-based on whatever rank-3 system Beukers 1987 Theorem 4
> actually uses (his construction is over $\mathbf Q(\sqrt5)$ with conjugate
> linear forms and is not the census companion of this recurrence).
> The loss is real: with $R=3$, $k=3$, $\mathrm{score}=-0.5939$, this host would
> have had deficit $+0.042$ in (K) and $+0.343$ in (D) — the only positive
> multi-period entry anywhere in the corpus.

### 2.4 The scale of the gate

Of the **227** rank-4 AESZ operators with $z$-degree $\ge3$ (i.e. $\ge2$
companions) in the CYCluster dump, **82 are discarded outright** because
$\lambda_1$ has non-unique maximal modulus **[computed]**, and $2$ more for
$|\lambda_2/\lambda_1|>0.97$. A further **85 of the remaining 143 are not
LCM-normalised at all** — their companions carry a geometric $p^{\beta n}$ factor
($\beta\approx1$, $p$ = the prime dividing the leading coefficient of $P_0$), so
by Lemma G of `CDT_NONCONGRUENCE.md` §2 they pay a further $\beta\log p$ and are
inadmissible in CDT's denominator type (6.0.9) as written. **58 survive**, all
with $k\in\{3,4\}$.

---

## 3. AESZ 207 — the one genuine multi-period host

### 3.1 Data **[verified exact unless marked]**

* Operator `4.4.38`, $\theta$-order $R=4$ (all five $P_i$ have $\deg_\theta=4$,
  $P_0=\theta^4$), $z$-degree $4$, so $r-1=3$ companions $B,C,D$.
* $\chi=(x-53248)^2(x^2+89344x-2^{24})$;
  $\lambda_1=-89531.389206720146782$ (simple), $53248$ (double),
  $187.389206720146782$; $|\lambda_2/\lambda_1|=0.5947411346$.
* Singular set $\{0\}\cup\{z_-,\ 1/53248,\ z_+\}\cup\{\infty\}$ with
  $z_\pm=(349\pm85\sqrt{17})/2^{17}$:
  $$\text{fold } z_-=-1.11692671\cdot10^{-5},\qquad
  \Sigma'=\{1.878004808\cdot10^{-5},\ 5.33648665\cdot10^{-3}\}.$$
* $k=4$ **sharp for all three companions** ($k=3$ first fails at $n=5,3,7$),
  $n\le1600$.
* Limits, $N=6000$ exact rationals, Cauchy-certified to $1354/1357/1361$ digits:
  $$\xi_B=-5.04554593441367088625425457970715161172\ldots\cdot10^{-4},$$
  $$\xi_C=+1.50199887161880282091966070356769386787\ldots\cdot10^{-8},\qquad
  \xi_D=-3.02619578142969433687387290057189731270\ldots\cdot10^{-13}.$$
  > **Correction to `MULTISLOPE_PROGRAM.md` §2.2**: the recorded mantissas of
  > $\xi_\infty(C)$ and $\xi_\infty(D)$ are right to $39$/$40$ digits but are
  > printed **without their exponents** $10^{-8}$ and $10^{-13}$. Confirmed by the
  > ratios $\xi_C/\xi_B=-2.97688077988600456\cdot10^{-5}$ and
  > $\xi_D/\xi_B=+5.99775687461135023\cdot10^{-10}$ — the same mantissas as the
  > recorded $-2976.88\ldots$, $+5997.76\ldots$, scaled by $10^{8}$, $10^{13}$.
* $\mathrm{prk}_\infty=3$ **[measured, 1354 digits]**: the LLL relation lattice of
  $\{1,\xi_B,\xi_C,\xi_D\}$ has all four reduced vectors at height
  $10^{337.3}$–$10^{337.8}$ against a spurious level $10^{338.5}$; no $2$- or
  $3$-subset relation; and over $\mathbf Q(\sqrt{17})$ (8-term `lindep`) height
  $10^{136.6}$ against spurious $10^{169.3}$ with residual only $10^{-958}$ —
  **no relation over either field**. $\mathrm{prk}_2=3$ likewise (19129 digits,
  `MULTISLOPE_PROGRAM.md` §2.2).

### 3.2 The conditional generator, the inventory, and $\tau$

$H=a_BB+a_CC+a_DD+a_0A$, fold-regular by §2.1, $k=4$, $\theta$-orbit
$\{H,\theta H,\theta^2H,\theta^3H\}$, so $m=5$ and
$\tau=4(1-1/25)=\mathbf{3.84}$.

**Pure inventory.** $\lambda_2=53248=2^{12}\cdot13\in\mathbf Z$, so the
polylogarithm module $\mathrm{Li}_j(53248\,z)=\sum(53248z)^n/n^j$ lives on
$\mathbf P^1\setminus\{0,1/53248,\infty\}\supset\Omega$ with integer numerators
and $e_i=j$. And that is *all* there is:

> **Proposition 3.1 [proved].** On the AESZ 207 host every unconditional
> $f\in\mathbf Q[\![z]\!]$ admissible for the bound branches only at $1/53248$.
> *Proof.* $z_\pm$ are conjugate over $\mathbf Q(\sqrt{17})$ and the singular
> locus of a minimal operator over $\mathbf Q(z)$ is Galois-stable
> (`CDT_NONCONGRUENCE.md` Proposition G), so any $f$ singular at $z_+$ is
> singular at $z_-$ too, i.e. at the fold — it is then not admissible on
> $\Omega$. $\square$

So the pure module is exactly the polylog orbit at the rational singularity.
Its $p$-adic slopes are $v_2(53248)=12$ and $v_{13}(53248)=1$, giving
$\sum_p\gamma_p\le\log53248=10.883$ in principle and, with CDT's $7{+}7$ split
($\gamma_p=\varsigma_p\log p\,(u-1)^2/m^2$, $u=7$, $m=14$),
$\gamma_2=1.528$, $\gamma_{13}=0.471$, total $\mathbf{2.00}$ **[computed]**.
*This is the largest adelic resource anywhere in the corpus* — and it is worth
$2$ nats against a deficit of $-12.9$.

### 3.3 Geometry and margin

$\#\Sigma'=2$, both on the positive real axis, so
$\Omega\supset\mathbf C\setminus[1/53248,\infty)$:

| template | $\log|\varphi'(0)|$ | $\mathrm{BC}$ | entry | margin | deficit |
|---|---|---|---|---|---|
| (K) Koebe slit, $\rho=4/53248$ | $-9.4964$ | $-9.4964$ | $-13.3364$ | $-57.186$ | $\mathbf{-14.296}$ |
| (C) hyperbolic ceiling $C=3.16827\cdot10^{-4}$ (exact, via $\lambda$) | $\le-8.0572$ | $\ge-8.0572$ | $\le-11.8972$ | $\le-51.429$ | $\le\mathbf{-12.857}$ |
| (D) Kodaira | **unavailable** ($\#\Sigma'=2$) | | | | |
| (S) symmetrised | **unavailable** ($z_+$ is not preserved by $w=sz/(z-s)$) | | | | |

$\mathrm{score}=\log(1/53248)-4=-14.883$; the (K) threshold at $m=5,k=4$ is
$k/m-\log4=-0.586$, so the deficit is $-14.30$, and the ceiling line is a
**rigorous upper bound over every admissible contour**.

### 3.4 The number-field remark, over $\mathbf Q(\sqrt{17})$ **[proved]**

The statement one would like is "$1,\xi_B,\xi_C,\xi_D$ are linearly independent
over $K=\mathbf Q(\sqrt{17})$" (the field of the conifolds; §3.1 excludes a
relation numerically to height $10^{136}$). By `CDT_FINDER.md` §3 the correctly
normalised budget over $K$ is
$$\frac1{[K:\mathbf Q]}\sum_{v\mid\infty}d_v\log|\varphi'(0)|_v
=\log\bigl(\text{const}\cdot|N(s)|^{1/[K:\mathbf Q]}\bigr),$$
with $\tau$ unchanged (the LCMs are rational integers). **Here the operator, the
row and all singular data are $\mathbf Q$-rational** — the pair $z_\pm$ is
Galois-stable as a *set*, and $\lambda_2=53248\in\mathbf Q$ — so both archimedean
places of $K$ see the same $\varphi$, $|N(s)|^{1/2}=|s|$, and

$$\text{margin}_{K}=\text{margin}_{\mathbf Q}=-51.43 .$$

*So on a $\mathbf Q$-rational row, linear independence over a real quadratic field
costs exactly nothing extra in the holonomy accounting.* (This is the honest,
and here vacuous, positive content of the number-field remark: the loss
`CDT_NONCONGRUENCE.md` §7.2 records for a Galois trace applies to rows whose
$\lambda_2$ is itself irrational, not to this one.)

---

## 4. The other candidate hosts

### 4.1 $\zeta(5)$ at level 16

$\lambda_1=2+4\sqrt2=7.65685$ (simple), $\lambda_2=-4$; fold $t=0.13060194$,
$\Sigma'=\{-\tfrac14,\,-0.27345908,\,-\tfrac12\}$ — all on one ray, so
$\rho_\Omega=4\cdot\tfrac14=1$ exactly and the hyperbolic ceiling is
$C=2.18844$ ($\log C=0.7832$). Minimal *joint* operator $(r,D)=(18,11)$, so
$R=11$ and $m=12$, $k=5$ sharp for $B$, $\tau=5(1-1/144)=4.9653$.

**But the row is not multi-period for this purpose.** Its higher companions
$X^{(1)},X^{(2)},X^{(3)}$ admit **no $k\le14$** (`MULTISLOPE_PROGRAM.md` §2.0):
forward iteration divides by a degree-11 $Q_0(n)$ whose prime support is
essentially random. They are genuine solutions with genuine $2$-adic limits, but
they are **not of CDT denominator type** and cannot enter the inventory. The only
admissible companion is $B$, $\xi_B=\tfrac{217}{1024}\zeta(5)$: a single period,
hence an irrationality target only.

entry(K) $=-4.9653$, margin $-59.58$, **deficit $-5.42$**; at the ceiling
$-4.63$. Headroom $\log C-\tfrac34k=0.783-3.75=-2.97$.

### 4.2 $\zeta(7)$: the level-60 programme and its level-12 parent

`ZETA7_LEVEL60.md` establishes there is no level-60 host at all. Its realisable
parent is the level-12 row, $\lambda_1=10$, $\lambda_2=6$, $k=7$ sharp,
$\Sigma'=\{-\tfrac16\}$ (the point $t=1/58$ is apparent). It *does* carry several
companions — one for each purified source $S$ in the 3-dimensional space
$\langle B_0,B_1,B_2\rangle$ — but
$$\xi_S=\tfrac{L(S,7)}{\ldots}\in\mathbf Q\cdot\zeta(7)\quad\text{for every }S,$$
so $\mathrm{span}_{\mathbf Q}\{1,\xi_{B_0},\xi_{B_1},\xi_{B_2}\}=\langle1,\zeta(7)\rangle$
and $\mathrm{prk}=1$: **several companions, one period**. A hypothesised relation
among them is just the irrationality of $\zeta(7)$. $\mathrm{score}=-8.79$,
$\log C=0.981$, headroom $0.981-5.25=-4.27$.

$\zeta(7)$ at level 24 ($\mathrm{Sym}^6$) is excluded by Proposition 2.1: its
characteristic polynomial is $q_5^{\,7}$, every root of multiplicity 7, and
`ROOT_ROWS.md` §5.4 already measured the polynomial ($n^{-1/2}$) convergence.

### 4.3 $\beta(4)$ at level 24

Excluded by Proposition 2.1 (§2.2). Counterfactually it would have had
$\lambda_2=4$, $k=4$, $m=5$, $\log C=\log4$, deficit $-3.41$ at the ceiling and
$-4.36$ with the Kodaira family at $r^*=0.62$ — i.e. *identical* to AESZ 388,
which is the same $(\lambda_2,k)$ with a simple dominant root.

### 4.4 The MUM survey **[computed]**

Sweeping the 297 rank-4 AESZ operators (`lattice/holonomy_lindep/survey_multiperiod.gp`):
227 have $z$-degree $\ge3$; after the multiplicity gate and the LCM gate, 58 are
admissible; **all 31 for which the period rank was determined have
$\mathrm{prk}\ge2$** — multi-periodicity is generic, not rare. Three genuine
$\mathbf Q$-relations were found (AESZ 102 and 304, one 6-term relation each, so
$\mathrm{prk}=6$; AESZ 348, $12\xi_1+1197\xi_2-2432\xi_3+11664\xi_4=0$,
$\mathrm{prk}=3$); in every case $1$ itself is outside the span. **No companion
limit was identified against any $\zeta$- or $L$-value** (rational multiples,
2- and 3-term combinations over a 9-constant library, 300-digit precision,
two-precision genuineness filter): the exclusion is roughly "no relation of
height below $10^{70}$".

---

## 5. The ranking, and why nothing works

$\mathrm{score}=\log(1/|\lambda_2|)-k$; deficit $=$ margin$/(m-1)$;
**headroom** $=\log C-\tfrac34k$ is a *necessary* condition for the architecture
to work with **any** inventory (the margin $(m-1)\log C-m\tau$ is positive for
large $m$ only if $\log C>\tau_\infty$, and $\tau_\infty\ge\tfrac34k$ for the best
inventory `CDT_FINDER.md` §4 scenario B allows, $u_j=m/2$).

| rank | row | $R$ | comp. | $k$ | $\mathrm{prk}_\infty$ | score | $\log\rho_{\rm Koebe}$ | $\log C$ | deficit (K) | deficit (best) | **headroom** |
|---|---|---|---|---|---|---|---|---|---|---|---|
| — | *calib.* Apéry $\zeta(3)$ | 3 | 1 | 3 | 1 | $+0.525$ | $4.912$ | $6.298$ | $+1.162$ | $+1.462$ (D) | $+4.048$ |
| — | *calib.* Apéry $\zeta(2)$ | 2 | 1 | 2 | 1 | $+0.406$ | $3.792$ | $5.179$ | $+1.126$ | $+1.245$ (D) | $+3.679$ |
| — | *calib.* Beukers $(136,10,16,4)$ | 2 | 1 | 2 | 1 | $+0.139$ | $3.525$ | $4.912$ | $+0.859$ | $+0.978$ (D) | $+3.412$ |
| — | *calib.* **CDT**, Zagier $\mathbf C$/$\Gamma_0(6)$ | 2 | 1 | 2 | **2** | $-2.000$ | $1.386$ | $2.773$ | $-1.280$ | $\mathbf{+0.0004}$ **(S)** | $+1.273$ |
| **1** | **AESZ 388** (4.3.18) | 4 | 2 | 4 | 2 | $-5.386$ | $0.000$ | $1.386$ | $-4.800$ | $\mathbf{-3.414}$ (C) | $-1.614$ |
| **2** | **AESZ 243** (4.5.46) | 4 | 4 | 4 | **4** | $-5.099$ | $-1.099$ | $1.256$ | $-5.899$ | $-3.544$ (C) | $-1.744$ |
| 3 | AESZ 313 (4.5.82) | 4 | 4 | 4 | 4 | $-5.792$ | $-1.792$ | $0.359$ | $-6.592$ | $-4.441$ (C) | $-2.641$ |
| 4 | $\zeta(5)$ level 16 | 11 | 1* | 5 | 1 | $-6.386$ | $0.000$ | $0.783$ | $-5.417$ | $-4.634$ (C) | $-2.967$ |
| 5 | AESZ 34 (4.3.1) | 4 | 2 | 4 | 2 | $-6.197$ | $-0.811$ | $0.146$ | $-5.611$ | $-4.654$ (C) | $-2.854$ |
| 6 | AESZ 380 (4.5.111) | 4 | 4 | 4 | 4 | $-6.303$ | $-2.303$ | $-0.430$ | $-7.103$ | $-5.230$ (C) | $-3.430$ |
| 7 | AESZ 104 (4.8.2) | 4 | 7 | 4 | **7** | $-6.197$ | $-2.197$ | $-0.663$ | $-6.997$ | $-5.463$ (C) | $-3.663$ |
| 8 | AESZ 410 (4.3.28) | 4 | 2 | 4 | 2 | $-7.466$ | $-2.079$ | $-0.693$ | $-6.879$ | $-5.493$ (C) | $-3.693$ |
| 9 | AESZ 392 (4.3.22) | 4 | 2 | **3** | 2 | $-7.682$ | $-3.296$ | $-1.910$ | $-6.896$ | $-5.510$ (C) | $-4.160$ |
| 10 | **AESZ 207** (4.4.38) | 4 | 3 | 4 | **3** | $-14.883$ | $-9.496$ | $-8.057$ | $-14.296$ | $-12.857$ (C) | $-11.057$ |
| — | $\zeta(7)$ level-12 parent | ? | 3 | 7 | 1 | $-8.792$ | $-0.405$ | $0.981$ | — | — | $-4.269$ |
| ✗ | $\beta(4)@24=\mathrm{Sym}^3\mathbf E$ | 4 | 5 | 4 | 3 | $-5.386$ | $0.000$ | $1.386$ | — | **no fold ($\mu=3$)** | $-1.614$ |
| ✗ | $X_1(5)\,\mathrm{Sym}^2$ | 3 | 3 | 3 | ? | $-0.594$ | $3.792$ | $5.179$ | — | **no fold ($\mu=2$)** | $+2.929$ |
| ✗ | $\zeta(7)$ level 24 $=\mathrm{Sym}^6$ | 7 | 2 | 7 | ? | $-8.574$ | $-1.575$ | $0.467$ | — | **no fold ($\mu=7$)** | $-4.783$ |

\* only one *admissible* companion; the others are not LCM-normalised.

**Reading.**

1. **No multi-period row has a positive margin.** The best, AESZ 388, is short by
   $3.41$ nats per function ($17.1$ in absolute margin at $m=5$).
2. **None can be rescued within the inventory ceiling this project uses.**
   Headroom is negative for every one of them, best $-1.61$. Headroom $>0$ is
   necessary for *any* contour once the inventory is capped at
   `CDT_FINDER.md` §4 scenario B ($u_j\le m/2$ denominator-free functions,
   already an upper-bound fantasy — CDT themselves have $u=(1,3)$ out of $14$).
   So this is not a search failure: adding pure functions up to that cap, Eichler
   integrals, symmetrisation or the adelic term cannot close a negative headroom.
   (For AESZ 207 the adelic term is worth $2.0$ nats against a headroom of
   $-11.06$.) Only an inventory *beyond* scenario B — asymptotically **all**
   functions denominator-free, which no host is known to supply — could evade the
   test, and then the Bost–Charles numerator becomes the binding constraint
   instead. **[estimated]**
3. **The obstruction is $k$, not multi-periodicity.** When exactly one
   singularity remains, $C=16|x_2|$ (Kodaira) and headroom $>0$ reads
   $$\tfrac34k+\log|\lambda_2|<\log16=2.7726 ,$$
   which with $|\lambda_2|\ge1$ forces $k\le3$; the only multi-period rows with
   $k=3$ (AESZ 392, and the two unnumbered `4.3.10`, `4.3.7`) have
   $|\lambda_2|\ge36$. When several singularities remain, $C$ is not bounded by
   $16|x_2|$ (as AESZ 207 shows: $C=3.168\cdot10^{-4}>16|x_2|=3.005\cdot10^{-4}$,
   because the second puncture is $284$ times further out) and must be computed;
   it was, exactly, for every row in the table, and the headroom is negative in
   every case. **CDT's own host satisfies the condition with room to spare**:
   $k=2$, $\lambda_2=1$, headroom $+1.27$.
4. **Multi-periodicity is not the scarce resource; a small $k$ with a small
   $\lambda_2$ is.** Two mechanisms produce several periods:
   *(i)* several Eisenstein classes on one rank-2 host — $c\le\#\text{cusps}-1$
   (`CDT_NONCONGRUENCE.md` §10.7) — which costs **nothing** in $k$, and is exactly
   what CDT use; *(ii)* higher companions of a higher-order recurrence, which
   requires weight $w\ge2$ hence $k=w+1\ge3$, and in practice $k=4$. Mechanism (i)
   is the only affordable one, and the search for it is the four-cusp/index-12
   search already carried out in `CDT_NONCONGRUENCE.md` §10.7 and
   `CDT_FINDER.md` §5.

**Candidate for audit: none.** No row here has a positive margin, so nothing is
flagged for audit. The nearest thing to a candidate — $X_1(5)\,\mathrm{Sym}^2$,
which *would* have had $+0.34$ — is disqualified by §2.3, and its disqualification
is the one result of this note that changes an existing ledger entry.

---

## 6. Task B — irrationality measures from the holonomy route

### 6.1 The formula, and the calibration

ICM Remark `withintegrationsremark`, eq. `withintegrations` (the quantitative
refinement of Theorem `hol_bound`, a special case of Theorem `theotrue` with
$K=\mathbf Q$, $S=\{v\}$):
$$m\ \le\ \frac{\mathrm{BC}}{\log|\varphi'(0)|-\tau(\mathbf b)-\tau^\sharp(\mathbf e)
-(1-\gamma)\bigl(\tfrac2\kappa-\tfrac{1-\gamma}{\kappa^2}\bigr)\log\tfrac1\rho},$$
where
* $\gamma m=\dim_{\mathbf Q(x)}\mathcal V(\varphi,\vec b)$ is the number of
  **unconditional** ("genuine") functions among the $m$; $E:=1-\gamma$ is the
  mean power-system exponent $\sum_\nu\nu m_\nu/m$;
* $\rho$ is the radius on which the **individual** $A_i,B_i$ converge after
  pullback — i.e. $\varphi(|z|<\rho)$ must sit inside the disc of convergence of
  $A$, which is bounded by the fold.

Writing $\Delta_0:=m(\log|\varphi'(0)|-\tau)-\mathrm{BC}$ for the *qualitative*
margin and $L:=\log(1/\rho)$, the inequality is contradicted as soon as
$\Delta_0\kappa^2-2mLE\kappa+mLE^2>0$, i.e. for
$$\boxed{\ \kappa\ \ge\ \kappa_*:=\frac{E\bigl(mL+\sqrt{mL(mL-\Delta_0)}\bigr)}{\Delta_0}
\ \approx\ \frac{2EmL}{\Delta_0}\ \ (\Delta_0\ll mL),\qquad
\mu_{\rm eff}(\eta)\le\kappa_* .\ }$$

**Calibration [computed, `measure.py`].** CDT's own $L(2,\chi_{-3})$: $m=14$,
$\gamma=\tfrac12$, $\log|\varphi'(0)|=\log(256\cdot0.6292232680)=5.081908$,
$\tau=\tfrac{27}{80}+\tfrac{191}{49}=4.235459$, $\mathrm{BC}=11.845$,
$\Delta_0=0.005288$, $\log(1/\rho)=\log11614=9.359967$. The formula returns
$$\kappa_*=24781.0,$$
i.e. **exactly CDT's printed $24781$** (ICM §6.1). The pipeline is calibrated.

### 6.2 Beukers' row and Apéry's $\zeta(2)$ row

Architecture (D) of `CDT_NONCONGRUENCE.md` (Kodaira $\varphi_r=x_2\lambda(rz)$),
inventory $\{1,H,\theta H\}$, so $m=3$, $\gamma=\tfrac13$, $E=\tfrac23$,
$\tau=16/9$; $\rho$ solves $\max_{|z|=\rho}|\varphi_r(z)|=|x_1|$, i.e.
$1/\rho=16r\,|\lambda_1/\lambda_2|$ to leading order (computed exactly from
$\max_{|w|=s}|\lambda(w)|$). $r$ is re-optimised **for $\kappa$**, which is not
the same $r$ as for the margin.

| row | architecture | $\log|\varphi'(0)|$ | $\mathrm{BC}$ | $\Delta_0$ | $\log\frac1\rho$ | $\kappa_*$ | classical |
|---|---|---|---|---|---|---|---|
| **Beukers** $(136,10,16,4)$ | (K) $m=2$ | $3.5255$ | $3.5255$ | $0.5255$ | $8.4377$ | $31.86$ | |
| | (K) $m=3$ | $3.5255$ | $3.5255$ | $1.7177$ | $8.4377$ | $19.31$ | |
| | **(D) $m=3$, $r^*=0.45$** | $4.1133$ | $5.0554$ | $1.9511$ | $9.0255$ | $\mathbf{18.16}$ | $\mu\le\mathbf{50.654}$ |
| **Apéry $\zeta(2)$** $(11,3,-1,0)$ | (K) $m=2$ | $3.7924$ | $3.7924$ | $0.7924$ | $6.2025$ | $15.40$ | |
| | (K) $m=3$ | $3.7924$ | $3.7924$ | $2.2514$ | $6.2025$ | $10.68$ | |
| | **(D) $m=3$, $r^*=0.40$** | $4.2624$ | $5.0071$ | $2.4467$ | $6.6725$ | $\mathbf{10.56}$ | Apéry $\mathbf{11.851}$, RV $\mathbf{5.095}$ |
| *control* **Apéry $\zeta(3)$** | (K) $m=4$ | $4.9118$ | $4.9118$ | $3.4854$ | $8.4377$ | $14.14$ | |
| | **(D) $m=4$, $r^*=0.48$** | $5.5641$ | $6.6435$ | $4.3630$ | $9.0900$ | $\mathbf{12.11}$ | Apéry $\mathbf{13.418}$, RV $\mathbf{5.514}$ |

### 6.3 Is it better or worse, and why

**Beukers' row: better, by a factor $2.79$** ($18.16$ against $50.654$).
**$\zeta(2)$: better than Apéry's own by $1.12$** ($10.56$ against $11.851$),
**worse than Rhin–Viola by $2.07$** ($10.56$ against $5.095$).
**$\zeta(3)$ (control): the same pattern** — $12.11$ against Apéry's $13.418$
($\times1.11$) and Rhin–Viola's $5.514$ ($\times2.20$).

The mechanism is transparent. Both routes are ratios; write
$\mathrm{score}=\log(1/|\lambda_2|)-k$ and $\log Q=k+\log\lambda_1$. Then
$$\mu_{\rm classical}-1=\frac{\log Q}{\mathrm{score}},\qquad
\kappa_{\rm holonomy}\approx\frac{2EmL}{\Delta_0}
=\frac{2EmL}{(m-1)(\mathrm{score}-\theta)},\quad \theta=-0.8385 .$$
The numerators are comparable in kind ($\log Q\approx7$ for Beukers,
$2EmL\approx36$), but the denominators are not: the classical one is the *bare*
score, the holonomic one is the score measured from the (D) **threshold**
$\theta=-0.83852$ — the $+0.84$ nats that `CDT_NONCONGRUENCE.md` §3.4 prices.
Hence

> **Observation 6.1.** The holonomy irrationality measure is worth most exactly
> where the elementary Apéry argument is weakest. As $\mathrm{score}\to0^+$ the
> classical measure diverges while $\kappa_*$ tends to the finite value
> $2EmL/((m-1)(-\theta))$; and on the entire window
> $$-0.83852<\mathrm{score}\le0$$
> — which contains the two Padé families and the $j=2,3$ rescalings of Apéry's
> $\zeta(2)$ row (`CDT_NONCONGRUENCE.md` §5.4) — the holonomy route supplies a
> **finite effective irrationality measure where the classical argument supplies
> no irrationality proof at all.**

Conversely, when the score is comfortable the two routes give nearly the same
number ($\times1.1$ in favour of holonomy for $\zeta(2)$ and $\zeta(3)$), and
both lose to the permutation-group/Rhin–Viola technology by a factor $\approx2$,
which optimises the *arithmetic* of a specific integral rather than the *count*
of functions. CDT's own caution ("in this form of the result, the extra
coefficient $2$ of $\mu(\eta)^{-1}$ means that, when the usual Apéry limit works,
Proposition `theomainreprise` does not recover the most basic standard bound")
is a statement about their $\gamma=\tfrac12$, $m=14$, $\Delta_0=0.005$
configuration; on a row where $\Delta_0$ is of order $1$ rather than $10^{-3}$
the same formula is competitive, and on a marginal row it wins outright.

**Caveats.** (i) Every measure above inherits the $\mathbf Q(x)$-linear
independence gap of $\{1,H,\theta H,\ldots\}$ (`CDT_NONCONGRUENCE.md` §8, item 1)
— CDT's Lemma 12.1.1 analogue, **[open]**. (ii) $\Delta(r)$ for the Kodaira
family is numerical, $\approx5\cdot10^{-6}$. (iii) $\kappa_*$ is
$\mu_{\rm eff}$ (an exponent with an explicit constant), which is the same notion
as the classical Apéry/Beukers/Rhin–Viola measures quoted, so the comparison is
like for like. (iv) The $r$-optimisation is on a grid of step $0.01$.

---

## 7. Honest ledger

**[proved]** Lemma 1.1 ($m=R+1$ from one relation, independent of the number of
periods); the $\tau$ formula $k(1-1/m^2)$; Proposition 2.1 (fold-regularity is
codimension $\mu$); Proposition 3.1 (the Galois obstruction pins AESZ 207's pure
module to the polylog orbit at $1/53248$); the number-field statement §3.4; the
necessary condition $\tfrac34k+\log|\lambda_2|<\log16$ of §5.3.

**[verified exact]** AESZ 207: $\chi$, $R=4$, $k=4$ sharp for all three
companions ($n\le1600$), the $\mathbf Q$-integrality of $A_n$;
$X_1(5)\,\mathrm{Sym}^2$: the minimal recurrence (order 4, degree 3),
$\chi=(x^2-11x-1)^2$, $k=3$ sharp; $\beta(4)@24$: $\chi=(x-8)^3(x-4)^3$, $k=4$
sharp for all five companions ($n\le1200$); the $z$-degree/multiplicity/LCM census
of the 297 AESZ operators.

**[measured]** All Apéry limits and period ranks: AESZ 207 ($1354$–$1361$ digits,
$\mathrm{prk}_\infty=3$, no relation below $10^{337}$ over $\mathbf Q$ or
$10^{136}$ over $\mathbf Q(\sqrt{17})$); $\beta(4)@24$ ($55$ digits by asymptotic
fit, $\mathrm{prk}_\infty=3$, two explicit relations); the MUM sweep's 31 rows.
All fold-regularity decay rates. The conformal radii (exact through the modular
$\lambda$-function; `conformal.py` is the one already calibrated in
`CDT_FINDER.md`). $\Delta(r)$ for the Kodaira family (reused, plus an extension
to $r\le0.85$; the quadrature degenerates at $r\ge0.86$, so $G(m)$ for $m\ge13$
is a lower bound).

**[estimated, flagged]** The (S)-architecture numbers quoted for comparison
transport CDT's contour loss $0.62922$ and $\mathrm{BC}$ shape term, as in
`CDT_FINDER.md` §8 estimate 1. The adelic $\gamma_p$ for AESZ 207 uses CDT's
$7{+}7$ split. The "best conceivable inventory" $\tau^\flat=\tfrac34\sigma_m$ is
`CDT_FINDER.md` §4 scenario B, an upper-bound fantasy — it is used only to make
the *negative* headroom conclusion stronger.

**[open]** (1) $\mathbf Q(x)$-linear independence of $\{1,H,\theta H,\dots\}$ on
every host — presumed by every margin and every measure here. (2) Whether
$\theta^{R-1}H$ is genuinely independent, i.e. whether $H$ is a cyclic vector for
the module; assumed. (3) Beukers 1987 Theorem 4's actual rank-3 system, now that
the $\mathrm{Sym}^2$ row is known not to be it (§2.3) — this is the concrete next
computation, and it is the only route by which a $k=3$ multi-period host with a
*simple* dominant root might still exist. (4) The identification of AESZ 207's
three periods, and of the 31 MUM rows' companion limits, remains open at height
$10^{70}$–$10^{337}$.

---

## 8. Reproduction

```
lattice/holonomy_lindep/measure.py            # Task B: kappa_* ; calibrates on CDT's 24781
lattice/holonomy_lindep/bound_ml.py           # Task A: tau, entry, margin, deficit, headroom, ranking
lattice/holonomy_lindep/delta_ext.py          # Delta(r) for the Kodaira family, r in [0.80,0.96]
lattice/holonomy_lindep/row_aesz207.gp        # AESZ 207: chi, k sharp, limits (N=6000), prk, fold test
lattice/holonomy_lindep/row_beta4_24.gp       # Sym^3(Zagier E): k sharp, asymptotic fit, prk, fold test
lattice/holonomy_lindep/x15a.gp, x15b.gp      # X_1(5) Sym^2: minimal recurrence, limits, fold test
lattice/holonomy_lindep/survey_multiperiod.gp # the 297-operator MUM sweep (multiplicity/LCM/prk)
lattice/holonomy_lindep/topgeom.gp            # singular sets of the top MUM rows
```
`python3 bound_ml.py` reproduces §5 and `python3 measure.py` reproduces §6
(both need `lattice/cdt_finder/conformal.py` and
`lattice/cdt_noncongruence/delta_table.json`). Logs sit beside the scripts.
