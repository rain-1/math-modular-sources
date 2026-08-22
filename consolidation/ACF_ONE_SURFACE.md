# Zagier's $\mathbf A,\mathbf C,\mathbf F$ are one system in three coordinates

*Claude (Opus 5), 2026-08-22. Scripts and logs: `lattice/acf_one_surface/`.
Companion to `HERFURTNER_CLASSIFICATION.md` §3.1/§7 (the three rows on
$I_6I_3I_2I_1$), `RIGIDITY_PROOF.md` (Theorem R and the level-12 cover),
`CONJ_D_PROOF.md` (the Calegari/Coleman input), `EULER_CRITERION.md`
(Theorem F and the $(\psi,\varphi)$ bookkeeping), `SLOPE_CENSUS.md`.*

---

## 0. Verdict first

The task was: prove the $p$-adic Apéry limits of $\mathbf A,\mathbf C,\mathbf F$
from a single **change of cusp** on $X_0(6)$. The outcome splits in two, and the
split is the interesting part.

* **Proved (algebraically, unconditionally).** Every second-order Apéry system
  admits two *cusp moves*, one for each characteristic root — a Möbius change of
  Hauptmodul together with an explicit **linear** gauge factor — and
  $$\boxed{\ t_{\mathbf A}=\frac{t_{\mathbf C}}{1-t_{\mathbf C}},\qquad
    t_{\mathbf F'}=\frac{t_{\mathbf C}}{1-9\,t_{\mathbf C}},\qquad
    F_{\mathbf A}=(1-t_{\mathbf C})F_{\mathbf C},\qquad
    F_{\mathbf F'}=(1-9\,t_{\mathbf C})F_{\mathbf C},\ }$$
  where $\mathbf F'$ is Zagier's $\mathbf F$ in the level-6 normalisation
  $u_n^{\mathbf F'}=(-1)^nu_n^{\mathbf F}$. So $\mathbf A,\mathbf C,\mathbf F$
  are literally **one** Picard–Fuchs system on $X_0(6)$ written in the three
  Hauptmoduls that put, respectively, the $I_2$, the $I_3$ and the $I_6$ fibre at
  $t=\infty$, the $I_1$ (MUM) fibre staying at $t=0$. Theorem 1, §2:
  proved symbolically in `04_cuspmove_proof.py`.

* **The relating map is *not* an Atkin–Lehner involution.** The group
  $\{1,W_2,W_3,W_6\}$ acts *simply transitively* on the four cusps of $X_0(6)$,
  so no nontrivial Atkin–Lehner element fixes the MUM cusp. The cusp move is a
  Möbius transformation of the $t$-line — a change of Hauptmodul, not an
  automorphism of the modular curve. A normaliser element does appear, but only
  to convert the level-6 presentation into Zagier's: $\gamma=\binom{2\ 1}{0\ 2}$,
  i.e. $\tau\mapsto\tau+\tfrac12$, normalises $\Gamma_0(12)$ (not $\Gamma_0(6)$),
  and $t_{\mathbf F}(\tau)=-t_{\mathbf F'}(\tau+\tfrac12)$,
  $\Phi_{\mathbf F}=-\gamma^*\Phi_{\mathbf F'}$.

* **The degree-2 cover of `RIGIDITY_PROOF.md` §3 is an artefact of Zagier's
  level-12 normalisation of $\mathbf F$.** In the level-6 normalisation
  $\mathbf F'$ the relation to $\mathbf C$ is a Möbius substitution, which is an
  *isometric automorphism of the relevant affinoid disc*. This lets one delete
  the Galois-descent step and, with it, the input **H3** ($w$-invariance) that
  `RIGIDITY_PROOF.md` flagged as 3-adic-only and whose complex proof was
  withdrawn. **Theorem 2** (§5) re-proves $\xi_3^{\mathbf F}=\tfrac54\xi_3^{\mathbf C}$
  from **H1** plus the C-side Calegari input alone.

* **The change of cusp does *not* transport anything to $\mathbf A$, and cannot.**
  The row transfers, but the *source* $\Phi=F\,Dt$ is not Möbius-covariant: it
  transforms as $\Phi\mapsto\Phi/(1-\lambda t)$, multiplication by a modular
  unit, and this moves $\Phi$ *between the two weight-3 Eisenstein families*:
  $$\Phi_{\mathbf C}=(1-8V_2)\,\mathcal S,\qquad
    \Phi_{\mathbf F'}=(1+V_2)\,\mathcal S,\qquad
    \Phi_{\mathbf A}=(1-V_2)\,\mathcal T,$$
  with $\mathcal S=E_{3,\chi_{-3},\mathbf 1}$ ($L=L(\chi_{-3},s)\zeta(s-2)$) but
  $\mathcal T=E_{3,\mathbf 1,\chi_{-3}}$ ($L=\zeta(s)L(\chi_{-3},s-2)$).
  Consequently $\mathbf A$ carries the period $\zeta(2)$ rather than
  $L(2,\chi_{-3})$, its alignment prime is $2$ and not $3$
  ($\sigma_p=v_p(d)$, $d_{\mathbf A}=-8$), and there
  $\xi_2^{\mathbf A}=\mathbf 0$ (verified to $2^{2392}$). **One proof covers
  $\mathbf C$ and $\mathbf F$; it does not cover $\mathbf A$, and the obstruction
  is exact, not a gap.**

* **Numerics.** $v_3\bigl(4\xi_3^{\mathbf F}-5\xi_3^{\mathbf C}\bigr)=1586$ at
  $N=800$ (the row's own Cauchy precision is $1581$/$1582$), i.e. the ratio
  $5/4$ holds to $1586$ 3-adic digits — five times the $\ge300$ requested and
  twice the $786$ of `RIGIDITY_PROOF.md` §5.

* **Census by-product.** The `SLOPE_CENSUS.md` entry "Zagier $\mathbf F$ vs
  Cooper $s_{18}$: growing but irregular, not a clean linear law with $r=5/4$"
  is a false negative of the test statistic used there. The clean test gives
  $v_3\bigl(4\xi^{\mathbf F}-5\xi^{s_{18}}\bigr)=589$ at $N=600$ — exactly the
  Cauchy precision of $s_{18}$ ($\sigma_3=1$, so $\approx N$, not $2N$). See §7.

---

## 1. Notation

Rows are written
$$(n+1)^2u_{n+1}=P(n)u_n-Q(n)u_{n-1},\qquad P(n)=a(n^2+n)+b,\quad Q(n)=d\,n^2,$$
$u_0=1$, $u_{-1}=0$; the triple is $(a,b,d)$. Zagier's six are
$\mathbf A(7,2,-8)$, $\mathbf B(9,3,27)$, $\mathbf C(10,3,9)$,
$\mathbf D(11,3,-1)$, $\mathbf E(12,4,32)$, $\mathbf F(17,6,72)$.
Characteristic roots $\lambda,\mu$: the roots of $\lambda^2-a\lambda+d$.
The Picard–Fuchs operator is $L=\theta^2-tP(\theta)+t^2Q(\theta+1)$,
equivalently
$$t^2R(t)\,y''+t\,S(t)\,y'+(-bt+dt^2)\,y=0,\qquad
R=1-at+dt^2,\quad S=1-2at+3dt^2. \tag{1.1}$$

Modular data: $q=e^{2\pi i\tau}$, $D=q\,d/dq$, $(V_kf)(\tau)=f(k\tau)$,
$\Phi_X=F_X\,Dt_X$, $\Theta_X=D^{-(w+1)}\Phi_X=D^{-2}\Phi_X$ ($w=1$),
$A^X(t_X)=F_X=\sum a_nt_X^n$, $B^X(t_X)=F_X\Theta_X=\sum b_nt_X^n$,
$R^X_\xi=B^X-\xi A^X$, $\xi_p^X=\lim_nb_n/a_n$ in $\mathbb Q_p$.

$$\mathcal S=E_{3,\chi_{-3},\mathbf 1}=\sum_{m\ge1}\Bigl(\sum_{k\mid m}\chi_{-3}(m/k)k^2\Bigr)q^m,
\qquad
\mathcal T=E_{3,\mathbf 1,\chi_{-3}}=c_0+\sum_{m\ge1}\Bigl(\sum_{k\mid m}\chi_{-3}(k)k^2\Bigr)q^m,$$
$\mathcal E=D^{-2}\mathcal S$, $\mathcal E_{\mathcal T}=D^{-2}\mathcal T$
(zero constant term in both cases),
$\Psi_0=(1-4V_2)\mathcal S$ (the class with vanishing Mellin factor,
`RIGIDITY_PROOF.md` §6), $\kappa=\tfrac12\zeta_3(2)$.

Level-6 eta data (`CONJ_D_PROOF.md` Table 1):
$$t_{\mathbf C}=\frac{\eta_1^4\eta_6^8}{\eta_2^8\eta_3^4},\qquad
F_{\mathbf C}=\frac{\eta_2^6\eta_3}{\eta_1^3\eta_6^2}.$$

---

## 2. Theorem 1: the cusp move

> **Theorem 1 (cusp move).** Let $(a,b,d)$ be a row as above with characteristic
> roots $\lambda,\mu$ ($a=\lambda+\mu$, $d=\lambda\mu$), let $y$ solve (1.1), and put
> $$s=\frac{t}{1-\lambda t}\quad\Bigl(\text{so } t=\frac{s}{1+\lambda s}\Bigr),
> \qquad z(s)=(1-\lambda t)\,y(t)=\frac{y(t)}{1+\lambda s}.$$
> Then $z$ solves the equation (1.1) of the row
> $$\boxed{\ (a,b,d)\ \longmapsto\ (a^\sharp,b^\sharp,d^\sharp)
>   =\bigl(\mu-2\lambda,\ \ b-\lambda,\ \ \lambda^2-\lambda\mu\bigr),\ }$$
> whose characteristic roots are $-\lambda$ and $\mu-\lambda$. Moreover
> $$t^\sharp=\frac{t}{1-\lambda t},\qquad
>   F^\sharp=(1-\lambda t)\,F,\qquad
>   \Phi^\sharp=F^\sharp\,Dt^\sharp=\frac{\Phi}{1-\lambda t}. \tag{2.1}$$

*Proof.* The substitution and gauge are an identity in the Weyl algebra; carried
out symbolically in `04_cuspmove_proof.py` (sympy), which returns
$\text{orig}/\text{new}=1$ and difference $0$ for generic $(\lambda,\mu,b)$.
The last clause of (2.1) is one line: $Dt^\sharp=Dt/(1-\lambda t)^2$, so
$F^\sharp Dt^\sharp=(1-\lambda t)F\cdot Dt/(1-\lambda t)^2=\Phi/(1-\lambda t)$. $\square$

**Reading.** The four singular points of (1.1) are $0$, $t_1=1/\lambda$,
$t_2=1/\mu$, $\infty$; $t\mapsto t/(1-\lambda t)$ is the unique Möbius map fixing
$0$ with derivative $1$ there that sends $t_1$ to $\infty$. So Theorem 1 is
exactly "**move the cusp $t_1$ to $t=\infty$, keep the MUM point at $t=0$ with
$t=q+O(q^2)$**". The normalisation $t=q+O(q^2)$ is what forces the multiplier to
be $1$ and the gauge to be *linear*.

**Consequences for Zagier's six** (`07_orbits.gp`, exact):

| row | $(a,b,d)$ | $\lambda,\mu$ | cusp moves |
|---|---|---|---|
| $\mathbf A$ | $(7,2,-8)$ | $-1,8$ | $(10,3,9)=\mathbf C$, $(-17,-6,72)=\mathbf F'$ |
| $\mathbf C$ | $(10,3,9)$ | $1,9$ | $(7,2,-8)=\mathbf A$, $(-17,-6,72)=\mathbf F'$ |
| $\mathbf F$ | $(17,6,72)$ | $8,9$ | $(-7,-2,-8)$, $(-10,-3,9)$ |
| $\mathbf E$ | $(12,4,32)$ | $4,8$ | $(0,0,-16)$ (degenerate), $(-12,-4,32)=\mathbf E$ |
| $\mathbf B$ | $(9,3,27)$ | $\tfrac{9\pm3\sqrt{-3}}2$ | not defined over $\mathbb Q$ |
| $\mathbf D$ | $(11,3,-1)$ | $\tfrac{11\pm5\sqrt5}2$ | not defined over $\mathbb Q$ |

Here $(a,b,d)\sim(-a,-b,d)$ under Zagier's own rescaling $t\mapsto t/\mu$ with
$\mu=-1$, i.e. $u_n\mapsto(-1)^nu_n$; we write $\mathbf F'$ for
$(-17,-6,72)$, so $u_n^{\mathbf F'}=(-1)^nu_n^{\mathbf F}$.

> **Corollary 1.1 (orbits).** Under the cusp move, Zagier's six rows fall into
> the orbits $\{\mathbf A,\mathbf C,\mathbf F\}$, $\{\mathbf B\}$,
> $\{\mathbf D\}$, $\{\mathbf E\}$. The move is defined over $\mathbb Q$ exactly
> when $a^2-4d$ is a square, which happens only for $\mathbf A,\mathbf C,\mathbf E,\mathbf F$.

This is the recurrence-level form of `HERFURTNER_CLASSIFICATION.md` §3.1
("$\mathbf A,\mathbf C,\mathbf F$ are three different rows on the same surface
$I_6I_3I_2I_1$ … differing only in which of the three remaining fibres is placed
at $\infty$"), and it makes that statement an identity rather than a
cross-ratio coincidence.

**Remark (the $\mathcal I=0$ entry of the census).** The third placement for
$\mathbf E$ is $(a,b,d)=(0,0,-16)$, $\mathcal I=a^2/d=0$. This row is *not*
non-existent: it is $u_{2m}=\binom{2m}{m}^2$, $u_{2m+1}=0$ — integral, but
degenerate: $P\equiv0$, and the characteristic roots are $\lambda,\mu=\pm4$, so
$|\lambda|=|\mu|$ and there is **no archimedean Apéry limit**. It is the
classical central-binomial-squared (Legendre) row in the variable $t^2$.
`HERFURTNER_CLASSIFICATION.md` §3.1's "no integral row exists there" should read
"no *non-degenerate* integral row"; verified in `07_orbits.gp`.

---

## 3. The modular picture on $X_0(6)$

The cusps of $\Gamma_0(6)$ are $\infty,1/2,1/3,0$ of widths $1,3,2,6$, carrying
the fibres $I_1,I_3,I_2,I_6$. Ligozat's formula (computed in
`06_units_and_level12.gp`) gives the divisors

| function | eta quotient | $\operatorname{ord}_0$ | $\operatorname{ord}_{1/2}$ | $\operatorname{ord}_{1/3}$ | $\operatorname{ord}_\infty$ |
|---|---|---|---|---|---|
| $t_{\mathbf C}$ | $\eta_1^4\eta_6^8/(\eta_2^8\eta_3^4)$ | $0$ | $-1$ | $0$ | $1$ |
| $1-t_{\mathbf C}$ | $\eta_1\eta_3^5/(\eta_2^5\eta_6)$ | $0$ | $-1$ | $1$ | $0$ |
| $1-9t_{\mathbf C}$ | $\eta_1^9\eta_6^3/(\eta_2^9\eta_3^3)$ | $1$ | $-1$ | $0$ | $0$ |
| $F_{\mathbf C}$ | $\eta_2^6\eta_3/(\eta_1^3\eta_6^2)$ | $0$ | $1$ | $0$ | $0$ |

(all four identities **EXACT** to $O(q^{260})$). So $t_{\mathbf C}$ has its zero
at the $I_1$ cusp $\infty$ and its pole at the $I_3$ cusp $1/2$; it takes the
value $1$ at the $I_2$ cusp $1/3$ and $1/9$ at the $I_6$ cusp $0$. Hence

$$t_{\mathbf A}=\frac{t_{\mathbf C}}{1-t_{\mathbf C}}\ \text{ puts }I_2\text{ at }\infty,\qquad
  t_{\mathbf F'}=\frac{t_{\mathbf C}}{1-9t_{\mathbf C}}\ \text{ puts }I_6\text{ at }\infty,$$

matching the invariants $\mathcal I=a^2/d=-\tfrac{49}8,\ \tfrac{100}9,\ \tfrac{289}{72}$
of `HERFURTNER_CLASSIFICATION.md` §3.1 for $I_2,I_3,I_6$ at $\infty$.

Because $\Gamma_0(6)$ has index $12$ and no elliptic points, a weight-$1$ form
has divisor degree $1$: $F_{\mathbf C}$ is the unique weight-$1$ form vanishing
at the cusp put at $t_{\mathbf C}=\infty$, and by Theorem 1 the same is true for
$F_{\mathbf A}=(1-t_{\mathbf C})F_{\mathbf C}$ (vanishes at $1/3$) and
$F_{\mathbf F'}=(1-9t_{\mathbf C})F_{\mathbf C}$ (vanishes at $0$).
Likewise, a weight-$3$ form has divisor degree $3$ and one computes from the
table that $\Phi_{\mathbf C}$ has **simple zeros at the three cusps other than
$1/2$ and no zero at $1/2$** — and again the same statement holds for
$\Phi_{\mathbf A}$ (no zero at $1/3$) and $\Phi_{\mathbf F'}$ (no zero at $0$).

**Verified $q$-expansion identities** (`01_cusp_move.gp`, `02_sources.gp`):

$$t_{\mathbf A}=\frac{t_{\mathbf C}}{1-t_{\mathbf C}},\quad
  t_{\mathbf F'}=\frac{t_{\mathbf C}}{1-9t_{\mathbf C}},\quad
  F_{\mathbf A}=(1-t_{\mathbf C})F_{\mathbf C},\quad
  F_{\mathbf F'}=(1-9t_{\mathbf C})F_{\mathbf C}\qquad\textbf{EXACT to }O(q^{241})$$

$$t_{\mathbf F'}(\tau)=-t_{\mathbf F}(\tau+\tfrac12),\qquad
  F_{\mathbf F'}(\tau)=F_{\mathbf F}(\tau+\tfrac12)\qquad\textbf{EXACT to }O(q^{241})$$

$$\Phi_{\mathbf A}=\frac{\Phi_{\mathbf C}}{1-t_{\mathbf C}},\qquad
  \Phi_{\mathbf F'}=\frac{\Phi_{\mathbf C}}{1-9t_{\mathbf C}},\qquad
  \Phi_{\mathbf F}(\tau)=-\Phi_{\mathbf F'}(\tau+\tfrac12)\qquad\textbf{EXACT to }O(q^{260})$$

**Why $\tau\mapsto\tau+\tfrac12$ and not an Atkin–Lehner element.**
$\begin{pmatrix}1&1/2\\0&1\end{pmatrix}$ normalises $\Gamma_0(N)$ iff $4\mid N$;
so it normalises $\Gamma_0(12)$ but not $\Gamma_0(6)$. This is precisely why
$t_{\mathbf F}$ is a level-12 eta quotient while $t_{\mathbf F'}$ is level 6, and
why `RIGIDITY_PROOF.md`, which works with Zagier's $t_{\mathbf F}$, sees a
degree-2 cover. As a cross-check, `06_units_and_level12.gp` re-verifies that
document's cover identities from the present eta data:
$$v^2-(9u^2-8u+1)v+u^2=0,\qquad t_{\mathbf F}=\frac{v-u}{9v-1},\qquad
u=t_{\mathbf C}(\tau),\ v=t_{\mathbf C}(2\tau)\qquad\textbf{EXACT to }O(q^{260}),$$
together with $t_{\mathbf F'}=-\bigl(\tfrac{v-u}{9v-1}\bigr)\!\mid_{\tau\mapsto\tau+1/2}$.
The Atkin–Lehner group $\{1,W_2,W_3,W_6\}$ of $X_0(6)$ acts simply transitively
on $\{\infty,1/2,1/3,0\}$ (each nontrivial $W$ is fixed-point-free on cusps), so
**no** Atkin–Lehner element fixes the MUM cusp: the A/C/F relation cannot be an
Atkin–Lehner involution, and is not one.

---

## 4. The sources: Eisenstein bookkeeping

All **EXACT** to $O(q^{260})$ (`02_sources.gp`):

$$\Phi_{\mathbf C}=(1-8V_2)\mathcal S,\qquad
  \Phi_{\mathbf F'}=(1+V_2)\mathcal S,\qquad
  \Phi_{\mathbf F}=(1+V_2)(1-8V_2)\mathcal S=(1-7V_2-8V_4)\mathcal S,$$
$$\Phi_{\mathbf A}=(1-V_2)\mathcal T .$$

Equivalently, as modular-function identities on $X_0(6)$ (also EXACT):
$$(1-9t_{\mathbf C})\,(1+V_2)\mathcal S=(1-8V_2)\mathcal S,\qquad
  (1-t_{\mathbf C})\,(1-V_2)\mathcal T=(1-8V_2)\mathcal S. \tag{4.1}$$

Since $D^{-2}V_k=k^{-2}V_kD^{-2}$ on zero-constant-term series, with
$\Theta_X=P_X(V_2/4)\mathcal E$ resp. $P_X(V_2/4)\mathcal E_{\mathcal T}$:

| $X$ | $\Phi_X$ | $P_X(Y)$ | $\Theta_X$ | $P_X(\tfrac14)$ | $L(\Phi_X,s)$ |
|---|---|---|---|---|---|
| $\mathbf C$ | $(1-8V_2)\mathcal S$ | $1-8Y$ | $(1-2V_2)\mathcal E$ | $-1$ | $P_{\mathbf C}\,L(\chi_{-3},s)\zeta(s-2)$ |
| $\mathbf F'$ | $(1+V_2)\mathcal S$ | $1+Y$ | $(1+\tfrac14V_2)\mathcal E$ | $\tfrac54$ | $P_{\mathbf F'}L(\chi_{-3},s)\zeta(s-2)$ |
| $\mathbf F$ | $(1+V_2)(1-8V_2)\mathcal S$ | $(1+Y)(1-8Y)$ | $(1+\tfrac14V_2)(1-2V_2)\mathcal E$ | $-\tfrac54$ | $P_{\mathbf F}\,L(\chi_{-3},s)\zeta(s-2)$ |
| $\mathbf A$ | $(1-V_2)\mathcal T$ | $1-Y$ | $(1-\tfrac14V_2)\mathcal E_{\mathcal T}$ | $\tfrac34$ | $P_{\mathbf A}\,\zeta(s)L(\chi_{-3},s-2)$ |
| $\Psi_0$ | $(1-4V_2)\mathcal S$ | $1-4Y$ | $(1-V_2)\mathcal E$ | $0$ | — |

(Here $Y=2^{-s}$ at $s=w+1=2$ is the same number as the $V_2$-eigenvalue
$\tfrac14$, which is why the archimedean and the $p$-adic rational factors below
are literally the same rational number.)

**The structural point.** $\{\mathcal S,V_2\mathcal S\}$ and
$\{\mathcal T,V_2\mathcal T\}$ span two *different* planes of the 4-dimensional
Eisenstein space $\mathcal M_3(\Gamma_0(6),\chi_{-3})$ (`03_companion.gp` §D:
the four series are linearly independent). The cusp move
$\Phi\mapsto\Phi/(1-\lambda t)$ multiplies by a modular unit and is **not**
Hecke-equivariant: from $\mathbf C$ it lands in the $\mathcal S$-plane for
$\lambda=9$ (giving $\mathbf F'$) but in the $\mathcal T$-plane for $\lambda=1$
(giving $\mathbf A$). That single fact is the whole difference between
$\mathbf F$ and $\mathbf A$ in what follows.

Inside the $\mathcal S$-plane, $\Theta$'s obey the exact linear relation
(`02_sources.gp` §D)
$$\boxed{\ \Theta_{\mathbf F'}=-\tfrac54\,\Theta_{\mathbf C}+\tfrac94\,\Theta_{\Psi_0}\ }
\qquad\textbf{EXACT to }O(q^{260}). \tag{4.2}$$

**Archimedean check.** $\xi_\infty^X=-P_X(2)\cdot\tfrac12L(2,\chi_{-3})$ for
$X\in\{\mathbf B,\mathbf C,\mathbf F\}$ gives
$r_{\mathbf C}=\tfrac12$, $r_{\mathbf F}=\tfrac58$
(`THEOREM_B_EXACT.md`), and $\xi_\infty^{\mathbf A}=\zeta(2)/4$ from the
$\mathcal T$-family — a different $L$-value, as the table predicts.

---

## 5. The companion, and Theorem 2

Multiply (4.2) by $F_{\mathbf F'}=(1-9t_{\mathbf C})F_{\mathbf C}$ and write
$\Xi:=F_{\mathbf C}\Theta_{\Psi_0}=\sum_n\psi_nt_{\mathbf C}^n$ for the
$\Psi_0$-companion on the C-curve (`RIGIDITY_PROOF.md` §6):

> **Proposition 2 (transfer identity, unconditional).** As formal power series,
> $$A^{\mathbf F'}=(1-9t_{\mathbf C})\,A^{\mathbf C},\qquad
>   B^{\mathbf F'}=(1-9t_{\mathbf C})\Bigl(-\tfrac54B^{\mathbf C}+\tfrac94\,\Xi\Bigr),$$
> hence for every $\xi$
> $$\boxed{\ B^{\mathbf F'}+\tfrac54\xi\,A^{\mathbf F'}
>   =(1-9t_{\mathbf C})\Bigl(-\tfrac54\,R^{\mathbf C}_\xi+\tfrac94\,\Xi\Bigr),
>   \qquad t_{\mathbf C}=\frac{t_{\mathbf F'}}{1+9t_{\mathbf F'}}.\ } \tag{5.1}$$

**VERIFIED EXACT** to $O(q^{260})$ in `03_companion.gp` §A; the $t$-expansions of
$A^{\mathbf F'},B^{\mathbf F'},A^{\mathbf C},B^{\mathbf C},A^{\mathbf A},B^{\mathbf A}$
extracted from the modular construction agree with the recurrence solutions
($u_0=1,u_1=b$ and $u_0=0,u_1=1$) to $n=60$ with residual $0$, and
$a_n^{\mathbf F'}=(-1)^na_n^{\mathbf F}$, $b_n^{\mathbf F'}=-(-1)^nb_n^{\mathbf F}$
exactly, so $\xi_p^{\mathbf F'}=-\xi_p^{\mathbf F}$.

Compare (5.1) with `RIGIDITY_PROOF.md` (4.1): there the right-hand side lives on
the level-12 cover and involves *both branches* $u,v$ of a degree-2 extension,
which is what forced Galois descent and the input **H3**. Here it lives on the
same curve, and the change of variable is a Möbius substitution.

> **Lemma 3 (the Möbius disc automorphism).** Fix $0<\delta<2$ and let
> $\mathbb D_\delta=\{|x|_3\le3^\delta\}$ with Tate algebra $T_\delta$ and
> sup-norm $\|\cdot\|$. Then $x\mapsto x/(1+9x)$ is an isometric automorphism of
> $\mathbb D_\delta$, and $1\pm9x$ and $(1\pm 9x)^{-1}$ are units of $T_\delta$ of
> norm $1$.
>
> *Proof.* $|9x|_3\le3^{\delta-2}<1$ on $\mathbb D_\delta$, so $|1+9x|_3=1$ and
> $|x/(1+9x)|_3=|x|_3$. $\square$

> **Theorem 2 (Conjecture D for $(\mathbf C,\mathbf F)$ at $p=3$, Möbius form).**
> Assume
> * **H1** $\ v_3(a_n^{\mathbf F})=O(\log n)$ (measured: $a_n\in\mathbb Z$ and
>   $\max_{n\le400}v_3(a_n^{\mathbf F})=11$, `RIGIDITY_PROOF.md` `radius.gp`);
> * **(ii)** $\xi:=\xi_3^{\mathbf C}$ exists and $R^{\mathbf C}_\xi$ converges on
>   $|t_{\mathbf C}|_3\le3^\delta$ for some $0<\delta<2$;
> * **(iii)** $\Xi$ converges on the same disc.
>
> Then $\ \xi_3^{\mathbf F'}=-\tfrac54\xi_3^{\mathbf C}$ and
> $\ \boxed{\xi_3^{\mathbf F}=\tfrac54\,\xi_3^{\mathbf C}}$.
>
> *Proof.* By Lemma 3 the substitution $t_{\mathbf C}=t_{\mathbf F'}/(1+9t_{\mathbf F'})$
> carries $T_\delta$ (in $t_{\mathbf F'}$) to $T_\delta$ (in $t_{\mathbf C}$)
> isometrically, so by (ii)+(iii) both $R^{\mathbf C}_\xi$ and $\Xi$, re-expanded
> in $t_{\mathbf F'}$, lie in $T_\delta$; and $1-9t_{\mathbf C}$ is a unit of
> $T_\delta$ of norm 1. Hence the right-hand side of (5.1) lies in $T_\delta$,
> i.e. $B^{\mathbf F'}-\bigl(-\tfrac54\xi\bigr)A^{\mathbf F'}$ converges on a disc
> of radius $3^\delta>1$. By **H1** and the uniqueness clause of Lemma 1 of
> `RIGIDITY_PROOF.md` §2 (equivalently `CONJ_D_PROOF.md` §1),
> $\xi_3^{\mathbf F'}=-\tfrac54\xi$. Finally
> $\xi_3^{\mathbf F}=-\xi_3^{\mathbf F'}$. $\square$

**Where (ii) and (iii) come from.** Both are the same Calegari/Coleman input
already isolated in `CONJ_D_PROOF.md` §3–§4: with
$E_{\kappa_0}=\kappa+\mathcal E$ the overconvergent weight-$\kappa_0$ Eisenstein
series ($\kappa_0=\langle\cdot\rangle^{-1}$, $\kappa=\tfrac12\zeta_3(2)$),
$$\Theta_{\mathbf C}=(1-2V_2)\mathcal E=(1-2V_2)E_{\kappa_0}+\kappa,\qquad
  \Theta_{\Psi_0}=(1-V_2)\mathcal E=(1-V_2)E_{\kappa_0},$$
because $V_2$ acts as the identity on constants and $P_{\mathbf C}(\tfrac14)=-1$,
$P_{\Psi_0}(\tfrac14)=0$. So $\Theta_{\mathbf C}-\kappa$ and
$\Theta_{\Psi_0}-0$ are overconvergent with zero constant term; (ii) and (iii)
are their $t_{\mathbf C}$-expansions on the ordinary locus. In particular the
general value formula
$$\xi_3\bigl[P(V_2)\mathcal S\bigr]=-P(\tfrac14)\cdot\tfrac12\zeta_3(2)$$
gives $\xi_3^{\mathbf C}=\tfrac12\zeta_3(2)$,
$\xi_3^{\mathbf F'}=-\tfrac54\cdot\tfrac12\zeta_3(2)$,
$\xi_3^{\mathbf F}=\tfrac58\zeta_3(2)$ — the values of `CONJ_D_PROOF.md`.

**What Theorem 2 buys over `RIGIDITY_PROOF.md` Theorem R.** It removes:
(a) the level-12 cover $\mathbb Q(u,v)/\mathbb Q(x)$ and the rank-2 Tate algebra
$\mathcal O=T_\delta[T]/(T^2-e_1T+e_2)$;
(b) Proposition 3 there ($w$-invariance) and hence the input **H3**, whose
archimedean proof was withdrawn in that document's review note and which
survived only as a 3-adic numerical check;
(c) the discriminant discussion $(1-8x)(1-6x)^2$.
It keeps **H1** and the C-side overconvergence. The mechanism is unchanged —
$\tfrac54=P_{\mathbf F}(\tfrac14)/P_{\mathbf C}(\tfrac14)$ is still the Mellin
eigenvalue ratio — but the *transfer of coordinates*, which was the hard step,
becomes Lemma 3.

**The $P(\tfrac14)=0$ class, re-measured** (`03_companion.gp` §C): with
$\psi_n=[t_{\mathbf C}^n]\Xi$,
$$v_3(\psi_n/a_n^{\mathbf C})=16,\,34,\,54,\,71,\,88,\,113\quad\text{at }n=10,20,\dots,60,$$
i.e. $2n-O(\log n)$: the $\Psi_0$-companion has $3$-adic Apéry limit exactly $0$
at the full slope, as `RIGIDITY_PROOF.md` §6 found.

---

## 6. Why $\mathbf A$ is not covered — and cannot be

> **Theorem 3.** The cusp move $\mathbf C\to\mathbf A$ does **not** transport the
> $p$-adic Apéry limit, for three independent and exact reasons.
>
> **(a) Different Eisenstein family.** $\Phi_{\mathbf A}=(1-V_2)\mathcal T$ while
> $\Phi_{\mathbf C},\Phi_{\mathbf F'}\in\langle\mathcal S,V_2\mathcal S\rangle$.
> Consequently no identity of the shape (4.2) exists: the three series
> $\Theta_{\mathbf C},\Theta_{\Psi_0},\Theta_{\mathbf A}$ have rank $3$ (verified,
> `03_companion.gp` §D), so $\Theta_{\mathbf A}\notin\langle\Theta_{\mathbf C},\Theta_{\Psi_0}\rangle$
> and there is no transfer identity (5.1) to run.
> The archimedean shadow is that $\mathbf A$'s period is $\zeta(2)/4$, not a
> multiple of $L(2,\chi_{-3})$.
>
> **(b) Different alignment prime.** $\sigma_p=v_p(d)$; $d_{\mathbf A}=-8$ gives
> $\sigma_2=3,\ \sigma_3=0$. Row $\mathbf A$ has **no $3$-adic Apéry limit at
> all**: $v_3\bigl(\xi_N-\xi_{N-1}\bigr)=0$ at $N=800$ (`05_padic.gp` §A).
>
> **(c) At its own prime the limit is $0$.** $\xi_2^{\mathbf A}=0$:
> $v_2(b_n^{\mathbf A}/a_n^{\mathbf A})=287,589,882,1191,1474,1785,2070,2392$ at
> $n=100,\dots,800$, i.e. $3n-O(\log n)$ (`05_padic.gp` §D). The reason is
> Theorem F of `EULER_CRITERION.md` §4.2(3): $P_{\mathbf A}(s)=1-2^{-s}$ *is* the
> Euler factor $\mathcal E_2$, so $Q=1$, and the co-divisor character
> $\varphi=\chi_{-3}$ is odd, so $\Lambda_2^{\mathbf A}=L_2(2,\omega)\equiv0$.
> Meanwhile $\xi_2^{\mathbf F}=\tfrac12L_2(2,\chi_{12})$ is a $2$-adic **unit**
> ($v_2(\xi_2^{\mathbf F})=0$, $v_2(\xi^{\mathbf A}-\xi^{\mathbf F})=0$), so
> $\mathbf A$ and $\mathbf F$ are not aligned at $p=2$ either.

So the honest headline is: **the change of cusp is a genuine identification of
the three rows, but only two of the three sources sit in one Hecke orbit.**
The one proof (Theorem 2) covers $\mathbf C$ and $\mathbf F$ at $p=3$; row
$\mathbf A$ is a *different* extension class attached to the *same* geometry, and
its $3$-adic slot is empty.

---

## 7. Alignment primes, and the census

$\sigma_p=v_p(d)$ (Proposition C with $\kappa_p=0$), computed exactly in
`08_census_check.gp`:

| row | $d$ | $\sigma_2$ | $\sigma_3$ | $\sigma_5$ | family $(\psi,\varphi)$ | $\xi_p$ |
|---|---|---|---|---|---|---|
| $\mathbf A$ | $-8$ | $3$ | $0$ | $0$ | $(\mathbf 1,\chi_{-3})$ | $\xi_2=0$ |
| $\mathbf B$ | $27$ | $0$ | $3$ | $0$ | $(\chi_{-3},\mathbf 1)$ | $\xi_3=\tfrac12\zeta_3(2)$ |
| $\mathbf C$ | $9$ | $0$ | $2$ | $0$ | $(\chi_{-3},\mathbf 1)$ | $\xi_3=\tfrac12\zeta_3(2)$ |
| $\mathbf D$ | $-1$ | $0$ | $0$ | $0$ | — | none |
| $\mathbf E$ | $32$ | $5$ | $0$ | $0$ | $(\chi_{-4},\mathbf 1)$ | $\xi_2=\tfrac12\zeta_2(2)$ |
| $\mathbf F$ | $72$ | $3$ | $2$ | $0$ | $(\chi_{-3},\mathbf 1)$ | $\xi_2=\tfrac12L_2(2,\chi_{12})$, $\xi_3=\tfrac58\zeta_3(2)$ |

Within the cusp-move orbit $\{\mathbf A,\mathbf C,\mathbf F\}$ the alignment
primes are therefore $\{2\}$, $\{3\}$, $\{2,3\}$: the orbit meets $p=3$ in
$\{\mathbf C,\mathbf F\}$ (where Theorem 2 applies) and $p=2$ in
$\{\mathbf A,\mathbf F\}$ (where it does not, by Theorem 3(c)).

**Measured cross-alignments** at $N=600$ (`08_census_check.gp`):

| test | $v_3$ | comment |
|---|---|---|
| $\xi^{\mathbf B}-\xi^{\mathbf C}$ | $1185$ | $\approx2N$ |
| $4\xi^{\mathbf F}-5\xi^{\mathbf C}$ | $1188$ | $\approx2N$ — Theorem 2 |
| $\xi^{s_{18}}-\xi^{\mathbf C}$ | $589$ | $\approx N$ ($\sigma_3(s_{18})=1$) |
| $4\xi^{\mathbf F}-5\xi^{s_{18}}$ | $\mathbf{589}$ | $\approx N$ — **new**, see below |

$s_{18}$'s own Cauchy precision at $N=600$ is $587$, so $589$ is *full
precision*. This settles the one open entry of `SLOPE_CENSUS.md` §1: the row
recorded there as "$s_{18}$ vs Zagier $\mathbf F$, $r=5/4$: growing but
irregular ($8$–$16$), not a clean linear law" is a false negative caused by the
test statistic $v_p(r\,a^X_nb^Y_n-r'a^Y_nb^X_n)$, which carries the extra factors
$a^X_na^Y_n$. On the clean statistic $v_3(4\xi^{\mathbf F}-5\xi^{s_{18}})$ the
alignment is perfect, at the smaller of the two slopes. So $\tfrac12\zeta_3(2)$
now has **four** aligned rows at $p=3$: $\mathbf B,\mathbf C,\mathbf F$ and
Cooper's $s_{18}$.

At $p=2$: $v_2(\xi^{\mathbf A})=1785$ at $N=600$ ($3N=1800$), while
$v_2(\xi^{\mathbf E})=v_2(\xi^{\mathbf E}-\xi^{\mathbf F})=-2$ — $\mathbf E$ and
$\mathbf F$ carry the distinct values $\tfrac12\zeta_2(2)$ and
$\tfrac12L_2(2,\chi_{12})$, so the $p=2$ half of the orbit gives **no**
alignment.

---

## 8. Verification log

All scripts in `lattice/acf_one_surface/`; `./run_all.sh` regenerates
`out/*.log`. PARI/GP 2.15.4, $q$-precision $260$, exact rational arithmetic
throughout.

| script | what | result |
|---|---|---|
| `lib.gp` | eta products, rows, $\mathcal S,\mathcal T$, $V_k$, $D^{\pm}$, $q\to-q$, $t$-expansion | — |
| `00_mirror.gp` | mirror maps $t_X(q)$ from the recurrence alone (no modularity) | $t_{\mathbf A}=t_{\mathbf C}/(1-t_{\mathbf C})$ to $O(q^{41})$; $t_{\mathbf F}\ne$ Möbius$(t_{\mathbf C})$ |
| `01_cusp_move.gp` | $t$'s, $F$'s, gauge factors, $\tau\mapsto\tau+\tfrac12$ | all **EXACT** to $O(q^{241})$ (row truncation) |
| `02_sources.gp` | $\Phi_X$ as Eisenstein combinations; (4.1); $\Theta_X$; (4.2) | all **EXACT** to $O(q^{260})$ |
| `03_companion.gp` | transfer identity (5.1); $t$-expansions vs recurrences; $\Xi$'s slope; rank test | (5.1) **EXACT** to $O(q^{260})$; rows match to $n=60$, residual $0$; rank $=3$ |
| `04_cuspmove_proof.py` | symbolic proof of Theorem 1 (sympy) | $\text{orig}/\text{new}=1$, difference $0$ |
| `05_padic.gp` | $3$-adic limits, $N=800$ | $v_3(4\xi^{\mathbf F}-5\xi^{\mathbf C})=\mathbf{1586}$; $v_2(\xi^{\mathbf A})=2392$ |
| `06_units_and_level12.gp` | eta quotients for $1-t_{\mathbf C}$, $1-9t_{\mathbf C}$; Ligozat divisors; level-12 cover cross-check | **EXACT** to $O(q^{260})$ |
| `07_orbits.gp` | cusp-move orbits of the six rows; the $(0,0,-16)$ row | as in §2 |
| `08_census_check.gp` | alignment primes; $p=3$ and $p=2$ clusters; $s_{18}$ | as in §7 |

Key numerics, $N=800$ (`out/05_padic.log`):

```
v_3(xi_C(N)-xi_C(N-1)) = 1581   [2N = 1600]      <- row's own Cauchy precision
v_3(xi_F(N)-xi_F(N-1)) = 1582   [2N = 1600]
v_3(4 xi_F - 5 xi_C)   = 1586                     <- the predicted ratio
v_3(xi_F - xi_C)       = -1     (control)
v_3(3 xi_F - 4 xi_C)   = -1     (control)
xi_3^C mod 3^40 : 3^-1 * 8386265965554334030      <- matches CONJ_D_PROOF §7
xi_3^F mod 3^40 : 3^-1 * 4403999727414453137
(5/4)xi_3^C     : 3^-1 * 4403999727414453137      <- identical
```

So the ratio $4\xi_3^{\mathbf F}=5\xi_3^{\mathbf C}$ is verified to **1586
3-adic digits**, against the requested $\ge300$, and the two controls confirm
that the agreement is specific to the factor $5/4$.

---

## 9. What this does **not** prove

* **It does not prove Conjecture D for $\mathbf C$.** Theorem 2 is a transfer:
  it converts the $\mathbf C$-statement into the $\mathbf F$-statement. The
  $\mathbf C$-statement itself still rests on the Coleman/Calegari
  overconvergence input of `CONJ_D_PROOF.md` §3–§4 (hypotheses (ii)–(iii)
  above), plus **H1** ($\kappa_3=0$), which is measured
  ($\max_{n\le400}v_3(a_n)\le13$) and not proved.

* **It does not cover $\mathbf A$, and the failure is structural, not a gap.**
  Row $\mathbf A$ has $\sigma_3=0$: there is *no* $3$-adic Apéry limit to relate.
  Its own alignment prime is $2$, and there $\xi_2^{\mathbf A}=0$ — the
  archimedean period $\zeta(2)/4$ has $2$-adic avatar $0$ because the Mellin
  polynomial $1-2^{-s}$ is exactly the Euler factor and the co-divisor character
  is odd. So the sentence "one proof of existence covers all three rows" in the
  task statement is **false as stated**; the correct statement is "covers
  $\mathbf C$ and $\mathbf F$". Nothing here says $\mathbf A$ is anomalous —
  its behaviour is exactly what `EULER_CRITERION.md` Theorem F predicts.

* **It does not make the A/C/F relation an Atkin–Lehner relation.** It is a
  change of Hauptmodul (a Möbius map of $\mathbb P^1_t$), which is *not* induced
  by any element of the normaliser of $\Gamma_0(6)$; see §3. The only normaliser
  element in the picture is $\tau\mapsto\tau+\tfrac12$ at level 12, and its only
  role is to convert $\mathbf F'$ into Zagier's $\mathbf F$.

* **It does not settle $\mathbf B$.** $t_{\mathbf B}$ is cubic over
  $\mathbb Q(t_{\mathbf F})$ (`RIGIDITY_PROOF.md` §7); $\mathbf B$'s
  characteristic roots $\tfrac12(9\pm3\sqrt{-3})$ are irrational, so
  $\mathbf B$ is not in the cusp-move orbit at all and gets nothing from
  Theorem 1. The open item (B2) — invariance of $\xi_p$ under a degree-$p$
  change of Apéry parametrisation — is untouched.

* **It does not settle the $p\mid d$ cases.** $\mathbf F$ at $p=2$
  ($\sigma_2=3$) is not reached: (5.1) needs $|9t|_p<1$, which is the statement
  $p=3$, and at $p=2$ the unit $1-9t_{\mathbf C}$ is no longer $\equiv1$ on the
  disc. Domb vs $\varepsilon$ at $p=2$ likewise.

* **The eta-quotient and Eisenstein identities are verified, not certified.**
  Everything marked EXACT holds to $O(q^{260})$ in exact rational arithmetic. A
  Sturm/Ligozat certification (bounds $\le24$ for level 6, $\le48$ for level 12)
  would upgrade them to proofs and is routine; it is not written out here. Only
  Theorem 1 is proved in the strict sense (symbolically, for generic
  $(\lambda,\mu,b)$).

* **The $1586$ digits are numerics.** They confirm $4\xi_3^{\mathbf F}=5\xi_3^{\mathbf C}$
  to that precision; the exact identity is what Theorem 2 delivers, conditionally.

---

## 10. Reproduction

```
cd lattice/acf_one_surface && ./run_all.sh      # writes out/*.log
```
Individual scripts: `timeout 560 gp -q 0X_name.gp </dev/null`, and
`python3 04_cuspmove_proof.py` (needs sympy). Total runtime a few minutes;
`05_padic.gp` at $N=800$ dominates.
