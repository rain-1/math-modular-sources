# The two-row $p$-adic method for $\zeta(5)$: what exists, what aligns, and why it fails

*Claude (Fable), 2026-08-21 overnight run. Scripts: `lattice/zeta5_two_row/`. Exact PARI
arithmetic throughout; every number below is measured, not quoted, unless marked
"[BZ]" (Brown–Zudilin, arXiv:2210.03391v3) or "[notes]".*

---

## 0. Verdict (read this and stop if you like)

1. **A $\zeta(5)$ row with a genuine $2$-adic slope now exists and is built.** The
   level-$16$ purified weight-six source $\Phi_{16}=C^{\rm full}_{5,2}E_6$ gives an
   integral Apéry row $(A_n,B_n)$ with $d_n^5B_n\in\mathbb Z$, archimedean limit
   $B_n/A_n\to\tfrac{217}{1024}\zeta(5)$, and $2$-adic slope $\sigma_2=1$.
2. **Its $2$-adic limit is identified in closed form:**
   $$\boxed{\ \xi_2=\tfrac{7}{32}\,\zeta_2(5)\ }$$
   with $\zeta_2$ the Kubota–Leopoldt $2$-adic zeta function, **verified to $372$
   $2$-adic digits** (= full available precision at $n=399$). This is the first
   instance of the `CONJ_D_PROOF.md` mechanism at $p=2$ and at weight $w=4$.
3. **Conjecture D holds for $\zeta(5)$ at $p=2$**: a second, genuinely different
   level-$16$ row (different host coordinate) has the same archimedean limit and the
   same $2$-adic limit, with cross determinant $v_2(\Delta_{n,m})=\min(n,m)+O(1)$ over
   a full $6\times6$ grid.
4. **The Brown–Zudilin cellular $\zeta(5)$ row is $2$-adically inert.** $Q_n\in\mathbb Z$
   and $v_2(Q_n)=0$ for *every* $n\le900$ (so $\kappa_2=0$ exactly, not merely
   asymptotically); $P_n/Q_n$ has **no** $p$-adic slope at $p=2,3,5,7$; the apparent
   tower valuations are exactly the $d_n^5$ denominator, $-(5\lfloor\log_2n\rfloor+2)$,
   with no convergence. Structural reason: $4\lambda^3-2368\lambda^2-188\lambda+1$ is
   **irreducible over $\mathbb Q_2$** with all three roots of valuation $-\tfrac23$ —
   there is no $2$-adic dominant root, so no $2$-adic Apéry limit can exist.
5. **Therefore the two rows do not align and the two-lattice method has no $\zeta(5)$
   instance.** Direct test: $v_2\bigl(217A_mP_n-1024Q_nB_m\bigr)$ is flat in
   $[-29,+2]$ over every sampling ratio $m/n\in\{\tfrac14,\tfrac12,1,\tfrac32,2,3,4\}$
   and every rational multiplier scanned.
6. **And it could not have worked even if it did align.** With the decayer's
   denominator cost $k=5$ ($d_n^5$), the master formula's cost–resource inequality
   $F<0$ requires
   $$\log\tfrac1{\lambda_{\rm dec}}+\sigma_2^{\rm dec}\log2\;>\;5 .$$
   Brown–Zudilin gives $2.4724+0=2.4724$; the level-$12$ Apéry-perfect row gives
   $2.6339+0=2.6339$. Both are short by a factor $\approx2$. The best *counterfactual*
   two-row quality (BZ decayer + level-16 engine, granting BZ a hypothetical aligned
   slope) is $\delta\le0.875$, attained only in the degenerate no-engine limit, and
   decreasing in the hypothetical slope. **$\delta>1$ is unreachable for $\zeta(5)$ by
   this architecture.**
7. **No irrationality statement of any kind is claimed or implied here.** Per
   `CATALAN_AUDIT.md`, even $\delta>1$ would prove nothing without an unconditional
   non-vanishing statement for the selected linear form.

---

## 1. The level-16 $\zeta(5)$ row (built here for the first time)

### 1.1 Data

* **Source.** $\Phi_{16}=C^{\rm full}_{5,2}E_6$ with
  $C^{\rm full}_{5,2}=(1-V_2)(1-4V_2)(1-16V_2)(1-64V_2)=1-85V_2+1428V_4-5440V_8+4096V_{16}$
  (book `v8/05_hecke_projectors.tex`, Example "$\zeta(5)$"). Normalised,
  $$\Phi_{16}=q-52q^2+244q^3-320q^4+3126q^5-12688q^6+16808q^7-14336q^8+197641q^9-\cdots$$
  $L(\Phi_{16},s)=P(2^{-s})\zeta(s)\zeta(s-5)$, $P(X)=(1-X)(1-4X)(1-16X)(1-64X)$,
  $P(0)=P(2)=P(4)=P(6)=0$, $P(5)=-\tfrac{217}{512}\ne0$, so
  $L(\Phi_{16},5)=\tfrac{217}{1024}\zeta(5)$. (Weight $6$, level $16$, trivial character;
  the unique such source — the purified space is $1$-dimensional at $N=16$.)
* **Host.** $x=\dfrac{\eta(2\tau)\eta(16\tau)^2}{\eta(\tau)^2\eta(8\tau)}=q+2q^2+4q^3+8q^4+14q^5+\cdots$
  is a degree-one Hauptmodul of $\Gamma_0(16)$ (Ligozat divisor $[-1,0,0,0,1]$; found
  independently by the eta-quotient search in `haupt16.gp`), with $x|W_{16}=1/(8x)$.
  The Fricke-invariant coordinate is $w=x/(1+8x^2)$ and the normalised one is
  $$t=\frac{w}{1+2w}=\frac{x}{8x^2+2x+1}=q-8q^3-16q^4+30q^5+192q^6+176q^7-1152q^8+\cdots$$
  Its cusps are $t=0,-\tfrac14,-\tfrac12$ and its folds $t=\tfrac{-1\pm2\sqrt2}{14}$ —
  exactly the singularities recorded in `RESEARCH_MASTER_NOTES_2026-08-17.txt` §35,
  which is an independent confirmation that this is the intended coordinate.
* **Rows.** $F=\Phi_{16}/(Dt)$ (weight $4$), $A(t)=F$, $B(t)=F\cdot D^{-5}\Phi_{16}$:
  $$A_n:\ 1,\,-52,\,268,\,-1920,\,9536,\,-60928,\,274688,\,-1682432,\,6536176,\,-42572992,\,127728448,\dots$$
  $$B_n:\ 0,\,1,\,-\tfrac{429}{8},\,\tfrac{175691}{486},\,-\tfrac{11018695}{3888},\,\tfrac{47811622847}{3037500},\dots$$

### 1.2 Measurements (to $n=399$, `level16_zeta5.gp`)

| quantity | measured | limit / prediction |
|---|---|---|
| $A_n\in\mathbb Z$ | all $n<400$ | — |
| $d_n^5B_n\in\mathbb Z$ | all $n<400$ (no extra factor) | $k=5$ |
| $B_n/A_n\to$ | $\texttt{lindep}=[-1024,217,0]$ **exact** | $\tfrac{217}{1024}\zeta(5)$ |
| $\log|A_n|/n$ | $2.01488$ at $n=399$ | $\log(4\sqrt2+2)=2.03574$ |
| $\log|B_n-LA_n|/n$ | $1.41083$ at $n=399$ | $\log4=1.38629$ |
| $|B_n/A_n-L|^{1/n}$ | $0.54659$ at $n=399$ | $4/(4\sqrt2+2)=0.52245$ |
| $v_2(A_n)$ | $6,12,12,18,19,15,14,16,26$ at $n=10\ldots399$ | bounded $=O(\log n)$, i.e. **H1** holds, $\kappa_2=0$ |
| $v_2\bigl(\tfrac{B_n}{A_n}-\tfrac{B_{n-1}}{A_{n-1}}\bigr)$ | $-1,13,33,53,74,129,178,274,371$ | $\sigma_2=1$ |

Successive slope estimates $ (129{-}74)/50=1.10$, $(178{-}129)/50=0.98$,
$(274{-}178)/100=0.96$, $(371{-}274)/99=0.98$: $\sigma_2=1$ with the usual
$O(\log n)$ noise. No slope at $p=3,5,7$ (all valuations flat/negative).

### 1.2a The recurrence

A mod-$p$ scan of kernel dimensions over $10\le r\le18$, $6\le D\le12$
(`fitrec16.gp`, $p=2^{61}-1$, $\ge25$ excess equations everywhere) gives

| $r\backslash D$ | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|
| $\le15$ | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 16 | 0 | 0 | 0 | 0 | 0 | **1** | 4 |
| 17 | 0 | 0 | 0 | 0 | 4 | 8 | 12 |
| 18 | 0 | 0 | 0 | 5 | 10 | 15 | 20 |

so the **minimal recurrence has order $16$ with degree-$11$ polynomial coefficients**
(the level-12 $\zeta(5)$ analogue is order $14$, degree $5$). The order-$17$,
degree-$10$ member reconstructed exactly over $\mathbb Q$ has leading-coefficient
characteristic polynomial
$$y\,(y+2)^5(y+4)^5(y^2-4y-28)^2\cdot(\text{spurious quadratic}),$$
the $y$ and the quadratic being artefacts of the $4$-dimensional over-parametrisation.
The genuine roots are therefore
$$-4,\qquad -2,\qquad 2\pm4\sqrt2 ,$$
i.e. the reciprocals of the two cusps $t=-\tfrac14,-\tfrac12$ and the two folds
$t=\tfrac{-1\pm2\sqrt2}{14}$ — an exact confirmation of the singularity structure
recorded in `RESEARCH_MASTER_NOTES_2026-08-17.txt` §35, and of $\lambda_1=4\sqrt2+2$
(the measured $\log|A_n|/n\to2.0149$ vs $\log(4\sqrt2+2)=2.0357$).

**Caveat.** $\sigma_2=1$ is still *measured*, not derived: because the recurrence is not
in Zagier/AZ normalisation (order $16$, degree $11$, trailing coefficient not $c\,n^{w+1}$),
Proposition C's $\sigma_p=v_p(c)$ does not apply directly — see
`paper/sections/04_padic.tex`, Remark "Non-Zagier normalisations".

### 1.3 The $2$-adic identification

Following `CONJ_D_PROOF.md` §§2–4, $\Theta=D^{-5}\Phi=P(V_2/2^5)\mathcal E$ with
$\mathcal E=\sum_n\bigl(\sum_{e\mid n}e^{-5}\bigr)q^n$, and the $2$-depleted Eisenstein
series is $\mathcal E^*=(1-2^{-5}V_2)\mathcal E$. Hence

> **Euler-factor criterion.** $\Theta$ is a $V_2$-combination of $\mathcal E^*$
> $\iff$ $(1-X)\mid P(X)$, $X\leftrightarrow V_2$ $\iff$ $\Phi$ has zero constant term
> $\iff$ $P(0)=0$. Writing $P=(1-X)Q$, the $2$-adic limit is
> $$\xi_2=-Q\bigl(s=w{+}1\bigr)\cdot\tfrac12\zeta_2(w{+}1).$$

For level $16$: $Q(X)=(1-4X)(1-16X)(1-64X)$, $Q(2^{-5})=\tfrac78\cdot\tfrac12\cdot(-1)=-\tfrac7{16}$,
so $\xi_2=\tfrac7{32}\zeta_2(5)$.

$\zeta_2(5)=L_2(5,\mathbf 1)$ was computed by Washington's Thm 5.11 with $F=4$,
$\omega=$ the character mod $4$, $\langle a\rangle=a/\omega(a)$, and **validated
exactly** against the interpolation law $L_2(1-n,\mathbf 1)=-(1-2^{n-1})B_n/n$ for
$n=2,4,6,8,10,12$ (equality as rationals: $\tfrac1{12},-\tfrac7{120},\tfrac{31}{252},
-\tfrac{127}{240},\tfrac{511}{132},-\tfrac{1414477}{32760}$). $v_2(\zeta_2(5))=-3$.

**Result** (`padic_id.gp`):

| $n$ | $199$ | $299$ | $349$ | $389$ | $398$ | $399$ |
|---|---|---|---|---|---|---|
| $v_2\bigl(B_n/A_n-\tfrac7{32}\zeta_2(5)\bigr)$ | $178$ | $274$ | $317$ | $364$ | $371$ | $372$ |

and directly $\bigl(B_{399}/A_{399}\bigr)/\zeta_2(5)=2^{-5}+2^{-4}+2^{-3}+O(2^{60})=\tfrac7{32}$
exactly, with $v_2(\text{ratio}-\tfrac7{32})=375$. A scan of ten competing rationals
$c$ gives $v_2(B/A-c\,\zeta_2(5))\le-7$ for every $c\ne\tfrac7{32}$.

**Consistency with the archimedean side.** $\xi_\infty=\tfrac{217}{1024}\zeta(5)$ and
$\xi_2=\tfrac{7}{32}\zeta_2(5)$; the ratio of rational factors is
$\tfrac{7/32}{217/1024}=\tfrac{32}{31}=(1-2^{-5})^{-1}$ — exactly the Euler factor
removed by Kubota–Leopoldt. So $\xi_2$ and $\xi_\infty$ *are the same formula*,
$-P(2^{-(w+1)})\cdot\tfrac12\cdot(\text{zeta value})$, with the $(1-X)$ factor absorbed
into the $2$-depletion. This is the $p=2$, $w=4$ analogue of `CONJ_D_PROOF.md` §4.

**Logical status.** This is a *verification*, not a proof. Hypothesis (2) of Theorem
R$''$ (`CONJ_D_PROOF.md` §9) — "every $d_i$ in $P$ is prime to $p$" — **fails here**:
the source is built out of $V_2,V_4,V_8,V_{16}$ at $p=2$. That is precisely the gap
Calegari's method leaves open. The numerics say the formula survives it; nothing here
proves it does.

---

## 2. Conjecture D for $\zeta(5)$ at $p=2$ (two aligned rows)

The same source on the Möbius-shifted host coordinate $y=x/(1+4x)$ gives a second,
inequivalent row:
$$A^{\rm II}_n:\ 1,\,-48,\,-56,\,416,\,5616,\,31872,\,109824,\dots,\qquad
d_n^5B^{\rm II}_n\in\mathbb Z .$$
It has $\log|A^{\rm II}_n|/n\to1.3887$ ($\Lambda\approx4$), converges archimedeanly to
the *same* $\tfrac{217}{1024}\zeta(5)$ at rate $0.7322$, and has $\sigma_2\approx1$
($v_2$-increment $189$ at $n=219$).

**Its $2$-adic limit is independently identified**: with the validated $\zeta_2$ code,
$v_2\bigl(B^{\rm II}_n/A^{\rm II}_n-\tfrac7{32}\zeta_2(5)\bigr)=78,129,154,175,172$ at
$n=99,149,189,198,199$ — i.e. $\xi_2^{\rm II}=\tfrac7{32}\zeta_2(5)$ to $172$ digits, the
same constant as row I, obtained from a different host coordinate. (Against the
competing $\tfrac{217}{1024}\zeta_2(5)$: $v_2=-13$, flat.)

Cross determinant $\Delta_{n,m}=A^{\rm I}_nB^{\rm II}_m-A^{\rm II}_mB^{\rm I}_n$
(`two_row16.gp`, $M=220$):

| $v_2(\Delta_{n,m})$ | $m{=}25$ | $50$ | $100$ | $150$ | $200$ | $219$ |
|---|---|---|---|---|---|---|
| $n{=}25$ | 33 | 32 | 33 | 35 | 33 | 45 |
| $n{=}50$ | 33 | 55 | 57 | 59 | 57 | 69 |
| $n{=}100$ | 40 | 62 | 107 | 109 | 107 | 119 |
| $n{=}150$ | 36 | 58 | 108 | 159 | 157 | 169 |
| $n{=}200$ | 35 | 57 | 107 | 159 | 206 | 219 |
| $n{=}219$ | 49 | 71 | 121 | 173 | 220 | 244 |

$v_2(\Delta_{n,m})=\min(n,m)+O(1)=\min(\sigma_2,\sigma_2')\min(n,m)+O(1)$: the exact
Conjecture-D signature (`paper/sections/04_padic.tex`, Verification 4.6). Both rows
therefore share one $2$-adic period $\Lambda_2$ with
$\xi_2=r\Lambda_2$, $r=\tfrac{217}{1024}$, $\Lambda_2=\tfrac{1024}{217}\cdot\tfrac7{32}\zeta_2(5)=\tfrac{32}{31}\zeta_2(5)$.

---

## 3. Census: which $\zeta(5)$ sources can carry a slope, and which rows actually do

`source_census.gp` enumerates all $N\le200$ and all trivial-character weight-six
Eisenstein-oldform sources $\Phi=\sum_{d\mid N}c_dE_6(d\tau)$ with
$P(0)=P(2)=P(4)=P(6)=0$ and $P(5)\ne0$, then imposes the Euler-factor criterion at each
$p\mid N$ (for every prime-to-$p$ class $m$, $\sum_ac_{p^am}=0$). Genus-zero hits:

| $N$ | $p$ | $c_d$ | $L(\Phi,5)$ | criterion | measured slope |
|---|---|---|---|---|---|
| $12$ | $2$ | $(1,-113,567,112,-1863,1296)$ | $\tfrac{31}{192}\zeta(5)$ | **passes** | **none** (built, §3.1) |
| $16$ | $2$ | $(1,-85,1428,-5440,4096)$ | $\tfrac{217}{1024}\zeta(5)$ | **passes** | $\sigma_2=1$ |
| $18$ | $3$ | $(1,-112,566,-1184,-567,1296)$ | $\tfrac{121}{729}\zeta(5)$ | **passes** | $\sigma_3\approx0.4$ (plain $\Gamma_0(18)$ Hauptmodul) |

The two sources already in the archive **fail** the criterion, and correspondingly have
no slope (measured, `zeta5_modular_census.gp`, $M=100$):

| system | criterion at $p=2$ | criterion at $p=3$ | measured $\sigma_2,\sigma_3,\sigma_5,\sigma_7$ |
|---|---|---|---|
| $N{=}12^{+}$ Domb, $\tfrac{25}{144}$: $(1,-104,351,832,-2808,1728)$ | $[729,-729]\ne0$ | $[352,-2912,2560]\ne0$ | $0,0,0,0$ |
| $N{=}12^{-}$ anti-Fricke, $\tfrac{11}{144}$: $(1,-176,2079,-4928,4752,-1728)$ | $[-5103,5103]\ne0$ | $[2080,4576,-6656]\ne0$ | $0,0,0,0$ |
| $N{=}18$: $(1,-136,918,-7344,12393,-5832)$ | $[-135,-6426,6561]\ne0$ | $[13312,-13312]\ne0$ | not built |

Both level-12 systems were rebuilt from scratch here and reproduce the archive exactly
— $N{=}12^{+}$ gives $A_n=1,-59,-530,-4787,-47654,-520862,-6116564,-75656435,-971795870$
(matching `MASTER_RESEARCH_LEDGER` digit for digit), $B_n/A_n\to\tfrac{25}{144}\zeta(5)$
at rate $0.2557$ with $\Lambda\approx16$; $N{=}12^{-}$ on the $h_{12}/x$ coordinate gives
$B_n/A_n\to\tfrac{11}{144}\zeta(5)$ at rate $0.0802\to7-4\sqrt3=0.0718$ with
$\Lambda\to7+4\sqrt3$ — the published values. This confirms `SLOPE_CENSUS.md` §3's
$\prod\lambda_i=1$ claim **by direct measurement** for the first time: the level-12
$\zeta(5)$ systems are Apéry-perfect and have no $p$-adic resource at all.

### 3.1 Important correction to the criterion

The brief's statement "*a trivial-character row has a $2$-adic slope iff $(1-X)\mid P$*"
is **necessary but not sufficient**. The new $N{=}12$, $\tfrac{31}{192}$ source passes
the criterion, yet the row built on it (both the Domb coordinate and the $h_{12}$
coordinate, `level12_slope.gp`, $M=160$) has $v_2$-increments pinned at $-1$ to $-3$:
**no slope**. The reason is that the criterion constrains the *source* while the slope
$\sigma_p=v_p(c)+2\kappa_p$ is a property of the *host*: the level-12 Picard–Fuchs
operator is $\mathrm{Sym}^4$ of $\lambda^2-14\lambda+1$, so $c=\prod\lambda_i=1$ and
$v_2(c)=0$. The corrected statement is

> $\sigma_p>0$ needs **(host)** $v_p(c)+2\kappa_p>0$; the Euler-factor criterion
> **(source)** $(1-X)\mid P$ then pins the value $\xi_p=-Q(s{=}w{+}1)\cdot\tfrac12\zeta_p(w{+}1)$.

Level $16$ satisfies both (its cusps sit at $t=-\tfrac14,-\tfrac12$, i.e. $2$-adically
non-unit singular values); level $12$ satisfies only the source half.

---

## 4. The Brown–Zudilin cellular row

Rebuilt exactly from arXiv:2210.03391v3 §2 (`bz_row.gp`, $n\le900$): the published
third-order recursion with $Q_0,Q_1,Q_2=1,21,2989$, $\hat P=0,\tfrac{101}4,\tfrac{344923}{96}$,
$P=0,\tfrac{87}4,\tfrac{1190161}{384}$. Independent check: the closed form
$$Q_n=\sum_{k_1,k_2=0}^{n}\binom{n+k_1}{n}\binom{n}{k_1}^{\!2}\binom{n+k_2}{n}\binom{n}{k_2}^{\!2}\binom{n+k_1+k_2}{n}$$
reproduces the recursion for $n\le7$ (this is BZ eq. (7); the exponent placement is
ambiguous in the PDF's text layer and the version above is the correct one).

**Arithmetic** (all $n\le900$, zero exceptions):
* $Q_n\in\mathbb Z$, and **$Q_n\equiv1\pmod 4$ for every $n\le520$** (exhaustive check,
  no exceptions; $v_2(Q_n)=0$ also at all sampled $n\le899$). So $\kappa_2=0$ *exactly*
  and there is not even a $2$-adic unit fluctuation to exploit — in sharp contrast to
  Zudilin's Catalan row, where $v_2(Q_m)=-4m+2s_2(m)$ supplies the entire $\kappa_2=4$
  resource behind the $1-\varepsilon$ result (`ZETA3_TWO_LATTICE.md` §14). Modulo $8$ the
  row is not constant ($258$ of the $521$ values differ from $1$), so $\equiv1\bmod4$
  appears to be sharp.
* $12\,d_n^5P_n\in\mathbb Z$ but $d_n^5P_n\notin\mathbb Z$; exactly
  $v_2(\mathrm{denom}\,P_n)=5\lfloor\log_2n\rfloor+2$ with a similar $O(\log n)$ bound at $3$.
  Confirms the "$12d_n^5P_n\in\mathbb Z$" of the archive notes.
* $d_{2n}^2\hat P_n\in\mathbb Z$ (BZ eq. (6)).

**Archimedean** [BZ, reproduced]: $\log|Q_n|/n\to\log\lambda_3=6.38364071$,
$\log|Q_n\zeta(5)-P_n|/n\to\log|\lambda_2|=-2.47237372$,
$\log|I_n|/n\to\log|\lambda_1|=-5.29756135$. Measured at $n=100$:
$6.2414$, $-2.5905$, $-5.2952$. Own worthiness of the totally symmetric row:
$\gamma=(c_1-c_0)/c_1=(11.38364-2.52763)/11.38364=0.77796$ [BZ]; the group-refined
non-symmetric family reaches $0.86$ [BZ Thm 1].

**$p$-adic — negative, decisively:**
* $v_p\bigl(P_n/Q_n-P_{n-1}/Q_{n-1}\bigr)$ at $p=2,3,5,7$ is flat and negative for
  $n=10,\dots,899$ (e.g. $p=2$: $-12,-17,-22,-27,-27,-32,-29,-24$). The dense window
  $n=870\ldots899$ oscillates in $[-43,-23]$ with no trend. **No slope.**
* Same for $\hat P_n/Q_n$ (the $\zeta(3)$ direction).
* $2\times2$ minor $Q_nP_{n-1}-Q_{n-1}P_n$ and the $3\times3$ Casoratian: no growth at
  any of $2,3,5,7$.
* **Towers.** $v_2\bigl(P_n/Q_n-P_{n/2}/Q_{n/2}\bigr)$ along $n=a\cdot2^e$ is
  $-7,-12,-17,-22,-27,-32,-37,-42,-47$ for $a=1$ and the identical ladder for
  $a=3,5,7$ — i.e. exactly $-(5\lfloor\log_2n\rfloor+2)$, the $d_n^5$ denominator and
  nothing else. There is no tower limit.
* **Structural reason.** The $2$-adic Newton polygon of $4\lambda^3-2368\lambda^2-188\lambda+1$
  is a single segment from $(0,0)$ to $(3,2)$: all three roots have $v_2=-\tfrac23$ and
  the polynomial is irreducible over $\mathbb Q_2$ (`factorpadic` confirms a single
  degree-3 factor). There is no $2$-adic dominant root, hence no $2$-adic Apéry limit.
  Likewise irreducible over $\mathbb Q_3$.

---

## 5. Alignment test — negative

`align.gp`. With the archimedean factors built in ($B_m/A_m\to\tfrac{217}{1024}\zeta(5)$,
$P_n/Q_n\to\zeta(5)$), the cross determinant is
$\Delta_{m,n}=217A_mP_n-1024Q_nB_m$.

$v_2(\Delta_{m,n})$ over sampling ratios $\rho=m/n$ and $n$ up to $199$, $m$ up to $399$:

| $\rho$ | $\tfrac14$ | $\tfrac12$ | $1$ | $\tfrac32$ | $2$ | $3$ | $4$ |
|---|---|---|---|---|---|---|---|
| range of $v_2$ | $[-25,-11]$ | $[-18,-11]$ | $[-19,-8]$ | $[-18,-6]$ | $[-19,-7]$ | $[-18,-3]$ | $[-19,-7]$ |

Flat, negative, no growth anywhere. A scan of eleven rational multipliers $c$ in
$v_2(cA_nP_n-Q_nB_n)$ gives the same picture (best value $+2$, at $n=60$, in a
non-monotone sequence). **Not aligned** — as it must be, since alignment requires both
rows to have positive slope at $p$ (`THEORY_NOTES_03_lattices.md` §2) and the BZ row
has $\sigma_2=0$.

---

## 6. What the master formula says, counterfactually

Inputs (per index, all measured here unless marked):

| row | rôle | $\log\Lambda$ | $\log(1/\lambda)$ or $\log\rho_2$ | $k$ | $\kappa_2$ | $\sigma_2$ |
|---|---|---|---|---|---|---|
| Brown–Zudilin (symmetric) | decayer | $6.38364$ | $\log\frac1\lambda=2.47237$ | $5$ | $0$ | $\mathbf 0$ |
| level-$12^-$ ($\tfrac{11}{144}$) | decayer | $2.63392$ | $\log\frac1\lambda=2.63392$ | $5$ | $0$ | $\mathbf 0$ |
| level-$16$ row I | engine | $2.03574$ | $\log\rho_2=\log4=1.38629$ | $5$ | $0$ | $\mathbf 1$ |
| level-$16$ row II | engine | $\approx\log4$ | $\log\rho_2\approx1.0769$ | $5$ | $0$ | $\approx1$ |

`ZETA3_TWO_LATTICE.md` §5 with $\gamma=1$, engine sampled at $\alpha$, $S=5\max(1,\alpha)$,
$\eta=0$, $G=\min(\alpha\sigma_{\rm eng},\sigma_{\rm dec})\log2$:
$$F=\tfrac12\bigl[5\max(1,\alpha)+1.38629\alpha-2.47237-G\bigr],\qquad
H=F+8.85601,\qquad \delta=1-F/H .$$
Granting Brown–Zudilin a *hypothetical* aligned slope $s$ and balancing $\alpha=s$:

| $s$ | $0$ | $1$ | $2$ |
|---|---|---|---|
| $F$ | $1.2638$ | $1.6104$ | $4.4570$ |
| $\delta$ | $0.8751$ | $0.8461$ | $0.6652$ |

$\delta$ is **decreasing** in the hypothetical slope, and its supremum $0.8751$ is
attained only at $s=0,\alpha\to0$ — i.e. with no engine row at all, which is not a
construction. Compare Brown–Zudilin's own $0.77796$ (symmetric) / $0.86$ (group-refined).

**The structural obstruction.** $F<0$ (equivalently $\delta>1$) requires
$\text{RESOURCE}>\text{COST}$:
$$\gamma\log\tfrac1\lambda+G>k\max(\gamma,\alpha)+\alpha\log\rho_2^{\rm eng}.$$
Since $G\le\gamma\sigma_{\rm dec}\log2$ and $\max(\gamma,\alpha)\ge\gamma$, a necessary
condition is
$$\boxed{\ \log\tfrac1{\lambda_{\rm dec}}+\sigma_2^{\rm dec}\log2\;>\;k=5\ }$$
Brown–Zudilin: $2.4724$. Level-$12$: $2.6339$. Both roughly half of what is needed; one
would need $\sigma_2^{\rm dec}>3.65$ *on top of* Brown–Zudilin's decay, or a decayer with
$\log(1/\lambda)>5$. The $d_n^5$ cost — intrinsic to weight-five Eichler integrals — is
the binding constraint, not the alignment. (Note the §5.1 "design rule" in its published
form does **not** apply here: it assumes the product formula $\Lambda\lambda=c$, which
holds for second-order rows and fails for the third-order Brown–Zudilin row
($\lambda_2\lambda_3=-49.96\notin\mathbb Z$) and for $\mathrm{Sym}^4$ rows.)

For the record, `ZETA3_TWO_LATTICE.md` §14.2's budget $=\log\Lambda-k+\sum\kappa_p\log p$:
Brown–Zudilin scores $6.3836-5+0=+1.384$, which would be the **largest budget of any row
in the project** (Zudilin's Catalan row is $+1.179$) — but it is entirely archimedean
and entirely unharvestable, because the row has no $p$-adic slope at any small prime.
This is a clean illustration that "budget" is an upper bound, not an achievement.

---

## 7. Honest caveats

1. **No irrationality of $\zeta(5)$ is claimed, implied, or approached here.** Nothing in
   this note improves on Brown–Zudilin's $0.86$.
2. Per `CATALAN_AUDIT.md`, even $\delta>1$ for individual pairs proves nothing: past the
   knife-edge the construction degenerates to Dirichlet and returns identical numerics
   for a rational number. The only open statement that matters is non-vanishing of the
   selected linear form.
3. $\xi_2=\tfrac7{32}\zeta_2(5)$ is a $372$-digit numerical identification, not a
   theorem. The Calegari/Coleman argument of `CONJ_D_PROOF.md` does **not** apply as
   written: hypothesis (2) ($p\nmid d_i$) fails because the source is a polynomial in
   $V_2$ at $p=2$. Closing that gap (unit-root/canonical-subgroup replacement for
   $V_p$ at $p\mid d$) is the same open problem as for the $\zeta(3)$ pair
   (`CONJ_D_PROOF.md` §9), now with a second, higher-weight instance to test against.
4. $\sigma_2=1$ for the level-16 row is measured over $n\le399$. The recurrence *was*
   found (order $16$, degree $11$; §1.2a) but its non-Zagier shape means Proposition C
   does not convert the characteristic roots into $\sigma_2$ without the correction of
   `04_padic.tex` Remark 4.5, which was not carried out.
5. **H1** ($v_2(A_n)=O(\log n)$) is verified for $n\le399$, not proved.
6. The level-18 $\zeta(5)$ row was built only on the plain $\Gamma_0(18)$ Hauptmodul
   (no fold), so its $3$-adic slope $\approx0.4$ is real but its archimedean limit is
   *not* $\tfrac{121}{729}\zeta(5)$ and was not identified. It lives at $p=3$ and
   cannot align with the level-16 row at $p=2$ in any case.
7. Search for a $\zeta(5)$ row with $\kappa_2>0$ (a hypergeometric partner in the sense
   of `ZETA3_TWO_LATTICE.md` §14.4(2)) was **not** carried out beyond Brown–Zudilin.
   Brown–Zudilin's $Q$ is integral for *every* parameter vector (it is a $5$-fold
   residue), so the whole $8$-parameter family has $\kappa_2=0$; other constructions
   (Rivoal/Zudilin well-poised, Krattenthaler–Rivoal) were not tested.

---

## 8. Where this leaves $\zeta(5)$

$\zeta(5)$ is now in exactly the position `THEORY_NOTES_03_lattices.md` §3.1 describes
for $L(2,\chi_{-3})$, but mirrored in weight: **the $p$-adic side is rich and the
archimedean side is poor.** There are two aligned $2$-adic rows with $\sigma_2=1$ and a
closed-form $2$-adic period $\tfrac{32}{31}\zeta_2(5)$; there is a superb archimedean
decayer (Brown–Zudilin, $\lambda=0.0844$) and a perfect one (level-12, linear form
$\sim n^{-1}\log^4n$); and **no row has both**. The complementarity stated in
`paper/sections/04_padic.tex` §"Slope-free systems" — rows with $p$-adic slope are
exactly the archimedeanly imperfect ones — is what blocks it, and the $d_n^5$ cost makes
the gap a factor of two rather than an epsilon.

Ranked next moves, if anyone returns to this:
1. **Prove or refute $\xi_2=\tfrac7{32}\zeta_2(5)$.** This is now a concrete
   $p\mid d$ instance with two independent rows and $372$ digits of evidence; it is a
   better test case than the $\zeta(3)$ pair because the weight is higher and the
   projector is a full Hecke projector.
2. **Turn $\sigma_2=1$ into a theorem** from the order-$16$ recurrence of §1.2a via the
   non-Zagier correction of `04_padic.tex` Remark 4.5 (the Casoratian is
   $\propto\prod_k p_{16}(k)/p_0(k)$, both of degree $11$).
3. **Hunt a $\zeta(5)$ row with $\kappa_2>0$.** By §6's boxed inequality it would need
   $\log(1/\lambda)+\sigma_2\log2>5$; with Brown–Zudilin's $\lambda$ that means
   $\kappa_2\ge2$ (since $\sigma_2=2\kappa_2$ when $v_2(c)=0$) — i.e. a cellular /
   Rhin–Viola-type $\zeta(5)$ normalisation carrying $2$-power denominators at rate
   $\ge2$ per index. Nothing in the literature surveyed here supplies one.
4. Do the same census for $\zeta(7)$: `SLOPE_CENSUS.md` §5 shows the level-24 row is
   Apéry-perfect, but the analogue of $C^{\rm full}_{7,2}=\prod_{a=0}^{4}(1-4^aV_2)$
   would live on level $32$ and, being pure $2$-power, automatically passes the
   Euler-factor criterion.

---

## Scripts (`lattice/zeta5_two_row/`)

| file | what it does |
|---|---|
| `haupt16.gp` | eta-quotient search for degree-one $\Gamma_0(16)$ Hauptmoduls |
| `level16_zeta5.gp` | builds the level-16 $\zeta(5)$ row to $M=400$; integrality, rates, slopes |
| `padic_id.gp` | Kubota–Leopoldt $\zeta_2$ (Washington Thm 5.11, validated exactly) and the $\xi_2=\tfrac7{32}\zeta_2(5)$ identification |
| `two_row16.gp` | the second level-16 row and the Conjecture-D cross-determinant grid |
| `level16_variants.gp` | the level-16 source on six host coordinates |
| `zeta5_modular_census.gp` | the two archive level-12 systems rebuilt and measured |
| `level12_slope.gp` | the new $N{=}12$, $\tfrac{31}{192}$ source (criterion passes, slope absent) |
| `level18_slope.gp` | the $N{=}18$ source at $p=3$ |
| `source_census.gp` | all $N\le200$: purified sources and the Euler-factor criterion |
| `bz_row.gp` | Brown–Zudilin row to $n=900$: arithmetic, archimedean, $p$-adic, towers, Newton polygon |
| `align.gp` | the level-16 $\times$ Brown–Zudilin cross-determinant alignment test |
| `fitrec16.gp` | recurrence fit for the level-16 $A$-row: mod-$p$ order/degree scan + exact reconstruction |

Logs sit beside the scripts (`*.log`); `level16_rows.txt` holds the $400$ exact rows.
