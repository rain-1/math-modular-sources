# Task 2 — the Γ₁(5) new period ξ at **both** real places of ℚ(√5)

**Verdict up front (item 4c), CORRECTED.** *An earlier draft of this file said the second
place fails unconditionally. That was wrong, and the error is recorded in §4(c) below: I
identified "the fold at $v_2$" with the singularity of $\mathcal L$ nearest $0$, which is
place-independent because $\mathcal L$ is defined over $\mathbb Q$. But which singularity must be
removed is fixed by the **pure module**, whose puncture $s=1/\lambda_2=t_2=-\varphi^5$ lies in $K$
and **does** move under $\sigma$. At $v_2$ the allowed puncture is $v_2(s)=+0.0902$ and the point
to be removed is $v_2(t_1)=-11.09$.*

**Corrected verdict.** The *source side is fully Galois-equivariant*: $\Phi_{\rm new}$ is
fold-regular at the $v_1$-fold $t_1$ and its conjugate $\Phi'_{\rm new}=\sigma(\Phi_{\rm new})$ is
fold-regular at the $v_2$-fold $t_2$ (§4(d) already showed this: the constant terms vanish at
cusp $0$ and at cusp $1/2$ respectively). The far-cusp periods exist and are (§4(e), two
independent methods, 174+ digits):
$$\pi_D=-\frac{11\,\zeta(2)}5,\qquad \pi'=\xi'=-\varphi^{-5}\mathrm{Im}L(2,\psi_4)-\mathrm{Re}L(2,\psi_4).$$
**The obstruction is therefore not analytic but arithmetic, and it lives entirely in the
hypothesis.** One $K$-linear relation
$a+b\frac{\zeta(2)}5+c\,\xi=0$ supplies a conditional function at $v_1$ **only**; the two-place
construction additionally needs the *independent* relation
$\sigma(a)-11\,\sigma(b)\frac{\zeta(2)}5+\sigma(c)\,\xi'=0$, and for every hypothesis actually
involving $\xi$ (i.e. $c\ne0$) these two are non-proportional $\mathbb Q$-linear functionals on
$\{1,\zeta(2),\mathrm{Re}L,\mathrm{Im}L\}\cdot\{1,\sqrt5\}$. Note also the **$-11$**: the
far-cusp period datum is *literally* $\sigma$-conjugate on the $\xi$-slot ($\pi'=\sigma(\xi)$) but
is $-11\times$ on the $\zeta(2)$-slot, so the second condition is not even the naive $\sigma$-image
of the first. Full statement in §4(c)/§4(f).

Scripts (all in this directory, all PARI/GP 2.15):
`01_build.gp` (exact construction, writes `build.txt`), `02_limits.gp`, `03_ode.gp`,
`04_ode_verify.gp`, `05_coeffs.gp`, `06_denoms.gp`, `07_foldreg.gp`, `08_outer.gp`,
`09_outerfit.gp`, `10_cusps.gp`, `11_cusps2.gp`, `12_cusptable.gp`, `13_haupt.gp`,
`14_rates.gp`, `15_checks.gp`, `16_summary.gp`.
All exact data (200 coefficients of $A$, $B_D$, $B_3$, $B_4$ as exact rationals) is in `build.txt`.

Throughout: $\varphi^5=\tfrac{11+5\sqrt5}{2}=11.0901699\ldots$, $\varphi^{-5}=\tfrac{5\sqrt5-11}{2}=0.0901699\ldots$,
$t_1=\varphi^{-5}$, $t_2=-\varphi^{5}$, $t_1t_2=-1$, $t_1+t_2=-11$, $t_1-t_2=5\sqrt5$;
$K=\mathbb Q(\sqrt5)$, $\sigma:\sqrt5\mapsto-\sqrt5$, so $\sigma(\varphi^5)=-\varphi^{-5}$;
$v_1$ is the place $\sqrt5\mapsto+2.236\ldots$, $v_2$ the place $\sqrt5\mapsto-2.236\ldots$.

---

## 0. Set-up and calibration  **[verified]**

$x=q\prod_{n\ge1}(1-q^n)^{5\left(\frac n5\right)}=q-5q^2+15q^3-30q^4+40q^5-26q^6-30q^7+125q^8-\cdots$,
$F=\sum_{n\ge0}A_nx^n$, $A_n=\sum_k\binom nk^2\binom{n+k}k=1,3,19,147,1251,11253,\ldots\in\mathbb Z$
(peeling $F$ back in the $x$-coordinate returns the $A_n$ exactly; recurrence
$(n+1)^2A_{n+1}=(11n^2+11n+3)A_n+n^2A_{n-1}$ checked for $n\le200$).

Real Eisenstein directions (rational $q$-coefficients), with
$P_1=E_3^{\mathbf 1,\psi_4}$, $P_2=E_3^{\mathbf 1,\bar\psi_4}$, $P_3=E_3^{\psi_4,\mathbf 1}$, $P_4=E_3^{\bar\psi_4,\mathbf 1}$,
$R_1=P_1+P_2$, $R_2=i(P_1-P_2)$, $R_3=P_3+P_4$, $R_4=i(P_3-P_4)$:

```
R1 = 2q + 2q^2 + 2q^3 - 30q^4 + 2q^5 + 74q^6 + 2q^7 - 30q^8 + ...
R2 =      -8q^2 +18q^3 -  8q^4       + 10q^6 - 98q^7 +120q^8 + ...
R3 = 2q + 8q^2 +18q^3 + 30q^4 + 50q^5 + 74q^6 + 98q^7 +120q^8 + ...
R4 =      -2q^2 + 2q^3 -  8q^4       - 10q^6 -  2q^7 - 30q^8 + ...
```

* $\Phi_D=\tfrac12R_1+R_2=q-7q^2+19q^3-23q^4+q^5+47q^6-97q^7+105q^8-62q^9+\cdots$
  and **$\Phi_D=F\cdot\theta_q x$ exactly** (checked to $O(q^{40})$, `15_checks.gp`) — the calibration
  against `hostscan/REPORT.md` §10.1 reproduces.
* $\Phi_{\rm new}=R_3+\varphi^5R_4=(1+i\varphi^5)P_3+(1-i\varphi^5)P_4$
  $=2q+(-3-5\sqrt5)q^2+(29+5\sqrt5)q^3+(-14-20\sqrt5)q^4+50q^5+\cdots\in\mathbb Z[\varphi]\llbracket q\rrbracket$.
* $\Phi_{\rm new}'=\sigma(\Phi_{\rm new})=R_3-\varphi^{-5}R_4$.

Companions: $B_\bullet$ = the $x$-expansion of $F\cdot D^{-2}\Phi_\bullet$, peeled with `peel2`.
Since $R_3,R_4$ are rational, $B_{\rm new}=B_3+\varphi^5B_4$ and $B_{\rm new}'=B_3-\varphi^{-5}B_4$
with $B_3,B_4\in\mathbb Q$; **all of Tasks 1–3 are therefore done in exact rational arithmetic**
to $n=200$.

---

## 1. Apéry limits and the first coefficients

### 1(a)  **[verified, 418 digits]**
$\displaystyle\lim_n \frac{B_{D,n}}{A_n}=\frac{\zeta(2)}5$: at $n=200$ the ratio agrees with
$\zeta(2)/5=0.328986813369645287294483033329205037843789980241359687547112\ldots$ to
**418 decimal digits** (`02_limits.gp`).

### 1(b)  **[verified, 417 digits]**
$\displaystyle\lim_n\frac{B_{{\rm new},n}}{A_n}=\xi
=0.655634188840656766330981413872399402411138459136702028463135\ldots$,
agreeing at $n=200$ with $\varphi^5\,{\rm Im}\,L(2,\psi_4)-{\rm Re}\,L(2,\psi_4)$ to
**417 decimal digits**. (${\rm Re}\,L(2,\psi_4)=0.958716122716883155\ldots$,
${\rm Im}\,L(2,\psi_4)=0.145565876785089590\ldots$)

`lindep` certificate at 200-digit precision (`07_foldreg.gp`):
$$\texttt{lindep}\big[\xi,\;1,\;\zeta(2),\;{\rm Re}L,\;{\rm Im}L,\;\sqrt5,\;\sqrt5{\rm Re}L,\;\sqrt5{\rm Im}L\big]
=[\,2,\;0,\;0,\;2,\;-11,\;0,\;0,\;-5\,]$$
i.e. $2\xi+2{\rm Re}L-(11+5\sqrt5){\rm Im}L=0$, i.e. $\xi=\varphi^5{\rm Im}L-{\rm Re}L$,
with **no** rational part, **no** $\zeta(2)$ admixture and no $\sqrt5$-only term.
This reproduces `hostscan/REPORT.md` §10.4 exactly.

For the conjugate value: $\texttt{lindep}[\xi',\ldots]=[-2,0,0,-2,11,0,0,-5]$, i.e.
$\xi'=\tfrac{11-5\sqrt5}2{\rm Im}L-{\rm Re}L=-\varphi^{-5}{\rm Im}L-{\rm Re}L
=-0.971841789638437582036104642602905464314\ldots$ — **but see §4(b): $\xi'$ is *not* an Apéry limit.**

### 1(c) First coefficients  **[exact]**

$$B_D:\quad 0,\;1,\;\frac{25}{4},\;\frac{1741}{36},\;\frac{6585}{16},\;\frac{13327519}{3600},\;
\frac{124308457}{3600},\;\frac{19427741063}{58800},\;\frac{2273486234953}{705600},\;\ldots$$

$B_{\rm new,n}=a_n+b_n\sqrt5$:

| $n$ | $a_n$ | $b_n$ |
|---|---|---|
| 0 | 0 | 0 |
| 1 | $2$ | $0$ |
| 2 | $61/4$ | $-5/4$ |
| 3 | $4733/36$ | $-565/36$ |
| 4 | $9625/8$ | $-685/4$ |
| 5 | $822319/72$ | $-65095/36$ |
| 6 | $8002129/72$ | $-75715/4$ |
| 7 | $1294849775/1176$ | $-38796015/196$ |
| 8 | $312067901387/28224$ | $-6501277035/3136$ |

$B_{\rm new}'$ has the same $a_n$ and the opposite $b_n$ (it is $\sigma(B_{\rm new})$).
Underlying rational data: $B_3=0,2,18,166,\frac{12639}8,\frac{123193}8,\frac{11000443}{72},\ldots$ and
$B_4=0,0,-\frac12,-\frac{113}{18},-\frac{137}2,-\frac{13019}{18},-\frac{15143}2,\ldots$;
$a_n=B_{3,n}+\tfrac{11}2B_{4,n}$, $b_n=\tfrac52B_{4,n}$.

---

## 2. Sharp denominators  **[verified, N = 200]**

Let $d_n={\rm lcm}(1,\ldots,n)$ and write $B_{{\rm new},n}=\tfrac{u_n+v_n\sqrt5}{2}$ with
$u_n=2B_{3,n}+11B_{4,n}$, $v_n=5B_{4,n}$; recall
$\mathbb Z[\varphi]=\{\tfrac{u+v\sqrt5}2:u,v\in\mathbb Z,\ u\equiv v\ (2)\}$. Then (`06_denoms.gp`):

* $A_n\in\mathbb Z$ for all $n\le200$.
* $d_n^2B_{D,n}\in\mathbb Z$ for **all** $1\le n\le200$, and $d_n^2u_n,\,d_n^2v_n\in\mathbb Z$
  with $d_n^2u_n\equiv d_n^2v_n\ (\mathrm{mod}\ 2)$ for **all** $1\le n\le200$; hence
  $d_n^2B_{{\rm new},n}\in\mathbb Z[\varphi]$ for all $n\le200$. The same holds separately for
  $B_{3,n}$ and $B_{4,n}$ (max exponent exactly 2 in both cases).
* **Sharpness.** $d_n^{1}B_{D,n}\notin\mathbb Z$ for **199 of the 200** values $n\le200$
  (the single exception is $n=1$), and $d_n^{1}B_{{\rm new},n}\notin\mathbb Z[\varphi]$ for
  **197 of the 200** values (exceptions $n=1,2,4$). The largest gap between consecutive
  "$d_n^1$ fails" indices is $1$ (for $B_D$) and $2$ (for $B_{\rm new}$) throughout $n\le200$.

**Statement supported:** the denominator type of both companions is exactly
$[1,\ldots,n]^2$ over $\mathbb Z$ resp. $\mathbb Z[\varphi]$, verified for $n\le200$ and sharp
at essentially every $n$ in that range. This is a numerical (exhaustive to $N=200$) statement,
not a proof for all $n$; it matches the $k=2$ measurement of `23b_knew.gp` (which went to $n=76$).

---

## 3. The Picard–Fuchs operator and the exact right-hand sides

### 3(a) The operator  **[proved / verified $n\le200$]**
From $(n+1)^2A_{n+1}=(11n^2+11n+3)A_n+n^2A_{n-1}$ one gets
$\theta^2y-x(11\theta^2+11\theta+3)y-x^2(\theta+1)^2y=0$, i.e. **the expected normalisation is
correct**:
$$\boxed{\;\mathcal L y \;=\; x(1-11x-x^2)\,y''\;+\;(1-22x-3x^2)\,y'\;-\;(3+x)\,y\;}$$
with $[x^{n-1}]\mathcal L y = n^2y_n-(11n^2-11n+3)y_{n-1}-(n-1)^2y_{n-2}$.
Singular points $0,\;t_1=\varphi^{-5},\;t_2=-\varphi^5,\;\infty$ (as $1-11x-x^2=-(x-t_1)(x-t_2)$),
all four regular singular; at $t_1$ and $t_2$ the exponents are $0,0$ (unipotent local
monodromy: the leading coefficient $p(x)=x(1-11x-x^2)$ and $q(x)=1-22x-3x^2$ satisfy
$q(t_i)=p'(t_i)$, checked exactly). $\mathcal L A=0$ verified for all $n\le200$.

### 3(b) The right-hand sides  **[verified: all 198 available coefficients, no fit residual]**
$$\boxed{\;\mathcal L B_D = 1,\qquad
\mathcal L B_3=\frac{2}{1-11x-x^2},\qquad
\mathcal L B_4=\frac{-2x}{1-11x-x^2}\;}$$
(the products $(\mathcal L B_\bullet)\cdot(1-11x-x^2)$ are the *exact polynomials*
$1-11x-x^2$, $2$, $-2x$; all 198 computable coefficients of the difference vanish, `04_ode_verify.gp`).

Factoring: since $2-2\varphi^5x=2\varphi^5(x-t_1)$ and $2+2\varphi^{-5}x=2\varphi^{-5}(x-t_2)$,
$$\boxed{\;\mathcal L B_{\rm new}\;=\;\frac{2}{1-x/t_2}\;=\;\frac{2}{1+\varphi^{-5}x},\qquad
\mathcal L B_{\rm new}'\;=\;\frac{2}{1-x/t_1}\;=\;\frac{2}{1-\varphi^{5}x}.\;}$$
**The pole of the inhomogeneity of the fold-regular source sits at $t_2$, and $\sigma$ moves it
to $t_1$.** This single line is the whole of item 4(c).

### 3(c) The conditional generator  **[proved, given 3(a),(b)]**
If $a+b\,\frac{\zeta(2)}5+c\,\xi=0$ with $a,b,c\in K$ and
$H=aA+bB_D+cB_{\rm new}=b\big(B_D-\tfrac{\zeta(2)}5A\big)+c\big(B_{\rm new}-\xi A\big)$, then
$$\boxed{\;x(1-11x-x^2)H''+(1-22x-3x^2)H'-(3+x)H\;=\;b+\frac{2c}{1+\varphi^{-5}x}\;=\;b+\frac{2c}{1-x/t_2}\;}$$
which is the exact structural analogue of CDT's Proposition `functionsH`
($x(1-x)(1-9x)y''+(1-20x+27x^2)y'+3(-1+3x)y=b+\frac c{1-x}$): a constant plus a simple pole at
the *non-fold* finite singular point ($x=1$ there, $x=t_2$ here). $H\in K\llbracket x\rrbracket$
with denominators of shape $d_n^2$ in $\mathbb Z[\varphi]$ (§2).

The $\sigma$-conjugate satisfies
$$\mathcal L H^\sigma \;=\;\sigma(b)+\frac{2\sigma(c)}{1-\varphi^5x}\;=\;\sigma(b)+\frac{2\sigma(c)}{1-x/t_1}.$$

### 3(d) Aside: the two outer directions individually  **[computed]**
$R_1$ and $R_2$ separately do **not** give a rational right-hand side (their
$\mathcal L B$ has radius of convergence $\varphi^{-5}$, ratios $\to11.09$, `09_outerfit.gp`).
The reason is §4(d): $E_3^{\mathbf 1,\psi_4}$ has a **non-zero constant term at $\infty$**
($-\tfrac{2+i}5$), which the Eichler integral $D^{-2}$ silently discards, destroying modularity.
The combination $\Phi_D=\tfrac12R_1+R_2$ is *precisely* the one whose $\infty$-constant term
cancels ($\mathrm{Re}\,c_0-2\,\mathrm{Im}\,c_0=-\tfrac25+\tfrac25=0$), which is why it, and only
it, has a clean rational RHS. $R_3,R_4$ have $c_0(\infty)=0$ automatically.

---

## 4. Fold-regularity and the second real place

### 4(a) Overconvergence at $t_1$ at the first place  **[verified, 200 coefficients]**
With $c_n=B_{D,n}-\frac{\zeta(2)}5A_n$ and $c_n=B_{{\rm new},n}-\xi A_n$ (`07_foldreg.gp`, `14_rates.gp`):

| quantity | $c_n/c_{n-1}$ at $n=200$ | $n\cdot(c_n/c_{n-1}+\varphi^{-5})$ |
|---|---|---|
| $B_D-\frac{\zeta(2)}5A$ | $-0.0897207202335\ldots$ | $0.08984$ (bounded) |
| $B_{\rm new}-\xi A$ | $-0.0897747710569\ldots$ | $0.07903$ (bounded) |

so $c_n/c_{n-1}=-\varphi^{-5}\big(1+O(1/n)\big)$ with $\varphi^{-5}=0.0901699437494742410\ldots$.
Both series are **overconvergent at $t_1=\varphi^{-5}$**: the decay rate is
$|t_2|^{-n}=\varphi^{-5n}$, not $|t_1|^{-n}$, and the **negative** ratio locates the
surviving singularity at $x=t_2=-11.0901699\ldots$, on the negative axis. ($|c_{200}|^{1/200}$
= $0.087208$ resp. $0.087984$, approaching $0.0901699$ from below at the expected $O(\log n/n)$ rate.)

Consequently $B_{{\rm new},n}/A_n-\xi=O(\varphi^{-10n})$, worth $10\log_{10}\varphi=2.0898$ digits
per step — matching the observed 417 digits at $n=200$.

### 4(b) The conjugate source, and a correction to the expected statement

$\Phi'_{\rm new}:=\sigma(\Phi_{\rm new})=R_3-\varphi^{-5}R_4$ (using $\sigma(\varphi^5)=-\varphi^{-5}$).
**Confirmed** (identity, and checked coefficientwise to $O(q^{40})$, `15_checks.gp`):
$$\Phi'_{\rm new}=(1-i\varphi^{-5})\,E_3^{\psi_4,\mathbf 1}+(1+i\varphi^{-5})\,E_3^{\bar\psi_4,\mathbf 1}.$$

**CORRECTION to the expected result.** $B'_{\rm new}$ has **no Apéry limit at all**:
$$\frac{B'_{{\rm new},n}}{A_n}\;\longrightarrow\;+\infty,\qquad
\frac{B'_{{\rm new},n}}{A_n}\;=\;\kappa\log n+O(1),\quad \kappa=0.504694\ldots$$
Measured $n\,(r_n-r_{n-1})$ = $0.5053716\ (n{=}200)$, $0.5055975\ (150)$, $0.5060496\ (100)$,
$0.5074073\ (50)$; Richardson $(a+b/n)$ extrapolation from any two of these gives
$\kappa=0.5046933\ldots$ consistently to 6 digits (`16_summary.gp`). The values themselves are
$r_{200}=3.212501$, $r_{150}=3.067505$, $r_{100}=2.863260$, $r_{50}=2.514608$.
So the number $\xi'=-\varphi^{-5}{\rm Im}L-{\rm Re}L$ is the correct **Galois conjugate of $\xi$**
(it is $\sigma$ applied to the $K$-coefficients of $\xi$'s closed form, certified by `lindep`
above), but it is **not** the limit of $B'_{{\rm new},n}/A_n$, and no companion on this host
converges to it.

*Why, exactly.* By §3(b), $\mathcal L B'_{\rm new}=2/(1-x/t_1)$ has a simple pole **at the
singular point $t_1$ itself**. Near $t_1$ the operator reads $c_1s\,y''+c_1y'+\cdots$ with
$s=x-t_1$, $c_1=-t_1(t_1-t_2)$; solving $(sy')'=r/(c_1s)$ with $r=-2t_1$ gives a particular
solution
$$y_p=\frac{r}{2c_1}\log^2 s=\frac{1}{t_1-t_2}\log^2(x-t_1)=\frac{1}{5\sqrt5}\log^2(x-t_1),$$
an **irremovable $\log^2$**: no homogeneous solution (span $\{1,\log\}$ locally) and no multiple
of $B_D$ (whose RHS $\equiv1$ is analytic everywhere) can cancel it. Hence for **every**
$\lambda,\mu\in\mathbb C$,
$$\big(B'_{{\rm new},n}-\lambda A_n-\mu B_{D,n}\big)\,t_1^{\,n}\,n \;=\;\frac{2}{5\sqrt5}\big(\log n+O(1)\big).$$
Verified numerically: with $\lambda=\xi'$, $g_nt_1^nn/\log n$ = $0.27953\ (n{=}200)$,
$0.28521\ (150)$, $0.29436\ (100)$, all three fitting
$\frac{2}{5\sqrt5}\big(1+2.98/\log n\big)$ with $\frac2{5\sqrt5}=0.1788854$ and the same constant
$2.98$ — a 3-point confirmation of the exact predicted coefficient $1/(5\sqrt5)$.
(Auxiliary: $A_n\sim\varphi^{5n+5/2}/(2\pi5^{1/4}n)$; Richardson gives
$\lim n\varphi^{-5n}A_n=0.3544426$ against $\varphi^{5/2}/(2\pi5^{1/4})=0.3544433$
**[verified, 6 digits]**, and $\kappa=\frac{2}{5\sqrt5}/0.3544433=0.504687$, matching the measured
$0.504693$ to 5 digits.)

### 4(c) **THE KEY TEST — verdict: FAILS, and is unrepairable**

**What actually changes between the two places.** The operator $\mathcal L$ has *rational*
coefficients, so its singular set $\{0,\;0.0901699\ldots,\;-11.0901699\ldots,\;\infty\}\subset\mathbb C$
is Galois-stable and is **the same at both places**; in particular the *fold* (the finite
non-zero singularity nearest $0$) is the complex number $0.0901699\ldots$ at both places.
What moves is the labelling by $K$-numbers, and with it the pole of the inhomogeneity:

| $K$-number | realisation at $v_1$ | realisation at $v_2$ |
|---|---|---|
| $t_1=\varphi^{-5}$ | $+0.0901699$ — **inner (the fold)** | $-11.0901699$ — outer |
| $t_2=-\varphi^{5}$ | $-11.0901699$ — outer | $+0.0901699$ — **inner (the fold)** |

$\mathcal L B_{\rm new}=2/(1-x/t_2)$: the pole is at the $K$-number $t_2$, which is **outer at
$v_1$** (harmless — this is fold-regularity) and **inner at $v_2$** (fatal).

**Numerical verification at the second embedding** (`07_foldreg.gp`, `14_rates.gp`).
Evaluating the $K$-coefficients of $B_{\rm new}$ under $v_2$ is the same as evaluating
$B'_{\rm new}=\sigma(B_{\rm new})$ under $v_1$. For $g_n=B'_{{\rm new},n}-\xi'A_n$:

$$g_n/g_{n-1}\big|_{n=200}=+11.0414635\ldots,\qquad n\,(g_n/g_{n-1}-\varphi^5)=-9.741\ \ (\text{bounded}),$$

so $g_n/g_{n-1}=+\varphi^{5}\big(1+O(1/n)\big)$ with $\varphi^5=11.0901699\ldots$. The rate is
$\varphi^{+5n}$, **not** $\varphi^{-5n}$, and the **positive** ratio locates the singularity at
$x=+\varphi^{-5}=t_1^{v_1}$ — i.e. exactly at the fold. The same holds for $\lambda=0$ and for
every other $\lambda$ tried; and by §4(b) it holds for every $(\lambda,\mu)$ as a theorem.
The radius of convergence of the conjugated conditional form at the second place is therefore
$\varphi^{-5}=0.0902$: it does not reach the fold, let alone overconverge past it to $|x|=\varphi^5$.

**The precise answer to the question posed.** Over $K$ the whole inner supply is
$c_3B_3+c_4B_4$, $c_3,c_4\in K$, with
$\mathcal L(c_3B_3+c_4B_4)=\dfrac{2(c_3-c_4x)}{1-11x-x^2}=\dfrac{-2(c_3-c_4x)}{(x-t_1)(x-t_2)}$.

* fold-regular at $v_1$ $\iff$ the numerator kills the $t_1$ pole $\iff c_3=c_4\,t_1=c_4\varphi^{-5}$
  (the line $\langle B_{\rm new}\rangle$, period $\xi$);
* fold-regular at $v_2$ $\iff$ $\sigma(c_3)=\sigma(c_4)t_1$ $\iff c_3=c_4\,\sigma(t_1)=c_4\,t_2=-c_4\varphi^{5}$
  (the line $\langle B'_{\rm new}\rangle$).

These two $K$-lines are distinct (their "slopes" $t_1$ and $t_2$ differ by $t_1-t_2=5\sqrt5\ne0$),
so
$$\boxed{\;c_3=c_4t_1\ \text{and}\ c_3=c_4t_2\ \Longrightarrow\ c_4(t_1-t_2)=0\ \Longrightarrow\ c_3=c_4=0.\;}$$

**Therefore: the second place's fold-regularity is *not* Galois-equivariant, and it is not merely
"an independent relation away" — it is impossible.** Two separate failures, in increasing order
of severity:

1. *(arithmetic)* From one hypothesised relation $a+b\frac{\zeta(2)}5+c\xi=0$ over $K$ one does
   **not** get $\sigma(a)+\sigma(b)\frac{\zeta(2)}5+\sigma(c)\xi'=0$: $\zeta(2)$, ${\rm Re}L$,
   ${\rm Im}L$ are not in $K$, so $\sigma$ does not act on the relation. A $K$-relation among
   $1,\zeta(2),{\rm Re}L,{\rm Im}L$ is *one* $\mathbb Q$-relation among the eight numbers
   $\{1,\zeta(2),{\rm Re}L,{\rm Im}L\}\cdot\{1,\sqrt5\}$; the conjugate relation is a genuinely
   second, independent one. So even the *hypothesis* would have to be doubled.
2. *(analysis — the fatal one)* Grant the doubled hypothesis anyway. The conjugate function
   $H^\sigma=\sigma(a)A+\sigma(b)B_D+\sigma(c)B'_{\rm new}$ satisfies
   $\mathcal L H^\sigma=\sigma(b)+2\sigma(c)/(1-x/t_1)$, whose pole is **at the fold**. Whenever
   $\sigma(c)\neq0$ (equivalently $c\neq0$, i.e. whenever the relation involves $\xi$ at all)
   $H^\sigma$ carries the $\log^2(x-t_1)$ of §4(b) and has radius of convergence $\varphi^{-5}$.
   The CDT architecture needs a conditional function holomorphic past the fold **at every
   archimedean place of $K$**; at $v_2$ no such function exists. This is not a bounded
   "second-place tax" that the $+0.0053$ margin might absorb — the object the argument requires
   does not exist at $v_2$.

**Consequence for the target.** The only $K$-rational sources fold-regular at *both* places are
those whose companion has no inhomogeneity pole at $t_1$ or $t_2$, i.e. (by §4(d)) whose constant
term vanishes at $\infty$, at cusp $0$ **and** at cusp $1/2$ — a **one-dimensional** space,
$\langle\Phi_D\rangle$, carrying the single period $\zeta(2)/5$. The mixed pair
$\{\zeta(2)/5,\;\xi\}$, which is the entire point of the §10.4 hit, is **2-dimensional at each
place separately but 1-dimensional at both places simultaneously**, and the surviving direction
is the old one. If the target was "run CDT at $\Gamma_1(5)$ with the new quartic period",
this kills it, and it kills it for a reason internal to the source space, before any entropy
budgeting.

*(Escapes considered and closed: the host $F$ is $\mathbb Q$-rational and the geometry
$\mathbf P^1\setminus\{0,t_1,t_2,\infty\}$ is defined over $\mathbb Q$, so there is no
"conjugate host" to move to at $v_2$; the module in the holonomy bound is $K$-rational, so its
$v_2$-realisation is forced to be $\sigma(H)$; and a hypothetical second companion converging
to $\xi'$ with $d_n^2$ denominators would have to be fold-regular at $v_1$ with period $\xi'$,
which §4(d) excludes — the fold-regular space at $v_1$ is $\langle\Phi_D,\Phi_{\rm new}\rangle$,
periods $\zeta(2)/5$ and $\xi$ only.)*

### 4(d) Cusp data — exact constant terms, and the dictionary  **[verified]**

$\Gamma_1(5)$ has four cusps $\infty=1/0$, $2/5$, $0=0/1$, $1/2$ (widths $1,1,5,5$), realised by
$\begin{psmallmatrix}1&0\\0&1\end{psmallmatrix},
\begin{psmallmatrix}2&1\\5&3\end{psmallmatrix},
\begin{psmallmatrix}0&-1\\1&0\end{psmallmatrix},
\begin{psmallmatrix}1&0\\2&1\end{psmallmatrix}$.
Constant terms of $F|_3\gamma$ computed with `mfinit`/`mfeisenstein`/`mfslashexpansion`
(`10_cusps.gp`–`12_cusptable.gp`).
**Note a PARI convention trap:** `mfeisenstein(k,C1,C2)` has
$c_n=\sum_{d|n}C_1(d)C_2(n/d)d^{k-1}$, i.e. the **opposite** order to `hostscan`'s
`eis(ch1,ch2)`; so $E_3^{\mathbf 1,\psi_4}=$`mfeisenstein(3,ps4,1)` and
$E_3^{\psi_4,\mathbf 1}=$`mfeisenstein(3,1,ps4)`. ($\psi_4$ = Conrey $\chi_5(2,\cdot)$, $\psi_4(2)=i$,
$\psi_4(-1)=-1$.)

Write $\rho:=\dfrac1{125\cdot5^{3/4}\varphi^{5/2}}=0.000718444762527336595907\ldots$
(so that $|c_0|=1/125$ at the width-5 cusps; PARI's expansions there are in $q^{1/5}$).

| cusp (width) | $E_3^{\mathbf 1,\psi_4}$ | $E_3^{\mathbf 1,\bar\psi_4}$ | $E_3^{\psi_4,\mathbf 1}$ | $E_3^{\bar\psi_4,\mathbf 1}$ | $\Phi_D$ | $\Phi_{\rm new}$ | $\Phi'_{\rm new}$ |
|---|---|---|---|---|---|---|---|
| $\infty$ (1) | $-\frac{2+i}5$ | $-\frac{2-i}5$ | $0$ | $0$ | $\mathbf 0$ | $\mathbf 0$ | $\mathbf 0$ |
| $2/5$ (1) | $\frac{-1+2i}5$ | $\frac{-1-2i}5$ | $0$ | $0$ | $\mathbf{-1}$ | $\mathbf 0$ | $\mathbf 0$ |
| $0$ (5) | $0$ | $0$ | $\rho(-1+i\varphi^5)$ | $\rho(1+i\varphi^5)$ | $\mathbf 0$ | $\mathbf 0$ | $\;2i\cdot5\sqrt5\,\rho=\frac{2i\,5^{3/4}}{125\,\varphi^{5/2}}$ |
| $1/2$ (5) | $0$ | $0$ | $i\rho(-1+i\varphi^5)$ | $-i\rho(1+i\varphi^5)$ | $\mathbf 0$ | $-2i(1{+}\varphi^{10})\rho=-\frac{2i\,5^{3/4}\varphi^{5/2}}{125}$ | $\mathbf 0$ |

(numerically $2\cdot5\sqrt5\rho=0.0160649132708981823896$,
$2(1+\varphi^{10})\rho=0.1781626183058570776311$; the vanishing entries are $0$ to $\ge96$ digits.
$c_0(E_3^{\mathbf 1,\psi_4},\infty)=-B_{3,\psi_4}/6=-(2+i)/5$ checks against the Bernoulli formula.
$1+\varphi^{10}=5\sqrt5\,\varphi^5$; the ratio ${\rm Im}/{\rm Re}$ of $c_0(E_3^{\psi_4,\mathbf 1},0)$
is $-\varphi^5$ to 78 digits; and $c_0(E_3^{\psi_4,\mathbf 1},1/2)=+i\,c_0(E_3^{\psi_4,\mathbf 1},0)$,
$c_0(E_3^{\bar\psi_4,\mathbf 1},1/2)=-i\,c_0(E_3^{\bar\psi_4,\mathbf 1},0)$, both exactly.)

The four constant-term functionals are independent on the 4-dimensional
$M_3^{\rm Eis}(\Gamma_1(5))$ (the matrix above is block-diagonal with both $2\times2$ blocks
non-singular, determinants $10i/25$ and $-2i\rho^2(-1-\varphi^{10})$), so the constant-term map
is an isomorphism onto $\mathbb C^4$.

**The cusp $\leftrightarrow$ $x$ dictionary, independently confirmed** by evaluating the
Hauptmodul directly (`13_haupt.gp`, $\tau\to$ cusp along a vertical line, $\ge38$ digits):
$$x(\infty)=0,\qquad x(0)=t_1=\varphi^{-5}=0.09016994374947424102293417182819058860,$$
$$x(1/2)=t_2=-\varphi^{5}=-11.0901699437\ldots,\qquad x(2/5)=\infty.$$

**Vanishing pattern $\leftrightarrow$ fold-regularity — a perfect bijection:**

| source | cusps with $c_0\ne0$ | image under $x$ | $\mathcal L B$ | fold-regular at |
|---|---|---|---|---|
| $\Phi_D$ | $2/5$ only | $x=\infty$ | $1$ (no finite pole) | **both places** |
| $\Phi_{\rm new}$ | $1/2$ only | $x=t_2$ | $2/(1-x/t_2)$ | $v_1$ only |
| $\Phi'_{\rm new}$ | $0$ only | $x=t_1$ | $2/(1-x/t_1)$ | $v_2$ only |

So: *a $K$-rational source is fold-regular at the place $v$ exactly when its constant term
vanishes at the cusp lying over the $v$-inner singularity* — cusp $0$ (over $t_1$) for $v_1$,
cusp $1/2$ (over $t_2$) for $v_2$. Being fold-regular at both places therefore forces
$c_0=0$ at three of the four cusps ($\infty$, $0$, $1/2$), which by the isomorphism above cuts
the 4-dimensional space down to the **line $\langle\Phi_D\rangle$**. This is the same conclusion
as §4(c), obtained by exact constant-term linear algebra rather than by the ODE — the two routes
agree, and they also agree with `eis/EIS_REPORT.md`'s "$\dim(\text{annihilated at }0,\infty)=2$"
(that census imposes vanishing at two cusps, hence gets 2; imposing the *second place* adds the
third cusp and drops it to 1).

---

## 5. Ledger

| claim | status |
|---|---|
| $F$, $x$, $\Phi_D=F\theta_qx$, $\Phi_{\rm new}$ rebuilt and calibrated against `hostscan` §10 | **[verified]** |
| $\lim B_{D,n}/A_n=\zeta(2)/5$ | **[verified, 418 digits]** |
| $\lim B_{{\rm new},n}/A_n=\xi=\varphi^5{\rm Im}L(2,\psi_4)-{\rm Re}L(2,\psi_4)$ | **[verified, 417 digits]** + `lindep` certificate $[2,0,0,2,-11,0,0,-5]$ |
| $d_n^2B_{D,n}\in\mathbb Z$, $d_n^2B_{{\rm new},n}\in\mathbb Z[\varphi]$, $n\le200$; sharp at 199/200 resp. 197/200 indices | **[verified, N=200]** |
| $\mathcal L=x(1{-}11x{-}x^2)\partial^2+(1{-}22x{-}3x^2)\partial-(3+x)$, $\mathcal LA=0$ | **[proved from the recurrence; verified $n\le200$]** |
| $\mathcal LB_D=1$, $\mathcal LB_3=\frac2{1-11x-x^2}$, $\mathcal LB_4=\frac{-2x}{1-11x-x^2}$ | **[verified, all 198 coefficients, exact]** |
| $\mathcal LB_{\rm new}=2/(1-x/t_2)$, $\mathcal LB'_{\rm new}=2/(1-x/t_1)$ | **[proved from the previous line]** |
| $B_D-\frac{\zeta(2)}5A$ and $B_{\rm new}-\xi A$ overconvergent at $t_1$, rate $-\varphi^{-5}$ | **[verified, 200 coefficients]** |
| $\Phi'_{\rm new}=(1-i\varphi^{-5})E_3^{\psi_4,\mathbf 1}+(1+i\varphi^{-5})E_3^{\bar\psi_4,\mathbf 1}$ | **[proved; identity]** — confirmed, no correction needed |
| $\xi'=-\varphi^{-5}{\rm Im}L-{\rm Re}L$ is $\sigma(\xi)$ | **[verified]** via `lindep` $[-2,0,0,-2,11,0,0,-5]$ |
| $\xi'$ is **not** an Apéry limit; $B'_{{\rm new},n}/A_n\sim0.50469\log n$ | **[verified numerically; proved from $\mathcal LB'_{\rm new}$]** — **corrects the expected statement** |
| $\log^2$ coefficient $=1/(t_1-t_2)=1/(5\sqrt5)$ | **[proved from the local analysis; verified to 3 significant figures at 3 values of $n$]** |
| second-place rate is $+\varphi^{5n}$ (singularity at the fold), not $\varphi^{-5n}$ | **[verified, 200 coefficients]** |
| no $K$-rational inner source is fold-regular at both places | **[proved]** (two ways: ODE pole, and cusp constant terms) |
| exact constant terms at all four cusps | **[verified, $\ge78$ digits; exact closed forms given]** |
| $x(\infty,0,1/2,2/5)=(0,\varphi^{-5},-\varphi^5,\infty)$ | **[verified, 38 digits]** |
| $A_n\sim\varphi^{5n+5/2}/(2\pi5^{1/4}n)$ | **[verified, 6 digits]** (auxiliary only) |

**Not settled / not attempted:** the $d_n^2$ denominator claim is exhaustive to $N=200$ only, not
proved for all $n$; the entropy/holonomy accounting itself was not recomputed (it is moot, since
the second-place object does not exist); nothing here bears on items 2, 3, 5 of `REPORT.md` §10.5.
