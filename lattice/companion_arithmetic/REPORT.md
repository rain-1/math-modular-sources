# Exact-arithmetic census of the companion sequences $b_n$ of the fifteen sporadic modular Apéry rows

*Working directory:* `/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/9a849c0a-95f8-4d19-b342-98033d0d9c03/scratchpad/companions`
*Source repository (read-only):* `/home/ubuntu/code/math-modular-sources`
*Tools:* PARI/GP 2.15.4, exact rational arithmetic throughout. Every table entry below is an exact
computation over the stated finite range; nothing is floating point.

**Evidence tags.** `[verified exact, range]` = exact rational / integer computation over that range,
no rounding anywhere. `[numerical]` = floating point. `[fit]` = a shape found by search and then
checked, not derived.

---

## 0. Set-up, and the data as actually used

### 0.1 Recurrences

$$\text{(R2)}\qquad (n+1)^2A_{n+1}=(an^2+an+b)A_n-c\,n^2A_{n-1},$$
$$\text{(R3)}\qquad (n+1)^3A_{n+1}=(2n+1)(an^2+an+b)A_n-n(cn^2+d)A_{n-1},$$
with the *row* $a_n$ given by $A_0=1,\ A_1=b$ and the *companion* $b_n$ by $B_0=0,\ B_1=1$.
All fifteen rows were regenerated from scratch and cross-checked against
`paper/sections/02_sources.tex`, `consolidation/ROW_LEDGER.md`,
`lattice/thmB_exact/common.gp`, `lattice/euler_criterion/rows.gp`,
`lattice/root_rows/03_sporadic.gp` and `lattice/padic_irrationality/cooper_scores.gp`.

| row | type | $(a,b,c,d)$ | $a_0,\dots,a_5$ | bad primes $p\mid c$ |
|---|---|---|---|---|
| $\mathbf A$ | R2 | $(7,2,-8)$ | $1,2,10,56,346,2252$ | $2$ |
| $\mathbf B$ | R2 | $(9,3,27)$ | $1,3,9,21,9,-297$ | $3$ |
| $\mathbf C$ | R2 | $(10,3,9)$ | $1,3,15,93,639,4653$ | $3$ |
| $\mathbf D$ | R2 | $(11,3,-1)$ | $1,3,19,147,1251,11253$ | — |
| $\mathbf E$ | R2 | $(12,4,32)$ | $1,4,20,112,676,4304$ | $2$ |
| $\mathbf F$ | R2 | $(17,6,72)$ | $1,6,42,312,2394,18756$ | $2,3$ |
| $\alpha$ (Domb) | R3 | $(10,4,64,0)$ | $1,4,28,256,2716,31504$ | $2$ |
| $\gamma$ (Apéry) | R3 | $(17,5,1,0)$ | $1,5,73,1445,33001,819005$ | — |
| $\delta$ | R3 | $(7,3,81,0)$ | $1,3,9,3,-279,-2997$ | $3$ |
| $\varepsilon$ ($=\mathbf T$) | R3 | $(12,4,16,0)$ | $1,4,40,544,8536,145504$ | $2$ |
| $\zeta$ | R3 | $(9,3,\mathbf{-27},0)$ | $1,3,27,309,4059,57753$ | $3$ |
| $\eta$ | R3 | $(11,5,125,0)$ | $1,5,35,275,2275,19255$ | $5$ |
| Cooper $s_7$ | R3 | $(13,4,-27,3)$ | $1,4,48,760,13840,273504$ | $3$ |
| Cooper $s_{10}$ | R3 | $(6,2,-64,4)$ | $1,2,18,164,1810,21252$ | $2$ |
| Cooper $s_{18}$ | R3 | $(14,6,192,-12)$ | $1,6,54,564,6390,76356$ | $2,3$ |

**Two corrections to the task statement.**

1. $\zeta$ has $c=\mathbf{-27}$, not $+27$. This is the value in `common.gp`
   (`["zeta", 3, 9, 3, -27]`), in `ROW_LEDGER.md` (row "AZ $(9,3,-27)$") and in
   `lattice/root_rows/03_sporadic.gp`. With $c=+27$ the (R3) sequence is **not integral** —
   $a_2=81/4$ — so it is not an Apéry-like row at all [verified exact, `chk_zeta.gp`].
   ($\varepsilon=(12,4,16)$ and $\eta=(11,5,125)$ are confirmed as guessed; $\varepsilon$ is the row
   called $\mathbf T$ in `ROW_LEDGER.md`.)
2. Cooper's rows in $(a,b,c,d)$ normalisation are
   $s_7=(13,4,-27,3)$, $s_{10}=(6,2,-64,4)$, $s_{18}=(14,6,192,-12)$, i.e.
   $$(n{+}1)^3A_{n+1}=(2n{+}1)(13n^2{+}13n{+}4)A_n+3n(9n^2{-}1)A_{n-1},\quad A_1=4,$$
   $$(n{+}1)^3A_{n+1}=2(2n{+}1)(3n^2{+}3n{+}1)A_n+4n(16n^2{-}1)A_{n-1},\quad A_1=2,$$
   $$(n{+}1)^3A_{n+1}=2(2n{+}1)(7n^2{+}7n{+}3)A_n-12n(16n^2{-}1)A_{n-1},\quad A_1=6.$$
   (Note $s_{10}=\sum_k\binom nk^4$: $1,2,18,164,1810$ — checked directly.)

All fifteen rows are integral, and none vanishes: $a_n\in\mathbf Z$, $a_n\neq0$
[verified exact, $0\le n\le3000$, all fifteen rows; the four rows with an exact digit law were
regenerated to $n\le10^4$].

### 0.2 The companion source

$\Phi=\sum_{m\ge1}c(m)q^m$ is Table `tab:sources`; $\Theta=D^{-r}\Phi$, $B=F\Theta$, so

$$b_n=\sum_{m=1}^{n}\frac{c(m)}{m^{r}}\,e_{n,m},\qquad e_{n,m}:=[t^n]\bigl(F(q)\,q^m\bigr)\in\mathbf Z .$$

The $c(m)$ were re-typed from `02_sources.tex` into `src.gp` and agree with `common.gp`.
For the parametrisation I use the identity $F=y_0\circ t$ (so $F$ *as a series in $t$* is literally
$\sum_n a_nt^n$), and $q(t)=t\exp\int(K^{-1}-1)\,dt/t$ with $K=P_2y_0^2$ ($r=2$), $K=y_0\sqrt{P_3}$ ($r=3$).

---

## 1. Task 1 — denominators

### 1.1 The exponent

**Result.** The least $k$ with $d_n^{\,k}b_n\in\mathbf Z$ for all $n\le3000$ is

$$k=2\ \ (\mathbf A,\dots,\mathbf F),\qquad k=3\ \ (\alpha,\gamma,\delta,\varepsilon,\zeta,\eta),
\qquad k=2\ \ (s_7,s_{10},s_{18}).$$

[verified exact, $1\le n\le3000$; script `t1a_den.gp`.] The Cooper value $k=2$ — one below the
order of the recurrence — is the *free integration*, and it is the only place in the census where
$k<r$. So $r$ is the correct exponent for the twelve Eisenstein rows and $2$ for Cooper's three,
exactly as stated in `ROW_LEDGER.md` Table 1 and in `cooper_scores.gp`.

### 1.2 $R_n$ versus $d_n^{\,r}$

With $R_n:=\operatorname{lcm}_{m\le n}\bigl(m^{r}/\gcd(m^{r},c(m))\bigr)$ one always has
$R_n\mid d_n^{\,r}$ (because $\operatorname{lcm}_{m\le n}m^r=d_n^{\,r}$).

* $R_n=d_n^{\,r}$ **identically** for the two Apéry-perfect rows $\mathbf D$ and $\gamma$ ($c=\pm1$):
  the source coefficients never cancel a needed prime power. [verified exact, $n\le3000$]
* For the other ten rows $R_n$ is a *proper* divisor of $d_n^{\,r}$; the quotient $d_n^{\,r}/R_n$
  reaches $7$–$14$ decimal digits by $n=3000$ (largest: $\mathbf F$, $14$ digits).

### 1.3 Is $\operatorname{den}(b_n)=R_n$?  **No** — and the failure is entirely at good primes

$\operatorname{den}(b_n)\mid R_n\mid d_n^{\,r}$ for every row and every $n\le3000$
[verified exact], but both divisibilities are strict for almost all $n$:

| row | $\#\{n\le3000:\operatorname{den}b_n=R_n\}$ | $\#\{n:\operatorname{den}b_n=d_n^{r}\}$ | max digits of $R_n/\operatorname{den}b_n$ | max single-prime defect exponent |
|---|---|---|---|---|
| $\mathbf A$ | 86 | 1 | 20 (at $n{=}2502$) | 8 |
| $\mathbf B$ | 176 | 2 | 18 | 6 |
| $\mathbf C$ | 68 | 2 | 19 | 5 |
| $\mathbf D$ | 24 | 24 | 29 | 7 |
| $\mathbf E$ | 71 | 1 | 19 | 6 |
| $\mathbf F$ | 87 | 1 | 20 | 6 |
| $\alpha$ | 91 | 1 | 18 | 5 |
| $\gamma$ | 24 | 24 | 24 | 7 |
| $\delta$ | 306 | 2 | 21 | 4 |
| $\varepsilon$ | 3 | 1 | 16 | 5 |
| $\zeta$ | 172 | 2 | 20 | 4 |
| $\eta$ | 34 | 4 | 23 | 6 |

[verified exact, $n\le3000$; `t1b_size.gp`, `t1c_defect.gp`.]
So $R_n$ is a *genericity upper bound*, not the denominator: the sum
$\sum_m c(m)m^{-r}e_{n,m}$ does cancel, at essentially every $n$. But the cancellation is
sharply localised:

> **Observation 1 (bad primes are exact).** For every row and every bad prime $p\mid c$,
> $$v_p\bigl(\operatorname{den}b_n\bigr)=v_p(R_n)\qquad(1\le n\le3000),$$
> with the single exception of $\varepsilon$ at $p=2$. Equivalently: **$b_n$ is $p$-integral at
> every bad prime**, except that $v_2(b^{\varepsilon}_n)=-1$ exactly when $n$ is a power of $2$
> (and $v_2(R_n)=1$ for all $n\ge2$, so $R_n$ over-predicts by one factor of $2$ at the other $n$).
> Aggregated over $n\le3000$, the total defect exponent at bad primes is $0$ for eleven rows and
> $2988\,(=3000-1-11)$ for $\varepsilon$ at $p=2$. [verified exact, $n\le3000$; `t1e_badden.gp`, `t2d_extra.gp` (d)]

All the remaining discrepancy lives at **good** primes; the good-prime defect totals run from
$7952$ ($\zeta$) to $22481$ ($\mathbf D$).

> **Observation 2 (mechanism at the large primes).** For $n/2<p\le n$ only $m=p$ contributes a $p$
> to the sum, and
> $$v_p(\operatorname{den}b_n)=\max\bigl(0,\ r-v_p\bigl(c(p)\,e_{n,p}\bigr)\bigr).$$
> [verified exact, all twelve Eisenstein rows, all $3\le n\le120$ and all $n/2<p\le n$:
> $859$ cases per row, no exception; `t1f_mech.gp`.]
> So $R_n$ fails precisely when $p\mid e_{n,p}$, i.e. the defect is a property of the *nome*
> $q(t)$, not of the source. This is the exact statement that replaces "$\operatorname{den}=R_n$".

The list of good primes that deviate is essentially all of them: e.g. row $\mathbf A$ first
deviates at $n=7$ ($p=5$), and by $n=3000$ every good prime $p\le n$ has deviated for some $n$.
Small-prime defects begin at
$\mathbf A$: $(n,p,e)=(7,5,1)$; $\mathbf B$: $(10,7,1)$; $\mathbf C$: $(7,5,1)$;
$\mathbf D$: $(4,3,2)$; $\mathbf E$: $(7,5,1)$; $\mathbf F$: $(9,7,1)$;
$\alpha$: $(6,3,2)$; $\gamma$: $(6,3,2)$; $\delta$: $(20,5,1)$;
$\varepsilon$: $(3,2,1)$; $\zeta$: $(14,7,1)$; $\eta$: $(6,3,3)$.

### 1.4 Cooper's rows

$s_7,s_{10},s_{18}$ have no Eisenstein source, so $R_n$ is not defined. What the free integration buys:
$d_n^{2}b_n\in\mathbf Z$ (not $d_n^3$) for all $n\le3000$ [verified exact], and this is sharp —
$d_n^{1}b_n\notin\mathbf Z$ already at $n=2$ ($b_2=45/4,\ 21/4,\ 51/4$ resp., $d_2=2$). At $p=2$ the
saving is exactly one power of $\operatorname{lcm}$ per level: see §3, where
$v_2(b_n)=s_2(n)-2\lfloor\log_2n\rfloor-1$ is proved sharp on the range for $s_{10}$ and $s_{18}$.

### 1.5 The companion formula itself

For all twelve Eisenstein rows, with $(t,F)$ built from the recurrence alone:

* $q(t)\in\mathbf Z[[t]]$ to $t^{60}$ [verified exact];
* $e_{n,m}=[t^n](F\,q^m)\in\mathbf Z$ for all $1\le m\le n\le60$ [verified exact];
* $$b_n=\sum_{m=1}^{n}\frac{c(m)}{m^{r}}\,e_{n,m}\qquad\text{exactly, as rationals, for all }n\le60,$$
  for every one of the twelve rows, with no exception. [verified exact, $n\le60$; `t1d_formula.gp`]

This is an independent re-derivation of Theorem A's twelve-family source table from the
recurrences: nothing but $(a,b,c)$ and the $c(m)$ of Table `tab:sources` enters.

---

## 2. Task 2 — valuation laws at the bad primes

Let $s_p(n)$ be the base-$p$ digit sum. Write $\mathrm{DIF}(n):=v_p(b_n)-v_p(a_n)=v_p(b_n/a_n)$.

### 2.1 $v_p(a_n)$

| row | $p$ | $v_p(c)$ | sharp $\lambda^{*}=\max_n\frac{v_p(a_n)}{s_p(n)}$ | least integer $\lambda_p$ with $v_p(a_n)\le\lambda_ps_p(n)$ | range of $v_p(a_n)-\lambda_ps_p(n)$ | exact law? |
|---|---|---|---|---|---|---|
| $\mathbf A$ | 2 | 3 | $24/7$ | 4 | $[-21,-3]$ | no |
| $\mathbf B$ | 3 | 3 | $1$ | 1 | $[0,0]$ | **yes** |
| $\mathbf C$ | 3 | 2 | $13/7$ | 2 | $[-14,-1]$ | no |
| $\mathbf E$ | 2 | 5 | $2$ | 2 | $[0,0]$ | **yes** |
| $\mathbf F$ | 2 | 3 | $10/3$ | 4 | $[-21,-3]$ | no |
| $\mathbf F$ | 3 | 2 | $13/7$ | 2 | $[-14,-1]$ | no |
| $\alpha$ | 2 | 6 | $9/2$ | 5 | $[-22,-2]$ | no |
| $\delta$ | 3 | 4 | $12/5$ | 3 | $[-21,-2]$ | no |
| $\varepsilon$ | 2 | 4 | $3$ | 3 | $[-1,0]$ | **yes, with a correction** |
| $\zeta$ | 3 | 3 | $27/14$ | 2 | $[-7,-1]$ | no |
| $\eta$ | 5 | 3 | $11/10$ | 2 | $[-25,-1]$ | no |
| $s_7$ | 3 | 3 | $1/2$ | 1 | $[-7,-1]$ | no |
| $s_{10}$ | 2 | 6 | $1$ | 1 | $[0,0]$ | **yes** |
| $s_{18}$ | 2 | 6 | $1$ | 1 | $[0,0]$ | **yes** |
| $s_{18}$ | 3 | 1 | $13/6$ | 3 | $[-21,-2]$ | no |

[verified exact, $1\le n\le3000$; `t2_val.gp`, `t3_laws.gp`.]
The $\lambda_p$ column reproduces the paper's list
(`04_padic_euler_draft.tex`, Verification (d): $4$ for $\mathbf A,\mathbf F$ at $2$; $5$ for
$\alpha$; $3$ for $\varepsilon$; $2$ for $\mathbf C,\mathbf F$ at $3$, $\zeta$, $\eta$; $3$ for
$\delta,s_{18}$) **exactly, with no exception**, and adds the sharp rational $\lambda^*$, which is
strictly smaller than $\lambda_p$ in every non-exact case (e.g. $\mathbf A$: $24/7=3.43$, not $4$;
$\eta$: $11/10$, not $2$). Adding a linear term or a $\lfloor\log_pn\rfloor$ term does not close
any of the non-exact cases: the best fit
$v_p(a_n)=An+Bs_p(n)+C$ over $|A|,|B|\le8$ always has $A=0$ and residual width
$\ge6$ (widths: $\mathbf A$ 16, $\mathbf C$ 10, $\mathbf F$ 15/10, $\alpha$ 15, $\delta$ 11,
$\zeta$ 6, $\eta$ 7, $s_7$ 6, $s_{18}@3$ 9). [fit + verified exact, $n\le3000$]

### 2.2 $v_p(b_n)-v_p(a_n)$ is constant, and equals $v_p(\xi_p)$

> **Observation 3.** For every $(\text{row},p)$ at which the Euler-factor criterion holds *with
> $\xi_p\neq0$*, $\mathrm{DIF}(n)$ is **constant on the whole range** $1\le n\le3000$ (at worst with
> one exceptional small $n$), and the constant is $v_p(\xi_p)$. The two rows with $\xi_p=0$
> ($\mathbf A$ at $2$, $\zeta$ at $3$) are treated below.

| row | $p$ | $\mathrm{DIF}(n)$ | $\#\{n\le3000:\mathrm{DIF}\neq\text{const}\}$ | $\xi_p$ (Theorem F / `EULER_CRITERION.md`) | $v_p(\xi_p)$ measured |
|---|---|---|---|---|---|
| $\mathbf B$ | 3 | $-1$ | $0$ | $\tfrac12\zeta_3(2)$ | $-1$ |
| $\mathbf C$ | 3 | $-1$ | $0$ | $\tfrac12\zeta_3(2)$ | $-1$ |
| $\mathbf E$ | 2 | $-2$ | $0$ | $\tfrac12\zeta_2(2)$ | $-2$ |
| $\mathbf F$ | 2 | $0$ | $1$ ($n=1$) | $\tfrac12L_2(2,\chi_{12})$ | $0$ |
| $\mathbf F$ | 3 | $-1$ | $0$ | $\tfrac58\zeta_3(2)$ | $-1$ |
| $\alpha$ | 2 | $-2$ | $0$ | $\tfrac13\zeta_2(3)$ | $-2$ |
| $\delta$ | 3 | $-1$ | $0$ | $\tfrac14\zeta_3(3)$ | $-1$ |
| $\varepsilon$ | 2 | $-4$ | $1$ ($n=1$) | $\tfrac14\zeta_2(3)$ | $-4$ |
| $\eta$ | 5 | $-1$ | $0$ | $\tfrac12\zeta_5(3)$ | $-1$ |
| $s_{18}$ | 3 | $-1$ | $1$ ($n=2$) | $\tfrac12\zeta_3(2)$ | $-1$ |

[verified exact, $1\le n\le3000$; `t2_val.gp`, `t2c_xi.gp`.]
Two internal consistency checks fall out for free and both hold:
$v_2(\zeta_2(3))=-2$ read from $\alpha$ ($\xi_2=\tfrac13\zeta_2(3)$, $v_2=-2$) agrees with
$\varepsilon$ ($\xi_2=\tfrac14\zeta_2(3)$, $v_2=-4$); and
$v_3(\zeta_3(2))=-1$ is read identically off $\mathbf B$, $\mathbf C$, $\mathbf F$ and $s_{18}$.

**The two zero-limit rows.** For $\mathbf A$ at $p=2$ and $\zeta$ at $p=3$ — exactly the two rows
whose co-divisor character $\varphi$ is non-trivial — $\mathrm{DIF}$ is *not* constant but grows:

$$\mathrm{DIF}(n)-3n\in[-54,6]\ (\mathbf A,\,p=2),\qquad
\mathrm{DIF}(n)-3n\in[-49,3]\ (\zeta,\,p=3),$$

i.e. $v_p(b_n/a_n)=\sigma_pn+O(1)$ with $\sigma_p=v_p(c)=3$, and $\xi_p=0$. Measured:
$v_2(\xi^{\mathbf A})=8972$ and $v_3(\xi^{\zeta})=8990$ at $N=3000$ (Cauchy precision $8954$, $8984$),
i.e. $\xi_p=0$ to $\ge p^{8954}$. [verified exact, $n\le3000$]

**The three non-convergent Cooper cases.** $s_7$ at $p=3$ and $s_{10},s_{18}$ at $p=2$ have
$p\mid c$ but no $p$-adic limit — $\mathrm{DIF}$ drifts to $-\infty$ logarithmically:
$$\mathrm{DIF}(n)+2\lfloor\log_pn\rfloor\in\begin{cases}\{-1\}&s_{10},s_{18}\text{ at }p=2\ (\text{exact}),\\ [0,1]&s_7\text{ at }p=3.\end{cases}$$
[verified exact, $n\le3000$.] So the "slope $0$ at $p=2,3$" recorded for Cooper's rows in
`ROW_LEDGER.md` / `SLOPE_CENSUS.md` is here made exact: the ratio $b_n/a_n$ loses exactly
$2$ digits per power of $p$, which is the free integration acting in reverse.

### 2.3 $\xi_p$: digits and Conjecture-D ratios

At $N=3000$, $\xi_p\approx b_N/a_N$:

| row | $p$ | $v_p(\xi)$ | Cauchy $v_p(\xi_N-\xi_{N-1})$ | $p^{-v}\xi \bmod p^8$ |
|---|---|---|---|---|
| $\mathbf A$ | 2 | $8972$ | $8954$ | $61$ |
| $\mathbf B$ | 3 | $-1$ | $8986$ | $2893$ |
| $\mathbf C$ | 3 | $-1$ | $5988$ | $2893$ |
| $\mathbf E$ | 2 | $-2$ | $14957$ | $219$ |
| $\mathbf F$ | 2 | $0$ | $8954$ | $179$ |
| $\mathbf F$ | 3 | $-1$ | $5988$ | $1976$ |
| $\alpha$ | 2 | $-2$ | $17932$ | $47$ |
| $\delta$ | 3 | $-1$ | $11984$ | $5770$ |
| $\varepsilon$ | 2 | $-4$ | $11940$ | $141$ |
| $\zeta$ | 3 | $8990$ | $8984$ | $2641$ |
| $\eta$ | 5 | $-1$ | $8968$ | $301776$ |
| $s_{18}$ | 3 | $-1$ | $2992$ | $2893$ |

Exact rational combinations, evaluated on the exact $b_N/a_N$:

$$v_3(\xi^{\mathbf B}-\xi^{\mathbf C})=5991,\quad
v_3(\xi^{\mathbf B}-\xi^{s_{18}})=2994,\quad
v_3(4\xi^{\mathbf F}-5\xi^{\mathbf B})=5991,\quad
v_2(3\xi^{\alpha}-4\xi^{\varepsilon})=11956 .$$

Each *exceeds* the Cauchy precision of the weaker row of the pair ($5988$, $2992$, $5988$, $11940$
respectively), so all four Conjecture-D ratios ($1$, $1$, $5/4$, $4/3$) are confirmed **to the full
precision the rows make available** [verified exact, $N=3000$].
$\mathbf B,\mathbf C,s_{18}$ give literally the same eight leading $3$-adic digits ($2893$).

---

## 3. Task 3 — the exact digit laws, and their consequences for $b_n$

Scanning all $15$ rows against all primes $p\le43$ for an exact law
$v_p(a_n)=\lambda s_p(n)+\mu$ over $1\le n\le3000$ turns up exactly **four** with $\lambda>0$
(the $\lambda=0$ hits are the different phenomenon "$p\nmid a_n$ for all $n\le3000$", which happens
for many good $p$ — e.g. $\mathbf B$ is never even, $\alpha$ is never divisible by $3$), plus one
near-law:

$$\boxed{\;v_3(a^{\mathbf B}_n)=s_3(n),\qquad v_2(a^{\mathbf E}_n)=2s_2(n),\qquad
v_2(a^{s_{10}}_n)=v_2(a^{s_{18}}_n)=s_2(n)\;}$$
$$v_2(a^{\varepsilon}_n)=3s_2(n)-[\,n\text{ odd}\,]\quad(\text{residual determined by }n\bmod 2).$$

[verified exact, $1\le n\le10^4$ for all five; `t3c_scan.gp`, `t3b_verify.gp`.]
The repo records only $\mathbf B$ at $3$, $\mathbf E$ at $2$ (`04_padic_euler_draft.tex`) and
$\varepsilon$ at $2$ (`05_two_row.tex` line 428, for the row it calls $\mathbf T$); the two
**Cooper laws are new here**, and they are *identical*: not merely the same law but the same
vector — $v_2(a^{s_{10}}_n)=v_2(a^{s_{18}}_n)$ **and** $v_2(b^{s_{10}}_n)=v_2(b^{s_{18}}_n)$ for
all $n\le3000$ [verified exact], matching the "identical vectors" already noticed for
$v_2(\operatorname{den}[q^m]\sqrt F)$ in `WEIGHT_DROP.md` §1.3.

**The implied laws for $b_n$, verified directly.** Combining with §2.2:

| row | $p$ | law for $b_n$ | verified |
|---|---|---|---|
| $\mathbf B$ | 3 | $v_3(b_n)=s_3(n)-1$ | exact, $1\le n\le6000$ |
| $\mathbf E$ | 2 | $v_2(b_n)=2s_2(n)-2$ | exact, $1\le n\le6000$ |
| $\varepsilon$ | 2 | $v_2(b_n)=3s_2(n)-[n\text{ odd}]-4$ | exact, $2\le n\le6000$ (fails only at $n=1$) |
| $s_{10}$ | 2 | $v_2(b_n)=s_2(n)-2\lfloor\log_2n\rfloor-1$ | exact, $1\le n\le6000$ |
| $s_{18}$ | 2 | $v_2(b_n)=s_2(n)-2\lfloor\log_2n\rfloor-1$ | exact, $1\le n\le6000$ |

[verified exact; `t3b_verify.gp`.] Three corollaries worth recording:

* $b^{\mathbf B}_n$ and $b^{\mathbf E}_n$ are $p$-integral for every $n$, with $v_p(b_n)=0$
  exactly on the $p$-power indices — the denominator of the companion is *prime to $c$*.
* $b^{\varepsilon}_n$ has $v_2=-1$ **iff** $n=2^k$ ($k\ge1$): $3s_2(n)-[n\text{ odd}]-4<0$
  forces $s_2(n)=1$, $n$ even. This is the single exception in Observation 1, now explained.
* For $s_{10},s_{18}$: $v_2(d_n^2b_n)=s_2(n)-1\ge0$ with equality exactly at $n=2^k$ — so
  $d_n^2b_n\in\mathbf Z_2$ is sharp at $p=2$ and $d_n^{1}b_n$ is not. The free integration is
  visible as the coefficient $2$ (not $3$) of $\lfloor\log_2n\rfloor$.

---

## 4. Task 4 — good primes

Here $\gamma$ ($c=1$) and $\mathbf D$ ($c=-1$) have **no** bad primes at all, so every $p$ is good.

### 4.1 Towers: the raw congruence fails, the $p^{w}$-twisted one holds with exponent $3$ per level

Write $w=r$ and $r_n=b_n/a_n$. Measured $v_p(r_{np^{k+1}}-r_{np^k})$ for $n=1,2,3$, $k=0,1,2$:

$$v_p\bigl(r_{np^{k+1}}-r_{np^k}\bigr)=-w(k+1)+v_p(r_n)\quad\text{(negative, decreasing)} .$$

So the *literal* tower congruence $r_{np^{k+1}}\equiv r_{np^k}\pmod{p^{e(k)}}$ **fails**, and fails by
exactly $p^{w(k+1)}$: e.g. $\gamma$ at $p=11$, $n=1$: $e=-3,-6,-9$; $\mathbf D$ at $p=13$, $n=1$:
$e=-2,-4,-6$. [verified exact, $n p^{k+1}\le6600$; `t4a_tower.gp`]

The correct normalisation is $\Lambda_k:=p^{wk}\,b_{np^k}/a_{np^k}$. (This is the tower object
$\Lambda_a=\lim_s\psi(p)^sp^{ws}b_{ap^s}/a_{ap^s}$ of `GOOD_PRIME_TOWERS.md` §2; the two agree
here because $\psi=\mathbf 1$ for both $\gamma$ and $\mathbf D$ — see §4.2 — so the $\psi(p)^s$
factor is absent.) Then

$$\Lambda_{k+1}\equiv\Lambda_k\pmod{p^{\,e(k)}\Lambda_k},\qquad
e(k)=\min\Bigl(v_p(\rho^a_k-1),\ v_p(p^{w}\rho^b_k-1)\Bigr),$$
$\rho^x_k=x_{np^{k+1}}/x_{np^k}$, and the measurement is

| cell | $e(0),e(1),e(2),\dots$ |
|---|---|
| $\gamma$, $p=5$: $n=1,3$ / $n=2$ | $2,5,8,11$ / $3,6,9,12$ |
| $\gamma$, $p=7$: $n=1,2$ / $n=3$ | $3,6,9$ / $2,5,8$ |
| $\gamma$, $p=11$: $n=1,2,3$ | $3,6,9$ |
| $\gamma$, $p=13$: $n=1,3$ / $n=2$ | $3,6,9$ / $2,5,8$ |
| $\mathbf D$, $p=5$: $n=1,3$ / $n=2$ | $3,6,9,12$ / $1,4,7,10$ |
| $\mathbf D$, $p=7$: $n=1,2$ / $n=3$ | $3,6,9$ / $1,4,7$ |
| $\mathbf D$, $p=11,13$: $n=1,2,3$ | $3,6,9$ |

[verified exact, $np^{k+1}\le6600$; `t4a2_tower.gp`.] So

> **$e(k)=3(k+1)$** in $18$ of the $24$ cells, and $e(k)=3(k+1)-|v_p(\Lambda_0)|$ in the other six
> (the cells where $\Lambda_0$ is not a $p$-adic unit). The gain is $\mathbf 3$ digits per level in
> **both** rows, although $w=3$ for $\gamma$ and $w=2$ for $\mathbf D$ — the exponent is the
> Beukers–Coster $3$, not $w$. Both constituents show it separately:
> $v_p(\rho^a_k-1)$ and $v_p(p^w\rho^b_k-1)$ each grow by exactly $3$ per level.

### 4.2 Lucas-type congruences

$N_n:=d_n^{\,r}b_n\in\mathbf Z$ does **not** satisfy $N_{np+m}\equiv N_nN_m\pmod p$: the hit rate is
$571/795$, $112/791$, $120/781$, $258/780$ for $\gamma$ at $p=5,7,11,13$ and similar for
$\mathbf D$ — i.e. essentially chance. Nor is there an ASD three-term law for $N_n$: the ratio
$N_{np}/N_n \bmod p$ wanders over $\mathbf F_p$ (e.g. $\gamma$, $p=13$: $12,-,8,8,1,-,12,12,5,\dots$),
so no $\alpha$ exists, and $v_p(N_{np}-N_n)$ is $0$ for most $n$.
[verified exact, $n\le800$; `t4b_lucas.gp`.] Two controls do hold:

* $a_{np+m}\equiv a_na_m\pmod p$ for **all** $15$ rows and all $p\in\{2,3,5,7,11,13\}$ and all
  $(n,m)$ with $np+m\le900$: $898/898$, $897/897$, $895/895$, $889/889$, $880/880$, $884/884$,
  no exception. (This is the Malik–Straub Lucas property; recomputed here as a control.)
* Beukers–Coster: $v_p(a_{np}-a_n)\ge3$ for every $n\le24$ and $p\in\{5,7,11,13\}$ in both rows,
  with equality most of the time.

**What does work — a Lucas law for the companion.** Put

$$\beta_n:=p^{\,k_p\lfloor\log_pn\rfloor}\,b_n\in\mathbf Z_p\quad(n\ge1),\qquad \beta_0:=0,$$

with $k_p=r$ ($=2$ or $3$) for the twelve Eisenstein rows and $k_p=2$ for Cooper's three (the
free-integration exponent). Then, at every **good** prime,

$$\boxed{\ \beta_{np+m}\;\equiv\;\psi(p)\,\beta_n\,a_m\pmod p,\qquad n\ge1,\ 0\le m<p,\ }$$

with $\psi$ the character of the *first* $L$-factor of $L(\Phi,s)$ (the $\psi$ of
`EULER_CRITERION.md` §1). Verified as a search over $u\in\mathbf F_p$ — the solution set is a
**single** value $u$ in every cell, and $u=\psi(p)$ in every cell:

| row | $\psi$ | $u$ at good $p=2,3,5,7,11,13,17,19,23,29,31,37,41,43$ |
|---|---|---|
| $\mathbf A,\ \alpha,\ \gamma,\ \delta,\ \varepsilon$ | $\mathbf 1$ | $1$ at every good $p$ |
| $\mathbf D$ | $\mathbf 1$ | $1$ at all fourteen |
| $\mathbf B,\mathbf C,\zeta$ | $\chi_{-3}$ | $1,-,4,1,10,1,16,1,22,28,1,1,40,1$ $=\chi_{-3}(p)\bmod p$ |
| $\mathbf F$ | $\chi_{-3}$ | $-,-,4,1,10,1,16,1,22,28,1,1,40,1$ (both $2,3$ bad) |
| $\mathbf E$ | $\chi_{-4}$ | $-,2,1,6,10,1,1,18,22,1,30,1,1,42$ $=\chi_{-4}(p)\bmod p$ |
| $\eta$ | $\chi_5$ | $1,2,-,6,1,12,16,1,22,1,1,36,1,42$ $=\chi_5(p)\bmod p$ |

[verified exact; $15$ rows $\times$ up to $14$ primes, all $(n,m)$ with $np+m\le2000$, $100\%$ of
$\approx1{,}500$–$2{,}000$ pairs per cell, unique solution $u$; `t4e_psi.gp`, `t4f_cooper.gp`,
`t4g_full.gp`.]

Two things follow.

1. **This is a mod-$p$ refinement of Observation `obs:tower`** of `04_padic_goodprimes.tex`:
   the tower statement $\rho^b_s\to\psi(p)p^{-w}$ is the case $m=0$, $n=ap^s$; the congruence above
   extends it from the towers $n=ap^s$ to *all* $n$, and to a genuine Lucas factorisation. Iterating
   over the base-$p$ digits $n=(n_j\cdots n_0)_p$ and using the $a$-Lucas law gives
   $$\beta_n\equiv\psi(p)^{\,j}\,\beta_{n_j}\prod_{i<j}a_{n_i}\pmod p,
   \qquad\text{so}\qquad \frac{\beta_n}{a_n}\equiv\psi(p)^{\,j}\frac{b_{n_j}}{a_{n_j}}\pmod p:$$
   modulo $p$, the normalised ratio depends only on the *leading* base-$p$ digit and the number of
   digits. **The iterated form was checked directly** and holds for every $1\le n\le2000$ in all
   $199$ (row, prime) cells listed above, with no exception [verified exact; `t4h_digit.gp`].
2. **It determines $\psi$ for Cooper's three rows**, whose sources are meromorphic and
   unidentified (Remark `rem:cooper`). The measurement gives, unambiguously,
   $$\psi_{s_7}=\mathbf 1,\qquad \psi_{s_{10}}=\mathbf 1,\qquad
   \psi_{s_{18}}=\chi_{-3},$$
   the last from $u=1,0,4,1,10,1,16,1,22,28,1,1,40,1$ at $p=2,3,5,7,\dots,43$: twelve unramified
   determinations at $5\le p\le43$, every one matching $\chi_{-3}(p)$, plus the ramified value
   $u=0=\chi_{-3}(3)$ at $p=3$ (at $p=2$ the two classes coincide, so that cell is uninformative).
   For $s_7$ and $s_{10}$ the normalisation $k_p=2$ works at *every* prime including their
   nominally bad ones ($3$ and $2$), giving $u=1$ there too — which is the exact form of the
   observation in `ROW_LEDGER.md` that these two rows carry **no** $p$-adic slope anywhere: at
   $p=3$ (resp. $2$) they behave like a good prime. This confirms as an exact computation
   the entries $(\psi,\varphi)=(\mathbf 1,\ast)$ for $s_7,s_{10}$ and the *conjectural*
   $(\chi_{-3},\mathbf 1)$ for $s_{18}$ recorded in `WEIGHT_DROP.md` §V2 and
   `EULER_CRITERION.md` §4.1 (where the $s_{18}$ line is flagged "conj.").

At the **bad** primes the same search on the unnormalised $b_n$ gives $u=1$ for
$\mathbf B,\mathbf C$ at $3$, $\mathbf E,\alpha$ at $2$, $\delta$ at $3$, $\eta$ at $5$, $\mathbf F$
at $3$ — the unipotent case — and $u=0$ for $\mathbf A$ at $2$ and $\zeta$ at $3$ — the two rows with
$\xi_p=0$. ($\mathbf F$ and $\varepsilon$ at $p=2$ also return $u=0$, but $p=2$ leaves only
$m\in\{0,1\}$ and the value is forced by $n=1$ alone, so this cell carries no information.)
This reproduces the three-case trichotomy of `obs:tower` mod $p$.

---

## 5. Task 5 — summary table

$k$ = least exponent with $d_n^{\,k}b_n\in\mathbf Z$; $\lambda_p$ = least integer with
$v_p(a_n)\le\lambda_ps_p(n)$; all entries [verified exact, $n\le3000$] unless a wider range is
stated.

| row | $r$ | $k$ | bad $p$ | $\operatorname{den}(b_n)=R_n$? | digit law for $a_n$ | $v_p(b_n)-v_p(a_n)$ | $v_p(\xi_p)$ | $\psi$ |
|---|---|---|---|---|---|---|---|---|
| $\mathbf A$ | 2 | 2 | $2$ | no (good $p$ only; exact at $p=2$) | $\lambda_2=4$, sharp $24/7$ | $3n+O(1)$, unbounded | $\xi_2=0$ ($\ge2^{8954}$) | $\mathbf 1$ |
| $\mathbf B$ | 2 | 2 | $3$ | no (good $p$ only; exact at $p=3$) | **$v_3(a_n)=s_3(n)$** ($n\le10^4$) | $-1$ (all $n$) | $-1$ | $\chi_{-3}$ |
| $\mathbf C$ | 2 | 2 | $3$ | no (good $p$ only; exact at $p=3$) | $\lambda_3=2$, sharp $13/7$ | $-1$ (all $n$) | $-1$ | $\chi_{-3}$ |
| $\mathbf D$ | 2 | 2 | — | no ($R_n=d_n^2$ identically) | — | — | — | $\mathbf 1$ |
| $\mathbf E$ | 2 | 2 | $2$ | no (good $p$ only; exact at $p=2$) | **$v_2(a_n)=2s_2(n)$** ($n\le10^4$) | $-2$ (all $n$) | $-2$ | $\chi_{-4}$ |
| $\mathbf F$ | 2 | 2 | $2,3$ | no (good $p$ only; exact at $2,3$) | $\lambda_2=4$ ($10/3$), $\lambda_3=2$ ($13/7$) | $0$ at $p{=}2$; $-1$ at $p{=}3$ | $0$; $-1$ | $\chi_{-3}$ |
| $\alpha$ Domb | 3 | 3 | $2$ | no (good $p$ only; exact at $p=2$) | $\lambda_2=5$, sharp $9/2$ | $-2$ (all $n$) | $-2$ | $\mathbf 1$ |
| $\gamma$ Apéry | 3 | 3 | — | no ($R_n=d_n^3$ identically) | — | — | — | $\mathbf 1$ |
| $\delta$ | 3 | 3 | $3$ | no (good $p$ only; exact at $p=3$) | $\lambda_3=3$, sharp $12/5$ | $-1$ (all $n$) | $-1$ | $\mathbf 1$ |
| $\varepsilon$ ($\mathbf T$) | 3 | 3 | $2$ | no; **and $R_n$ over-predicts $2^1$** except at $n=2^k$ | **$v_2(a_n)=3s_2(n)-[n\text{ odd}]$** ($n\le10^4$) | $-4$ ($n\ge2$) | $-4$ | $\mathbf 1$ |
| $\zeta$ | 3 | 3 | $3$ | no (good $p$ only; exact at $p=3$) | $\lambda_3=2$, sharp $27/14$ | $3n+O(1)$, unbounded | $\xi_3=0$ ($\ge3^{8984}$) | $\chi_{-3}$ |
| $\eta$ | 3 | 3 | $5$ | no (good $p$ only; exact at $p=5$) | $\lambda_5=2$, sharp $11/10$ | $-1$ (all $n$) | $-1$ | $\chi_5$ |
| $s_7$ | 3 | **2** | $3$ | n/a (no Eisenstein source) | $\lambda_3=1$, sharp $1/2$ | $-2\lfloor\log_3n\rfloor+[0,1]$ | no limit | $\mathbf 1$ (confirmed) |
| $s_{10}$ | 3 | **2** | $2$ | n/a | **$v_2(a_n)=s_2(n)$** ($n\le10^4$) | $-2\lfloor\log_2n\rfloor-1$ (exact) | no limit | $\mathbf 1$ (confirmed) |
| $s_{18}$ | 3 | **2** | $2,3$ | n/a | **$v_2(a_n)=s_2(n)$** ($n\le10^4$); $\lambda_3=3$ ($13/6$) | $p{=}2$: $-2\lfloor\log_2n\rfloor-1$; $p{=}3$: $-1$ | $p{=}2$: none; $p{=}3$: $-1$ | $\chi_{-3}$ (conj. confirmed) |

---

## 6. What is new here relative to the repository

1. **$R_n$ is not the denominator.** The refined prediction is an upper bound only, and the
   discrepancy is entirely at good primes (Observation 1) with the exact mechanism $p\mid e_{n,p}$
   (Observation 2). The two "Apéry-perfect" rows $\mathbf D,\gamma$ are exactly the rows where
   $R_n$ collapses back to $d_n^{\,r}$.
2. **Two new exact digit laws**: $v_2(a^{s_{10}}_n)=v_2(a^{s_{18}}_n)=s_2(n)$, together with the
   identity of the two full valuation vectors for $a$ *and* $b$.
3. **Exact companion laws** for all five rows that have an exact $a$-law, including
   $v_2(b_n)=s_2(n)-2\lfloor\log_2n\rfloor-1$ for $s_{10},s_{18}$, which is the free integration
   made exact at $p=2$.
4. **A Lucas congruence for the companion**, $\beta_{np+m}\equiv\psi(p)\beta_na_m\pmod p$, which
   upgrades the tower Observation from towers to all $n$ — and which **determines $\psi$ for
   Cooper's three rows** ($\mathbf 1,\mathbf 1,\chi_{-3}$), confirming the conjectural $s_{18}$
   entry of `EULER_CRITERION.md` §4.1 by an exact computation that never touches the (meromorphic,
   unidentified) source.
5. The good-prime tower gain is **$3$ digits per level in both $\gamma$ and $\mathbf D$**, i.e.
   the Beukers–Coster exponent, not $w$; the raw tower congruence on $b_n/a_n$ fails by exactly
   $p^{w(k+1)}$.

## 7. Open / not settled here

* No exact law was found for $v_p(a_n)$ at the ten remaining (row, bad prime) cells; the residuals
  $v_p(a_n)-\lambda_ps_p(n)$ have width $6$–$16$ and are not functions of $n\bmod p^j$, $j\le4$.
  Whether a carry-counting law exists is open.
* The $\beta$-Lucas law is a congruence mod $p$ and no better: $\min_{n,m}
  v_p(\beta_{np+m}-\psi(p)\beta_na_m)=1$ in $187$ of the $199$ cells and $2$ in the other twelve
  ($\mathbf A@3$, $\mathbf C@2$, $\mathbf D@5$, $\alpha@3$, $\gamma@2$, $\gamma@3$, $\zeta@2$,
  $\eta@3$, $s_7@3$, $s_{10}@2$, $s_{18}@2$, $s_{18}@3$), so it is not a supercongruence
  [verified exact, $np+m\le2000$].
  The right mod-$p^2$ correction term is not identified.
* $e(k)=3(k+1)-|v_p(\Lambda_0)|$ in the six non-unit cells is a fitted description of $24$ cells,
  not a derivation.
* Observation 2 was checked to $n\le120$ (it needs the $e_{n,m}$ matrix); pushing it to $n\le3000$
  is a matter of series arithmetic time, not of new ideas.

## 8. Files

| file | purpose |
|---|---|
| `lib.gp` | the fifteen rows, exact $(a_n,b_n)$ generator, integer-numerator generator |
| `src.gp` | the twelve source coefficient functions $c(m)$ from Table `tab:sources` |
| `t0_check.gp`, `t0_int.gp`, `chk_zeta.gp` | integrality / first-terms sanity checks ($n\le3000$); the $c=+27$ non-integrality check |
| `t1a_den.gp`, `t1b_size.gp`, `t1c_defect.gp` | Task 1: exponent $k$, $R_n$ vs $d_n^r$ vs $\operatorname{den}$, defect anatomy |
| `t1d_formula.gp` | Task 1: the companion formula $b_n=\sum c(m)m^{-r}e_{n,m}$, $n\le60$ |
| `t1e_badden.gp`, `t1f_mech.gp`, `t1g_sharp.gp` | Task 1: bad-prime exactness, large-prime mechanism, sharpness |
| `t2_val.gp` (`t2_val.log`), `t2c_xi.gp`, `t2d_extra.gp` | Task 2: valuations, $\xi_p$, Conjecture-D ratios |
| `t3_laws.gp` (`t3_laws.log`), `t3b_verify.gp`, `t3c_scan.gp` | Task 3: digit-law search and verification to $n\le10^4$ |
| `t4a_tower.gp`, `t4a2_tower.gp` | Task 4: tower congruences |
| `t4b_lucas.gp`, `t4c_lucas2.gp`, `t4d_all.gp`, `t4e_psi.gp`, `t4f_cooper.gp`, `t4g_full.gp`, `t4h_digit.gp` | Task 4: Lucas / ASD searches, the $\psi(p)$ law, its iterated digit form |
