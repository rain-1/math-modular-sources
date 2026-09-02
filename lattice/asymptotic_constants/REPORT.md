# Asymptotic constants of modular Apéry sequences: a computational census

*Agent report for `consolidation/ASYMPTOTIC_CONSTANTS.md` §4.  Working directory
`scratchpad/asymK/`; PARI/GP 2.15.4.  Every claim is tagged
`[exact]` (identity of power series / rationals, no floating point),
`[proved]` (short argument given here), or `[numerical, d digits]`.
Definitions taken from the repository (row recurrences, the two $\zeta(5)$
sources and the $\zeta(7)$ source, and the three hosts of §3) are cited where
used; every **number** below is recomputed here, and the one repository value
used as a check is marked `[repo]`.*

---

## 0. Verdict

Ten rows were run.  **Seven are elementary, three are not**, and the split is
*not* the one the note's headline suggests.  Two mechanisms act independently:

| | $F(\tau_c)\neq0$ (i.e. $i^k\varepsilon=+1$) | $F(\tau_c)=0$ (i.e. $i^k\varepsilon=-1$) |
|---|---|---|
| **$k=2$** | $K\in\overline{\mathbf Q}\,\pi^{-3/2}$ — **elementary** (7 rows) | (does not occur here) |
| **$k>2$** | $K\in\overline{\mathbf Q}\,\Omega^{k-2}\pi^{-3/2}$ (2 rows) | $K\in\overline{\mathbf Q}\,\Omega^{k}\pi^{-1/2}$ (1 row) |

with $\Omega^2$ the weight-two Chowla–Selberg unit of the CM field of $\tau_c$.
So the correct statement is

> **$K$ is elementary iff the companion has weight $2$ *and* does not vanish at
> the fold.**  Vanishing at the fold costs $\Omega^{2}$; weight above $2$ costs
> $\Omega^{k-2}$.  Both costs are visible in the census.

For the seven weight-two rows a *single closed formula* covers all of them
(§2.2), with no case distinction and no modular input beyond the row's own
characteristic roots:
$$\boxed{\;K=\frac{\sqrt N}{2\pi^{3/2}}\sqrt{\frac{\lambda_1}{\lambda_1-\lambda_2}}\;}$$
$\lambda_1,\lambda_2$ the characteristic roots of the recurrence.  In particular
$K^2\pi^3=\tfrac{N\lambda_1}{4(\lambda_1-\lambda_2)}$ has degree $\le2$ over
$\mathbf Q$.  Apéry's classical
$K=(1+\sqrt2)^2/(2^{9/4}\pi^{3/2})$ is the case $N=6$: **confirmed, exponents
unchanged** [numerical, 73 digits].

The three Cooper rows $s_7,s_{10},s_{18}$ **do** have a Fricke fold as their
dominant singularity (this had to be checked), and
they fall under the same formula, giving the three prettiest values in the
census,
$K=\tfrac{3\sqrt3}{4}\pi^{-3/2}$, $\sqrt2\,\pi^{-3/2}$, $3\sqrt2\,\pi^{-3/2}$.

---

## 1. Method

### 1.1 Building $(x,F)$ from the recurrence alone

For a row given by a three-term recurrence
$p_3(n)a_{n+1}=p_2(n)a_n+p_1(n)a_{n-1}$, $p_3(n)=(n+1)^3$, $a_0=1$, $a_{-1}=0$
(MUM: triple indicial root $0$), deform $n\mapsto n+\epsilon$, work in
$\mathbf Q[\epsilon]/(\epsilon^2)$, and set
$$y_0=\sum a_n(0)\,x^n,\qquad g=\sum a_n'(0)\,x^n,\qquad q=x\,e^{g/y_0},\qquad
x(q)=q\text{-reversion},\qquad F(q)=y_0(x(q)).$$
This is the canonical mirror map; it uses **no modular input**.  Implemented in
`lib.gp` (`mirror`).  For Apéry's row it returns
$x(q)=q-12q^2+66q^3-\dots=(\eta_1\eta_6/\eta_2\eta_3)^{12}$ and
$F=1+5q+13q^2+\dots=\eta_2^7\eta_3^7/(\eta_1^5\eta_6^5)$ [exact, checked to
$O(q^{300})$], which is the standard sanity check that the construction is the
right one.

### 1.2 Locating the fold

Newton on $x'(q)=0$ along $(0,1)$, started anywhere near $e^{-2\pi/\sqrt N}$.
For every one of the seven weight-two rows the root is
$$q_c=e^{-2\pi/\sqrt N},\qquad x(q_c)=x_+=1/\lambda_1,\qquad D^2x(\tau_c)\neq0 .$$
From the recurrence-built series this is [numerical, $\ge200$ digits] for six
rows and $\ge103$ digits for $s_{18}$ (whose $x(q)$ has a pole inside
$|q|<1$ — at $u=-7\pm4\sqrt3$ — so its $q$-expansion converges only for
$|q|\lesssim1/4$, barely outside $q_c=0.2274$).  Re-run from the eta quotient
$u$ of §2.1, whose coefficients grow subexponentially, all seven give
$u(\tau_c)=1/\sqrt C$ and $x(\tau_c)=x_+$ to $\ge400$ digits, $s_{18}$ included.
So the dominant singularity is a simple fold at the Fricke point in all seven
cases; in particular **Cooper's $s_7,s_{10},s_{18}$ qualify** — a point I did
not find recorded anywhere in the repository.

### 1.3 Two independent evaluations of $K$

* **Fold formula.**  $K^2=x_+\,(DF(\tau_c))^2/\bigl(2\pi(-D^2x(\tau_c))\bigr)$,
  evaluated from the $q$-expansions at $q_c$ ($\ge600$ terms for the sporadic
  rows, $1200$ from the eta quotients; tails $<10^{-340}$).
  Sign: along $q\uparrow q_c$ one has $q-q_c=-\sqrt{2x_+/(-x'')}\sqrt{1-x/x_+}$,
  so $h(x_+)=-F'(q_c)\sqrt{2x_+/(-x'')}$ and
  $a_n\sim +\frac{DF(\tau_c)}{2\sqrt\pi}\sqrt{\frac{2x_+}{-D^2x}}\,x_+^{-n}n^{-3/2}$;
  with $DF(\tau_c)>0$ all seven rows have $a_n>0$, as they must.
* **Richardson.**  $a_n$ from the exact recurrence to $n=6000$
  (integrality verified at $n=6000$ for all seven), then
  $K_n=|a_n|x_+^{\,n}n^{3/2}$ extrapolated in $1/n$ with $20$ nodes
  $n=6000,5850,\dots,3150$ at $400$-digit working precision.

The two agree to $46$–$73$ digits (table in §2.3); the limiting error is the
divergence of the $1/n$ asymptotic series, except for $s_{18}$ where the second
fold ($\lambda_2=12$, ratio $(3/4)^n$) also enters.

---

## 2. The seven weight-two rows

### 2.1 A uniform structure (this is the main structural finding)

Let $\lambda_1>|\lambda_2|$ be the characteristic roots and put
$$B=\tfrac{\lambda_1+\lambda_2}{2},\qquad C=\Bigl(\tfrac{\lambda_1-\lambda_2}{4}\Bigr)^2,
\qquad\text{so }\lambda_{1,2}=B\pm2\sqrt C .$$

> **Structure Theorem (verified).** There is a normalised degree-one Hauptmodul
> $u=q+O(q^2)$, an eta quotient with integral $q$-expansion, such that
> $$x=\frac{u}{1+Bu+Cu^2},\qquad F=\frac{Du}{u}=D\log u,\qquad
> P_3(x):=1-2Bx+(B^2-4C)x^2=\Bigl(\frac{1-Cu^2}{1+Bu+Cu^2}\Bigr)^{2},$$
> and $u\big|W_N=\dfrac{1}{Cu}$.

* The three displayed series identities hold **[exact, to $O(q^{300})$]** for all
  seven rows (`main3.gp`); $u$ has integral coefficients in every case.
* $u\,|\,W_N=1/(Cu)$ is **[proved]**: for $u=\prod_d\eta_d^{r_d}$ with
  $\sum_d r_d=0$, $\eta(-1/(N\tau))^{\ }$ gives
  $u\circ W_N=\bigl(\prod_d(N/d)^{r_d/2}\bigr)\prod_d\eta_{N/d}^{r_d}$, and for
  each row $\prod_d\eta_{N/d}^{r_d}=1/u$ and $\prod_d(N/d)^{r_d}=1/C^2$
  **exactly as rationals** (`alconst.gp`).  Numerically
  $u(W_N\tau)u(\tau)=1/C$ to $\ge90$ digits at a generic $\tau$.
* **Consequence [proved]:** $F=D\log u$ with $\log u\circ W_N=-\log u+\log(1/C)$
  gives $F(-1/(N\tau))=-N\tau^2F(\tau)$, i.e.
  $$\varepsilon_F=-1,\qquad i^k\varepsilon_F=i^2(-1)=+1 .$$
  So all seven are **Fricke-odd of weight two**: Lemma 2.1 case (a).  Measured:
  $F|_2W_N/F=-1$ to $\ge100$ digits, $x\circ W_N=x$ to $\ge100$ digits
  (`fricke.gp`), and $DF(\tau_c)=\frac{2\sqrt N}{4\pi}F(\tau_c)$ to
  $10^{-404}$ (`main4.gp`).

The eta quotients (found by exact linear algebra on $\log(u/q)$, `etaid.gp`;
each verified against the recurrence-built $u$ to $O(q^{60})$ and, through
$(x,F)$, against the recurrence to $n=40$ **[exact]**):

| row | $N$ | $u$ | $B$ | $C$ | $\lambda_1$ | $\lambda_2$ |
|---|---|---|---|---|---|---|
| $\alpha$ (Domb, $(10,4,64)$) | 12 | $\bigl(\eta_2\eta_3\eta_{12}/(\eta_1\eta_4\eta_6)\bigr)^4$ | 10 | 9 | 16 | 4 |
| $\gamma$ (Apéry $\zeta(3)$, $(17,5,1)$) | 6 | $\eta_2\eta_6^5/(\eta_1^5\eta_3)$ | 17 | 72 | $17+12\sqrt2$ | $17-12\sqrt2$ |
| $\varepsilon$ ($(12,4,16)$) | 8 | $\eta_2^2\eta_8^4/(\eta_1^4\eta_4^2)$ | 12 | 32 | $12+8\sqrt2$ | $12-8\sqrt2$ |
| $\zeta$ ($(9,3,-27)$) | 9 | $(\eta_9/\eta_1)^3$ | 9 | 27 | $9+6\sqrt3$ | $9-6\sqrt3$ |
| Cooper $s_7$ | 7 | $(\eta_7/\eta_1)^4$ | 13 | 49 | 27 | $-1$ |
| Cooper $s_{10}=\sum_k\binom nk^4$ | 10 | $\bigl(\eta_5\eta_{10}/(\eta_1\eta_2)\bigr)^2$ | 6 | 25 | 16 | $-4$ |
| Cooper $s_{18}$ | 18 | $\bigl(\eta_2\eta_3^2\eta_{18}/(\eta_1\eta_6^2\eta_9)\bigr)^6$ | 14 | 1 | 16 | 12 |

*(Note on $\alpha$: the level-$6$ eta quotient
$(\eta_2\eta_6/\eta_1\eta_3)^6=q+6q^2+21q^3+\cdots$ of `SPORADIC_SEARCH.md` is the
$q\mapsto-q$ twist of the coordinate used here:
$x_\alpha(-q)=-(\eta_2\eta_6/\eta_1\eta_3)^6(q)$ **[exact, $O(q^{9})$]**, i.e. it
parametrises the $(-1)^na_n$ normalisation of the same row.  The untwisted Domb
row $a_n=1,4,28,256,2716,\dots$ lives at level $12$, and that is where the fold
sits ($q_c=e^{-2\pi/\sqrt{12}}$), as `03_fricke_fold.gp` already recorded.)*

First terms $a_0,\dots,a_6$ (from the eta-quotient $(x,F)$, agreeing with the
recurrence for $n\le40$ **[exact]**): $\alpha$: $1,4,28,256,2716,31504,387136$;
$\gamma$: $1,5,73,1445,33001,819005,21460825$;
$\varepsilon$: $1,4,40,544,8536,145504,2618176$;
$\zeta$: $1,3,27,309,4059,57753,866349$;
$s_7$: $1,4,48,760,13840,273504,5703096$;
$s_{10}$: $1,2,18,164,1810,21252,263844$;
$s_{18}$: $1,6,54,564,6390,76356,948276$.

### 2.2 The closed formula

With $r(u)=u/(1+Bu+Cu^2)$ and $u_c:=u(\tau_c)=1/\sqrt C$ (the fixed point of
$u\mapsto1/(Cu)$) one gets **[proved from §2.1]**
$$\Bigl(\tfrac{F}{Du}\Bigr)(\tau_c)=\frac1{u_c}=\sqrt C,\qquad
x_+=r(u_c)=\frac1{B+2\sqrt C}=\frac1{\lambda_1},\qquad
x''(u_c)=r''(u_c)=-\frac{2C^{3/2}}{\lambda_1^{2}},\qquad
\frac{2x_+}{-r''(u_c)}=\frac{\lambda_1}{C^{3/2}},$$
so that
$$K=\frac{k\sqrt N}{8\pi^{3/2}}\Bigl|\tfrac{F}{Du}(\tau_c)\Bigr|
\sqrt{\frac{2x_+}{-r''(u_c)}}
=\frac{\sqrt N}{4\pi^{3/2}}\sqrt{\frac{\lambda_1}{\sqrt C}}
=\frac{\sqrt N}{2\pi^{3/2}}\sqrt{\frac{\lambda_1}{\lambda_1-\lambda_2}},
\qquad K^2\pi^3=\frac{N\lambda_1}{4(\lambda_1-\lambda_2)} .$$

The three pieces requested separately:

| row | $(F/Du)(\tau_c)=\sqrt C$ | $x''(h_c)=r''(u_c)$ | $2x_+/(-x''(h_c))$ |
|---|---|---|---|
| $\alpha$ | $3$ | $-27/128=-0.2109375$ | $16/27$ |
| $\gamma$ | $6\sqrt2$ | $-864\sqrt2/(577+408\sqrt2)=-1.0588227343\ldots$ | $(17+12\sqrt2)/(72\sqrt{72})$ |
| $\varepsilon$ | $4\sqrt2$ | $-256\sqrt2/(272+192\sqrt2)=-0.6660889655\ldots$ | $(12+8\sqrt2)/(32\sqrt{32})$ |
| $\zeta$ | $3\sqrt3$ | $-162\sqrt3/(189+108\sqrt3)=-0.7461339179\ldots$ | $(9+6\sqrt3)/(27\sqrt{27})$ |
| $s_7$ | $7$ | $-686/729=-0.9410150892\ldots$ | $27/343$ |
| $s_{10}$ | $5$ | $-250/256=-0.9765625$ | $16/125$ |
| $s_{18}$ | $1$ | $-1/128=-0.0078125$ | $16$ |

all confirmed against the $q$-series values to $\ge200$ digits (`main4.gp`).

### 2.3 Census table, weight two

$k=2$, $\varepsilon=-1$, $i^k\varepsilon=+1$, $F(\tau_c)\neq0$ for **all seven**.

| row | $N$ | $x_+$ | $K$ (30 digits) | closed form | Richardson vs closed form |
|---|---|---|---|---|---|
| $\alpha$ | 12 | $1/16$ | $0.359174244250333123378163967255$ | $\dfrac{2}{\pi^{3/2}}$ | $5.0\cdot10^{-60}$ |
| $\gamma$ | 6 | $17-12\sqrt2$ | $0.220043767112643037850689759810$ | $\dfrac{(1+\sqrt2)^2}{2^{9/4}\pi^{3/2}}$ | $4.6\cdot10^{-73}$ |
| $\varepsilon$ | 8 | $\tfrac{3-2\sqrt2}{4}$ | $0.257797308917838939757311495404$ | $\dfrac{\sqrt{4+3\sqrt2}}{2\,\pi^{3/2}}$ | $1.7\cdot10^{-68}$ |
| $\zeta$ | 9 | $\tfrac{2\sqrt3-3}{9}$ | $0.260201758994440992303217903862$ | $\dfrac{3(\sqrt2+\sqrt6)}{8\,\pi^{3/2}}$ | $3.0\cdot10^{-69}$ |
| $s_7$ | 7 | $1/27$ | $0.233290514929399000627453098124$ | $\dfrac{3\sqrt3}{4\,\pi^{3/2}}$ | $7.6\cdot10^{-71}$ |
| $s_{10}$ | 10 | $1/16$ | $0.253974543736963879143053219739$ | $\dfrac{\sqrt2}{\pi^{3/2}}$ | $1.7\cdot10^{-67}$ |
| $s_{18}$ | 18 | $1/16$ | $0.761923631210891637429159659216$ | $\dfrac{3\sqrt2}{\pi^{3/2}}$ | $3.8\cdot10^{-46}$ |

$K^2\pi^3$ values: $4$; $\tfrac{24+17\sqrt2}{32}$; $\tfrac{4+3\sqrt2}{4}$;
$\tfrac{9(2+\sqrt3)}{16}$; $\tfrac{27}{16}$; $2$; $18$.  Minimal polynomials of
$K\pi^{3/2}$ over $\mathbf Q$, in the same order:
$$x-2;\quad 512x^4-768x^2-1;\quad 8x^4-16x^2-1;\quad 256x^4-576x^2+81;\quad
16x^2-27;\quad x^2-2;\quad x^2-18 .$$
Every one is returned by `algdep(K\pi^{3/2},4)` on the Richardson value at
$45$-digit precision **[numerical, confirmed]**; the degree is $4$ exactly when
$\lambda_1\notin\mathbf Q$.

$\gamma$: the classical Apéry asymptotic $a_n\sim(1+\sqrt2)^{4n+2}/(2^{9/4}\pi^{3/2}n^{3/2})$
is **confirmed exactly as stated**: $x_+^{-n}=(1+\sqrt2)^{4n}$ and
$K=(1+\sqrt2)^2/2^{9/4}\cdot\pi^{-3/2}$ [numerical, 73 digits].

---

## 3. The three higher-weight rows

Here $F$ (or the companion $G$) has weight $k>2$, $x=r(h)$ for a Hauptmodul $h$
of $\Gamma_0(N)$, and
$$\frac{F(\tau_c)}{\sqrt{-D^2x(\tau_c)}}=\frac{(F/(Dh)^{k/2})(\tau_c)\,(Dh)^{k/2}}{|Dh|\sqrt{-r''(h_c)}}
\in\overline{\mathbf Q}\cdot(Dh(\tau_c))^{k/2-1},$$
so $K\in\overline{\mathbf Q}\,\Omega^{k-2}\pi^{-3/2}$ when $F(\tau_c)\neq0$, and
$K\in\overline{\mathbf Q}\,\Omega^{k}\pi^{-1/2}$ when $F(\tau_c)=0$
($\Omega^2:=$ the weight-two CM unit at $\tau_c$; $\Gamma(1/3)^6/\pi^4$ for
$\mathbf Q(\sqrt{-3})$, $\Gamma(1/4)^4/\pi^3$ for $\mathbf Q(i)$).

### 3.1 Level-$12$ $\zeta(5)$ — the control, CM case, reproduced

Data rebuilt from scratch: $h_{12}=\eta_1^3\eta_4\eta_6^2/(\eta_2^2\eta_3\eta_{12}^3)$,
$v:=1/h_{12}=q+3q^2+7q^3+\cdots$, $x=v/(1+7v+12v^2)$ **[exact, $O(q^{903})$]**
— i.e. this system is the $(B,C)=(7,12)$ instance of the §2.1 shape, with
$\lambda_{1,2}=7\pm4\sqrt3$; $E=(-9E_4(3\tau)+16E_4(4\tau))/7$,
$Q(x)=143x^3+189x^2+21x+7$, $G=7E/Q(x)=\sum c_nx^n$, $k=4$.

* $E|_4W_{12}=-E$ **[numerical, 230 digits]**, so $i^4\varepsilon=-1$ and
  $E(\tau_c)=G(\tau_c)=0$ **[numerical, $<10^{-230}$]** — the CM case.
* New exact structure: with $R:=v^2G/(Dv)^2$ (a modular function of $v$),
  $$R(v)=\frac{(1+3v)(1+4v)(1-12v^2)}{(1+6v)^2(1+2v)^2}\qquad[\text{exact},\ O(q^{502})],$$
  so $R(v_c)=0$ at $v_c=1/\sqrt{12}$ *visibly* (the factor $1-12v^2$), and
  $R'(v_c)=-3/2$ **[exact]**.  Hence $DG(\tau_c)=-18\,(Dv(\tau_c))^3$ and
  $$K=\frac{9\,(Dv(\tau_c))^2}{\sqrt\pi}\sqrt{\frac{7+4\sqrt3}{24\sqrt3}} .$$
* $\bigl(Dv(\tau_c)/\Omega^2\bigr)^3=\dfrac{3(45+26\sqrt3)}{2^{17}}$
  **[numerical, 156 digits]** ($\Omega^2=\Gamma(1/3)^6/\pi^4$); this is the
  repo's $(Dh/W^2)^3=-81(45+26\sqrt3)/2048$ transported through $v=1/h$,
  $h_c^2=12$.
* **Closed form (equivalent to, and simpler than, the published one):**
  $$K=9\Bigl(\frac{3(45+26\sqrt3)}{2^{17}}\Bigr)^{2/3}\sqrt{\frac{7+4\sqrt3}{24\sqrt3}}\;
  \frac{\Gamma(1/3)^{12}}{\pi^{17/2}} = 0.685371848490534405951430044881\ldots$$
  Agrees with the fold formula to **$\ge148$ digits** and with the published
  value `[repo: ZETA5_K_CLOSED_FORM.md]` to all $50$ published digits.
* Independent Richardson: $c_n$ regenerated from the $q$-series, checked against
  the order-$14$ recurrence for $n=14\ldots40$ **[exact, all residuals $0$]**,
  run to $n=4200$ ($c_{4200}\in\mathbf Z$, $\log_{10}|c_{4200}|=4798.76$),
  $20$-point extrapolation: agrees with the closed form to
  **$1.2\cdot10^{-62}$**.  $c_n<0$ for large $n$.
* $K\pi^{17/2}/\Gamma(1/3)^{12}=0.0843620623288518148561713194047\ldots$

### 3.2 Level-$16$ $\zeta(5)$ — weight $4$, Fricke-**even**, $F(\tau_c)\neq0$

$x=\eta_2\eta_{16}^2/(\eta_1^2\eta_8)$ (a genuine degree-one Hauptmodul of
$\Gamma_0(16)$; Ligozat: simple zero at $\infty$, simple pole at $0$),
$t=x/(8x^2+2x+1)$ — again the §2.1 shape with $(B,C)=(2,8)$,
$\lambda_{1,2}=2\pm4\sqrt2$, $t_+=(2\sqrt2-1)/14$, $\tau_c=i/4$,
$q_c=e^{-\pi/2}$; $\Phi_{16}=-\frac1{504}\bigl(E_6-85E_6(2\tau)+1428E_6(4\tau)
-5440E_6(8\tau)+4096E_6(16\tau)\bigr)$, $F=\Phi_{16}/Dt$ (weight $4$),
$A_n=[t^n]F=1,-52,268,-1920,9536,\dots$ (integral to $n=900$).

* $\Phi_{16}|_6W_{16}=+\Phi_{16}$ and $F|_4W_{16}=+F$
  **[numerical, 230 digits]** — and **[proved]** from
  $E_k(d\tau)|_kW_N=(N/d^2)^{k/2}E_k((N/d)\tau)$, which maps the coefficient
  vector $(1,-85,1428,-5440,4096)$ to itself.
* $\varepsilon_F=+1$, $i^4\varepsilon=+1$: $F(\tau_c)=-5.8337830202\ldots\neq0$
  and $DF(\tau_c)=\frac{4\cdot4}{4\pi}F(\tau_c)$ to $10^{-230}$.
  **Not** the CM mechanism — yet $K$ is still not elementary, because $k=4$.
* $R:=x^2F/(Dx)^2$ is a rational function of $x$ of degree $12$
  **[exact fit, $168$ over-determined equations]**, with
  $$R(x_c)=\frac{7083\sqrt2-10021}{2}=-2.0626688556838846688\ldots \quad(x_c=\tfrac1{2\sqrt2}),$$
  and $Dx(\tau_c)=\dfrac{2+\sqrt2}{32}\cdot\dfrac{\Gamma(1/4)^4}{\pi^3}$
  **[numerical, 213 digits]**.
* **Closed form:**
  $$K=\frac{(10021-7083\sqrt2)(2+\sqrt2)\sqrt{4+\sqrt2}}{16}\cdot
  \frac{\Gamma(1/4)^4}{\pi^{9/2}}=2.04997125634010629357510092802\ldots,$$
  $K\pi^{9/2}/\Gamma(1/4)^4=2.04832043898164671480429450103\ldots$
  Agrees with the fold formula to $\ge150$ digits; Richardson ($A_n$ by series
  reversion to $n=900$, $18$ nodes down to $n=305$) agrees to
  **$1.0\cdot10^{-41}$**.  $A_n<0$ for large $n$.

### 3.3 Level-$12$ $\zeta(7)$ parent — weight $6$, Fricke-odd, $F(\tau_c)\neq0$

$x=\eta_4^2\eta_{12}^2/(\eta_1^2\eta_3^2)$, $t=x/(16x^2+2x+1)$,
$(B,C)=(2,16)$, $\lambda_1=10$, $\lambda_2=-6$, $t_+=1/10$, $\tau_c=i/\sqrt{12}$;
$\Phi_{12}^-=\frac1{480}\sum_{d\mid12}c_dE_8(d\tau)$,
$c=(1,-572,11583,-36608,46332,-20736)$, $F=\Phi^-_{12}/Dt$ (weight $6$),
$A_n=1,-443,13816,-120335,1042756,\dots$ (integral to $n=900$).

* **New fact:** this host is *not* a Hauptmodul — Ligozat gives $x$ two simple
  zeros ($\infty$ and the cusp $1/4$) and two simple poles (cusps $0$ and
  $1/3$), i.e. $\deg x=2$ on $X_0(12)$, which is why $F/(Dx)^3$ admits no
  rational expression in $x$ (searched to degree $26$).  In the *genuine*
  Hauptmodul $v=1/h_{12}$ of §3.1 one has
  $$x=\frac{v(1+3v)}{1+4v},\qquad
  t=\frac{v(12v^2+7v+1)}{(12v^2+(5+\sqrt3)v+1)(12v^2+(5-\sqrt3)v+1)}\quad[\text{exact},\ O(q^{953})],$$
  and $v^3F/(Dv)^3$ *is* rational in $v$, of degree $20$ **[exact fit]**.
  This corrects the tacit assumption that the level-12 $\zeta(7)$ host is a
  Hauptmodul, and explains the repo's "no recurrence of order $\le20$".
* $\Phi^-_{12}|_8W_{12}=-\Phi^-_{12}$, $F|_6W_{12}=-F$
  **[numerical, 230 digits; proved from the $E_8$ oldform rule]**.
  $i^6\varepsilon=(-1)(-1)=+1$, so $F(\tau_c)=66.7635091135\ldots\neq0$ and
  $DF(\tau_c)=\frac{6\sqrt{12}}{4\pi}F(\tau_c)$ to $10^{-228}$.
* Exact CM ratios at $\tau_c$ (same point as §3.1, so the same $Dv$):
  $$\frac{F(\tau_c)}{(Dv)^3}=2025\bigl(490730\sqrt3-849969\bigr),\qquad
  \frac{D^2t(\tau_c)}{(Dv)^2}=\frac{1152\sqrt3-2016}{25},\qquad
  \frac{Dx(\tau_c)}{\Omega^2}=\frac{3^{3/2}}{2^{17/3}},$$
  the first by `lindep` (stable at $80$, $120$, $200$ digits), the second and
  third by `algdep` degree $2$ / $6$ and then verified to $198$/$213$ digits.
* **Closed form:**
  $$K=\frac{3\sqrt3}{2}\cdot2025\bigl(490730\sqrt3-849969\bigr)
  \Bigl(\frac{3(45+26\sqrt3)}{2^{17}}\Bigr)^{2/3}
  \sqrt{\frac{5}{2016-1152\sqrt3}}\;\frac{\Gamma(1/3)^{12}}{\pi^{19/2}}
  =31.7215208279640907727003835703\ldots$$
  $K\pi^{19/2}/\Gamma(1/3)^{12}=12.2666168265985972583025946036\ldots$
  Agrees with the fold formula to $\ge148$ digits; Richardson ($n\le900$, $18$
  nodes) to **$9.7\cdot10^{-42}$**.

---

## 4. Master table

| row | $N$ | $k$ | $\varepsilon$ | $i^k\varepsilon$ | $F(\tau_c)=0$? | $K$ (30 digits) | closed form | verified |
|---|---|---|---|---|---|---|---|---|
| $\alpha$ Domb | 12 | 2 | $-1$ | $+1$ | no | $0.359174244250333123378163967255$ | $2\pi^{-3/2}$ | 60 d |
| $\gamma$ Apéry | 6 | 2 | $-1$ | $+1$ | no | $0.220043767112643037850689759810$ | $(1{+}\sqrt2)^2 2^{-9/4}\pi^{-3/2}$ | 73 d |
| $\varepsilon$ | 8 | 2 | $-1$ | $+1$ | no | $0.257797308917838939757311495404$ | $\tfrac12\sqrt{4{+}3\sqrt2}\,\pi^{-3/2}$ | 68 d |
| $\zeta$ | 9 | 2 | $-1$ | $+1$ | no | $0.260201758994440992303217903862$ | $\tfrac38(\sqrt2{+}\sqrt6)\pi^{-3/2}$ | 69 d |
| Cooper $s_7$ | 7 | 2 | $-1$ | $+1$ | no | $0.233290514929399000627453098124$ | $\tfrac34\sqrt3\,\pi^{-3/2}$ | 71 d |
| Cooper $s_{10}$ | 10 | 2 | $-1$ | $+1$ | no | $0.253974543736963879143053219739$ | $\sqrt2\,\pi^{-3/2}$ | 67 d |
| Cooper $s_{18}$ | 18 | 2 | $-1$ | $+1$ | no | $0.761923631210891637429159659216$ | $3\sqrt2\,\pi^{-3/2}$ | 46 d |
| $\zeta(5)$ level 12 | 12 | 4 | $-1$ | $-1$ | **yes** | $0.685371848490534405951430044881$ | $\overline{\mathbf Q}\,\Gamma(1/3)^{12}\pi^{-17/2}$ | 62 d |
| $\zeta(5)$ level 16 | 16 | 4 | $+1$ | $+1$ | no | $2.04997125634010629357510092802$ | $\overline{\mathbf Q}\,\Gamma(1/4)^{4}\pi^{-9/2}$ | 41 d |
| $\zeta(7)$ level 12 | 12 | 6 | $-1$ | $+1$ | no | $31.7215208279640907727003835703$ | $\overline{\mathbf Q}\,\Gamma(1/3)^{12}\pi^{-19/2}$ | 42 d |

("verified" = digits of agreement between the Richardson extrapolation of the
sequence and the closed form / fold formula.  The fold formula and the closed
form agree to $\ge148$ digits throughout; the Richardson figure is the binding
one.)

---

## 5. What this changes in `ASYMPTOTIC_CONSTANTS.md`

1. **Theorem 3.1(a) needs its $k=2$ clause promoted from a "moreover" to the
   defining hypothesis.**  The census contains two rows with $i^k\varepsilon=+1$
   whose $K$ is *not* elementary (level-16 $\zeta(5)$, level-12 $\zeta(7)$): the
   surviving transcendental is $\Omega^{k-2}$, not $\Omega^{k}$.  The clean
   statement is $K\in\overline{\mathbf Q}\cdot\Omega^{k-2}\pi^{-3/2}$ in case (a)
   and $K\in\overline{\mathbf Q}\cdot\Omega^{k}\pi^{-1/2}$ in case (b), the two
   agreeing at $k=2$ up to which power of $\pi$ is displayed.  The title claim
   is exactly right *within* the weight-two world (the four real-Fricke-fold
   third-order sporadic rows $\alpha,\gamma,\varepsilon,\zeta$ plus all three
   Cooper rows).
2. **Remark 3.2 is confirmed and can be sharpened.**  Not only is Apéry's
   $(1+\sqrt2)^2/2^{9/4}$ equal to $\frac{\sqrt6}{4}|F/Dh|(\tau_c)\sqrt{2x_+/(-r'')}$
   — the whole weight-two family collapses to
   $K=\frac{\sqrt N}{2\pi^{3/2}}\sqrt{\lambda_1/(\lambda_1-\lambda_2)}$, because
   $(F/Du)(\tau_c)=\sqrt C$ and $2x_+/(-r''(u_c))=\lambda_1 C^{-3/2}$ *always*.
   $K$ therefore depends on the row only through $(N,\lambda_1,\lambda_2)$ —
   no modular data at all beyond the level.
3. **Remark 3.3's claim "a fold at the Fricke point with $F(\tau_c)\neq0$ forces
   $\varepsilon_F=-1$" is confirmed and now has a one-line proof**: $F=D\log u$
   and $u|W_N=1/(Cu)$.  The identity $F=D\log u$ (equivalently
   $\sum a_nx^n=D\log u$ with $x=u/(1+Bu+Cu^2)$) is itself worth recording; it
   is what makes the weight-two case elementary.
4. **Cooper's three rows belong in the table.**  Their dominant singularity is a
   Fricke fold ($q_c=e^{-2\pi/\sqrt N}$ to $\ge200$ digits), they satisfy the same
   structure theorem, and they give the three prettiest constants in the census.
5. **Correction to the level-12 $\zeta(7)$ description.**  Its host
   $x=\eta_4^2\eta_{12}^2/(\eta_1^2\eta_3^2)$ has degree $2$ on $X_0(12)$; the
   Hauptmodul is $v=1/h_{12}$, the *same* one that carries the level-12
   $\zeta(5)$ system, with $x=v(1+3v)/(1+4v)$.  Both $\zeta(5)$ and $\zeta(7)$
   at level $12$ therefore live on one coordinate and share one CM period
   $Dv(i/\sqrt{12})$; they differ only in the companion's Fricke sign, hence in
   $\pi^{-17/2}$ vs $\pi^{-19/2}$.

---

## 6. Reproduction

All scripts in this directory; `gp -q <file>`.

| script | what |
|---|---|
| `lib.gp` | mirror map / exact sequence from a MUM order-3 recurrence |
| `main1.gp` | fold location, $q_c=e^{-2\pi/\sqrt N}$, first fold-formula values, 7 rows |
| `main3.gp` | the three exact series identities $x=r(u)$, $F=Du/u$, $P_3=((1-Cu^2)/(1+Bu+Cu^2))^2$ |
| `etaid.gp` | eta-quotient identification of $u$ (exact linear algebra) |
| `alconst.gp` | $\prod_d(N/d)^{r_d}=1/C^2$ exactly |
| `fricke.gp` | $u|W_N=1/(Cu)$, $F|_2W_N=-F$, $x\circ W_N=x$ numerically |
| `main4.gp` | $a_n$ from $(x,F)$ vs recurrence to $n=40$; recurrence to $n=6000$; Richardson; fold formula; closed form; all seven rows |
| `cm1.gp` | level-12 $\zeta(5)$ control from scratch |
| `cm2.gp`,`cm4.gp` | level-16 $\zeta(5)$ |
| `cm3.gp`,`cm6.gp`,`cm7.gp` | level-12 $\zeta(7)$ |
| `cm5.gp` | exact evaluation of the companion rational functions at the fold |
| `cm8.gp` | level-12 $\zeta(5)$: order-14 recurrence to $n=4200$ + Richardson |
| `fits.gp`,`fits2.gp` | rational-function fits of $F/(Dh)^{k/2}$ |
| `final.gp`,`table.gp` | closed-form assembly and the master table |

Runtimes: everything under 10 minutes per script on this machine; `main4.gp`
(seven rows, $n\le6000$ exact, 400-digit Richardson) is the longest.
