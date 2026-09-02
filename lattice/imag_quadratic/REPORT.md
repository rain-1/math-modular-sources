# Imaginary quadratic fields in the arithmetic-holonomy architecture:
# the number-field tax is zero, and there is nothing to spend it on

*Working directory `lattice/imag_quadratic/`. All scripts and raw outputs listed in §10.
Conventions: `CDT_FINDER.md` §§1–6 ($\tau$, ceil, entryC, entryR, margin),
`NUMBER_FIELD_HOLONOMY.md` (the number-field budget is the **average** over archimedean places),
`CUSP_PERIODS.md` §3 (fold-regularity = vanishing of the constant term at $\infty$ and at the
near cusp), `lattice/hostscan/REPORT.md` (the census of genus-zero four-point hosts over $\mathbf Q$).
Tags: **[exact]** = closed form / exact arithmetic; **[verified]** = high-precision or exact
machine computation in this task; **[cited]**. No irrationality claim is made anywhere.*

---

## 0. Verdict

| claim | verdict |
|---|---|
| **Task 0. For $K$ imaginary quadratic the arithmetic holonomy bound is *literally* the rational bound.** One archimedean place, complex, $d_v=2=n_K$, so $m\le \frac{d_vM_v}{d_v\log\lvert\varphi_v'(0)\rvert-n_K\tau}=\frac{M}{\log\lvert\varphi'(0)\rvert-\tau}$. There is no second place, no hypothesis tax, and fold-regularity at the conjugate embedding is complex conjugation, hence automatic. Derived from `NUMBER_FIELD_HOLONOMY.md` Thm 2.1 in §1 below, with the meaning of $M_v$, $\log\lvert\varphi_v'(0)\rvert$ and the Liouville step at a complex place spelled out. | **[exact]** §1 |
| **But an imaginary quadratic field buys a new Diophantine statement only if the period is *not real*,** because $K\cap\mathbf R=\mathbf Q$: for real $\Pi_1,\dots,\Pi_r$, $K$-linear independence of $1,\Pi_1,\dots,\Pi_r$ is *equivalent* to $\mathbf Q$-linear independence (take real and imaginary parts). | **[exact]** §2.1 |
| **A non-real period forces a nebentypus $\varepsilon$ with $\varepsilon\ne\bar\varepsilon$, i.e. $\operatorname{ord}\varepsilon\ge3$.** And a source with nebentypus $\varepsilon$ lives on $\Gamma_{H'}$ with $H'\subseteq\ker\varepsilon$; for it to be a section over the host curve one needs $H'\cdot\{\pm1\}=\bar H$. On any $\Gamma_0(N)$-type host ($\bar H=(\mathbf Z/N)^\times$) this forces $\operatorname{ord}\varepsilon\le2$. | **[exact]** §3.1 |
| **Complete census of intermediate hosts.** Over *all* subgroups $H\le(\mathbf Z/N)^\times$ and all $N\le60$ there are exactly **25** genus-zero curves $X_H$, of which exactly **six** have four special points: $X_0(5),X_0(6),X_0(7),X_0(8),X_0(9),X_1(5)$ — an independent re-derivation of `hostscan` §2, now over the full intermediate lattice and not only $\Gamma_0$. | **[verified]** §3.2 |
| **Of those six, exactly one admits a nebentypus of order $>2$: $X_1(5)$ with the odd quartic $\psi_4$ (values in $\mathbf Q(i)$).** The other five admit only characters of order $\le2$. Among all 25 genus-zero curves the ones admitting order $\ge3$ are $X_1(5)$ (4 special points), $X_1(7)$ (6), $X_1(9)$ (8), $X_1(10)$ (8), two level-13 curves (12 and 8), a level-16 curve (10) and a level-25 curve (12). | **[verified]** §3.3 |
| **And on $X_1(5)$ the fold-regular subspace misses the $\psi_4$-component entirely.** $\dim M_3^{\rm Eis}(\Gamma_1(5),\psi_4)=2$, spanned by the inner $E_3^{\psi_4,\mathbf 1}$ ($a_0=0$ at $\infty$, $\ne0$ at the near cusp $0$) and the outer $E_3^{\mathbf 1,\psi_4}$ ($a_0\ne0$ at $\infty$, $=0$ at $0$); the two constant-term functionals therefore have **rank 2** on it and $\dim(\text{fold-regular}\cap\psi_4)=\mathbf 0$. (On the full 4-dimensional space the rank is 2 and the fold-regular subspace is 2-dimensional — `hostscan` §10.3 — but it is the conjugation-stable diagonal, whose interesting line $\Phi_{\rm new}$ is defined over the **real** quadratic field $\mathbf Q(\sqrt5)$ and pays the $1.46$-nat second-place tax of `GAMMA15_CLOSURE.md`.) | **[verified, exact]** §4 |
| **The supply is not the obstruction; the geometry is.** Fold-regular Eisenstein directions with coefficients in $\mathbf Q(i)$ or $\mathbf Q(\zeta_3)$ and a genuinely non-real period exist from level **7** upward ($N=7,9,10,13,14,15,16,\dots$: exact PARI ranks for all $N\le30$ and both weights). Every level that carries them has $\ge6$ special points; every four-point host has only real characters. The two lists are disjoint. | **[verified]** §5 |
| **Route (b), hosts over $K$.** The only four-point genus-zero curves with a complex-conjugate pair of special points are $X_0(5)$ (two $\mathbf Z/2$ points at $\mu=-11\pm2i$), $X_0(7)$ (two $\mathbf Z/3$ points at $\mu=-\tfrac{13}2\pm\tfrac{3\sqrt3}2 i$) and $X_0(9)$ (the conjugate cusps $1/3,2/3$ at $\mu=-\tfrac92\mp\tfrac{3\sqrt3}2i$) — all three $\mu$-configurations recomputed here to 80 digits and recognised exactly. | **[verified]** §6.1 |
| **Placing the pole of $x$ at a cusp always gives $\lambda_1=\bar\lambda_2$: dead.** $\lvert\lambda_1\rvert=\lvert\lambda_2\rvert$, no dominant characteristic root, no Apéry limit, no geometric convergence. This kills the three configurations `hostscan` §5 records as "complex fold" ($\Gamma_0(5)$ pole 0, $\Gamma_0(7)$ pole 0, Zagier B). | **[exact]** §6.2 |
| **Placing the pole at one member of the conjugate pair gives a genuinely $K$-rational coordinate — three of them, and all three fail.** $X_0(5)/\mathbf Q(i)$: $\lambda=(11-2i,\,-4i)$, $\lambda_2^{\rm norm}=4$, minimal recurrence **5 terms** (not CDT-shape). $X_0(7)/\mathbf Q(\sqrt{-3})$: $\lambda=(\tfrac{13-3\sqrt{-3}}2,\,-3\sqrt{-3})$, $\lambda_2^{\rm norm}=3\sqrt3$, minimal recurrence **4 terms** (weight one) / **5 terms** (weight two). $X_0(9)/\mathbf Q(\sqrt{-3})$: pole at the cusp $1/3$, $\lambda=(\tfrac{9+3\sqrt{-3}}2,\,3\sqrt{-3})$, **genuinely three-term**, but $\lvert\lambda_1\rvert=\lvert\lambda_2\rvert=3\sqrt3$ — the triangle $\{0,\mu,\bar\mu\}$ is equilateral — so again no Apéry limit. | **[verified]** §6.3–6.4 |
| **The unit condition over an imaginary quadratic field is "$\lambda_2$ is a root of unity".** $\mathcal O_K^\times=\mu_K$, so $\lambda_2^{\rm norm}=1\iff\lambda_2\in\{\pm1,\pm i,\pm\zeta_6^{\,j}\}$, and then $\lvert\lambda_1\rvert=\lvert N(\lambda_1)\rvert^{1/2}\ge\sqrt2$. No four-point configuration realises it. (`hostscan` §7's off-shape $N=18$ configuration, $\lambda_2$ a primitive sixth root of unity, is the near miss: $\mathbf Q(\sqrt{-3})$, unit, but **seven** singular points.) | **[exact]** §2.2, §6.5 |
| **Atkin–Lehner and Fricke quotients add nothing.** They are quotients of $\Gamma_0(N)$-type curves, so $\operatorname{ord}\varepsilon\le2$ by the rule above; the only intermediate curve admitting order $\ge3$ that survives an Atkin–Lehner involution at all is $X_1(10)+W_2$ ($W_5$ and $W_{10}$ send $\psi_4\mapsto\bar\psi_4$), and it has $\ge5$ special points. Weight-four sources with even characters of order $3,4,6$ need trivial nebentypus, hence $E_4^{\psi,\bar\psi}$ with $\operatorname{cond}(\psi)^2\mid N$, i.e. $N\ge25$ (order 4) or $N\ge49$ (orders 3, 6). | **[exact]** §7 |
| **Best margin over everything imaginary quadratic:** the only live scores are $-18.0$ ($X_0(5)$, and only if a $k=2$ row existed there, which it does not: the honest entry is $-42.4$) and $-21.4$ ($X_0(7)$). The prize that does not exist — an Apéry-perfect imaginary-quadratic host with a non-real period — would have scored **entryC $+1.310$, entryR $+0.846$, margin $+0.005$: CDT's own numbers, with no tax at all.** | **[verified]** §8 |

**One sentence.** *Over an imaginary quadratic field the arithmetic holonomy bound costs
nothing extra — one complex place, $d_v=2=n_K$, the two factors cancel and the inequality is
letter-for-letter the rational one — but the only thing such a field could buy is a **non-real**
period, that needs a nebentypus of order at least three, a nebentypus of order at least three
needs a group $\Gamma_{H'}$ with $H'\subseteq\ker\varepsilon$ and therefore a curve with at
least six special points at every level $N\le60$ except $X_1(5)$, and on $X_1(5)$ the two
constant-term conditions have rank two on the two-dimensional $\psi_4$-component and annihilate
it; the three genuinely $\mathbf Q(i)$- and $\mathbf Q(\sqrt{-3})$-rational four-point
coordinates that do exist (poles at the conjugate special points of $X_0(5)$, $X_0(7)$, $X_0(9)$)
are respectively not three-term, not three-term, and three-term with two singularities at the
same distance.*

---

## 1. Task 0 — the bound over an imaginary quadratic field  **[exact]**

Let $K$ be imaginary quadratic, $n_K=[K:\mathbf Q]=2$. $K$ has **one** archimedean place $v$,
and it is complex: $K_v=\mathbf C$, $d_v=[K_v:\mathbf R]=2$, $\sum_vd_v=2=n_K$. The two
embeddings $\sigma,\bar\sigma\colon K\hookrightarrow\mathbf C$ induce the *same* absolute value
and are the two embeddings counted by $d_v=2$.

Run `NUMBER_FIELD_HOLONOMY.md` Thm 2.1 verbatim.

* **Data.** $f_1,\dots,f_m\in K[[x]]$, $K(x)$-linearly independent, with $L(n)a_{i,n}\in\mathcal O_K$
  for a *rational* lcm-type denominator $L(n)$, $\log L(n)\sim\sigma n$ (this is the modular
  situation: $d_n^{\,k}B_n$ is an algebraic integer by Beukers' argument, and $d_n\in\mathbf Z$).
* **What $M_v$ and $\log\lvert\varphi_v'(0)\rvert$ mean at a complex place.** One holomorphic
  $\varphi=\varphi_v\colon(\mathbf D,0)\to(\mathbf C,0)$ with $f_i^{(v)}\circ\varphi$ meromorphic
  on $\mathbf D$, where $f^{(v)}$ means $\sigma$ applied to the coefficients; $M_v=\max_{|z|=1}\log^+|\varphi(z)|$
  and $\log|\varphi'(0)|$ use the **ordinary** complex absolute value, the local degree $d_v=2$
  appearing as an explicit factor, not inside the absolute value.
  *Subtlety, and it is the whole point:* the conjugate embedding $\bar\sigma$ gives the functions
  $\bar f_i(\bar x)$, i.e. the complex conjugates. It is the *same place*: one imposes holomorphy
  of $f_i^{(v)}\circ\varphi$ once, and the conjugate condition
  $z\mapsto\overline{f_i^{(v)}(\varphi(\bar z))}$ is then automatic, because the conjugate of a
  function holomorphic on $\mathbf D$ is holomorphic on $\overline{\mathbf D}=\mathbf D$. The
  disc, the contour and the Bost–Charles integral are **shared**; there is nothing to choose
  twice and nothing to pay twice. In the real quadratic case one has genuinely two discs and
  the deeper of the two conjugate geometries is a tax (`GAMMA15_CLOSURE.md` §3.1: $-1.46$ nats
  on $\Gamma_1(5)$); here that mechanism does not exist.
* **Siegel.** $mD$ unknowns in $\mathcal O_K$ = $2mD$ rational unknowns; $T$ equations over $K$ =
  $2T$ rational equations; the ratio, hence the admissible $T=(1-\varepsilon)mD$, is unchanged.
* **Liouville.** $\beta=[x^{n_0}]F\in K^\times$ with $L(n_0)\beta$ (times a height factor) a nonzero
  algebraic integer, so $|N_{K/\mathbf Q}(L(n_0)\beta)|\ge1$. At a complex place
  $N_{K/\mathbf Q}(\beta)=\sigma(\beta)\overline{\sigma(\beta)}=|\sigma(\beta)|^2$, i.e.
  $\sum_vd_v\log|\sigma_v(\beta)|=2\log|\sigma(\beta)|\ge-2\log L(n_0)\ge-n_K\sigma n_0-o(n_0)$.
* **Cauchy.** $\log|\sigma(\beta)|+n_0\log|\varphi'(0)|\le DM_v+o(D)$; multiply by $d_v=2$.

Both sides carry the factor $2$:
$$n_0\bigl(2\log|\varphi'(0)|-2\sigma\bigr)\le 2DM_v+o(D)
\qquad\Longrightarrow\qquad
\boxed{\;m\;\le\;\frac{2M_v}{2\log|\varphi'(0)|-2\sigma}\;=\;\frac{M_v}{\log|\varphi'(0)|-\sigma}\;}$$
which is the $\mathbf Q$-statement of `CDT_UNPACKED.md` §1 with no modification. The same happens
in the refined (Bost–Charles / slopes) form, Thm 2.2: the finite places above $p$ contribute
$\sum_{v\mid p}d_v\log|L(n)|_v=n_K\log|L(n)|_p=2\log|L(n)|_p$, again matching the factor $2$ in
the archimedean terms. **Over an imaginary quadratic field the denominator rate $\tau$, the
ceiling $\log(256/\lambda_2^{\rm norm})$ and the numerator are exactly CDT's**, with
$\lambda_2^{\rm norm}=|N(\lambda_2)|^{1/2}=|\lambda_2|$, the modulus at the unique place.
`10_score.out` reproduces the three entries side by side:

```
Q, CDT's own host                                     L-tau = +0.84645   m <= 13.994   margin(14) = +0.0053
Q(i) host, same geometry (one complex place, d_v=2)   L-tau = +0.84645   m <= 13.994   margin(14) = +0.0053
Q(sqrt5) Gamma_1(5) (two real places, averaged)       L-tau = +0.66956   m <= 16.154   margin(14) = -1.4425
```

---

## 2. Two elementary facts that decide the whole question

### 2.1 $K\cap\mathbf R=\mathbf Q$: a real period gains nothing  **[exact]**

Let $\Pi_1,\dots,\Pi_r\in\mathbf R$ and $K=\mathbf Q(\sqrt{-d})$. A relation
$a_0+\sum a_j\Pi_j=0$ with $a_j=\alpha_j+\beta_j\sqrt{-d}\in K$ not all zero splits into the two
rational relations $\alpha_0+\sum\alpha_j\Pi_j=0$ and $\beta_0+\sum\beta_j\Pi_j=0$, at least one
of them nontrivial; conversely a rational relation is a $K$-relation. Hence

> **$K$-linear independence of $1,\Pi_1,\dots,\Pi_r$ over an imaginary quadratic $K$ is
> *equivalent* to $\mathbf Q$-linear independence, whenever the $\Pi_j$ are real.**

Contrast the real quadratic case, where $K$-independence is genuinely stronger and is what
Beukers' $8\zeta(3)-5\sqrt5L(3,\chi_5)$ and the $\Gamma_1(5)$ target (T) are about. So an
imaginary quadratic system is Diophantically new **iff its period is not real**, and then the
statement it would prove is of the shape $L(2,\psi_4)\notin\mathbf Q(i)$ or
$L(3,\psi_3)\notin\mathbf Q(\zeta_3)$ for a non-real character $\psi$ — statements that are
open and are not implied by any real-period result.

### 2.2 The unit condition becomes "root of unity"  **[exact]**

$\lambda_2$ is an algebraic integer (integrality of the pure module $\mathrm{Li}_j(x/s)$,
`CDT_FINDER.md` §3), and `hostscan` §1: Apéry-perfect $\iff|N(\lambda_2)|=1\iff\lambda_2\in\mathcal O_K^\times$.
For $K$ imaginary quadratic $\mathcal O_K^\times=\mu_K$ is finite, so

> Apéry-perfect over an imaginary quadratic $K$ $\iff$ $\lambda_2$ is a **root of unity**
> ($\pm1$; $\pm i$ in $\mathbf Q(i)$; $\pm\zeta_6^{\,j}$ in $\mathbf Q(\sqrt{-3})$),

and then $|\lambda_2|=1$ exactly, ceiling $=\log256=5.545$, and $|\lambda_1|=|N(\lambda_1)|^{1/2}\ge\sqrt2$
so the convergence rate is at worst $2^{-1/2}$. This is far more rigid than the real quadratic
condition ($\lambda_1\lambda_2=\pm1$ with a fundamental unit of infinite order, which is how
Apéry's $17\pm12\sqrt2$ and Zagier D's $\varphi^{\pm5}$ arise).

---

## 3. Where a non-real nebentypus can live

### 3.1 The rule  **[exact]**

Let the host curve be $X_{\bar H}$, $\bar H=H\cdot\{\pm1\}$. A source $\Phi$ (and the weight-$w$
form $F$ carrying the row) with nebentypus $\varepsilon$ is modular on $\Gamma_{H'}$ with
$H'\subseteq\ker\varepsilon$; for $F$ to be a section of a line bundle on the host curve — i.e.
for $A_n=[x^n]F$ to be defined at all — one needs $\Gamma_{H'}$ to have the *same image in*
$\mathrm{PSL}_2(\mathbf Z)$ as the host group, i.e.
$$H'\cdot\{\pm1\}=\bar H .$$
Consequently the available nebentypus characters are exactly the characters of
$(\mathbf Z/N)^\times/H'$ for such $H'$. For a $\Gamma_0(N)$-type host ($\bar H=(\mathbf Z/N)^\times$)
this forces $(\mathbf Z/N)^\times/\ker\varepsilon$ to be killed by $\{\pm1\}$, i.e.
$\operatorname{ord}\varepsilon\le2$: **only real characters, only real Eisenstein coefficients,
only real periods.** Odd weight is possible (that is how CDT's $\chi_{-3}$ and Zagier E's
$\chi_{-4}$ occur) but the character is quadratic.

### 3.2 The complete four-point census over the intermediate lattice  **[verified]** (`01_groups.py`)

For every $N\le60$ and every subgroup $H\le(\mathbf Z/N)^\times$ with $-1\in H$, the cosets
$\Gamma_H\backslash\mathrm{SL}_2(\mathbf Z)$ are the primitive pairs $(c,d)\in(\mathbf Z/N)^2$
modulo $(c,d)\sim(hc,hd)$; right multiplication by $S$, $T$, $ST$ gives
$\nu_2,\nu_\infty,\nu_3$ and $g=1+\mu/12-\nu_2/4-\nu_3/3-\nu_\infty/2$. Result: **25 genus-zero
curves**, of which **six** have $\nu_2+\nu_3+\nu_\infty=4$:

| curve | $H$ | $\mu$ | $\nu_2$ | $\nu_3$ | cusps | #special |
|---|---|---|---|---|---|---|
| $X_1(5)$ | $\{1,4\}$ | 12 | 0 | 0 | 4 | **4** |
| $X_0(5)$ | all | 6 | 2 | 0 | 2 | **4** |
| $X_0(6)$ | all | 12 | 0 | 0 | 4 | **4** |
| $X_0(7)$ | all | 8 | 0 | 2 | 2 | **4** |
| $X_0(8)$ | all | 12 | 0 | 0 | 4 | **4** |
| $X_0(9)$ | all | 12 | 0 | 0 | 4 | **4** |

This reproduces `hostscan` §2 and strengthens it: the sweep there covered $\Gamma_0(N)$ and
$\Gamma_1(5)$; here the *entire* lattice of intermediate groups $\Gamma_H$ at $N\le60$ is
enumerated and nothing else is four-point. (Full list of all 25 in `01_groups.out`; e.g.
$X_1(7)$ has 6 cusps, $X_1(9)$ and $X_1(10)$ have 8, the two level-13 curves have 12 and 8.)

### 3.3 Character availability  **[verified]** (`03_chars.py`)

For every one of the 25 curves, all $H'$ with $H'\{\pm1\}=\bar H$ and all characters of
$(\mathbf Z/N)^\times/H'$:

| curve | #special | character orders available | non-real possible? |
|---|---|---|---|
| **$X_1(5)$** | **4** | $1,2,\mathbf 4$ | **yes: $\psi_4$, $\mathbf Q(i)$** |
| $X_0(5)$ | 4 | $1$ | no |
| $X_0(6)$ | 4 | $1,2$ ($\chi_{-3}$) | no |
| $X_0(7)$ | 4 | $1,2$ ($\chi_{-7}$) | no |
| $X_0(8)$ | 4 | $1,2$ ($\chi_{-4}$ or $\chi_{-8}$) | no |
| $X_0(9)$ | 4 | $1,2$ ($\chi_{-3}$ mod 9) | no |
| $X_1(7)$ | 6 | $1,2,3,6$ | yes, $\mathbf Q(\zeta_3)$ |
| $X_1(9)$ | 8 | $1,2,3,6$ | yes, $\mathbf Q(\zeta_3)$ |
| $X_1(10)$ | 8 | $1,2,4$ | yes, $\mathbf Q(i)$ |
| $X_{\{1,5,8,12\}}(13)$ | 12 | $1,3$ | yes, $\mathbf Q(\zeta_3)$ |
| $X_{\{1,3,4,9,10,12\}}(13)$ | 8 | $1,2,4$ | yes, $\mathbf Q(i)$ |
| $X_{\{1,7,9,15\}}(16)$ | 10 | $1,2,4$ | yes, $\mathbf Q(i)$ |
| $X_{\{1,4,6,\dots\}}(25)$ | 12 | $1,2,4$ | yes, $\mathbf Q(i)$ |
| (the remaining 12 curves) | 3–10 | $\le2$ | no |

**The two columns "four special points" and "non-real character available" intersect in exactly
one row, $X_1(5)$.**

---

## 4. $X_1(5)$: the one candidate, and why it is empty  **[verified, exact]** (`04_g15.gp`)

Mod 5 the characters are $\mathbf 1$, $\chi_5$ (even, quadratic) and the two odd quartic
$\psi_4,\bar\psi_4$ with $\psi_4(2)=i$; weight three needs odd nebentypus, so a $\mathbf Q(i)$-rational
weight-three source must lie in the $\psi_4$-component (a line inside a single non-real nebentypus
component is the only way a *fold-regular direction* can be non-real without being a $K$-combination
of $\mathbf Q$-rational ones, since the fold-regular subspace itself is defined over $\mathbf Q$).

PARI, exact arithmetic in $\mathbf Q(\zeta_{24})$ for the cusp-0 expansion:

```
dim M_3^Eis(Gamma_1(5)) = 2 + 2 = 4
rank of (a_0 at oo, a_0 at cusp 0) on the full 4-dimensional space = 2
dim(fold-regular) = 2                        [= hostscan sec.10.3]

dim M_3^Eis(Gamma_1(5), psi_4) = 2
basis 1 : 0 + q + (4+i)q^2 + (9-i)q^3 + ...      a_0(oo) = 0 ,  a_0(cusp 0) != 0
basis 2 : (-2-i)/5 + q + (1+4i)q^2 + ...          a_0(oo) != 0 , a_0(cusp 0) = 0
rank = 2 ,  dim(fold-regular cap psi_4-component) = 0
```

Structurally: basis 1 is the **inner** direction $E_3^{\psi_4,\mathbf 1}$ and basis 2 the **outer**
$E_3^{\mathbf 1,\psi_4}$; $a_0(\infty)$ vanishes exactly on inner directions and $a_0(\text{cusp }0)$
exactly on outer ones (`CUSP_PERIODS.md` Thm 3.2 and its Fricke consequence). With $N=\operatorname{cond}\psi_4=5$
there is **no oldform shift** ($\tau(N/f)=1$) to build $(1-d^{k}V_d)E$ from, and **no mixed
direction** $E^{\psi_1,\psi_2}$ with both $\psi_i\ne\mathbf 1$ (that needs $f_1f_2\mid N$ with
$f_1,f_2>1$, impossible at $N=5$). So the component is exactly one inner plus one outer direction
and the two conditions annihilate it. This is precisely why the interesting fold-regular line on
$\Gamma_1(5)$ is the *diagonal* $\Phi_{\rm new}=(1+i\varphi^5)E_3^{\psi_4,\mathbf1}+(1-i\varphi^5)E_3^{\bar\psi_4,\mathbf1}$
of `hostscan` §10, whose coefficients are real (in $\mathbf Z[\varphi]$) and whose field of
definition is the **real** quadratic $\mathbf Q(\sqrt5)$ — hence two real places, hence the tax
that `GAMMA15_CLOSURE.md` measures as $-1.44$ nats.

At weight four ($r=3$) on $\Gamma_1(5)$ the nebentypus must be even, i.e. $\mathbf 1$ or $\chi_5$:
everything is $\mathbf Q$-rational (Beukers), no imaginary quadratic direction exists
(`hostscan` §10.6, recomputed here as part of `03_chars.out`).

---

## 5. The supply: where imaginary-quadratic fold-regular sources do live  **[verified]** (`09_census.gp`)

Exact PARI ranks of the pair (constant term at $\infty$, constant term at the cusp $0$) on the
$\varepsilon$-component of $M_k^{\rm Eis}(\Gamma_1(N))$, for every $N\le30$ and every $\varepsilon$
of order $3$, $4$ or $6$ (the orders whose values generate an imaginary quadratic field:
$\mathbf Q(\zeta_3)$ for $3$ and $6$, $\mathbf Q(i)$ for $4$); $k=3$ for $\varepsilon$ odd
(period $L(2,\psi)$), $k=4$ for $\varepsilon$ even (period $L(3,\psi)$). Excerpt:

| $N$ | $k$ | $\operatorname{ord}\varepsilon$ | $\dim(\varepsilon\text{-comp})$ | rank | $\dim$ fold-regular | $K$ | host geometry |
|---|---|---|---|---|---|---|---|
| **5** | 3 | 4 | 2 | 2 | **0** | $\mathbf Q(i)$ | **four-point ($X_1(5)$)** |
| 7 | 4 | 3 | 3 | 2 | 1 | $\mathbf Q(\zeta_3)$ | $X_1(7)$: 6 special points |
| 7 | 3 | 6 | 2 | 2 | 0 | $\mathbf Q(\zeta_3)$ | $X_1(7)$: 6 |
| 9 | 3 | 6 | 3 | 2 | 1 | $\mathbf Q(\zeta_3)$ | $X_1(9)$: 8 |
| 9 | 4 | 3 | 4 | 2 | 2 | $\mathbf Q(\zeta_3)$ | $X_1(9)$: 8 |
| 10 | 3 | 4 | 5 | 2 | 3 | $\mathbf Q(i)$ | $X_1(10)$: 8 |
| 13 | 3 | 4 | 4 | 2 | 2 | $\mathbf Q(i)$ | 8 |
| 14 | 3 | 6 | 6 | 2 | 4 | $\mathbf Q(\zeta_3)$ | $X_1(14)$: genus 1 |
| 15 | 3 | 4 | 6 | 2 | 4 | $\mathbf Q(i)$ | genus 1 |
| 16 | 3 | 4 | 5 | 2 | 3 | $\mathbf Q(i)$ | 10 |
| 18 | 4 | 3 | 11 | 2 | 9 | $\mathbf Q(\zeta_3)$ | genus 0, $>4$ |
| $\dots$ | | | | | | | (full table `09_census.out`, 78 rows) |

The mechanisms are the two CDT use at level 6. (i) An **oldform shift**: the inner direction
$E_k^{\psi,\mathbf 1}$ has $a_0=0$ at $\infty$ automatically, and its cusp-$0$ constant term is a
nonzero multiple of $\sum_dc_d\,d^{-k}$ (`CUSP_PERIODS.md` Thm 3.2(i)), so
$(1-d^{\,k}V_d)E_k^{\psi,\mathbf 1}$ is fold-regular as soon as $\operatorname{cond}\psi\cdot d\mid N$
with $d>1$. (ii) A **mixed direction** $E_k^{\psi_1,\psi_2}$ with both $\psi_i\ne\mathbf 1$, which has
$a_0=0$ at both cusps automatically and needs $f_1f_2\mid N$.

The smallest instance of (i) with a non-real character is exhibited explicitly in
`11_level10.out` **[verified]**: at $N=10$, $k=3$,
$$\Phi\;=\;E_3^{\psi_4,\mathbf 1}(\tau)-8\,E_3^{\psi_4,\mathbf 1}(2\tau)
\;=\;q+(-4+i)q^2+(9-i)q^3+(-17-4i)q^4+25q^5+\cdots\in\mathbf Z[i][[q]],$$
with $a_0=0$ at $\infty$ and at the cusp $0$ (both checked exactly in $\mathbf Q(\zeta_{24})$),
lying in the $\psi_4$-component of $M_3^{\rm Eis}(\Gamma_1(10))$, and with near-cusp period
$-\tfrac12L(2,\psi_4)$ where
$$L(2,\psi_4)=0.9587161227168831553919364293\ldots+0.1455658767850895904617045118\ldots\,i$$
is **not real**. A CDT-shape host carrying this source would prove $L(2,\psi_4)\notin\mathbf Q(i)$
at no number-field cost. At $N=5$ the same construction is unavailable because
$\operatorname{cond}\psi_4=5=N$ leaves no shift $d>1$ — that is the entire content of §4.

So the sources are plentiful, they carry genuinely non-real $L$-values, and their coefficient
fields are exactly $\mathbf Q(i)$ and $\mathbf Q(\zeta_3)$. **The obstruction is entirely
geometric: every level carrying one has at least six special points; the four-point locus stops
at $N=9$.** This is the same "supply versus geometry" verdict as `hostscan` §10.5–10.6, one
character-order further out.

---

## 6. Route (b): hosts defined over an imaginary quadratic field

### 6.1 Exact $\mu$-configurations  **[verified, 80 digits]** (`02_geom.gp`)

With a reference Hauptmodul $h=q+O(q^2)$ and $\mu:=1/h$, evaluated at $\gamma(60i)$ for the cusps
and directly at the elliptic fixed points, all values recognised exactly by `bestappr`:

| curve | reference $h$ | $\mu$ of the three non-$\infty$ special points |
|---|---|---|
| $X_0(5)$ | $(\eta_5/\eta_1)^6$ | $0$ (cusp 0), $\;-11+2i$, $\;-11-2i$ (the two $\mathbf Z/2$ points $\tau=(\mp2+i)/5$) |
| $X_0(6)$ | $\eta_2\eta_6^5/(\eta_1^5\eta_3)$ | $0,\;-9,\;-8$ |
| $X_0(7)$ | $(\eta_7/\eta_1)^4$ | $0$, $\;-\tfrac{13}2\pm\tfrac{3\sqrt3}2 i$ (the two $\mathbf Z/3$ points $\tau=(-5+i\sqrt3)/14,(-9+i\sqrt3)/14$) |
| $X_0(8)$ | $\eta_2^2\eta_8^4/(\eta_1^4\eta_4^2)$ | $0,\;-8,\;-4$ |
| $X_0(9)$ | $(\eta_9/\eta_1)^3$ | $0$, $\;-\tfrac92\mp\tfrac{3\sqrt3}2i$ (the **conjugate cusps** $1/3$, $2/3$) |
| $X_1(5)$ | Zagier D's $x$ | $0,\;\varphi^5,\;-\varphi^{-5}$ (exact from the recurrence $t^2-11t-1$) |

$X_0(6)$, $X_0(8)$ and $X_1(5)$ have all special points rational (resp. real quadratic): no
imaginary quadratic coordinate exists on them. The three curves with a conjugate pair are
$X_0(5)$, $X_0(7)$, $X_0(9)$.

### 6.2 Pole at the cusp $0$: $\lambda_1=\bar\lambda_2$, dead  **[exact]**

`hostscan` §3: the admissible normalisations are $\lambda_i=\mu_i-\mu_j$, one per special point
$P_j$ carrying the pole of $x$. Taking $P_j$ = the cusp $0$ ($\mu_j=0$) leaves
$\lambda=(\mu,\bar\mu)$: a complex-conjugate pair, so
$$|\lambda_1|=|\lambda_2|.$$
Then the three-term recurrence has two characteristic roots of equal modulus, $A_n$ oscillates
($A_n\sim2\mathrm{Re}(c\lambda^n)$ up to polynomial factors), no dominant root exists, the
companion ratio $B_n/A_n$ has no limit and **no linear form tends to zero**. This is the honest
reason the three rows `hostscan` §5 marks "complex fold" — $\Gamma_0(5)$ pole 0 ($\lambda_2^{\rm norm}=11.18$),
$\Gamma_0(7)$ pole 0 ($7$), $\Gamma_0(9)$ pole 0 = Zagier B ($5.196$) — carry no period: they are
not merely badly scored, they are geometrically inert. *A host whose near singularity has a
complex-conjugate partner at the same distance is excluded.*

### 6.3 Pole at one member of the conjugate pair: three $K$-rational hosts  **[verified]** (`05_rows.gp`)

Taking $P_j$ = one of the conjugate pair gives $\lambda=(-\mu,\;\bar\mu-\mu)$, a genuinely
$K$-rational configuration with $K=\mathbf Q(\mu)$ imaginary quadratic:

| host | $K$ | $\lambda_1$ | $\lambda_2$ | $\lvert\lambda_1\rvert$ | $\lambda_2^{\rm norm}=\lvert\lambda_2\rvert$ | unit? | minimal recurrence | CDT-shape? |
|---|---|---|---|---|---|---|---|---|
| $X_0(5)$, pole at $\mathbf Z/2$ point, $F\in M_2(\Gamma_0(5))$ | $\mathbf Q(i)$ | $11-2i$ | $-4i$ | $\sqrt{125}=11.180$ | $4$ | no ($N=16$) | **5 terms**, deg 4 | **no** |
| $X_0(7)$, pole at $\mathbf Z/3$ point, $F\in M_1(\Gamma_1(7),\chi_{-7})$ | $\mathbf Q(\sqrt{-3})$ | $\tfrac{13-3\sqrt{-3}}2$ | $-3\sqrt{-3}$ | $7$ | $3\sqrt3=5.196$ | no ($N=27$) | **4 terms**, deg 2 | **no** |
| $X_0(7)$, same pole, $F\in M_2(\Gamma_0(7))$ | $\mathbf Q(\sqrt{-3})$ | " | " | $7$ | $5.196$ | no | **5 terms**, deg 4 | **no** |
| $X_0(9)$, pole at the cusp $1/3$, $F\in M_1(\Gamma_1(9),\chi_{-3})$ | $\mathbf Q(\sqrt{-3})$ | $\tfrac{9+3\sqrt{-3}}2$ | $3\sqrt{-3}$ | $3\sqrt3$ | $3\sqrt3=5.196$ | no ($N=27$) | **3 terms**, deg 4 | **yes, but $\lvert\lambda_1\rvert=\lvert\lambda_2\rvert$** |

Calibration of the same machinery on $\Gamma_0(6)$: pole at the cusp $1/2$ reproduces
**Zagier C = CDT's own row $1,3,15,93,639,4653$** with a three-term recurrence of degree 2 and
characteristic polynomial $y^2-10y+9$ (roots $9,1$); pole at $1/3$ gives Zagier A
($3y^2-21y-24$, roots $8,-1$); pole at $0$ gives Zagier F ($y^2+17y+72$, roots $-9,-8$).
$X_0(9)$ pole at cusp $0$ reproduces Zagier B ($y^2+9y+27$). **[verified]**

### 6.4 Why the elliptic placements are not CDT-shape, and why $X_0(9)$'s is dead  **[verified/exact]**

*Elliptic points.* When a **finite** singularity of the coordinate is an elliptic point of order
$e\ge2$ of the covering group, the Wronskian $W=F^2/(dx/d\tau)$ of the local system acquires a
fractional order there ($dx/d\tau$ vanishes to order $(e-1)/e$ in the $x$-coordinate), the
leading coefficient of the Fuchsian operator picks up extra factors of $(x-x_e)$, and the
recurrence lengthens. Measured exactly: $X_0(5)$ pole at cusp $0$ has characteristic polynomial
$y^4+44y^3+734y^2+5500y+15625=(y^2+22y+125)^2$ — the two singularities *doubled*; $X_0(7)$ pole
at cusp $0$ gives $(y^2+13y+49)^2$. This is the same phenomenon as `hostscan` §7's off-shape
configurations, and it is why the CDT-shape hosts of the census are exactly those whose finite
singularities are cusps (weight-one family) or Atkin–Lehner folds absorbed by a $W_N$-eigen $F$
(Fricke family) — never bare elliptic points of the covering group.

*The equilateral obstruction.* For a conjugate pair $\{\mu,\bar\mu\}$ and the third special
point at $\mu=0$, the placement at $\mu$ gives $\lambda=(-\mu,\;-2i\,\mathrm{Im}\,\mu)$, so
$$|\lambda_1|=|\lambda_2|\iff|\mu|=|\mu-\bar\mu|\iff\{0,\mu,\bar\mu\}\text{ is equilateral}
\iff\arg\mu=\pm150^\circ .$$
$X_0(9)$: $\mu=-\tfrac92-\tfrac{3\sqrt3}2i$, $\arg\mu=210^\circ$ — equilateral, $|\lambda_1|=|\lambda_2|=3\sqrt3$,
dead. $X_0(5)$ ($\sqrt{125}$ against $4$) and $X_0(7)$ ($7$ against $3\sqrt3$) are non-degenerate
— and are the two that fail the three-term test instead. **All three conjugate-pair curves are
excluded, each for a different reason, and no fourth one exists at $N\le60$.**

*And even had one survived,* its period would have been real: on $X_0(5)$, $X_0(7)$, $X_0(9)$ the
nebentypus is of order $\le2$ (§3.3), the Eisenstein source has rational coefficients, and by
`CUSP_PERIODS.md` Thm 3.1 the period at the cusp $0$ ($\zeta=e(0)=1$, all $w_\varphi$ real) is
real — the slots are $\zeta(3)$, $L(2,\chi_{-7})$, $L(3,\chi_{-3})$. By §2.1 nothing is gained.

### 6.5 Root-of-unity $\lambda_2$  **[exact]**

By §2.2 an Apéry-perfect imaginary quadratic host needs $\lambda_2\in\mu_K$. None of the
configurations above has it ($N(\lambda_2)=16$, $27$, $27$). The nearest thing in the whole
project is `hostscan` §7's $N=18$, $C=6$, $B=5$ configuration, whose extra characteristic roots
are a *primitive sixth root of unity, doubled* — a genuine $\mathcal O_{\mathbf Q(\sqrt{-3})}^\times$
unit far singularity — but it sits on a **seven**-singular-point geometry with a six-term,
degree-5 minimal recurrence, so the descent has nothing to remove there and the conditional
function would have to be regular at seven points.

---

## 7. Atkin–Lehner, Fricke and weight-four sources  **[exact]**

* **Atkin–Lehner quotients of $\Gamma_0$-type curves** ($\bar H=(\mathbf Z/N)^\times$) inherit the
  rule of §3.1: $\operatorname{ord}\varepsilon\le2$, whatever the quotient. (A finer, and
  independently sufficient, constraint: for $f\in M_k(N,\varepsilon)$, $f|W_Q\in M_k(N,\bar\varepsilon_Q\varepsilon_{N/Q})$,
  so the Fricke involution $W_N$ requires $\varepsilon=\bar\varepsilon$, i.e. quadratic; a partial
  $W_Q$ requires $\varepsilon_Q$ quadratic. And for **odd** weight $W_Q$ carries a factor $Q^{k/2}$,
  so its eigenvectors live over $F(\sqrt Q)$ — a *real* quadratic extension, which would
  reintroduce a second archimedean place and its tax, exactly the mechanism §1 says imaginary
  quadratic fields avoid.)
* **Atkin–Lehner quotients of the non-$\Gamma_0$ curves that do admit order $\ge3$.** $W_5$ on
  $X_1(5)$, $W_7$ on $X_1(7)$, $W_9$ on $X_1(9)$, $W_{13}$, $W_{16}$, $W_{25}$ are all Fricke
  involutions and send $\varepsilon\mapsto\bar\varepsilon$, killing every order-$\ge3$ character.
  The single survivor is $X_1(10)+W_2$ (since $\varepsilon_2=\mathbf 1$ for $\varepsilon=\psi_4$),
  and it is excluded on geometry: $\mu'=18$, $g=0$, cone orders $\in\{2,4\}$, so
  $A=\mu'/6+2-\nu'_\infty=5-\nu'_\infty$ and $\#\text{special}=\nu'_\infty+\#\text{cone}\ge\nu'_\infty+\tfrac43A=(20-\nu'_\infty)/3\ge5$
  for every admissible $\nu'_\infty\le5$. **Five or more special points, never four.**
* **Weight-four sources ($r=3$) with even characters of order 3, 4, 6 on Fricke-type hosts.** On a
  Fricke host the parametrising form $F=\mathcal D\log u$ has trivial character, so the source must
  too; a trivial-nebentypus Eisenstein direction is $E_4^{\psi,\bar\psi}(d\tau)$, which needs
  $\operatorname{cond}(\psi)^2d\mid N$. Order 4 has $\operatorname{cond}\psi\ge5$, hence $N\ge25$;
  orders 3 and 6 have $\operatorname{cond}\psi\ge7$, hence $N\ge49$. No four-point Atkin–Lehner
  quotient exists at $N=25,49,50$ or $81$ (`hostscan` §2's list at $N\le60$ is
  $5,6,7,8,9,10,11,12,13,14,15,17,18,20,21,22,30$). **Empty.**

---

## 8. Scoring  **[verified]** (`10_score.py`, using `lattice/cdt_finder/cdt_bound.py`)

CDT's inventory ($m=14$, $u=(1,3)$ for $k=2$, $(1,3,5)$ for $k=3$, $b_j=2$, $\sum e_i=6$;
$u_1=1$ forced by `INVENTORY_BOUND.md`), $\tau=4.2355$ ($k=2$) / $5.980$ ($k=3$),
ceil $=\log(256/\lambda_2^{\rm norm})$, entryR $=$ ceil $+\log0.62922-\tau$,
margin $=14\,\text{entryR}-(11.845+\log s)$.

```
host                                                 K                    l2n  k    tau   ceil  entryC  entryR   margin  shape live
CDT / Zagier C   Gamma_0(6) pole 1/2                 Q                 1.0000  2  4.235  5.545  +1.310  +0.846   +0.005  yes   yes
Zagier D  Gamma_1(5) pole 2/5                        Q(sqrt5)          1.0000  2  4.235  5.545  +1.310  +0.846   +0.005  yes   yes   (-1.44 after the 2nd real place)
HYPOTHETICAL imag. quad. host, lam2 a root of unity  Q(i) / Q(zeta_3)  1.0000  2  4.235  5.545  +1.310  +0.846   +0.005  yes   yes   *** does not exist ***
X_0(5) pole at ell_2   (weight 2, k=3)               Q(i)              4.0000  3  5.980  4.159  -1.821  -2.285  -42.445  no    yes
X_0(5) pole at ell_2   (if k=2 were reached)         Q(i)              4.0000  2  4.235  4.159  -0.077  -0.540  -18.017  no    yes
X_0(7) pole at ell_3   (weight 1, chi_-7, k=2)       Q(sqrt-3)         5.1962  2  4.235  3.897  -0.338  -0.801  -21.418  no    yes
X_0(7) pole at ell_3   (weight 2, k=3)               Q(sqrt-3)         5.1962  3  5.980  3.897  -2.083  -2.546  -45.846  no    yes
X_0(9) pole at cusp 1/3 (weight 1)                   Q(sqrt-3)         5.1962  2  4.235  3.897  -0.338  -0.801  -21.418  yes   NO   |lam_1|=|lam_2|
X_0(5)/X_0(7)/X_0(9) pole at cusp 0                  lam in imag.quad.   --    --    --     --      --      --       --   --    NO   lam_1 = conj(lam_2)
Zagier E  Gamma_0(8)  (Catalan, reference)           Q                 4.0000  2  4.235  4.159  -0.077  -0.540  -18.017  yes   yes
```

**Best margin over every imaginary-quadratic candidate: $-18.0$**, and that entry is itself
hypothetical (it assumes a $k=2$ row on $X_0(5)$, which does not exist — the honest $X_0(5)$
number is $-42.4$). The best *live, CDT-shape* imaginary-quadratic configuration has no Apéry
limit at all. Meanwhile the thing that is missing — an Apéry-perfect host over $\mathbf Q(i)$ or
$\mathbf Q(\zeta_3)$ with a non-real period — would have scored **exactly CDT's $+0.005$ with no
number-field correction whatsoever**, and would have proved a statement of the form
$L(2,\psi_4)\notin\mathbf Q(i)$. That is the precise size of what §§3–6 rule out.

---

## 9. Coverage — honest ledger

**Covered exhaustively and exactly.**
1. All subgroups $H\le(\mathbf Z/N)^\times$, all $N\le60$: genus, $\nu_2,\nu_3,\nu_\infty$, special
   point count (`01_groups.py`). Complete; 25 genus-zero curves, six four-point.
2. All nebentypus characters available on each of those 25 curves, with orders and parities
   (`03_chars.py`). Complete.
3. The $\mu$-configuration of every four-point curve, exactly (`02_geom.gp`, 80 digits,
   recognised by `bestappr`), hence all $3$ pole placements per curve, hence all
   $K$-rational coordinates for $K$ imaginary quadratic.
4. The minimal recurrence (number of terms, degree, characteristic polynomial) of the natural row
   for every such placement, with a calibration that reproduces Zagier A/C/F and B (`05_rows.gp`).
5. The fold-regular rank on the $\psi_4$-component of $M_3^{\rm Eis}(\Gamma_1(5))$, exactly
   (`04_g15.gp`), and on every order-$\{3,4,6\}$ component at $N\le30$, both weights
   (`09_census.gp`, `08_verify.gp`).

**Covered by argument, not by machine.**
* The Atkin–Lehner branch (§7): the character rule of §3.1 applies to any quotient of a
  $\Gamma_0$-type curve, and the six remaining cases are handled by the Fricke rule
  $\varepsilon\mapsto\bar\varepsilon$ plus one orbifold-defect inequality for $X_1(10)+W_2$. The
  four-point Atkin–Lehner list itself is **[cited]** from `hostscan` §2 and was not recomputed.
* The elliptic-point lemma of §6.4 is stated as the Wronskian/exponent mechanism and *verified*
  case by case (the perfect-square characteristic polynomials), not proved in general.
* CDT's contour loss $0.62922$ and Bost–Charles numerator $11.845$ are transported unchanged, as
  everywhere in this project; no contour is designed here.

**Not covered — explicit gaps.**
* **$N>60$.** The four-point census, the character census and the supply census all stop at
  $N\le60$ ($N\le30$ for the PARI supply ranks). Nothing suggests a four-point genus-zero curve
  exists at larger $N$ (the index grows and the special points with it), but it is not proved here.
* **Non-Eisenstein sources.** Cusp-form sources over an imaginary quadratic field were not swept;
  by §3.1 they face the same nebentypus obstruction, and by `CUSPFORM_SOURCES.md` the periods
  $L(f,r)$ on the four-point hosts are real in any case.
* **Weight $\ge5$ / $r\ge4$.** Only $r=2$ (weight-three sources) and $r=3$ (weight-four) were
  swept, matching the CDT architecture; higher $r$ costs $\approx1.7$ nats each and is hopeless
  from the table of §8.
* **Higher-degree CM fields.** $\mathbf Q(\zeta_5)$, $\mathbf Q(\zeta_8)$ and their kin are
  totally imaginary of degree 4: two complex places, so the bound averages *two* geometries and
  the no-tax statement of §1 does **not** extend to them unless the two places see the same host
  data. §7 notes where such a field would appear (odd-weight Atkin–Lehner eigenvectors); it was
  not pursued because the character rule closes those cases first.
* No independence check of the would-be function inventory, no contour, and no claim of
  irrationality anywhere.

---

## 10. Scripts and data

| file | what it does |
|---|---|
| `01_groups.py` | all $H\le(\mathbf Z/N)^\times$, $N\le60$: $\mu,\nu_2,\nu_3,\nu_\infty$, genus, #special $\to$ `01_groups.out`, `01_groups.json` (25 genus-zero curves, 6 four-point) |
| `02_geom.gp` | exact $\mu$-configuration of every four-point curve (cusps at $\gamma(60i)$, elliptic points directly, 80 digits) $\to$ `02_geom.out` |
| `03_chars.py` | nebentypus characters available on each of the 25 genus-zero curves, with order and parity $\to$ `03_chars.out` |
| `04_g15.gp` | the decisive computation: $\dim(\text{fold-regular}\cap\psi_4\text{-component})=0$ on $\Gamma_1(5)$, and $=2$ on the full space $\to$ `04_g15.out` |
| `05_rows.gp` | the rows on every pole placement of the four-point curves, minimal recurrence and characteristic polynomial; calibration on Zagier A/C/F/B $\to$ `05_rows.out` |
| `06_wt1.gp` | first pass at the weight-one rows on $X_0(7)$, $X_0(9)$ (superseded by `05_rows.gp`) $\to$ `06_wt1.out` |
| `07_supply.gp` | closed-form (direction-counting) version of the supply census $\to$ `07_supply.out` |
| `08_verify.gp` | PARI check of the two constant-term rules on a sample of $(N,k,\varepsilon)$ $\to$ `08_verify.out` |
| `09_census.gp` | exact PARI supply census, $N\le30$, both weights, all $\varepsilon$ of order 3, 4, 6 $\to$ `09_census.out` |
| `11_level10.gp` | the explicit first $\mathbf Q(i)$-rational fold-regular source, level 10 weight 3, and its non-real period $L(2,\psi_4)$ $\to$ `11_level10.out` |
| `10_score.py` | CDT scoring of every candidate + the no-tax comparison $\to$ `10_score.out`, `10_scored.json` |
