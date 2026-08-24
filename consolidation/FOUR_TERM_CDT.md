# The honest CDT entry and margin for the four-term Catalan hosts

*Claude (Opus 5), 2026-08-24.  Scripts, logs and data: `lattice/four_term_cdt/`.
Sequel to `FOUR_TERM_DEEP.md` §6 (the six mixed-exponent five-point rows and their
period matrices), `CATALAN_OBSTRUCTION.md` §3 (evasion (E1)), `CDT_FINDER.md` §§1–4
(the scoring conventions), `ADELIC_HOLONOMY.md` §2.6 (Theorem A, the adelic entry
condition), `CDT_NONCONGRUENCE.md` §10 (the $c$-ladder) and
`CATALAN_THREE_PERIOD.md` §1 (the fold-regularity test).*

*Tags: **[proved]** = exact argument; **[verified]** = exact / high-precision
computation over a stated range; **[measured]** = numerical at stated precision;
**[estimate]** = an input transported from CDT's own host, flagged at each use.*

**No irrationality claim is made anywhere in this document.**

---

## 0. Verdict

These are the first Catalan-class hosts *outside* the hypotheses of the width law of
`CATALAN_OBSTRUCTION.md` — they are not covers of $\mathbf P^1\setminus\{0,1,\infty\}$,
their singular points are irrational, and two of the five are order-$2$ orbifold points
rather than cusps.  So the ceiling $\log(16|t_2|)$ had to be replaced by an honest
numerical Fuchsian uniformisation.  It was.

| claim | verdict |
|---|---|
| The post-hypothesis geometry of every one of the six rows is a **four**-point orbifold on $\mathbf P^1$ (cusp at $0$, cusp at $\infty$, plus two more), never a three-punctured sphere; its conformal size at the cusp $0$ is computed here by numerical Fuchsian uniformisation, validated to $10$–$14$ digits against five exactly known cases | **[verified]** §5 |
| **Fold-regularity holds on all six.**  $H=B-\xi A$ has its radius of convergence enlarged from $\vert t_{\rm fold}\vert$ to $\vert t_2^{\rm post}\vert$, by $0.069$–$2.634$ nats depending on the row | **[measured]**, $n\le260$, §3 |
| **$c=1$: there is exactly one conditional generator, and the three-period hypothesis buys nothing.** The admissible companions form the affine line $B+\mathbf QA$, so fold-regularity is *equivalent* to "$\xi\in\mathbf Q$" for that row's own $\xi$; and the further inhomogeneous solutions $L y=t^k$ ($k=2,3$) have fold constants that are **not** in $\mathbf Q+\mathbf Q\zeta(2)+\mathbf QG$ (resp. $+\mathbf QL(2,\chi_{-3})$) at height $10^{12}$ over $218$–$220$ digits, on all six rows | **[proved]+[verified]** §4 |
| Consequently **only row 1 ($\xi=\frac14G$) is activated by the hypothesis $G\in\mathbf Q$**; rows 2 and 3 need the specific relations $b/c=-\frac38$, $b/c=-\frac34$ in $a+b\zeta(2)+cG=0$, which a generic relation does not supply | **[proved]** §4.2 |
| Denominators: $a_n\in\mathbf Z$ and $\mathrm{den}(b_n)\mid[1,\dots,n]^2$ **sharply**, with no $n^e$ improvement, on all seven rows ($n\le120$) | **[verified]** §6 |
| A **denominator-free algebraic pure function exists on every row**: $\sqrt{1-\sigma t+\pi t^2}\in\mathbf Z[\![t]\!]$, the square root of the reciprocal of the quadratic factor of the characteristic cubic.  It is branched exactly at the two order-$2$ points — hence at the *fold* on rows 1,2,3,5,6, so it costs the whole gain of the hypothesis there.  On rows 4 and 7 the fold is the rational cusp and the function is free.  $\sqrt{1-rt}\in\mathbf Z[\![t]\!]$ iff $4\mid r$ (rows 1, 4, 7) | **[verified]** §6.2 |
| **Symmetrisation is unavailable.**  Only rows 1, 4, 7 admit an involution of $\mathbf P^1$ fixing $0$ and permuting the singular set; on 4 and 7 the fold *is* the involution centre (identified with $\infty$ in the quotient, so the hypothesis cannot delete it), and on row 1 the fold's involution partner survives, leaving a four-point quotient of size exactly $4$ — entry $-2.849$, *worse* than the unsymmetrised $-1.914$ | **[proved]+[verified]** §7 |
| **A $2$-adic gain exists on row 1 and nowhere else**: $v_2(a_n)=v_2(b_n)=\tfrac n2+O(1)$, so the *conditional* function is $2$-adically overconvergent with $R_2=2^{1/2}$, worth $+0.322$ nats at $m=14$.  This does not contradict Calegari: the mechanism $\mathrm{slope}(H)=\mathrm{slope}(A)$ still holds, it is $\mathrm{slope}(A)$ that is $\tfrac12$ rather than $0$.  All other rows have slope $0$ at every $p\le13$ | **[measured]** §8 |
| **Entry at the hard ceiling.**  $-1.914$ (R1), $-0.703$ (R2, R3), $-4.081$ (R4), $+0.008$ (R5, R6), against level $8$'s $-0.077$ | **[verified]+[estimate]** §9 |
| **No four-term host carrying $G$ beats $-0.077$.**  The two that do — R5, R6, by $+0.084$ — carry $L(2,\chi_{-3})$, where the target is already a theorem (CDT).  The row whose hypothesis is literally $G\in\mathbf Q$ is $1.84$ nats *behind* level 8 | **[verified]** §9 |
| Every margin is $\le-15$ at $m=14$ with CDT's own contour loss transported; with that loss even R5/R6 fail entry ($-0.456$) | **[estimate]** §9 |

**One sentence.**  *The mixed-exponent four-term hosts realise evasion (E1) — irrational
conjugate singular points, no width law — and the geometry does improve where the
period is a mixed form, but the improvement is spent on the wrong rows: the only row
whose hypothesis is $G\in\mathbf Q$ has its fold and its second singularity almost on
top of each other ($|t_{\rm fold}/t_2|=0.83$), so the hypothesis buys $0.19$ nats where
level $8$ buys $0.69$ and CDT buy $2.20$, and after the honest four-point uniformisation
its entry is $-1.91$ against level $8$'s $-0.077$.*

---

## 1. The rows and the exact geometry

All six live in the mixed class $(\rho_p,\rho_r;M,j_1,j_2)=(-\tfrac12,0;1,0,0)$,
$R(n)=Cn^2$, of `FOUR_TERM_DEEP.md` §6.4.  Write
$$(n+1)^2u_{n+1}=P(n)u_n-Q(n)u_{n-1}+R(n)u_{n-2},\qquad
  \chi(\lambda)=\lambda^3-a\lambda^2+d\lambda-g=(\lambda-r)(\lambda^2-\sigma\lambda+\pi),$$
$g=C$, $t_i=1/\lambda_i$.  The Picard–Fuchs operator is
$$L=t\,\mathcal R(t)\,\partial_t^2+\mathcal S(t)\,\partial_t+\mathcal V(t),\qquad
  \mathcal R=1-at+dt^2-gt^3,$$
of **order two** with **five** regular singular points $0,t_1,t_2,t_3,\infty$;
$LA=0$ and $LB=t$ for the companion ($b_0=0$, $b_1=1$).  `01_geometry.py`.

### 1.1 Local exponents and Jordan structure  **[proved]**

| point | exponents | difference | monodromy | Kodaira | Jordan block |
|---|---|---|---|---|---|
| $t=0$ | $0,0$ | $0$ | unipotent (log) | $I_n$ | $2$ |
| $t=1/r$ (the rational root) | $0,0$ | $0$ | unipotent (log) | $I_n$ | $2$ |
| $t=1/\lambda_\pm$ (the conjugate pair) | $0,-\tfrac12$ | $\tfrac12$ | semisimple, order $2$ in $\mathrm{PSL}_2$ | $III$ | $1$ |
| $t=\infty$ | $2,2$ | $0$ | unipotent (log) | $I_n$ | $2$ |

A double indicial root at a regular singular point of a second-order operator forces a
logarithm, so the three $I_n$ points are genuinely logarithmic and **no apparent
singularity is possible**; the Fuchs relation $\sum=3$ checks:
$(0{+}0)+(0{+}0)+2\cdot(0-\tfrac12)+(2{+}2)=3$.  The characteristic cubics have
discriminants $2048,\,32768,\,32768,\,768,\,192,\,49152$ — all non-zero, so **the three
finite singular points are simple and distinct**, and the "fold-regularity codimension
$=$ Jordan block size" criterion of `CDT_NONCONGRUENCE.md` §10 applies with block size
$1$ at the order-$2$ points (one linear condition) and $2$ at the cusps.

### 1.2 The table  **[verified]**

$|t|$ increasing; the **fold** is the singular point nearest $0$, the one that governs
$\lim b_n/a_n$.

| row | $(a,c,d,f,C)_r$ | $\xi=\lim b_n/a_n$ | $t_{\rm fold}$ | $t_2^{\rm post}$ | $t_3^{\rm post}$ | score |
|---|---|---|---|---|---|---|
| **R1** | $(16,8,48,0,-128)_{8}$ | $\tfrac14G$ | $\tfrac{\sqrt2-1}4=0.103553$ (orb$_2$) | $\tfrac18$ (cusp) | $-\tfrac{1+\sqrt2}4=-0.603553$ (orb$_2$) | $-4.0794$ |
| **R2** | $(14,8,28,4,8)_{2}$ | $\tfrac12G-\tfrac3{16}\zeta(2)$ | $\tfrac32-\sqrt2=0.085786$ (orb$_2$) | $\tfrac12$ (cusp) | $\tfrac32+\sqrt2=2.914214$ (orb$_2$) | $-2.6931$ |
| **R3** | $(6,4,-32,-8,32)_{-2}$ | $\tfrac38\zeta(2)-\tfrac12G$ | $\tfrac{\sqrt2-1}4=0.103553$ (orb$_2$) | $-\tfrac12$ (cusp) | $-0.603553$ (orb$_2$) | $-2.6931$ |
| **R4** | $(16,8,68,8,32)_{8}$ | $\tfrac1{32}(2\zeta(2)+15L)$ | $\tfrac18$ (**cusp**) | $1-\tfrac{\sqrt3}2=0.133975$ (orb$_2$) | $1+\tfrac{\sqrt3}2=1.866025$ (orb$_2$) | $-4.0101$ |
| **R5** | $(17,10,32,8,16)_{1}$ | $\tfrac1{16}(15L-6\zeta(2))$ | $\tfrac12-\tfrac{\sqrt3}4=0.066987$ (orb$_2$) | $\tfrac12+\tfrac{\sqrt3}4=0.933013$ (orb$_2$) | $1$ (cusp) | $-2.0693$ |
| **R6** | $(13,8,-13,-1,-1)_{-1}$ | $\tfrac18(2\zeta(2)-3L)$ | $7-4\sqrt3=0.071797$ (orb$_2$) | $-1$ (cusp) | $7+4\sqrt3=13.928203$ (orb$_2$) | $-2.0000$ |
| *R7* | $(8,4,32,8,64)_{4}$ | *no arch. limit* | $\tfrac14$ (cusp) | $\tfrac{1\pm i\sqrt3}8$, $|t|=\tfrac14$ | — | $-3.3863$ |

$L=L(2,\chi_{-3})$; score $=\log|t_2^{\rm post}|-k$ with $k=2$.  R7 has all three
$|t_i|=\tfrac14$ (no dominant root, hence no Apéry limit) and is carried through the
computations as a control only.

> **The structural point that decides everything below.**  On R1 and R4 the fold and the
> next singularity nearly coincide, $|t_{\rm fold}/t_2^{\rm post}|=0.828$ and $0.933$; on
> R5 and R6 they are far apart, $0.0718$ and $0.0718$.  The hypothesis buys
> $\log|t_2^{\rm post}/t_{\rm fold}|$ nats, and **R1 — the only row whose hypothesis is
> $G\in\mathbf Q$ — is the worst of the six in exactly this respect.**

### 1.3 Periods re-verified

The full period matrices of `FOUR_TERM_DEEP.md` §6.5 were recomputed independently with
the Frobenius + Taylor-continuation code of `lattice/four_term_deep/11_foldmix.py` (the
machinery of `lattice/k3_period`), at two working precisions, in
`lattice/four_term_cdt/02_periods.sh`.  Agreement between the two precisions is
$134$–$136$ digits at all $21$ of the three finite singular points of the seven rows, and
the fold constants reproduce the identifications of `FOUR_TERM_DEEP.md` §6.4–6.5
exactly — e.g. R1 $\xi(t_3)=0.228991398544304753763650878733\ldots=\tfrac14G$ with
$\operatorname{Im}\xi=O(10^{-394})$, R2 $\xi(t_3)=0.149557659554567050688723913720\ldots
=\tfrac12G-\tfrac3{16}\zeta(2)$.  **[verified]**

The fold constant at the fold is real on all six rows.  At the other two finite points it
is complex whenever the continuation from $0$ has to cross a singularity, with
$\operatorname{Im}\xi\in\mathbf Q\zeta(2)$ resp. $\mathbf Q\pi^2\sqrt3$ (R2:
$\tfrac9{10}\zeta(2)$ at $t=\tfrac32+\sqrt2$ and $\tfrac38\zeta(2)$ at $t=\tfrac12$;
R4: $\tfrac5{96}\pi^2\sqrt3$ and $\tfrac1{96}\pi^2\sqrt3$; R5: two non-zero
imaginary parts).  On R1 and R3 the third point lies on the ray *opposite* the fold and
its constant is real as well:
$$\text{R1: }\xi(-\tfrac{1+\sqrt2}4)=\tfrac18(2G-3\zeta(2)),\qquad
  \text{R3: }\xi(-\tfrac12)=-\tfrac12G-\tfrac38\zeta(2)$$
(both **[verified]** to $135$–$136$ agreeing digits).  These are *different* elements of
$\mathbf QG+\mathbf Q\zeta(2)$ from the row's own $\xi$, so a single $H=B-\xi A$ can
never be regular at two singular points at once; deleting two would need a
two-dimensional conditional space, i.e. $c=2$, which §4 rules out.  The fold is therefore
forced to be the nearest singular point.

---

## 2. What the hypothesis removes, and what it does not

$H=bB-aA$ is regular at the fold under $\xi=a/b$.  It is *not* regular at the
Galois-conjugate order-$2$ point: Proposition G of `CDT_NONCONGRUENCE.md` §10.2
constrains the singular locus of the **minimal operator**, not the branching of a single
solution, and the level-$5$ host is the precedent — there too the conditional function is
regular at $x_1$ and branched at its conjugate $x_2$.  The numbers confirm it: the
radius of $H$ stops exactly at $|t_2^{\rm post}|$ (§3), not beyond.

So the post-hypothesis domain is
$$\Omega=\mathbf P^1\setminus\{0,\ t_2^{\rm post},\ t_3^{\rm post},\ \infty\}
  \quad\text{as an \emph{orbifold}: order-2 cone points where the exponent difference is }\tfrac12 .$$
This is a **four**-point orbifold, never CDT's three-punctured sphere.  The rigorous
monotonicity ceiling (fill one of the two finite points in, Schwarz–Pick) is
$16\min(|t_2^{\rm post}|,|t_3^{\rm post}|)=16|t_2^{\rm post}|$, and that is what the
"score $+\log16$" arithmetic of `FOUR_TERM_DEEP.md` amounts to.  It is an over-estimate
by $0.41$–$2.56$ nats (§5).

---

## 3. Fold-regularity: the radius-enlargement test  **[measured]**

`03_foldreg.py`.  $a_n,b_n$ exactly over $\mathbf Q$ to $n=260$; $\xi$ at $260$ digits
from its closed form; $h_n=b_n-\xi a_n$; the reported quantity is $\log|h_n|/n$, whose
limit is $-\log|t|$ for the radius $|t|$.

| row | target $-\log|t_{\rm fold}|$ | target $-\log|t_2^{\rm post}|$ | $\log|h_n|/n$ at $n=65,130,195$ | control ($\xi$ perturbed by $10^{-8}$) at $n=260$ | **enlargement** |
|---|---|---|---|---|---|
| R1 | $+2.2677$ | $+2.0794$ | $1.99997,\ 2.03458,\ 2.04750$ | $2.18198$ | $+0.1882$ |
| R2 | $+2.4559$ | $+0.6931$ | $0.60365,\ 0.64312,\ 0.65773$ | $2.36626$ | $+1.7627$ |
| R3 | $+2.2677$ | $+0.6931$ | $0.62577,\ 0.65397,\ 0.66492$ | $2.17718$ | $+1.5745$ |
| R4 | $+2.0794$ | $+2.0101$ | $1.96982,\ 1.98738,\ 1.99387$ | $1.99013$ | $+0.0693$ |
| R5 | $+2.7033$ | $+0.0693$ | $0.02953,\ 0.04672,\ 0.05320$ | $2.61247$ | $+2.6339$ |
| R6 | $+2.6339$ | $0$ | $-0.09795,\ -0.05425,\ -0.03823$ | $2.54282$ | $+2.6339$ |

Every row converges to the *post-fold* rate and every control reverts to the *fold*
rate.  So $B-\xi A$ is genuinely fold-regular at an **order-$2$ orbifold point** on rows
1, 2, 3, 5, 6 (exponents $0,-\tfrac12$: the hypothesis kills the
$(t-t_{\rm fold})^{-1/2}$ component, one linear condition, Jordan block size $1$) and at
a **logarithmic cusp** on row 4 (Jordan block size $2$, but the inhomogeneity is regular
there and the condition is again one-dimensional).  `CATALAN_AL_HOSTS.md` §5's
conclusion that "on a quotient whose fold is an order-$2$ orbifold point there is no
Eichler period and no conditional function" is **false on these hosts**: there is one,
and it is exactly as good a conditional function as a cusp fold gives.

---

## 4. How many conditional generators?  $c=1$.

### 4.1 The count  **[proved]**

The holonomic rank is $2$ and $LB=t$, so the $\mathbf Q(x)$-module generated by one
conditional $H$ is $\{H,\theta H\}$ together with the constant $1$: $\theta^2H$ is a
$\mathbf Q(x)$-combination of $H$, $\theta H$ and $1$.  **The $\theta$-orbit of one
generator has size $2$**, against CDT's $4$ derivatives $+3$ integrals $=7$
(`CATALAN_THREE_PERIOD.md` §4.1); adding the two cheap integrations
$\int H\,dx/x$, $\int\theta H\,dx/x$ (cost $e_i=1$, i.e. $\tau^\sharp$ only) gives at
most $4$ conditional functions here.

The admissible power-series companions with rational coefficients form the affine line
$B+\mathbf QA$ (the second solution of $L$ at $0$ is logarithmic and is not a power
series).  The fold constant is $\xi(B+\alpha A)=\xi-\alpha$.  Hence
$$\boxed{\ \exists\ \text{a fold-regular rational combination}\iff\xi\in\mathbf Q\ }$$
— one condition, one generator, $c=1$.

The $c$-ladder of `CDT_NONCONGRUENCE.md` §10.7 gives $c\le\#\text{cusps}-1$; these hosts
have **three** cusps ($0$, $1/r$, $\infty$; the two order-$2$ points are not cusps), so
$c\le2$.  The realised value is $1$.

### 4.2 The three-period hypothesis buys nothing  **[verified]**

The natural candidates for a second direction are the further inhomogeneous solutions
$y^{(k)}$ of $Ly=t^k$ — in recurrence terms $b_0=\dots=b_{k-1}=0$, $b_k=1$ — which are
new $\mathbf Q[\![x]\!]$ functions with the same singular set, exactly analogous to
CDT's $H_A,H_B,H_C$.  A fold-regular rational combination
$\alpha y^{(1)}+\beta y^{(2)}+\gamma A$ exists iff $1,\xi^{(1)},\xi^{(2)}$ are
$\mathbf Q$-linearly dependent.  `09_second.py` ($n=1200$, $220$ digits, PSLQ at height
$10^{12}$):

| row | $\xi^{(1)}$ | identified | $\xi^{(2)}$ | in $\mathbf Q+\mathbf Q\zeta(2)+\mathbf Q(G\text{ or }L)$? |
|---|---|---|---|---|
| R1 | $0.2289913985\ldots$ | $\tfrac14G$ ✅ | $0.0352787486\ldots$ | **no** (219 digits, `09b`) |
| R2 | $0.1495576596\ldots$ | $\tfrac12G-\tfrac3{16}\zeta(2)$ ✅ | $0.0168853302\ldots$ | **no** (220 digits) |
| R3 | $0.1588674780\ldots$ | $\tfrac38\zeta(2)-\tfrac12G$ ✅ | $0.0205750248\ldots$ | **no** (219 digits) |
| R5 | $0.1156207370\ldots$ | $\tfrac1{16}(15L-6\zeta(2))$ ✅ | $0.0101221354\ldots$ | **no** (219 digits) |
| R6 | $0.1182451119\ldots$ | $\tfrac18(2\zeta(2)-3L)$ ✅ | $0.0108898665\ldots$ | **no** (219 digits) |
| R4 | $0.4690438852\ldots$ | $\tfrac1{32}(2\zeta(2)+15L)$ ✅ | $0.1084842159\ldots$ | **no** (218 digits, `09b`) |

($\xi^{(3)}$ likewise.  R1 and R4, whose ratios $|t_{\rm fold}/t_2^{\rm post}|$ are
$0.828$ and $0.933$, converge too slowly at $n=1200$; `09b_second_slow.py` re-runs them at
$n=9000$ and confirms the same negative at **$218$–$219$ digits**, so the verdict is
uniform over all six rows.)  The recovery of $\xi^{(1)}$ in every case is the self-test.

**Consequence.**  On these hosts the hypothesis that activates a conditional function is
"$\xi\in\mathbf Q$" for that row's own $\xi$, and nothing weaker:

* $\xi_{\rm R1}=\tfrac14G\in\mathbf Q\iff G\in\mathbf Q$.
* $\xi_{\rm R2}=\tfrac12G-\tfrac3{16}\zeta(2)\in\mathbf Q$: under a relation
  $a+b\zeta(2)+cG=0$ this holds iff $b/c=-\tfrac38$; for R3 iff $b/c=-\tfrac34$.
  A *generic* relation among $1,\zeta(2),G$ activates none of the three rows.
* So **under $G\in\mathbf Q$ alone: one conditional generator, on R1 only.  Under the
  three-period hypothesis: still one generator, and only on the row whose $\xi$-vector
  happens to be proportional to the relation.**

This is the four-term analogue of `CDT_NONCONGRUENCE.md` §10.7's "no admissible second
fold-regular class", reached by a different route: the second class *exists as a
function* ($y^{(2)}$), but its period is a new unidentified constant, so it carries no
hypothesis one wants to contradict.

---

## 5. The conformal size: numerical Fuchsian uniformisation

### 5.1 What has to be computed, and why the width law does not apply

`CATALAN_OBSTRUCTION.md` proves the width law only for covers of
$\mathbf P^1\setminus\{0,1,\infty\}$.  These hosts are not such covers: their singular
points are irrational, two of them are order-$2$ cone points, and the post-hypothesis
geometry has four special points.  The CDT quantity $\log|\varphi'(0)|$ is the **cusp
size** of the host orbifold: if $\Gamma$ uniformises $\Omega$ and the cusp over $x=0$ is
normalised to width $1$ ($\tau\mapsto\tau+1$, $q=e^{2\pi i\tau}$), then
$x=\varphi(q)=cq+O(q^2)$ and $|\varphi'(0)|=|c|$.  For
$\mathbf P^1\setminus\{0,s,\infty\}$ this is $16|s|$ (modular $\lambda$); for a
four-point orbifold there is no closed form and the accessory parameter must be solved
numerically.

### 5.2 The method (`05_uniformise.py`)  **[verified]**

* Projective structure $u''+Q(z)u=0$,
  $Q=\sum_j\bigl[\tfrac{1-\theta_j^2}{4(z-z_j)^2}+\tfrac{B_j}{z-z_j}\bigr]$ with
  $\theta=0$ at a cusp and $\theta=\tfrac12$ at an order-$2$ cone point; regularity at
  $\infty$ ($Q\sim\tfrac1{4z^2}$) gives two linear conditions on the three $B_j$,
  leaving **one real accessory parameter**.
* Frobenius basis at the cusp $0$ (exponents $\tfrac12,\tfrac12$, logarithmic):
  $u_1=z^{1/2}v_1$, $u_2=u_1\log z+z^{1/2}w$ with $v_1(0)=1$, $w(0)=0$, so
  $W:=u_2/u_1=\log z+O(z)$, $\tau_0:=W/(2\pi i)$ satisfies $\tau_0\mapsto\tau_0+1$
  around $z=0$, and $e^{W}=z(1+O(z))$.
* Monodromy by Taylor-step continuation of $u''=-Qu$ along straight lassos from a base
  point placed **off the real axis** (a detour would conjugate a generator and destroy
  the test below).
* The uniformising coordinate is $\tau=\tau_0+\beta$ for a constant $\beta$, and the
  monodromy group is Fuchsian iff **every generator preserves one horizontal line**
  $\operatorname{Im}\tau_0=h$.  That is one real equation in the accessory parameter
  (plus $\operatorname{Im}c=0$ per generator as a check).  Then
  $$\boxed{\ |\varphi'(0)|=e^{-2\pi h}\ }$$
  since $q=e^{2\pi i\tau}=e^{2\pi i\beta}e^{W}$.
* The equation has a *ladder* of real-monodromy solutions with increasing size, only one
  of which is the uniformisation; the correct one is the **largest root not exceeding
  the rigorous monotonicity ceiling** $16\min_j|z_j|$.  This selection rule is validated
  below.

### 5.3 Validation  **[verified]**

| configuration | computed | exact | source of the exact value |
|---|---|---|---|
| $\{0,1,\infty\}$ all cusps | $15.999999999999977$ | $16$ | modular $\lambda$ |
| $\{0,1,-1,\infty\}$ all cusps | $4.000000000$, accessory $=\tfrac14$ | $4$ | $z\mapsto z^2$ cover of $\{0,1,\infty\}$ |
| $\{0,\tfrac19,1,\infty\}$ all cusps ($X_0(6)$, CDT's own host) | $1.000000000$, accessory $=-\tfrac5{16}$ | $1$ | $x=q\prod\frac{(1-q^n)^4(1-q^{6n})^8}{(1-q^{2n})^8(1-q^{3n})^4}=q+O(q^2)$ at the width-$1$ cusp |
| $\{0,\infty$ cusps$;\ 1$ order-$2\}$ | $63.99999999999873$ | $64$ | $\Gamma_0(2)$: $1/|v(\tfrac{1+i}2)|$ with $v=(\eta(2\tau)/\eta(\tau))^{24}$, $v(\tfrac{1+i}2)=-\tfrac1{64}$ |
| $\{0,\infty$ cusps$;\ \pm\tfrac12$ order-$2\}$ | $3.999999999999983$ | $4$ | $z\mapsto z^2$ cover of the previous, at scale $\tfrac12$ |

The fourth line also **reproduces CDT's symmetrised ceiling from scratch**: their
quotient is $\{0,\infty$ cusps; $4s$ order-$2\}$, size $64\cdot4s=256s$. ✅

### 5.4 The sizes  **[verified]**, `06_sizes.py`

| row | post-hypothesis orbifold | **size** $|\varphi'(0)|$ | $\log$ | ceiling $16|t_2^{\rm post}|$ | loss to the 4th point |
|---|---|---|---|---|---|
| R1 | $\{0,\infty$ cusps$;\ \tfrac18$ cusp$;\ -0.603553$ orb$_2\}$ | $1.050903594$ | $+0.04965$ | $2.000000$ | $0.6435$ |
| R2 | $\{0,\infty$ cusps$;\ \tfrac12$ cusp$;\ 2.914214$ orb$_2\}$ | $4.871667446$ | $+1.58344$ | $8.000000$ | $0.4960$ |
| R3 | $\{0,\infty$ cusps$;\ -\tfrac12$ cusp$;\ -0.603553$ orb$_2\}$ | $4.871667446$ | $+1.58344$ | $8.000000$ | $0.4960$ |
| R4 | $\{0,\infty$ cusps$;\ 0.133975,\ 1.866025$ orb$_2\}$ | $0.166183685$ | $-1.79466$ | $2.143594$ | $2.5571$ |
| R5 | $\{0,\infty$ cusps$;\ 0.933013$ orb$_2;\ 1$ cusp$\}$ | $9.918076446$ | $+2.29436$ | $14.928203$ | $0.4089$ |
| R6 | $\{0,\infty$ cusps$;\ -1$ cusp$;\ 13.928203$ orb$_2\}$ | $9.918076446$ | $+2.29436$ | $16.000000$ | $0.4782$ |

Two internal consistency checks fell out.  **R2 and R3 have exactly the same size**, and
so do **R5 and R6**: the map $x\mapsto x/(x-1)$ sends the normalised cross-ratio
$\mu=3+2\sqrt2$ (R2) to $\mu/(\mu-1)=\tfrac{1+\sqrt2}2$ (R3) and $\mu=8-4\sqrt3$ (R5) to
$-(7+4\sqrt3)$ (R6), so these are the *same orbifold in different coordinates*, and in
each pair the cusp sits at the same distance $|t|$ from $0$ — a coincidence the
computation reproduces to $10$ digits without being told.

R7's two outer points are a complex-conjugate pair, for which the "common horizontal
line" residual degenerates ($h_1\equiv h_2$ by conjugation symmetry); its size was not
computed.  R7 has no Apéry limit and is not needed.

---

## 6. Denominators and pure functions

### 6.1 The array $\mathbf b$  **[verified]**, `04_denom.py`, $n\le120$

For every one of the seven rows: $a_n\in\mathbf Z$; the least $k$ with
$[1,\dots,n]^k b_n\in\mathbf Z$ is $\mathbf 2$, **sharply**; and no pair $(k,e)$ with
$k<2$ works, i.e. $n^e[1..n]^kb_n\in\mathbf Z$ first at $(k,e)=(2,0)$.  So
$$\mathbf b=\bigl(\underbrace{0,\dots}_{u_1},1\bigr)\times\bigl(\underbrace{0,\dots}_{u_2},1\bigr),
\qquad \sigma_m=2,\qquad \text{no geometric part},$$
in the CDT normalisation of `CDT_FINDER.md` §1 — the four-term echo of $k=w+1=2$.  (For
comparison, CDT's own host is symmetrised, $b_j=2$, $\sigma_m=4$.)

### 6.2 Pure and algebraic functions  **[verified]**

The three finite singular points are **not** one Galois orbit — $1/r$ is rational and the
other two are conjugate — so the Proposition-G obstruction of `CDT_NONCONGRUENCE.md`
§10.2 does not apply a priori, and indeed:

* **the polylogarithm module** $\mathrm{Li}_j(rx)=\sum r^nx^n/n^j$ lives on
  $\mathbf P^1\setminus\{0,1/r,\infty\}$, has integral numerators ($r\in\mathbf Z$) and
  denominator type $n^{j}$ only ($b_{ij}=0$).  This is the exact analogue of CDT's pure
  module and is why the transported inventory is not absurd;
* $$\sqrt{1-\sigma t+\pi t^2}\in\mathbf Z[\![t]\!]\quad\text{on all seven rows}$$
  ($1-\sigma t+\pi t^2$ is the reciprocal of the quadratic factor of $\chi$; in each case
  it has the shape $1-4w$ with $w\in\mathbf Z[t]$, or is a square mod $2$).  Heads:
  R1/R3 $1-4t-16t^2-64t^3-\cdots$, R2 $1-6t-16t^2-96t^3-\cdots$,
  R6 $1-7t-24t^2-168t^3-\cdots$.  This is a **denominator-free** ($\sigma=0$) algebraic
  pure function — precisely the object that raises $u_j$.  **But it is branched exactly
  at the two order-$2$ points**, one of which is the *fold* on rows 1,2,3,5,6, so
  including it forces $\Omega$ to keep the fold and destroys the whole gain of the
  hypothesis (a re-run of the level-$5$ $\sqrt{1-44t-16t^2}$ instance of Proposition G,
  now with the loophole *visible*: on rows 4 and 7, whose fold is the rational cusp, the
  function is free);
* $\sqrt{1-rt}\in\mathbf Z[\![t]\!]$ iff $4\mid r$, i.e. on R1 ($r=8$), R4 ($r=8$) and
  R7 ($r=4$).  On R1 it is branched only at the *outer* cusp $\tfrac18$ and is therefore
  **admissible and denominator-free** — one extra free pure function.  It is worth
  $\Delta\tau^\flat=-\bigl[(2u_1{+}1)+(2u_2{+}1)\bigr]/m^2=-0.051$ at $m=14$: real, and negligible.

### 6.3 What a better inventory would be worth

$\tau^\flat=\sigma_m-\frac1{m^2}\sum_ju_j^2b_j$ with $\sigma_m=2$, $b=(1,1)$:
CDT's profile $u=(1,3)$ gives $\tau^\flat=1.948980$; the best conceivable profile
$u_j=m/2$ gives $\tau^\flat=1.5$.  With $\tau^\sharp=0.3375$ transported,
$$\tau_{\rm A}=2.286480\ \ (\text{CDT inventory}),\qquad
  \tau_{\rm B}=1.837500\ \ (\text{best conceivable}),$$
against $\tau=4.235459$ for the *symmetrised* level-$8$ and CDT hosts.  The entire
inventory question is worth at most $0.449$ nats here.

---

## 7. Symmetrisation is unavailable  **[proved]+[verified]**, `10_symm.py`

CDT's $+\log16$ comes from descending along $w(x)=sx/(x-s)$ (fixes $0$, swaps
$s\leftrightarrow\infty$), $y=x^2/(x-s)$: the cusp at $0$ acquires half its width and the
size goes $16s\mapsto 64\cdot(4s)=256s$; the price is $y\sim x^2$, i.e.
$[1..n]\mapsto[1..2n]$ and $\tau:2.2865\mapsto4.2355$.

Enumerating all involutions of $\mathbf P^1$ that fix $0$ and permute the singular set
type-preservingly:

| row | pre-hypothesis (5 points) | post-hypothesis (4 points) |
|---|---|---|
| R1 | $w(x)=\tfrac18x/(x-\tfrac18)$ | none |
| R2, R3, R5, R6 | **none** | none |
| R4, R7 | $w(x)=\tfrac18x/(x-\tfrac18)$, $\tfrac14x/(x-\tfrac14)$ | none |

* On **R4 and R7 the fold *is* the involution centre** $s$, which the quotient identifies
  with the cusp at $\infty$; the hypothesis cannot delete it.  Symmetrisation is not
  available.
* On **R1** the involution swaps the two order-$2$ points $0.103553\leftrightarrow
  -0.603553$ and $\tfrac18\leftrightarrow\infty$.  Its quotient has a cusp at $y=0$, a
  cusp at $y=\infty$, an order-$2$ branch point at $y=4s=\tfrac12$, and the merged
  order-$2$ pair at $y=-\tfrac12$.  $\mathrm{Sym}^+H=H(x)+H(w(x))$ is regular at the fold
  but **not** at the fold's involution partner, so $y=-\tfrac12$ survives: the symmetrised
  post-hypothesis geometry is $\{0,\infty$ cusps; $\pm\tfrac12$ order-$2\}$, whose cusp
  size is **exactly $4$** (the $z\mapsto z^2$ cover of $\{0,\infty;1$ order-$2\}$ of size
  $64$, at scale $\tfrac12$; confirmed numerically to $14$ digits).  Entry
  $=\log4-4.235459=\mathbf{-2.849}$, against the unsymmetrised $-1.914$.
  **Symmetrisation makes R1 worse by $0.935$ nats.**

That is the whole of level $8$'s advantage, stated structurally: level $8$'s
post-hypothesis geometry is a *three*-punctured sphere and its involution is compatible
with the fold's deletion, so it collects the full $\log16$; no four-term host does.

---

## 8. The $p$-adic places  **[measured]**, `07_padic.py`, $n\le160$, $p\le13$

Theorem A of `ADELIC_HOLONOMY.md` §2.6 needs $|c_{i,n}|_p\le Cn^AR_p^{-n}$ for **every**
member of the inventory, so $R_p=\min_i p^{\,\mathrm{slope}_i}$.

| row | slope of $a_n$ at $2$ | of $b_n$ | of $b_n-\tfrac37a_n$ | $R_2$ | slopes at $p=3,5,7,11,13$ | polylog module $\mathrm{Li}_j(rx)$ |
|---|---|---|---|---|---|---|
| **R1** | $+0.50108$ | $+0.50108$ | $+0.50108$ | $\mathbf{2^{1/2}}$ | $0$ | $2^{3}$ |
| R2 | $+0.00108$ | $+0.00108$ | $+0.00108$ | $1$ | $0$ | $2^{1}$ |
| R3 | $+0.00108$ | $+0.00108$ | $+0.00108$ | $1$ | $0$ | $2^{1}$ |
| R4 | $-0.00156$ | $-0.00156$ | $-0.00156$ | $1$ | $0$ | $2^{3}$ |
| R5 | $-0.01608$ | $-0.01608$ | $-0.01608$ | $1$ | $0$ | — |
| R6 | $+0.00707$ | $-0.02918$ | $-0.02918$ | $1$ | $0$ | — |
| R7 | $-0.01310$ | $-0.01310$ | $-0.01310$ | $1$ | $0$ | $2^{2}$ |

(the $O(10^{-2})$ entries are the $-2\log_2 n$ of the $[1..n]^2$ denominator over a
finite window; the slope set is $\{2\}$ or empty, exactly as `FOUR_TERM_DEEP.md` §6.8
found for $\sigma_p$.)

> **Row 1's $2$-adic gain is real, and it is the first one anywhere in the Catalan
> ledger.**  $v_2(a_n)=v_2(b_n)=\tfrac n2+O(1)$, so the *conditional* function inherits
> slope $\tfrac12$: $R_2=2^{1/2}$, worth $(1-\tfrac1m)\log2^{1/2}=+0.322$ nats at
> $m=14$.  This does **not** contradict `CATALAN_OBSTRUCTION.md` (2.3): the Calegari
> mechanism says $\mathrm{slope}(H)=\mathrm{slope}(A)$ because $\xi_2-\xi\ne0$
> ($\xi_2=\tfrac14\zeta_2(2)$, irrational), and that is exactly what is observed — it is
> $\mathrm{slope}(A)$ that is $\tfrac12$ instead of $0$.  The gain is also not a
> normalisation artefact: rescaling $x\mapsto\rho x$ moves $\log|\varphi'(0)|$ by
> $\log\rho$ and $\sum_p\log R_p$ by $-\log\rho$, and the entry functional carries the
> weight $1-\tfrac1m$ on the second, so $x$ as given is within $0.05$ nats of optimal.
> On every other row the $2$-adic place contributes nothing.

---

## 9. Entry and margin  `08_entry.py`

Conventions exactly those of `CDT_FINDER.md` §§1–4 and `ADELIC_HOLONOMY.md` §2.6:
$$\text{entry}=\log|\varphi'(0)|+\Bigl(1-\tfrac1m\Bigr)\sum_p\log R_p-\tau(\mathbf b;\mathbf e),
\qquad \text{margin}=m\cdot\text{entry}-\mathrm{BC}(\varphi).$$
$m=14$, $\tau^\sharp=27/80$ and the shape term $\mathrm{BC}-\log|\varphi'(0)|=6.763$ are
**transported from CDT's own host** and are the only *estimated* inputs; the sizes, the
denominators and the $p$-adic slopes are computed here.  "ceiling" uses the exact
orbifold size (a genuine supremum, not attained); "contour" applies CDT's own realised
loss $0.62922$ on top.

### 9.1 At the hard uniformisation ceiling

| row | $\xi$ | $\log$ size | adelic | $\tau_{\rm A}$ | **entry (ceiling)** | entry, best inventory $\tau_{\rm B}$ | vs level 8 |
|---|---|---|---|---|---|---|---|
| **R1** | $\tfrac14G$ | $+0.0497$ | $+0.3225$ | $2.2865$ | $\mathbf{-1.9143}$ | $-1.4653$ | $-1.838$ |
| **R2** | $\tfrac12G-\tfrac3{16}\zeta(2)$ | $+1.5834$ | $0$ | $2.2865$ | $\mathbf{-0.7030}$ | $-0.2541$ | $-0.627$ |
| **R3** | $\tfrac38\zeta(2)-\tfrac12G$ | $+1.5834$ | $0$ | $2.2865$ | $\mathbf{-0.7030}$ | $-0.2541$ | $-0.627$ |
| **R4** | $\tfrac1{32}(2\zeta(2)+15L)$ | $-1.7947$ | $0$ | $2.2865$ | $\mathbf{-4.0811}$ | $-3.6322$ | $-4.005$ |
| **R5** | $\tfrac1{16}(15L-6\zeta(2))$ | $+2.2944$ | $0$ | $2.2865$ | $\mathbf{+0.0079}$ | $+0.4569$ | $+0.084$ |
| **R6** | $\tfrac18(2\zeta(2)-3L)$ | $+2.2944$ | $0$ | $2.2865$ | $\mathbf{+0.0079}$ | $+0.4569$ | $+0.084$ |
| R1, symmetrised | $\tfrac14G$ | $+1.3863$ | $0$ | $4.2355$ | $-2.8492$ | — | $-2.773$ |

Benchmarks in the same conventions: **level 8 (Zagier $\mathbf E$), symmetrised,
$\mathbf{-0.0766}$**; level 8 unsymmetrised $-0.9002$; level 16 symmetrised $-1.4629$;
CDT's $X_0(6)$ symmetrised $+1.3097$.

### 9.2 With a realised contour, and the margins  **[estimate]**

| row | entry (CDT contour loss transported) | $\mathrm{BC}$ | **margin at $m=14$** |
|---|---|---|---|
| R1 | $-2.3774$ | $6.349$ | $-39.64$ |
| R2, R3 | $-1.1661$ | $7.883$ | $-24.21$ |
| R4 | $-4.5442$ | $4.505$ | $-68.13$ |
| R5, R6 | $-0.4552$ | $8.594$ | $-14.97$ |

Even R5 and R6, the only rows that clear entry at the ceiling, fail it as soon as any
contour loss is paid; and no margin is within $15$ of zero.

### 9.3 The answer to the question as asked

> **Does any four-term host beat level 8's $-0.077$?**
>
> **Yes — two of them, by $+0.084$ nats each: R5 and R6.  Both carry
> $L(2,\chi_{-3})$, where linear independence of $1,\zeta(2),L(2,\chi_{-3})$ is already
> CDT's theorem.  No host carrying $G$ does.**  The best Catalan four-term entry is
> $-0.703$ (R2 and R3, whose hypothesis is a specific mixed relation, not $G\in\mathbf Q$),
> $0.627$ nats behind level 8; and R1 — the *only* row activated by $G\in\mathbf Q$
> itself, even with its $+0.322$ of genuine $2$-adic overconvergence — is at $-1.914$,
> $1.838$ nats behind.

**The margin story, in one accounting.**  Take level $8$ ($-0.0766$) and walk to R1:

| step | $\Delta$ |
|---|---|
| level 8, symmetrised, at the ceiling | $-0.077$ |
| lose the symmetrisation ($\tau:4.2355\to2.2865$, size $64\to4$) | $-0.823$ |
| the fold moves from $|t|=\tfrac18$ to $|t|=0.1036$ and the outer point from $\tfrac14$ to $\tfrac18$: net ceiling $4\to2$ | $-0.693$ |
| the *fourth* special point (the outer order-$2$ point at $-0.6036$), which level $8$ does not have | $-0.644$ |
| the $2$-adic place, $R_2=2^{1/2}$ — the only credit | $+0.322$ |
| **R1** | $\mathbf{-1.914}$ |

The four-term world therefore does what `CATALAN_OBSTRUCTION.md` (E1) hoped in one
respect — the width law's hypothesis is gone, and the *score* $\log|t_2|-k$ improves by
$0.693$ on R2/R3 — but it pays for it three times over: no symmetrisation, an extra
puncture, and (on the row that matters) a fold that is barely inside its neighbour.

### 9.4 Against the $c$-ladder

`CDT_NONCONGRUENCE.md` §10.7 prices a host by $c$, the number of periods in one
hypothesised relation.  Here $c\le\#\text{cusps}-1=2$ and the realised value is
$c=1$ (§4).  A host with $c=3$ needs four cusps; these hosts have three, because two of
the five special points are order-$2$ orbifold points — the very feature that makes them
evade the width law is the feature that caps $c$ at $2$.  There is no configuration of
this class with more cusps: the class $(-\tfrac12,0;1,0,0)$ has exactly one rational
characteristic root by construction.

---

## 10. Status ledger

**[proved]**  §1.1 (exponents, Jordan blocks, Fuchs relation, no apparent
singularity); §2 (why Proposition G does not force the conjugate point out); §4.1
($c=1$, the $\theta$-orbit has size $2$, fold-regular $\iff\xi\in\mathbf Q$); §4.2 (which
relation activates which row); §7 (the involution enumeration and why the symmetrisation
fails on every row).

**[verified]**  §1.2–1.3 (geometry exactly, periods re-verified to $134$–$136$ agreeing
digits by an independent continuation); §3 (radius enlargement, $n\le260$); §4.2
($\xi^{(2)},\xi^{(3)}$ unidentified on all six rows, $218$–$220$ digits, height $10^{12}$); §5.3 (five
exact validations of the uniformiser, including $X_0(6)$ and the reconstruction of CDT's
$256s$); §5.4 (the six sizes, and the two coordinate-coincidence cross-checks); §6
(denominators sharp $k=2$, $n\le120$; the integral algebraic functions); §7 (the
symmetrised R1 size $=4$).

**[measured]**  §8 ($p$-adic slopes, $n\le160$, $p\le13$).

**[estimate]**  $m=14$ and CDT's $u=(1,3)$, $\tau^\sharp=27/80$ and the shape term
$6.763$ are transported from CDT's own host, exactly as in `CDT_FINDER.md` §8; the
honest module count here is smaller ($4$ conditional functions, not $7$), which moves
$\tau$ by at most $0.03$ and the margins by a few units.  No contour was designed and no
Bost–Charles integral was computed on any of these orbifolds.

**[open]**
1. The ladder-selection rule for the accessory parameter ("largest real-monodromy root
   not exceeding the monotonicity ceiling") is validated on four exactly known cases but
   not proved; the standard discriminator is an oscillation count and was not
   implemented.  A wrong choice could only make the sizes *smaller*, hence the entries
   *worse*.
2. R7 (complex-conjugate outer pair) was not uniformised; the residual degenerates by
   conjugation symmetry.  R7 has no Apéry limit.
3. The identification of $\xi^{(2)}$, $\xi^{(3)}$ — five new unidentified constants, one
   per row, attached to the second and third inhomogeneous solutions of the same
   operator.  They are the four-term analogue of the $\xi^*$ of `FOUR_TERM_DEEP.md` §6.7.
4. Whether the *pure* inventory transported from CDT is realisable here at all: no
   $\mathbf Q(x)$-independence check was run on $\{1,\sqrt{1-rt},\mathrm{Li}_j(rx),
   H,\theta H,\dots\}$.

---

## 11. Reproduction

```
cd /home/ubuntu/code/math-modular-sources/lattice/four_term_cdt
python3 01_geometry.py                 # exact geometry -> out/geometry.json
./02_periods.sh                        # period matrices re-verified -> out/02_periods.log
NMAX=260 python3 03_foldreg.py         # radius enlargement -> out/foldreg.json
NMAX=120 python3 04_denom.py           # denominators, algebraic functions -> out/denom.json
python3 05_uniformise.py test          # the five exact validations of the uniformiser
python3 06_sizes.py                    # the conformal sizes -> out/sizes.json
NMAX=160 python3 07_padic.py           # p-adic slopes -> out/padic.json
python3 08_entry.py                    # entry and margin -> out/entry.json
NMAX=1200 python3 09_second.py         # is there a second conditional direction?
NMAX=9000 ROWS=R1,R4 OUTJSON=second_slow.json python3 09b_second_slow.py
python3 10_symm.py                     # the involution / symmetrisation analysis
```
PARI/GP is not needed here (all the exact arithmetic is over $\mathbf Q$ in
`fractions`); `05_uniformise.py` uses numpy/`cmath` in double precision and mpmath is
used only for the periods and the PSLQ searches.
