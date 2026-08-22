# The adelic arithmetic holonomy bound, and what it does for Catalan

*Fable (Opus 5), 2026-08-22.  Scripts: `lattice/adelic_holonomy/`
(`adelic_bound.py`, `final_table.py`, `pure_2adic.py`, `rowE_sources.gp`,
`hosts_padic.gp`).  Calibrated on and reusing `lattice/cdt_finder/`
(`consolidation/CDT_FINDER.md`).  Literature read in source: CDT,
**The unbounded denominators conjecture**, arXiv:2109.09040 §2; CDT,
**The linear independence of 1, ζ(2), L(2,χ₋₃)**, arXiv:2408.15403 §§2,5,6,7,8.
Tags: [proved] [verified] [open].  **No irrationality claim is made anywhere.***

---

## 0. Verdict

| claim | verdict |
|---|---|
| The $p$-adic gain can be inserted into CDT's proof; it is a **one-line change to the Liouville step**, not a change to the analytic machinery | **[proved]**, §2 |
| The gain is $\gamma_p=\log p\bigl[\varsigma_{\min}(1-\tfrac1m)+\tfrac2m\int_0^1\delta^{\downarrow}(t)\max(mt-1,0)\,dt\bigr]$ — a rearrangement functional of the $p$-adic slopes, exactly mirroring CDT's $\tau^\flat$ with the order reversed | **[proved]**, §2.4 |
| The resulting bound is **exactly invariant** under rescaling the coordinate by a rational number (product formula), verified to $10^{-15}$ | **[proved]** + **[verified]** `adelic_bound.py` |
| $\gamma=0$ reproduces CDT's own numbers to the last digit ($\tau=16603/3920$, bound $13.9938$, margin $+0.0053$) | **[verified]** |
| **Level-8 Catalan host: the adelic gain flips the entry test.** Archimedean-only entry $-0.0766$ at the uniformisation ceiling; with the (unconditional) $2$-adic slope $2$ of the pure polylogarithm module, entry $=+0.1780$ | **[verified]**, §4.2 |
| …but the **margin stays at $-7.97$** (ceiling) / $-14.45$ (CDT's realised contour), and with CDT's realised contour the entry test still fails ($-0.2852$) | **[verified]** |
| **The doubly-small module on the level-8 host is empty.** The target-zero class $(1-4V_2)E$ has $\xi_\infty=\xi_2=0$, is $2$-adically overconvergent of slope $3$, but its companion has **non-trivial monodromy at the fold** $t=1/8$: $B_0(t)/A(t)\sim0.207\log n$, so it is not single-valued on the CDT domain | **[verified]**, §4.3 — *this is the load-bearing negative* |
| **Apéry-perfect hosts have zero adelic gain at every prime.** $\sigma_p=v_p(c)$ and $c=\pm1$; and $\lambda_2$ a unit kills the pure module's gain too | **[proved]** (Cor. 3.2 of `CRYSTAL_THEOREM_F.md`) + **[verified]** |
| Hence the $X_1(5)\,\mathrm{Sym}^2$ host ($L(3,\chi_5)$) gets **no help at all**: its deficit stays $-0.435$ (ceiling) / $-0.898$ (contour) | **[verified]** |
| **Structural ceiling.** $\log(256/\lambda_2)+\sigma_p\log p\le\log(256\lambda_1)$, with equality iff $c$ is a power of $p$ | **[proved]**, §5 |

One sentence: *the adelic bound is a real theorem and it is worth about $+0.25$ to $+0.9$ of entry budget on the Catalan host — enough to pass CDT's entry test at the uniformisation ceiling for the first time, and not remotely enough to close the $-8$ margin; and the one architectural device that would close it, an unconditional doubly-small orbit, does not exist on that host.*

---

## 1. What had to be corrected first

The task, and the GPT-era `book/v8/07a_adelic_layerC.tex`, both assume CDT's
proof runs through a **Gelfond determinant** $\Delta$, to be bounded above by
Cauchy estimates and below by $1/\mathrm{den}(\Delta)$.  **There is no
determinant in either CDT paper**, and CDT say explicitly why they avoid one
(§2.2 of 2408.15403: Cramer's rule gives a bound exponential in $M\alpha\asymp\alpha^2$
instead of $\alpha$).  The architecture is:

1. a **Thue–Siegel box principle** producing a nonzero auxiliary
   $F(\mathbf x)=\sum_{\mathbf{i,k}}c_{\mathbf{i,k}}\mathbf x^{\mathbf k}\prod_s f_{i_s}(x_s)$,
   $c\in\mathbf Z$, vanishing to order $\alpha=mdD/2$ at $\mathbf 0$;
2. **one rational number** $\beta$ — the coefficient of a lowest-order monomial
   $\beta\mathbf x^{\mathbf n}$;
3. a **Liouville lower bound** $|\beta|\ge1/\mathrm{den}(\beta)\ge e^{-\alpha\tau(\mathbf b;\mathbf e)}$
   against a **Cauchy upper bound**
   $\log|\beta|\le\frac\alpha m\!\int_0^1\!2t\,g^*+\ \cdots-\alpha\log|\varphi'(0)|$.

The good news is that this makes the $p$-adic insertion *cleaner* than the
determinant version, because $\beta$ is a single rational number:

> **if a nonzero rational $\beta$ has $\mathrm{den}(\beta)\le e^{\tau\alpha}$ and
> $v_p(\beta)\ge\gamma\alpha/\log p$, then $|\beta|\ge e^{(\gamma-\tau)\alpha}$.**

(Write $\beta=N/D$ in lowest terms; $v_p(D)\ge0$ forces $v_p(N)\ge v_p(\beta)$,
so $|N|\ge p^{v_p(\beta)}$.)  No product formula gymnastics are needed; the
statement *is* the product formula for the two-element set $\{p,\infty\}$.

`book/v8`'s Theorem "Frobenius-filtered holonomy" is, read carefully, **not**
about $p$-adic overconvergence at all: its $\tau_F=\sum_jb_j(1-u_j^2/m^2)$ is
verbatim CDT's own $\tau^\flat(\mathbf b)$ (6.0.4), i.e. the *denominator*
functional.  It is a correct citation with a misleading name; it contains no
$p$-adic input, and the "$\kappa_p$" discussion around it is decoration.  The
$p$-adic input is the new $\gamma_p$ below.

---

## 2. The theorem

### 2.1 Hypotheses

Fix $m,r\ge1$, an integration vector $\mathbf e\in\mathbf N^m$, and an $m\times r$
array $\mathbf b=(b_{i,j})$ of nonnegative reals whose columns have CDT's **step
shape** $0=b_{1,j}=\dots=b_{u_j,j}<b_{u_j+1,j}=\dots=b_{m,j}=:b_j$.  Put
$\sigma_i=\sum_jb_{i,j}$ (nondecreasing in $i$) and

$$\tau^\flat(\mathbf b)=\frac1{m^2}\sum_{i=1}^m(2i-1)\sigma_i=\sigma_m-\frac1{m^2}\sum_ju_j^2b_j,
\qquad \tau^\sharp(\mathbf e)\ \text{as in CDT Def. 6.0.1/6.0.5}.$$

Let $S$ be a finite set of primes.  For each $p\in S$ let
$\varsigma^{(p)}=(\varsigma^{(p)}_1,\dots,\varsigma^{(p)}_m)\in\mathbf R^m$ be a
vector of **$p$-adic slopes**.

Suppose $f_1,\dots,f_m\in\mathbf Q[\![x]\!]$ are $\mathbf Q(x)$-linearly
independent and holonomic, of CDT denominator type
$$f_i(x)=a_{i,0}+\sum_{n\ge1}a_{i,n}\frac{x^n}{n^{e_i}\,[1,\dots,b_{i,1}n]\cdots[1,\dots,b_{i,r}n]},
\qquad a_{i,n}\in\mathbf Z,$$
and in addition satisfy, for every $p\in S$, writing $c_{i,n}$ for the actual
coefficient of $x^n$,

$$(\ast_p)\qquad v_p(c_{i,n})\ \ge\ \varsigma^{(p)}_i\,n-o(n)\qquad(n\to\infty),$$

equivalently: $f_i$ converges and is bounded on the $p$-adic disc
$|x|_p\le p^{\varsigma_i^{(p)}}$ (an *overconvergence* hypothesis when
$\varsigma_i^{(p)}>0$, a *geometric-denominator* hypothesis when
$\varsigma_i^{(p)}<0$).

Let $\varphi:(\overline{\mathbf D},0)\to(\mathbf C,0)$ be holomorphic with every
$\varphi^*f_i$ meromorphic on $|z|<1$.

### 2.2 Statement

> **Theorem (adelic holonomy bound).**  With
> $$\boxed{\ \tau_{\mathrm{ad}}\ :=\ \tau^\flat(\mathbf b)+\tau^\sharp(\mathbf e)\ -\ \sum_{p\in S}\gamma_p\ },\qquad
> \gamma_p:=\log p\left[\varsigma_{\min}^{(p)}\Bigl(1-\frac1m\Bigr)+\frac2m\int_0^1\delta^{\downarrow}_p(t)\max(mt-1,0)\,dt\right],$$
> where $\varsigma^{(p)}_{\min}=\min_i\varsigma^{(p)}_i$, $\delta_{p,i}=\varsigma^{(p)}_i-\varsigma^{(p)}_{\min}$
> and $\delta_p^{\downarrow}$ is the **non-increasing** rearrangement of
> $(\delta_{p,i})$ read as a step function on $[0,1]$ with steps of width $1/m$:
> if $\log|\varphi'(0)|>\tau_{\mathrm{ad}}$ then
> $$\boxed{\ m\ \le\ \frac{\displaystyle\iint_{\mathbf T^2}\log|\varphi(z)-\varphi(w)|\,\mu(z)\mu(w)}{\log|\varphi'(0)|-\tau_{\mathrm{ad}}}\ }$$
> and likewise for all of CDT's refined numerators (the rearrangement integral
> 6.0.15, the multi-$\varphi$ form 6.0.2, the convexity forms of §7).  If the
> $f_i$ are not assumed holonomic a priori, the entry condition must be
> strengthened to $\log|\varphi'(0)|>\max\{\sigma_m-\sum_p\gamma_p,\ \tau_{\mathrm{ad}}\}$.

Two special cases:
* **uniform profile** ($\varsigma_i\equiv\varsigma$): $\gamma_p=\varsigma\log p\,(1-1/m)$;
* **step profile** ($u$ functions of slope $\varsigma$, the rest $0$):
  $\gamma_p=\varsigma\log p\,(u-1)^2/m^2$.

### 2.3 Proof

Run CDT's proof of Theorem 6.0.2 (arXiv:2408.15403 §6) verbatim.  Every step —
the Siegel lemma §6.2, the statistics of the lowest-order exponent
(Cor. `F stats`), the Vandermonde damping, the numerical integration, the
archimedean Cauchy bound §6.3, the denominator computation §6.4 — is untouched:
none of them sees $(\ast_p)$.  The single change is in the **arithmetic bound on
$\beta$**, CDT's §6.4 endgame.

CDT prove $\mathrm{den}(\beta)\le\exp(\alpha\tau^\flat(\mathbf b)+\alpha\tau^\sharp(\mathbf e)+o(\alpha))$.
We add: for each $p\in S$,
$$v_p(\beta)\ \ge\ \frac{\gamma_p}{\log p}\,\alpha-o(\alpha). \tag{$\dagger$}$$
Granting $(\dagger)$, the Liouville step becomes
$$\log|\beta|\ \ge\ \sum_{p\in S}\gamma_p\,\alpha-\alpha\bigl(\tau^\flat+\tau^\sharp\bigr)+o(\alpha)
= -\alpha\,\tau_{\mathrm{ad}}+o(\alpha),$$
by the elementary remark of §1 applied one prime at a time (the primes of $S$
are distinct, so the $p$-parts of the numerator multiply).  Comparing with
CDT's unchanged analytic upper bound and letting $\alpha\to\infty$, $d\to\infty$,
$\epsilon\to0$ gives the theorem.

*Proof of $(\dagger)$.*  Write $\beta=\sum_{\mathbf{i,k}}c_{\mathbf{i,k}}\prod_{s=1}^d c_{i_s,\,n_s-k_s}$,
a $\mathbf Z$-combination of $\exp(o(\alpha))$ terms; $v_p(\beta)$ is at least the
minimum of $\sum_sv_p(c_{i_s,n_s-k_s})$ over the terms.  By $(\ast_p)$ each term
has
$$\sum_s v_p\bigl(c_{i_s,n_s-k_s}\bigr)\ \ge\ \sum_s\varsigma_{i_s}(n_s-k_s)-o(\alpha),
\qquad n_s\ge k_s\ \ (\text{else the term is absent}).$$
Split $\varsigma_i=\varsigma_{\min}+\delta_i$.

*(a) The uniform part.*  $\sum_s\varsigma_{\min}(n_s-k_s)=\varsigma_{\min}\bigl(|\mathbf n|-|\mathbf k|\bigr)$.
By Cor. `F stats`, $|\mathbf n|=\alpha(1+O(\epsilon))$ and $|\mathbf k|\le\alpha/m\,(1+O(\epsilon))$
(the $\mathbf k$ live in $[0,D]^d$ with $\alpha=mdD/2$), so this contributes
$\varsigma_{\min}\alpha(1-1/m)+O(\epsilon\alpha)$.  *(This is exactly the $p$-adic
Gauss-norm/maximum-modulus estimate $|\beta|_p\le\|F\|_{R}R^{-|\mathbf n|}$ at
$R=p^{\varsigma_{\min}}$, and it is what makes the whole bound scale-covariant.)*

*(b) The residual part.*  The species index $\mathbf i$ is constrained to
$V_m^d$ (each species occurs exactly $d/m$ times: CDT's §6.1.3 "arithmetic
sharpening", the same hypothesis that produces $\tau^\flat$).  Hence
$\sum_s\delta_{i_s}(n_s-k_s)$ is minimised by giving the largest $\delta$ to the
slots with the smallest $n_s-k_s$.  We know from `F stats` that, sorted
increasingly, $n_{\pi(j)}=mD\,(j/d)+O(\epsilon D)$, and $k_s\le D$ always; so
$n_s-k_s\ge\max\{mD(j/d)-D,0\}$ for the slot of rank $j$.  Therefore
$$\min_{\mathbf{i}\in V_m^d}\sum_s\delta_{i_s}(n_s-k_s)\ \ge\
dD\int_0^1\delta^{\downarrow}(t)\max(mt-1,0)\,dt-O(\epsilon\alpha)
=\frac{2\alpha}m\int_0^1\delta^{\downarrow}(t)\max(mt-1,0)\,dt-O(\epsilon\alpha),$$
using $dD=2\alpha/m$.  Adding (a) and (b) and dividing by $\log p$ gives $(\dagger)$. $\square$

**Where $(\ast_p)$ actually bites.**  The hypothesis is on the *coefficients*, not
on $f_i$ as a function on a $p$-adic domain — but the two are the same statement,
because the LCM denominators are harmless at a **fixed** prime:
$v_p\bigl([1,\dots,bn]\bigr)=\lfloor\log(bn)/\log p\rfloor=O(\log n)=o(n)$.
*This is the clean answer to "the entries have both a denominator part and an
overconvergence part": at a fixed $p$ the LCM part contributes $o(n)$ and drops
out, so $\varsigma^{(p)}_i$ may be read off the integral numerators $a_{i,n}$ or
off the true coefficients indifferently.*  What does **not** separate is a
**geometric** denominator $A^{n}$ with $p\mid A$: that is a genuine negative
slope $-v_p(A)$ and the theorem accounts for it automatically (and correctly:
setting all $\varsigma_i=-v$ gives $\gamma_p=-v\log p\,(1-1/m)$, i.e. $\tau$
increases, which is the right cost).

### 2.4 Consistency checks (all machine-verified, `adelic_bound.py`)

1. **$\gamma=0$ is CDT.**  $m=14$, $\mathbf b$ with $u=(1,3)$, $b=(2,2)$,
   $\mathbf e$ CDT's: $\tau^\flat=191/49$, $\tau^\sharp=27/80$,
   $\tau=16603/3920=4.235459$, $|\varphi'(0)|=161.0812$, bound $13.9938$,
   margin $+0.0053$.  Identical to `CDT_FINDER.md` §2.
2. **Scale covariance = the product formula.**  Replacing the coordinate
   $x\rightsquigarrow x/p^{v}$ multiplies $|\varphi'(0)|$ and the
   Bost–Charles integral by $p^{-v}$ and shifts every $\varsigma_i$ by $+v$.
   The signed margin $m(\log|\varphi'(0)|-\tau_{\mathrm{ad}})-\mathrm{BC}$ is
   **exactly** unchanged: verified to $1.8\cdot10^{-15}$ at $v=3$, $m=21$.
   Note this works *only* because part (a) of the proof carries the weight
   $(1-1/m)$, the same weight as the archimedean $\log|\varphi'(0)|-\mathrm{BC}/m$;
   a term-by-term-only bound (weight $(1-1/m)^2$) would break it.  The check is
   therefore a genuine test of the exponents, not a tautology.
3. **No false contradictions.**  Applied to the unconditional pure
   polylogarithm module alone on the level-8 host (7 functions, slope $2$,
   type $[1..2n]^2$, $u=(1,3)$), the bound gives $m\le9.46$: consistent, with
   room.  (It becomes violated at $m=10$, i.e. the theorem now *asserts* that
   that module has dimension $\le9$.)

### 2.5 What is lossy

Part (b) uses only $k_s\le D$; a sharper joint law for $(n_s,k_s)$ would improve
the step-profile weight from $(u-1)^2/m^2$ towards $(u/m)^2$.  The improvement
is at most $O(1/m)$ and is not what decides anything below.  Part (a) is sharp.
CDT's own Remark `genbcomplicated` (the column-shape hypothesis on $\mathbf b$
is *necessary* for their Lemma "single step valuationwise") has no analogue
here: the $\gamma$-side optimisation is a plain rearrangement inequality with no
shape hypothesis, because it is a *minimum* over a product of a monotone
sequence with an arbitrary vector.

---

## 3. Mixed profiles: what the Catalan hypothesis actually produces

Under the hypothesis $G=\mathrm{Catalan}\in\mathbf Q$, on a host $(\Gamma,t,F)$
carrying a $\chi_{-4}$ Eisenstein class $\Phi$ with companion $B$, the three
kinds of function have **genuinely different local profiles**:

| orbit | archimedean | $2$-adic | rational coefficients? |
|---|---|---|---|
| pure $\mathrm{Li}_j(x/s)$, $s=1/\lambda_2$ | on $\mathbf P^1\setminus\{0,s,\infty\}$ | slope $v_2(\lambda_2)$ in $x$ | yes, unconditionally |
| conditional $H=bB-aA$, $a+b\xi_\infty=0$ | **trivial monodromy at the fold** $t_1$ — this is what the hypothesis buys | **slope $0$**: $2$-adically $H\approx(b\xi_2-a)A$ and $\xi_2\ne\xi_\infty$ *as numbers* | yes, conditionally |
| doubly small $B_i-rB_j$, $r\in\mathbf Q$, $r=\xi_\infty^{(i)}/\xi_\infty^{(j)}$ | trivial monodromy at $t_1$ (both conditional parts are regular there) | slope $>0$ (both companions converge to the *same* $\xi_2$ up to the same rational $r$, by Theorem F / Conjecture D) | yes, **unconditionally** |

The middle row is the subtlety the task names, and it is real and unavoidable:
$\xi_2=\tfrac12\zeta_2(2)$ is irrational (Calegari 2005), so $b\xi_2-a\ne0$ and
the conditional function is $2$-adically of slope exactly $0$.  **The hypothesis
$G\in\mathbf Q$ buys archimedean smallness and buys nothing $2$-adically.**

The theorem of §2 is stated for exactly this situation: each $f_i$ carries its
own $\varsigma_i^{(p)}$, and they enter through a rearrangement functional, just
as CDT's per-function denominator types $b_{i,j}$ enter through $\tau^\flat$.
The two rearrangements run in **opposite directions** — denominators are
maximised (weight $2t$ against the increasing rearrangement), $p$-adic slopes
are minimised (weight $\max(mt-1,0)$ against the decreasing one) — and they are
optimised **independently**, which is legitimate because each is a valid bound
over all species assignments separately.

**The coordinate.**  CDT's normaliser descent $w(x)=sx/(x-s)$,
$y=x+w(x)=x^2/(x-s)$ changes the slopes: for $|x|_p<|s|_p$ one has
$|y|_p=|x|_p^2/|s|_p$, so
$$\varsigma_y\ \ge\ \max\{0,\ 2\varsigma_x-2v_p(\lambda_2)\},$$
the $0$ coming from the integrality already assumed in the denominator type.
All numbers below are in the $y$-coordinate, which is the one in which the
finder's $|\varphi'(0)|\le256\,s$ and $\mathrm{BC}=11.845+\log s$ are stated.

---

## 4. Catalan, level 8 (Zagier E)

### 4.1 The host and its classes (`rowE_sources.gp`, exact)

$(a,b,c)=(12,4,32)$; singular $t$ are the roots of $1-12t+32t^2=(1-4t)(1-8t)$,
so $\lambda_1=8$ (fold, $t_1=1/8$), $\lambda_2=4$ ($t_2=1/4$), $s=1/\lambda_2=1/4$;
$k=2$ sharp ($d_n^2b_n\in\mathbf Z$, $d_n b_n\notin\mathbf Z$, checked $n\le31$).
Source $\Phi_{\mathbf E}=(1-8V_2)E$ with $E=E_3^{\chi_{-4},\mathbf 1}$;
$\xi_\infty=G/2$, $\xi_2=\tfrac12\zeta_2(2)$.

The $\chi_{-4}$ **inner** Eisenstein space on the host is $2$-dimensional,
$\langle E,V_2E\rangle$.  For $\Phi=P(V)E$, Theorem F and Theorem B$^*$ give
$$\xi_\infty=-\tfrac12P(2)\,G,\qquad \xi_2=-\tfrac12P(2)\,\zeta_2(2)\quad(\mathcal E_2\equiv1\ \text{since}\ \chi_{-4}(2)=0),$$
so *every* class has $\xi_\infty:\xi_2=G:\zeta_2(2)$ — rigidity — and there is a
unique **target-zero** line $P(2)=0$, i.e. $\Phi_0=(1-4V_2)E$.  The outer
orientation $E_3^{\mathbf 1,\chi_{-4}}$ carries $\zeta(2)$ archimedean and
$\xi_2=0$ (the "$\pi^2$ Tate defect" of the archive) and is not a Catalan
direction.

I rebuilt $t(q)=q-4q^2+12q^3-32q^4+78q^5-\cdots$ and
$F(q)=1+4q+4q^2+4q^4+8q^5+\cdots$ from $\Phi_{\mathbf E}=F\,\mathcal D_qt$ by
series reversion, and verified $B_{\Phi_{\mathbf E}}(t)=\sum b_nt^n$ **equals the
census companion exactly** ($b=0,1,7,404/9,2603/9,\dots$, 8 coefficients checked)
— so the companion construction $B_\Psi=$ ($t$-expansion of $F\cdot\mathcal D^{-2}\Psi$)
is correct.

### 4.2 The gain that is real

| profile | $2$-adic slope in $x$ | in $y$ | evidence |
|---|---|---|---|
| pure $\mathrm{Li}_j(4x)$ orbit | $2=v_2(\lambda_2)$ | $2$ | `pure_2adic.py`: $v_2$ of the $y$-coefficients of $\mathrm{Sym}^{\pm}\mathrm{Li}_j(4x)$, $j\le4$, $k\le12$, is $2k-O(\log k)$ |
| conditional orbit | $0$ | $0$ | $v_2(a_n^{\mathbf E})=2s_2(n)$, $d_n^2b_n\in\mathbf Z$, $\xi_2\notin\mathbf Q$ |
| $B_{\mathbf E}-\xi_2A$ (irrational coefficients) | $5$ | $8$ | `hosts_padic.gp`: $v_2(\xi_{n+1}-\xi_n)/n=4.70,4.80,4.858$ at $n=60,110,155\to5=v_2(32)$ |

With CDT's architecture ($m=14=7$ pure $+\,7$ conditional), only the pure orbit
contributes, giving $\gamma_2=+0.2546$:

| accounting | $\tau$ | $\gamma_2$ | **entry** | **margin** |
|---|---|---|---|---|
| archimedean only, ceiling $\log(256s)=4.158883$ | $4.235459$ | $0$ | $\mathbf{-0.0766}$ | $-11.531$ |
| archimedean only, CDT's realised contour $3.695614$ | $4.235459$ | $0$ | $-0.5398$ | $-18.017$ |
| **adelic, ceiling** | $4.235459$ | $+0.2546$ | $\mathbf{+0.1780}$ | $-7.966$ |
| adelic, CDT's realised contour | $4.235459$ | $+0.2546$ | $-0.2852$ | $-14.452$ |

$\mathrm{BC}=11.845+\log s=10.458706$ throughout.  **This is the first time the
Catalan host passes CDT's entry condition.**  It is also the whole of the good
news: the margin is still $-7.97$, i.e. one needs roughly $14+8/0.18\approx58$
functions, or a much better pure inventory, to reach $0$.

For reference, with the "best conceivable" inventory ($u_1=u_2=m/2$, i.e. half
the functions free of both LCM layers — nothing on this host is known to
realise it) the adelic margin is $+4.61$ at $m=14$ against $+1.04$ archimedean-only,
and $+13.56$ vs $+7.94$ at $m=20$.  The adelic gain is worth about $3.5$ units of
margin at $m=14$ and $5.6$ at $m=20$ — real, but an order of magnitude smaller
than the inventory lever.

### 4.3 The load-bearing negative: no doubly-small orbit on this host

The doubly-small construction of §3 requires **two** classes whose companions
are *each* regular at the fold.  On a host with a fixed $(t,F)$ there is only
one canonical source ($\Phi=F\mathcal D_qt$), and the measurement is decisive.
Writing $B_c:=B_{(1-cV_2)E}$ (so $B_8=B_{\mathbf E}$, $B_4=B_0$ is the
target-zero class):

| class | $P(2)$ | $\xi_\infty$ predicted | measured $b_n/a_n$ at $n=10,20,30$ |
|---|---|---|---|
| $B_8$ (canonical) | $-1$ | $G/2=0.4579828$ | $0.4574377,\ 0.4579822,\ 0.45798280$ — converges ✅ |
| $B_4$ (target zero) | $0$ | $0$ | $0.8227,\ 0.9680,\ 1.0498$ — **diverges like $0.207\log n$** |
| $B_2$ | $1/2$ | $G/4=0.2289914$ | $1.0053,\ 1.2230,\ 1.3456$ — **diverges like $0.311\log n$** |

The divergence rate is exactly affine in $c$ (slope $0.0518\,(8-c)$, fitted to
three digits at two independent pairs of $n$), so it is the generic behaviour and
the canonical class is the unique zero of it.  Interpretation: $H_\xi=F(\Theta-\xi)$
is a weight-zero modular function on $X_0(8)$; the constant $\xi$ can only fix
its behaviour at the cusp $t=0$, and regularity at the *second* cusp $t=1/8$ is a
property of $\Phi$ alone.  Only the canonical $\Phi$ has it.  Consequently
$B_0$ has **non-trivial local monodromy at $t_1=1/8$**, is not single-valued on
CDT's domain $\mathbf P^1\setminus\{0,4s,\infty\}$ (whose universal cover is the
$\lambda$-uniformisation used for $\varphi$), and **cannot be used**, however
good its $2$-adic behaviour.  (For the record, that behaviour is good:
$v_2$ of the $t$-coefficients of $B_0$ is $3n-O(\log n)$ — measured
$v_2(b^0_n)/n=2.90,\,2.96,\,3.06$ at $n=30,27,31$ — i.e. slope $3$ in $x$,
$4$ in $y$; it is genuinely $2$-adically overconvergent, just archimedean-ly
multivalued.)

Also relevant, and independent of the above: even if such functions existed, the
theorem of §2 **caps how many there can be**, because they are unconditional.
On this host, with no denominator-free members, the cap is
$$\#\{\text{doubly small}\}\ \le\ 7\ (\varsigma_y=2),\quad 5\ (3),\quad 4\ (4),\quad 3\ (5),\quad 2\ (8).$$
So the $d=7$ rows of the hypothetical table (the only ones with positive margin,
$+2.58$ resp. $+7.33$) are **self-contradictory**: they assert more doubly-small
functions than the same theorem permits.  Honest hypotheticals are $d\le4$
($\varsigma_y=4$: margin $-2.86$) and $d\le2$ ($\varsigma_y=8$: margin $-5.61$).

Where to look instead: a host on which **two** distinct $\chi_{-4}$ classes have
fold-regular companions.  That is exactly a host with two Apéry rows sharing
$(t,F)$ — equivalently a second canonical source — or a level-16/level-144 host
where the "period-killed" classes are engineered to vanish at the fold cusp.
This is a well-posed, purely archimedean, finite computation on
$M_3(\Gamma_0(16),\chi_{-4})$ and $M_3(\Gamma_0(144))$, and it is **the** thing to
compute next.  It was not done here.

---

## 5. Where the gain comes from, and its hard ceiling

For an (R2)/(R3) row the singular $t$ are the roots of $P_r$, of product $1/c$,
and $\sigma_p=v_p(c)$ (Cor. 3.2 of `CRYSTAL_THEOREM_F.md`, proved from the exact
Casoratian).  The pure module $\mathrm{Li}_j(x/s)=\sum(\lambda_2x)^n/n^j$ has
$p$-adic slope $v_p(\lambda_2)$.  Since $\lambda_1\lambda_2=c$,

$$\boxed{\ \underbrace{\log\frac{256}{\lambda_2}}_{\text{archimedean ceiling}}\ +\ \underbrace{\sigma_p\log p}_{\text{$p$-adic ceiling}}\ \le\ \log\bigl(256\,\lambda_1\bigr)\ }$$

with equality iff $c$ is a power of $p$.  Verified:

| host | $\lambda_1$ | $\lambda_2$ | $c$ | ceiling | max gain | total | $\log(256\lambda_1)$ |
|---|---|---|---|---|---|---|---|
| Zagier E (Catalan) | $8$ | $4$ | $32$ | $4.1589$ | $3.4657$ | $7.6246$ | $7.6246$ |
| Zagier C (CDT) | $9$ | $1$ | $9$ | $5.5452$ | $2.1972$ | $7.7424$ | $7.7424$ |
| Apéry $\zeta(3)$ | $1$ | $1$ | $1$ | $5.5452$ | $0$ | $5.5452$ | $5.5452$ |
| $X_1(5)\mathrm{Sym}^2$ | $\varphi^5$ | unit | $-1$ | $5.5452$ | $0$ | $5.5452$ | $7.9512$ |
| Domb | $16$ | $4$ | $64$ | $4.1589$ | $4.1589$ | $8.3178$ | $8.3178$ |

**The archimedean and $p$-adic budgets are in direct competition.**  The
Apéry-perfect / Fricke-palindromic hosts, which are exactly the ones that
maximise the archimedean ceiling ($\lambda_2=1$ or a unit), have $c=\pm1$ and
therefore **zero $p$-adic slope at every prime, on every row and on the pure
module** — verified directly on Zagier D $(11,3,-1)$, the weight-one layer of the
$X_1(5)$ host: $v_p(b_{n+1}/a_{n+1}-b_n/a_n)\le0$ at $p=2,3,5,11$, $n=150$.
Hence:

> **$X_1(5)\ \mathrm{Sym}^2$ ($1,\zeta(3),L(3,\chi_5)$ over $\mathbf Q(\sqrt5)$) gets
> no adelic help whatever.  Its entry deficit stays $-0.435$ at the ceiling and
> $-0.898$ with CDT's contour, margin $-24.4$; $\gamma_p\equiv0$.**

The coordinator's suggestion that the $\chi_5$ classes on that host carry a
$5$-adic slope $3$ conflates it with the $\eta$ row $(11,5,125)$, which lives on
**level 20**, has $\lambda_2\approx11.18$ (census) and $c=125$, hence
$\sigma_5=v_5(125)=3$ and a gain of $3\log5=4.83$ — bought at the price of an
archimedean ceiling $\log(256/11.18)=3.13$; the total is $\log(256\lambda_1)$ again.  The two hosts cannot be combined:
the bound needs one common $\varphi$ and one common coordinate.

---

## 6. What remains

1. **The theorem's remaining softness.**  (i) The $O(1/m)$ loss in part (b) of
   the proof; harmless.  (ii) The $\sigma_m$-vs-$\tau$ entry threshold in the
   *non*-holonomic case: CDT's Remark `converse remark` constructs a continuum
   of counterexamples when $\log|\varphi'(0)|\le\sigma_m$, and the $p$-adic
   analogue of that construction (choices constrained mod $p^{\varsigma n}$) has
   not been worked out.  All our $f_i$ are holonomic (Picard–Fuchs), so the
   relaxed entry condition applies and this is not load-bearing.  **[open, minor]**
2. **Number fields.**  CDT's Remark `BCboundK` (§8.5) gives the several-real-places
   version; the $p$-adic gain fits it with $\gamma_p$ multiplied by $[K:\mathbf Q]$
   in the denominator, matching `CDT_FINDER.md` §3's normalisation.  Not written
   out.  **[open]**
3. **The Catalan search.**  Find a host carrying **two** $\chi_{-4}$ Eisenstein
   classes with fold-regular companions (level 16, level 144, or an Atkin–Lehner
   quotient with $\deg t>1$).  This is the only route to a doubly-small orbit
   and it is a finite exact computation.  **[open, the next task]**
4. **Functional independence.**  Every margin above presumes $\mathbf Q(y)$-linear
   independence of the $m$ functions (CDT's Lemma 12.1.1 analogue).  Nothing here
   proves it; the finder's rank check exists only for CDT's own $14$.  **[open]**
5. **Continuation proofs.**  The transported contour loss $0.62922$ and the
   Bost–Charles shape term $6.763$ are *estimates* carried over from CDT's
   level-6 contour (`CDT_FINDER.md` §8, estimate 1).  The inner singularity
   $y(t_1)$ sits at a different hyperbolic depth on the level-8 host, and the
   entry test at the ceiling now passes by only $+0.178$ — so the contour is no
   longer a second-order question, it is the question.  **[open]**

**No irrationality claim.**  Nothing above proves anything about Catalan's
constant, $\zeta(3)$, or $L(3,\chi_5)$.  The one new *unconditional* assertion is
the dimension cap of §4.3 (at most $4$ functions on the level-8 host that are
simultaneously fold-regular, of type $[1..2n]^2$, and $2$-adically overconvergent
of slope $4$), and it is a consequence of the theorem of §2, whose proof is a
modification of CDT's that I believe to be correct but which has had no
independent check.

---

## 7. Correction (2026-08-22): the conditional orbit's $y$-slope is $-2$, not $0$

*Added by Fable (Opus 5) after `consolidation/CATALAN_THREE_PERIOD.md` §5.0 found
the same defect on the level-16 host.  Scripts: `lattice/catalan_three_period/`
(`20_level8_slopes.gp`, `21_level8_margin.py`; logs `20run.log`, `21run.log`).
**The tables of §4.2 are superseded, not deleted**: the old rows are reproduced
below marked SUPERSEDED so that the change can be audited.*

### 7.1 What was wrong

§3 asserts, for the normaliser descent $w(x)=sx/(x-s)$, $y=x+w(x)=x^2/(x-s)$,
$$\varsigma_y\ \ge\ \max\{0,\ 2\varsigma_x-2v_p(\lambda_2)\},$$
"the $0$ coming from the integrality already assumed in the denominator type".
**That floor is not available here.**  $x$ and $w(x)$ are the two roots of
$X^2-yX+sy=0$, i.e. $x=\tfrac12\bigl(y-\sqrt{y^2-4sy}\bigr)$; with $s=\tfrac14$
this is $x=\tfrac12\bigl(y-\sqrt{y^2-y}\bigr)$, whose coefficients carry $2$-power
**denominators**.  So $\mathrm{Sym}^{\pm}$ of an integral $x$-series need not be
integral in $y$, and the CDT denominator hypothesis
($a_{i,n}\in\mathbf Z$ over $n^{e}[1..b_1n][1..b_2n]$) is *not* satisfied by the
functions in this $y$ — the extra $2$-power denominator is a genuine geometric
denominator, i.e. a genuine negative slope, which the theorem of §2 accounts for
correctly and at full weight $(1-1/m)$ through $\varsigma_{\min}$.

### 7.2 Measured slopes on the level-8 host (`20run.log`)

$t$-series to $n=300$; $y$-expansions extracted exactly over $\mathbf Q$ by the
`to_y` algorithm of `lattice/cdt_finder/indep_check2.py` (no $\sqrt y$ needed).
$\sigma(t)=t/(4t-1)$, $y=4t^2/(4t-1)$.

| object | $v_2(d_k)/k$ at $k=30,60,147$ | $\varsigma_y$ | assumed in §4.2 |
|---|---|---|---|
| $\mathrm{Sym}^+\mathrm{Li}_2(4t)$ | $1.800,\ 1.867,\ 1.973$ | $+2$ | $2$ ✅ |
| $\mathrm{Sym}^+\mathrm{Li}_3(4t)$ | $1.633,\ 1.767,\ 1.918$ | $+2$ | $2$ ✅ |
| $\mathrm{Sym}^-\mathrm{Li}_2(4t)/(t-\sigma t)$ | $1.767,\ 1.867,\ 1.932$ | $+2$ | $2$ ✅ |
| $\mathrm{Sym}^-\mathrm{Li}_3(4t)/(t-\sigma t)$ | $1.600,\ 1.767,\ 1.878$ | $+2$ | $2$ ✅ |
| $\mathrm{Sym}^+A$ (the host) | $-1.700,\ -1.850,\ -1.939$ | $\mathbf{-2}$ | — |
| $\mathrm{Sym}^+B_{\mathbf E}$ (raw companion) | $-1.767,\ -1.883,\ -1.952$ | $\mathbf{-2}$ | — |
| $\mathrm{Sym}^+(B_{\mathbf E}+\mu A)$, $\mu=\tfrac73,\tfrac18,-\tfrac12$ (**conditional**) | $-1.767,\ -1.883,\ -1.952$ | $\mathbf{-2}$ | $\mathbf 0$ ❌ |
| $\mathrm{Sym}^-(B_{\mathbf E}+\mu A)/(t-\sigma t)$ | $-1.733,\ -1.867,\ -1.932$ | $\mathbf{-2}$ | — |
| $\mathrm{Sym}^+B_0$, $\Phi_0=(1-4V_2)E$ | $2.900,\ 2.950,\ 2.993$ | $+3$ | $4$ ❌ |

**[verified]**  $\mathrm{Sym}^+1$ and $\mathrm{Sym}^+\mathrm{Li}_1(4t)$ are
identically zero beyond low order — $(1-4t)(1-4\sigma t)=1$, so the $\log$ cancels.

The **$x$-slopes of §4.2 were all correct** (measured in $t$: $A$ $0.027$,
$B_{\mathbf E}$ $0.020$, $\mathrm{Li}_2(4t)$ $1.987$, $B_0$ $2.990$ at $n=300$,
i.e. $0,0,2,3$).  The error is entirely in the descent: the correct rule on this
host is $\varsigma_y=2\varsigma_x-2$ (matching `CATALAN_TWO_CLASSES.md` §3's
$\varsigma_y=2(\varsigma_x-1)$ for $v_2(s)=-2$), with **no floor at $0$** — which
gives $2\cdot0-2=-2$ for the conditional orbit.  ($B_0$ measures $+3$ where the
rule would give $+4$; $B_0$ is unusable anyway, §4.3.)

### 7.3 The corrected §4.2 margin table (`21run.log`)

$m=14=7$ pure $+7$ conditional, $u=(1,3)$, $b=(2,2)$, CDT's $\mathbf e$,
$\tau=4.235459$, $\mathrm{BC}=10.458706$, ceiling $4.158883$, realised $3.695614$.

| accounting | $\gamma_2$ | **entry** | **margin** |
|---|---|---|---|
| archimedean only, ceiling | $0$ | $-0.0766$ | $-11.531$ |
| archimedean only, realised contour | $0$ | $-0.5398$ | $-18.017$ |
| ~~adelic (pure 2, cond 0), ceiling~~ **SUPERSEDED** | $+0.2546$ | $+0.1780$ | $-7.966$ |
| ~~adelic (pure 2, cond 0), realised~~ **SUPERSEDED** | $+0.2546$ | $-0.2852$ | $-14.452$ |
| **adelic (pure $2$, cond $-2$), ceiling** | $\mathbf{-0.7780}$ | $\mathbf{-0.8546}$ | $\mathbf{-22.423}$ |
| **adelic (pure $2$, cond $-2$), realised** | $\mathbf{-0.7780}$ | $\mathbf{-1.3179}$ | $\mathbf{-28.909}$ |

> **The headline claim of §0 and §4.2 — "the adelic gain flips the entry test on
> the level-8 Catalan host" — is withdrawn.**  With the measured slopes
> $\gamma_2$ is **negative** on the $m=14$ inventory, the entry condition
> $\log|\varphi'(0)|>\tau_{\mathrm{ad}}$ **fails** ($-0.855$ at the ceiling), and
> the adelic bound gives no bound at all there.  The level-8 Catalan host never
> passed CDT's entry test.

The adelic gain is not worthless — it is worthless *in the presence of the
conditional functions*.  Optimising over sub-inventories at the ceiling
(`21run.log`):

| sub-inventory | $m$ | $\tau$ | $\gamma_2$ | entry | margin |
|---|---|---|---|---|---|
| best containing a conditional function: $7$ pure $+1$ | $8$ | $4.1172$ | $+0.3466$ | $+0.3883$ | $-7.353$ |
| best unconditional: the pure orbit alone | $7$ | $4.0782$ | $+1.1883$ | $+1.2689$ | $-1.576$ |

So on the level-8 host the best configuration that can prove anything sits at
$-7.35$, and the unconditional pure orbit alone sits at $-1.58$ (consistent, as
it must be — a positive margin there would be a false theorem; with this
$\mathbf e$ the theorem asserts the pure module has $\mathbf Q(y)$-dimension
$\le8$).

### 7.4 What this changes elsewhere in this document

* §0, row "**Level-8 Catalan host: the adelic gain flips the entry test**":
  **withdrawn** (§7.3).
* §0, row "…but the margin stays at $-7.97$ / $-14.45$": the correct figures for
  that inventory are $-22.42$ / $-28.91$.
* §4.2's slope table: the $y$-column entries $0$ (conditional) and $4$ ($B_0$)
  are replaced by $-2$ and $+3$; the $x$-column is unchanged and confirmed.
* §4.3's dimension cap is *unaffected* (it concerns the hypothetical doubly-small
  orbit, whose slope was measured in $x$ and whose $y$-slope enters only the
  hypothetical rows; those rows should be recomputed with $\varsigma_y=2\varsigma_x-2$
  and no floor).
* §5's structural ceiling $\log(256/\lambda_2)+\sigma_p\log p\le\log(256\lambda_1)$
  is unaffected (it is about $\varsigma_x$, not the descent).
* The $X_1(5)\,\mathrm{Sym}^2$ verdict is unaffected ($\gamma_p\equiv0$ there).

**Still no irrationality claim, and now one fewer positive claim.**
