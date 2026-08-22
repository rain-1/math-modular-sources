# The fold cusp of the level-16 Catalan host, a second doubly-small function, and the honest conditional count

*Fable (Opus 5), 2026-08-22.  Scripts and logs: `lattice/catalan_three_period/`
(`10_periods.gp`, `11_a0.gp`, `12_exact.gp`, `13_outer_dbl.gp`, `14_symm.gp`,
`15_rank.gp`, `16_margin.py`, `17_conditions.gp`, `18_slopes.gp`, `19_margin3.py`,
`20_level8_slopes.gp`, `21_level8_margin.py`; logs `10run.log`–`21run.log`).
Revises `consolidation/CATALAN_TWO_CLASSES.md` §1.1, §2.3, §4, §6 and
`consolidation/ADELIC_HOLONOMY.md` §4.3; conventions from
`consolidation/CDT_FINDER.md` §§1–3 and `consolidation/CDT_NONCONGRUENCE.md` §10.
Tags **[proved]** = derived here from the transformation law; **[verified]** =
exact or high-precision computation in this task; **[estimate]** = transported
input, flagged at each use.  **No irrationality claim is made anywhere.**

---

## 0. Verdict

| claim | verdict |
|---|---|
| The **fold cusp** of the level-16 host ($x=-1/4$, $\lambda_1=4$) is the cusp $\mathfrak a_1=\tfrac12$ of $\Gamma_0(16)$, of width $4$, with stabiliser generator $\gamma_1=\begin{pmatrix}-7&4\\-16&9\end{pmatrix}$ | **[verified]** §1.1–1.3 |
| Fold-regularity of the *power-series* companion $B_\Phi=F\,\mathcal D^{-2}\Phi_+$ is **two** conditions, not one: $a_0(\Phi)=0$ (vanishing of the constant term of $\Phi$) **and** $-16\alpha-8\beta=0$ on the Eichler period polynomial $r_{\gamma_1}=\alpha+\beta\tau$; then $\xi=-\beta/16$ | **[proved]** §1.2 |
| On the **inner** space $\langle E,V_2E,V_4E\rangle$ the first condition is vacuous and the second is exactly $8c_2+c_4=0$, with $\xi=-\tfrac12P(2)G$ — reproducing `CATALAN_TWO_CLASSES.md` §2.3 from the period polynomial | **[verified]** §1.3 |
| On the **outer** space $\langle T,V_2T,V_4T\rangle$ the *second* condition is **identically satisfied** ($\beta=-2\alpha$ always) and the *first* is the binding one: $a_0=-\tfrac14P(0)$, so fold-regular $\iff P(0)=0$. This is where §1.1's counting rule fails: the right count is still $d-1=2$, but for a completely different reason | **[proved]+[verified]** §1.4, §3 |
| On the outer fold-regular plane $\ \boxed{\xi_\infty=\tfrac32\zeta(2)\bigl(P(1)-P(2)\bigr)}\ $ — **not** Theorem B\*'s $\tfrac12P(2)\zeta(2)$. It reproduces the two anomalous measured limits of `CATALAN_TWO_CLASSES.md` §4 to $29$ digits | **[verified]** §1.4 |
| **A second, outer doubly-small function exists:** $\Phi_0^{\rm out}=(1+3V_2-4V_4)T$, $\xi_\infty=\xi_2=0$. Radius $\sqrt2/4>1/4$, rational coefficients, $d_n^2B\in\mathbf Z$ sharp, $2$-adic slope $1$ in $x$ and $0$ in $y$. **`CATALAN_TWO_CLASSES.md` §4 and §2.4 are wrong**: the doubly-small orbit has dimension $2$, not $1$ | **[verified]** §2 |
| §4 of `CATALAN_TWO_CLASSES.md` missed it because it tested each outer class against Theorem B\*'s prediction $\tfrac12P(2)\zeta(2)$, which is the wrong $\xi$ on this host; with the wrong $\xi$ every $|c_n|^{1/n}$ necessarily returns $\lambda_1=4$ | **[verified]** §3.1 |
| CDT's $m=14=7+7$: the $7$ conditional functions are the **$\theta$-orbit of one** symmetrised conditional function $\mathcal G=\mathrm{Sym}^+H$, $\{\mathcal G,\mathcal G',\mathcal G'',\mathcal G''',\int\mathcal G,\int\mathcal G/y,\int\mathcal G/y^2\}$. The three periods $1,\zeta(2),L(2,\chi_{-3})$ contribute **one** generator, not three | **[verified from `lattice/cdt_finder/indep_check2.py`]** §4.1 |
| Consequently the **one-period hypothesis $G\in\mathbf Q$ and the three-period hypothesis $a+b\zeta(2)+cG=0$ give exactly the same function count on this host** ($2$ unconditional $+1$ conditional generator). The stronger conclusion is free | **[proved]** §4.2 |
| The $21$ members of the three small orbits are $\mathbf Q(y)$-linearly independent (rank $21/21$, $42/42$, $63/63$ at $\deg\le0,1,2$) | **[verified]** §4.3 |
| **The same defect is present on the level-8 host and withdraws a headline claim there:** measured $\varsigma_y=-2$ for the level-8 conditional orbit, so $\gamma_2=-0.778$, the entry test **fails** ($-0.855$ at the ceiling) and the $m=14$ margin is $-22.42$ (ceiling) / $-28.91$ (transported contour), not $-7.97$/$-14.45$. `ADELIC_HOLONOMY.md` §7 records the correction | **[verified]** §5.1, `ADELIC_HOLONOMY.md` §7 |
| **The conditional functions' $2$-adic slope in $y$ is $-2$, not $0$.** $y=4x^2/(4x+1)$ is not an integral coordinate ($x(v)=\tfrac12(v^2+v\sqrt{1+v^2})$ carries $2$-power denominators), so the "$\max\{0,\cdot\}$ floor" used in `ADELIC_HOLONOMY.md` §3 and `CATALAN_TWO_CLASSES.md` §6 does not hold here. Measured: $\mathrm{Sym}^+A$ and $\mathrm{Sym}^+(B_E+\mu A)$ both $\to-2$ | **[verified]** §5.0 |
| **Corrected margins at the ceiling.** With the previously *assumed* profile: previous $-7.966$/$-7.551$, best inventory containing a conditional function $-2.006$, best *unconditional* inventory $-0.737$. With the *measured* profile: previous configurations $-22.42$/$-22.43$, best with a conditional function $\mathbf{-5.795}$, best unconditional $-0.737$ (unchanged) | **[verified]**, §5 |
| **The conditional function is a net cost either way.** Adding one to the best unconditional inventory costs $1.27$ nats on the assumed profile and $5.06$ nats on the measured one. On this host $\gamma_2$ is the only thing holding the entry test open, and every conditional function dilutes it | **[verified]** §5.3 |
| The $-0.737$ unconditional row is a **consistency alarm**: with the transported pure inventory the adelic bound is within $0.74$ nats of asserting something false. That is a bound on how optimistic the transported pure inventory can be | **[verified]** §5.2 |

One sentence: *the level-16 host has two unconditional doubly-small functions, not
one — the missing one is the outer class $(1+3V_2-4V_4)T$, and it was missed
because Theorem B\* does not compute the period at this host's fold cusp; the
three-period hypothesis is free but buys nothing; and once the counting and the
$2$-adic slopes are done honestly the conditional function is a net liability, so
the margin story ends at $-2.0$ (assumed slopes) or $-5.8$ (measured slopes)
rather than $-7.6$, with no route to $0$ in either case.*

---

## 1. The fold cusp and the exact criterion

### 1.1 Which cusp is the fold

$\Gamma_0(16)$ has six cusps, $c\in\{1,2,4,4',8,16\}$ (denominators), of widths
$16,4,1,1,1,1$.  From the Ligozat orders of $x=\eta_2\eta_{16}^2/(\eta_1^2\eta_8)$,
$x$ has its only pole at the cusp $0$ and its only zero at $\infty$; the four
remaining cusps carry the four finite non-zero values
$\{-\tfrac12,-\tfrac14,\tfrac{-1\pm i}4\}$ of `CATALAN_TWO_CLASSES.md` §2.1.  The
Galois-conjugate pair $\{\tfrac14,\tfrac34\}$ must carry $\tfrac{-1\pm i}4$, and
the degree-$2$ covering $X_0(16)\to X_0(8)$ (unramified over $\infty_8$, ramified
over $0_8$ and $(\tfrac12)_8$) forces $x(\tfrac18)=-\tfrac12$ and
$$\boxed{\ x\bigl(\tfrac12\bigr)=-\tfrac14=x_1\quad\text{(the fold)},\qquad
h=4,\qquad \gamma_1=\sigma T^4\sigma^{-1}=\begin{pmatrix}-7&4\\-16&9\end{pmatrix},
\ \ \sigma=\begin{pmatrix}1&0\\2&1\end{pmatrix}. }$$
This identification is **independently confirmed**: the period polynomials were
computed at *all five* non-$\infty$ cusps (`10run.log`), and only
$\mathfrak a_1=\tfrac12$ produces the fold-regularity condition $8c_2+c_4=0$ that
`CATALAN_TWO_CLASSES.md` §2.3 verified numerically on twenty classes:

| cusp | $x$-value | width | $\gamma$ | condition $\ c\,\alpha-(d-1)\beta=0$ | condition on $(c_1,c_2,c_4)$, inner |
|---|---|---|---|---|---|
| $0$ | $\infty$ | $16$ | $\begin{pmatrix}1&0\\-16&1\end{pmatrix}$ | $\alpha=0$ | $64c_1+8c_2+c_4=0$ |
| $\mathbf{1/2}$ | $\mathbf{-1/4}$ | $4$ | $\begin{pmatrix}-7&4\\-16&9\end{pmatrix}$ | $2\alpha+\beta=0$ | $\mathbf{8c_2+c_4=0}$ ✅ |
| $1/4$ | $\tfrac{-1+i}4$ | $1$ | $\begin{pmatrix}-3&1\\-16&5\end{pmatrix}$ | $4\alpha+\beta=0$ | $c_4=0$ |
| $3/4$ | $\tfrac{-1-i}4$ | $1$ | $\begin{pmatrix}-11&9\\-16&13\end{pmatrix}$ | $4\alpha+3\beta=0$ | $c_4=0$ |
| $1/8$ | $-1/2$ | $1$ | $\begin{pmatrix}-7&1\\-64&9\end{pmatrix}$ | $8\alpha+\beta=0$ | *vacuous* |

**[verified]** — `17run.log` recovers the condition coefficients exactly as
$(64,8,1)$, $(0,8,1)$, $(0,0,1)$, $(0,0,1)$ and *vacuous*.  Each cusp cuts out a
different plane, and only the cusp $\tfrac12$ cuts out the one whose members were
measured to have the enlarged radius $\sqrt2/4$.  (The last row says every inner
class is already single-valued at $x=-\tfrac12$, which is consistent: that
singularity sits at $|x|=\tfrac12$, outside the disc of radius $\sqrt2/4$ that the
fold-regular companions actually reach.)

Note the inversion relative to level 8, where the fold **is** the cusp $0$ and the
condition **is** $\alpha=0$; that is the source of the phrasing "vanishing of the
constant coefficient" in `CATALAN_TWO_CLASSES.md` §1.1, which is correct at level
8 and wrong at level 16.

### 1.2 The criterion, derived

Let $\Phi\in M_3(\Gamma_0(16),\chi_{-4})$ with constant term $a_0$, and let
$\Theta=\mathcal D^{-2}\Phi_+$ be the project's Eichler integral (positive part
only — this is the object in $B_\Phi=F\Theta$; adding the constant term would add
$-2\pi^2a_0\tau^2F$, which is not a power series in $x$).  Put
$\widehat\Theta=\Theta-2\pi^2a_0\tau^2$, so that $\mathcal D^2\widehat\Theta=\Phi$
exactly.  Bol's identity gives, for $\gamma=\begin{pmatrix}a&b\\c&d\end{pmatrix}\in\Gamma_0(16)$,
$$(c\tau+d)\,\widehat\Theta(\gamma\tau)=\chi_{-4}(d)\bigl(\widehat\Theta(\tau)+\alpha+\beta\tau\bigr),$$
$\alpha+\beta\tau=r_\gamma$ the Eichler period polynomial.  Hence for the actual
$\Theta$,
$$R(\tau):=\chi_{-4}(d)^{-1}(c\tau+d)\Theta(\gamma\tau)-\Theta(\tau)
=\alpha+\beta\tau+a_0\,g_\gamma(\tau),\qquad
g_\gamma(\tau)=2\pi^2\Bigl[\tfrac{(a\tau+b)^2}{\chi_{-4}(d)(c\tau+d)}-\tau^2\Bigr].$$
$g_\gamma$ has a **pole** at $\tau=-d/c$ of residue $2\pi^2a_0/c^2\ne0$ whenever
$a_0\ne0$ (the numerator of $g_\gamma$ takes the value $1/c^2$ there, by
$ad-bc=1$).  With $F|_1\gamma=\chi_{-4}(d)F$ and $\chi_{-4}(d)^2=1$,
$$H_\xi:=F\,(\Theta-\xi)\ \Longrightarrow\
H_\xi\circ\gamma-H_\xi=F\cdot\bigl[R(\tau)-\xi\bigl(\chi_{-4}(d)(c\tau+d)-1\bigr)\bigr].$$

> **Proposition 1 [proved].**  $H_\xi$ is invariant under $\gamma_1$ — i.e.
> $B_\Phi-\xi A$ is single-valued at the fold — **iff**
> $$a_0(\Phi)=0\qquad\text{and}\qquad -16\alpha-8\beta=0,\qquad\text{and then}\qquad \xi=-\tfrac{\beta}{16}=\tfrac{\alpha}8 .$$
> If $a_0=0$ but $-16\alpha-8\beta\ne0$, then in the cusp-local frame
> $r_{\gamma_1}|_{-1}\sigma=(2\alpha+\beta)\tau'+\alpha$ and $B_\Phi-\xi A$ acquires a
> $\log^2$ at $x_1$ for **every** $\xi$; the coefficient $2\alpha+\beta$ of that
> $\log^2$ is independent of $\xi$.  This is the "one linear condition" of
> `CATALAN_TWO_CLASSES.md` §1.1, correctly located: it is the vanishing of the
> **leading** (degree-$1$) coefficient of the period polynomial in the cusp frame,
> and $\xi$ is then read off the constant coefficient.  If $a_0\ne0$ no $\xi$
> works at all, because $R$ is not a polynomial.

### 1.3 The computation, and the cross-validation on the inner orientation

`10_periods.gp`, `11_a0.gp`, `12_exact.gp`: $\Theta$ and $\Theta'$ are summed to
$130$ digits at the symmetric point $\tau_0=\tfrac{9}{16}+\tfrac i{16}$ of the
$\gamma_1$-horocycle (both $\tau_0$ and $\gamma_1\tau_0$ have $\operatorname{Im}=1/|c|$),
and $(\alpha,\beta,a_0)$ are solved from $R(\tau_0)$, $R'(\tau_0)$, $R(\tau_1)$ at a
second point.  Agreement between the two points is $\sim10^{-131}$.  At the fold
cusp $\mathfrak a_1=\tfrac12$ (`12run.log`), with $G$ = Catalan:

| $\Phi$ | $a_0$ | $\alpha$ | $\beta$ | $-16\alpha-8\beta$ | $\xi=-\beta/16$ |
|---|---|---|---|---|---|
| $E$ | $0$ | $-4G$ | $8G$ | $0$ | $-\tfrac12G$ |
| $V_2E$ | $0$ | $-G-\tfrac5{16}i\pi^2$ | $2G+\tfrac12i\pi^2$ | $i\pi^2$ | $-\tfrac18G-\tfrac1{32}i\pi^2$ |
| $V_4E$ | $0$ | $-\tfrac14G-\tfrac5{128}i\pi^2$ | $\tfrac12G+\tfrac1{16}i\pi^2$ | $\tfrac18i\pi^2$ | $-\tfrac1{32}G-\tfrac1{256}i\pi^2$ |
| $T$ | $-\tfrac14$ | $4\zeta(2)$ | $-8\zeta(2)$ | $0$ | $\tfrac12\zeta(2)$ |
| $V_2T$ | $-\tfrac14$ | $7\zeta(2)$ | $-14\zeta(2)$ | $0$ | $\tfrac78\zeta(2)$ |
| $V_4T$ | $-\tfrac14$ | $\tfrac{25}4\zeta(2)$ | $-\tfrac{25}2\zeta(2)$ | $0$ | $\tfrac{25}{32}\zeta(2)$ |

(all entries recognised by `bestappr` from $130$-digit values; $a_0=-\tfrac14$ is
$L(-2,\chi_{-4})/2$, as it must be.)  Writing $\Phi=(c_1+c_2V_2+c_4V_4)E$,
$P(s)=c_1+c_22^{-s}+c_44^{-s}$, $\ell=8c_2+c_4$:
$$\alpha=-4G\,P(2)-\tfrac{5i\pi^2}{128}\ell,\qquad
\beta=8G\,P(2)+\tfrac{i\pi^2}{16}\ell,\qquad
-16\alpha-8\beta=\tfrac{i\pi^2}8\,\ell .$$

> **The inner cross-validation.**  Fold-regular $\iff\ell=8c_2+c_4=0$, and then
> $\xi=-\tfrac12P(2)G$.  This is *exactly* `CATALAN_TWO_CLASSES.md` §2.3, now
> derived from the period polynomial rather than fitted to twenty numerical
> experiments.  **[verified]**

### 1.4 The outer orientation

For $\Phi=(c_1+c_2V_2+c_4V_4)T$, $T=E_3^{\mathbf1,\chi_{-4}}$, the same table gives
$$\boxed{\ \alpha=4\zeta(2)\,Q,\qquad \beta=-8\zeta(2)\,Q,\qquad
Q:=P(0)+3P(1)-3P(2),\qquad a_0=-\tfrac14P(0).\ }$$
So $\beta=-2\alpha$ **identically**: the functional $-16\alpha-8\beta$, which cuts
out fold-regularity on the inner space, is **identically zero** on the outer
space.  The binding condition is the other one:
$$\boxed{\ \Phi\ \text{outer fold-regular}\iff P(0)=0,\qquad\text{and then}\quad
\xi_\infty=\tfrac{\zeta(2)}2Q=\tfrac32\zeta(2)\bigl(P(1)-P(2)\bigr).\ }$$

Two independent confirmations from the *old* numerics (`catalan_two_classes/03run.log`),
which had been read as evidence against the outer orientation:

| class | $P(0)$ | predicted $\xi_\infty$ | measured $b_n/a_n$ at $n=218$ | Theorem B\*'s $\tfrac12P(2)\zeta(2)$ |
|---|---|---|---|---|
| $T-V_2T$ | $0$ | $-\tfrac38\zeta(2)=-0.61685027506808491368$ | $-0.61685027506808491368$ (29 digits) | $+\tfrac38\zeta(2)$ — **wrong sign** |
| $T-5V_2T+4V_4T$ | $0$ | $-\tfrac34\zeta(2)=-1.2337005501361698274$ | $-1.2337005501361698274$ (29 digits) | $0$ — **wrong value** |

Theorem B\* computes $L(\Phi,2)$, which is the period at the cusp $0$; on this
host the fold is the cusp $\tfrac12$, and the correct functional mixes $P(1)$ and
$P(2)$ — as it must, a period polynomial at a general cusp seeing both critical
values.  §4 of `CATALAN_TWO_CLASSES.md` compared measurements against $\tfrac12P(2)\zeta(2)$
and, finding disagreement, concluded "the outer orientation contributes nothing".
The disagreement was in the prediction.

---

## 2. The outer doubly-small function

### 2.1 Definition

$$\boxed{\ \Phi_0^{\rm out}=(1+3V_2-4V_4)\,T=T+3V_2T-4V_4T,\qquad
P(0)=0,\quad P(1)=P(2)=\tfrac32,\quad \xi_\infty=0.\ }$$
Equivalently $\Phi_0^{\rm out}=2(T-V_2T)-(T-5V_2T+4V_4T)$, and its period
polynomial at the fold cusp **vanishes identically**: $\alpha=\beta=0$ and
$a_0=0$.  Its companion begins
$$B_{\Phi_0^{\rm out}}(x)=x-x^2-\tfrac{44}9x^3+\tfrac{88}3x^4-\tfrac{20416}{225}x^5
+\tfrac{2704}{15}x^6-\tfrac{1642304}{11025}x^7-\tfrac{787712}{1575}x^8
+\tfrac{260982784}{99225}x^9-\cdots$$

### 2.2 Verification (`13run.log`, `14run.log`; $N=220$ terms, $n=218$)

| property | $\Phi_0^{\rm in}=(1+4V_2-32V_4)E$ | $\Phi_0^{\rm out}=(1+3V_2-4V_4)T$ |
|---|---|---|
| rational coefficients | exact | exact |
| $b_n/a_n$ at $n=218$ | $1.03\cdot10^{-31}$ | $9.19\cdot10^{-31}$ |
| $|b_n|^{1/n}$ at $n=60,120,218$ | $2.688,\,2.744,\,2.748$ | $2.649,\,2.722,\,2.775$ |
| limit | $\to2\sqrt2=2.8284$: radius $\sqrt2/4$ ✅ | $\to2\sqrt2$: radius $\sqrt2/4$ ✅ |
| $v_2(b_n)/n$ at $n=60,120,218$ | $1.467,\,1.483,\,1.495\to\tfrac32$ | $0.983,\,1.000,\,0.991\to1$ |
| $\xi_2$ (since $v_2(a_n)=O(1)$) | $0$ ✅ | $0$ ✅ |
| $d_n^kB\in\mathbf Z$ | $k=2$ yes, $k=1$ fails at $n=3$ | $k=2$ yes, $k=1$ fails at $n=3$ |
| $2$-adic slope in $y$ | $0.95,\,0.975,\,1.000\to1$ | $0.050,\,0.050,\,-0.010\to0$ |
| $y$-radius $|d_k|^{1/k}$ | $1.80,\,1.87,\,1.93\to2$ | $1.70,\,1.81,\,1.90\to2$ |

Both are singular at $|y|=1/2$, as §5.3 of `CATALAN_TWO_CLASSES.md` predicts for
everything that is not the pure module.  The $y$-slope obeys the same rule
$\varsigma_y=2(\varsigma_x-1)$: $\varsigma_x=\tfrac32\Rightarrow\varsigma_y=1$,
$\varsigma_x=1\Rightarrow\varsigma_y=0$.

$\xi_2=0$ for the outer orientation is *forced* (Theorem F: $\varphi\ne\mathbf1$),
and here it is also *measured*: $v_2(b_n)\sim n$ against $v_2(a_n)=O(1)$
($v_2(a_n)=8,6$ at $n=60,100$).

### 2.3 What it changes

* **The doubly-small orbit on this host has dimension $2$, not $1$.**
  `CATALAN_TWO_CLASSES.md` §0 and §2.4 ("the doubly-small orbit has dimension
  exactly 1", "the outer orientation contributes no target-zero fold-regular
  direction") are **retracted**.
* $B_{\Phi_0^{\rm out}}$ is a **new unconditional object**: a rational power series
  of radius $\sqrt2/4$ on a host of radius $1/4$, of denominator type $d_n^2$
  (sharp), $2$-adically overconvergent of slope $1$ in $x$.  It is
  $\mathbf Q(y)$-independent of $B_{\Phi_0^{\rm in}}$, of $A$, and of the
  conditional function (§4.3).
* It does **not** change any irrationality status.  Its $y$-slope is $0$ — two
  below the pure module's and one below $\Phi_0^{\rm in}$'s — so adding its orbit
  *dilutes* the adelic gain; §5 shows the full $m=28$ inventory is margin-worse
  than the $m=14$ one that omits it.

### 2.4 The full fold-regular picture on the $6$-dimensional space

$$\underbrace{\{8c_2^{\rm in}+c_4^{\rm in}=0\}}_{2\text{-dim, period }\in\mathbf Q\,G}
\ \oplus\ \underbrace{\{c_1^{\rm out}+c_2^{\rm out}+c_4^{\rm out}=0\}}_{2\text{-dim, period }\in\mathbf Q\,\zeta(2)}
\qquad(\dim=4)$$
with $\xi_\infty(\Phi)=-\tfrac12P_{\rm in}(2)\,G+\tfrac32\bigl(P_{\rm out}(1)-P_{\rm out}(2)\bigr)\zeta(2)$.
The two conditions $\xi_\infty=0$ (one per period) cut out the $2$-dimensional
doubly-small plane $\langle\Phi_0^{\rm in},\Phi_0^{\rm out}\rangle$.  There are **no
genuinely mixed** extra fold-regular directions: the fold-regularity functionals
separate, one living only on the inner coordinates and one only on the outer.

---

## 3. Where the counting rule of `CATALAN_TWO_CLASSES.md` §1.1 fails

> **Theorem-shaped claim (the outer obstruction) [proved].**
> Let $\Phi\in M_3(\Gamma_0(N),\chi)$ with constant term $a_0$, let $\mathfrak a$ be a
> cusp with $\mathfrak a\ne\infty$ (so its stabiliser generator has $c\ne0$), and let
> $B_\Phi=F\cdot\mathcal D^{-2}\Phi_+$ be the power-series companion.  Then the
> transformation defect $R$ of $\mathcal D^{-2}\Phi_+$ under that generator is a
> polynomial of degree $\le1$ **iff $a_0=0$**; otherwise it is a rational function
> with a simple pole at $\tau=-d/c$, of residue $2\pi^2a_0/c^2$.  Consequently no
> constant $\xi$ makes $B_\Phi-\xi A$ single-valued at $\mathfrak a$, and $\Phi$ is not
> fold-regular there.  *(The obstruction is exactly the $\tau^2F$ term one would
> have to add to $B_\Phi$ to restore modularity, and $\tau^2F$ is not a power
> series in $x$: it carries $\log^2x$.)*

So the counting rule "fold-regularity is one linear condition, hence the
fold-regular subspace has dimension $\ge d-1$" is:

* **correct as stated on the inner orientation** ($a_0\equiv0$ there, and the one
  condition is $-16\alpha-8\beta=0$, i.e. $8c_2+c_4=0$);
* **numerically correct but mechanically wrong on the outer orientation**: the
  condition $-16\alpha-8\beta=0$ is vacuous, and the one condition that actually
  cuts the space is $a_0=0$, i.e. $P(0)=0$.  Both give $\dim=d-1=2$; the
  coincidence of dimensions is why the rule was never caught.
* The rule's second half — "$\xi$ is determined by $\beta$" — is correct in both
  cases, but the *value* of $\xi$ is $-\beta/16$, which equals Theorem B\*'s
  $L(\Phi,2)$ on the inner orientation and does **not** on the outer.

### 3.1 Why §4 of `CATALAN_TWO_CLASSES.md` missed the class

`catalan_two_classes/03_level16.gp` computes $|c_n|^{1/n}$ for $c_n=b_n-\xi a_n$
with $\xi$ set to Theorem B\*'s $\tfrac12P(2)\zeta(2)$.  For a fold-regular class
with true period $\xi^\ast\ne\xi$, $c_n=(\xi^\ast-\xi)a_n+O(\lambda_2^n)$, so
$|c_n|^{1/n}\to\lambda_1=4$ *by construction*.  Every outer diagnostic in §4 was
therefore a test of Theorem B\*, not of fold-regularity.  With the correct $\xi$
(`13run.log`) the same classes give $|c_n|^{1/n}=2.775\to2\sqrt2$:

| class | $P(0)$ | $|c_n|^{1/n}$ at $n=218$, correct $\xi$ | verdict |
|---|---|---|---|
| $T$, $V_2T$, $V_4T$ | $1$ | $3.815\to4$ | not fold-regular ✅ (constant term) |
| $V_2T-8V_4T$ | $-7$ | $3.850\to4$ | not fold-regular ✅ |
| $T-V_2T$ | $0$ | $2.7753\to2\sqrt2$ | **fold-regular** |
| $T-5V_2T+4V_4T$ | $0$ | $2.7753\to2\sqrt2$ | **fold-regular** |
| $T+V_2T-2V_4T$, $2T+V_2T-3V_4T$ | $0$ | $2.7753,\,2.7842$ | **fold-regular** |
| $\Phi_0^{\rm out}=T+3V_2T-4V_4T$ | $0$ | $2.7753$, and $\xi=0$ | **doubly small** |

---

## 4. The honest conditional count

### 4.1 What CDT's "$7$" actually is

Reading `lattice/cdt_finder/indep_check2.py`, which rebuilds all fourteen of
CDT's functions and verifies their $\mathbf Q(y)$-independence (rank $14,28,42,56,70,84$
at $\deg P_i\le0,\dots,5$), the inventory is

* **pure ($7$):** $B_1=1$; $B_2,B_3,B_5,B_6$ (factorial series);
  $B_4=\mathrm{Sym}^-\mathrm{Li}_2$; $B_7=\int B_4\,dy/y$;
* **conditional ($7$):** $\mathcal G,\ \mathcal G',\ \mathcal G'',\ \mathcal G''',\
  \int\mathcal G\,dy,\ \int\mathcal G\,dy/y,\ \int\mathcal G\,dy/y^2$,

where $\mathcal G=\mathrm{Sym}^+H$ and $H=aH_A+bH_B+cH_C$ is the **single**
function fixed by the hypothesised relation $a+b\zeta(2)+cL(2,\chi_{-3})=0$.
So the multiplication is
$$\boxed{\ 7=\#\{\text{$\mathbf Q(y)$-independent members of the $\theta$-orbit of one small function}\}
=\underbrace{4}_{\text{$\partial^j$, }j\le3}+\underbrace{3}_{\int\cdot\,y^{-i},\ i\le2},\ }$$
**not** a factor $3$ coming from the three periods.  The $\mathbf Z/2$ of the
normaliser descent does not double the conditional count either: it is used once,
to push $H$ down from $\mathbf P^1\setminus\{0,s,\infty\}$ to the $y$-line (raising
the ceiling $16|s|\to256|s|$ and doubling the denominators $b_j:1\to2$), and once
more on the *pure* side, where both $\mathrm{Sym}^+$ and
$\mathrm{Sym}^-/(x-\sigma x)$ are used ($B_4$).  The integration profile
$\mathbf e$ with $\sum e_i=6$, $\max e_i=1$ is exactly $\{B_5,B_6,B_7\}$ plus the
three conditional integrals.

$H_A,H_B,H_C$ are a **basis of the inhomogeneous solution space**, i.e. the
$c=3$ of `CDT_NONCONGRUENCE.md` §10 counts the *parameters available to satisfy
the fold condition*, not the number of usable functions.  Only their hypothesis-fixed
combination is admissible.

### 4.2 The level-16 inventory, honestly

The host carries exactly **two non-trivial** periods.  The nebentypus is forced: $F$ has
character $\chi_{-4}$, so $\Phi=F\mathcal D x$ must lie in $M_3(\Gamma_0(16),\chi_{-4})$,
and the only pairs $(\psi,\varphi)$ with $\psi\varphi=\chi_{-4}$ and
$\operatorname{cond}\psi\cdot\operatorname{cond}\varphi\mid16$ are
$(\chi_{-4},\mathbf1)$ and $(\mathbf1,\chi_{-4})$.  So the period set is
$\{1,\zeta(2),G\}$ and $c=3$ — well below `CDT_NONCONGRUENCE.md` §10's cusp bound
$c\le\#\text{cusps}-1=5$.  *The binding constraint here is the character, not the
cusp count.*

Let $V$ be the $4$-dimensional fold-regular space of §2.4 and consider
$H=B_\Phi+\mu A$ with $\Phi\in V$ (rational coordinates) and $\mu\in\mathbf Q$; its
period is $\mu+r(\Phi)\,G+s(\Phi)\,\zeta(2)$, where
$r(\Phi)=-\tfrac12P_{\rm in}(2)$ and $s(\Phi)=\tfrac32(P_{\rm out}(1)-P_{\rm out}(2))$
are jointly surjective onto $\mathbf Q^2$ with $2$-dimensional kernel.  "$H$ small"
means that period vanishes; $H$ is *usable* only if its coefficients are rational,
which is automatic here.

| hypothesis | condition on $(\mu,r,s)$ | small space | unconditional part | **new conditional generators** |
|---|---|---|---|---|
| none | $\mu=r=s=0$ | $2$ | $2$ | $0$ |
| **(a)** $G\in\mathbf Q$ | $s=0$, $\mu=-rG$ | $3$ | $2$ | $\mathbf 1$ |
| **(b)** $a+b\zeta(2)+cG=0$, $c\ne0$ | $cs=br$, $\mu=ra/c$ | $3$ | $2$ | $\mathbf 1$ |

(In (b), $c\ne0$ is automatic since $1,\zeta(2)$ are $\mathbf Q$-independent —
$\pi^2\notin\mathbf Q$.  Under (a) the conditional generator is purely inner;
under (b) it is a genuinely mixed inner$+$outer class.)

> **Proposition 2 [proved].**  On the level-16 Catalan host the one-period
> hypothesis $G\in\mathbf Q$ and the three-period hypothesis
> $a+b\zeta(2)+cG=0$ yield the **same** number of small functions:
> $2$ unconditional $+1$ conditional generator.  The three-period statement — whose
> contradiction proves that $1,\zeta(2),G$ are $\mathbf Q$-linearly independent, and
> in particular $G\notin\mathbf Q$ — is therefore **free**, and is the statement one
> should target.

Applying CDT's $\times7$ orbit multiplier:
$$m=\underbrace{7}_{\text{pure [estimate]}}+\underbrace{7}_{\Phi_0^{\rm in}}
+\underbrace{7}_{\Phi_0^{\rm out}}+\underbrace{7}_{\text{conditional}}=28,
\qquad\text{of which } 21 \text{ are unconditional.}$$

### 4.3 Independence, verified (`15run.log`)

Working in $y=4x^2/(4x+1)$, with $g_1=\mathrm{Sym}^+B_{\Phi_0^{\rm in}}$,
$g_2=\mathrm{Sym}^+B_{\Phi_0^{\rm out}}$, $g_3=\mathrm{Sym}^+(B_E+\tfrac\gamma2A)$
($\gamma\in\mathbf Q$ generic, standing for $G$):
$$g_1=y-2y^2+\tfrac{32}9y^3-\tfrac{56}9y^4+\cdots,\qquad
g_2=\tfrac12y-y^2+\tfrac{31}{18}y^3-\tfrac{26}9y^4+\cdots$$

| set | $\deg\le0$ | $1$ | $2$ | $3$ |
|---|---|---|---|---|
| $\{g_1,g_2,g_3\}$ | $3/3$ | $6/6$ | $9/9$ | $12/12$ (and $15/15$, $18/18$ at $\deg\le4,5$) |
| $\{g_1,g_2,g_3,\mathrm{Sym}^+A\}$ | $4/4$ | $8/8$ | $12/12$ | $16/16$ |
| the three $7$-member orbits ($21$) | $21/21$ | $42/42$ | $63/63$ | $83/84$ |
| $+$ a candidate pure set ($28$) | $27/28$ | $54/56$ | $81/84$ | — |

**[verified]**, ranks mod $2^{61}-1$, series to $y^{88}$.  The three small
generators and their orbits are $\mathbf Q(y)$-independent (the $83/84$ at
$\deg\le3$ is at the edge of the available series length and is not resolved).
The single deficiency at $28$ is inside the *pure* block: of the seven candidates
$\{1,\mathrm{Sym}^\pm\mathrm{Li}_j(-4x)\}_{j=1,2,3}$ only six are independent.
**No pure module on this host has been constructed; the $7$ pure functions with
$u=(1,3)$ remain CDT's, transported. [estimate]**

---

## 5. Margins (`16run.log`, `18run.log`–`21run.log`)

$\mathrm{BC}=11.845+\log\tfrac14=10.458706$; ceiling $\log(256|s|)=4.158883$;
"realised" $=$ ceiling $+\log0.62922=3.695614$ **[estimate, and see §5.4]**;
"rigorous" $=\log(\pi^2/\Gamma(3/4)^4)=1.476336$ (`CATALAN_TWO_CLASSES.md` §5.3).

### 5.0 The $2$-adic slope profile, measured rather than assumed

`ADELIC_HOLONOMY.md` §3 takes $\varsigma_y\ge\max\{0,2\varsigma_x-2v_p(\lambda_2)\}$,
"the $0$ coming from the integrality already assumed in the denominator type",
and `CATALAN_TWO_CLASSES.md` §6 uses that floor to give the conditional orbit
$\varsigma_y=0$.  **The floor does not hold on this host.**  The symmetrised
coordinate $y=4x^2/(4x+1)$ is not integral: inverting $y=v^2$ gives
$x(v)=\tfrac12\bigl(v^2+v\sqrt{1+v^2}\bigr)=\tfrac v2+\tfrac{v^2}2+\tfrac{v^3}4-\cdots$,
whose coefficients carry $2$-power **denominators**, so $\mathrm{Sym}^+$ of an
integral $x$-series need not be integral in $y$.  Measured (`18run.log`, $200$
terms, $k\le96$):

| object | $v_2(d_k)/k$ at $k=20,40,96$ | $\varsigma_y$ |
|---|---|---|
| $\mathrm{Sym}^+\mathrm{Li}_2(-4x)$ | $1.700,\ 1.800,\ 1.875$ | $+2$ |
| $\mathrm{Sym}^+\mathrm{Li}_3(-4x)$ | $1.450,\ 1.650,\ 1.802$ | $+2$ |
| $\mathrm{Sym}^+B_{\Phi_0^{\rm in}}$ | $0.950,\ 0.975,\ 0.990$ | $+1$ |
| $\mathrm{Sym}^+B_{\Phi_0^{\rm out}}$ | $0.050,\ 0.050,\ 0.042$ | $\ \ 0$ |
| $\mathrm{Sym}^+A$ (the host) | $-1.750,\ -1.875,\ -1.948$ | $\mathbf{-2}$ |
| $\mathrm{Sym}^+B_E$ (raw companion) | $-1.850,\ -1.925,\ -1.969$ | $\mathbf{-2}$ |
| $\mathrm{Sym}^+(B_E+\mu A)$, $\mu\in\mathbf Q$ (conditional) | $-1.850,\ -1.925,\ -1.969$ | $\mathbf{-2}$ |

**[verified]**  Only relative slopes matter: the bound is exactly covariant under
$y\mapsto y/4$, which adds $+2$ to every slope and subtracts $\log4$ from both
$\log|\varphi'(0)|$ and $\mathrm{BC}$ (checked to $7\cdot10^{-15}$, `19run.log`).
In the integral coordinate $Y=y/4$ the profile is
$$(\text{pure},\ \Phi_0^{\rm in}\text{-orbit},\ \Phi_0^{\rm out}\text{-orbit},\ \text{conditional})
=(4,\ 3,\ 2,\ 0),$$
against the $(2,1,0,0)$ assumed hitherto.  Below, **P1** is the assumed profile
and **P2** the measured one; both are reported, because P1 is what the earlier
documents used.

*A structural aside*: $\mathrm{Sym}^+B_{T-V_2T}$ vanishes **identically** (checked
to $x^{58}$) — that outer companion is $\sigma$-anti-invariant, so an outer
conditional function has to be carried by $\mathrm{Sym}^-/(x-\sigma x)$.
$\mathrm{Sym}^+B_{\Phi_0^{\rm out}}=2y^2-8y^3+16y^4-\tfrac{1312}9y^6+\cdots$ is not
degenerate, which is what §4.3's rank check needs.

### 5.1 The table (ceiling)

| inventory | $m$ | $\tau^\flat$ | $\tau^\sharp$ | $\tau$ | $\gamma_2$ (P1) | margin (P1) | $\gamma_2$ (P2) | **margin (P2)** |
|---|---|---|---|---|---|---|---|---|
| *previous:* $7$ pure $+7$ conditional | $14$ | $3.8980$ | $0.3375$ | $4.2355$ | $+0.2546$ | $-7.966$ | $-0.7780$ | $-22.423$ |
| *previous:* $+1$ doubly-small | $15$ | $3.9111$ | $0.3158$ | $4.2269$ | $+0.2619$ | $-7.551$ | $-0.7301$ | $-22.431$ |
| (a)/(b) full: $7$ pure $+3$ small orbits | $28$ | $3.9745$ | $0.2177$ | $4.1922$ | $+0.1812$ | $-6.316$ | $-0.4482$ | $-23.942$ |
| **best containing a conditional function** ($7$ pure $+7$ $\Phi_0^{\rm in}$ $+1$) | $15$ | $3.9111$ | $0.3158$ | $4.2269$ | $+0.6315$ | $\mathbf{-2.006}$ | $+0.3789$ | $\mathbf{-5.795}$ |
| **best unconditional** ($7$ pure $+7$ $\Phi_0^{\rm in}$) | $14$ | $3.8980$ | $0.3375$ | $4.2355$ | $+0.7709$ | $-0.737$ | $+0.7709$ | $-0.737$ |
| **level 8** (Zagier E), $7$ pure $+7$ conditional — `ADELIC_HOLONOMY.md` §4.2 | $14$ | $3.8980$ | $0.3375$ | $4.2355$ | $+0.2546$ | $-7.966$ | $-0.7780$ | $-22.423$ |
| **level 8**, best containing a conditional fn ($7$ pure $+1$) | $8$ | $3.6875$ | $0.4297$ | $4.1172$ | $+0.7798$ | $-3.887$ | $+0.3466$ | $-7.353$ |
| **level 8**, best unconditional (pure orbit alone) | $7$ | $3.5918$ | $0.4864$ | $4.0782$ | $+1.1883$ | $-1.576$ | $+1.1883$ | $-1.576$ |

(the "best unconditional" rows contain no slope-$(-2)$ member, so P1 and P2 agree
there.  The level-8 rows are computed in §7 of `ADELIC_HOLONOMY.md` from the
*measured* level-8 slopes, `20run.log`/`21run.log`: $\varsigma_y=+2$ for the whole
pure orbit — confirming `pure_2adic.py` — and $\varsigma_y=-2$ for
$\mathrm{Sym}^+A$, for $\mathrm{Sym}^+B_{\mathbf E}$ and for
$\mathrm{Sym}^+(B_{\mathbf E}+\mu A)$, i.e. for the conditional orbit.  Level 8 and
level 16 share ceiling, $\mathrm{BC}$ and $\tau$, so the two $m=14$ rows coincide
exactly.  **This withdraws `ADELIC_HOLONOMY.md`'s headline claim that the adelic
gain flips the entry test on the level-8 Catalan host**: with the measured slopes
$\gamma_2=-0.778$ and the entry condition fails by $0.855$ at the ceiling; on the
transported contour that level-8 row is $-28.909$.)

With the transported contour the first three level-16 rows become $-14.452,\,-14.500,\,-19.288$ (P1)
and $-28.909,\,-29.380,\,-36.914$ (P2); with the rigorous thrice-punctured ceiling
the entry test fails outright in every row and there is no bound at all — the
$-2.76$ contour issue of `CATALAN_TWO_CLASSES.md` §5.3 is untouched and remains
decisive.

**Hypotheses (a) and (b) give identical rows** (Proposition 2): they are the same
computation, so there is nothing to put side by side.

### 5.2 The consistency alarm

The best **unconditional** inventory sits at $-0.737$, on either profile.  A
positive margin there would be a contradiction derived from no hypothesis at all,
i.e. a false theorem.  So the adelic bound of `ADELIC_HOLONOMY.md` §2 already
constrains the transported inputs: the pure inventory ($7$ functions of slope $2$
with $u=(1,3)$) together with the ceiling contour is within $0.74$ nats of being
**inconsistent**.  Either the pure inventory is smaller than CDT's on this host,
or the ceiling is not approachable, or both.  This sharpens the qualitative
"self-contradictory $d=7$ rows" observation of `ADELIC_HOLONOMY.md` §4.3 into a
number.

### 5.3 Why the conditional function is a liability

Adding one conditional function to the best unconditional inventory moves the
margin from $-0.737$ to $-2.006$ (P1) or $-5.795$ (P2): a cost of $1.27$ resp.
$5.06$ nats.  The mechanism is structural:

* every conditional function has $2$-adic slope $0$ in $x$ (`ADELIC_HOLONOMY.md`
  §3: $\xi_2=\tfrac12\zeta_2(2)\ne\xi_\infty$ as numbers, and $\zeta_2(2)$ is
  irrational by Calegari 2005), hence $\varsigma_y=2\varsigma_x-2=-2$ by the
  level-16 descent rule — it is the **worst** member of the inventory, and
  $\varsigma_{\min}$ enters $\gamma_p$ with the large weight $1-1/m$;
* on this host $\gamma_2$ is the *only* thing that keeps the entry test open at
  all (archimedean-only entry is $-0.0766$);
* the outer doubly-small orbit is also dilutive ($\varsigma_y=0$ against the
  inner orbit's $1$ and the pure module's $2$): the full $m=28$ row is worse than
  the $m=14$ unconditional row on both profiles.

> **Consequence.** On the level-16 Catalan host the CDT architecture cannot close,
> and the obstruction is no longer just the size of the deficit: the best
> configuration that could in principle prove anything is $-2.0$ (assumed slopes)
> or $-5.8$ (measured slopes), and every device that raises $m$ lowers $\gamma_2$
> faster.  The only remaining levers are (i) a genuinely denominator-free pure
> inventory (which raises $\tau^\flat$'s $u_j$, not $m$), and (ii) a better
> contour.

### 5.4 What did **not** change

The archimedean geometry: ceiling $\log 64=4.158883$, identical to level 8;
extra singularity at relative radius $\tfrac12$; the transported contour loss
$0.62922$ still not defensible; the rigorous thrice-punctured ceiling still
$1.476$, i.e. entry $\approx-2.7$ and no bound.  Nothing in this task bears on
`CATALAN_TWO_CLASSES.md` §5.3, which remains **the** decisive open input.

---

## 6. Ledger

**[proved]** Proposition 1 (§1.2, the two-part fold-regularity criterion and the
constant-term obstruction); the outer obstruction theorem (§3); Proposition 2
(§4.2, one-period $=$ three-period on this host); the character-theoretic
determination of the period set $\{1,\zeta(2),G\}$ (§4.2).

**[verified]** The fold cusp $\mathfrak a_1=\tfrac12$ and $\gamma_1$ (§1.1, §1.3);
the whole measured $2$-adic slope profile in $y$, including $\varsigma_y=-2$ for
the host and for the conditional functions (§5.0);
all six period polynomials at all five non-$\infty$ cusps, to $130$ digits, with
$(\alpha,\beta,a_0)$ recognised exactly (§1.3–1.4); the inner cross-validation
$\ell=8c_2+c_4$ and $\xi=-\tfrac12P(2)G$; the outer functionals $Q$ and
$\xi=\tfrac32\zeta(2)(P(1)-P(2))$, matching the archived measurements to $29$
digits; the existence, radius, denominators, and $2$-adic slopes of
$\Phi_0^{\rm out}$ (§2.2); the $\mathbf Q(y)$-ranks (§4.3); every margin in §5.

**[estimate]** The pure inventory ($m_{\rm pure}=7$, $u=(1,3)$, $\mathbf e$
profile) is CDT's, transported — nothing on this host realises it, and §4.3 shows
the obvious candidate set has only six independent members.  The $\times7$ orbit
multiplier is CDT's, transported, but is *verified* here for the level-16 small
functions (rank $21/21$).  The "realised contour" rows use CDT's $0.62922$, which
`CATALAN_TWO_CLASSES.md` §5.3 already showed is not transportable.

**[open]**
1. The contour (`CATALAN_TWO_CLASSES.md` §8 item 1) — unchanged and still decisive.
2. A pure module on this host, with its denominator type.  §4.3's $27/28$ says
   the naive polylogarithm candidates are one short.  **[open]**
3. Whether the $21$ small functions are independent at $\deg\le3$ (the rank came
   out $83/84$ at the edge of the series length).  **[open, minor]**
4. ~~Whether `ADELIC_HOLONOMY.md` §4.2's level-8 table has the same defect.~~
   **Done** (`20_level8_slopes.gp`, `21_level8_margin.py`): it does.  The
   level-8 conditional orbit measures $\varsigma_y=-2$, its $m=14$ margin is
   $-22.423$ at the ceiling and $-28.909$ on the transported contour, and its
   entry test fails by $0.855$.  Recorded as §7 of `ADELIC_HOLONOMY.md`, dated,
   with the superseded rows kept.  **[closed]**
5. The general shape of $\xi_\infty$ at a non-zero cusp — §1.4 exhibits
   $\tfrac32\zeta(2)(P(1)-P(2))$ on one host and one cusp; the general
   "Theorem B\* at an arbitrary cusp" is not formulated.  **[open]**
6. Whether $\Phi_0^{\rm out}$ has a level-8 shadow (it cannot: level-8 sources are
   ramified at the fold, `CATALAN_TWO_CLASSES.md` §1.3), and whether the
   two doubly-small functions plus the host span a rank-$4$ $\mathbf Q(y)$-module
   with an interesting Frobenius structure.  **[open]**

**No irrationality claim.**  The one new unconditional object is
$B_{\Phi_0^{\rm out}}\in\mathbf Q[[x]]$ of §2: radius $\sqrt2/4$ against the host's
$1/4$, $2$-adic slope $1$ in $x$ and $0$ in $y$, denominator type $d_n^2$ (sharp).
It doubles the doubly-small orbit and, by §5.3, makes the margin worse.
