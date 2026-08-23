# Catalan from the dilogarithm module at fourth roots of unity

*Fable, 2026-08-23.  Scripts: `lattice/catalan_mu4/` (`calib.py`, `mu4_series.py`,
`mu4_geom.py`, `mu4_bc.py`, `mu4_tau.py`, `mu4_indep.py`, `mu4_table.py`).
Sources: CDT, "The linear independence of $1,\zeta(2),L(2,\chi_{-3})$"
(`papers/cdt/cdt2/L2chi.tex`) — §7 (pure functions on $\mathbf P^1\setminus\{0,1,\infty\}$),
§9 (products of two logarithms), Remark `BCboundK` (number fields);
our `CDT_FINDER.md`, `CDT_UNPACKED.md`, `ADELIC_HOLONOMY.md` §2.6,
`CATALAN_OBSTRUCTION.md`.*

**No theorem is claimed.**  §§1–4 are a construction (the analytic statements are
proved; the arithmetic ones are measured to $n=160$).  §§5–8 are a *pricing*: the
conformal numbers are exact, the $\tau$'s are exact rationals from CDT's own
functional, the Bost–Charles numerators are computed numerically for explicit
concentric contours.  §9 is the verdict.

---

## 0. Verdict in one line

The polylogarithm module at $\mu_4$ gives a genuine CDT architecture for Catalan —
a pure module of **ten** $\mathbb Q(x)$-independent functions with unusually good
denominators, a fold at $x=\pm1$, and a conditional function that is fold-regular
under $G\in\mathbb Q$ — but it **fails the entry test by $-0.56$ nats at the hard
conformal ceiling**, because the Catalan class forces a degree-$2$ cover of
$\mathbf P^1\setminus\{0,1,\infty\}$ ramified at the cusp, and the $\lambda$-template's
conformal size on a degree-$d$ such cover is $16^{1/d}$, not $16$.  That is
**worse than the level-8 modular host** (entry $-0.077$ symmetrised) by $0.49$
nats and **equal to it** unsymmetrised ($-0.55$).  The loss is not repairable by
symmetrising ($x\mapsto-x$ is the only available involution and it is a net
$-0.56$), by going up the tower ($x^4$: ceiling $2$, entry $\approx-1.25$), or by
passing to $\mathbb Q(i)$ (exactly neutral, §8).

---

## 1. Calibration (`calib.py`) — all three CDT examples reproduced exactly

| example | $\tau^\flat$ | $\tau^\sharp$ | $\tau$ | CDT |
|---|---|---|---|---|
| Thm A, $1,\zeta(2),L(2,\chi_{-3})$, $m=14$ | $191/49$ | $27/80$ | $16603/3920=4.235459$ | ✅ identical |
| Thm `logsmain`, $y$-line, $m=17$ | $1136/289$ | $78419/242760$ | $1032659/242760=4.253827$ | ✅ identical |
| §7 remark, $x$-line $n[1..n]^2$, $m=17$ | $558/289$ | $83711/242760$ | $552431/242760=2.275626$ | ✅ identical |

For the third: $\log|\varphi'(0)|=\log16+\log0.9163768=2.685261$, entry $+0.409635$,
and CDT's quotient $22.7527$ forces $\mathrm{BC}=9.3203$.  The Bost–Charles
integrator of `mu4_bc.py` reproduces $\mathrm{BC}(\rho z)=\log\rho$ to $8$ decimals,
including CDT's own $\rho=161.081157$ $\Rightarrow$ $5.08190114$.  **[verified]**

Signed margin convention throughout (`CDT_FINDER.md`):
$$\textbf{margin}:=m\bigl(\log|\varphi'(0)|-\tau(\mathbf b;\mathbf e)\bigr)-\mathrm{BC}(\varphi),
\qquad \textbf{entry}:=\log|\varphi'(0)|-\tau(\mathbf b;\mathbf e).$$

---

## 2. The Catalan function and its host

Set
$$f(x)\ :=\ \sum_{k\ge0}\frac{(-1)^kx^{2k+1}}{(2k+1)^2}\ =\ \frac{\operatorname{Li}_2(ix)-\operatorname{Li}_2(-ix)}{2i}\ =\ \operatorname{Ti}_2(x)\ \in\ \mathbb Q[\![x]\!].$$

**Facts.**  **[proved / verified numerically]**

1. $f$ is holonomic on $\mathbf P^1\setminus\{0,i,-i,\infty\}$, with $\theta^2f=x/(1+x^2)$
   ($\theta=x\,d/dx$); the annihilating operator $\theta\circ\frac{1+x^2}{x}\circ\theta^2$
   has order $3$ with solution space $\{1,\log x,f\}$.
2. **$x=\pm1$ are regular points of $f$**, and $f(1)=G$, $f(-1)=-G$.
3. At the singular points, $f(\pm i)=\pm i\pi^2/8$ (checked to 20 digits).
4. Denominator type **$n^2$ exactly**: $n^2f_n\in\{0,\pm1\}$.  This is the single best
   arithmetic feature of the architecture: the Catalan function is $\mathrm{lcm}$-free.

So the geometry is *inverted* relative to a modular Apéry row: $G$ is a value at a
**regular** point of the generating function, and the singular points carry $\pi^2$.
The fold has to be created, not found.

### The pure module (`mu4_series.py`, measured to $n=159$)

Ten $\mathbb Q(x)$-linearly independent $G$-functions on $\mathbf P^1\setminus\{0,\pm i,\infty\}$,
holomorphic at $0$, with the measured minimal types $n^{e}[1..b_1n][1..b_2n]$:

| $i$ | function | $e$ | $(b_1,b_2)$ | $\sigma_i$ |
|---|---|---|---|---|
| 1 | $1$ | 0 | — | 0 |
| 2 | $\arctan x=\sum(-1)^kx^{2k+1}/(2k+1)$ | 1 | — | 0 |
| 3 | $\log(1+x^2)$ | 1 | — | 0 |
| 4 | $\operatorname{Li}_2(-x^2)$ | 2 | — | 0 |
| 5 | $f(x)=\operatorname{Ti}_2(x)$ | 2 | — | 0 |
| 6 | $\log^2(1+x^2)$ | 1 | $(\tfrac12,0)$ | $\tfrac12$ |
| 7 | $M(x)=\int_0^x\frac{\log(1+t^2)}{1+t^2}dt=\sum(-1)^{m-1}H_m\frac{x^{2m+1}}{2m+1}$ | 1 | $(\tfrac12,0)$ | $\tfrac12$ |
| 8 | $\arctan^2x$ | 1 | $(1,0)$ | 1 |
| 9 | $\arctan x\cdot\log(1+x^2)$ | 1 | $(1,0)$ | 1 |
| 10 | $J(-x^2)$, CDT's fifth function $\frac1{\sqrt{1+u}}\int_0^u\frac{\log(1+t)}{t\sqrt{1+t}}dt$ at $u=x^2$ | 0 | $(\tfrac12,\tfrac12)$ | 1 |

Rows 1–9 are the $\mathbb Q$-form of the weight-$\le2$ multiple polylogarithm ring on
the $4$-punctured sphere: $\dim=1+2+6=9$, the six weight-$2$ words being three
products of logarithms, $\operatorname{Li}_2(\pm ix)$ and the antisymmetric
$\omega_{+}\omega_{-}-\omega_{-}\omega_{+}$ (a representative of which is $M$).
Row 10 is the $\sqrt{1+x^2}$-twisted "mystery" function of CDT §7 (their
Remark `mystery three`), transported.  **[measured]** $\mathbb Q(x)$-independence:
rank $10(D+1)$ for $\deg P_i\le D$, $D\le5$, mod $2^{61}-1$, series order $200$
(`mu4_indep.py`).  **[verified numerically, not proved]**

**Five of the ten carry no $\mathrm{lcm}$ layer at all.**  CDT have one ($B_1=1$).
This is the architecture's one real arithmetic advantage, and it is worth
$\approx0.33$ nats of $\tau$ (§5).

---

## 3. The conditional function: where the hypothesis enters

Put
$$H(x)\ :=\ \int_0^x\frac{2f(t)}{1-t^2}\,dt\ =\ \int_0^{x^2}\frac{g(v)}{1-v}\,dv,\qquad
g(u):=\sum_{n\ge0}\frac{(-1)^nu^n}{(2n+1)^2},\ \ g(1)=G,$$
and $A(x):=-\log(1-x^2)$.

**Fold mechanism.**  **[proved]**  Near $t=1$ the integrand is
$\frac{2f(t)}{1-t^2}=\frac{f(1)}{1-t}+O(1)$, so
$$H(x)=-G\log(1-x^2)+(\text{holomorphic at }x=\pm1),$$
i.e. $\Delta_{x=1}H=-2\pi i\,G$ while $\Delta_{x=1}A=-2\pi i$.  Hence
$$\boxed{\ H-G\!\cdot\!A\ \text{is holomorphic at }x=\pm1\ }$$
and under the hypothesis $G=a/b$ the series
$$H_{\mathrm{cond}}:=b\,H-a\,A\ \in\ \mathbb Q[\![x]\!]$$
has rational coefficients **and** is fold-regular.  This is exactly the mechanism of
CDT §9 ($H_B-\tfrac12L\,H_A$ overconvergent past $x=1/9$) and of their Beukers-type
$H(a,x)-\eta_aA(a,x)$: the *hypothesised period is the period of the extension class
at the fold*.

**[verified numerically]** With $H=\sum_nH_nx^{2n}$ one measures
$$H_n-\frac G n\ =\ -\frac{1}{8n^3}\bigl(1+o(1)\bigr)\quad(\text{sign alternating}),$$
summable, so $\sum(H_n-G/n)$ converges (partial sum $0.0735539721$ at $n<160$),
whereas $\sum H_n$ diverges like $G\log n$.

**Denominators.**  **[measured, $n\le159$]** $H$ has type $n\,[1..n]^2$ on the
$x$-line — CDT's own conditional type.  $A$ has type $n$.

**How much hypothesis is needed.**  **[proved]**  The standard $\operatorname{Li}_2$
monodromy $\Delta_{z=1}\operatorname{Li}_2(z)=-2\pi i\log z$ gives, for the branch of
$g$ continued once around $u=-1$,
$$(\Delta_{-1}g)(1)\ =\ \pi\log(-i)\ =\ -\,\tfrac{i\pi^2}{2}$$
(checked numerically).  Therefore the **conjugates** of the loop around the fold act on
$H_{\mathrm{cond}}$ by $-2\pi i\,b\,\kappa\,i\pi^2\in\mathbb Q^\times\pi^3$, non-zero.
Consequently:

* under $G\in\mathbb Q$ alone, $H_{\mathrm{cond}}$ is fold-regular **on the principal
  branch only** — CDT's `stacky overconvergence` situation, with
  $\Sigma^0=\{1,-1\}$ and the requirement $\#\varphi^{-1}(\pm1)=1$;
* full descent to $\mathbf P^1\setminus\{0,\pm i,\infty\}$ needs the whole conjugacy
  class killed, i.e. the **three-period hypothesis** that $1,\pi^2,G$ are
  $\mathbb Q$-linearly dependent — CDT's own architecture, and the stronger theorem.

Either way the analytic requirement on $\varphi$ is the same as CDT's, and it is
*cheap here*: the fold preimages are far apart (§4).

**Size of the conditional orbit.**  **[measured]**  Unlike CDT's, the derivatives of
$H_{\mathrm{cond}}$ are **not new**: $H'=\dfrac{2f}{1-x^2}$ is a $\mathbb Q(x)$-multiple
of the pure function $f$, and $H''\in\mathbb Q(x)\{f,\arctan x\}$.  `mu4_indep.py`
finds exactly $2$ relations of degree $3$ in
$\{\text{pure}_{10},H,H',H''\}$ and none in the pure module alone.  The conditional
orbit is therefore $\{H_{\mathrm{cond}}\}$ together with *added integrations*
$\int x^jH\,dx$; the first few are independent (rank $m(D+1)$ up to $D=1$; a
defect $1$ appears at $D=2$ once three integrations are used).  Realistically
$m\approx12$–$15$, not CDT's $14$ from $7+7$.  **This turns out not to matter: the
architecture fails at entry, and entry is monotone decreasing in $m$ here** (§6).

---

## 4. The three hosts and their conformal ceilings (`mu4_geom.py`)

The Catalan class lives on the tower of covers of $\mathbf P^1\setminus\{0,1,\infty\}$,
all totally ramified over the two cusps $0,\infty$:

$$\underbrace{\mathbf P^1\setminus\{0,\pm1,\pm i,\infty\}}_{x^4}\ \longrightarrow\
\underbrace{\mathbf P^1\setminus\{0,\pm i,\infty\}}_{x^2=u,\ \textbf{host A}}\ \longrightarrow\
\underbrace{\mathbf P^1_u\setminus\{0,-1,\infty\}}_{\textbf{host B}}
\ \stackrel{w(u)=\frac{-u}{u+1}}{\longrightarrow}\
\underbrace{\mathbf P^1_v,\ v=\tfrac{u^2}{u+1}}_{\textbf{host C}}$$

**The $\lambda$-template on a degree-$d$ cover.**  **[proved]**  CDT's template
$\varphi:(\mathbb D,0)\to\mathbb C\setminus\{\text{outer}\}$ with $\varphi^{-1}(0)=\{0\}$ is
the quotient of the universal cover by the parabolic at the cusp $x=0$.  On the
degree-$d$ cover the cusp width multiplies by $d$, so the local coordinate is
$\zeta=z^{1/d}$ and
$$x=\bigl(\lambda(\zeta^{d})\bigr)^{1/d}=16^{1/d}\,\zeta+O(\zeta^{d+1}),
\qquad\boxed{\ |\varphi'(0)|_{\max}=16^{1/d}\ }$$
Verified numerically for $d=2$: $|\varphi_4(\zeta)/\zeta|\to 4.000000000$
($\varphi_4(\zeta)=i\sqrt{\lambda(\zeta^2)}$).

| host | coordinate | outer singularity | $\varphi$ ceiling | $\log$ |
|---|---|---|---|---|
| $x^4$ | $x$, $\mu_4$-full | $\pm1,\pm i$ | $16^{1/4}=2$ | $0.6931$ |
| **A** | $x$ | $\pm i$ | $16^{1/2}=4$ | $1.3863$ |
| **B** | $u=x^2$ | $-1$ | $16$ | $2.7726$ |
| **C** | $v=u^2/(u+1)$ | — ($\mathbb Z/2$ at $v=-4$) | $256$ | $5.5452$ |

(Schwarz–Pick check without the $\varphi^{-1}(0)=\{0\}$ condition: the conformal
radius of $\mathbb C\setminus\{\pm i\}$ at $0$ is $4.376879230453$, so $4$ is
consistent and the extra condition costs $\log(4.3769/4)=0.090$.)

**Symmetrisations.**  **[proved]**  The Möbius stabiliser of $\{0,\infty,i,-i\}$ that
fixes the base point $x=0$ is $\{1,\,x\mapsto-x\}$: the *only* descent available on
host A is $x\mapsto x^2$, i.e. host B.  There is no involution of host A swapping
an outer singularity with $\infty$ (the analogue of CDT's $w(x)=x/(x-1)$), because
$\{i,-i,\infty\}$ admits no such fixed-point pattern.  Host B *does* admit CDT's
involution $w(u)=-u/(u+1)$, giving host C.

**Fold preimages.**  Folds are $x=\pm1\Leftrightarrow u=1\Leftrightarrow\lambda=-1$,
i.e. the $\Gamma(2)$-orbit of $\tau=i-1$, with $\operatorname{Im}\gamma\tau=1/((d-c)^2+c^2)$:

| host | nearest preimage | second-nearest | admissible concentric radius |
|---|---|---|---|
| A ($\zeta$) | $0.20788$ (one for $x=+1$, one for $x=-1$) | $0.73040$ | $r<0.7304$ |
| B ($z$) | $0.043214$ | $0.533488$ | $r<0.5335$ |
| C ($Q$) | $0.0018674$ | $0.284648$ | $r<0.2846$ |

(For comparison, CDT's fold $y=-1/72$ has its nearest preimage at $|Q|\approx5.4\cdot10^{-5}$
and its second at $|Q|\approx0.40$; they buy $|\psi'(0)|=0.62922$ with a 4-slit + lune
Riemann map.  Our folds are *better placed than CDT's* relative to the disc, which is
why simple concentric contours already lose only $\log0.7304=-0.314$.)

---

## 5. $\tau(\mathbf b;\mathbf e)$ (`mu4_tau.py`)

Two consistent descriptions of the same module, since $n^{e}\mid[1..n]^{e}$:

* **"measured"** — use the sharp $n^{e}$ types.  Then $\max_ie_i=2$ and CDT's
  $\tau^\sharp$ functional charges for it.
* **"relaxed"** — declare every $n^{e}$ as $[1..n]^{e}$.  Then $\mathbf e=0$,
  $\tau^\sharp=0$, but the $u_j$ drop.

On host A the relaxed description is the better of the two (the $\tau^\sharp$ price of
$\max e_i=2$ exceeds the $\tau^\flat$ gain):

| host | description | $c$ | $m$ | $\sigma_m$ | $\tau^\flat$ | $\tau^\sharp$ | $\tau$ | entry @ ceiling |
|---|---|---|---|---|---|---|---|---|
| A | measured | 4 | 14 | 2 | 1.4592 | 0.7144 | 2.1735 | $-0.7872$ |
| A | **relaxed** | 4 | 14 | 2 | 1.9490 | 0 | **1.9490** | $\mathbf{-0.5627}$ |
| A | relaxed | 7 | 17 | 2 | 1.9654 | 0 | 1.9654 | $-0.5791$ |
| A | relaxed | 10 | 20 | 2 | 1.9750 | 0 | 1.9750 | $-0.5887$ |
| B | relaxed | 4 | 14 | 4 | 3.8980 | 0 | 3.8980 | $-1.1254$ |
| B | relaxed | 7 | 17 | 4 | 3.9308 | 0 | 3.9308 | $-1.1582$ |
| C | relaxed | 4 | 14 | 8 | 7.7959 | 0 | 7.7959 | $-2.2507$ |
| C | measured | 4 | 14 | 8 | 7.0816 | 0.6020 | 7.6837 | $-2.1385$ |

Two structural readings:

1. **Entry is monotone decreasing in $m$** on every one of these hosts, because every
   added conditional function carries the full $\sigma_m=2$ (resp. $4$, $8$) and
   dilutes the $u_j$.  So "find more functions" cannot fix entry here — only a better
   *fraction* of $\mathrm{lcm}$-free functions can.
2. **Each descent step costs $-0.56$.**  Doubling the denominators doubles
   $\tau^\flat$ ($\approx+1.95$) while the ceiling gains only $\log4=1.386$.  On CDT's
   own host the same descent gains $\log16=2.773$ against the same $+1.96$, i.e.
   $+0.81$ — because their involution *deletes* a puncture (swaps the outer
   singularity with $\infty$), whereas $x\mapsto-x$ merely *merges* two.
   **This is the precise reason the $\mu_4$ architecture cannot be symmetrised into
   profit.**

**How good would the pure module have to be?**  With $\sigma_m=2$ and $b_1=b_2=1$,
$\tau^\flat=2-2(u/m)^2$, so entry at the ceiling $\log4$ needs
$$\boxed{\ u_1/m=u_2/m\ >\ 0.5544\ }$$
i.e. **strictly more than 55.4% of all $m$ functions must carry no $\mathrm{lcm}$
layer at all.**  We have $5$ out of $10$ pure and $0$ out of $c$ conditional.  Since
the $\mathrm{lcm}$-free functions are exactly the depth-one $\mu_4$ polylogarithms
$\sum_{n\equiv r(4)}x^n/n^k$ — two per weight over $\mathbb Q$ on host A — raising the
fraction means raising the weight, which raises $\max_ie_i$ and hence $\tau^\sharp$
faster than it lowers $\tau^\flat$.  Even the fantasy $u_1=u_2=m/2$ gives
$\tau^\flat=1.5>\log4=1.3863$: **host A cannot pass entry with a half-free module,
at any $m$.**

---

## 6. Bost–Charles numerators and margins (`mu4_bc.py`, `mu4_table.py`)

$\mathrm{BC}(\varphi)=\iint_{\mathbb T^2}\log|\varphi(z)-\varphi(w)|$ computed as
$\iint\log\bigl|\frac{\varphi(z)-\varphi(w)}{z-w}\bigr|$ (the diagonal handled by
$\log|\varphi'|$), $N=1024$–$2048$ nodes.  Sanity: $\mathrm{BC}(\rho z)=\log\rho$ to
$8$ decimals for $\rho=0.5,3,161.081157$.

| host | contour | $\log|\varphi'(0)|$ | $\mathrm{BC}$ | $\mathrm{BC}-\log|\varphi'(0)|$ |
|---|---|---|---|---|
| A | $\varphi_4$ on $\lvert\zeta\rvert<0.4$ | $+0.47000$ | $+0.47000$ | $0$ (univalent) |
| A | $\lvert\zeta\rvert<0.60$ | $+0.87547$ | $+1.16627$ | $+0.29080$ |
| A | $\lvert\zeta\rvert<0.7304$ | $+1.07158$ | $+1.75885$ | $+0.68726$ |
| B | $-\lambda$ on $\lvert z\rvert<0.5335$ | $+2.14434$ | $+3.51776$ | $+1.37$ |
| C | $-h$ on $\lvert Q\rvert<0.2846$ | $+4.28874$ | $+7.03620$ | $+2.75$ |

Best margins over the admissible concentric contours, with the relaxed inventory:

| host | $m$ | $\tau$ | best $r$ | $\log|\varphi'(0)|$ | entry | $\mathrm{BC}$ | **margin** |
|---|---|---|---|---|---|---|---|
| **A** | 14 | 1.9490 | $0.730$ | $+1.0721$ | $\mathbf{-0.8768}$ | $1.7635$ | $\mathbf{-14.04}$ |
| A | 17 | 1.9654 | $0.730$ | $+1.0721$ | $-0.8933$ | $1.7635$ | $-16.95$ |
| B | 14 | 3.8980 | $0.533$ | $+2.1443$ | $-1.7537$ | $3.5178$ | $-28.07$ |
| C | 14 | 7.7959 | $0.285$ | $+4.2887$ | $-3.5072$ | $7.0362$ | $-56.14$ |

Entry at the *hard ceiling* (no contour loss at all) is $-0.5627$ / $-1.1254$ /
$-2.2507$ for A / B / C.  Since entry $<0$, the bound is vacuous and the margin
column is bookkeeping only.

---

## 7. The $x^4$ host is worse still

$\mathbf P^1\setminus\{0,\pm1,\pm i,\infty\}\to\mathbf P^1\setminus\{0,1,\infty\}$,
$x\mapsto x^4$, is an unramified degree-$4$ cover of the *punctured* surfaces.  The
pure module gains: the $\mathrm{lcm}$-free depth-one functions become
$\sum_{n\equiv r\,(4)}x^n/n^k$, four per weight ($\log(1-x^4)$, $\operatorname{Li}_2(x^4)$,
$\arctan$, $f$, … all reappear as $\mathbb Q$-combinations).  But $\varphi$ must now
avoid all four of $\pm1,\pm i$, and by §4 the ceiling is $16^{1/4}=2$,
$\log2=0.6931$, against $\sigma_m=2$ unchanged.  Entry needs $u/m>0.808$: out of
reach.  **The fold cannot live on this host either** — the whole point of host A is
that $\pm1$ are *not* punctures of the pure module, which is what makes the fold
mechanism free.

---

## 8. The number-field variant ($K=\mathbb Q(i)$) is exactly neutral

CDT's Remark `BCboundK`: for $f_i\in K[\![x]\!]$ with $\mathcal O_K$-numerators and
*rational-integer* denominator layers, entry reads
$\frac1{[K:\mathbb Q]}\sum_{\sigma}\log|\varphi'_\sigma(0)|>\sigma_m$ and
$$m\ \le\ \frac{\sum_\sigma\iint_{\mathbb T^2}\log|\varphi_\sigma(z)-\varphi_\sigma(w)|}
{\bigl(\sum_\sigma\log|\varphi'_\sigma(0)|\bigr)-[K:\mathbb Q]\bigl(\overline\tau(\mathbf b)+\tau^\sharp(\mathbf e)\bigr)}.$$

For $K=\mathbb Q(i)$ there is **one complex place**, i.e. **two** embeddings
$\sigma,\bar\sigma$ in the sum (the normalisation is per embedding, not per place).
Our host, its singular set $\{0,\pm i,\infty\}$ and the whole module are stable under
complex conjugation, so the optimal pair is conjugation-symmetric,
$\varphi_{\bar\sigma}=\overline{\varphi_\sigma}$, whence
$\sum_\sigma\log|\varphi'_\sigma(0)|=2\log|\varphi'(0)|$,
$\sum_\sigma\mathrm{BC}_\sigma=2\,\mathrm{BC}$, and the inequality is **identical** to
the $\mathbb Q$ one.  **[proved]**

Nor does $m$ increase: the $\mathbb Q(i)$-module is $\mathbb Q$-module $\otimes_{\mathbb Q}\mathbb Q(i)$
(the host and the local system are defined over $\mathbb Q$), so
$\dim_{K(x)}=\dim_{\mathbb Q(x)}$.  Concretely $\operatorname{Li}_2(ix)$ has
$\mathbb Z[i]$-numerators over $n^2$ — the *same* profile as its two $\mathbb Q$-forms
$\tfrac12\operatorname{Li}_2(-x^2)$ and $f$.  **Gaussian coefficients buy nothing here.**
(Contrast `CDT_FINDER.md` §3/§6, where the gain would come from a *real* quadratic
field with two real places and conjugate singular sets of different sizes.)

---

## 9. Comparison with the modular hosts, and verdict

Like-for-like, the same functional applied to Catalan's modular host (Zagier row E,
level 8, $\lambda_2=4$; pure module = the five CDT functions on
$\mathbf P^1\setminus\{0,\tfrac14,\infty\}$):

| architecture | ceiling | $\sigma_m$ | $\tau$ | **entry @ ceiling** |
|---|---|---|---|---|
| CDT's own $L(2,\chi_{-3})$, $\mathbf P^1\setminus\{0,1,\infty\}$, $n[1..n]^2$ | $16$ | 2 | $2.2756$ | $+0.4970$ |
| CDT's own $L(2,\chi_{-3})$, symmetrised $Y_0(2)$ | $256$ | 4 | $4.2355$ | $+1.3097$ |
| **$\mu_4$ host A ($x$-line), $m=14$** | $\mathbf 4$ | 2 | $1.9490$ | $\mathbf{-0.5627}$ |
| Catalan level 8, unsymmetrised | $16/4=4$ | 2 | $1.9383$ | $-0.5520$ |
| Catalan level 8, symmetrised ($Y_0(2)$) | $256/4=64$ | 4 | $4.2355$ | $-0.0766$ |
| $\mu_4$ host B ($u$-line) | $16$ | 4 | $3.8980$ | $-1.1254$ |
| $\mu_4$ host C ($v$-line) | $256$ | 8 | $7.7959$ | $-2.2507$ |
| $\mu_4$, $x^4$ host | $2$ | 2 | $\approx1.95$ | $\approx-1.25$ |

**Answer to the question asked.**  The polylogarithm-at-$\mu_4$ architecture is

* **equal** to the *unsymmetrised* modular host ($-0.563$ vs $-0.552$) — and this is
  not an accident: both have ceiling exactly $4$, one because $16/\lambda_2=16/4$,
  the other because $16^{1/2}$.  Two different mechanisms, the same number;
* **worse** than the *symmetrised* modular host by $\mathbf{0.486}$ nats, i.e. the
  best Catalan host remains level 8 at $-0.077$;
* **worse than CDT's own host by $\log4=1.386$ of conformal size**, of which the
  richer pure module (five $\mathrm{lcm}$-free functions instead of one) returns
  $0.327$, leaving the $-1.06$ that separates $+0.497$ from $-0.563$.

**Why, in one sentence.**  Catalan's period needs the character $\chi_{-4}$, hence the
fourth roots of unity, hence a degree-$2$ cover of $\mathbf P^1\setminus\{0,1,\infty\}$
totally ramified at the cusp where the $G$-function is expanded; the $\lambda$-template's
conformal size on such a cover is $16^{1/d}$, and the denominators do **not** improve
by the compensating factor, because the cover is unramified over the *outer* puncture.
The modular hosts lose the same $\log 4$ for a different reason ($|t_2|\le\tfrac14$, a
motivic invariant — `CATALAN_OBSTRUCTION.md` (2.2)), but they can recover $0.81$ of it
by an Atkin–Lehner descent that deletes a puncture, and the $\mu_4$ host has no such
involution.

**What would evade, in this world.**  Exactly one thing: a pure module on
$\mathbf P^1\setminus\{0,\pm i,\infty\}$ of type $[1..n]^2$ in which more than $55.4\%$
of the functions are $\mathrm{lcm}$-free — and §5 shows that the only source of such
functions (depth-one $\mu_4$ polylogarithms) is exhausted at two per weight and is
taxed by $\tau^\sharp$ as the weight grows.  Equivalently: one needs a *degree-1*
realisation of the $\chi_{-4}$ class, i.e. a host where $G$ is a period without the
Gaussian cover — which is `CATALAN_OBSTRUCTION.md` (E3) again.

**Not done.** No slit/gobble contour was designed for host A (only concentric discs);
no convexity-improved numerator was used; the $\mathbb Q(x)$-independence of the ten
pure functions is checked numerically to degree $5$, not proved; the exact size of the
conditional orbit (added integrations) was measured only to $D=2$; the adelic
($p$-adic) twist of `ADELIC_HOLONOMY.md` §2.6 was **not** applied — but note that the
$2$-adic obstruction of `CATALAN_OBSTRUCTION.md` (2.3) transfers verbatim: under
$G=a/b$ the conditional function's $2$-adic avatar is governed by $\xi_2$, which is
irrational, so no $\gamma_p$ can be harvested from the fold here either.
