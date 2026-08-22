# The CDT arithmetic holonomy bound on second-order (non-congruence) Apéry rows

*Claude (Opus 5), 2026-08-22. Scripts and logs: `lattice/cdt_noncongruence/`.
Calibrated on and reusing `lattice/cdt_finder/` (`CDT_FINDER.md`) and
`lattice/adelic_holonomy/` (`ADELIC_HOLONOMY.md`). Rows and scores from
`NONCONGRUENCE_SCAN.md`, `SQRT_APERY.md`, `ROOT_ROWS.md`. Literature: CDT,
arXiv:2408.15403v2 (Thm 6.0.2 / 7.0.1) and arXiv:2510.04156.
Tags: **[proved]**, **[computed]**, **[estimated]**, **[open]**.
**No irrationality claim is made anywhere in this note.***

---

## 0. Verdict

The programme was: *"the CDT holonomy bound counts functions against monodromy and is
insensitive to $\lambda$ in the way the score is; a score $-1.67$ non-congruence row may be
much closer in CDT units."* The answer is **half right, and the half that is right is
worth $+0.84$ nats, not $+1.67$.**

| claim | verdict |
|---|---|
| A geometric denominator $\lambda^n$ is handled by CDT's framework in exactly one way: rescale $t\mapsto x=t/\lambda$. Its cost is **exactly $\log\lambda$ of entry budget** and $(m-1)\log\lambda$ of margin — this *is* the $-\log\lambda$ of Theorem N1 | **[proved]** §2, and **[computed]** to machine precision as an identity with the adelic uniform-slope treatment (`geom_denom.py`) |
| For a second-order row the whole admissible inventory from one $2$-term relation is $\{1,H,\theta H\}$: $m=3$, the true dimension. $\theta^2H$ is $\mathbf Q(x)$-dependent | **[proved]** §1.3 |
| The **rigorous univalent (Koebe)** form of the bound gives a contradiction iff $\operatorname{score}>-0.71963$; the **multivalent (Kodaira) map $\varphi_r=x_2\lambda(rz)$**, with its Bost–Charles numerator computed numerically, improves this to $\operatorname{score}>\mathbf{-0.83852}$ at $r^\ast=0.46$ | **[computed]** §3, `optimum.py`, `bc_multivalent.py` |
| So **the CDT route is worth $+0.8385$ nats over the elementary Beukers criterion** on any second-order row — and *only* that much, unless the host admits CDT's symmetrisation | **[computed]** |
| **CDT's own symmetrised architecture requires $\lambda_2^{\rm norm}=|N(\lambda_2)|^{1/[K:\mathbf Q]}=1$** (the Apéry-perfect / unit-fold condition), and its threshold is $\operatorname{score}_{\rm norm}>-2.00041$ — met *exactly* at $\lambda_2^{\rm norm}=1$, margin $+0.0053$. This is the structural reason `CDT_FINDER.md` §4 found no new host | **[proved]** (the rationality of the descent) $+$ **[computed]** (the threshold) §3.3 |
| **Beukers' row is re-proved by the CDT route**, at $m=2$ (bound $1.7406<2$) and at $m=3$ (bound $2.1705<3$), deficit $+0.9777$ — but it gets **nothing** from CDT's own architecture ($\lambda_2^{\rm norm}=4$, margin $-18.02$). The two routes see different things | **[computed]** §5 |
| **The $\Gamma_0(5)+5$ non-congruence row's deficit is $-0.8349$, not $-1.6734$** — the CDT route halves it. It still fails, by a factor $2.305$ in $|t_2|$ | **[computed]** §5 |
| The entry condition for the level-$5$ row **fails** at the univalent map ($-0.0649$) and **passes** at the Kodaira map ($+0.5449$): the theorem does not even apply without the multivalent contour | **[computed]** §5 |
| $\sum_p\gamma_p=0$ for every row in the $(K)/(D)$ architectures: those inventories contain no function of positive $p$-adic slope. The adelic gain lives only in $(S)$, i.e. only on hosts with $\lambda_2\in\mathbf Z$ | **[proved]** §6 |
| **New host at CDT parity:** the square root of Cooper's $s_7$ row, $\lambda_1=27$, $\lambda_2=-1$, $k=2$ — margin $+0.0053$, tying CDT. Its period is $L(f_7,2)$, a CM weight-three critical value, so (as with the other $\lambda_2^{\rm norm}=1$ hosts) no new theorem follows | **[computed]** §5, §7 |
| **No second fold-regular class exists on $(\Gamma_0(5)+5,t)$**, and the obstruction is a dimension count: $c=\dim\{\text{usable weight-2 forms}\}\le\#\text{cusps}-1$, and $\Gamma_0(5)+5$ has **one** cusp. The nebentypus route costs $1.693$ against a gain of $0.585$; levels $10,15,20,25$ return only $E_{2,5}$ itself; the gauge family gives either $\mathbf Q(x)$-multiples, or $k'=10$, or doubled roots | **[proved]** $+$ **[computed]** §10 |
| No index-$9$/index-$10$ non-congruence host group contributes a second-order row with real $\lambda_2<1$; and the newly opened window $-0.8385<\operatorname{score}\le0$ contains, in every box the scan searched, only the two Padé families and the rescalings $(22,6,-4,0)$, $(33,9,-9,0)$ of Apéry's $\zeta(2)$ row | **[computed]** §5.4, `index910.py` |

One sentence: *the CDT bound is not insensitive to $\lambda$ — it pays exactly the same
$\log\lambda$ — but it is worth a flat $+0.84$ nats of archimedean budget over the
elementary criterion on every second-order row, and it is worth an extra $\approx1.2$ nats
on top of that only on hosts whose fold is a unit, which is precisely the class the
elementary criterion already handles.*

---

## 1. The CDT bound for a second-order Apéry row

### 1.1 The row, the hypothesis, and the conditional function

Work in the **integral coordinate** $x$: the row is
$$(n+1)^2a_{n+1}=(\alpha n^2+\beta n+\gamma)a_n-(\delta n^2+\varepsilon n+\zeta)a_{n-1},
\qquad a_0=1,\ a_1=\gamma,\ a_n\in\mathbf Z,$$
with companion $b_n$ ($b_0=0$, $b_1=1$) satisfying $d_n^2b_n\in\mathbf Z$ ($k=2$, sharp for
every row below), $\lambda_{1,2}$ the roots of $L^2-\alpha L+\delta$ ($|\lambda_1|\ge|\lambda_2|$),
$\delta=\lambda_1\lambda_2$, and singular points of the Picard–Fuchs operator at
$$x_1=1/\lambda_1,\qquad x_2=1/\lambda_2,\qquad\text{plus }0,\infty .$$
Put $A=\sum a_nx^n$, $B=\sum b_nx^n$, $\xi=\lim b_n/a_n$.

**Hypothesis to be contradicted:** $\xi=p/q\in\mathbf Q$. Set
$$H:=qB-pA=-q\sum_{n\ge0}r_nx^n,\qquad r_n:=a_n\xi-b_n .$$
Since $|r_n|^{1/n}\to|\lambda_2|$ (Poincaré; the Casoratian is non-vanishing for every row
here), $H$ has radius of convergence $|x_2|>|x_1|$: **$H$ is holomorphic at the inner
singular point $x_1$.** That is the entire content of the hypothesis — it removes $x_1$
from the branch locus of the module. Without it the module contains $A$, which branches at
$x_1$, and the domain shrinks from $|x_2|$ to $|x_1|$.

*This is Calegari's "codimension-one holomorphy at $\alpha$" and CDT's fold-regularity, for
a second-order row.*

### 1.2 Denominator type

$H$ has $d_n^2H_n\in\mathbf Z$, i.e. CDT type (6.0.9) with $r=2$ layers of rate $b_j=1$ and
$e=0$; $\theta H=xH'$ has coefficients $nH_n$, same type; the constant $1$ is
denominator-free. So $\sigma_m=2$ and the step data is $u_1=u_2=1$, $b_1=b_2=1$
— *one* denominator-free function, exactly CDT's $u_1=1$, but without their $u_2=3$.

Note the contrast with CDT's own array ($b_j=2$, $\sigma_m=4$): the doubling there is
caused by the normaliser descent, which we do **not** perform in architectures (K)/(D).

### 1.3 The inventory is $m=3$, and that is the true dimension

$H$ satisfies the second-order **inhomogeneous** equation obtained from the row's operator,
$$P_2(x)\,\theta^2H+P_1(x)\,\theta H+P_0(x)\,H=c(x)\cdot 1,\qquad P_i,c\in\mathbf Q[x],$$
so $\{1,H,\theta H,\theta^2H\}$ carries a $\mathbf Q(x)$-relation. Hence

> **the admissible inventory from one $2$-term relation is exactly $\{1,H,\theta H\}$, $m=3$**,

which is the "true dimension of the holonomic module ($=2$, or $3$ with the $1$)" of the
brief. Enlarging it requires either Eichler integrals $\theta^{-j}H$ (which cost $e_i=j$,
hence $\tau^\sharp$), or a genuinely second fold-regular class (§7). Both are priced in §3.

*Caveat.* $\mathbf Q(x)$-linear independence of $\{1,H,\theta H\}$ is a **hypothesis of
CDT's theorem** (their Lemma 12.1.1), not proved here. It amounts to $H\notin\mathbf Q(x)$
and $H'/H\notin\mathbf Q(x)$. **[open]**

### 1.4 The domain

All three functions are holomorphic at $x=0$ and their only finite branch point is $x_2$
(monodromy at $x_1$ is trivial by §1.1; there are no other singular points). Two admissible
families of $\varphi:(\overline{\mathbf D},0)\to(\mathbf C,0)$:

* **(K) univalent.** $\Omega=\mathbf C\setminus(\text{the ray from }x_2\text{ away from }0)$
  is simply connected, contains $x_1$, and every $f_i$ is single-valued holomorphic on it.
  By Koebe's $1/4$ theorem $\rho(\Omega)=4|x_2|$ **exactly**, and this is the largest
  conformal radius of any simply connected domain omitting $x_2$. For univalent $\varphi$,
  $$\mathrm{BC}(\varphi)=\iint_{\mathbf T^2}\log|\varphi(z)-\varphi(w)|\,d\mu\,d\mu=\log|\varphi'(0)|$$
  (Grunsky: $\log\frac{\varphi(z)-\varphi(w)}{z-w}=\sum b_{mn}z^mw^n$ with $b_{00}=\log\varphi'(0)$).
  **No contour has to be designed; this architecture is rigorous end-to-end.**
* **(D) multivalent (André/Kodaira).** $\varphi_r(z)=x_2\,\lambda(rz)$, $r<1$, with
  $\lambda=16w-128w^2+\cdots$ the modular $\lambda$-function in the nome $w=e^{i\pi\tau}$.
  Then $\lambda(w)=16w\,u(w)$ with $u$ non-vanishing on $\mathbf D$ and $\lambda\ne1$ on
  $\mathbf D$, so a loop in $\mathbf D$ winds **once around $0$ and zero times around $x_2$**;
  since every $f_i$ is holomorphic at $0$, every $\varphi_r^*f_i$ is single-valued
  holomorphic on $\mathbf D$. Here $|\varphi_r'(0)|=16r|x_2|\to16|x_2|$ (Kodaira's ceiling)
  but $\mathrm{BC}$ grows; the deficit
  $$\Delta(r):=\mathrm{BC}(\varphi_r)-\log|\varphi_r'(0)|$$
  is computed numerically in `bc_multivalent.py` (§3.2).

---

## 2. Geometric denominators: how CDT's framework absorbs $\lambda^n$

The modular rows of `NONCONGRUENCE_SCAN.md` are $a_n=\lambda^n[t^n]\sqrt F$ with
$\lambda\in\{1,2,4\}$: in the **modular** coordinate $t$ the coefficients $c_n=[t^n]\sqrt F$
carry unbounded powers of $2$. A geometric denominator $\lambda^{-n}$ is **not** of CDT's LCM
type (6.0.9), so $\sqrt F(t)$ is inadmissible as an $f_i$ in the coordinate $t$. There are
exactly two ways to fix it, and they agree.

> **Lemma G [proved].** Let $f(t)=\sum c_nt^n$ with $\lambda^nc_n$ of CDT type, $\lambda\in\mathbf Z_{>0}$.
>
> **(A) Rescale.** Put $x=t/\lambda$; then $F(x):=f(\lambda x)$ *is* of CDT type. A map
> $\varphi_t:\mathbf D\to\Omega_t$ is the same as $\varphi_x=\varphi_t/\lambda$, and both
> $\log|\varphi'(0)|$ and $\mathrm{BC}$ are homogeneous of degree $1$ in $\varphi$:
> $$\log|\varphi_x'(0)|=\log|\varphi_t'(0)|-\log\lambda,\qquad \mathrm{BC}(\varphi_x)=\mathrm{BC}(\varphi_t)-\log\lambda,$$
> while $\tau(\mathbf b;\mathbf e)$ is unchanged. Hence
> $$\textbf{entry}\ \rightsquigarrow\ \textbf{entry}-\log\lambda,\qquad
> \textbf{margin}\ \rightsquigarrow\ \textbf{margin}-(m-1)\log\lambda .$$
>
> **(B) Keep $t$, declare a negative $p$-adic slope.** In the coordinate $t$ every $f_i$ has
> $v_p(c_n)\ge -v_p(\lambda)n$: the **uniform** slope profile $\varsigma_i\equiv-v_p(\lambda)$.
> The adelic bound (`ADELIC_HOLONOMY.md` Thm 2.2) then gives
> $\gamma_p=-v_p(\lambda)\log p\,(1-\tfrac1m)$, i.e.
> $\tau_{\rm ad}=\tau+(1-\tfrac1m)\log\lambda$, hence the same
> $\textbf{margin}-(m-1)\log\lambda$.
>
> **(A) and (B) give identical margins.** Verified to $0$ (exactly, $\le5\cdot10^{-15}$)
> for $\lambda\in\{2,4,8\}$ and $m\in\{2,3,5,14,30\}$ in `geom_denom.log`.

**Consequence — this *is* N1's $-\log\lambda$.** Per function, and in the entry condition,
a geometric denominator costs exactly $\log\lambda$; in the margin it costs
$(1-\tfrac1m)\log\lambda$ per function, $\to\log\lambda$ as $m\to\infty$. Writing
$\lambda_2=\lambda/|t_2|$, the entry budget $\log(C|x_2|)=\log C+\log|t_2|-\log\lambda$
carries the $-\log\lambda$ of
$$\operatorname{score}=\log|t_2|-k-\log\lambda=\log(1/|\lambda_2|)-k .$$
**So the CDT bound is *not* insensitive to $\lambda$: it pays exactly the same price.** What
it *is* insensitive to is the split of $\lambda_2$ into $\lambda$ and $|t_2|$ — only the
product matters — which is why the two routes can still differ, and by exactly the amounts
computed in §3.

*(Sanity: the scan's rows $(22,6,-4,0)$, $(33,9,-9,0)$, $(44,12,-16,0)$ are Apéry's $\zeta(2)$
row rescaled by $j=2,3,4$; their scores drop by $\log j$ and, in §5, so do their CDT
deficits — by exactly $\log j$.)*

---

## 3. Three architectures and their exact thresholds

Because everything here is scale-covariant, each architecture reduces to a single number:
the **threshold on the elementary score** at which it starts to give a contradiction. Write
$\operatorname{score}=\log(1/|\lambda_2|)-2=\log|x_2|-2$; then in (K)/(D)
$$\textbf{margin}=m\bigl(\log|\varphi'(0)|-\tau\bigr)-\mathrm{BC}=(m-1)\bigl(\operatorname{score}-\text{threshold}\bigr).$$

### 3.1 (K) univalent

$\log|\varphi'(0)|=\mathrm{BC}=\log(4|x_2|)$, so a contradiction needs
$\log(4|x_2|)>\tau\cdot\frac{m}{m-1}$.

| $m$ | inventory | $\tau^\flat$ | $\tau^\sharp$ | $\tau$ | threshold on score |
|---|---|---|---|---|---|
| $2$ | $\{1,H\}$ | $3/2$ | $0$ | $1.50000$ | $-0.38629$ |
| **3** | $\{1,H,\theta H\}$ | $16/9$ | $0$ | $\mathbf{1.77778}$ | $\mathbf{-0.71963}$ |
| $4$ | $+\ \theta^{-1}H$ | $30/16$ | $0.43750$ | $2.31250$ | $-0.30296$ |
| $5$ | $+\ \theta^{-2}H$ | $1.92000$ | $0.91000$ | $2.83000$ | $+0.15121$ |

Eichler integrals make it **worse**: $\tau^\sharp$ grows faster than $\tau^\flat$ falls. The
optimum is $m=3$, i.e. exactly the true dimension of the module.

### 3.2 (D) multivalent Kodaira map — the Bost–Charles numerator, computed

$\Delta(r)=\mathrm{BC}(\varphi_r)-\log|\varphi_r'(0)|$, computed by removing the diagonal
singularity analytically ($\iint_{\mathbf T^2}\log|z-w|=0$) and trapezoidal quadrature on
$N$ points per circle; converged to $\approx5$ decimals ($N=128,192,256,384$ agree to
$5\cdot10^{-6}$). $\Delta(r)=0$ for $r\le0.20$ — a genuine check, since $\lambda$ is
univalent there and (K)'s Grunsky identity must hold.

| $r$ | $0.20$ | $0.25$ | $0.30$ | $0.40$ | $0.46$ | $0.53$ | $0.60$ | $0.70$ | $0.80$ |
|---|---|---|---|---|---|---|---|---|---|
| $\Delta(r)$ | $0.00000$ | $0.11705$ | $0.32225$ | $0.74469$ | $0.98175$ | $1.35369$ | $1.75586$ | $2.43881$ | $3.36078$ |

Margin $=(m-1)\log(4|x_2|)-m\tau+G(m)$ with
$G(m)=\max\bigl(0,\max_r[(m-1)\log(4r)-\Delta(r)]\bigr)$ the gain over Koebe:

| $m$ | $2$ | **3** | $4$ | $5$ | $7$ | $9$ | $14$ |
|---|---|---|---|---|---|---|---|
| $r^\ast$ | — | $\mathbf{0.46}$ | $0.53$ | $0.62$ | $0.73$ | $0.78$ | $0.80$ |
| $G(m)$ | $0$ | $\mathbf{0.23779}$ | $0.90055$ | $1.74975$ | $3.74664$ | $5.95222$ | $11.76018$ |

With the (§1.3) inventory the optimum is again $m=3$:
$$\boxed{\ \textbf{(D) threshold}\ =\ \frac{3\cdot\frac{16}9-0.23779}{2}-\log4-2\ =\ \mathbf{-0.83852}\ }$$
(with a pure polylogarithm module, available only when $\lambda_2\in\mathbf Z$, $m=4$ gives
$-1.10315$; see §3.3).

### 3.3 (S) CDT's own symmetrised architecture — and when it is available

CDT descend along the normaliser involution $w(x)=\frac{sx}{x-s}$, $y=x+w(x)=\frac{x^2}{x-s}$,
$s=1/\lambda_2$: the ceiling rises from $16|s|$ to $256|s|$, the denominators double
($b_j:1\to2$, $\sigma_m:2\to4$), and one gets $m=14$ with $u=(1,3)$, $\tau=16603/3920=4.235459$.

> **Proposition S [proved].** The descent is defined over $\mathbf Q$ iff $s=1/\lambda_2\in\mathbf Q$;
> and integrality of the pure module $\mathrm{Li}_j(x/s)=\sum\lambda_2^nx^n/n^j$ then forces
> $\lambda_2\in\mathbf Z$, hence $|\lambda_2|\ge1$ and $\operatorname{score}\le-2$. Over
> $K=\mathbf Q(\lambda_2)$ the normalisation of `CDT_FINDER.md` §3 replaces $|\lambda_2|$ by
> $$\lambda_2^{\rm norm}=|N_{K/\mathbf Q}(\lambda_2)|^{1/[K:\mathbf Q]}=|\delta|^{1/2}\ \ (\text{quadratic case}),$$
> and $\lambda_2^{\rm norm}\ge1$ always, since $\delta=\lambda_1\lambda_2\in\mathbf Z\setminus\{0\}$.

**[computed]** The threshold in (S), with CDT's realised contour loss $0.62922$ and
$\mathrm{BC}=11.845+\log s$:
$$\text{margin}=13\log s+14\bigl(\log256+\log0.62922-\tau\bigr)-11.845
=13\bigl(\operatorname{score}_{\rm norm}+2.00041\bigr),$$
$$\boxed{\ \textbf{(S) threshold}=-2.00041\ \ (\text{hard ceiling: }-2.49931),\qquad
\operatorname{score}_{\rm norm}:=-\log\lambda_2^{\rm norm}-2\ }$$
i.e. **(S) works iff $\lambda_2^{\rm norm}\le1.00041$, i.e. iff $\lambda_2^{\rm norm}=1$** —
the Apéry-perfect / Fricke-palindromic condition $|N(\lambda_2)|=1$, met with margin
$+0.0053$. *This is the structural explanation of `CDT_FINDER.md` §4: CDT sit exactly on
their own threshold, and there is no room at all for a host with $\lambda_2^{\rm norm}\ge2$.*

Asymptotically (many independent functions, $u_j=1$) the two architectures reach
$$\text{(D)}:\ \operatorname{score}>2-\log16-2=-2.77259,\qquad
\text{(S)}:\ \operatorname{score}_{\rm norm}>4-\log256-2=-3.54518 .$$

### 3.4 Summary of thresholds

| architecture | requires | threshold | worth over elementary |
|---|---|---|---|
| elementary (Beukers) | — | $\operatorname{score}>0$ | — |
| (K) univalent, $m=2$ | — | $>-0.38629$ | $+0.386$ |
| **(K) univalent, $m=3$** | — | $>-0.71963$ | $+0.720$ |
| **(D) Kodaira $\varphi_{0.46}$, $m=3$** | — | $\mathbf{>-0.83852}$ | $\mathbf{+0.839}$ |
| (D$+$) $+$ pure module, $m=4$ | $\lambda_2\in\mathbf Z$ | $>-1.10315$ | $+1.103$ |
| (S) CDT symmetrised, $m=14$ | $\lambda_2^{\rm norm}=1$ | $\operatorname{score}_{\rm norm}>-2.00041$ | — |
| (D), $c\to\infty$ fold-regular classes | $c$ classes | $\to-2.77259$ | $+2.773$ |
| (S), $m\to\infty$ | $\lambda_2\in\mathbf Q$ | $\operatorname{score}_{\rm norm}\to-3.54518$ | — |

---

## 4. The hosts (exact data, re-verified)

`verify_rows.gp` rebuilds every row to $n=160$ and checks $a_n\in\mathbf Z$, $k=2$ sharp
($d_nb_n\notin\mathbf Z$), the characteristic roots, $|a_n\xi-b_n|^{1/n}\to|\lambda_2|$, and
$\xi$ to $30$ digits (agreeing with `NONCONGRUENCE_SCAN.md` and `SQRT_APERY.md` in every case).

| host | $(\alpha,\gamma,\delta,\zeta)$ | $\lambda_1$ | $\lambda_2$ | $\delta$ | $\lambda$ | $\lambda_2^{\rm norm}$ | $|x_2|$ | score | period |
|---|---|---|---|---|---|---|---|---|---|
| **Beukers $\sqrt{\text{Apéry}}$** $\Gamma_0(6)+6$ | $(136,10,16,4)$ | $4(17{+}12\sqrt2)$ | $4(17{-}12\sqrt2)$ | $16$ | $4$ | $4$ | $8.49264$ | $\mathbf{+0.1392}$ | $L(\Psi,2)=0.10018744\ldots$ |
| **level-5 Fricke** $\Gamma_0(5)+5$ | $(88,6,-64,-12)$ | $44{+}20\sqrt5$ | $44{-}20\sqrt5$ | $-64$ | $2$ | $8$ | $1.38627$ | $-1.6734$ | $0.164306701064\ldots$ **unid.** |
| $\sqrt T$ (level 8) | $(24,2,16,4)$ | $12{+}8\sqrt2$ | $12{-}8\sqrt2$ | $16$ | $1$ | $4$ | $1.45711$ | $-1.6235$ | $L(f_8,2)$ (CM) |
| $\sqrt{\text{Domb}}$ (level 6) | $(20,2,64,16)$ | $16$ | $4$ | $64$ | $1$ | $4$ | $0.25$ | $-3.3863$ | $L(f_{12},2)$ (CM) |
| $\sqrt{s_7}$ (level 7) | $(26,2,-27,-6)$ | $27$ | $-1$ | $-27$ | $1$ | $\mathbf1$ | $1$ | $-2.0000$ | $L(f_7,2)$ (CM) |
| $\sqrt{s_{10}}$ (level 10) | $(24,2,-256,-60)$ | $32$ | $-8$ | $-256$ | $2$ | $8$ | $0.125$ | $-4.0794$ | $0.31692536\ldots$ **unid.** |
| $\sqrt{s_{18}}$ (level 18) | $(56,6,768,180)$ | $32$ | $24$ | $768$ | $2$ | $24$ | $0.04167$ | $-5.1781$ | $0.48423755\ldots$ **unid.** |
| $\sqrt{\mathrm{AZ}(9,3,-27)}$ | $(72,6,-432,-108)$ | $36{+}24\sqrt3$ | $36{-}24\sqrt3$ | $-432$ | $4$ | $20.785$ | $0.17956$ | $-3.7173$ | $0.14551448\ldots$ **unid.** |
| *calib.* Apéry $\zeta(2)$ $\Gamma_1(5)$ | $(11,3,-1,0)$ | $\tfrac{11+5\sqrt5}2$ | $\tfrac{11-5\sqrt5}2$ | $-1$ | $1$ | $\mathbf1$ | $11.09017$ | $+0.4061$ | $\zeta(2)/5$ **proved** |
| *calib.* Zagier C $=$ CDT | $(10,3,9,0)$ | $9$ | $1$ | $9$ | $1$ | $\mathbf1$ | $1$ | $-2.0000$ | $L(2,\chi_{-3})/2$ **proved (CDT)** |
| *calib.* Zagier E $=$ Catalan | $(12,4,32,0)$ | $8$ | $4$ | $32$ | $1$ | $4$ | $0.25$ | $-3.3863$ | $G$ **[open]** |

Two structural remarks that the table makes visible.

1. **$\lambda_2^{\rm norm}=1$ happens in two ways.** Either $\lambda_2=\pm1$ rational
   (then $\operatorname{score}=-2$ exactly: $\sqrt{s_7}$, CDT's own host), or $\lambda_2$ is
   a unit of a real quadratic order, $|\delta|=1$ (then $\operatorname{score}$ can be
   positive: Apéry's $\zeta(2)$ row, $\delta=-1$). **The Apéry-perfect class is exactly the
   intersection of the two routes** — and it is the only place where both work.
2. **Beukers' row has $|\delta|=16$, the level-5 row $|\delta|=64$.** That is *why* they are
   invisible to CDT's architecture: their fold is not a unit. In `NONCONGRUENCE_SCAN.md` §6
   the same fact appears as "level $6$ wins because $t_1t_2=1$"; in CDT units it appears as
   $\lambda_2^{\rm norm}=|\delta|^{1/2}=4$ resp. $8$, i.e. budgets $\log64$ and $\log32$
   instead of $\log256$.

---

## 5. Margins

### 5.1 Entry conditions, in `CDT_FINDER.md` units

`entries.log`. $(K)$: $\log|\varphi'(0)|=\log(4|x_2|)$, $\tau=16/9$. $(D)$: $+\log(4r^\ast)=+0.60977$.
$(S)$: ceiling $\log(256/\lambda_2^{\rm norm})$, $\tau=4.235459$, entryR $=$ entryC $-0.46327$.

| host | $\log 4|x_2|$ | **(K) entry** | **(D) entry** | ceil$_S$ | **(S) entryC** | **(S) entryR** | $\gamma_2$ | (S) adC | (S) adR |
|---|---|---|---|---|---|---|---|---|---|
| Beukers | $3.5255$ | $+1.7477$ | $+2.3575$ | $4.1589$ | $-0.0766$ | $-0.5398$ | $0$ | $-0.0766$ | $-0.5398$ |
| **level-5** | $1.7129$ | $\mathbf{-0.0649}$ | $\mathbf{+0.5449}$ | $3.4657$ | $-0.7697$ | $-1.2330$ | $0$ | $-0.7697$ | $-1.2330$ |
| $\sqrt T$ | $1.7627$ | $-0.0150$ | $+0.5947$ | $4.1589$ | $-0.0766$ | $-0.5398$ | $0$ | $-0.0766$ | $-0.5398$ |
| $\sqrt{\text{Domb}}$ | $0.0000$ | $-1.7778$ | $-1.1680$ | $4.1589$ | $-0.0766$ | $-0.5398$ | $+0.2546$ | $+0.1780$ | $-0.2852$ |
| $\sqrt{s_7}$ | $1.3863$ | $-0.3915$ | $+0.2183$ | $5.5452$ | $+1.3097$ | $+0.8464$ | $0$ | $+1.3097$ | $+0.8464$ |
| $\sqrt{s_{10}}$ | $-0.6931$ | $-2.4709$ | $-1.8612$ | $3.4657$ | $-0.7697$ | $-1.2330$ | $+0.3819$ | $-0.3878$ | $-0.8511$ |
| $\sqrt{s_{18}}$ | $-1.7918$ | $-3.5695$ | $-2.9598$ | $2.3671$ | $-1.8683$ | $-2.3316$ | $+0.5837$ | $-1.2846$ | $-1.7479$ |
| $\sqrt{\mathrm{AZ}}$ | $-0.3310$ | $-2.1087$ | $-1.4990$ | $2.5110$ | $-1.7245$ | $-2.1878$ | $0$ | $-1.7245$ | $-2.1878$ |
| Apéry $\zeta(2)$ | $3.7924$ | $+2.0146$ | $+2.6243$ | $5.5452$ | $+1.3097$ | $+0.8464$ | $0$ | $+1.3097$ | $+0.8464$ |
| CDT host | $1.3863$ | $-0.3915$ | $+0.2183$ | $5.5452$ | $+1.3097$ | $+0.8464$ | $0$ | $+1.3097$ | $+0.8464$ |
| Catalan | $0.0000$ | $-1.7778$ | $-1.1680$ | $4.1589$ | $-0.0766$ | $-0.5398$ | $+0.2546$ | $+0.1780$ | $-0.2852$ |

(The Catalan row reproduces `ADELIC_HOLONOMY.md` §4.2 exactly: $-0.0766$, $-0.5398$,
$+0.1780$, $-0.2852$.)

**Note the level-5 line.** With the univalent map the entry condition **fails**
($-0.0649$): CDT's theorem does not apply at all. With the Kodaira map it **passes**
($+0.5449$). *The multivalent contour is not an optimisation here, it is the difference
between having a theorem and not having one.*

### 5.2 Margins and dimension bounds

| host | (K) $m=2$: $m\le$ | (D) $m=3$: entry / BC / $m\le$ | contradiction? |
|---|---|---|---|
| **Beukers** | $\mathbf{1.7406}<2$ ✅ | $+2.3575$ / $5.1170$ / $\mathbf{2.1705}<3$ ✅ | **yes, both** |
| Apéry $\zeta(2)$ | $1.6543<2$ ✅ | $+2.6243$ / $5.3839$ / $2.0515<3$ ✅ | yes, both |
| level-5 | $8.0452$ | $+0.5449$ / $3.3044$ / $6.0643$ | no |
| $\sqrt T$ | $6.7089$ | $+0.5947$ / $3.3543$ / $5.6399$ | no |
| $\sqrt{s_7}$ | entry fails | $+0.2183$ / $2.9778$ / $13.642$ | no (but **yes in (S)**) |
| Catalan, $\sqrt{\text{Domb}}$, $\sqrt{s_{10}}$, $\sqrt{s_{18}}$, $\sqrt{\mathrm{AZ}}$ | entry fails | entry fails | no |

### 5.3 The ranking table — deficits in score units

$\text{deficit}=\operatorname{score}-\text{threshold}$; positive $=$ contradiction. The
$(S)$ columns score the host by $\operatorname{score}_{\rm norm}=-\log\lambda_2^{\rm norm}-2$.
The margin in `CDT_FINDER.md` units is $(m-1)\times$deficit in (K)/(D) and $13\times$deficit in (S).

| host | score | $\operatorname{score}_{\rm norm}$ | (K) | **(D)** | (D$+$) | (S) real | (S) ceil | (S)$+$ad real | (S)$+$ad ceil | **best** |
|---|---|---|---|---|---|---|---|---|---|---|
| **Beukers** | $+0.1392$ | $-3.3863$ | $+0.8588$ | $\mathbf{+0.9777}$ | — | $-1.3859$ | $-0.8870$ | $-1.3859$ | $-0.8870$ | **(D) $+0.978$ ✅** |
| Apéry $\zeta(2)$ | $+0.4061$ | $-2.0000$ | $+1.1257$ | $+1.2446$ | — | $+0.0004$ | $+0.4993$ | $+0.0004$ | $+0.4993$ | (D) $+1.245$ ✅ |
| $\sqrt{s_7}$ | $-2.0000$ | $-2.0000$ | $-1.2804$ | $-1.1615$ | $-0.8969$ | $\mathbf{+0.0004}$ | $+0.4993$ | $+0.0004$ | $+0.4993$ | **(S) $+0.0004$ ✅** |
| CDT host | $-2.0000$ | $-2.0000$ | $-1.2804$ | $-1.1615$ | $-0.8969$ | $+0.0004$ | $+0.4993$ | $+0.0004$ | $+0.4993$ | (S) $+0.0004$ ✅ |
| **Catalan** | $-3.3863$ | $-3.3863$ | $-2.6667$ | $-2.5478$ | $-2.2831$ | $-1.3859$ | $-0.8870$ | $-1.1117$ | $\mathbf{-0.6128}$ | (S)$+$ad ceil $-0.613$ |
| $\sqrt{\text{Domb}}$ | $-3.3863$ | $-3.3863$ | $-2.6667$ | $-2.5478$ | $-2.2831$ | $-1.3859$ | $-0.8870$ | $-1.1117$ | $-0.6128$ | $-0.613$ |
| $\sqrt T$ | $-1.6235$ | $-3.3863$ | $-0.9039$ | $\mathbf{-0.7850}$ | — | $-1.3859$ | $-0.8870$ | $-1.3859$ | $-0.8870$ | (D) $-0.785$ |
| **level-5** | $-1.6734$ | $-4.0794$ | $-0.9538$ | $\mathbf{-0.8349}$ | — | $-2.0790$ | $-1.5801$ | $-2.0790$ | $-1.5801$ | (D) $-0.835$ |
| $\sqrt{s_{10}}$ | $-4.0794$ | $-4.0794$ | $-3.3598$ | $-3.2409$ | $-2.9763$ | $-2.0790$ | $-1.5801$ | $-1.6677$ | $-1.1688$ | $-1.169$ |
| $\sqrt{s_{18}}$ | $-5.1781$ | $-5.1781$ | $-4.4584$ | $-4.3395$ | $-4.0749$ | $-3.1776$ | $-2.6787$ | $-2.5490$ | $-2.0501$ | $-2.050$ |
| $\sqrt{\mathrm{AZ}}$ | $-3.7173$ | $-5.0342$ | $-2.9976$ | $-2.8787$ | — | $-3.0338$ | $-2.5349$ | $-3.0338$ | $-2.5349$ | $-2.535$ |

*Calibration: the (S) column reproduces `CDT_FINDER.md` exactly — margin
$13\times(+0.0004)=+0.0053$ at $\lambda_2^{\rm norm}=1$, $13\times(-1.3859)=-18.02$ at
$\lambda_2^{\rm norm}=4$, $13\times(-2.0790)=-27.03$ at $8$.*

### 5.4 Index-9/10 groups, and the window the CDT route opens

**[computed]** (`index910.py`). By `NONCONGRUENCE_SCAN.md` Thm N3 and §3, the twelve
non-congruence four-special-point groups of index $7,9,10$ contribute, at best, the
$e=\infty$ (two-cusp, "Zagier") recurrence class, which Scan A searched exhaustively.
Re-reading every scored output for rows with **real** $|\lambda_2|<1$, $k\ge2$ and
non-degenerate Casoratian, and restricting to the modular window $\alpha\le150$,
$|\delta|\le100$ (which contains every census row), the complete list in the newly opened
window $-0.8385<\operatorname{score}\le0$ is

| score | $(\alpha,\gamma,\delta,\zeta)$ | identification |
|---|---|---|
| $+0.4061$ | $(11,3,-1,0)$ | Apéry's $\zeta(2)$ row |
| $+0.1392$ | $(136,10,16,4)$ | Beukers' row |
| $+0.0138$, $-0.1297$, $-0.2973$, $-0.4990$, $-0.7524$ | $\gamma=-\tfrac54\alpha$, $\delta=16$, $\zeta=140$ | Legendre/Padé $\log$, Hecke variant |
| $-0.2871$ | $(22,6,-4,0)$ | Apéry's $\zeta(2)$ row rescaled by $2$ |
| $-0.4579$, $-0.7927$ | $(126,-84,-27,-84)$, $(90,-60,-27,-84)$ | Padé $\arctan$ family |
| $-0.6926$ | $(33,9,-9,0)$ | Apéry's $\zeta(2)$ row rescaled by $3$ |

**No index-$9$/$10$ non-congruence host contributes anything**, and the window contains no
new modular row: it contains the two Padé families and the rescalings of Apéry's own row.
(The rescalings are a good consistency check on Lemma G: the CDT route survives $j=2,3$ —
deficits $+0.552$, $+0.146$ — and dies at $j=4$ ($-0.141$), exactly $\log j$ per step, and it
re-proves $\zeta(2)\notin\mathbf Q$ each time.)

---

## 6. The adelic term $\sum_p\gamma_p$

**[proved]** In architectures (K)/(D) the inventory is $\{1,H,\theta H\}$ and

* the constant has slope $0$;
* $H=qB-pA$ has $p$-adic slope $0$ at every $p$: $p$-adically $H\approx(q\xi_p-p)A$ with
  $\xi_p$ the $p$-adic Apéry limit, and $\xi_p\ne\xi_\infty$ *as numbers*, so the leading
  behaviour does not cancel (this is exactly the middle row of `ADELIC_HOLONOMY.md` §3, and
  it is why the hypothesis "buys archimedean smallness and buys nothing $p$-adically");
* $\theta H$ likewise.

Hence $\gamma_p=0$ and $\tau_{\rm ad}=\tau$: **the adelic refinement is worth exactly nothing
in the architectures that reach the near-miss hosts.** The gain lives only in (S), where the
pure polylogarithm orbit $\mathrm{Li}_j(\lambda_2x)$ has slope $v_p(\lambda_2)$ — which needs
$\lambda_2\in\mathbf Z$, i.e. $\operatorname{score}\le-2$. With CDT's $7+7$ split
($\gamma_p=\varsigma\log p\,(u-1)^2/m^2$, $u=7$, $m=14$):

| host | $\lambda_2$ | slopes | $\sum_p\gamma_p$ | (S) margin | $\to$ adelic |
|---|---|---|---|---|---|
| Catalan, $\sqrt{\text{Domb}}$ | $4$ | $v_2=2$ | $+0.2546$ | $-18.017$ | $-14.452$ |
| $\sqrt{s_{10}}$ | $-8$ | $v_2=3$ | $+0.3819$ | $-27.027$ | $-21.680$ |
| $\sqrt{s_{18}}$ | $24$ | $v_2=3,v_3=1$ | $+0.5837$ | $-41.309$ | $-33.137$ |
| Beukers, level-5, $\sqrt T$, $\sqrt{\mathrm{AZ}}$ | irrational | — | $0$ | — | — |
| $\sqrt{s_7}$, CDT, Apéry $\zeta(2)$ | unit | — | $0$ | $+0.005$ | $+0.005$ |

This is `ADELIC_HOLONOMY.md` §5's structural ceiling
$\log(256/\lambda_2)+\sigma_p\log p\le\log(256\lambda_1)$ again: the archimedean and $p$-adic
budgets are in direct competition, and the hosts that maximise one have zero of the other.

---

## 7. Ranking, and the single ingredient that would close each

### 7.1 Ranked by best deficit (open periods in bold)

| rank | deficit | architecture | host | period |
|---|---|---|---|---|
| — | $+1.2446$ | (D) | Apéry $\zeta(2)$ row | $\zeta(2)/5$ — already proved |
| — | $+0.9777$ | (D) | **Beukers' row** | $L(\Psi,2)$ — proved (Beukers 1987 Thm 3) |
| — | $+0.0004$ | (S) | $\sqrt{s_7}$, CDT's own host | $L(f_7,2)$ (CM), $L(2,\chi_{-3})/2$ |
| **1** | $-0.6128$ | (S)$+$adelic, ceiling | **Catalan** (Zagier E) | $\mathbf{G}$ |
| 2 | $-0.6128$ | (S)$+$adelic, ceiling | $\sqrt{\text{Domb}}$ | $L(f_{12},2)$ (CM) |
| 3 | $-0.7850$ | (D) | $\sqrt T$ | $L(f_8,2)$ (CM) |
| **4** | $-0.8349$ | (D) | **level-5 $\Gamma_0(5)+5$** | $\mathbf{0.164306701064\ldots}$ **unidentified** |
| 5 | $-1.1688$ | (S)$+$adelic, ceiling | $\sqrt{s_{10}}$ | $0.31692536\ldots$ unidentified |
| 6 | $-2.0501$ | (S)$+$adelic, ceiling | $\sqrt{s_{18}}$ | $0.48423755\ldots$ unidentified |
| 7 | $-2.5349$ | (S), ceiling | $\sqrt{\mathrm{AZ}(9,3,-27)}$ | $0.14551448\ldots$ unidentified |

**The nearest host to a CDT proof is still Catalan** ($-0.613$ at the hard ceiling with the
$2$-adic slope, $-1.112$ with CDT's realised contour), with the level-$5$ non-congruence row
**second at $-0.835$** and $\sqrt T$ between them at $-0.785$. The three CM rows
($\sqrt T$, $\sqrt{\text{Domb}}$, $\sqrt{s_7}$) carry critical values of CM weight-three
newforms, which are algebraic multiples of CM periods and therefore not open problems; the
genuinely open targets are **Catalan** and the four unidentified non-congruence periods.

### 7.2 What each ingredient is worth

**[computed]** (`ingredients.py`). The lever with by far the largest value is **more
fold-regular conditional classes**: $c$ of them give $m=2c+1$ functions with $e_i\equiv0$,
so $\tau^\flat=2-2/m^2\to2$ while $G(m)$ grows.

| $c$ (fold-regular classes) | $1$ | $2$ | $3$ | $4$ | $5$ | $\to\infty$ |
|---|---|---|---|---|---|---|
| $m=2c+1$ | $3$ | $5$ | $7$ | $9$ | $11$ | — |
| threshold on score | $-0.8385$ | $-1.4237$ | $-1.7250$ | $-1.9081$ | $-2.0316$ | $-2.7726$ |
| gain | — | $+0.585$ | $+0.887$ | $+1.070$ | $+1.193$ | $+1.934$ |

Reading this against §5.3:

* **level-5 $\Gamma_0(5)+5$** (deficit $-0.8349$): **two** further fold-regular classes
  ($c=3$) close it, with $+0.0516$ to spare. One further class leaves $-0.2497$.
* **$\sqrt T$** (deficit $-0.7850$): $c=3$ closes it with $+0.1015$; $c=2$ leaves $-0.1998$.
* **Catalan** (deficit $-0.6128$ in (S)): unreachable in (D) at any $c$ (its score
  $-3.3863$ is below the asymptote $-2.7726$); it needs the *symmetrised* architecture, and
  there the missing ingredient is the one `ADELIC_HOLONOMY.md` §4.3 already identified — a
  **second $\chi_{-4}$ class with a fold-regular companion** — plus the contour, which at
  $+0.178$ of entry is now the binding constraint.
* A **$p$-adic slope** is worth $\gamma_p\approx0.25$–$0.58$, i.e. $\tfrac13$ to $\tfrac12$
  of one extra class, and is **available only where the elementary route is already dead**.
* A **Galois trace to a number field** is a strict **loss** for every decaying row: it
  replaces $\log|x_2|$ by $-\log\lambda_2^{\rm norm}$, costing
  $$\text{Beukers } 3.5255,\qquad \text{level-5 } 2.4061,\qquad \sqrt T\ 1.7627 .$$
  Buying a pure polylogarithm module over $K$ therefore never pays on these hosts.

**Where a second fold-regular class could come from.** For a second-order row the
hypothesis $a+b\xi=0$ is $2$-term, so the conditional space is $1$-dimensional — this is the
structural handicap against CDT, whose $3$-term relation $a+b\zeta(2)+cL(2,\chi_{-3})=0$
gives $\dim=3$. The two honest routes to $\dim\ge2$ are (i) a **host carrying two periods in
one relation** (two weight-two forms $F,F'$ on the *same* $(\Gamma,t)$ — by
`NONCONGRUENCE_SCAN.md` Lemma 1.1 the singular set, hence the whole archimedean side, is
independent of $F$, so this costs nothing geometrically), and (ii) CDT's own device of an
inhomogeneous ODE with several free constants. On $\Gamma_0(5)+5$ route (i) is a finite exact
computation in $M_2(\Gamma_0(5))^{\pm}$ and its index-two multiplier covers; **it was not
done here**, and it is the single thing to compute next for the level-$5$ row.

### 7.3 Comparison with the two benchmarks

| target | project's quoted deficit | in the units of this note |
|---|---|---|
| **Catalan** | margin $-7.55$ (level 16, two classes) / $-7.97$ (level 8, adelic, ceiling) / $-14.45$ (realised contour) | deficit $-0.613$ (ceiling) / $-1.112$ (contour); $=$ margin$/13$ |
| **$L(3,\chi_5)$**, $X_1(5)\,\mathrm{Sym}^2$ | entryC $-0.435$, entryR $-0.898$, margin $-24.4$ | threshold for $k=3$ ($\sigma_m=6$, $\tau=5.980357$) is $\operatorname{score}_{\rm norm}>-0.12129$; the host has $\lambda_2^{\rm norm}=1$, $\operatorname{score}_{\rm norm}=-2$, deficit $\mathbf{-1.8787}$ ($=-24.42/13$) |
| **level-5 non-congruence** | score $-1.6734$ | deficit $\mathbf{-0.8349}$ |
| **Beukers** | score $+0.1392$ | deficit $\mathbf{+0.9777}$ ✅ |

So on the common scale: **level-5 ($-0.835$) is closer than $L(3,\chi_5)$ ($-1.879$) and
farther than Catalan ($-0.613$)**. The ordering is *not* that of the raw scores
($-1.673$, $-2.000$, $-3.386$), because the three targets live in three different
architectures whose thresholds are

| target | architecture | threshold | worth over the elementary criterion |
|---|---|---|---|
| level-5 | (D), $m=3$, $k=2$ | $-0.8385$ | $+0.839$ |
| $L(3,\chi_5)$ | (S), $m=14$, $k=3$ | $-0.1213$ | $+0.121$ |
| Catalan | (S)$+$adelic at the ceiling, $m=14$, $k=2$ | $-2.7735$ | $+2.774$ |

— i.e. the extra Eichler integration of the $L(3,\chi_5)$ host ($k=3$, $\sigma_m=6$) costs
almost the whole value of CDT's architecture, which is exactly `CDT_FINDER.md` §6's verdict
read in these units.

---

## 8. Honest ledger

**Computed exactly here.** All $\tau^\flat,\tau^\sharp$ (via `cdt_finder/cdt_bound.py`,
recalibrated to CDT's $191/49$, $27/80$, $16603/3920$, margin $+0.0053$); every row's
$a_n\in\mathbf Z$, $k=2$ sharp, $\lambda_{1,2}$, $\delta$, $|x_2|$, $\xi$ to $30$ digits
($n\le160$, `verify_rows.gp`); the exact identity of the two treatments of a geometric
denominator (§2, to $5\cdot10^{-15}$); the thresholds of §3; the ranking of §5.3 and §7.1;
the index-9/10 re-scan of §5.4.

**Computed numerically.** $\Delta(r)=\mathrm{BC}(\varphi_r)-\log|\varphi_r'(0)|$ for the
Kodaira family, by trapezoidal quadrature after analytic removal of the diagonal
singularity; converged to $\approx5\cdot10^{-6}$ over $N=128\ldots384$, and equal to $0$ for
$r\le0.20$ where $\lambda$ is univalent (an independent check against the Grunsky identity).
The optimum $r^\ast$ is on a grid of step $0.01$; $r^\ast$ is capped at $0.80$ for
$m\ge11$, so the large-$m$ entries of the §7.2 table are lower bounds.

**Cited, not recomputed.** CDT's Bost–Charles integral $11.845$ and their contour loss
$0.62922$; the rows' provenance and periods (`SQRT_APERY.md`, `NONCONGRUENCE_SCAN.md`,
`ROOT_ROWS.md`); the adelic theorem of `ADELIC_HOLONOMY.md` §2 (used only for the
consistency check of §2 and the $\gamma_p$ of §6).

**Estimated, and flagged.** (i) Architecture (S) transports CDT's contour loss and
Bost–Charles *shape* term to other hosts — `CDT_FINDER.md` §8 estimate 1, unchanged here.
(ii) Architecture (S)'s pure inventory $u=(1,3)$, $m=14$ is CDT's, transported; nothing here
identifies pure functions on any other host. (iii) The number-field normalisation of
Proposition S is `CDT_FINDER.md` §3's, which that note flags as its one genuinely uncertain
input; the *rationality* half of Proposition S (the descent needs $s\in\mathbf Q$) is
unconditional.

**Open / not done.**
1. **$\mathbf Q(x)$-linear independence** of $\{1,H,\theta H\}$ on any host. Every margin in
   (K)/(D) presumes it. This is the analogue of CDT's Lemma 12.1.1. **[open]**
2. The **$m=3$ Bost–Charles optimum is over one one-parameter family** ($\varphi_r$). CDT's
   own contour is a $4$-slit-plus-lune refinement of the same shape; a comparable refinement
   here would improve $-0.83852$ by an unknown amount, plausibly a few hundredths to a tenth.
   **[open]**
3. **A second fold-regular class on $\Gamma_0(5)+5$** — the finite exact computation of
   §7.2. **[open, the next task]**
4. **The identification of the level-5 period** $0.164306701064\ldots$ remains open
   (`NONCONGRUENCE_SCAN.md` §5 excluded all weight-three newforms of level $M\le120$ with
   $5\mid M$).
5. The (D) architecture at $c\ge2$ classes assumes the extra classes share the *same*
   $(\Gamma,t)$, hence the same $x_2$; on a host where they do not, the common domain
   shrinks to the smaller $|x_2|$ and the gain is lost. **[open]**

**No irrationality claim.** The two positive lines of §5.2 (Beukers' row, Apéry's $\zeta(2)$
row) re-derive theorems that are already proved by the elementary route; they are presented
as a calibration of the CDT machinery, and they inherit gaps 1 and 2 above, so they are not
offered as independent proofs.

---

## 9. Reproduction

```
lattice/cdt_noncongruence/hosts_nc.py        # the rows, lam_1, lam_2, delta, lam_2^norm, score
lattice/cdt_noncongruence/verify_rows.gp     # PARI: a_n in Z, k=2 sharp, roots, xi (n<=160)
lattice/cdt_noncongruence/geom_denom.py      # Lemma G: rescaling == uniform p-adic slope
lattice/cdt_noncongruence/arch_k.py          # architectures (K)/(D): inventories and tau
lattice/cdt_noncongruence/bc_multivalent.py  # Delta(r) for phi_r = x_2 lambda(rz)
lattice/cdt_noncongruence/delta_table.py     # fine Delta(r) grid + convergence check
lattice/cdt_noncongruence/optimum.py         # optimise (m, inventory, r) -> thresholds
lattice/cdt_noncongruence/entries.py         # entry conditions in CDT_FINDER units
lattice/cdt_noncongruence/table_nc.py        # margins, deficits, final ranking
lattice/cdt_noncongruence/ingredients.py     # what each extra ingredient is worth
lattice/cdt_noncongruence/index910.py        # index-9/10 groups; the newly opened window
lattice/cdt_noncongruence/20_dims.gp         # PARI: dim M_2 at level 5 and above (Sec.10.3)
lattice/cdt_noncongruence/20_second_class.gp # the level-5 nebentypus row and the s-scan (Sec.10.4)
lattice/cdt_noncongruence/21_higher_levels.gp# M_2(Gamma_0(10/15/20/25)) on the same t-line (Sec.10.5)
lattice/cdt_noncongruence/22_gauge_family.gp # F*P(t)^b, sqrt(P), and the identification of the hits
lattice/cdt_noncongruence/23_costs.py        # the cost of every route to a second class (Sec.10.4-10.6)
```
Logs of every run are beside the scripts (`*.log`). All of §§3–7 is reproduced by
`python3 table_nc.py` after `python3 delta_table.py`.

---

## 10. Is there a second fold-regular class on $(\Gamma_0(5)+5,\ t)$?

*Added after §9. Scripts: `20_dims.gp`, `20_second_class.gp`, `21_higher_levels.gp`,
`22_gauge_family.gp`, `23_costs.py` (logs beside them).*
**The answer is no, for a reason that is structural rather than a search failure**, and
the minimal modification that would supply extra classes moves the host off the row.

### 10.1 What counts as a fold-regular class in architecture (D)

A function $f$ is admissible for architecture (D) iff

1. $f\in\mathbf Q[\![x]\!]$ (rational coefficients — conditionally is allowed),
2. $f$ is holomorphic at $x=0$ and its Taylor coefficients are of CDT type
   $a_n/(n^{e}[1..b_1n]\cdots)$ with $\sigma=\sum_jb_j$ as small as possible,
3. **$f$ is single-valued on $\varphi_r(\mathbf D)=\mathbf C\setminus\{0,x_2\}$, with
   trivial monodromy at $0$ and at $x_1$** — this is *fold-regularity*, and it is what
   the hypothesis buys,
4. the $f$'s are $\mathbf Q(x)$-linearly independent.

Condition 3 does **not** require the $f$'s to solve the same operator $L_1$; it requires
only that their branch locus be $\{x_2,\infty\}$. This is why
`NONCONGRUENCE_SCAN.md` **Lemma 1.1 makes a second weight-two form "free"**: the finite
singular points of the operator attached to $\sqrt{F'}$ are the values $t(P)$ at the
special points of $\Gamma$ and **do not depend on $F'$** — so a second form on the same
$(\Gamma,t)$ carries no geometric cost at all: same $x_2$, same $\Omega$, same $\varphi$,
same $\mathrm{BC}$. The whole cost of a second class is arithmetic, in three places:
$\lambda'$ (Lemma G), $k'$ (the value of $\sigma$), and *the hypothesis*.

**The hypothesis is the binding constraint.** Fold-regularity at $x_1$ is *one linear
condition* on the two-dimensional local solution space at $x_1$; a hypothesised
$\mathbf Q$-linear relation with $\ell$ terms supplies exactly $\ell-1$ such conditions
(CDT's $H_A,H_B,H_C$ are the $\ell-1=2\,{+}\,1$ basis functions of the solution space of
their inhomogeneous ODE with the three free constants $a,b,c$ of
$a+b\zeta(2)+cL(2,\chi_{-3})=0$). A second-order row carries **one** period, so the
relation is $2$-term and the conditional space is $1$-dimensional:
$$\boxed{\ c\ =\ \#\{\text{periods on the host in one relation}\}\ =\ \dim\{\text{usable weight-2 forms with }F(0)=1\}\ }$$
Two rows with two *separate* hypotheses would only prove "not both rational", which is not
the target.

### 10.2 The Galois obstruction: no unconditional class exists

> **Proposition G [proved].** Let $f\in\mathbf Q[\![x]\!]$ be holonomic and let $L$ be its
> minimal operator over $\mathbf Q(x)$. The singular locus of $L$ is stable under
> $\mathrm{Gal}(\overline{\mathbf Q}/\mathbf Q)$. On the level-$5$ host
> $x_1,x_2$ are conjugate over $\mathbf Q(\sqrt5)$, so any $f$ whose operator is singular at
> $x_2$ is singular at $x_1$ as well. Hence **every** admissible $f$ owes its regularity at
> $x_1$ to a codimension-one condition on the solution, i.e. to the hypothesis: there is no
> pure (unconditional) module on this host, of polylogarithmic or of any other type.

A striking instance. Put $P(t)=1-44t-16t^2=\prod_i(1-t/t_i)$, the symbol polynomial. Then
$$\sqrt{P(t)}=\sqrt{1-4(11t+4t^2)}\in\mathbf Z[\![t]\!]$$
(the binomial coefficients $\binom{1/2}{k}4^k=2(-1)^{k-1}\mathrm C_{k-1}$ are integers) —
a **denominator-free** ($\sigma=0$, $\lambda=1$) holonomic function on the host, of order
$2$ with characteristic roots $\tfrac12\lambda_1=44.3607$, $\tfrac12\lambda_2=-0.36068$
(verified, `22_gauge_family.log`). It is exactly the function a pure module would want —
**and it branches at $x_1$ as well as at $x_2$**, because $P$ is irreducible over
$\mathbf Q$. Proposition G in one line.

### 10.3 The dimension count on this host

**[computed]** (`20_dims.log`, PARI `mfdim`).

| space | $\dim M_2$ | $\dim S_2$ | usable ($F(0)=1$, divisor on special points) |
|---|---|---|---|
| $M_2(\Gamma_0(5))$ | $\mathbf1$ | $0$ | $\mathbf1$: $E_{2,5}$ |
| $M_2(\Gamma_0(5),\chi_5)$ (the quadratic, even character) | $2$ | $0$ | see §10.4 |
| $M_2(\Gamma_1(5))=M_2(\Gamma_0(5))\oplus M_2(\Gamma_0(5),\chi_5)$ | $3$ | $0$ | — |
| $M_2(\Gamma_0(10))$, $M_2(\Gamma_0(15))$, $M_2(\Gamma_0(20))$, $M_2(\Gamma_0(25))$ | $3,4,6,5$ | $0,1,1,0$ | see §10.5 |

Two elementary facts fix the shape of the answer immediately.

* **$\Gamma_0(5)+5$ has exactly one cusp**, so every weight-zero modular unit with cuspidal
  divisor has divisor $0$ and is constant: *the family $F\mapsto FU$ of the original brief is
  empty on this host.*
* On $\Gamma_0(5)$ (two cusps) the unit is $u=(\eta_5/\eta_1)^6$, but a weight-two form has
  divisor of degree $\mu/6=1$ and $\operatorname{ord}_\infty(E_{2,5}u^j)=j$, so $j=0$:
  $E_{2,5}$ is the **unique** weight-two form with $\operatorname{ord}_\infty=0$,
  holomorphic, cuspidal divisor. This is $\dim M_2(\Gamma_0(5))=1$ read geometrically.

### 10.4 Route (a): the nebentypus space $M_2(\Gamma_0(5),\chi_5)$

**[computed]** (`20_second_class.log`). The space is $2$-dimensional and Eisenstein; only
one basis direction has $a_0\ne0$, giving
$$E^{\mathbf1,\chi_5}=1-5q+5q^2+10q^3-15q^4-5q^5-10q^6+30q^7+25q^8-\cdots$$
Its row on the $t$-line:
$$\lambda'=4,\qquad\text{minimal recurrence order }\mathbf6,\ \deg 4,$$
$$\text{singular }t\ \in\ \Bigl\{\underbrace{t_1,\ t_2}_{\text{the host's}},\ \underbrace{\tfrac{-6+\sqrt{39}}{12}=0.0204165,\ \tfrac{-6-\sqrt{39}}{12}=-1.0204165}_{\textbf{new}}\Bigr\}.$$
The two new points are a **new** quadratic irrationality, $\mathbf Q(\sqrt{39})$: they are
the images of the zeros of $E^{\mathbf1,\chi_5}$ in $\mathfrak H$, where $\sqrt{F'}$
branches. This is exactly `NONCONGRUENCE_SCAN.md` §1's "non-eta $F$ is strictly worse",
made explicit. **Contour cost:**

| item | value |
|---|---|
| $\lambda'=4$ against $\lambda=2$: the common integral coordinate must be $x=t/4$ | $\log2=0.69315$ |
| the new branch point sits at $|t|=1.0204165<|t_2|=2.7725425$: the common $\Omega$ shrinks | $\log\frac{2.7725425}{1.0204165}=0.99955$ |
| **total archimedean cost** | $\mathbf{1.69270}$ |
| gain from going $c=1\to c=2$ | $+0.58521$ |
| **net** | $\mathbf{-1.10749}$ |

Threshold with this class: $+0.2690$ (against the baseline $-0.83852$); the level-5
deficit worsens from $-0.8349$ to $-1.9424$. The whole one-parameter family
$E^{\mathbf1,\chi_5}+sE^{\chi_5,\mathbf1}$ was scanned over $34$ rationals $s$: **no member
admits any recurrence of order $\le6$ and degree $\le4$** — generically the zeros in
$\mathfrak H$ move and multiply. Route (a) is dead.

### 10.5 Route (b): the same $t$-line at levels $10,15,20,25$

**[computed]** (`21_higher_levels.log`, `22_gauge_family.log`). For every rational basis
element of $M_2(\Gamma_0(N))$, $N=10,15,20,25$, with $a_0\ne0$, normalised to $F'(0)=1$:

| level | $\dim M_2$ | elements with $a_0\ne0$ | order-$2$ with the host's singular set | identification |
|---|---|---|---|---|
| $10$ | $3$ | $3$ | $1$ | **$=E_{2,5}$** |
| $15$ | $4$ | $3$ | $1$ | **$=E_{2,5}$** |
| $20$ | $6$ | $5$ | $1$ | **$=E_{2,5}$** |
| $25$ | $5$ | $2$ | $1$ | **$=E_{2,5}$** |

(verified as an equality of $q$-expansions to $q^{20}$). Every other element either has
non-$2$-power denominators in $[t^n]\sqrt{F'}$ (so it is not a Theorem-R1 row at all) or
admits no recurrence of order $\le6$, degree $\le4$ — again the zeros in $\mathfrak H$.
**Raising the level on the same $t$-line supplies literally nothing new**: the map
$M_2(\Gamma_0(5))\hookrightarrow M_2(\Gamma_0(5m))$ is the only source, and it is the same
one-dimensional space. Route (b) is dead — and it never even reaches the stage where the
extra cusps of $\Gamma_0(5m)$ would have to be paid for.

### 10.6 Route (c): the same operator, gauged — and the $\mathrm{Sym}^w$ tower

**[proved]** Two weight-one forms $g,g'$ on the same $\Gamma$ (up to multiplier) differ by
a weight-zero function, and on a genus-zero $\Gamma$ every weight-zero modular function is a
rational function of $t$. So $g'=g\sqrt{R_0(t)}$, $R_0\in\mathbf Q(t)$; and for $\sqrt{R_0}$
to branch only inside $\{0,t_1,t_2,\infty\}$ one needs $R_0=P(t)^b$. Hence the *complete*
list of same-singular-set deformations is $F'=E_{2,5}\,P(t)^b$, $b\in\mathbf Z_{\ge0}$.

**[computed]** (`22_gauge_family.log`):

| $b$ | order | char. roots | $k'$ | verdict |
|---|---|---|---|---|
| even | $2$ | $\lambda_1,\lambda_2$ | — | $g'=gP^{b/2}$, so $A'=P(4x)^{b/2}A$: **$\mathbf Q(x)$-dependent**, no new function |
| $1$ | $2$, $\deg4$ | $\lambda_1,\lambda_2$ | $\mathbf{10}$ | genuinely new, but $\sigma=10$ instead of $2$ |
| $3$ | $\mathbf4$, $\deg2$ | $\lambda_1,\lambda_1,\lambda_2,\lambda_2$ | $2$ | **doubled roots**: `ROOT_ROWS.md` (H4), the linear form converges polynomially, no fold-regular companion |

The $b=1$ function *is* admissible, and it is priced exactly: with columns
$u_1=u_2=1$ (two layers of rate 1) and $u_3=\dots=u_{10}=3$ (eight further layers carried
by that one function), $m=4$ gives $\tau^\flat=5.375$ and threshold $\mathbf{+3.4802}$ —
worse than the baseline by $4.32$. Route (c) is dead.

**The $\mathrm{Sym}^w$ tower.** Solutions of $\mathrm{Sym}^wL_1$ are $g^{w-j}(g\!\int\!W)^j$;
their companions have $k=w+1$, i.e. $\sigma=w+1\ge3$. Adding the $w=2$ pair
($\sigma=3$) to $\{1,H,\theta H\}$ gives $m=5$, $\tau^\flat=2.56$, threshold
$\mathbf{-0.6237}$ — worse than the baseline by $0.215$; and the $\mathrm{Sym}^2$
conditional function needs its **own** rationality hypothesis (its period is the parent
row's, a different number). Route (d) is dead.

### 10.7 Verdict and the minimal modification of the host

> **No admissible second fold-regular class exists on $(\Gamma_0(5)+5,\ t)$.** The level-5
> row stays at $c=1$, threshold $-0.83852$, deficit $\mathbf{-0.8349}$.

The obstruction is a dimension count, and it is sharp: by §10.1, $c$ equals the number of
usable weight-two forms with $F(0)=1$ on the host, which for a genus-zero group with $s$
cusps is at most $\dim\mathrm{Eis}_2=s-1$; and
$$\Gamma_0(5)+5\ \text{has}\ s=1 .$$
**So the minimal modification is a host with more cusps.** Combining this with the scan's
own classification:

* By `NONCONGRUENCE_SCAN.md` Theorem N2 a three-term (second-order) row forces exactly
  **four special points**, so $s\le4$ and $c\le3$.
* By Theorem N3, for $\lambda_2$ to be an irrational conjugate (the only way
  $|\lambda_2|<1$ with a monic integral symbol) the two finite singular points $t_1,t_2$
  must carry **equal** exponent differences — so the elliptic point, if any, must be the one
  $t$ sends to $\infty$.
* Hence the admissible signatures are exactly the rows of the scan's §3 table with
  $\ge3$ cusps: $(e_2,e_3,\text{cusps})=(1,0,3)$ at index $9$, $(0,1,3)$ at index $10$
  ($c\le2$), and $(0,0,4)$ at index $12$ ($c\le3$).

| target $c$ | needed signature | groups | what is actually there |
|---|---|---|---|
| $2$ | $4$ special points, $3$ cusps | index $9$: $(0;2;3)$, index $10$: $(0;3;3)$ — four non-congruence groups each, plus their congruence analogues | **no row with real $|\lambda_2|<1$ in any box the scan searched** (§5.4), and a non-congruence Hauptmodul pays an extra $\log\mu$ on top of $\log\lambda$ |
| $3$ | $4$ cusps | index $12$: **Beauville's six**, all congruence — Zagier's six sporadic rows | thresholds $-1.7250$; the only one of the six above it is **Apéry's $\zeta(2)$ row** ($+0.406$, already proved). Zagier A and C sit at $-2.000$, still $0.275$ short; Catalan at $-3.386$ |

Two consistency checks that this reading is the right one.

1. **CDT's own host is the $c=3$ case.** $\Gamma_0(6)$ has four cusps, $\dim\mathrm{Eis}_2=3$,
   three periods $1,\zeta(2),L(2,\chi_{-3})$ in one relation, and CDT's $H_A,H_B,H_C$ are
   exactly the $c=3$ conditional basis. Our $c$-ladder therefore reproduces their function
   count from the geometry of the host.
2. **$\Gamma_0(6)+6$, Beukers' host, has $s=2$ cusp orbits** ($\infty\!\sim\!0$,
   $\tfrac12\!\sim\!\tfrac13$ under $W_6$), so $c\le1$ there too — and it does not need more:
   its deficit is already $+0.9777$.

**The cost of moving the level-5 row to a $3$- or $4$-cusp host is that the row does not
exist there.** By Lemma 1.1 the archimedean data $(t_1,t_2)$ is a property of $(\Gamma,t)$,
so changing $\Gamma$ changes $|t_2|$ and hence the score; the scan already searched every
such host exhaustively and found nothing with $|\lambda_2|<1$ beyond Apéry's row, Beukers'
row, $\sqrt T$ and this one. In the honest summary: **the level-5 row is short by $0.835$
nats, and the one lever that could supply $0.585$ of it ($c=2$) is unavailable on its host
for a dimension reason ($s=1$ cusp), while the lever that would supply $0.887$ ($c=3$)
requires four cusps and therefore a different row altogether.**

### 10.8 Ledger for §10

**[proved]** §10.1 (what admissibility means, and $c=\ell-1$); §10.2 (Proposition G, the
Galois obstruction, and the $\sqrt{P}$ instance); the two geometric facts of §10.3;
the classification $F'=E_{2,5}P(t)^b$ of §10.6.

**[computed]** all dimensions (`mfdim`); the level-5 nebentypus row (order $6$, $\lambda'=4$,
new singular points $(-6\pm\sqrt{39})/12$); the $34$-member $s$-scan; the levels
$10,15,20,25$ sweep and the identification of every order-$2$ hit as $E_{2,5}$ itself;
$k'=10$ at $b=1$ and the doubled roots at $b=3$; every threshold in §10.4–10.7.

**[open]** (i) $q$-precision is $80$ and the recurrence search is bounded by order $6$,
degree $4$ — a class whose row satisfies only a larger recurrence would be missed, though by
§10.1 such a class is useless anyway (its conditional space is not one-dimensional).
(ii) The claim "$c$ periods $\Rightarrow$ $c$ fold-regular conditional functions" is CDT's
mechanism transported; it is not proved that the resulting functions are
$\mathbf Q(x)$-linearly independent (the same gap as §8, item 1).
(iii) Whether the level-5 period $0.164306701064\ldots$ lies in $\mathbf Q+\mathbf Q\xi'$
for some second period on a *related* host was not investigated; there is no second period
on this host to test against.
