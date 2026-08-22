# A $p$-adic holonomy census — running CDT's $\zeta_2(5)$ recipe on every near miss

*Claude (Opus 5), 2026-08-22.  Scripts and logs: `lattice/padic_holonomy/`.
Companion to `PADIC_IRRATIONALITY_CENSUS.md` (which scores the same cells by
Calegari's **elementary** criterion $\theta_p>1$), `CDT_FINDER.md` (the same
holonomy bound applied at the archimedean place only), `CDT_UNPACKED.md`,
`ZUDILIN_2ADIC.md`.*

*Sources.  F. Calegari, "Irrationality of certain $p$-adic periods for small $p$",
math/0408214 (`papers/calegari_0408214.txt`); F. Calegari, V. Dimitrov, Y. Tang,
**Arithmetic holonomy bounds and their applications** (ICM survey),
`papers/cdt/cdt/icm/ICM.tex`, §6.2 lines 1063–1201 — this is the source of the
recipe; and **The linear independence of $1,\zeta(2),L(2,\chi_{-3})$**,
arXiv:2408.15403 (`papers/cdt/cdt2/L2chi.tex`) for Theorem 2.5.1/6.0.2/7.0.1,
the rearrangement integral and the lune contours of Appendix A.  For the **crude
(sup-norm) multi-place** holonomicity theorem, with $p$-adic places included:
V. Dimitrov's LNT lecture notes `papers/dimitrov_LNT_zeta2_5.txt`, proved in
V. Dimitrov, arXiv:1912.12545 (the Schinzel–Zassenhaus paper).*

**Nothing here is a claimed theorem.**  §9 is the honest computed/estimated/assumed
ledger; §8 flags the one value that would be new *if* the numerics and the
quoted multi-place form of the bound both hold; **§11 is a full certification of
that one value**, with the residual caveats listed in §11.8.

---

## 0. Verdict

We ran CDT's ICM §6.2 recipe — Buzzard's $p$-adic radius $R_p$ at the finite
place, the Kontsevich–Zagier order-$(2k+1)$ inhomogeneous ODE for the function
inventory, and an optimised archimedean template — on every $p$-adic near miss of
`PADIC_IRRATIONALITY_CENSUS.md`, plus the two known calibrations.

* **The recipe reduces to a single scalar inequality**
  $$\mathrm{cost}(\Omega)=\mathrm{BC}(x\circ\psi)-m\log|\psi'(0)|
  \;<\;d\Bigl(\log R_p-m+\frac1m\Bigr)=\mathrm{budget}(p,k),\qquad m=d+1,\ d=2k+1,$$
  with *all* the arithmetic in the budget and *all* the analysis in the cost.
* **Calibration is exact.**  CDT's two contours for $\zeta_2(5)$ reproduce their
  printed numerators $2.13322$ and $3.92881$ (to $5\times10^{-4}$) and their
  bounds $4.43206$, $4.48866$ (to $5\times10^{-5}$), and their irrationality
  measures $22.0724$, $19.7439$.  In the process we identify which integral each
  of their numbers is: **circle $\to$ Bost–Charles, lune $\to$ rearrangement**.
* **The admissible search space is larger than "pick a domain $\Omega$".**  Since
  $H,H',\dots$ pull back to functions meromorphic on the *whole* $q$-disc, any
  holomorphic $\psi:\mathbb D\to\mathbb D$ with $\psi(0)=0$ works — $\psi$ need not
  be injective.  Optimising over the exact space (via outer functions: the only
  degree of freedom is the boundary modulus $u\le0$, and
  $\log|\psi'(0)|=\bar u$) beats CDT's own contours everywhere: for $\zeta_2(5)$
  the margin goes $3.697\ (\Omega_{\rm circ})\to4.123\ (\Omega_{\rm lune})\to
  \mathbf{4.188}$, improving their irrationality measure from $19.744$ to
  $\mathbf{19.44}$.
* **$\boxed{\zeta_5(3)}$ — the primary target — gives a contradiction**, and §11
  hardens every input of it.  With $m=4$, $R_5=5^3$ (measured from the
  coefficients to $n=2000$, §11.1), $\tau(\mathbf b)=\frac{45}{16}$ (proved,
  §11.2) and an explicit template $\psi(z)=z\exp\bigl(\sum_{k\le40}c_kz^k\bigr)$
  lying in $|q|\le0.8363$:
  $$m\;\le\;3.9716148\ \ (\text{rearrangement, }N=2^{24}\text{, a-priori error }6.8\cdot10^{-4})
  \qquad\text{resp.}\qquad 3.9499804\ \ (\text{Bost--Charles})\qquad <\;4=m .$$
  Certified margin $\ge+0.04214$ out of a budget of $3.23494$ ($1.3\%$);
  irrationality measure $\kappa\le687$ (RE) resp. $390$ (BC).  **The claim needs
  the *refined* multi-place bound**: the published *crude* (sup-norm) form of
  Dimitrov, which suffices for $\zeta_2(5)$ ($5.5556<6$, reproduced), gives only
  $m\le\mathbf{5.7143}$ here — §11.8 decomposes the $1.59$ gap.  **§11.9(e) closes
  the citation gap** for the Bost–Charles numerator: the adelic form with the sharp
  weight $1-\frac1m$ is proved there by an $R_p$-Gauss twist of the evaluation
  lattice in CDT's own slopes proof (L2chi §6.3), the sharpness being confirmed by
  the fact that at $m=2$ it *is* the Borel–Dwork criterion $\rho\prod_pR_p>1$.  So
  the operative statement is the Bost–Charles row — bound $\mathbf{3.9499804}$,
  margin $+0.0742553$, $\kappa\le390$ — and every ingredient is published except
  one new elementary lemma (§11.9(e), Lemma B).  Calegari's
  elementary criterion fails here ($\theta_5=0.8918$) and Beukers 2008 reaches only
  $\zeta_5(3)\pm L_5(3,\chi_5)$, so **this would be a new theorem** — subject to
  the caveats of §11.8, above all the fact that the multi-place form of the bound
  is a citation to CDT's companion paper and not verified here.
* **A host correction.**  The $3$-adic value of the Zagier $\mathbf B$ row is
  $\frac12L_3(2,\chi_{-3})$ ($=\frac12\zeta_3(2)$ in Teichmüller normalisation).
  On the $X_0(9)$ host of the census the bound fails (margin $-0.390$), but the
  same value lives on $X_0(3)$ — $\chi_{-3}$ has conductor $3$ — where
  $R_3=3^6$ rather than $3^3$, the budget jumps from $1.258$ to $7.850$, and the
  contradiction is easy (margin $+4.127$).  The value is already known irrational
  (Beukers 2008), but the census's identification of the *host* was the wrong one.
* **Everything else fails, and fails structurally.**  $\zeta_7(3)$ ($-2.05$),
  $\zeta_3(5)$ ($-2.37$), $\zeta_2(7)$ ($-6.46$) miss by margins no contour can
  recover; and for $\zeta_5(5)$ (and every $\zeta_{13}(\cdot)$) the *entry
  condition* $\log|\varphi'(0)|+\log R_p>\tau(\mathbf b)$ fails for **every**
  admissible template, since $|\varphi'(0)|\le1$ by Schwarz and
  $\log R_5=4.8283<\tau(\mathbf b)=4.8611$.  The deficit is $0.033$.
* **Margins table** — see §5; the one-line version is

  | | $\zeta_2(3)$ | $\zeta_3(3)$ | $L_2(2,\chi_{-4})$ | $\zeta_2(5)$ | $L_3(2,\chi_{-3})^{X_0(3)}$ | $\zeta_5(3)$ | $L_3(2,\chi_{-3})^{X_0(9)}$ | $\zeta_7(3)$ | $\zeta_3(5)$ | $\zeta_2(7)$ | $\zeta_5(5)$ |
  |---|---|---|---|---|---|---|---|---|---|---|---|
  | $m$ | 4 | 4 | 3 | 6 | 3 | **4** | 3 | 4 | 6 | 8 | 6 |
  | bound | 2.276 | 2.691 | 2.028 | 4.439 | 1.871 | **3.9716** | 3.343 | 6.822 | 8.063 | 15.595 | — |
  | margin | $+7.394$ | $+3.795$ | $+2.879$ | $+4.188$ | $+4.127$ | $\mathbf{+0.042}$ | $-0.390$ | $-2.048$ | $-2.367$ | $-6.462$ | n/a |

---

## 1. The recipe, unpacked for general $(p,k)$

### 1.1 The curve, the coordinate, and $R_p$

Let $p-1 \mid 24$ so that $X_0(p)$ has genus $0$ ($p=2,3,5,7,13$), and put
$$x \;=\; \Bigl(\frac{\Delta(p\tau)}{\Delta(\tau)}\Bigr)^{1/(p-1)}
\;=\; q\prod_{n\ge1}\Bigl(\frac{1-q^{pn}}{1-q^{n}}\Bigr)^{24/(p-1)}
\;=\;\Bigl(\frac{\eta(p\tau)}{\eta(\tau)}\Bigr)^{24/(p-1)}\in q+q^2\mathbb Z[\![q]\!].$$
This is a hauptmodul: $x=0$ at the cusp $\infty$, $x=\infty$ at the cusp $0$.
The Fricke involution $w_p:\tau\mapsto-1/(p\tau)$ sends
$$x\longmapsto p^{-12/(p-1)}/x ,$$
so it exchanges the two components $|x|_p\le1$ and $|x|_p\ge p^{12/(p-1)}$ of
the ordinary locus, the supersingular locus being the complementary annulus.  By
Buzzard [Thm 5.2] the finite-slope overconvergent form $E^*_{-2k}$ extends across
the *whole* supersingular locus, so
$$\boxed{\;R_p \;=\; p^{12/(p-1)}\;}\qquad
(p=2:\,2^{12},\quad 3:\,3^{6},\quad 5:\,5^{3},\quad 7:\,7^{2}),$$
and $R_\ell=1$ for every $\ell\ne p$.  For the two **level-$p^2$** rows we use
the same construction with
$$x=\Bigl(\frac{\Delta(p^2\tau)}{\Delta(\tau)}\Bigr)^{1/(p^2-1)},\qquad
R_p=p^{24/(p^2-1)}\qquad(p=2:\,2^{8}=X_1(4);\quad p=3:\,3^{3}=X_0(9)).$$

### 1.2 The function inventory, and $\tau(\mathbf b)$

With $E^*_{2k}=E_{2k}(\tau)-p^{2k-1}E_{2k}(p\tau)$ and
$E'_{-2k}=\sum_{n\ge1}\bigl(\sum_{d\mid n,\,p\nmid d}d^{-2k-1}\bigr)q^n$
(the constant term $\zeta_p(2k+1)/2$ omitted), set
$$H \;=\; E^*_{2k}\bigl(E'_{-2k}+\eta\bigr)\;=\;\sum_{n\ge0}(a_n+\eta\,b_n)\,x^n .$$
$H$ is an overconvergent modular *function* precisely when $\eta=\zeta_p(2k+1)/2$.
Assuming for contradiction that $\zeta_p(2k+1)\in\mathbb Q$ gives $H\in\mathbb Q[\![x]\!]$
of denominator type $[1,\dots,n]^{2k+1}$ (§2 certifies this and its sharpness).
By Kontsevich–Zagier, $H$ satisfies an **inhomogeneous** linear ODE over $\mathbb Q(x)$
of exact minimal order $2k+1$, so
$$\{\,1,\;H,\;H',\;\dots,\;H^{(2k)}\,\}\qquad\text{are }m=2k+2\text{ functions,}$$
$\mathbb Q(x)$-linearly independent (§2 verifies the rank), all of denominator type
$[1,\dots,n]^{2k+1}$ except the first.  So $r=1$, $\mathbf b=(0,d,\dots,d)$ with
$d=2k+1$, $u_1=1$, and
$$\tau(\mathbf b)=\sigma_m-\frac{u_1^2b_1}{m^2}=d\Bigl(1-\frac1{m^2}\Bigr),
\qquad m=d+1 .$$
For the **character rows** — $\chi_{-4}$ on $X_1(4)$ at $p=2$ (Calegari §4), and
$\chi_{-3}$ at $p=3$, which can be run either on $X_0(9)$ (the level-$p^2$
continuation, $R_3=3^3$) or on $X_0(3)$ (possible because $\mathrm{cond}(\chi_{-3})=3$,
$R_3=3^6$; see §5.4) — the weight-one Eisenstein pair $F_1,F'_{-1}$ replaces
$E^*_{2k},E'_{-2k}$, the denominator type is $[1,\dots,n]^2$ and the minimal order
is $2$, so $d=2$, $m=3$.

### 1.3 The bound, and the single inequality it reduces to

CDT's multi-place holonomy bound (ICM §6.2, eq. (6.2)), with the $p$-adic
template $\varphi_p(z)=R_p^{-1}z$ and an archimedean template
$\varphi=\varphi_\infty=x\circ\psi$, reads
$$m \;\le\; \frac{\mathrm{BC}(\varphi)+\sum_v \log R_v}
{\log|\varphi'(0)|+\sum_v\log R_v-\tau(\mathbf b)},\qquad
\mathrm{BC}(\varphi)=\iint_{\mathbb T^2}\log|\varphi(z)-\varphi(w)|\,\mathrm d\mu\,\mathrm d\mu .$$
Writing $L:=\log R_p$ (the only nontrivial term in $\sum_v\log R_v$) and
$\rho:=|\varphi'(0)|$, the contradiction "$\text{bound}<m$" is **exactly**
$$\boxed{\;\mathrm{cost}\;:=\;\mathrm{BC}(\varphi)-m\log\rho\;<\;d\Bigl(L-m+\frac1m\Bigr)=:\;\mathrm{budget}\;}$$
— a clean separation: *budget* is pure arithmetic ($p$ and $k$ only), *cost* is
pure complex analysis (the contour only).  Equivalently
$$\text{bound}=\frac{\mathrm{BC}+L}{\log\rho+L-\tau(\mathbf b)},\qquad
\text{margin}=m\bigl(\log\rho+L-\tau(\mathbf b)\bigr)-(\mathrm{BC}+L)=\mathrm{budget}-\mathrm{cost} .$$
Note $\text{margin}>0\iff\text{bound}<m$, but *minimising the bound* and
*maximising the margin* are different optimisations; the margin is the one that
controls the irrationality measure (§7), so it is the one we optimise, and we
report both.

**Entry condition.**  $\log\rho+L>\tau(\mathbf b)$.  Since $\varphi=x\circ\psi$ with
$\psi:\mathbb D\to\mathbb D$ and $x=q+O(q^2)$, Schwarz gives $\rho=|\psi'(0)|\le1$,
so a *necessary* condition on $(p,k)$ alone is
$$L\;>\;\tau(\mathbf b)=d\Bigl(1-\frac1{m^2}\Bigr).$$

### 1.4 The two numerators

CDT prove the bound both with the Bost–Charles double integral (Thm 7.0.1) and,
by a fully independent high-dimensional argument, with the **rearrangement
integral** (Thm "elementary form"),
$$\mathrm{RE}(\varphi)=\int_0^1 2t\,\bigl(\log|\varphi(e^{2\pi it})|\bigr)^{*}\mathrm dt
=\iint_{[0,1]^2}\max\bigl(g(s),g(t)\bigr)\,\mathrm ds\,\mathrm dt ,\qquad
g(t)=\log|\varphi(e^{2\pi it})| ,$$
$(\cdot)^*$ the **increasing** rearrangement.  Nazarov's argument
(L2chi §rmk\_Nazarov) gives $\mathrm{BC}\le\mathrm{RE}$ always.  We therefore
compute **both**: RE is a one-dimensional integral, numerically far more stable,
and gives the conservative verdict; BC is the sharper one CDT actually quote.

### 1.5 What an admissible archimedean template is

$H(x(q))$, $H'(x(q))$, … are **meromorphic on the whole disc $|q|<1$**: $E^*_{2k}$
is holomorphic on $\mathfrak H$, the $q$-series of $E'_{-2k}$ has bounded
coefficients hence converges on $|q|<1$, and $\mathrm d/\mathrm dx=(x'(q))^{-1}\mathrm d/\mathrm dq$
only introduces poles at the finitely many $\Gamma_0(p)$-orbits of critical points
of $x$.  Hence

> $\varphi=x\circ\psi$ is admissible for **every** holomorphic
> $\psi:\mathbb D\to\mathbb D(0,1)$ with $\psi(0)=0$.

$\psi$ need **not** be injective and $\Omega=\psi(\mathbb D)$ need not be a Jordan
domain.  This is a genuine enlargement of the search space over "choose a
simply-connected $\Omega\subset\mathbb D$ and take its Riemann map", and it is what
we exploit in §4.

---

## 2. The arithmetic inputs, certified from the coefficients

All of §2 is **exact rational arithmetic** in PARI/GP
(`lattice/padic_holonomy/cdt_ab.gp`, `cdt_task1.gp`, `cdt_task23.gp`,
`cdt_task4.gp`, `cdt_task5.gp`; exact $a_n,b_n$ for $n\le60$ in
`data_X0_p_k*.txt`).  Nothing below is assumed from the papers.

### 2.1 Normalisation (a trap)

We use $E^*_{2k}=E_{2k}(\tau)-p^{2k-1}E_{2k}(p\tau)$ **undivided**, i.e. with
constant term $1-p^{2k-1}$, and $\eta=-\lim_n a_n/b_n$.  Rescaling $E^*_{2k}$ to
be $1+O(q)$ (Calegari's normalisation, harmless at $(p,k)=(2,1)$ where the factor
is $-1$) **destroys $b_n\in\mathbb Z$** in every other case: the denominators of
$b_n$ become $7,31,2,13,4,31,6$ for $(2,2),(2,3),(3,1),(3,2),(5,1),(5,2),(7,1)$
respectively.  Calegari's §3 displayed formula $H=\sum(a_n-b_n\eta)$ and his §3
table disagree by a sign (his §4 table is consistent); $\eta=-\lim a_n/b_n$ is
the convergent choice, confirmed independently against Kubota–Leopoldt below.

**Calibration.**  With the extra global $-1$, $(p,k)=(2,1)$ reproduces Calegari's
printed $b_n:1,24,-552,19392,-810024,37210944,-1815620160$ and
$a_n:0,1,1,-\frac{8072}{27},\frac{160841}{9},-\frac{1088512616}{1125},
\frac{175310024408}{3375}$ (his printed $a_6$ is a duplication typo for $a_5$).
The $X_1(4)$ run reproduces his $a_n$ verbatim, his $b_n$ up to sign, and his
$2(-a_6/b_6)=783269/13060350$ exactly.

### 2.2 $\tau(H)=2k+1$ is exact and sharp ($n\le60$, exact arithmetic)

| case | $b_n\in\mathbb Z$ | $[1..n]^{2k+1}a_n\in\mathbb Z$ | least $n$ with $[1..n]^{2k}a_n\notin\mathbb Z$ | offending $\ell$ |
|---|---|---|---|---|
| $X_0(2)$, $k=1,2,3$ | yes | yes | $3,3,3$ | $3$ |
| $X_0(3)$, $k=1,2$ | yes | yes | $4,4$ | $2$ |
| $X_0(5)$, $k=1,2$ | yes | yes | $3,3$ | $3$ |
| $X_0(7)$, $k=1$ | yes | yes | $4$ | $2$ |
| $X_1(4)$, $\chi_{-4}$ | yes | yes ($d=2$) | $3$ | $3$ |
| $X_0(9)$, $\chi_{-3}$ | yes | yes ($d=2$) | $2$ | $2$ |
| $X_0(3)$, $\chi_{-3}$ | yes | yes ($d=2$) | $2$ | $2$ |

So $\tau(\mathbf b)$ is *measured*, not assumed.  **Extra structural fact** (exact,
apparently stronger than anything stated in either paper): the primes that ever
obstruct the smaller exponent are all $\ell\le n$ **except $\ell=p$**, so
$$\bigl(\text{prime-to-}p\ \mathrm{lcm}\{1,\dots,n\}\bigr)^{2k+1}a_n\in\mathbb Z .$$
(Immediate from $E'_{-2k}$ summing only over $p\nmid d$, but not stated in either
paper.)  It gives **no asymptotic gain**, however: deleting the $p$-part removes
only $\lfloor\log_pn\rfloor\log p=O(\log n)$ from
$\log\mathrm{lcm}\{1,\dots,n\}=n+o(n)$, so the rate is still $2k+1$ nats per
index and $\tau(\mathbf b)$ is unchanged.  It is recorded because it is the exact
statement, and because it closes off one obvious hope for the near misses.

### 2.3 The identification of $\eta$

Exact rational arithmetic against (generalised) Bernoulli numbers: taking
$2K\equiv-2k \bmod (p-1)p^N$, the valuation $v_p\bigl(2\eta-(-B_{2K}/2K)\bigr)$
grows linearly in $N$ in every case (the opposite sign stays bounded).  Hence
$$2\eta=\zeta_p(2k+1)\ \text{(all eight $X_0(p)$ cases)},\qquad
2\eta=L_p(2,\chi)\ \text{(character cases)} .$$
For $X_1(4)$, $\chi=\chi_{-4}$: this is Calegari's $L_2(2,\chi_{-4})$, the
$2$-adic avatar of Catalan's constant $G=L(2,\chi_{-4})$, interpolating
$L(-2k,\chi_{-4})=E_{2k}/2$ (Euler numbers).  For $X_0(9)$, $\chi=\chi_{-3}$:
$$\boxed{\ \text{the $X_0(9)$ (Zagier }\mathbf B\text{) 3-adic Apéry limit is }
\tfrac12 L_3(2,\chi_{-3})\ }$$
— the Kubota–Leopoldt $3$-adic $L$-function of the odd quadratic character of
conductor $3$ at $s=2$, i.e. the $3$-adic avatar of $L(2,\chi_{-3})$, interpolating
$L(-2k,\chi_{-3})$.  In the Teichmüller normalisation used in
`PADIC_IRRATIONALITY_CENSUS.md` §3 this is the number written $\zeta_3(2)$
(at $p=3$ the Teichmüller character *is* $\omega=\chi_{-3}$).  A level-$3$ control
run ($X_0(3)$ hauptmodul, same $F_1,F'_{-1}$) converges to the **same** value.
No $3$-stabilisation is applied and none is needed: $\chi_{-3}(d)=0$ whenever
$3\mid d$, so the form is automatically $p$-stabilised, exactly as for $\chi_{-4}$
at $p=2$.  (Using $G_1-G_1(3\tau)$ instead would select the critical,
infinite-slope stabilisation — not Calegari's construction.)

### 2.4 $R_p$ measured from the coefficients

Exact $v_p(a_n+\eta b_n)$ with $\eta=-a_{200}/b_{200}$; the table is the measured
slope $s_n=v_p(a_n+\eta b_n)/n$ against the Buzzard prediction $\log_p R_p$:

| case | target | $n=20$ | $40$ | $60$ | $80$ | $100$ | $120$ |
|---|---|---|---|---|---|---|---|
| $X_0(2)$ $k=1$ | 12 | 11.100 | 11.475 | 11.600 | 11.700 | 11.790 | 11.775 |
| $X_0(2)$ $k=2$ | 12 | 10.850 | 11.300 | 11.417 | 11.588 | 11.710 | 11.667 |
| $X_0(2)$ $k=3$ | 12 | 10.250 | 10.950 | 11.250 | 11.388 | 11.540 | 11.567 |
| $X_0(3)$ $k=1$ | 6 | 5.550 | 5.775 | 5.800 | 5.813 | 5.920 | 5.917 |
| $X_0(3)$ $k=2$ | 6 | 5.400 | 5.575 | 5.717 | 5.763 | 5.800 | 5.817 |
| $X_0(5)$ $k=1$ | 3 | 2.650 | 2.775 | 2.867 | 2.900 | 2.900 | 2.933 |
| $X_0(5)$ $k=2$ | 3 | 2.550 | 2.725 | 2.817 | 2.863 | 2.860 | 2.900 |
| $X_0(7)$ $k=1$ | 2 | 1.800 | 1.850 | 1.883 | 1.913 | 1.950 | 1.967 |
| $X_1(4)$ | 8 | 7.400 | 7.650 | 7.733 | 7.800 | 7.860 | 7.850 |
| $X_0(9)$, $\chi_{-3}$ | 3 | 2.800 | 2.900 | 2.883 | 2.938 | 2.960 | 2.942 |
| $X_0(3)$, $\chi_{-3}$ | 6 | 5.650 | 5.825 | 5.850 | 5.863 | 5.930 | 5.925 |

The *deficit* $\rho n-v_p(a_n+\eta b_n)$ stays $O(\log n)$ throughout (e.g. for
$X_0(2)$, $k=1$: $18,21,24,24,21,27$), which is exactly the shape
$v_p\ge-\frac{12}{p-1}n+O(\log n)$ that pins $R_p=p^{12/(p-1)}$ (resp. $2^8$,
$3^3$).  The bound we use is the *lower* bound from Buzzard, which is all the
holonomy bound needs; the measurement confirms there is no slack to harvest.

### 2.5 Archimedean radius and the branch points (a correction to the folklore)

Measured by $(|b_{2n}|/|b_n|)^{1/n}$ at $N=400$ (immune to phase oscillation),
with the singularity located **exactly** by factoring the relevant modular
polynomial:

| case | $1/R_{\mathrm{arch}}$ | $n=80$ | $160$ | $200$ | nearest singularity (exact) |
|---|---|---|---|---|---|
| $X_0(2)$ | $2^6=64$ | 63.174 | 63.586 | 63.668 | order-2 elliptic point, $x=-1/64$ |
| $X_0(3)$ | $3^3=27$ | 26.614 | 26.806 | 26.845 | **order-3** elliptic point, $x=-1/27$ |
| $X_0(5)$ | $5^{3/2}=11.1803$ | 11.061 | 11.036 | 11.220 | conjugate pair, roots of $125x^2+22x+1$, $\lvert x\rvert =5^{-3/2}$ |
| $X_0(7)$ | $7$ | 6.800 | 7.007 | 6.883 | conjugate pair, roots of $49x^2+13x+1$, $\lvert x\rvert =1/7$ |
| $X_1(4)$ | $2^4=16$ | 15.812 | 15.905 | 15.924 | **cusp** $1/2$, $z=-1/16$ |
| $X_0(9)$, $\chi_{-3}$ | $3^{3/2}=5.19615$ | 5.151 | 5.174 | 5.178 | **cusps** $1/3,2/3$: $u=3^{-3/2}e^{\pm5\pi i/6}$ |
| $X_0(3)$, $\chi_{-3}$ | $3^3=27$ | 26.697 | 26.847 | 26.877 | elliptic point, $u=-1/27$ |

Two corrections worth recording.  (i) $\Gamma_0(3)$ has **no** order-2 elliptic
point ($\nu_2=0$); its branch point sits over $j=0$ and has order 3 — Calegari's
"$1/2+\sqrt{-3}/6$" is the order-3 point.  (ii) $X_0(9)$ has **no** elliptic
points at all; its branch points are the two *cusps* $1/3,2/3$ — the exact
analogue of Calegari's $z=-2^{-4}$ at the cusp $1/2$ of $X_1(4)$.  Their argument
$\pm5\pi/6$ explains the period-$12$ oscillation of $|b_n/b_{n-1}|$ there.
For $X_0(7)$ the *other* factor $2401x^2+245x+1$ has a real root at
$|x|\approx0.00426$, much closer to $0$, but those points have $e=3$ and are
therefore **not** branch points of $x$.

*These archimedean radii do not enter the holonomy bound at all* (only meromorphy
of $\varphi^*f_i$ on $\mathbb D$ does — §1.5).  They are what governs Calegari's
elementary criterion, and they are recorded here because they are the exact
quantity by which the two methods differ.

### 2.6 $\mathbb Q(x)$-independence (rank mod $2^{61}-1$)

With $\eta$ specialised to a rational (5/7; $\eta=0$ gives identical tables), the
matrix of $\{x^jf_i\}_{i\le m,\,j\le D}$ has **full rank $m(D+1)$** in every case
for every $D$ scanned, and the *extended* inventory
$\{1,H,\dots,H^{(2k+1)}\}$ first drops rank by exactly $1$ at a definite $D$:

| case | $m$ | base inventory | first deficiency of the extended inventory |
|---|---|---|---|
| $X_0(2)$ $k=1$ | 4 | full, $D\le6$ | $D=5$, rank $29/30$ |
| $X_0(2)$ $k=2$ | 6 | full, $D\le12$ | $D=11$, rank $83/84$ |
| $X_0(2)$ $k=3$ | 8 | full, $D\le21$ | $D=20$, rank $188/189$ |
| $X_0(3)$ $k=1$ | 4 | full, $D\le6$ | $D=5$, rank $29/30$ |
| $X_0(3)$ $k=2$ | 6 | full, $D\le15$ | $D=14$, rank $104/105$ |
| $X_0(5)$ $k=1$ | 4 | full, $D\le9$ | $D=8$, rank $44/45$ |
| $X_0(5)$ $k=2$ | 6 | full, $D\le19$ | $D=18$, rank $132/133$ |
| $X_0(7)$ $k=1$ | 4 | full, $D\le9$ | $D=8$, rank $44/45$ |
| $X_1(4)$, $\chi_{-4}$ | 3 | full, $D\le4$ | $D=3$, rank $15/16$ |
| $X_0(9)$, $\chi_{-3}$ | 3 | full, $D\le6$ | $D=5$, rank $23/24$ |
| $X_0(3)$, $\chi_{-3}$ | 3 | full, $D\le4$ | $D=3$, rank $15/16$ |

The drop is exactly $1$ at the first deficient $D$ and exactly $2$ at $D+1$: a
single minimal relation, i.e. an inhomogeneous ODE of order exactly $2k+1$
(resp. $2$).  This is independent confirmation of the Kontsevich–Zagier
minimal-order statement, and it is the hypothesis of the holonomy bound that
would otherwise have to be taken on faith.

---

## 3. Calibration: CDT's own $\zeta_2(5)$ contours reproduced

`lattice/padic_holonomy/calib.py`, `calib.log`.  $p=2$, $k=2$: $m=6$, $d=5$,
$L=12\log2=8.3177662$, $\tau(\mathbf b)=\frac{35}{36}\cdot5=\frac{175}{36}=4.8611111$,
**budget** $=5(L-6+\tfrac16)=12.4221645$.

| contour | $\lvert \varphi'(0)\rvert $ | our BC | our RE | CDT's printed numerator | our bound | CDT's bound |
|---|---|---|---|---|---|---|
| $\Omega_{\mathrm{circ}}=\{\,\lvert z+\frac25\rvert \le\frac35\,\}$, $\psi(z)=\frac{z}{2z+3}$ | $1/3$ | **2.133686** | 2.205440 | $2.13322$ | $4.432257$ | $4.43206$ |
| $\Omega_{\mathrm{lune}}=\frac23\cdot\bigl(-h(-z,\frac52)\bigr)$ (L2chi (A.1.1)) | $14/29$ | 3.920949 | **3.929038** | $3.92881$ | $4.488612$ | $4.48866$ |

Both of CDT's printed numerators are reproduced, and the identification of *which*
integral each is falls out unambiguously: **their circle number is the
Bost–Charles double integral** (agreement $4.7\times10^{-4}$; our RE for the same
contour is $2.2054$, nowhere near) and **their lune number is the rearrangement
integral** (agreement $2.3\times10^{-4}$; our BC is $3.9209$).  Reproduced
bounds agree to $5\times10^{-5}$, i.e. to CDT's own printed precision.

Convergence: the lune contour lies in $|q|\le\frac23$ and both integrals are
stable to $10^{-6}$ from $N=4000$ quadrature nodes upward.  The circle contour is
internally tangent to $|q|=1$ at $q=-1$; there $x\to0$, but the contour maps
under $\gamma=\binom{1\;-1}{2\;-1}\in\Gamma_0(2)$ to a curve asymptotic to the
horizontal line $\operatorname{Im}\tau=\pi/6$, along which $\log|x|$ oscillates
**periodically and without limit** in a band $[-4.23,-2.42]$.  Equidistribution
makes the quadrature converge anyway, but only to $\sim10^{-4}$; this is the
source of the residual disagreement with CDT and is why every contour we
*optimise* below is kept strictly inside $|q|\le1-\delta$.

Margins and irrationality measures (§7): $\Omega_{\rm circ}$ gives margin
$+3.6968$ and $\kappa=22.071$ (CDT: $22.0724$); $\Omega_{\rm lune}$ gives margin
$+4.1237$ (RE) and $\kappa=19.746$ (CDT: $19.7439$).  Note the two orderings
disagree: the lune has the *larger* bound ($4.4886>4.4323$) but the *larger*
margin, and it is the margin that improves $\kappa$ — exactly as CDT report.

---

## 4. Optimising the archimedean template

### 4.1 The exact search space

By §1.5 the admissible $\varphi$ are exactly $x\circ\psi$ for holomorphic
$\psi:\mathbb D\to\mathbb D$, $\psi(0)=0$.  Write $\psi(z)=z\,G(z)$ with
$G:\mathbb D\to\overline{\mathbb D}$.  For a **prescribed boundary modulus**
$|G|=e^{u}$ on $|z|=1$ ($u\le0$), Jensen gives
$\log|\psi'(0)|=\log|G(0)|\le\frac1{2\pi}\int u$, with equality iff $G$ is the
**outer** function of $u$; and then the boundary contour is
$$\psi(e^{2\pi it})=\exp\bigl(2\pi it+u(t)+i\,\widetilde u(t)\bigr),\qquad
\widetilde u=\text{conjugate function of }u .$$
So *the whole optimisation is over the single real function $u\le0$*, with
$\log|\varphi'(0)|=\bar u$ exactly.  Since $|x(\bar q)|=|x(q)|$, the optimum can
be taken even in $t$, and we parametrise by a cosine polynomial.  The final
templates are therefore fully explicit:
$$\boxed{\ \psi(z)=z\exp\Bigl(c_0+\sum_{k=1}^{K}c_kz^k\Bigr),\quad c_k\in\mathbb R,\qquad
\log|\varphi'(0)|=c_0,\qquad
\log|\varphi(e^{2\pi it})|\!=\!\log\bigl|x(\psi(e^{2\pi it}))\bigr| ,\ }$$
admissible as soon as $\max_t\bigl(c_0+\sum_kc_k\cos2\pi kt\bigr)\le0$ (checked
on a $2\times10^5$-point grid; $c_0+\sum_k|c_k|\le0$ is the cheap sufficient
test).  This family strictly contains discs, off-centre discs, tangent discs,
CDT's lunes and gobbles, and every slit domain.

### 4.2 Evaluating $x$ anywhere in the $q$-disc

The naive $\eta$-product is useless near $|q|=1$.  We evaluate
$$\log|x_N(\tau)|=\frac1e\bigl(\log|\Delta(N\tau)|-\log|\Delta(\tau)|\bigr),
\qquad e=24/s,$$
with $\log|\Delta|$ computed by $\mathrm{SL}_2(\mathbb Z)$-reduction
($\log|\Delta(\tau)|=\log|\Delta(\gamma\tau)|-12\log|c\tau+d|$).  This is
branch-free, has no convergence problem anywhere in $\mathfrak H$, and agrees
with the direct product to $3\times10^{-15}$ wherever the latter converges
(`haupt.py::log_absDelta`, `logabs_x`).  Complex values of $x$ (needed only for
BC) come from a $\Gamma_0(N)$-reduction implemented by Gauss/Lagrange reduction of
the lattice $\mathbb Z\!\cdot\!N\tau+\mathbb Z$ followed by a small enumeration to
enforce $\gcd(c,d)=1$ (`haupt.py::reduce_g0N`).  Unit tests: $x(\frac{1+i}2)=-2^{-6}$,
$x_{X_0(3)}=-3^{-3}$ at the order-3 point, $|x_{X_0(5)}|=5^{-3/2}$ at
$\frac{2+i}5$, $z_{X_1(4)}(\frac12)=-2^{-4}$, $|u_{X_0(9)}(\frac13)|=3^{-3/2}$,
$\Gamma_0(2)$-invariance to $10^{-16}$, and the exact identity
$\int_{|q|=r}\log|x|\,\mathrm d\mu=\log r$.

### 4.3 Where the contour may and may not touch $|q|=1$

$|x|\to\infty$ at every $\Gamma_0(p)$-image of the cusp $0$ — the points
$q=e^{2\pi ia/c}$ with $p\nmid c$, which are **dense** in $|q|=1$.  A contour may
therefore meet $|q|=1$ only in a measure-zero set, and only at images of the
cusp $\infty$ (where $x=0$), i.e. $q=\zeta_p^{\,j}$ at depth $c=p$.

*A trap worth recording* (§11.8(b)): a disc internally **tangent at an image of
the cusp $0$** has *bounded a.e. boundary values*, because the tangential
approach lands on a horocycle at finite height — for $\Omega=D(0.17,0.83)$ at
$p=5$ the sampled boundary maximum of $\log|x|$ is $-0.027$.  It is nevertheless
inadmissible: $x\circ\psi$ is **unbounded on $\mathbb D$** (the point $q=0.9$ is
interior with $\log|x|=70.1$; a grid over $\Omega$ reaches $6271$), i.e. it has
singular-inner-function behaviour and is not holomorphic on $\overline{\mathbb D}$.
Any numerical search must carry a maximum-principle guard.  For a disc
internally tangent at such a point with radius $b$, the contour maps under the
$\Gamma_0(p)$ element killing that cusp to a curve asymptotic to the horizontal
line
$$\operatorname{Im}\tau=\frac{\pi(1-b)}{p^{2}b},$$
so the *cost of tangency grows like $p^{2}$*: at $p=2$, $b=\frac35$ this is
$\pi/6$ and $|x|$ stays in $[0.0145,0.089]$ (CDT's circle); at $p=5$, $b=\frac35$
it is $0.0838$, where $|x|$ on that line reaches $\sim4\times10^{4}$.  This is why
tangency is decisive for CDT at $p=2$ and **useless for every $p\ge3$** — a fact
the numerics confirms independently (§5, the `T*` rows are never best).
For $X_1(4)$ and $X_0(9)$ the *shallow* cusps $\frac12$ resp. $\frac13,\frac23$
(depth $c=2$ resp. $3$) also have bounded $x$ and are admissible tangency points.

### 4.4 The two numerators, numerically

$\mathrm{RE}$ is a $1$-dimensional integral: sample $g$ at $N$ equispaced $t$,
sort, and take $\sum_i g_{(i)}(2i-1)/N^2$.  $\mathrm{BC}$ is the $N^2$ sum with
the diagonal handled by $\log|\varphi(z)-\varphi(w)|=\log|z-w|+\log\bigl|\frac{\varphi(z)-\varphi(w)}{z-w}\bigr|$
and $\iint\log|z-w|=0$ exactly.  All optimisations were run against $\mathrm{RE}$
(the conservative, rigorous-side numerator) at $N=1024$ and every reported number
was recomputed at $N=4096,\dots,32768$; only $N$-converged values are quoted.

---

## 5. The census

All contours below are the **free-template** optima of §4.1 (scripts `freeopt.py`,
`certify.py`, `assemble.py`; results in `FINAL.json`, logs `free_*.log`).  Every
number is $N$-converged: RE and BC agree to $\le10^{-7}$ between $N=4096$ and
$N=32768$ quadrature nodes, and every contour stays in $|q|\le0.95$.

| target | host | $m$ | $R_p$ | $\log R_p$ | $\tau(\mathbf b)$ | budget | best $\log\lvert \varphi'(0)\rvert $ | RE | BC | bound (RE) | bound (BC) | margin (RE) | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| $\zeta_2(3)$ | $X_0(2)$, $k{=}1$ | 4 | $2^{12}$ | 8.31777 | 2.81250 | 13.70330 | -1.215476 | 1.447778 | 1.443172 | **2.276462** | 2.275389 | +7.39361 | **contradiction** |
| $\zeta_3(3)$ | $X_0(3)$, $k{=}1$ | 4 | $3^{6}$ | 6.59167 | 2.81250 | 8.52502 | -0.880279 | 1.208721 | 1.199899 | **2.690817** | 2.687774 | +3.79518 | **contradiction** |
| $L_2(2,\chi_{-4})$ | $X_1(4)$ | 3 | $2^{8}$ | 5.54518 | 1.77778 | 5.75702 | -0.804653 | 0.463751 | 0.386521 | **2.028161** | 2.002094 | +2.87931 | **contradiction** |
| $\zeta_2(5)$ | $X_0(2)$, $k{=}2$ | 6 | $2^{12}$ | 8.31777 | 4.86111 | 12.42216 | -0.774053 | 3.590146 | 3.580998 | **4.438942** | 4.435532 | +4.18770 | **contradiction** |
| $L_3(2,\chi_{-3})$ | $X_0(3)$, $\chi_{-3}$ | 3 | $3^{6}$ | 6.59167 | 1.77778 | 7.85001 | -1.157963 | 0.249344 | 0.234164 | **1.871210** | 1.867058 | +4.12678 | **contradiction** |
| $\zeta_5(3)$ | $X_0(5)$, $k{=}1$ | 4 | $5^{3}$ | 4.82831 | 2.81250 | 3.23494 | -0.531289 | 1.066965 | 1.035529 | **3.971615** | 3.949980 | +0.04214 | **contradiction** (certified, §11) |
| $L_3(2,\chi_{-3})$ | $X_0(9)$, $\chi_{-3}$ | 3 | $3^{3}$ | 3.29584 | 1.77778 | 1.25834 | -0.380379 | 0.507585 | 0.320415 | **3.343139** | 3.178620 | -0.39038 | no contradiction |
| $\zeta_7(3)$ | $X_0(7)$, $k{=}1$ | 4 | $7^{2}$ | 3.89182 | 2.81250 | 0.42546 | -0.353806 | 1.057868 | 1.000917 | **6.822317** | 6.743819 | -2.04763 | no contradiction |
| $\zeta_3(5)$ | $X_0(3)$, $k{=}2$ | 6 | $3^{6}$ | 6.59167 | 4.86111 | 3.79170 | -0.583663 | 2.656254 | 2.642173 | **8.063418** | 8.051140 | -2.36653 | no contradiction |
| $\zeta_2(7)$ | $X_0(2)$, $k{=}3$ | 8 | $2^{12}$ | 8.31777 | 6.89062 | 3.09936 | -0.576370 | 4.949915 | 4.937999 | **15.594892** | 15.580886 | -6.46151 | no contradiction |
| $\zeta_5(5)$ | $X_0(5)$, $k{=}2$ | 6 | $5^{3}$ | 4.82831 | 4.86111 | $-5.02510$ | — | — | — | — | — | — | **method inapplicable** |

Read the verdict column against $m$: a contradiction is exactly `bound < m`.

### 5.1 $\zeta_5(5)$: the method cannot even start

For $(p,k)=(5,2)$ one has $L=3\log5=4.828314$ and
$\tau(\mathbf b)=\frac{35}{36}\cdot5=4.861111>L$.  Since $\varphi=x\circ\psi$ with
$\psi:\mathbb D\to\mathbb D$, $\psi(0)=0$ forces $|\varphi'(0)|\le1$ (Schwarz), the
denominator $\log|\varphi'(0)|+L-\tau(\mathbf b)\le-0.0328$ is **negative for every
admissible template**.  This is a structural obstruction, not a numerical one: no
choice of contour can help, and the deficit $\tau-L=0.0328$ is small enough to be
tantalising.  Closing it would need $R_5>5^3$ (excluded by the measured
coefficient slopes of §2.4), or a better denominator array (excluded: §2.2 shows
$2k+1$ is sharp and the prime-to-$p$ refinement gains only $O(\log n)$), or a
larger function inventory with a smaller $\tau$ — i.e. extra *unconditional*
functions of the same denominator type, in CDT's level-$6$ style.  That is the
one live route, and it is not attempted here (§9 A4).

### 5.2 The primary target $\zeta_5(3)$

$p=5$, $k=1$: $m=4$, $d=3$, $R_5=5^3$, $L=3\log5=4.8283137$,
$\tau(\mathbf b)=3(1-\tfrac1{16})=\tfrac{45}{16}=2.8125$, **budget** $=3(L-4+\tfrac14)=3.2349412$.

**(a) A simple contour that already suffices.**  Take
$\Omega=\tfrac{13}{20}\cdot L(c)^{-}$, i.e. $\psi=\mathrm{SCALE}(0.65)\circ\mathrm{BITE}(0,c)$
with $c=4.149353$ — the disc of radius $0.65$ with a single CDT lune-bite of
parameter $c$ taken out towards $q=+1$ (the cusp $0$, where $|x|$ blows up).
Then $|\varphi'(0)|=0.65\cdot\frac{c^2-1}{c^2+1}=0.5786$ and (`precise5.log`,
stable from $N=8000$ to $N=32000$)
$$\mathrm{RE}=1.0455734,\quad \text{bound}=3.9992789<4,\qquad
\mathrm{BC}=0.9813230,\quad\text{bound}=3.9555336<4 .$$

**(b) The optimised template.**  With $\psi(z)=z\exp(c_0+\sum_{k\le8}c_kz^k)$ and
$$c=(-0.5335130,\,-0.2311512,\,-0.1257210,\,-0.0056622,\,-0.0269403,\,-0.0354941,\,0.0566503,\,0.0199393,\,-0.0138686)$$
(so $|\varphi'(0)|=e^{c_0}=0.586529$, $\max_t\log|\psi|=-0.2629$, i.e. the contour
lies in $|q|\le0.7688$), the values are stable to $10^{-7}$ from $N=4096$ to
$N=65536$:
$$\mathrm{RE}=1.0658811,\ \ \text{bound}=3.9763826,\ \ \text{margin}=+0.0350080;\qquad
\mathrm{BC}=1.0286380,\ \ \text{bound}=3.9512575,\ \ \text{margin}=+0.0722511 .$$
Extending the degree to $40$ and polishing directly in the $c_k$ (`pol_05_K40.json`,
§11.7) gives the best found,
$$\boxed{\ \mathrm{RE}=1.0669650,\quad m\le 3.9716148<4=m,\quad\text{margin}\ge+0.0421385\ }$$
and $\mathrm{BC}=1.0355293$, bound $3.9499804$, margin $+0.0742553$.  Ten
independent optimisations (degrees $8$–$40$, three parametrisations, three
boundary constraints) all land in $\text{margin}\in[+0.0394,+0.0428]$, so the
family has converged.  **§11 certifies this**: the error term is a-priori, the
denominator type is proved rather than measured, and $R_5=5^3$ is read off the
coefficients to $n=2000$.

Irrationality measure (§7): $\kappa\le687$ from RE, $\kappa\le390$ from BC.  These
are large because the margin is small; the qualitative conclusion is the point.

### 5.3 The near misses, and why they fail

The failure is always the **budget**, i.e. arithmetic, never the contour.  Since
$\mathrm{budget}=d\bigl(L-m+\tfrac1m\bigr)$ with $m=d+1$, and $L=\frac{12\log p}{p-1}$
resp. $\frac{24\log p}{p^2-1}$:

| target | $L$ | $m$ | budget | best cost found | deficit |
|---|---|---|---|---|---|
| $\zeta_5(3)$ | $4.8283$ | 4 | $3.2349$ | $3.1941$ | $\mathbf{-0.0408}$ (**succeeds**) |
| $X_0(9)$ model of $L_3(2,\chi_{-3})$ | $3.2958$ | 3 | $1.2583$ | $1.6487$ | $+0.3904$ |
| $\zeta_7(3)$ | $3.8918$ | 4 | $0.4255$ | $2.4731$ | $+2.0476$ |
| $\zeta_3(5)$ | $6.5917$ | 6 | $3.7917$ | $6.1582$ | $+2.3665$ |
| $\zeta_2(7)$ | $8.3178$ | 8 | $3.0994$ | $9.5609$ | $+6.4615$ |
| $\zeta_5(5)$ | $4.8283$ | 6 | $-5.0251$ | — | entry condition fails |

The structural reason is transparent from the boxed inequality: the budget is
$d\,L$ minus $d(m-\frac1m)\approx d\,m$.  Raising $k$ costs $\sim(2k+1)(2k+2)$ while
$L$ is fixed, and raising $p$ shrinks $L$ like $12\log p/p$.  So the *only* window
is small $k$ and small $p$, and the census exhausts it.

The cost side is essentially saturated: over the whole free family the cost is
remarkably insensitive to the host (all optimal contours have
$\max|q|\in[0.62,0.95]$ and $\log|\varphi'(0)|\in[-1.22,-0.35]$), so no amount of
further contour engineering can move a deficit of $+0.39$, let alone $+2.0$.

### 5.4 The $X_0(9)$ vs $X_0(3)$ lesson (a correction to `PADIC_IRRATIONALITY_CENSUS.md` §6.1)

`PADIC_IRRATIONALITY_CENSUS.md` identifies Zagier $\mathbf B$ at $p=3$ with the
$X_0(p^2)$ continuation of Calegari's $X_1(4)$ construction, and reads off
$S_3=\frac{12\log3}{9-1}-2=-0.35208$: a near miss.  That is correct **for that
host**, and the holonomy bound fails there too (margin $-0.390$).  But the value
$L_3(2,\chi_{-3})$ does **not** need level $9$: $\chi_{-3}$ has conductor $3$, so
$F_1F_{-1}$ is already a modular function on $\Gamma_0(3)$ (weight $0$, character
$\chi_{-3}^2=1$), the hauptmodul is $x=(\Delta(3\tau)/\Delta(\tau))^{1/2}$, and

$$R_3 = 3^{6}\quad\text{not}\quad 3^{3}$$

— verified by the measured coefficient slopes (§2.4, the last row: $\to6$, not
$3$), with $d=2$ still sharp and the rank still full.  The budget jumps from
$1.2583$ to $\mathbf{7.8500}$ and the contradiction is easy (margin $+4.127$,
bound $1.871$).  So the level-$9$ framing was simply the wrong host for that
number; the same descent is **not** available at $p=2$, where $\chi_{-4}$ has
conductor $4$ and $X_1(4)$ is already minimal.

---

## 6. The contours, described numerically

Every optimal template is
$\psi(z)=z\exp\bigl(c_0+\sum_{k\ge1}c_kz^k\bigr)$ with real $c_k$, so the region
is described completely by $\rho(\theta)=|\psi(e^{i\theta})|=\exp\bigl(c_0+\sum_kc_k\cos k\theta\bigr)$;
the coefficient vectors are in `lattice/padic_holonomy/FINAL.json`.  Below,
$\theta/2\pi$ runs over $j/12$, and we give $\rho$ and $\log|x(\psi(e^{i\theta}))|$.
$\theta=0$ points at the cusp $0$ (where $x=\infty$); the contour always pulls in
there and bulges out towards the cusp-$\infty$ directions.

| $\theta/2\pi$ | 0 | 1/12 | 2/12 | 3/12 | 4/12 | 5/12 | 6/12 |
|---|---|---|---|---|---|---|---|
| **$\zeta_2(3)$** $\rho$ | 0.1668 | 0.1739 | 0.1995 | 0.2715 | 0.5462 | 0.5144 | 0.4046 |
| $\log\lvert x\rvert $ | 2.70 | 2.70 | 2.71 | 2.69 | $-3.40$ | $-9.84$ | $-10.91$ |
| **$\zeta_2(5)$** $\rho$ | 0.2886 | 0.3053 | 0.3787 | 0.6125 | 0.6241 | 0.6004 | 0.4960 |
| $\log\lvert x\rvert $ | 7.57 | 7.53 | 7.44 | $-2.77$ | $-3.80$ | $-10.33$ | $-14.07$ |
| **$\zeta_2(7)$** $\rho$ | 0.3780 | 0.4081 | 0.6028 | 0.6143 | 0.6602 | 0.6854 | 0.5643 |
| $\log\lvert x\rvert $ | 11.97 | 12.01 | 8.50 | $-4.75$ | $-2.91$ | $-9.20$ | $-17.25$ |
| **$\zeta_3(3)$** $\rho$ | 0.2615 | 0.2767 | 0.3447 | 0.5610 | 0.4598 | 0.5009 | 0.6351 |
| $\log\lvert x\rvert $ | 3.22 | 3.22 | 3.22 | $-2.26$ | $-5.53$ | $-5.56$ | 0.65 |
| **$\zeta_3(5)$** $\rho$ | 0.3841 | 0.4159 | 0.5950 | 0.6998 | 0.5427 | 0.6347 | 0.6873 |
| $\log\lvert x\rvert $ | 7.16 | 7.15 | 3.42 | $-3.92$ | $-7.13$ | $-6.41$ | 2.17 |
| **$\zeta_5(3)$** $\rho$ | 0.4046 | 0.4416 | 0.6082 | 0.6077 | 0.7619 | 0.6443 | 0.6822 |
| $\log\lvert x\rvert $ | 3.90 | 3.89 | $-0.53$ | $-3.15$ | $-1.43$ | $-3.49$ | 0.36 |
| **$\zeta_7(3)$** $\rho$ | 0.5027 | 0.5663 | 0.6804 | 0.8380 | 0.7811 | 0.7696 | 0.7856 |
| $\log\lvert x\rvert $ | 4.31 | 4.26 | $-1.68$ | $-1.93$ | $-1.76$ | $-3.21$ | 1.92 |
| **$L_2(2,\chi_{-4})$** $\rho$ | 0.2601 | 0.3312 | 0.5224 | 0.9365 | 0.5224 | 0.3312 | 0.2601 |
| $\log\lvert x\rvert $ | 1.78 | 1.77 | $-3.75$ | $-2.77$ | $-3.75$ | 1.77 | 1.78 |
| **$L_3(2,\chi_{-3})^{X_0(3)}$** $\rho$ | 0.1848 | 0.2260 | 0.4385 | 0.5979 | 0.4385 | 0.2260 | 0.1848 |
| $\log\lvert x\rvert $ | 1.20 | 1.20 | $-4.77$ | $-0.17$ | $-4.77$ | 1.20 | 1.20 |
| **$L_3(2,\chi_{-3})^{X_0(9)}$** $\rho$ | 0.4302 | 0.6595 | 0.9348 | 0.7918 | 0.9348 | 0.6595 | 0.4302 |
| $\log\lvert x\rvert $ | 1.89 | 0.41 | $-1.65$ | 1.43 | $-1.65$ | 0.41 | 1.89 |

(the profiles are even in $\theta$, so $j/12$ for $j=7,\dots,11$ mirrors
$j=5,\dots,1$).  Three qualitative features are stable across all hosts:

1. **$\rho$ is minimal at $\theta=0$** (the cusp $0$) and the minimum deepens as
   the target gets harder — $0.167$ at $\zeta_2(3)$, $0.405$ at $\zeta_5(3)$,
   $0.503$ at $\zeta_7(3)$.  The pole at $q=1$ is by far the dominant cost; the
   next ones, at $q=-1$ ($c=2$) and $q=e^{\pm2\pi i/3}$ ($c=3$), are weaker by
   $c^{2}$.
2. $\log|x|$ is **flat and positive** on the pulled-in arc near $\theta=0$: the
   optimiser equalises the boundary modulus exactly where the pole is, which is
   the signature of a level-set (topographic) contour, cf. CDT's Figure in
   L2chi §A.
3. **Tangency to $|q|=1$ is never used for $p\ge3$**: every structured run with a
   tangent disc (`res_*_RE.json`, the `T*` rows) is beaten by an interior
   contour, exactly as predicted by the $\pi(1-b)/p^2b$ horocycle-height formula
   of §4.3.  At $p=2$ tangency is competitive but still not optimal: CDT's
   $\Omega_{\rm circ}$ gives margin $+3.697$ where a free interior contour gives
   $+7.394$ for $\zeta_2(3)$ and $+4.188$ for $\zeta_2(5)$.

---

## 7. Irrationality measures

CDT's measure refinement (ICM Prop. `theomainreprise` / eq. `weakermeasure`) is,
in the multi-place form used for $\zeta_2(5)$,
$$m\;\le\;\frac{\mathrm{BC}+L}
{\log\rho+L-\tau(\mathbf b)-L\cdot\bigl(2(1-\gamma)\kappa^{-1}-(1-\gamma)^2\kappa^{-2}\bigr)},
\qquad \gamma=\frac{\dim_{\mathbb Q(x)}\mathcal V(\varphi,\vec b)}{m}=\frac1m$$
(only the constant function $1$ is unconditionally in $\mathbb Q[\![x]\!]$), with
$\rho_p=R_p^{-1}$ so that $\log(1/\rho_p)=L$.  Setting the right-hand side equal
to $m$ and solving the resulting quadratic in $\kappa^{-1}$ gives the closed form
we use throughout:
$$\boxed{\ \kappa\;=\;\frac{1-\gamma}{1-\sqrt{\,1-\dfrac{\text{margin}}{mL}\,}}\ ,
\qquad \gamma=\frac1m ,\quad \text{margin}=m\bigl(\log\rho+L-\tau(\mathbf b)\bigr)-(\mathrm{BC}+L).}$$
**Check against CDT.**  Their circle ($\text{margin}=3.69727$, $m=6$,
$L=12\log2$) gives $\kappa=22.071$ — they print $22.0724$.  Their lune
($\text{margin}=4.12345$) gives $\kappa=19.746$ — they print $19.7439$.  The
formula is therefore the right one, and the residual $10^{-3}$ is their rounding
of the numerator.

$\kappa$ is monotone decreasing in the margin, so the margin — not the bound — is
the right optimisation target, and the two disagree (§3).  For the six cells where
the contradiction holds:

| value | margin (RE) | $\kappa$ (RE) | margin (BC) | $\kappa$ (BC) | previous best |
|---|---|---|---|---|---|
| $\zeta_2(3)$ | $+7.3936$ | $6.4$ | $+7.3982$ | $6.3$ | Calegari's $\theta=1.16188$ gives $\theta/(\theta-1)=7.18$ |
| $\zeta_3(3)$ | $+3.7952$ | $10.0$ | $+3.8040$ | $10.0$ | Calegari's $\theta=1.04699$ gives $22.28$ |
| $L_2(2,\chi_{-4})$ | $+2.8793$ | $7.4$ | $+2.9565$ | $7.2$ | Calegari's $\theta=1.16188$ gives $7.18$ |
| $\zeta_2(5)$ | $+4.1877$ | $19.4$ | $+4.1968$ | $19.4$ | CDT print $\kappa<19.7439<20$ |
| $L_3(2,\chi_{-3})$ on $X_0(3)$ | $+4.1268$ | $6.0$ | $+4.1420$ | $6.0$ | — |
| $\zeta_5(3)$ | $+0.0421$ | $687.2$ | $+0.0743$ | $389.8$ | **no previous value** |

So this route improves CDT's own $\zeta_2(5)$ measure from $19.744$ to $19.39$,
and Calegari's $\zeta_3(3)$ measure from $22.3$ to $10.0$; it is slightly worse
than his for the $2$-adic Catalan constant ($7.2$ vs $7.18$), where his elementary
argument is already very efficient.  For $\zeta_5(3)$ the measure is large
precisely because the margin is $1.3\%$ of the budget.

---

## 8. Novelty — what would be new, and against what it must be checked

**Not an audit.** These are flags for one.  The three references that must be
checked, per the brief, are

* **Calegari 2005** (math/0408214) — proves $\zeta_2(3)$, $\zeta_3(3)$,
  $L_2(2,\chi_{-4})$; explicitly *fails* at $\zeta_2(5)$ ($\theta=0.90816$),
  $\zeta_5(3)$ ($\theta=0.89179$), $\zeta_7(3)$, $\zeta_{13}(3)$.
* **Beukers 2008**, "Irrationality of some $p$-adic $L$-values" — per
  `PADIC_IRRATIONALITY_CENSUS.md` §6 reaches $\zeta_2(3)$, $\zeta_2(2)$,
  $\zeta_3(2)$, $\zeta_3(3)$; at $p=5$ only the *combination*
  $\zeta_5(3)\pm L_5(3,\chi_5)$, **not** $\zeta_5(3)$ alone.
* **Lai–Sprang–Zudilin 2025** — a new proof of $\zeta_2(5)\notin\mathbb Q$ via a
  Calabi–Yau local system (CDT ICM, remark after Thm `old2adic`).  **Whether they
  also cover $\zeta_5(3)$ or $\zeta_p(3)$ for $p=5,7$ has not been checked here.**

Status of each value reached by this recipe:

| value | verdict here | already known? |
|---|---|---|
| $\zeta_2(3)$ | contradiction, large margin | Calegari 2005 Thm 3.3; Beukers 2008 — **not new** |
| $\zeta_3(3)$ | contradiction, large margin | Calegari 2005 Thm 3.4 — **not new** |
| $L_2(2,\chi_{-4})$ ($=\zeta_2(2)$, $2$-adic Catalan) | contradiction, large margin | Calegari 2005 Thm 4.2; Beukers 2008 Cor 7.3 — **not new** |
| $\zeta_2(5)$ | contradiction, margin $+4.19$, improves CDT's $\kappa$ | CDT ICM 2025; Lai–Sprang–Zudilin 2025 — **not new**, this is the calibration |
| $L_3(2,\chi_{-3})$ ($=\zeta_3(2)$) on the **level-3** model | contradiction, margin $+4.13$ | Beukers 2008 Cor 7.3 — **not new** (but the level-9/Zagier-$\mathbf B$ model *fails*, see §5) |
| $\boxed{\zeta_5(3)}$ | **contradiction**, certified margin $\ge+0.0421$ (RE) / $+0.0743$ (BC), §11 | Calegari fails; Beukers reaches only $\zeta_5(3)\pm L_5(3,\chi_5)$; **would be NEW** — flag for audit against LSZ 2025 |
| $\zeta_2(7)$, $\zeta_3(5)$, $\zeta_7(3)$ | **no contradiction** (margins $-6.46$, $-2.35$, $-2.04$) | open, and stays open by this route |
| $\zeta_5(5)$ | **method inapplicable**: the entry condition fails for *every* template | open |

**The one candidate new theorem is $\zeta_5(3)\notin\mathbb Q$**, i.e. the $5$-adic
zeta value at $3$, equivalently $\zeta_5(3)=L_5(3,\mathrm{id})$ in Calegari's
notation.  Its margin is small ($+0.0421$ out of a budget of $3.2349$, i.e.
$1.3\%$), so §9 and especially §11.8 must be read before any claim is made.

---

## 9. Computed / estimated / assumed — the honest ledger

### (E) Exact, certified

* Every $a_n,b_n$ for $n\le60$ in all eleven cases; $b_n\in\mathbb Z$;
  $[1..n]^{2k+1}a_n\in\mathbb Z$ and the sharpness of $2k+1$ (least failing $n$
  and the offending prime); the prime-to-$p$ refinement (§2.2).
* $\tau(\mathbf b)=d(1-1/m^2)$ from that data; the exact rationals
  $\frac{175}{36},\frac{45}{16},\frac{16}{9}$.
* The identification $2\eta=\zeta_p(2k+1)$ resp. $L_p(2,\chi)$, by exact rational
  arithmetic against (generalised) Bernoulli numbers (§2.3).
* $\mathbb Q(x)$-independence of $\{1,H,\dots,H^{(2k)}\}$ (rank mod $2^{61}-1$,
  full for every scanned $D$), and the exact minimal ODE order $2k+1$ (§2.6).
* The exact algebraic location of all archimedean branch points (§2.5).

### (C) Numerically computed, converged, but not interval-certified

* $\log|x|$ on every contour: cross-validated against the direct $\eta$-product
  ($3\times10^{-15}$), against $\Gamma_0(N)$-invariance ($10^{-16}$), and against
  the exact identity $\int_{|q|=r}\log|x|\,d\mu=\log r$.
* **Independent end-to-end recomputation of the $\zeta_5(3)$ contour.**  The
  contour lies well inside $|q|<1$, so the direct $\eta$-product converges and
  $\log|x|$ was recomputed at $40$ digits in `mpmath` with **no** modular
  reduction of any kind: agreement $4.9\times10^{-15}$ at every sampled point (on
  the degree-$16$ template, where the recomputed bound is $3.972717235$ with
  margin $+0.040441262$).  The identity $\log|\psi'(0)|=c_0=\overline u$ holds to
  machine precision.  §11.5 repeats the exercise on the final degree-$40$
  template with an a-priori error bound.
* $\mathrm{RE}$: the quadrature is exact to $10^{-8}$ on a closed-form test
  ($g=\cos2\pi t$, $\mathrm{RE}=4/\pi^2$), and every reported value is stable to
  $10^{-7}$ over $N=4096\dots65536$.
* $\mathrm{BC}$: exact on linear $\varphi$; but on a *degenerate* $\varphi$ (e.g.
  $z\mapsto z^2$, where $\varphi(z)=\varphi(-z)$ creates an antidiagonal
  $\log$-singularity) the error is $O(1/N)$ — $-0.013$ at $N=4000$.  None of our
  optimal contours shows this (all converge to $10^{-7}$), but **BC is the
  numerically fragile numerator and RE is the robust one**; all verdicts in §5
  are stated on RE.
* $R_p$ is *measured* from coefficient valuations (§2.4) only as a consistency
  check; the bound uses the *proved* lower bound from Buzzard's theorem, which is
  what the holonomy bound needs.

### (A) Assumed / not verified here — read before claiming anything

**A1. The multi-place form of the bound.**  *(Now settled as a citation gap —
see §11.4, which locates exactly what is and is not in L2chi.)*  The inequality actually used,
$$m\le\frac{\mathrm{BC}+\sum_v\log R_v}{\log|\varphi'(0)|+\sum_v\log R_v-\tau(\mathbf b)},$$
is quoted from CDT's ICM survey §6.2 eq. (6.2) (their $\zeta_2(5)$ proof), where
it is attributed to [L2chi, Thm 2.5.1] with the $p$-adic template
$\varphi_p(z)=R_p^{-1}z$.  The copy of L2chi in `papers/cdt/cdt2/L2chi.tex`
states Theorems 2.5.1/6.0.2/7.0.1 over $\mathbb Q$ at the archimedean place only;
we did **not** locate the multi-place statement (their §15.7) in that file.  We
take the ICM form at face value — it is calibrated exactly against their own
$\zeta_2(5)$ numbers in §3, which is strong evidence the transcription is right,
but it is an input, not a verification.

**A2. The rearrangement form of the multi-place bound.**  *(See §11.4.)*  Nazarov's
$\mathrm{BC}\le\mathrm{RE}$ is a theorem, and CDT prove a holonomy bound with the
rearrangement integral in place of BC with the *same* $\tau(\mathbf b)$
(Thm `main:elementary form`).  We assume the same substitution is legitimate in
the multi-place form.  If it is not, the $\zeta_5(3)$ verdict rests on the BC
computation alone (bound $3.9504<4$, margin $+0.0735$), which is still a
contradiction but is the numerically softer of the two.

**A3. Optimality of $\mathbf b$.**  We use CDT's array $\mathbf b=(0,d,\dots,d)$,
which §2.2 certifies as sharp.  The prime-to-$p$ refinement found there gives no
asymptotic improvement ($O(\log n)$), so $\tau(\mathbf b)$ is genuinely optimal
for this inventory.  What is *not* excluded is a better inventory (A4) or the
refined $\tau^{\sharp}(\mathbf e)$ machinery of L2chi §6 applied to a different
presentation of the same functions.

**A4. The function inventory.**  Only $\{1,H,\dots,H^{(2k)}\}$.  Adding
unconditional ("pure") functions of the same denominator type, as CDT do at
level $6$, would change both $m$ and $u_1$ and could rescue a near miss; none
were sought here.  Cf. `CDT_FINDER.md` §7.

**A5. Contour optimality, and a numerical trap.**  §4.1 is the exact search
space, but the optimisation is a finite-dimensional (cosine-polynomial,
$K\le24$) Nelder–Mead/Powell search from random restarts.  It gives a *lower*
bound on what is achievable: our margins can only improve.

The trap: the optimiser is run at $N=1024$–$2048$ quadrature nodes and is
therefore free to cheat by pushing the contour to $|q|=1-10^{-3}$ at a point that
is **not** a cusp, where $\log|x|\sim2\pi/(pc^{2}\operatorname{Im}\tau)$ has an
enormous spike of very small width — invisible to a coarse grid.  Such optima are
spurious, and they are easy to spot: they have $\max|q|\approx0.999$ and
$\max\log|x|$ in the thousands, and their value is *not* $N$-stable.  Two of our
restarts found them — $\zeta_7(3)$: apparent cost $0.826$ at $N=1024$, true cost
$3.955$ at $N=8192$; $\zeta_5(3)$: apparent cost $1.623$ ($\text{margin}+1.61$) at
$N=2048$, true cost $12.07$ at $N=4096$, $7.29$ at $N=8192$, with
$\max|q|=0.99867$ and $\max\log|x|=5484$.  **Every reported number was therefore recomputed at
$N=4096,8192,16384,32768$**, and only contours whose values are stable there (all
have $\max|q|\le0.95$) are quoted; the $\zeta_5(3)$ optimum is stable to $10^{-7}$
over $N=4096\dots65536$ and has $\max|q|=0.785$.  With that filter, four
independent searches for $\zeta_5(3)$ ($K=8,12,16,24$; $N=1024,2048$) all land at
margin $+0.040\pm0.002$, so the family is converged; a genuinely different family
could still do better.

**A6. Admissibility of the explicit templates.**  $\psi(z)=z e^{P(z)}$ maps
$\mathbb D$ into $\mathbb D$ iff $\max_{|z|=1}\operatorname{Re}P\le0$; we check
this on a $2\times10^5$-point grid with slack $\ge0.23$ (for $\zeta_5(3)$), not by
interval arithmetic.  Given the slack this is not a serious gap, but a formal
write-up would use the sufficient condition $c_0+\sum_{k\ge1}|c_k|\le0$ (which our
optimum narrowly fails: $-0.5335+0.5705=+0.037$) or a certified trigonometric
maximum.

**Not done.**  No literature audit (§8 flags what one must check).  No
interval-arithmetic certification of RE or BC.  No use of CDT's convexity/
multi-radius improvements of the numerator (§7–8 of L2chi), which improved their
level-$6$ numerator by $\approx2.5\%$ and would add roughly $+0.03$ to the
$\zeta_5(3)$ margin and $\approx+0.2$ to the $\zeta_2(5)$ margin.  No attempt at
$p=13$ (where $L=\log 13=2.565<\tau=2.8125$, so the entry condition fails
outright, exactly as for $\zeta_5(5)$).

---

## 10. Scripts

`lattice/padic_holonomy/`:

| file | what it does |
|---|---|
| `cdt_ab.gp`, `cdt_audit.gp`, `cdt_task1.gp`, `cdt_task23.gp` | exact $a_n,b_n$; integrality and sharpness of $2k+1$ (§2.1–2.3) |
| `cdt_task4.gp`, `cdt_task4b.gp`, `cdt_branch.gp`, `cdt_x09cusp.gp` | $p$-adic slopes ($R_p$), archimedean growth, exact branch points (§2.4–2.5) |
| `cdt_task5.gp`, `cdt_task5b.gp` | rank / $\mathbb Q(x)$-independence tables (§2.6) |
| `cdt_zetacheck.gp`, `cdt_Lcheck.gp`, `cdt_sharp.gp` | identification of $\eta$ against Kubota–Leopoldt; obstructing primes |
| `data_X0_p_k*.txt` | the exact coefficient data |
| `haupt.py` | hauptmoduln of $X_0(p),X_0(p^2),X_1(4)$ anywhere in $\lvert q\rvert <1$: $\Gamma_0(N)$-reduction by Gauss lattice reduction, and the branch-free $\log\lvert x\rvert $ via $\log\lvert \Delta\rvert $ and $\mathrm{SL}_2(\mathbb Z)$ |
| `family.py` | composable structured templates: SCALE, OFF, TANG, BITE (CDT lunes) |
| `outer.py` | the exact free parametrisation $\psi=z\cdot(\text{outer of }u)$; RE and BC quadratures |
| `bcint.py` | earlier structured-family versions of the same integrals |
| `targets.py` | the eleven cells with $d,m,\log R_p,\tau(\mathbf b)$, budget |
| `optimise.py`, `driver.py`, `scan1.py` | structured-family scans (`res_*.json`, `drv_*.log`) |
| `freeopt.py` | free-template optimisation (`free_*.json`, `free_*.log`) |
| `certify.py` | conversion to the explicit $\psi(z)=ze^{P(z)}$ form and admissibility check |
| `precise5.py` | high-precision 2-parameter optimum for $\zeta_5(3)$ (`precise5.log`) |
| `calib.py` | reproduction of CDT's two $\zeta_2(5)$ contours (`calib.log`) |
| `assemble.py` | final table, $\kappa$, contour profiles (`FINAL.json`, `assemble.log`) |
| `freeopt2.py`, `polish.py` | re-optimisation with the vectorised evaluator and the two-resolution anti-aliasing guard (§11.7) |
| `freeopt3.py` | the non-even (sine-terms) variant, §11.7 |
| `certnum.py`, `certrig.py` | §11.5: Lipschitz bound, $N$-refinement to $2^{24}$, and the fully a-priori certificate |
| `crude.py`…`crude5.py` | §11.8: the crude (sup-norm) multi-place bound — Dimitrov calibration, maximum-principle guard, and the family scans |
| — | §11.9: the algebraic identity between the CDT ICM form and `ADELIC_HOLONOMY.md` §2 (uniform slope), verified symbolically |
| `ap_lib.gp`, `ap_run2000.gp`, `ap_xcheck.gp` | §11.1: $a_n,b_n,e_n$ mod $5^{6150}$ to $n=2000$ by peeling, cross-checked against exact rationals to $n=600$ |
| `ap_kl.gp` | §11.1: $2\eta=\zeta_5(3)$ against Kubota–Leopoldt/Bernoulli |
| `ap_rank.gp`, `ap_rank2.gp`, `ap_ode.gp` | §11.3: rank to $D=120$ and the explicit order-3 ODE |
| `ap_lip.gp` | §11.5: $\max\lvert q\,d\log x/dq\rvert $ on circles |
| `ap_out_2000.txt` | the valuation data $\bigl(n,v_5(a_n),v_5(b_n),v_5(e_n),v_5(h_n)\bigr)$, $n\le2000$ |

---

## 11. Certification of the $\zeta_5(3)$ verdict

The margin in §5.2 is $1.3\%$ of the budget, so this section hardens every input
separately.  Everything below refers to $(p,k)=(5,1)$: $m=4$, $d=2k+1=3$,
$R_5=5^3$, $\tau(\mathbf b)=45/16$, $L=3\log5$.

### 11.1 $R_5=5^3$ read off the coefficients, to $n=2000$

**What is actually needed.**  The bound needs a *lower* bound $R_5\ge125$ for the
$5$-adic radius of $H(x)=\sum_n h_nx^n$, $h_n=a_n+\eta b_n$ — equivalently
$\liminf_n v_5(h_n)/n\ge3$.

**The computation** (`ap_lib.gp`, `ap_run2000.gp`; **exact**).  Since
$a_n,b_n,e_n\in\mathbb Z_5$ (§11.2), the whole construction can be run in
$\mathbb Z/5^{M}$ with $M=3N+150=6150$; residues are then exact reductions of the
true values, with no loss.  $x(q)$ is built from the pentagonal Euler product and
the re-expansion in $x$ is done by *peeling* (repeatedly take the constant term
and divide by $x$) rather than by series reversion.  Validation: the residues
agree with the exact-rational data of §2 for $n\le40$, and with a fully
independent exact-rational recomputation (`serreverse`+`subst` over $\mathbb Q$)
**coefficient for coefficient for all $n\le600$**; and $\eta$ computed from
$N=500$ and from $N=2000$ agree to $5^{1486}$, exactly the predicted accuracy of
the $N=500$ value.

**$\eta$.**  $v_5(b_N)=1$, so $\eta:=-a_N/b_N$ at $N=2000$ is pinned mod
$5^{6149}$:
$$\eta=1+2\cdot5+3\cdot5^2+3\cdot5^3+5^4+3\cdot5^5+4\cdot5^6+3\cdot5^7+2\cdot5^8+\cdots,\qquad v_5(\eta)=0 .$$

**The measurement** (`ap_out_2000.txt`).  With $h_n=a_n+\eta b_n$ and
$d_n:=v_5(h_n)-3n$ (only $n\le0.8N=1600$ is quoted; beyond that the answer is
contaminated by the choice of $\eta$):

| $n$ | 50 | 100 | 200 | 400 | 600 | 800 | 1000 | 1200 | 1400 | 1600 |
|---|---|---|---|---|---|---|---|---|---|---|
| $v_5(h_n)$ | 141 | 290 | 588 | 1189 | 1789 | 2387 | 2985 | 3585 | 4187 | 4785 |
| $v_5(h_n)/n$ | 2.8200 | 2.9000 | 2.9400 | 2.9725 | 2.9817 | 2.9838 | 2.9850 | 2.9875 | 2.9907 | 2.9906 |
| $d_n$ | $-9$ | $-10$ | $-12$ | $-11$ | $-11$ | $-13$ | $-15$ | $-15$ | $-13$ | $-15$ |

* $v_5(h_n)/n\to3$ cleanly, reaching $2.99063$ at $n=1600$.
* Least squares $d_n=A+B\log n$ gives $A=-0.974$, $B=-1.728$ over
  $50\le n\le1600$ ($A=-1.463$, $B=-1.654$ over $2\le n\le1600$).  Windowed
  minima ($-7$ on $[2,25]$, $-10$ on $[51,100]$, $-12$ on $[101,200]$, $-14$ on
  $[401,800]$, $-15$ on $[801,1600]$) show genuine $-C\log n$ behaviour with
  $C\approx1.8$, not a bounded constant — exactly the $O(\log n)$ deficit the
  theory predicts.
* $\max_n d_n=-3$ (at $n=2$): $d_n$ is **never positive**, so the radius is
  *exactly* $125$ and there is no slack to harvest.  ($\min d_n=-15$ at $n=937$.)

**The two factors separately — an instructive correction.**  One might expect
$E^*_{2}(x)$ to have radius $1$ and $E^*_{-2}(x)$ radius $125$.  **That is false.**
The data ($n\le2000$, exact) gives
$$v_5\bigl([x^n]E^*_2\bigr)/n\to0,\qquad v_5\bigl([x^n]E'_{-2}\bigr)/n\to0$$
— both $\le7$ in absolute terms over the whole range, both fitting $\approx0.27\log n$,
and $v_5(a_n)=v_5(b_n)$ for **every** $1\le n\le2000$.  So *both* naively
$x$-expanded factors have $5$-adic radius exactly $1$, and the radius $125$ is a
property of the **product alone**, produced by cancellation between $\eta b_n$ and
$\sum_{i<n}b_ie_{n-i}$.  There is no contradiction with non-archimedean
multiplicativity (the Gauss norm at $R=1$ is multiplicative and consistent;
"Newton polygon of $fg$ = sum of polygons" fails for series with unbounded
slopes).  The structural reason is that $E^*_2$ has weight $2$ and $E^*_{-2}$
weight $-2$: their bare $q\mapsto x$ re-expansions are not functions on the curve
and carry no intrinsic radius — only the weight-$0$ product does.  *This is why
the criterion is intrinsically about the pair $(a_n,b_n)$ and not about either
series alone.*

**$2\eta=\zeta_5(3)$, independently of the $a_n/b_n$** (`ap_kl.gp`, exact
rational Bernoulli numbers).  Using $\zeta_5(3)=\lim\zeta^*(1-2j)$ over
$2j\equiv-2\pmod{4\cdot5^r}$, $\zeta^*(1-2j)=-(1-5^{2j-1})B_{2j}/(2j)$:

| $r$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|
| $2j$ | 18 | 98 | 498 | 2498 | 12498 | 62498 | 312498 |
| $v_5\bigl(2\eta-\zeta^*(1-2j)\bigr)$ | 2 | 3 | 4 | 5 | 6 | 7 | 8 |

Agreement is exactly $r+1$ digits at every $r$ — the Kubota–Leopoldt continuity
rate — confirming $2\eta=\zeta_5(3)$ to $8$ $5$-adic digits, with no input from
the modular construction.

**Which theorem supplies $R_5\ge5^3$.**  K. Buzzard, *Analytic continuation of
overconvergent eigenforms*, J. Amer. Math. Soc. **16** (2003), 29–55,
**Theorem 5.2** (and Corollary 6.2 for the level-$p^2$ variants), applied exactly
as in Calegari's Lemma 3.1: $E^*_{-2k}$ is an overconvergent eigenform of level
$\Gamma_0(p)$ of **finite slope** — indeed slope $0$, since it is the weight-$(-2k)$
member of the *ordinary* Eisenstein family ($U_p$-eigenvalue $1$) — and Buzzard's
theorem continues such a form **across the entire supersingular annulus**, not
merely part of it (Calegari: *"such sections extend far into the supersingular
annuli.  In particular, by Theorem 5.2 of Buzzard, it extends entirely over the
supersingular annuli"*).  Together with the Fricke description of the ordinary
locus (§1.1) this is precisely the disc $|x|_5<5^{12/(p-1)}=5^3$.  Two caveats we
do not close: (i) Buzzard's paper is written for classical integer weights, and
the weight $-2k$ specialisation is via Coleman's theory / the eigencurve — this
is what Calegari invokes and we follow him; (ii) Calegari states the $p=2,3$ cases
and asserts in his §3.2 that "the same technique can also be applied to other $p$
when $X_0(p)$ has genus zero", which is the $p=5$ statement we use.  **The
measurement above is what makes this safe in practice: it shows the radius is
exactly $5^3$, so the cited lower bound is attained and not merely plausible.**

### 11.2 The denominator type: a proof, not a measurement

> **Proposition.**  Let $p$ be a prime with $p-1\mid24$, let
> $x=(\Delta(p\tau)/\Delta(\tau))^{1/(p-1)}$, and for $1\le k\le3$ put
> $$E^*_{2k}=E_{2k}(\tau)-p^{2k-1}E_{2k}(p\tau),\qquad
> E'_{-2k}=\sum_{n\ge1}\Bigl(\sum_{d\mid n,\ p\nmid d}d^{-2k-1}\Bigr)q^n,$$
> $H_\eta=E^*_{2k}(E'_{-2k}+\eta)=\sum_n(a_n+\eta b_n)x^n$.  Let
> $$D'_n:=\operatorname{lcm}\{\,d:1\le d\le n,\ p\nmid d\,\}$$
> be the **prime-to-$p$** lcm.  Then $b_n\in\mathbb Z$ and
> $(D'_n)^{2k+1}a_n\in\mathbb Z$ for every $n$; a fortiori
> $[1,\dots,n]^{2k+1}a_n\in\mathbb Z$.

*Proof.*  (i) $24/(p-1)\in\mathbb Z$, so
$x=q\prod_{n\ge1}\bigl((1-q^{pn})/(1-q^{n})\bigr)^{24/(p-1)}\in q+q^2\mathbb Z[\![q]\!]$.
A power series $q+\sum_{m\ge2}\alpha_mq^m$ with $\alpha_m\in\mathbb Z$ has compositional
inverse $q(x)=x+\sum_{m\ge2}\beta_mx^m$ with $\beta_m\in\mathbb Z$ (induction: $\beta_m$
is a universal polynomial with integer coefficients in $\alpha_2,\dots,\alpha_m$).

(ii) $E_{2k}=1-\frac{4k}{B_{2k}}\sum_{n\ge1}\sigma_{2k-1}(n)q^n$ and
$\frac{4k}{B_{2k}}=24,-240,504$ for $2k=2,4,6$, all integers; hence
$E_{2k}\in\mathbb Z[\![q]\!]$ and $E^*_{2k}\in\mathbb Z[\![q]\!]$ for $k\le3$.
Composing with $q(x)\in x+x^2\mathbb Z[\![x]\!]$ gives $E^*_{2k}\in\mathbb Z[\![x]\!]$,
i.e. $b_n\in\mathbb Z$.  (For $k\ge6$, $4k/B_{2k}\notin\mathbb Z$ and one rescales by a
*fixed* rational, which is harmless; it does not arise in this census.)

(iii) The $q^n$-coefficient of $E'_{-2k}$ is $c_n=\sum_{d\mid n,\,p\nmid d}d^{-2k-1}$.
Every $d$ occurring satisfies $1\le d\le n$ and $p\nmid d$, so $d\mid D'_n$ and
$(D'_n)^{2k+1}c_n\in\mathbb Z$.

(iv) Substituting $q=q(x)$: $E'_{-2k}(x)=\sum_{n\ge1}c_n\,q(x)^n$ with
$q(x)^n\in x^n\mathbb Z[\![x]\!]$, so $[x^N]E'_{-2k}(x)$ is a $\mathbb Z$-linear
combination of $c_1,\dots,c_N$, hence lies in $(D'_N)^{-(2k+1)}\mathbb Z$ because
$D'_n\mid D'_N$ for $n\le N$.  Multiplying by $E^*_{2k}(x)\in\mathbb Z[\![x]\!]$
preserves this, giving $(D'_N)^{2k+1}a_N\in\mathbb Z$. $\square$

This is exactly the Eichler-integral/$q$-expansion mechanism of the Theorem-B\*
integrality proofs elsewhere in this project: the only source of denominators is
the $d^{-(2k+1)}$ in the $(2k+1)$-fold formal primitive, and the change of
variable $q\mapsto x$ is integral in both directions.

**Sharpness** is the exact computation of §2.2: for $(p,k)=(5,1)$ the exponent $3$
cannot be lowered — $[1,\dots,n]^{2}a_n\notin\mathbb Z$ already at $n=3$, with
$\ell=3$ the obstructing prime.

**Use the Euler operator, not $d/dx$.**  The $x^{N}$-coefficient of $H^{(j)}$ is
$\frac{(N+j)!}{N!}h_{N+j}$, whose denominator is controlled by $D'_{N+j}$ rather
than $D'_N$ — a discrepancy of size $e^{O(\log N)}$, harmless for the *rate* but
outside the literal shape $a_{i,n}\in\mathbb Z$ demanded by the theorem.  Replacing
$d/dx$ by the Euler operator $\delta=x\,\tfrac{d}{dx}$ removes the issue exactly:
$$\delta^{j}H=\sum_n n^{j}h_n x^{n},\qquad (D'_n)^{2k+1}\,n^{j}h_n\in\mathbb Z .$$
Since $\delta^{j}=x^{j}\frac{d^{j}}{dx^{j}}+(\text{lower order, coefficients in }\mathbb Q[x])$
with invertible leading coefficient $x^j$,
$$\operatorname{span}_{\mathbb Q(x)}\{1,H,\delta H,\dots,\delta^{2k}H\}
=\operatorname{span}_{\mathbb Q(x)}\{1,H,H',\dots,H^{(2k)}\},$$
so the inventory, its size $m$, and its $\mathbb Q(x)$-independence are unchanged,
while the denominator type is now literally $[1,\dots,n]^{2k+1}$ for every member.
$\varphi^*\delta^jH$ is still meromorphic on $\mathbb D$ ($\delta=x(dq/dx)\,d/dq$).

**The array $\mathbf b$ and the $u_j$.**  Take $r=2k+1$ and
$$b_{1,j}=0\ (j=1,\dots,r),\qquad b_{i,j}=1\ (2\le i\le m,\ j=1,\dots,r),$$
i.e. denominators $[1,\dots,1\cdot n]^{\,r}=[1,\dots,n]^{2k+1}$ on the nose.  Every
column has the step shape with $u_j=1$ and $b_j=1$, so $\sigma_1=0$,
$\sigma_i=2k+1$ for $i\ge2$, and
$$\tau(\mathbf b)=\frac1{m^2}\sum_{i=1}^m(2i-1)\sigma_i
=\frac{2k+1}{m^2}\sum_{i=2}^{m}(2i-1)=\frac{2k+1}{m^2}(m^2-1)=(2k+1)\Bigl(1-\frac1{m^2}\Bigr),$$
equivalently $\sigma_m-\frac1{m^2}\sum_ju_j^2b_j=(2k+1)-\frac{2k+1}{m^2}$.
For $(p,k)=(5,1)$: $m=4$, $r=3$, $u=(1,1,1)$, $b=(1,1,1)$, $\sigma_m=3$,
$$\boxed{\ \tau(\mathbf b)=3\Bigl(1-\tfrac1{16}\Bigr)=\tfrac{45}{16}=2.8125\ \text{exactly}.}$$
CDT's presentation ($r=1$, $\mathbf b=(0,2k+1,\dots,2k+1)$, denominators
$[1,\dots,(2k+1)n]$) gives the *same* $\tau$ because only the growth rate enters;
the presentation above matches the actual denominators exactly and is the one we
use.  The prime-to-$p$ refinement gives no asymptotic gain (§2.2), so
$\tau(\mathbf b)=45/16$ is optimal for this inventory.

### 11.3 $\mathbb Q(x)$-independence of $\{1,H,\delta H,\delta^2H\}$

> **Lemma.**  Suppose $H$ satisfies an inhomogeneous linear ODE over $\mathbb Q(x)$
> of **exact minimal order** $N$ — i.e. there are $P_0,\dots,P_N,R\in\mathbb Q(x)$
> with $P_N\ne0$ and $\sum_{i\le N}P_iH^{(i)}=R$, and no such relation with a
> smaller order.  Then $1,H,H',\dots,H^{(N-1)}$ are $\mathbb Q(x)$-linearly
> independent.

*Proof.*  Let $\lambda_{-1}\cdot1+\sum_{i=0}^{N-1}\lambda_iH^{(i)}=0$ with
$\lambda_\bullet\in\mathbb Q(x)$.  If some $\lambda_i\ne0$ with $i\ge0$, take
$j=\max\{i\ge0:\lambda_i\ne0\}\le N-1$; then
$\sum_{i\le j}\lambda_iH^{(i)}=-\lambda_{-1}$ is an inhomogeneous ODE for $H$ of
order $j<N$ with nonzero leading coefficient — contradiction.  Hence all
$\lambda_i=0$ for $i\ge0$, and then $\lambda_{-1}=0$. $\square$

So the whole of the independence hypothesis reduces to **minimal order $=2k+1$**.
That is Kontsevich–Zagier (the proof of Fact 1 in their §2.3), quoted by CDT for
$\zeta_2(5)$.  The structural reason, in our case $k=1$, $p=5$:

* $\theta^{3}E'_{-2}=E^{\mathrm{evil}}_{4}=E_4(\tau)-E_4(5\tau)$, a *classical*
  holomorphic modular form of weight $4$ on $\Gamma_0(5)$
  ($\theta=q\,d/dq$).  So $E'_{-2}$ is an Eichler integral: a $3$-fold formal
  primitive of a weight-$4$ form.
* Consequently the monodromy of $E'_{-2}$ around the $x$-line adds, for each
  $\gamma\in\Gamma_0(5)$, a period polynomial of degree $\le2$ in $\tau$; after
  multiplication by the weight-$2$ form $E^*_{2}$, the monodromy orbit of $H$ is
  contained in the affine space $H+\operatorname{span}\{E^*_2,E^*_2\tau,E^*_2\tau^2\}$.
* The three functions $E^*_2\tau^{\,i}$ ($i=0,1,2$) are exactly a basis of
  solutions of $\operatorname{Sym}^2$ of the Picard–Fuchs (order-$2$) operator of
  $X_0(5)$, an order-$3$ operator over $\mathbb Q(x)$; hence $H$ satisfies an
  inhomogeneous ODE of order $\le3$.
* *Exactness* (order not $\le2$) is the non-degeneracy of the Eichler cocycle of
  $E^{\mathrm{evil}}_4$, i.e. that its period polynomial spans the full space of
  degree-$\le2$ polynomials — this is the Kontsevich–Zagier input, and it is what
  we take on citation.

**Independent verification** (exact, rank over $\mathbb F_q$, $q=2^{61}-1$;
`ap_rank.gp`, `ap_rank2.gp`, `ap_ode.gp`).  For $\eta=5/7$ and for $\eta=0$ alike:

| family | $D$ | required rank | actual | deficiency |
|---|---|---|---|---|
| $\{1,H,H',H''\}$ | $0,1,\dots,120$ | $4(D+1)$, up to $484$ | $4(D+1)$ | $0$ **everywhere** |
| $\{1,H,H',H'',H'''\}$ | $\le7$ | $5(D+1)$ | full | $0$ |
| $\{1,H,H',H'',H'''\}$ | $8$ | $45$ | $44$ | $\mathbf 1$ |
| $\{1,H,H',H'',H'''\}$ | $9,10,\dots$ | $5(D+1)$ | $-$ | $2,3,\dots$ (grows by exactly $1$) |

So: (a) **no $\mathbb Q(x)$-relation among $\{1,H,H',H''\}$ with
$\deg P_i\le120$** — far beyond any plausible degree; (b) an inhomogeneous ODE of
order exactly $3$ *does* exist, its kernel at $D=8$ is $1$-dimensional, and the
polynomial coefficients of $(1,H,H',H'',H''')$ in it have degrees
$(5,4,6,7,8)$ — the $H'''$ coefficient and the inhomogeneous term are both
nonzero, so it is a genuine order-$3$ inhomogeneous relation; (c) the deficiency
incrementing by exactly $1$ per step confirms that relation is unique.

This is a check, not a substitute for Kontsevich–Zagier: a rank certificate to
degree $D$ cannot by itself exclude a relation of larger degree.  But the
juxtaposition is striking — the minimal order-$3$ relation has degree $8$, while
a hypothetical order-$\le2$ relation would have to have degree $>120$.

### 11.4 The multi-place bound: exactly what is being used, and where it lives

The inequality used is
$$m\;\le\;\frac{\displaystyle\sum_{v}\iint_{\mathbb T^2}\log|\varphi_v(z)-\varphi_v(w)|_v\,\mathrm d\mu\,\mathrm d\mu}
{\displaystyle\sum_{v}\log|\varphi_v'(0)|_v-\tau(\mathbf b)}
\;=\;\frac{\mathrm{BC}(\varphi_\infty)+\sum_p\log R_p}{\log|\varphi_\infty'(0)|+\sum_p\log R_p-\tau(\mathbf b)},$$
where the second equality is the specialisation $\varphi_p(z)=R_p^{-1}z$: then
$\log|\varphi_p'(0)|_p=\log R_p$, and
$$\iint_{\mathbb T^2}\log|\varphi_p(z)-\varphi_p(w)|_p=\log R_p+\iint_{\mathbb T^2}\log|z-w|_p=\log R_p+0 ,$$
so **$\sum_p\log R_p$ enters numerator and denominator once each** — which is
precisely CDT's ICM eq. (6.2) for $\zeta_2(5)$, and is confirmed empirically: with
$\mathrm{BC}=2.13322$, $|\varphi'(0)|=1/3$, $\tau=175/36$ and $\log R_2=12\log2$
it reproduces their $4.43206$ to $5\times10^{-5}$ (§3).  Any other placement of
$\log R_2$ fails to reproduce it.

**Where the statement lives.**  We searched `papers/cdt/cdt2/L2chi.tex`
(arXiv:2408.15403) for the $p$-adic-place version and it is **not there**:

* Theorems 2.5.1 / 6.0.2 / 7.0.1 (`basic main`, `main:elementary form`,
  `main:BC form`) are stated over $\mathbb Q$ with a **single archimedean**
  $\varphi$.
* Remark `BCboundK` (§ after Thm `main:BC form`) gives the number-field
  generalisation, with $\varphi_\sigma$ one per **complex embedding**
  $\sigma:K\hookrightarrow\mathbb C$:
  $$m\le\frac{\sum_\sigma\iint\log|\varphi_\sigma(z)-\varphi_\sigma(w)|}
  {\bigl(\sum_\sigma\log|\varphi_\sigma'(0)|\bigr)-[K:\mathbb Q]\bigl(\tau(\mathbf b)+\tau^\sharp(\mathbf e)\bigr)} .$$
  This is the exact shape displayed above — a sum over places upstairs and
  downstairs — but only over **archimedean** places.
* L2chi says so explicitly in its concluding section: *"in this paper, we have
  emphasized the archimedean place as special when it comes to overconvergence.
  In [CDT24] we plan to write our holonomicity bound in a more general Arakelov
  adelic form over a global field."*  Here **[CDT24] $=$ F. Calegari,
  V. Dimitrov, Y. Tang, _Arithmetic holonomy bounds and the irrationality of the
  $2$-adic $\zeta(5)$_, 2024** — the companion treatise, which is **not in this
  repository**.

> **Therefore:** the inequality we use is stated here as *"the CDT ICM §6.2 form,
> the one they themselves use for $\zeta_2(5)$"*, proved in the companion paper
> [CDT24], and **not verified against a proof available to us**.  *(Superseded by
> §11.9(e): for the Bost–Charles numerator the inequality is now proved, from
> published ingredients plus one new elementary lemma.)*  §11.8 completes
> the picture: the *crude* (sup-norm) multi-place bound **is** published and
> available (Dimitrov, arXiv:1912.12545, stated in his LNT notes with $p$-adic
> places included), and it suffices for $\zeta_2(5)$ — but **not** for
> $\zeta_5(3)$, which needs the refined form.  Its shape is
> forced by, and consistent with, L2chi Remark `BCboundK`, and its numerical
> transcription is confirmed by reproducing their $4.43206$; but this remains the
> single load-bearing citation of the whole census.  The same caveat applies a
> fortiori to the **rearrangement** variant of the multi-place bound (L2chi's
> Thm `main:elementary form` is the archimedean rearrangement bound; we assume
> the multi-place version substitutes $\mathrm{RE}$ for $\mathrm{BC}$ with the
> same $\tau(\mathbf b)$, which is legitimate at least in the direction
> $\mathrm{BC}\le\mathrm{RE}$ of Nazarov's inequality).

### 11.5 The archimedean numerics, with an explicit error bound

**The certified template.**  $\psi(z)=z\exp\bigl(c_0+\sum_{k=1}^{40}c_kz^k\bigr)$,
$c_k\in\mathbb R$ (`pol_05_K40.json`):

$$\begin{aligned}
c_{0\ldots5}&=(-0.531289158,\,-0.225379074,\,-0.127598877,\,-0.001824369,\,-0.027735780,\,-0.041418343)\\
c_{6\ldots11}&=(+0.051444137,\,+0.026087140,\,-0.010853938,\,+0.001314954,\,-0.003547849,\,-0.023294126)\\
c_{12\ldots17}&=(+0.000125751,\,+0.020158067,\,-0.001454440,\,-0.004262361,\,+0.007019143,\,-0.000198726)\\
c_{18\ldots23}&=(-0.007684408,\,-0.006339902,\,+0.005397468,\,+0.001260337,\,+0.000134428,\,+0.014051878)\\
c_{24\ldots29}&=(-0.011391873,\,-0.005386029,\,-0.000202917,\,+0.001464054,\,-0.000123972,\,+0.008670335)\\
c_{30\ldots35}&=(-0.000051833,\,-0.006877094,\,-0.001450246,\,-0.001060183,\,+0.002930635,\,+0.002103211)\\
c_{36\ldots40}&=(+0.001750273,\,-0.002516914,\,-0.002244428,\,+0.000524421,\,+0.001281738)
\end{aligned}$$

**(a) Admissibility.**  On $|z|=1$, $\log|\psi|=u(t):=c_0+\sum_kc_k\cos2\pi kt$, a
trigonometric polynomial with exact Lipschitz constant $2\pi\sum_kk|c_k|=2\pi\cdot4.3589$.
Sampling on $2^{22}$ points and adding $2\pi\sum k|c_k|/2^{23}$ gives
$$\max_t u(t)\;\le\;-0.17880 .$$
Hence $|\psi(z)|\le e^{u}\le e^{-0.1788}<1$ on $|z|=1$, so $\psi(\mathbb D)\subset\mathbb D$
by the maximum principle, and $\varphi=x\circ\psi$ is admissible (§1.5).  The
slack is $0.179$, four orders of magnitude above any plausible evaluation error;
the cheaper sufficient test $c_0+\sum_{k\ge1}|c_k|\le0$ is *not* satisfied
($=+0.127$), so the grid-plus-Lipschitz certificate is the operative one.

**(b) $\log|\varphi'(0)|$ is exact.**  $x=q+O(q^2)$ so $x'(0)=1$, and
$\psi'(0)=e^{c_0}$; hence
$$\log|\varphi'(0)|=c_0=-0.531289158096\quad\text{exactly (a decimal we chose).}$$

**(c) A Lipschitz bound for $g(t)=\log|x(\psi(e^{2\pi it}))|$.**  With
$\Lambda(q):=q\frac{d}{dq}\log x=1+6\sum_{n\ge1}n\Bigl(\frac{q^n}{1-q^n}-\frac{5q^{5n}}{1-q^{5n}}\Bigr)$
(from $x=q\prod((1-q^{5n})/(1-q^n))^6$; verified against $d g/d\theta=-\operatorname{Im}\Lambda$
on circles), and $q=z\,e^{P(z)}$ giving $\frac{dq}{dt}=2\pi i\,q\,(1+zP'(z))$,
$$g'(t)=\operatorname{Re}\bigl[\,2\pi i\,\Lambda(q)\,(1+zP'(z))\,\bigr],\qquad
|g'|\le 2\pi\max|\Lambda|\cdot\Bigl(1+\sum_kk|c_k|\Bigr).$$
We use **two** bounds for $\max|\Lambda|$:

* an **a-priori** one, from the triangle inequality on the Lambert series
  (`ap_lip.gp` verifies the $q$-series identity
  $\Lambda=-\tfrac14(E_2(q)-5E_2(q^5))=1+6q+18q^2+24q^3+\cdots$ to $O(q^{61})$):
  $$\max_{|q|\le r}|\Lambda(q)|\;\le\;M(r):=1+6\sum_{n\ge1}n\Bigl(\frac{r^n}{1-r^n}+\frac{5r^{5n}}{1-r^{5n}}\Bigr),
  \qquad M(0.836265)=339.3697 ,$$
  giving $\operatorname{Lip}(g)\le2\pi\cdot339.3697\cdot(1+4.358941)=\mathbf{11427.0}$
  — **no sampled quantity enters this**;
* the **true** $\max_t|g'(t)|$ on the contour, sampled on $2^{20}$ points:
  $\operatorname{Lip}(g)=259.68$.  (It is $44\times$ smaller because the contour
  reaches $|q|=0.836$ only at angles far from $\theta=0$, where $\Lambda$ peaks:
  $\max|\Lambda|$ *on the contour* is only $15.25$, against $339$ on the full
  circle of that radius.)

**(d) The quadrature bound.**  $\mathrm{RE}_N=\frac1{N^2}\sum_{i,j}\max(g_i,g_j)$
(identically equal to the sorted formula $\sum_ig_{(i)}(2i-1)/N^2$) is a
left-endpoint Riemann sum on the torus of
$F(s,t)=\max(g(s),g(t))$, which is $\operatorname{Lip}(g)$-Lipschitz in the
sup-metric.  Hence
$$\bigl|\mathrm{RE}-\mathrm{RE}_N\bigr|\;\le\;\frac{\operatorname{Lip}(g)}{N} .$$
$N$-refinement (`certnum.py`):

| $N$ | $\mathrm{RE}_N$ | $\mathrm{RE}_N-\mathrm{RE}_{N/2}$ | $\operatorname{Lip}/N$ | bound$_N$ |
|---|---|---|---|---|
| $2^{12}$ | 1.066964686390 | — | $6.34\cdot10^{-2}$ | 3.971155821 |
| $2^{13}$ | 1.066964916829 | $2.30\cdot10^{-7}$ | $3.17\cdot10^{-2}$ | 3.971155976 |
| $2^{14}$ | 1.066964976643 | $5.98\cdot10^{-8}$ | $1.59\cdot10^{-2}$ | 3.971156016 |
| $2^{15}$ | 1.066964991126 | $1.45\cdot10^{-8}$ | $7.93\cdot10^{-3}$ | 3.971156026 |
| $2^{16}$ | 1.066964994868 | $3.74\cdot10^{-9}$ | $3.96\cdot10^{-3}$ | 3.971156028 |
| $2^{17}$ | 1.066964995801 | $9.33\cdot10^{-10}$ | $1.98\cdot10^{-3}$ | 3.971156029 |
| $2^{18}$ | 1.066964996037 | $2.36\cdot10^{-10}$ | $9.91\cdot10^{-4}$ | 3.971156029 |
| $2^{19}$ | 1.066964996096 | $5.86\cdot10^{-11}$ | $4.95\cdot10^{-4}$ | 3.971156029 |
| $2^{20}$ | 1.066964996111 | $1.46\cdot10^{-11}$ | $2.48\cdot10^{-4}$ | 3.971156029 |

The successive differences fall by a factor of exactly $4$ per doubling — the
$O(N^{-2})$ behaviour expected for a smooth periodic integrand with a
measure-zero crease — so Richardson gives
$\mathrm{RE}=1.06696499612\pm10^{-11}$.  We do **not** use that.

**(e) The fully a-priori certificate** (`certrig.py`), using only the majorant
$\operatorname{Lip}(g)\le11427.0$ of (c), i.e. with *no sampled quantity in the
error bound at all*:

| $N$ | $\mathrm{RE}_N$ | $\varepsilon=\operatorname{Lip}/N$ | bound $\le$ | margin $\ge$ |
|---|---|---|---|---|
| $2^{20}$ | 1.066964996111 | $1.090\cdot10^{-2}$ | 3.978496850 | $+0.031921955$ |
| $2^{21}$ | 1.066964996114 | $5.449\cdot10^{-3}$ | 3.974826439 | $+0.037370769$ |
| $2^{22}$ | 1.066964996115 | $2.724\cdot10^{-3}$ | 3.972991234 | $+0.040095176$ |
| $2^{23}$ | 1.066964996115 | $1.362\cdot10^{-3}$ | 3.972073632 | $+0.041457380$ |
| $2^{24}$ | 1.066964996116 | $6.811\cdot10^{-4}$ | **3.971614831** | $\mathbf{+0.042138482}$ |

Note the verdict already holds at $N=2^{20}$ with the crudest a-priori chain
(bound $3.9785<4$).  Using instead the sampled $\operatorname{Lip}(g)=259.68$ at
$N=2^{20}$ gives $\varepsilon=2.476\times10^{-4}$ and bound $\le3.9713229$.

**(f) The certified inequality.**  With $L=3\log5=4.828313737$ (exact),
$\tau(\mathbf b)=45/16$ (exact, §11.2), $c_0=-0.531289158096$ (exact),
$$\boxed{\ m\;\le\;\frac{\mathrm{RE}+L}{c_0+L-\tau(\mathbf b)}
\;\le\;\frac{1.067646097+4.828313737}{1.484524579}\;=\;3.9716148\;<\;4\;=\;m\ }$$
$$\textbf{margin}\;=\;m\bigl(c_0+L-\tau(\mathbf b)\bigr)-(\mathrm{RE}+L)\;\ge\;+0.0421385 ,$$
using the **a-priori** chain at $N=2^{24}$ (the sampled-Lipschitz variant at
$N=2^{20}$ gives the slightly better $3.9713229$ and $+0.0425719$).
With the (sharper, $N$-converged to $10^{-8}$) Bost–Charles numerator
$\mathrm{BC}=1.03552928$ the same template gives bound $=3.9499804$ and margin
$=+0.0742553$.

**Irrationality measure.**  $\kappa=(1-\gamma)/\bigl(1-\sqrt{1-\text{margin}/(mL)}\bigr)$
with $\gamma=1/4$: the certified RE margin $+0.0421385$ gives $\kappa\le687.0$;
the (uncertified but $N$-converged) BC margin $+0.0742553$ gives $\kappa\le389.8$.

**What "certified" does and does not mean here.**  (i) The error bound
$\varepsilon=\operatorname{Lip}(g)/N$ is a genuine a-priori quadrature bound, and
in the boxed statement $\operatorname{Lip}(g)$ is itself a-priori (a
triangle-inequality majorant on $|q|\le\max|q|$), so **no sampled quantity enters
the error term**.  (ii) $\max_t u\le-0.1788$ is a $2^{22}$-point grid maximum
plus the exact Lipschitz constant of a trigonometric polynomial — rigorous up to
float64.  (iii) All arithmetic is float64; the quantities are $O(1)$ with
$\sim10^{-11}$ round-off, seven orders below $\varepsilon$, and the evaluation of
$\log|x|$ was cross-checked against a $40$-digit `mpmath` recomputation with no
modular reduction (agreement $4.9\times10^{-15}$).  (iv) What is *not* done: a
genuine interval-arithmetic enclosure of $\log|\Delta|$, the
$\mathrm{SL}_2(\mathbb Z)$ reduction, the FFT and the sort.  Given that the
verdict survives even the crudest chain at $N=2^{20}$ by a factor of $4$ in
$\varepsilon$, this is a formality rather than a risk — but it is not done.

### 11.6 Sensitivity: what the verdict actually hinges on

The margin is linear in each input, so the sensitivity is exact.  With
$\text{margin}=m(c_0+L-\tau)-(\mathrm{RE}+L)$, $m=4$, $L=4.82831$, $\tau=45/16$:

| perturbation | effect on the margin | tolerance |
|---|---|---|
| $\mathrm{RE}$ up by $\varepsilon$ | $-\varepsilon$ | survives to $\varepsilon=+0.042$, i.e. a $4\%$ error in RE |
| $\log R_5$ down by $\varepsilon$ | $-(m-1)\varepsilon=-3\varepsilon$ | survives to $\varepsilon=0.014$, i.e. $R_5\ge5^{2.991}$ |
| $\tau(\mathbf b)$ up by $\varepsilon$ | $-4\varepsilon$ | survives to $\varepsilon=0.0106$: the exponent $2k+1=3$ must be right to $0.35\%$ |
| $\log\lvert\varphi'(0)\rvert$ down by $\varepsilon$ | $-4\varepsilon$ | survives to $\varepsilon=0.0106$ |
| $m\to3$ (an order-$2$ inhomogeneous ODE) | $\tau\to\tfrac83$, margin $\to-1.004$ | **fails** |
| $m\to5$ (one extra independent function, same type) | $\tau\to\tfrac{72}{25}$, margin $\to+1.185$ | comfortable |
| $m\to6$ | $\tau\to\tfrac{35}{12}$, margin $\to+2.380$ | comfortable |

So there are exactly three ways the verdict could fail, in decreasing order of
worry:

1. **$m=4$ wrong downward** — $\{1,H,\delta H,\delta^2H\}$ $\mathbb Q(x)$-dependent,
   i.e. minimal inhomogeneous ODE of order $\le2$.  Excluded by Kontsevich–Zagier
   (§11.3) and, to the degrees scanned, by the rank certificate.  An error in the
   *other* direction only strengthens the conclusion: each extra independent
   function of the same denominator type is worth about $+1.2$ of margin.
2. **The multi-place form of the bound** (§11.4) — cited, not verified here.
3. **The numerics** — bounded explicitly in §11.5, at $\sim0.3\%$ of the margin.

The archimedean template is the *least* critical input: it is a lower bound by
construction (a better contour only helps), and even a plain disc $|q|\le r$ at
the optimal $r=0.5965$ already gives a positive margin, $+0.0053$.

### 11.7 Improving the margin, and the one large lever not pulled

Independent re-optimisations, all with the anti-aliasing guard of §9 A5 (the
objective is evaluated at both $N$ and $2N$ and the *larger* cost used):

| run | family | $K$ | $N$ | cost | margin |
|---|---|---|---|---|---|
| `precise5` | $\mathrm{SCALE}(r)\circ\mathrm{BITE}(0,c)$, 2 parameters | — | 32000 | $3.233882$ | $+0.001060$ |
| `free_05` | even cosine profile | 14 | 1024 | $3.194500$ | $+0.040441$ |
| `f2_05_a` | even, guarded, $u_{\min}=10^{-3}$ | 16 | 4096 | $3.194244$ | $+0.040697$ |
| `f2_05_e` | even, guarded, $u_{\min}=10^{-4}$ ($\lvert q\rvert\le0.9999$ allowed) | 24 | 8192 | $3.193693$ | $+0.041249$ |
| `f2_05_d` | even, guarded, $u_{\min}=10^{-2}$ ($\lvert q\rvert\le0.99$) | 24 | 8192 | $3.194328$ | $+0.040614$ |
| `f3_05_a` | **non-even** (sine terms allowed) | 12 | 4096 | $3.195493$ | $+0.039449$ |
| `pol_05_K16` | direct polish in the $c_k$ | 16 | 16384 | $3.193916$ | $+0.041025$ |
| `pol_05_K24` | " | 24 | 16384 | $3.193078$ | $+0.041863$ |
| `pol_05_K32` | " | 32 | 16384 | $3.192493$ | $+0.042448$ |
| `pol_05_K40` | " | 40 | 16384 | $3.192703$ | $+0.042239$ |

Three conclusions.  (i) **Non-even profiles do not help** — expected, since
$|x(\bar q)|=|x(q)|$; the optimum is genuinely even.  (ii) **Letting the contour
approach $|q|=1$ does not help**: $u_{\min}=10^{-2}$ ($|q|\le0.99$) and
$10^{-4}$ ($|q|\le0.9999$) give the same answer as $10^{-3}$, and the optimum has
$\max|q|=0.785$ — the constraint is *inactive*, the optimum is an interior
stationary point.  Pushing to $|q|\le0.9$ or beyond is therefore not a lever, and
the coordinator's alternative is ruled out numerically.  (iii) Raising the degree
$K$ from $8$ to $40$ buys about $+0.007$ of margin in total and is saturating; the
best stable value is $\text{margin}=+0.0424$ at $K=32$.

**The one large lever not pulled** is CDT's *convexity / multi-radius* improvement
of the numerator (L2chi §§7–8, Thms `main:BC conv discrete`, `main:BC fullconv`),
which replaces $\mathrm{BC}(\varphi)$ by a strictly smaller quantity built from
$\varphi(rz)$ at several radii.  On their own level-$6$ application it moved the
numerator by $2.3\%$ ($13.9938\to13.621$).  Here the numerator is
$\mathrm{RE}+L=5.895$, of which the *contour-dependent* part is
$\mathrm{RE}=1.067$; a $20\%$ reduction of that part would add $\approx+0.21$ to
the margin, a fivefold improvement.  We did not implement it: doing so correctly
needs the multi-place version of those theorems, which is subject to the same
citation gap as §11.4.

Finally, the **Bost–Charles numerator is the sharper of the two** and is the one
CDT actually use: with it the margin is $+0.073$ rather than $+0.042$ and the
bound $3.950$ rather than $3.972$.  The headline uses the conservative
rearrangement numerator.

### 11.8 The crude (sup-norm) multi-place bound, and the exact citation chain

#### (a) The two multi-place statements, and where each lives

There are **two** multi-place holonomicity bounds in play, and the $\zeta_5(3)$
claim rests on the second, not the first.

**(I) The crude form — Dimitrov.**  Stated explicitly in V. Dimitrov's LNT
lecture notes (`papers/dimitrov_LNT_zeta2_5.txt`, "The arithmetic holonomicity
theorem"), over a number field $K$ with a radius $R_v>0$ and a normalised
holomorphic $x_v(z)=z+O(z^2)$ at **every** place $v\in M_K$ (with $x_v(z)=z$ for
all but finitely many $v$), and $f(x_v(z))$ the germ of a $v$-adic meromorphic
function on $D(0,R_v)$:

> **Theorem.**  If $\sum_{v\in M_K}\log R_v>\tau_K(f)$ then $f$ is holonomic, and
> with $S_v:=\sup_{|z|_v=R_v}|x_v(z)|_v$ there is an inhomogeneous equation
> $L(f)\in K[x]$ with
> $$r\;\le\;\frac{\sum_v\log^{+}S_v}{\sum_v\log R_v-\tau_K(f)} .$$

Proved in **V. Dimitrov, arXiv:1912.12545** (the Schinzel–Zassenhaus paper), of
which the $K=\mathbb Q$, $R_p=1$, $x_p(z)=z$ case is the "gist of
Schinzel–Zassenhaus" specialisation.  $p$-adic places are included from the
outset.  **This is a published, available statement.**

**(II) The refined form — CDT.**  The same shape with $\log^{+}S_\infty$ replaced
by the Bost–Charles double integral (or the rearrangement integral) and
$\tau_K(f)=2k+1$ replaced by the array-refined $\tau(\mathbf b)$:
$$m\;\le\;\frac{\mathrm{BC}(\varphi_\infty)+\sum_p\log R_p}
{\log|\varphi_\infty'(0)|+\sum_p\log R_p-\tau(\mathbf b)} .$$
This is CDT's ICM survey **arXiv:2510.04156, §6.2 eq. (6.2)** (their `fullforce`),
used by them exactly there for $\zeta_2(5)$.  The ICM survey **states and uses it
but does not restate or prove it**; it attributes the archimedean ingredients to
arXiv:2408.15403 (`L2chi`) and the whole package to the companion treatise
**[CDT24] = Calegari–Dimitrov–Tang, _Arithmetic holonomy bounds and the
irrationality of the $2$-adic $\zeta(5)$_ (2024)**.  As established in §11.4,
`L2chi` contains only the archimedean and the multi-*embedding* number-field
forms and explicitly defers the adelic form to [CDT24].

**Calibration of (I).**  Dimitrov's own $\zeta_2(5)$ application: $R_2=2^{12}$,
$x_2(z)=z$; $R_\infty=1/5$ with
$x_\infty(z)=x\bigl(z/(1+3z)\bigr)$, i.e. $\Omega=B=\{|q+\tfrac3{16}|\le\tfrac5{16}\}$;
$\tau=5$, $m=6$.  We reproduce
$$\sup_{\partial B}|x|=|x(1/8)|=3.23160,\qquad
r\le\frac{12\log2+\log3.2316}{12\log2+\log(1/5)-5}=\mathbf{5.555573}<6$$
(his printed "$<9.5/1.7=5.58$"), i.e. crude cost $10.8296<$ crude budget
$5L-30=11.5888$ — a contradiction with margin $+0.759$.  So the **crude form
alone suffices for $\zeta_2(5)$**, which is why the notes can present the whole
proof without the refined machinery.

#### (b) The crude form applied to $\zeta_5(3)$: it fails

Here $m=4$, $\tau=2k+1=3$ (not $45/16$), $L=\log R_5=3\log5$.  Contradiction
requires
$$\log^{+}S_\infty-4\log R_\infty\;<\;(m-1)L-m\tau=3L-12=\mathbf{2.4849412}.$$

**A trap that must be excluded first.**  The crude form needs
$S=\sup_{\overline{\mathbb D}}|x\circ\psi|$ — the $H^\infty$ norm — not merely the
sampled boundary values.  A disc internally tangent to $|q|=1$ at an image of the
**cusp $0$** looks spectacular and is worthless: for $\Omega=D(0.17,0.83)$
(tangent at $q=+1$) the boundary samples give $\max\log|x|=-0.027$ while the
interior contains $q=0.9$ with $\log|x(0.9)|=70.1$, and a grid over $\Omega$
reaches $\log|x|=6271.3$.  The reason is that the tangential approach to the cusp
$0$ lands on a horocycle at finite height, so the a.e. boundary values stay
bounded, while $x\circ\psi$ blows up radially: $x\circ\psi$ is unbounded on
$\mathbb D$ with bounded boundary values — singular-inner-function behaviour — so
it is not holomorphic on $\overline{\mathbb D}$ and violates the hypothesis of
*both* forms.  **Tangency to $|q|=1$ is admissible only at images of the cusp
$\infty$** (where $x\to0$), as in CDT's $\Omega_{\rm circ}$; §4.3 above is
correct, and every scan below carries a maximum-principle guard that rejects such
contours (`crude3.py`).

**The optimisation** (`crude2.py`, `crude3.py`, `crude4.py`, `crude5.py`).  Best
crude cost over each family, all with the guard and all $N$-stable:

| family | best | $\log R$ | $\log^{+}S$ | crude cost | crude bound |
|---|---|---|---|---|---|
| centred disc $\lvert q\rvert\le r$ | $r=0.240$ | $-1.42712$ | $+0.68043$ | $6.38889$ | $9.24$ |
| off-centre disc $D(-a,b)$ | $a=0.23,b=0.44$ | $-1.14015$ | $+0.19245$ | $4.75303$ | $6.84$ |
| off-centre disc $D(+a,b)$ | $a=0.01,b=0.24$ | $-1.42885$ | $+0.84695$ | $6.56237$ | $9.42$ |
| tangent at $\zeta_5^{\,j}$ (cusp $\infty$) | $j=2,b=0.565$ | $-1.46929$ | $+8.57280$ | $14.44996$ | $27.7$ |
| free even profile ($K=10$) | — | $-0.88494$ | $+0.63742$ | $4.17718$ | $5.79381$ |
| gobble (two bites) | $r{=}0.7,c_1{=}2,c_2{=}20$ | $-0.87250$ | $+0.57061$ | $4.06061$ | $5.75$ |
| **lune $\mathrm{SCALE}(r)\circ\mathrm{BITE}(0,c)$** | $r=0.680698$, $c=1.948571$ | $\mathbf{-0.924090}$ | $\mathbf{+0.338722}$ | $\mathbf{4.035083}$ | $\mathbf{5.714335}$ |

$$\boxed{\ \text{best crude bound for }\zeta_5(3):\quad m\;\le\;5.714335\;\not<\;4=m\ }$$

crude cost $4.035083$ against a crude budget of $2.484941$: **it fails by
$1.550$**.  (Conformal radius $R_\infty=0.39689$, $S_\infty=1.40315$; the contour
is the disc $|q|\le0.680698$ with a single lune-bite of parameter $c=1.948571$
taken out towards $q=+1$, reaching $\max|q|=0.68070$.)  To succeed one would need
$S$ smaller by a factor $e^{1.55}=4.7$ at the same $R$, or $R$ larger by
$e^{0.39}=1.48$ at the same $S$ — neither is remotely available: the free-template
optimisation over all admissible $\psi$ does *worse* ($4.177$) than the
two-parameter lune, so the sup-norm optimum is essentially attained.

#### (c) Exactly which refinement carries the claim

Both numerator and denominator refinements are needed.  Evaluating all four
combinations on **each** contour (each optimised for its own objective):

| contour | numerator | $\tau$ | bound | margin |
|---|---|---|---|---|
| crude-optimal ($\log R=-0.92409$) | $\log^{+}S=0.33872$ | $3$ | $5.714331$ | $-1.550138$ |
| " | $\log^{+}S$ | $45/16$ | $4.732913$ | $-0.800138$ |
| " | $\mathrm{RE}=-0.10362$ | $3$ | $5.225144$ | $-1.107804$ |
| " | $\mathrm{RE}$ | $45/16$ | $4.327742$ | $-0.357804$ |
| refined-optimal ($\log R=-0.53129$) | $\log^{+}S=3.95870$ | $3$ | $6.774750$ | $-3.598918$ |
| " | $\log^{+}S$ | $45/16$ | $5.919078$ | $-2.848918$ |
| " | $\mathrm{RE}=1.06696$ | $3$ | $4.545233$ | $-0.707180$ |
| **"** | $\mathbf{\mathrm{RE}}$ | $\mathbf{45/16}$ | $\mathbf{3.971156}$ | $\mathbf{+0.042820}$ |

* The $\tau$ refinement $3\to\frac{45}{16}$ is worth exactly
  $m\cdot\frac3{16}=0.75$ of margin, on any contour.
* The numerator refinement $\log^{+}S\to\mathrm{RE}$ is worth $2.892$ on the
  refined-optimal contour ($3.9587\to1.0670$) — it is the **decisive** one.
* Neither alone suffices: with the crude numerator and the refined $\tau$ the best
  margin is $-0.800$; with the refined numerator and $\tau=3$ it is $-0.707$.

$$\textbf{Therefore the }\zeta_5(3)\textbf{ claim rests specifically on the refined
multi-place bound}$$
— the form with a Bost–Charles/rearrangement archimedean numerator *and* the
array-refined $\tau(\mathbf b)$, at $p$-adic as well as archimedean places.  That
is CDT ICM arXiv:2510.04156 §6.2 eq. (6.2), **used there but not restated or
proved there**, and attributed to the companion [CDT24].  Dimitrov's published
crude form (arXiv:1912.12545) — which is enough for $\zeta_2(5)$ — is **not**
enough for $\zeta_5(3)$: it gives $m\le5.714$ where $m=4$ is needed.

**The single published statement the proof needs** is therefore the refined
multi-place (adelic) holonomy bound of [CDT24].  Until that is in hand, the
$\zeta_5(3)$ computation should be read as: *"the CDT ICM §6.2 inequality, applied
to $X_0(5)$ with $k=1$, is violated"* — a conditional statement whose one
unverified hypothesis is that inequality itself.

### 11.9 The dependency chain, and the theorem as it would be written

#### (a) The algebraic reduction, verified

The refined multi-place form is *identically* an assertion about a single
weighted combination:
$$\frac{\mathrm{BC}+\sum_p\log R_p}{\log|\varphi'(0)|+\sum_p\log R_p-\tau(\mathbf b)}\;\ge\;m
\iff
m\Bigl(\log|\varphi'(0)|-\tau(\mathbf b)+\bigl(1-\tfrac1m\bigr)\!\sum_p\log R_p\Bigr)\;\le\;\mathrm{BC}.$$
(Cross-multiplication; verified symbolically, difference $\equiv0$.)  So the
$p$-adic input enters as a single additive gain
$$\gamma_p=\Bigl(1-\frac1m\Bigr)\log R_p$$
subtracted from $\tau(\mathbf b)$ — which is exactly the **uniform-slope case**
of the theorem of `consolidation/ADELIC_HOLONOMY.md` §2
($\gamma_p=\varsigma\log p\,(1-1/m)$ with $R_p=p^{\varsigma}$, all $m$ functions
sharing the slope).  Our $m$ functions do share it: $1$ trivially, and
$\delta^jH$ for $j=0,1,2$ because $v_5$ of the $x^n$-coefficient of $\delta^jH$
is $v_5(n^jh_n)\ge v_5(h_n)\ge3n-O(\log n)$ (§11.1).

**Why the weight is exactly $1-\frac1m$.**  It is forced by the product formula:
rescaling the coordinate $x\rightsquigarrow x/p^{v}$ sends
$\log|\varphi'(0)|\mapsto\log|\varphi'(0)|-v\log p$, $\mathrm{BC}\mapsto\mathrm{BC}-v\log p$,
$\log R_p\mapsto\log R_p+v\log p$, and leaves $\tau(\mathbf b)$ alone; the signed
margin then shifts by $v\log p\,(mw-m+1)$, which vanishes **iff** $w=1-\frac1m$.
Any other weight would make the bound depend on an arbitrary choice of
coordinate scale.

#### (b) Critical reading of `ADELIC_HOLONOMY.md` §2.3 — a gap, recorded there

The proof there is "CDT's archimedean proof verbatim + one modified Liouville
step": $\mathrm{den}(\beta)\le e^{\alpha\tau}$ is CDT's, and one adds
$v_p(\beta)\ge\gamma_p\alpha/\log p$.  The hypothesis is exactly what we need —
$(\ast_p)$: $v_p(c_{i,n})\ge\varsigma_i n-o(n)$ on the *actual* coefficients,
with $\mathbb Z$ auxiliary polynomials — and §2.3 correctly observes that the
$\mathrm{lcm}$ denominators are invisible at a fixed $p$
($v_p[1,\dots,bn]=O(\log n)=o(n)$), so the $\tau$-rearrangement and the $p$-adic
gain **do not interact**: they are additive.  That part is sound.

What is **not** sound is the weight.  Part (a) asserts $|\mathbf k|\le\alpha/m$,
whereas the stated facts ($\mathbf k\in[0,D]^d$, $\alpha=mdD/2$) give only
$|\mathbf k|\le dD=2\alpha/m$.  The two estimates actually available are
$$\text{term-by-term: } w=\Bigl(1-\frac1m\Bigr)^{2},\qquad
p\text{-adic Gauss norm: } w=1-\frac2m,$$
the former being the better, and §2.4 note 2 itself observes that
$(1-\frac1m)^2$ *breaks* scale covariance.  A **dated note (2026-08-22)** has been
added to `ADELIC_HOLONOMY.md` §2 recording this, downgrading its status to
"[proved for $w=(1-1/m)^2$]".

This matters here quantitatively: the contradiction needs
$$w\;>\;\frac{\mathrm{RE}/m+\tau(\mathbf b)-\log|\varphi'(0)|}{\log R_5}
=\frac{0.266741+2.8125+0.531289}{4.828314}=\mathbf{0.747783}$$
(and $0.747818$ if the certified $\mathrm{RE}+\varepsilon$ is used),
against $1-\frac1m=0.75$: a slack of only $0.00218$ in $w$, i.e. $0.29\%$.  With
$w=(1-\frac1m)^2=0.5625$ the denominator is $-0.6279<0$ (margin $-3.58$) and with
$w=1-\frac2m=0.5$ it is $-0.9296$ (margin $-4.79$): the method does not even
start.  **The $\zeta_5(3)$ result depends on the exact product-formula weight to
within $0.3\%$**, so it cannot be rescued by any partial version of the adelic
step.

#### (c) The exact dependency chain

| # | input | status |
|---|---|---|
| 1 | **Archimedean holonomy bound over $\mathbb Q$** with the array-refined $\tau(\mathbf b)$: BC numerator = CDT L2chi Thm 2.5.1 (`basic main`) / Thm 7.0.1 (`main:BC form`); rearrangement numerator = Thm 6.0.2 (`main:elementary form`).  *(We checked: both numerator variants are in the published theorems; `basic main` and `main:BC form` carry $\iint\log\lvert \varphi(z)-\varphi(w)\rvert $, `main:elementary form` carries the piecewise rearrangement integral.)* | **published**, arXiv:2408.15403 |
| 2 | **Adelic architecture** — the same bound with radii $R_v$ at all places, equivalently the weight $1-\frac1m$ on $\sum_p\log R_p$ — for the **crude (sup-norm)** numerator | **published**, Dimitrov arXiv:1912.12545 (stated with $p$-adic places in his LNT notes) |
| 3 | **1 + 2 combined**: adelic architecture *with* the refined numerator and $\tau(\mathbf b)$ | **used but not restated or proved** in CDT ICM arXiv:2510.04156 §6.2 eq. (6.2); attributed to the companion **[CDT24]**, *Arithmetic holonomy bounds and the irrationality of the $2$-adic $\zeta(5)$* — not in this repository |
| 3′ | the Liouville-modification route to 3 (`ADELIC_HOLONOMY.md` §2) | gives only $w=(1-\frac1m)^2$, **insufficient here** (§11.9(b)) |
| 4 | **Kontsevich–Zagier**: $H$ satisfies an inhomogeneous ODE over $\mathbb Q(x)$ of exact minimal order $2k+1=3$ | **published**, KZ *Periods* §2.3 Fact 1; corroborated by our rank certificate to $D=120$ (§11.3) |
| 5 | **Buzzard**: $R_5=5^{3}$, i.e. the finite-slope overconvergent $E^*_{-2}$ continues across the whole supersingular annulus of $X_0(5)$ | **published**, JAMS 16 (2003) Thm 5.2, applied as in Calegari 2005 §3; *measured* to be exactly $5^3$ (§11.1) |
| 6 | denominator type $\tau(\mathbf b)=45/16$; $\mathbb Q(x)$-independence; the certified archimedean numerics | **this document**, §§11.2, 11.3, 11.5 |

$$\boxed{\ \text{The single unavailable input is row 3.}\ }$$

#### (d) The theorem, as it would appear

> **Theorem.**  Let $\zeta_5(3)\in\mathbb Q_5$ be the Kubota–Leopoldt $5$-adic zeta
> value, i.e. the unique element of $\mathbb Q_5$ for which
> $$E^*_{-2}=\frac{\zeta_5(3)}{2}+\sum_{n\ge1}\Bigl(\sum_{d\mid n,\,5\nmid d}d^{-3}\Bigr)q^n$$
> is an overconvergent $5$-adic modular form of weight $-2$ and level $\Gamma_0(5)$.
> Then $\zeta_5(3)\notin\mathbb Q$.  Moreover, for all $p/q\in\mathbb Q$ with
> $\max(|p|,|q|)$ sufficiently large,
> $$\Bigl|\zeta_5(3)-\frac pq\Bigr|_5\;>\;\frac1{\max(|p|,|q|)^{688}} .$$
> *(With the Bost–Charles numerator in place of the rearrangement numerator the
> exponent improves to $390$.)*
>
> **Conditional on** the refined multi-place holonomy bound of [CDT24]
> (row 3 above), in the form used by Calegari–Dimitrov–Tang, ICM survey
> arXiv:2510.04156, §6.2 eq. (6.2).
>
> *(Superseded by §11.9(e), which proves that bound for the Bost–Charles numerator
> from published ingredients plus one elementary lemma; the statement then becomes
> unconditional, with exponent $390$ in place of $688$.)*

**Proof sketch.**  Let $x=(\Delta(5\tau)/\Delta(\tau))^{1/4}=q\prod((1-q^{5n})/(1-q^n))^6$
be the hauptmodul of $X_0(5)$, $E^*_2=E_2(\tau)-5E_2(5\tau)$, and
$H=E^*_2\cdot E^*_{-2}=\sum_nh_nx^n$.  Suppose $\zeta_5(3)\in\mathbb Q$.

1. *Arithmetic.*  $h_n\in\mathbb Q$ with $[1,\dots,n]^{3}h_n\in\mathbb Z$ up to a
   fixed denominator (§11.2, proved: the only source of denominators is $d^{-3}$
   in the third formal primitive, and $q\leftrightarrow x$ is integral both ways).
   Hence the array $\mathbf b$ with $r=3$, $b_{i,j}=1$ for $i\ge2$, $u_j=1$, and
   $\tau(\mathbf b)=3(1-\tfrac1{16})=\tfrac{45}{16}$.
2. *$5$-adic.*  $E^*_{-2}$ is overconvergent of slope $0$, so by Buzzard it
   continues over the entire supersingular annulus of $X_0(5)$; with the Fricke
   description of the ordinary locus this is the disc $|x|_5<5^{3}$, so
   $R_5=5^3$ and $R_\ell=1$ for $\ell\ne5$.
3. *The inventory.*  By Kontsevich–Zagier the minimal inhomogeneous ODE for $H$
   over $\mathbb Q(x)$ has order exactly $3$, hence $1,H,\delta H,\delta^2H$
   ($\delta=x\,d/dx$) are $m=4$ $\mathbb Q(x)$-linearly independent elements of
   $\mathcal V(\varphi_5,\varphi_\infty;\mathbf b)$, all of the same denominator
   type and all with meromorphic pullback.
4. *The template.*  $H(x(q))$ and its $\delta$-derivatives are meromorphic on the
   whole disc $|q|<1$, so $\varphi=x\circ\psi$ is admissible for every
   holomorphic $\psi:\mathbb D\to\mathbb D$ with $\psi(0)=0$.  Take
   $\psi(z)=z\exp\bigl(\sum_{k=0}^{40}c_kz^k\bigr)$ with the $c_k$ of §11.5;
   then $\max_{|z|=1}|\psi|\le e^{-0.1788}<1$, $\log|\varphi'(0)|=c_0$, and the
   contour lies in $|q|\le0.8363$.
5. *The bound.*  Apply the multi-place bound with $\varphi_5(z)=5^{-3}z$:
   $$4=m\;\le\;\frac{\mathrm{RE}(\varphi)+3\log5}{c_0+3\log5-\tfrac{45}{16}}\;\le\;3.9716148\;<\;4 ,$$
   a contradiction.  The measure follows from the same inequality with CDT's
   $\kappa$-refinement, $\gamma=1/m=1/4$, $\rho_5=5^{-3}$. $\square$

**The numerical certificate.**

| quantity | value | provenance |
|---|---|---|
| $m$ | $4$ | KZ minimal order $3$ + §11.3 Lemma; rank full to $D=120$ |
| $\tau(\mathbf b)$ | $45/16=2.8125$ | §11.2, proved and sharp |
| $\log R_5$ | $3\log5=4.828313737$ | Buzzard Thm 5.2; measured slope $2.99063$ at $n=1600$, deficit $O(\log n)$ |
| $c_0=\log\lvert\varphi'(0)\rvert$ | $-0.531289158096$ | exact ($x'(0)=1$, $\psi'(0)=e^{c_0}$) |
| $\max_t\log\lvert\psi\rvert$ | $\le-0.1788$ | $2^{22}$ grid + exact trig-polynomial Lipschitz constant |
| $\max\lvert q\rvert$ on the contour | $0.836265271$ | — |
| $\operatorname{Lip}(g)$ | $\le11427.0$ | a-priori: $2\pi M(r)(1+\sum k\lvert c_k\rvert)$ |
| $\mathrm{RE}$ | $1.066964996\pm6.81\cdot10^{-4}$ | $N=2^{24}$, error $=\operatorname{Lip}/N$ |
| $\mathrm{BC}$ | $1.03552928$ | $N$-converged to $10^{-8}$ |
| **bound** | $\le3.9716148$ (RE) / $3.9499804$ (BC) | vs $m=4$ |
| **margin** | $\ge+0.0421385$ (RE) / $+0.0742553$ (BC) | |
| $\kappa$ | $\le687.2$ (RE) / $389.8$ (BC) | §7 formula, $\gamma=1/4$ |
| crude (sup-norm) form | $m\le5.714335$ — **fails** | §11.8; needs the refined form |
| required adelic weight | $w>0.747818$; product formula gives $w=0.75$ (slack $0.0022$, i.e. $0.29\%$) | §11.9(b) |

#### (e) Update (2026-08-22): the citation gap is closed for the Bost–Charles numerator

Row 3 of the chain above — the adelic form with the refined numerator — has now
been **proved** for the Bost–Charles numerator, by the single-variable slopes
method.  The argument is written out in `consolidation/ADELIC_HOLONOMY.md` §2.6
("Proof of the BC-numerator form with weight $1-\frac1m$ via the slopes method");
here is what it does and why it works where the Thue–Siegel route (§2.3 there,
§11.9(b) here) does not.

**The mechanism.**  Run CDT's own proof of L2chi Thm `main:BC form` (§6.3)
with one change: at each $p\in S$ replace the $\mathbb Z$-structures by their
$R_p$-Gauss counterparts,
$$\Lambda_D:=\bigoplus_{j=0}^{D}\mathbb Z\,(R_px)^{-j}\ \supset\ \mathbb Z[1/x]_{\le D},
\qquad G^{(n)}:=R_p^{\,n-D}\,\mathbb Z\cdot x^{n-D},$$
leaving every **archimedean** norm untouched.  Then, since $h_v$ depends only on
the norms at $v$ and the real vector spaces and their Euclidean norms are
unchanged, *all* of CDT's archimedean estimates hold verbatim, and three
elementary computations account for the twist:

| ingredient | contribution to the slopes inequality | source |
|---|---|---|
| $\widehat{\deg}\,\overline E_D$ gains the index $\prod_pR_p^{D(D+1)/2}$ per summand | $+\frac m2D^2\sum_p\log R_p$ | elementary |
| $\widehat\mu_{\max}(\overline{G^{(n)}})=-(n-D)\log R_p$, summed with $\sum_{n\in\mathcal V_D}n=\frac{m^2}2D^2+o(D^2)$ | $-\bigl(\frac{m^2}2-m\bigr)D^2\sum_p\log R_p$ | L2chi Thm `KolchinSolved` + $\#\mathcal V_D=m(D+1)$ — the *same* input CDT already use for the archimedean sum |
| **Lemma B** (the $p$-adic operator norm): $h_p(\psi_D^{(n)})\le\log(Cn^A)$ | $o(D^2)$ | new, one paragraph |

Net: $-\bigl(\frac{m^2}2-m\bigr)+\frac m2=-\frac{m(m-1)}2$, i.e. **weight
$1-\frac1m$**.  Assembling with L2chi (BC-sum-arch-ht), (BC-sum-finite-ht),
(BCardegED) and optimising over $\xi,y_h$ exactly as CDT do gives
$$m\bigl(\log|\varphi'(0)|-\tau(\mathbf b;\mathbf e)\bigr)+(m-1)\sum_p\log R_p\ \le\ \mathrm{BC}(\varphi),$$
verified symbolically (the reduction's two sides differ by $0$ identically).

**Lemma B.**  With $\|\sum_{i,j}q_{ij}f_ix^{-j}\|_{E,p}=\max_{i,j}|q_{ij}|_pR_p^{-j}$
and $\|c\,x^{n-D}\|_{G,p}=|c|_pR_p^{\,n-D}$, and with
$|c_{i,k}|_p\le Cn^{A}R_p^{-k}$ on the *actual* coefficients, the ultrametric
inequality gives
$$|c_n(Q)|_p=\Bigl|\sum_{i,j}q_{ij}c_{i,n-D+j}\Bigr|_p\le Cn^{A}\max_{i,j}|q_{ij}|_pR_p^{-(n-D+j)}=Cn^{A}R_p^{-(n-D)}\|Q\|_{E,p},$$
so $\|\psi_D^{(n)}[Q]\|_{G,p}\le Cn^{A}\|Q\|_{E,p}$; and since $c_n(Q)$ depends only
on the class $[Q]$, the bound passes to the quotient norm.  $\square$

**Why the single-variable route succeeds.**  Both routes pay "$-D$ per slot" for
the pole of the auxiliary polynomials at $x=0$ — in §11.9(b) it appears as
$\sum_s(n_s-k_s)\ge\alpha(1-\frac1m)^2$, here as
$\widehat\mu_{\max}(\overline{G^{(n)}})=-(n-D)\log R_p$.  Only the slopes method
**pays it back**, through the $+\frac m2D^2\log R_p$ that the twisted lattice adds
to $\widehat{\deg}\,\overline E_D$; the Thue–Siegel route never forms the
arithmetic degree of the auxiliary lattice and so cannot.

**Sharpness check.**  For $m=2$, $\mathbf b=0$ and univalent $\varphi(z)=\rho z$
(so $\mathrm{BC}=\log\rho$), the theorem reads: $\rho\prod_pR_p>1$ forces
$f\in\mathbb Q(x)$ — the **Borel–Dwork rationality criterion**, exactly.  The
weight $(1-\frac1m)^2$ would give only $\rho\prod_pR_p^{1/2}>1$.  So $1-\frac1m$
is sharp and the Thue–Siegel weight is genuinely lossy.

**Consequence: use the BC numerator.**  Our two numerators give

| numerator | multi-place status | bound | margin | $\kappa$ |
|---|---|---|---|---|
| Bost–Charles | **proved** with weight $1-\frac1m$ (`ADELIC_HOLONOMY.md` §2.6) | $3.9499804$ | $+0.0742553$ | $\le389.8$ |
| rearrangement | only $(1-\frac1m)^2$ — **insufficient** ($w$ needed $>0.7478$) | $3.9716148$ | $+0.0421385$ | $\le687.2$ |

so the $\zeta_5(3)$ verdict should be read off the **Bost–Charles** row, whose
margin is the larger of the two anyway.  The revised dependency chain is:

| # | input | status |
|---|---|---|
| 1 | Bost's slopes inequality (adelic by construction, $h(\psi)=\sum_vh_v$) | **published**, Bost, *Foliations*, Prop. 4.6 |
| 2 | Bost–Charles: $(\overline{\mathcal O(1)}_\varphi\cdot\overline{\mathcal O(1)}_\varphi)=\iint\log\lvert \varphi(z)-\varphi(w)\rvert $, and arithmetic Hilbert–Samuel — used at $\infty$ only | **published**, Bost–Charles Thm 5.4.1/Prop. 5.4.2; Zhang, Gillet–Soulé, Bost |
| 3 | CDT's archimedean and finite-place evaluation-height estimates and $\tau(\mathbf b;\mathbf e)$ (L2chi (BCardegED), (BC-sum-arch-ht), (BC-sum-finite-ht), Thm `KolchinSolved`) | **published**, arXiv:2408.15403 §6.3 |
| 4 | the $R_p$-Gauss twist + Lemma B | **ours**, `ADELIC_HOLONOMY.md` §2.6 — one paragraph, elementary |
| 5 | Kontsevich–Zagier: minimal inhomogeneous order $2k+1=3$ | **published**, KZ *Periods* §2.3; rank certificate to $D=120$ here |
| 6 | Buzzard Thm 5.2 $\Rightarrow R_5=5^3$ | **published**, JAMS 16 (2003); measured exactly here (§11.1) |
| 7 | $\tau(\mathbf b)=45/16$, $\mathbb Q(x)$-independence, the certified numerics | **this document**, §§11.2, 11.3, 11.5 |

**Nothing in the chain is now unavailable.**  Row 4 is new — written today, not
refereed — and it is the only step that is not a citation.  Its hypotheses are met
here: $R_5=5^3$ is an integral power of $5$ (needed only so that $R_5^{-j}\in\mathbb Q$);
$(\ast_5)$ holds with a polynomial factor, since $v_5$ of the $x^n$-coefficient of
$\delta^jH$ is $v_5(n^jh_n)\ge v_5(h_n)\ge3n-O(\log n)$ by §11.1; and the
$5$-part of the lcm denominators is $O(\log n)$, so it is absorbed in the
$o(D^2)$.

**Restated theorem.**  With row 4 accepted, the statement of §11.9(d) becomes
**unconditional**, with the Bost–Charles exponent:

> **Theorem.**  $\zeta_5(3)\notin\mathbb Q$; and for all $p/q\in\mathbb Q$ with
> $\max(|p|,|q|)$ sufficiently large,
> $\bigl|\zeta_5(3)-\tfrac pq\bigr|_5>\max(|p|,|q|)^{-390}$.

**The remaining gap, stated exactly.**  There is no longer a missing *published*
statement; what remains is that step 4 is our own and unreviewed.  Concretely, the
three assertions a referee must check are: (i) that replacing the $\mathbb Z$-structure
of $\Gamma(\mathcal X,\mathcal L^{\otimes D})$ by $\Lambda_D$ leaves every
archimedean estimate of L2chi §6.3 untouched — it does, because $E_{D,\mathbb R}$
and its Euclidean norm are unchanged and $h_v$ depends only on the $v$-adic norms;
(ii) the index computation $[\Lambda_D:\mathbb Z[1/x]_{\le D}]=R_p^{D(D+1)/2}$ and
$\widehat\mu_{\max}(\overline{G^{(n)}})=-(n-D)\log R_p$; (iii) Lemma B.  We regard
(i)–(iii) as routine, but they have not been checked by anyone else, and the
$\zeta_5(3)$ margin ($1.5\%$ of budget in the BC row) leaves no room for an error
in the *weight*.

### 11.10 The certified statement, and the caveats that remain open

**Certified** (float64 arithmetic; every error term a-priori):

| input | value | how |
|---|---|---|
| $m$ | $4$ | Kontsevich–Zagier minimal order $2k+1=3$ + §11.3 Lemma; rank certificate to $D=120$ |
| $\tau(\mathbf b)$ | $45/16=2.8125$ exactly | §11.2 Proposition (proof) + sharpness by exact computation |
| $\log R_5$ | $3\log5=4.828313737$ | Buzzard Thm 5.2 (cited); *measured* to be exactly $3$ at $n\le1600$ |
| $\log\lvert\varphi'(0)\rvert$ | $c_0=-0.531289158096$ exactly | $x'(0)=1$, $\psi'(0)=e^{c_0}$ |
| admissibility of $\psi$ | $\max_t\log\lvert\psi\rvert\le-0.1788$ | $2^{22}$ grid + exact Lipschitz constant of a trig. polynomial |
| $\mathrm{RE}$ | $1.066964996\pm6.9\cdot10^{-4}$ | $N=2^{24}$ + a-priori $\operatorname{Lip}(g)\le11427$ |
| crude (sup) form | best $m\le5.7143$ — **fails** | §11.8; the claim needs the refined form |
| $\mathrm{BC}$ | $1.03552928$ ($N$-converged to $10^{-8}$) | the numerator whose multi-place form is proved (§11.9(e)) |

$$\boxed{\ \text{bound}\;\le\;3.9716148\;<\;4\;=\;m,\qquad \text{margin}\;\ge\;+0.0421385,\qquad \kappa\;\le\;687.2\ }$$

with the **rearrangement** numerator (a-priori error bound, but its multi-place
form is only proved with weight $(1-\frac1m)^2$), and

$$\boxed{\ \text{bound}\;=\;3.9499804\;<\;4,\qquad \text{margin}\;=\;+0.0742553,\qquad \kappa\;\le\;389.8\ }$$

with the **Bost–Charles** numerator, whose multi-place form is proved in
§11.9(e) — this is the operative one.

**Caveats that remain open**, in order of importance:

1. **The refined multi-place bound — no longer a missing citation, but step 4 of
   the chain is ours.**  §11.9(e) proves the adelic form *with the Bost–Charles
   numerator* and the sharp weight $1-\frac1m$, by an $R_p$-Gauss twist of the
   evaluation lattice in CDT's own §6.3 slopes proof.  Everything else in the
   chain is published (Bost's slopes inequality; Bost–Charles's self-intersection
   formula and arithmetic Hilbert–Samuel; CDT's evaluation-height estimates and
   $\tau(\mathbf b)$).  What remains is that the twist and Lemma B are **new and
   unreviewed** — see §11.9(e) for the three assertions a referee must check.  The
   *published crude* form (Dimitrov arXiv:1912.12545) fails here ($m\le5.7143$,
   §11.8), and the rearrangement-numerator route gives only $(1-\frac1m)^2$, which
   is insufficient — so the verdict is read off the **Bost–Charles** row.

2. **The rearrangement-numerator multi-place form remains unproved** with the sharp
   weight (only $(1-\frac1m)^2$, insufficient).  The headline therefore uses the
   Bost–Charles numerator, whose margin is the larger ($+0.0743$) but whose
   quadrature is the softer: $O(1/N)$ rather than $O(1/N^2)$ where $\varphi$ is
   non-injective on the boundary.  On this contour $\mathrm{BC}$ is $N$-stable to
   $10^{-8}$ over $N=8192\dots32768$, so this is not a practical worry; but unlike
   $\mathrm{RE}$ it does not yet carry an a-priori error bound.
3. **Kontsevich–Zagier minimal order $2k+1$** is cited, not proved here.  §11.3
   reduces the whole independence hypothesis to it and gives the structural
   reason; the rank certificate excludes any relation of degree $\le120$ while the
   genuine order-$3$ relation has degree $8$.
4. **Buzzard's theorem at weight $-2k$ and $p=5$** — Calegari's §3.2 asserts the
   $p=5,7,13$ cases without repeating the argument; we follow him.  The
   coefficient measurement of §11.1 shows the radius is exactly $5^3$, so nothing
   is being over-claimed, but the *proof* that it is at least $5^3$ is his.
5. **No interval arithmetic.**  The float64 chain has $\sim10^{-11}$ round-off
   against an error budget of $6.8\times10^{-4}$; a certified enclosure of
   $\log|\Delta|$, the $\mathrm{SL}_2(\mathbb Z)$ reduction, the FFT and the sort
   has not been implemented.
6. **Novelty is unaudited** (§8): the check against Lai–Sprang–Zudilin 2025 has
   not been made.

**Not a caveat:** the contour.  It is a lower bound by construction, and a plain
disc already gives a positive margin.
