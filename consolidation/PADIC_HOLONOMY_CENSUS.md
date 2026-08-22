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
the rearrangement integral and the lune contours of Appendix A.*

**Nothing here is a claimed theorem.**  §9 is the honest computed/estimated/assumed
ledger; §8 flags the one value that would be new *if* the numerics and the
quoted multi-place form of the bound both hold.

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
* **$\boxed{\zeta_5(3)}$ — the primary target — gives a contradiction.**  With
  $m=4$, $R_5=5^3$, $\tau(\mathbf b)=\frac{45}{16}$ and an explicit template
  $\psi(z)=z\exp\bigl(\sum_{k\le16}c_kz^k\bigr)$ lying in $|q|\le0.786$:
  $$m\;\le\;3.9727173\ \ (\text{rearrangement})\qquad\text{resp.}\qquad
  3.9504220\ \ (\text{Bost--Charles})\qquad <\;4=m .$$
  Margin $+0.0404$ out of a budget of $3.2349$ ($1.3\%$); irrationality measure
  $\kappa\approx716$ (RE) resp. $394$ (BC).  Calegari's elementary criterion
  fails here ($\theta_5=0.8918$) and Beukers 2008 reaches only
  $\zeta_5(3)\pm L_5(3,\chi_5)$, so **this would be a new theorem** — subject to
  the ledger in §9, in particular the multi-place form of the bound (A1/A2) and
  the novelty audit of §8.
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
  | bound | 2.276 | 2.691 | 2.028 | 4.439 | 1.871 | **3.9727** | 3.343 | 6.822 | 8.063 | 15.595 | — |
  | margin | $+7.394$ | $+3.795$ | $+2.879$ | $+4.188$ | $+4.127$ | $\mathbf{+0.040}$ | $-0.390$ | $-2.048$ | $-2.367$ | $-6.462$ | n/a |

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
cusp $\infty$ (where $x=0$), i.e. $q=\zeta_p^{\,j}$ at depth $c=p$.  For a disc
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
| $\zeta_5(3)$ | $X_0(5)$, $k{=}1$ | 4 | $5^{3}$ | 4.82831 | 2.81250 | 3.23494 | -0.533513 | 1.060448 | 1.027399 | **3.972717** | 3.950422 | +0.04044 | **contradiction** |
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
Extending to degree $16$ (`cert_05.json`) gives the best found,
$$\boxed{\ \mathrm{RE}=1.0604479,\quad m\le 3.9727173<4=m,\quad\text{margin}=+0.0404412\ }$$
and $\mathrm{BC}=1.0273996$, bound $3.9504220$, margin $+0.0734895$.  Independent
re-optimisations at $K=24$, $N=2048$ reproduce $\text{margin}\in[+0.0398,+0.0414]$,
so the optimum of this family is $\approx+0.041$ and the search has converged.

Irrationality measure (§7): $\kappa=716$ from RE, $\kappa=394$ from BC.  These are
large because the margin is small; the qualitative conclusion is the point.

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
| $\zeta_5(3)$ | $+0.0404$ | $716.0$ | $+0.0735$ | $393.8$ | **no previous value** |

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
| $\boxed{\zeta_5(3)}$ | **contradiction**, margin $+0.041$ (RE) / $+0.073$ (BC) | Calegari fails; Beukers reaches only $\zeta_5(3)\pm L_5(3,\chi_5)$; **would be NEW** — flag for audit against LSZ 2025 |
| $\zeta_2(7)$, $\zeta_3(5)$, $\zeta_7(3)$ | **no contradiction** (margins $-6.46$, $-2.35$, $-2.04$) | open, and stays open by this route |
| $\zeta_5(5)$ | **method inapplicable**: the entry condition fails for *every* template | open |

**The one candidate new theorem is $\zeta_5(3)\notin\mathbb Q$**, i.e. the $5$-adic
zeta value at $3$, equivalently $\zeta_5(3)=L_5(3,\mathrm{id})$ in Calegari's
notation.  Its margin is small ($+0.041$ out of a budget of $3.235$, i.e. $1.3\%$),
so §9 must be read before any claim is made.

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
* **Independent end-to-end recomputation of the $\zeta_5(3)$ contour.**  With the
  contour lying in $|q|\le0.7852$, the direct $\eta$-product converges, so
  $\log|x|$ was recomputed at $40$ digits in `mpmath` with **no** modular
  reduction of any kind: agreement $4.9\times10^{-15}$ at every sampled point, and
  the recomputed bound is $3.972717235$ with margin $+0.040441262$.  The identity
  $\log|\psi'(0)|=c_0=\overline u$ holds to machine precision.
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

**A1. The multi-place form of the bound.**  The inequality actually used,
$$m\le\frac{\mathrm{BC}+\sum_v\log R_v}{\log|\varphi'(0)|+\sum_v\log R_v-\tau(\mathbf b)},$$
is quoted from CDT's ICM survey §6.2 eq. (6.2) (their $\zeta_2(5)$ proof), where
it is attributed to [L2chi, Thm 2.5.1] with the $p$-adic template
$\varphi_p(z)=R_p^{-1}z$.  The copy of L2chi in `papers/cdt/cdt2/L2chi.tex`
states Theorems 2.5.1/6.0.2/7.0.1 over $\mathbb Q$ at the archimedean place only;
we did **not** locate the multi-place statement (their §15.7) in that file.  We
take the ICM form at face value — it is calibrated exactly against their own
$\zeta_2(5)$ numbers in §3, which is strong evidence the transcription is right,
but it is an input, not a verification.

**A2. The rearrangement form of the multi-place bound.**  Nazarov's
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
