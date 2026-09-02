# Split denominators in two variables: the companion cannot be split

*Working note, 2026-09-02. All computations exact (Python `Fraction`; mod-$p$ linear algebra
for recurrence fitting); scripts `lib2v.py`, `libx.py`, `01_pullback.py` ... `22_firstint.py` and outputs `out*.txt` in this
directory. Framework: `consolidation/TWO_VARIABLE_HOLONOMY.md` (naive bound
$m\le 2M_1M_2/(e_1e_2)$, multiplicativity of tensor inventories, and Section 4 "what a gainful
system would need"), `CDT_UNPACKED.md` Section 1, `INVENTORY_BOUND.md`,
`COMPANION_ARITHMETIC.md` Section 1.*

## 0. Summary

The question: does a **non-tensor** two-variable extension of an Apery-like row admit a
companion with **split denominators** $[1..a]^{r_1}[1..b]^{r_2}$, $r_1+r_2=r$?  The answer is
**no**, and the obstruction is structural, not numerical.

1. **The natural decouplings are pullbacks, not tensor products.**  Every separated decoupling
   $c_{a,b}=\sum_k f(a,k)g(b,k)$ whose $k$-generating functions are geometric is a *pullback*
   $C(x,y)=R(x,y)h(u(x,y))$ with $u=\alpha(x)\beta(y)$, of a **lower-weight** one-variable
   system.  For Apery's $\zeta(3)$ row,
   $$c_{a,b}=\sum_k\binom ak\binom bk\binom{a+k}k\binom{b+k}k
     =\sum_k\binom{2k}k^2\binom{a+k}{2k}\binom{b+k}{2k},\qquad
     C=\frac{{}_2F_1(\tfrac12,\tfrac12;1;16u)}{(1-x)(1-y)},\quad
     u=\frac{xy}{(1-x)^2(1-y)^2}.$$
   Rank $2$ (not $4$); singular locus the **irreducible** curve $16xy=(1-x)^2(1-y)^2$ plus the
   polar lines $x=1$, $y=1$ -- genuinely coupled, not a union of horizontal and vertical lines,
   not a tensor product.  Weight has dropped from $3$ to $1$: the $\zeta(3)$-ness survives only
   in the Hadamard diagonal.

2. **A two-variable companion with the period constant along the whole fold curve exists** --
   the object `TWO_VARIABLE_HOLONOMY.md` Section 4 asks for -- **but its denominators are of
   MAX type.**  The power-series solution space of the two-variable system is exactly
   $2$-dimensional, spanned by the row $c$ and one companion $s$ with
   $s_{a,b}/c_{a,b}\to\zeta(3)-1$ along *every* ray $a=\lambda b$, at the same exponential rate
   as Apery's own forms; and
   $$\operatorname{den}(s_{a,b})=\operatorname{lcm}(1..\max(a,b))^{3}$$
   exactly (up to bounded factors), so the only uniform split is $(r_1,r_2)=(3,3)$: **the rates
   sum to $2r$, not $r$.**  Since the solution space is $2$-dimensional and the row is integral,
   *every* companion of the system has this denominator type: it is an invariant of the system,
   not a choice of splitting.

3. **The split (MIN-type) candidate exists but diverges.**  The unique pullback series with
   diagonal equal to the one-variable companion,
   $d_{a,b}=\sum_k g_k\binom{a+k}{2k}\binom{b+k}{2k}$ with $d_{n,n}=b_n$, has the *optimal*
   denominator $\operatorname{lcm}(1..\min(a,b))^{3}$ -- split for every $(r_1,r_2)$ with
   $r_1+r_2=3$ -- but $g_k$ is not holonomic and grows like $(k!)^{\gamma}$: $h(u)=\sum g_ku^k$
   is a divergent Gevrey series, $d_{a,b}/c_{a,b}$ blows up off the diagonal, and it satisfies
   neither recurrence.

A sweep over ten decouplings of four rows ($\zeta(3)$, $\zeta(2)$, Zagier's Catalan row E,
Cooper's $s_{10}$) confirms this uniformly: in every case, every solution whose ray limit has a
nonzero period coefficient is MAX-type of rate $r$ or worse; the only solutions with split (or
better) denominators are "small" harmonic solutions whose ratio to the row tends to $0$.

Together: **denominator-splitting and period-constancy are exactly complementary.**  The
$(3,0)$-solution (no $b$-denominators at all) has fixed-$b$ Apery limit $\zeta(3)-H^{(3)}_b$;
to make the period constant one must add $H^{(3)}_b c_{a,b}$, and that costs exactly
$\operatorname{lcm}(1..b)^3$ -- one full lcm layer, precisely what the split was to save.  Call
this the **tail obstruction**.

On the conformal side there is a matching conservation law.  For a pullback host along a
*product* $u=\alpha(x)\beta(y)$, **Landau's theorem** gives the sharp ceiling
$$\ell_1+\ell_2\ \le\ \log 16 = 2.7726,\qquad \ell_k=\log|\varphi_k'(0)|,$$
whereas in one variable the ceiling is $\log(16|t_{\rm outer}|)$ with $t_{\rm outer}$ the outer
singularity (for Apery $t_2=(1+\sqrt2)^4=33.97$, so $\ell\le6.298$).  Entry in two variables
with max-type denominators needs $\ell_1+\ell_2>2r=6$: it fails by $3.23$ nats, where the
one-variable test passes by $3.30$.

## 1. The two-variable rows (Task 1)

### 1.1 Separated decouplings and the pullback lemma

Write $c_{a,b}=\sum_kf(a,k)g(b,k)$, so $C(x,y)=\sum_kF_k(x)G_k(y)$.

**Lemma 1.1 (pullback).** If $F_k(x)=F(x)\alpha(x)^k$ and $G_k(y)=G(y)\beta(y)^k$ then
$C(x,y)=F(x)G(y)h(u)$ with $u=\alpha(x)\beta(y)$: the two-variable system is the pullback of
the one-variable system of $h$ along the dominant map $u$, twisted by a rational prefactor --
same rank, hence *sub*-multiplicative.

The classical Apery-type binomial sums satisfy the hypothesis after
$\binom ak\binom{a+k}k=\binom{a+k}{2k}\binom{2k}k$ and
$\sum_a\binom{a+k}{2k}x^a=x^k(1-x)^{-2k-1}$, $\sum_a\binom ak x^a=x^k(1-x)^{-k-1}$,
$\sum_b\binom{2b-2k}{b-k}y^b=y^k(1-4y)^{-1/2}$.  Verified coefficientwise for $0\le a,b\le25$
(`01_pullback.py`):

| row | decoupling $c_{a,b}$ | closed form | $u$ | rank |
|---|---|---|---|---|
| $\zeta(3)$ | $\sum_k\binom ak\binom bk\binom{a+k}k\binom{b+k}k$ | $\frac{{}_2F_1(\frac12,\frac12;1;16u)}{(1-x)(1-y)}$ | $\frac{xy}{(1-x)^2(1-y)^2}$ | 2 |
| $\zeta(2)$ | $\sum_k\binom ak\binom{2k}k\binom{b+k}{2k}$ | $\frac{1}{(1-x)(1-y)\sqrt{1-4v}}$ | $v=\frac{xy}{(1-x)(1-y)^2}$ | 1 (algebraic) |
| row E | $\sum_k\binom ak\binom{2k}k\binom{2b-2k}{b-k}$ | $\frac{1}{(1-x)\sqrt{1-4y}\sqrt{1-4w}}$ | $w=\frac{xy}{1-x}$ | 1 (algebraic) |

Diagonals are the rows ($c_{n,n}=a_n$).  Note the weight drop: $\zeta(3)$'s weight-3 rank-3 MUM
system becomes the weight-1 rank-2 elliptic system $\theta^2=4u(2\theta+1)^2$; the $\zeta(2)$
and Catalan rows become **algebraic**.  The Apery structure lives in the Hadamard diagonal,
not in the two-variable function.

**Tensor test.**  $C$ does not factor ($c_{a,0}=c_{0,b}=1$, $c_{1,1}=5$); the rank is $2<2\cdot2$;
the singular locus is not a union of horizontal and vertical lines.  These are non-tensor,
genuinely coupled systems -- the class `TWO_VARIABLE_HOLONOMY.md` Section 4 identifies as the
only hope, and the class in which the rank is sub-multiplicative.

### 1.2 Recurrences and singular locus (zeta(3) D1)

Exact minimal recurrences over $\mathbf Q$ (`04_exactrec.py`: mod-$p$ nullspace at three primes,
CRT, rational reconstruction), order 2, coefficient degree 3:
$$(a+2)^3c_{a+2,b}-\bigl[2a^3+9a^2+(15+4b+4b^2)a+9+6b+6b^2\bigr]c_{a+1,b}+(a+1)^3c_{a,b}=0,$$
and the mirror in $b$.  With $a=n-1$, $B=b(b+1)$:
$$(n+1)^3u_{n+1}=\bigl[2n^3+3n^2+(3+4B)n+1+2B\bigr]u_n-n^3u_{n-1}.\tag{1.2}$$
A one-parameter deformation of Apery's recurrence with the same $\lambda_1\lambda_2=1$
("Apery-perfect") shape; both characteristic roots are $1$ because $c_{a,b}$ is, for fixed $b$,
a polynomial in $a$ of degree $2b$.  The leading and trailing coefficients depend on $a$
**only** -- this is what forces the denominator dichotomy of Section 2.

Singular locus: $\Gamma=\{16xy=(1-x)^2(1-y)^2\}$ together with $x=1$, $y=1$.  $\Gamma$ is
**irreducible**: as a quadratic in $x$ its discriminant is $64y(1+y)^2$, whose square root
$8(1+y)\sqrt y$ is not rational in $y$; the two points of $\Gamma$ over a fixed $y$ are Galois
conjugate over $\mathbf Q(y)$.  On the diagonal they are $3\pm2\sqrt2=(\sqrt2\pm1)^2$, whose
squares are Apery's $t_{1,2}=(\sqrt2\mp1)^4$: **the one-variable fold and outer singularities
are the two branches of one irreducible two-variable curve**, and the polyradii obey
$\rho_1\rho_2\le(3-2\sqrt2)^2=(\sqrt2-1)^4=t_1$ -- the product of the two-variable radii is the
one-variable radius (Hadamard).

## 2. The two-variable companions (Task 2)

### 2.1 The solution space is two-dimensional; the companion is forced

Parametrising power-series solutions of the pair of recurrences by
$(t_{00},t_{10},t_{01},t_{11})$ and imposing the $a$-recurrence for $b\ge2$ as linear conditions
(`08_solspace.py`) gives a **2-dimensional** consistent space, spanned by

* $c$: data $(1,1,1,5)$ -- the row;
* $s$: data $(-1,0,0,1)$, with $s_{a,b}/c_{a,b}\to\zeta(3)-1$;

equivalently $\tilde d:=c+s$ (data $(0,1,1,6)$) has $\tilde d_{a,b}/c_{a,b}\to\zeta(3)$.  Both
are symmetric in $a\leftrightarrow b$.  Two independent constructions of the same object:

* the "second solution in $a$" $d^{(R)}$ of (1.2) with $d^{(R)}_{0,b}=0$, $d^{(R)}_{1,b}=1$ has
  $d^{(R)}_{a,0}=H^{(3)}_a$ exactly, denominators $\operatorname{lcm}(1..a)^3$ with **no $b$
  denominators at all** (rates $(3,0)$), and fixed-$b$ Apery limit
  $$\xi(b)=\lim_{a\to\infty}d^{(R)}_{a,b}/c_{a,b}=\zeta(3)-H^{(3)}_b=\sum_{m>b}m^{-3}$$
  (verified to 20 digits for $b=0,\dots,8$, `06_secondsol.py`).  It does **not** satisfy the
  $b$-recurrence.
* $\tilde d=d^{(R)}+H^{(3)}_b\,c$ satisfies **both** recurrences exactly (residual $0$ on the
  tested grid) and has constant period $\zeta(3)$.

### 2.2 The tail obstruction

> The $b$-denominator-free solution has period $\zeta(3)-H^{(3)}_b$; the correction that makes
> the period constant is exactly $H^{(3)}_b$, whose denominator is $\operatorname{lcm}(1..b)^3$.

Measured denominators of $s$ on $0\le a,b\le34$ (`09_denom_law.py`):

* $\operatorname{den}(s_{a,b})\mid\operatorname{lcm}(1..\max(a,b))^3$ -- **true**;
  $\mid\operatorname{lcm}(1..\max(a,b))^2$ -- false; the quotient
  $\operatorname{lcm}(1..\max)^3/\operatorname{den}$ takes values $1,5,9,23,125,\dots$: the
  max-type law is **sharp**, not merely an upper bound;
* $\mid\operatorname{lcm}(1..\min(a,b))^3$ -- **false**;
* for $a<b$ no power of $\operatorname{lcm}(1..a)$ alone clears the denominator (primes in
  $(a,b]$ occur), and symmetrically; hence the only uniform split is $(3,3)$.

Because the solution space is 2-dimensional and $c$ is integral, every solution $\alpha c+\beta s$
with $\beta\ne0$ has this type.  **The denominator type of the companion is an invariant of the
two-variable system.**

**Theorem 2.1 (closed form of the companion; proves the max-type bound).**  Verified exactly on
$0\le a,b\le30$ (`20_identity.py`):
$$s_{a,b}\;=\;\underbrace{d^{(R)}_{a,b}}_{\text{den}\,\mid\,\operatorname{lcm}(1..a)^3,\;
\text{period }\zeta(3)-H^{(3)}_b}\;+\;\underbrace{\bigl(H^{(3)}_b-1\bigr)}_{\text{den}\,\mid\,
\operatorname{lcm}(1..b)^3}\,c_{a,b}.$$
Since $c$ is integral and $\operatorname{den}(d^{(R)}_{a,b})\mid\operatorname{lcm}(1..a)^3$
(checked for $a\le26$, all $b$), this **proves**
$\operatorname{den}(s_{a,b})\mid\operatorname{lcm}(1..a)^3\vee\operatorname{lcm}(1..b)^3
=\operatorname{lcm}(1..\max(a,b))^3$.  It is sharp: for a prime $p$ with $a<p\le b$ and
$p^2>b$, the first term is a $p$-unit and $v_p(H^{(3)}_b)=-3$ generically, so
$v_p(\operatorname{den}s_{a,b})=3$ -- confirmed in $1099$ of $1200$ such cells (the $101$
exceptions are the ones where $p\mid c_{a,b}$ or a Wolstenholme-type cancellation occurs).

The identity is the tail obstruction in one line: the split part carries the *variable* period
$\zeta(3)-H^{(3)}_b$, and the rational correction that makes the period constant is the whole of
the $\operatorname{lcm}(1..b)^3$ layer.

*Mechanism.*  The $a$-recursion divides by $(a+2)^3$ at every step, the $b$-recursion by
$(b+2)^3$; reaching $(a,b)$ needs steps in both directions, so primes up to $\max(a,b)$ enter
to the third power.  For the row they all cancel (integrality); for any other solution they do
not.  The saving over the naive product $\operatorname{lcm}(1..a)^3\operatorname{lcm}(1..b)^3$
is real (primes $\le\min$ occur to power 3, not 6) but it is the wrong saving: it does not
split the rate.

### 2.3 The MIN-type candidate: optimal denominators, divergent

Solving the unitriangular system $b_n=\sum_k g_k\binom{n+k}{2k}^2$ gives $g_1=6$,
$g_2=135/4$, $g_3=6095/9$, $g_4=-492065/288$, ... with
$\operatorname{den}(g_k)\mid\operatorname{lcm}(1..k)^3$ for all $k\le30$ (`01_pullback.py`).
Hence
$$d^{T}_{a,b}=\sum_{k\le\min(a,b)}g_k\binom{a+k}{2k}\binom{b+k}{2k},\qquad
\operatorname{den}(d^T_{a,b})\mid\operatorname{lcm}(1..\min(a,b))^3,$$
verified on the whole grid $a,b\le34$: a **min-type** denominator, which divides
$\operatorname{lcm}(1..a)^{r_1}\operatorname{lcm}(1..b)^{r_2}$ for *every* $r_1+r_2=3$ -- the
best conceivable split -- with $d^T_{n,n}=b_n$ exactly.  But $d^T$ is inadmissible:

* $\log|g_k|/(k\log k)\to\approx 2$: $g_k$ grows like $(k!)^2$, so $\sum g_ku^k$ has radius $0$;
* no linear recurrence of order $\le4$, degree $\le5$ annihilates $g_k$: not holonomic;
* $d^T$ fails the $a$-recurrence at 1122 of 1155 grid points, residuals with no rational pattern;
* $\log(|d^T_{a,b}|/c_{a,b})$ is $0.184=\log\zeta(3)$ on $a=b$ but $12.9$ at $(10,20)$ and
  $28.4$ at $(15,25)$: the ray limits blow up.

The same holds for every row (`13_mintype_growth.py`): $\log|g_k|\sim\gamma k\log k$ with
$\gamma\approx1.95$ for $\zeta(3)$, $\zeta(2)$, $s_{10}$ and $\gamma\approx0.99$ for row E --
always divergent, though Catalan's row E is only Gevrey-1.

### 2.4 Harmonic-sum splittings all fail the system

The arithmetic of a split period is *not* the obstruction: the truncated MZV
$Z(a,b)=\sum_{a\ge m>k\ge1,\,k\le b}1/(m^2k)\to\zeta(2,1)=\zeta(3)$ has
$\operatorname{den}Z(a,b)\mid\operatorname{lcm}(1..a)^2\operatorname{lcm}(1..b)$ -- a genuine
$(2,1)$ split (checked for $a,b\le25$).  But none of the natural harmonic-split candidates
satisfies the system (`14_harmonic_splits.py`): $c\,H^{(3)}_a$, $c\cdot\frac12(H^{(3)}_a+H^{(3)}_b)$,
$c\,Z(a,b)$ and $\sum_k\binom{2k}k^2\binom{a+k}{2k}\binom{b+k}{2k}H^{(3)}_k$ all violate the
$a$-recurrence at every one of 169 tested points.  This is forced: the solution space is
2-dimensional and the companion is $s$.

## 3. The two-variable count (Task 3)

### 3.1 Region calculus for the three denominator geometries

Redoing `TWO_VARIABLE_HOLONOMY.md` Prop 2.1 with a general Liouville weight
$w(a,b)=\ell_1a+\ell_2b-\log\operatorname{den}$ and $|R_T|=\kappa T^2$ gives $m\le4\kappa M_1M_2$
(optimise $D_2/D_1=M_1/M_2$).  With $e_k=\ell_k-r_k$:

| denominator type | entry condition | bound $m\le$ |
|---|---|---|
| product $[1..a]^{r_1}[1..b]^{r_2}$ | $\ell_1>r_1$, $\ell_2>r_2$ | $2M_1M_2/(e_1e_2)$ (Prop 2.1) |
| min $[1..\min(a,b)]^{r}$ | $\ell_1,\ell_2>0$ and $\ell_1+\ell_2>r$ | $\frac{2M_1M_2(\ell_1+\ell_2)}{(\ell_1+\ell_2-r)\ell_1\ell_2}$ |
| max $[1..\max(a,b)]^{r}$ | $\ell_1>r$ **and** $\ell_2>r$ | $\frac{2M_1M_2}{\ell_1+\ell_2-r}\left(\frac1{\ell_1-r}+\frac1{\ell_2-r}\right)$ |

Min-type is the only geometry that weakens entry: it asks $\ell_1+\ell_2>r$ where one variable
asks $\ell>r$.  Max-type strengthens it to two separate conditions, i.e. $\ell_1+\ell_2>2r$.  In
the symmetric case the max-type bound is $B\cdot\frac{4M}{2\ell-r}$ with $B=M/(\ell-r)$ the
one-variable bound at the same data, and $\frac{4M}{2\ell-r}>2$ always (Cauchy: $M\ge\ell$):
with equal data the two-variable count is more than twice as expensive.

### 3.2 The conformal ceiling: Landau's theorem forces $\ell_1+\ell_2\le\log16$

On a pullback host the pair $(\varphi_1,\varphi_2)$ must make $h(u)$ single-valued, i.e.
$v=u\circ(\varphi_1,\varphi_2)$ must omit the branch point $u=1/16$.  With
$\alpha(x)=x/(1-x)^2$, $\beta(y)=y/(1-y)^2$ (Koebe, $\alpha'(0)=\beta'(0)=1$),
$A=(\alpha\circ\varphi_1)(\mathbf D)$, $B=(\beta\circ\varphi_2)(\mathbf D)$, the constraint is
$1/16\notin A\cdot B$.  Both $A,B$ are open connected sets containing $0$, so each contains a
disc $D(0,\rho)$, and $D(0,\rho_A)D(0,\rho_B)=D(0,\rho_A\rho_B)$ must miss $1/16$:
$\rho_A\rho_B\le1/16$.  **Landau's theorem** (sharp constant $1/16$, extremal the modular
function): a holomorphic $f$ on $\mathbf D$ with $f(0)=0$ omitting $w_0$ has
$|w_0|\ge|f'(0)|/16$, i.e. its image contains $D(0,|f'(0)|/16)$.  Hence
$\rho_A\ge e^{\ell_1}/16$, $\rho_B\ge e^{\ell_2}/16$, and
$$\ell_1+\ell_2\ \le\ \log16 = 2.7726,$$
sharp.  In one variable the same theorem gives $\ell\le\log(16|t_{\rm outer}|)$ -- CDT's
$16|s|$ -- and $t_{\rm outer}$ can be *large*.  This is the conformal conservation law:
because $u$ is a **product**, the two variables share one Landau budget, and the
"far outer singularity" advantage of the one-variable host is destroyed by the decoupling.

### 3.3 The numbers for Apery's row

| | one variable (Apery row) | two variables (D1 host) |
|---|---|---|
| rate(s) | $\tau=3$ | $\sigma_1=\sigma_2=3$ (max type), total 6 |
| conformal ceiling | $\ell\le\log(16(1+\sqrt2)^4)=6.298$ | $\ell_1+\ell_2\le\log16=2.773$ |
| entry margin | $+3.298$ (passes) | $-3.227$ (**fails**) |
| inventory | row + companion, 1 conditional generator | the same (rank-2 pullback: no new pure functions) |
| bound | $M/3.298$ | none (entry fails); granting entry, $>2\times$ |

In general on such a host the ceiling is $\ell_1+\ell_2\le\log(256|u_0|)$ with $u_0$ the branch
point, so the entry margins are $\log(256|u_0|)-r$ (min type) and $\log(256|u_0|)-2r$ (max type).
For the $\zeta(3)$ host $u_0=1/16$: $-0.227$ and $-3.227$, both failing.  For the rate-2 hosts
$u_0=1/4$ and the margins are $+2.159$ and $+0.159$ (Section 6) -- i.e. **a convergent min-type
companion would pass entry comfortably for $\zeta(2)$, Catalan's row E and Cooper's $s_{10}$.**
That is exactly why the question was worth asking, and exactly where Section 2.3 bites: the
min-type companion is a divergent Gevrey series for all four rows, so the passing entry test has
no functions to apply itself to.

### 3.4 Inventory: two variables buy nothing here

The system is a pullback $\pi^*\mathcal M_u$ of a rank-2 one-variable system; pullback along a
dominant map preserves rank and $\mathbf Q(x,y)$-linear independence, so the pure inventory is
exactly the $u$-line inventory: no super-multiplicativity.  Two remarks:

* **Rate-0 directions are useless.**  If $\sigma_2=0$ and the $y$-radius exceeds 1 then, by
  Polya-Carlson applied to each $F_a(y)=\sum_bL_1(a)c_{a,b}y^b\in\mathbf Z[[y]]$, every
  admissible $f$ is a polynomial in $y$ over $\mathbf Q[[x]]$; writing $f_i=\sum_jg_{ij}(x)y^j$,
  the $f_i$ lie in the $\mathbf Q(x,y)$-span of a $\mathbf Q(x)$-basis of the $g_{ij}$, so
  $m\le\dim_{\mathbf Q(x)}\langle g_{ij}\rangle$ -- the one-variable inventory -- while the bound
  is at least $2B_1B_2\ge2B_1$.  This kills the $(3,0)$-solution $d^{(R)}$ despite its ideal
  denominators.
* The conditional function exists and *is* fold-regular in the two-variable sense:
  $H_{a,b}=s_{a,b}-(\zeta(3)-1)c_{a,b}$ has $\log|H|/(a+b)\to-\log(3+2\sqrt2)=-1.7627$ on the
  diagonal (measured $-1.835$ at $(40,40)$, converging from below), i.e. $H$ is regular across
  the fold branch of $\Gamma$ and singular on the outer branch -- the exact analogue of the
  one-variable $\log|b_n-\zeta(3)a_n|/n\to-\log(17+12\sqrt2)$ (measured $-3.670$ at $n=40$).
  The mechanism is intact; the denominators and the conformal ceiling are what fail.

## 4. Literature (Task 4)

*From knowledge; flagged where less certain.*

* **Beukers' integrals** (1979).  The double integral for $\zeta(2)$ and triple for $\zeta(3)$
  represent the *one-variable* forms; the extra variables are integration variables, not
  variables of a holonomic system, and the denominator $\operatorname{lcm}(1..n)^r$ comes out of
  the single index.
* **Apery numbers as diagonals** (Straub, *Multivariate Apery numbers and supercongruences of
  rational functions*, 2014; Bostan-Lairez-Salvy).  $A_n$ is the diagonal of a rational function
  in 4 variables, and Straub's multivariate Apery numbers $A(n_1,\dots,n_4)$ with
  $A(n,n,n,n)=A_n$ are multi-index extensions of the *row*.  I know of no treatment of a
  multivariate **companion**.  The $c_{a,b}$ above is a two-index specialisation of that family.
* **Appell/Horn systems.**  $F_1$ has rank $3<4$; $F_4$ has rank 4 but by **Bailey's formula**
  is a product of two ${}_2F_1$'s in transformed variables, singular locus the quadric
  $(1-x-y)^2=4xy$.  The host found here is exactly of that type
  ($C=R\cdot{}_2F_1(\frac12,\frac12;1;16\alpha(x)\beta(y))$, quadric $16xy=(1-x)^2(1-y)^2$).
  "Genuinely coupled but a pullback along a product" is the generic outcome for two-variable
  hypergeometrics, and it is why the rank is sub-multiplicative, as
  `TWO_VARIABLE_HOLONOMY.md` Section 4 anticipated.  Integrality of Taylor coefficients of
  Appell functions with rational parameters is governed by Landau-type Gamma-quotient criteria
  (Delaygue; Delaygue-Rivoal-Roques in one variable, with multivariate versions); the integral
  cases are the mirror-map/diagonal ones, integral in *all* indices, i.e. rate $(0,0)$ for the
  row -- consistent with what is found here.
* **Deformed Apery families and Apery limits.**  Recurrence (1.2) is a one-parameter deformation
  of Apery's with $\lambda_1\lambda_2=1$; the fixed-$b$ Apery limits $\zeta(3)-H^{(3)}_b$ are the
  tails of $\zeta(3)$.  This is the pattern of Golyshev-Zagier (*Interpolated Apery numbers,
  quasiperiods of modular forms, and motivic gamma functions*) and Zagier's Apery-limit
  formalism: limits of such families are periods, and members of a family have limits differing
  by rationals.  The tail obstruction is the arithmetic shadow of that.
* **Several-variable transcendence machinery.**  Schneider-Lang in several variables (Bombieri,
  *Algebraic values of meromorphic maps*, 1970) has the same $d^d/d!$ loss as
  `TWO_VARIABLE_HOLONOMY.md` Prop 2.1; Bost's slope method and Bost-Charles are developed on
  curves, and I am not aware of a higher-dimensional arithmetic holonomy theorem of CDT type.
* I am **not aware of any two-variable Apery-type system in the literature with a companion
  having split denominators**, nor of any statement that such a thing is impossible.

## 5. Verdict (Task 5)

**No.**  For Apery's $\zeta(3)$ row -- and for the $\zeta(2)$ row, Zagier's Catalan row E and
Cooper's $s_{10}$ -- the natural non-tensor two-variable extensions do **not** admit a companion
with split denominators.  Precisely:

1. The extensions are genuinely coupled but *degenerate*: pullbacks of a lower-weight
   one-variable system along a multiplicatively separated map $u=\alpha(x)\beta(y)$.  Rank is
   sub-multiplicative, so the pure inventory is at best the one-variable one.
2. The system has a 2-dimensional space of power-series solutions.  Its companion **does** have a
   constant fold period along every ray -- the super-multiplicative ingredient Section 4 of
   `TWO_VARIABLE_HOLONOMY.md` asks for -- but its denominator is
   $\operatorname{lcm}(1..\max(a,b))^r$, rates $(r,r)$ summing to $2r$: worse than the
   one-variable $r$, and the type is an invariant of the system.
3. The min-type (perfectly split) candidate exists as a formal series with the correct diagonal,
   but its hidden-variable coefficients grow like $(k!)^\gamma$: divergent, non-holonomic, no ray
   limits.
4. The obstruction is the **tail obstruction**: a solution with no $b$-denominators has fixed-$b$
   period $\zeta(3)-H^{(3)}_b$, and the rational correction making the period constant costs
   exactly $\operatorname{lcm}(1..b)^r$.  Splitting and period-constancy are complementary, unit
   for unit.
5. The conformal side conspires: Landau gives $\ell_1+\ell_2\le\log16$ on any pullback host along
   a product, against a required $\ell_1+\ell_2>2r$; the one-variable ceiling
   $\log(16|t_{\rm outer}|)$ is unavailable because the decoupling puts fold and outer
   singularity on one irreducible curve and bounds only the *product* of the radii.

The sweep of Section 7 makes this uniform: over ten decouplings of the four rows, **every**
solution of **every** system whose ray limit has a nonzero period coefficient is of MAX type of
rate $r$, or worse (outside the lcm class altogether when the recurrence's leading coefficient
picks up apparent singularities depending on the other variable).  The systems that are *not*
pullbacks along a separated map -- so that the Landau ceiling of Section 3.2 does not apply --
are exactly the ones whose companions are worst arithmetically.  The two halves of the
obstruction cover each other.

What remains open, in decreasing order of plausibility:

1. **The sharpened problem.**  Is there any holonomic two-variable object with
   $\operatorname{lcm}(1..\min(a,b))^2$ denominators and a fold period constant along the fold
   curve?  For the rate-2 hosts such a companion would pass entry by $2.159$ nats -- including
   Catalan's row E, where the one-variable test fails by $0.614$.  The formal candidate exists
   with the right diagonal and the right denominators and is Gevrey-1 divergent; some resummation
   or a different diagonal-preserving extension is the only visible route.
2. **A non-product, non-pullback host.**  All ten systems here satisfy a first-order fibration
   PDE, i.e. all are "one-variable systems in disguise".  A two-variable system with no first
   integral at all (holonomic rank genuinely $>$ the rank of any one-variable pullback) has not
   been exhibited from an Apéry-like row; its inventory could in principle be
   super-multiplicative.  Straub's four-index multivariate Apéry numbers with three or four
   independent indices are the natural place to look.
3. **The conformal extremal problem** raised by Section 3.2 in general: maximise
   $|\alpha'(0)\beta'(0)|$ over holomorphic $\alpha,\beta:\mathbf D\to\mathbf C$ with
   $\alpha(0)=\beta(0)=0$ and $u_0\notin\alpha(\mathbf D)\beta(\mathbf D)$ -- solved here (the
   answer is $256|u_0|$, by Landau) only because the constraint is multiplicative.  For a
   non-separated singular curve the analogous sharp constant is not known to me.

## 6. Landau ceilings and entry margins, host by host (`19_ceilings.py`)

One variable: $\varphi$ omits the OUTER singularity $t_2$, so $\ell\le\log(16|t_2|)$ (Landau;
this is CDT's $16|s|$, and it reproduces the ledger's numbers exactly -- Catalan's $-0.61$ and
CDT's own host $\log16-2=+0.77$).  Two variables on a pullback host $u=\alpha(x)\beta(y)$ with
$\alpha'(0)=\beta'(0)=1$ and branch point $u_0$: $\ell_1+\ell_2\le\log(256|u_0|)$, plus a
one-variable Landau bound for any branch point lying in $x$ or $y$ alone.

| row | one variable: $\ell\le$, $\tau$, margin | two variables: $\ell_1+\ell_2\le$ | min-type margin ($>r$) | max-type margin ($>2r$) |
|---|---|---|---|---|
| $\zeta(3)$ | $6.298$, $3$, $+3.298$ | $\log16=2.773$ | $-0.227$ | $-3.227$ |
| $\zeta(2)$ | $5.179$, $2$, $+3.179$ | $\log64=4.159$ | $+2.159$ | $+0.159$ |
| row E (Catalan) | $1.386$, $2$, $\mathbf{-0.614}$ | $\log64=4.159$ | $+2.159$ | $+0.159$, but see below |
| Cooper $s_{10}$ | $1.386$, $2$, $-0.614$ | (as row E class) | | |

Two refinements for row E, whose decoupling
$c_{a,b}=\sum_k\binom ak\binom{2k}k\binom{2b-2k}{b-k}$ is asymmetric ($w=xy/(1-x)$,
$\alpha=x/(1-x)$, $\beta=y$):

* $y=\tfrac14$ is a branch point of $\sqrt{1-4y}$, in the variable $y$ alone, so Landau gives
  $\ell_2\le\log4=1.386$.  Max-type entry needs $\ell_2>2$: it **fails by
  $2-\log4=0.614$, exactly the one-variable Catalan deficit.**  The two-variable count
  inherits the obstruction unchanged.
* Min-type entry needs only $\ell_1+\ell_2>2$ against a ceiling $\log64=4.159$: it would
  **pass by $2.159$ nats** -- a swing of $\log16$ relative to the one-variable $-0.614$, which
  is exactly the extra Kodaira factor bought by the second variable.

So for Catalan the entire question reduces to one point: **is there a convergent, holonomic
two-variable object with $\operatorname{lcm}(1..\min(a,b))^2$ denominators and a constant fold
period?**  The natural candidate exists formally, has the right diagonal and the right
denominators, and diverges: for row E its hidden-variable coefficients satisfy
$\log|g_k|/k=\gamma\log k+O(1)$ with $\gamma\to0.993$ measured at $k=180$
(`18_rowE_growth.py`), i.e. $g_k\sim k!\,C^{-k}$: Gevrey-1, radius $0$.  (Milder than
$\zeta(3)$'s Gevrey-2, and Borel-summable, but a Borel sum is not analytic at the origin, so it
is not a power series on a polydisc and the count cannot use it.)

## 7. Sweep over ten decouplings (companion computation)

Recurrences fitted mod $p$ over all $1\le a,b\le28$, re-verified exactly on $0\le a,b\le30$
(0 exceptions in 20 recurrences); order minimality checked by ruling out order-2 relations up to
coefficient degree 12; exact coefficients by CRT over 5 primes + rational reconstruction.  The
solution space was computed as the *simultaneous* grid nullspace of both recurrences (robust to
vanishing leading coefficients), dimension stable on grids $N=6,\dots,16$, basis verified exactly
and extended to $N=32$ ($N=64$ for E_D1).  Denominator type measured
normalisation-robustly: type $(r_1,r_2)$ iff $\operatorname{den}(t)/\gcd(\operatorname{den},
\operatorname{lcm}(1..a)^{r_1}\operatorname{lcm}(1..b)^{r_2})$ stays *bounded* as the grid grows.

| case | $r$ | dim | companion denominator type | ray limit $t/c$ | split? |
|---|---|---|---|---|---|
| z3_D1 $\sum\binom ak\binom bk\binom{a+k}k\binom{b+k}k$ | 3 | 2 | **(3,3)**, MAX $r=3$ | $\zeta(3)$, all 7 rays, 20 digits | no |
| z3_D2 $\sum\binom ak^2\binom{b+k}k^2$ | 3 | 3 | **not lcm-type at all** (excess over $\operatorname{lcm}^4\operatorname{lcm}^4$ grows to $4.4\cdot10^{14}$) | $\zeta(3)/2$, all rays | no |
| z3_D3 $\sum\binom ak^2\binom{a+k}k\binom{b+k}k$ | 3 | 3 | **(3,3)**, MAX $r=3$ | $\tfrac27\zeta(3)$, all rays | no |
| z3_D4 $\sum\binom ak\binom bk\binom{a+k}k^2$ | 3 | 3 | **(3,3)**, MAX $r=3$ | $2\zeta(3)$ for $\lambda\le4$; **diverges at $\lambda=5$** | no |
| z2_D1 $\sum\binom ak\binom{2k}k\binom{b+k}{2k}$ | 2 | 2 | **(2,2)**, MAX $r=2$ | $\zeta(2)/2$, all rays | no |
| z2_D2 $\sum\binom ak\binom bk\binom{a+k}k$ | 2 | 2 | **(2,2)**, MAX $r=2$ | $\zeta(2)$, all rays | no |
| z2_D3 $\sum\binom ak^2\binom{b+k}k$ | 2 | 2 | — **no companion exists**: the period functional vanishes on the whole space | limits $0$ and $1$, rational | n/a |
| E_D1 $\sum\binom ak\binom{2k}k\binom{2b-2k}{b-k}$ | 2 | 3 | **not lcm-type at all** (same unbounded excess) | $-G$, $2G$, $1+2G$ on the diagonal and $\lambda\le1$ (11+ digits at $(64,64)$); for $\lambda>1$ drifts or diverges | no |
| s10_D1 $\sum\binom ak^2\binom bk^2$ | 2 | 3 | **(2,2)**, MAX $r=2$ | $\zeta(2)$, $-\tfrac32\zeta(2)$, all rays | no |
| s10_D2 $\sum\binom ak^3\binom bk$ | 2 | 3 | **(2,2)**, MAX $r=2$ | $\zeta(2)/2$, $-\zeta(2)$, all rays | no |

Every primitive direction with $|x_i|\le4$ in each solution space was scanned (24 for dim 2, 289
for dim 3) and classified.  **In all ten cases no solution with a nonzero period coefficient has
a split denominator**: the type is exactly MAX of rate $r$, or worse.

Only z3_D1, z2_D1 and z2_D2 have *both* leading coefficients free of the other variable (the
clean situation of Section 1.2); for the rest the leading coefficient of one or both recurrences
depends on the other variable, and the apparent singularities it introduces (e.g.
$9a^2-3ab+21a-2b^2-6b+12$ for s10_D1, $2a-2b+5$ for E_D1) contribute primes $p>\max(a,b)$ to the
companion's denominator, taking it *outside* the lcm class altogether — worse than max type.

**The decoy.**  In seven of the ten cases the solution space does contain a distinguished
solution whose denominators are split and in fact *better* than split — type $(0,1)$, $(1,0)$ or
$(1,1)$, MAX rate $1$ — but its ratio to the row tends to $0$ along every ray.  Structurally it
is a harmonic-number solution ($t(0,b)=H_b$, supported on a half-plane $b>a$ or $a>b$); it is the
two-variable analogue of the *small* solution, not a companion, and adding any multiple of it
never repairs a companion's denominators (it is $p$-integral exactly at the primes where the
companion has poles).

**Two ancillary facts.**  (i) Ray-independence of the Apéry limit — which holds perfectly for
z3_D1 — is *not* general: z3_D4's ratio diverges at $\lambda=5$, and E_D1's holds only for
$\lambda\le1$.  (ii) z2_D3 shows a decoupling can fail even to contain the companion: the second
solution of that system vanishes identically on the diagonal.

**Closing the last door.**  z3_D2, z3_D3, z3_D4, s10_D1, s10_D2 are *not* pullbacks along a
separated map: they do satisfy a first-order fibration PDE
$A(x,y)\theta_xC+B(x,y)\theta_yC+E(x,y)C=0$ (`21_fibration.py`, `22_firstint.py`) — for s10_D1,
$3A=4x^2y-3x^2-4xy^2-2xy+6x+y^2+2y-3$, $3B=-3A(y,x)$, $E=(x-y)(2xy-x-y+1)$ — but the first
integral of $dx/(xA)=dy/(yB)$ is *not* of the form $\alpha(x)\beta(y)$, so the Landau ceiling of
Section 3.2 does not apply to them.  (s10_D1 is the Legendre generating function
$\sum_kt^kP_k(X)P_k(Y)$, $X=\frac{1+x}{1-x}$, $Y=\frac{1+y}{1-y}$ — an elliptic integral in a
genuinely non-separated modulus.)  This was the one place where the conformal argument could have
been evaded — and it is closed on the arithmetic side instead: their companions are still
MAX-type of rate $r$, or not lcm-type at all.
