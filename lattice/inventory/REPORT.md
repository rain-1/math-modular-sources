# The quantitative inventory bound: how large can $u_2$ and the conditional orbit be, and what does that buy?

*Fable, 2026-09-02. Working directory `lattice/inventory/` (this directory). Sequel to
`consolidation/INVENTORY_BOUND.md` (which proves $u_1=1$) and `consolidation/NUMBER_FIELD_HOLONOMY.md`
(which fixes the number-field normalisation as the average over places). All CDT equation numbers
are those of arXiv:2408.15403v2.*

---

## 0. Summary — the four results

**(R1) CDT's fourteen functions are the complete supply.** On the descent orbifold
$\mathbf P^1_y\setminus\{0,\infty\}$ with a $\mathbf Z/2$ point at $y=4s$, the number of
$\mathbf Q(y)$-independent admissible functions with denominator exponent $\max_ie_i\le 1$ is
exactly $14$, and CDT use all of them, function for function:

| block | measured supply | CDT's functions |
|---|---|---|
| $0$ lcm layers, $e=0$ | $1$ | $B_1=1$ |
| $1$ layer, $e=0$ | $1$ | $B_2$ ($=\mathbf Q(y)$-equivalent to $F_1$) |
| $1$ layer, $e=1$ | $1$ | $B_3=F_2=2\arcsin^2(\sqrt y/2)$ |
| $2$ layers, $e=0$ | $2$ pure $+\,4$ conditional | $B_4,B_5;\ \mathcal G,\mathcal G',\mathcal G'',\mathcal G'''$ |
| $2$ layers, $e=1$ | $2$ pure $+\,3$ conditional | $B_6=F_3,B_7;\ \int\mathcal G,\int\mathcal G\frac{dy}y,\int\mathcal G\frac{dy}{y^2}$ |
| **total** | **14** | **14** |

Both conditional counts are *saturation points*, not choices: $\mathcal G''''$ and
$\int\mathcal G\,dy/y^3$ are $\mathbf Q(y)$-dependent on the rest (§3).

**(R2) The single-layer module is exactly the arcsine tower, one new function per unit of $e$.**
$N(1,e)=e+2$: the functions with at most one $[1,\dots,2n]$ layer and exponent $\le e$ span
$\mathbf Q(y)\{1,F_1,F_2,\dots,F_{e+1}\}$, $F_j:=\sum_{n\ge1}y^n/(n^j\binom{2n}n)$ (rescaled
$y\mapsto y/s$). So $u_2=1+a$ costs the exponent profile $(0,1,\dots,a-1)$ — and that price is
what defeats every host.

**(R3) The entry test is *not* the Catalan obstruction.** With CDT's transported $m=14$ inventory
Catalan's level-8 host misses entry by $-0.0765$ at the hard ceiling (`CDT_FINDER.md` §4A). With the
*best* inventory it passes entry comfortably: $+1.2700$ at the ceiling, $+0.8067$ at CDT's realised
contour. What fails is the **margin**: the best achievable is $\boxed{-8.04}$, i.e. the method is
short by $8$ nats, or by a factor $\approx4.3$ in $m$ ($m\le12.97$ is available; the largest
inventory compatible with entry has $m=3$).

**(R4) The same for $X_1(5)\,\mathrm{Sym}^2$:** best entry $+1.5452$ (ceiling) / $+1.0819$ (CDT's
contour); best margin $\boxed{-6.99}$. The sharp thresholds are: Catalan needs the single-layer
supply to be $\nu=4$ per exponent where it is $1$; $X_1(5)$ needs $\nu=2$, **or** twelve
two-layer $e=0$ pure functions where four exist.

---

## 1. Calibration: CDT's numbers reproduced exactly (`tau.py`, `calib.py`)

An independent exact-rational implementation of Definition 6.0.1 and (6.0.4)–(6.0.6):

$$\tau^\flat(\mathbf b)=\frac1{m^2}\sum_i(2i-1)\sigma_i=\sigma_m-\frac1{m^2}\sum_ju_j^2b_j,\qquad
\tau^\sharp(\mathbf e)=\frac2{m^2}\min_{\xi\in[0,m]}\Bigl\{\xi\sum_ie_i+(\max_ie_i)\,I^m_\xi(\xi)\Bigr\}.$$

| quantity | CDT | computed here |
|---|---|---|
| $\tau^\flat$ ($m=14$, $u=(1,3)$, $b=(2,2)$) | $191/49$ | $\mathbf{191/49}$ (both forms agree) |
| $I_2^{14}(2)$ | $21.075$ | $\mathbf{843/40}$ **exactly** |
| $\tau^\sharp$, $\mathbf e=(0,0,1;0^6;1^5)$ | $27/80$, minimiser $[2,\tfrac{13}6]$ | $\mathbf{27/80}$, and the objective is exactly $27/80$ at $\xi=2,\ \tfrac{25}{12},\ \tfrac{13}6$ and larger at $\xi=\tfrac94$ |
| $\tau(\mathbf b;\mathbf e)$ | $16603/3920=4.235459$ | $\mathbf{16603/3920}$ |
| $\log|\varphi'(0)|=\log(256\cdot0.62922\ldots)$ | $5.081908$ | $5.081908$ |
| bound $\mathrm{BC}/(\log|\varphi'(0)|-\tau)$ | $13.9938$ | $\mathbf{13.9938}$ |
| signed margin | — | $+0.0053$ |

*The exact value $I_2^{14}(2)=843/40$ appears not to be printed by CDT.*

**Contour convention.** CDT's realised loss is the factor $|\psi'(0)|=0.6292232680$, i.e.
$\log|\varphi'(0)|=\text{ceiling}+\log 0.62922=\text{ceiling}-0.46327$. (The task brief wrote
"ceiling $-\,0.62922$"; the finder's `LOSS` and CDT's $5.081908$ both use $\log0.62922=-0.46327$,
which is what is used below. Reading it as a further $-0.62922$ would lower every "entry@contour"
column by $0.16595$ and every margin by $0.166\,m$ — it changes no verdict.)

---

## 2. The single-layer module (`t2_arcsine.py`, `t2_pareto.py`, `t2_polylog.py`, `t2_indep.py`, `t2_words.py`, `t2_supply.py`, `t2_b5b.py`)

Setting: the descent $y=x^2/(x-s)$, $s=1/\lambda_2$, involution $w(x)=sx/(x-s)$, punctures $y=0,\infty$,
order-two point $y=4s$, principal radius $R=4|s|$. All computations below are done at $s=1$; the
transported host only rescales $y\mapsto y/s$ with $\lambda_2=1/s$ an algebraic integer, which multiplies
the $n$-th coefficient by $\lambda_2^n\in\mathcal O_K$ and changes **no** denominator type
(checked independently for $\lambda_2=1,2,4,8$).

### 2.1 The arcsine tower and its sharp type  [VERIFIED to $n=200$; the divisibility is a theorem]

$$F_e(y):=\sum_{n\ge1}\frac{y^n}{n^e\binom{2n}n},\qquad \theta F_e=F_{e-1},\quad F_2=2\arcsin^2(\sqrt y/2).$$

* $\binom{2n}n\mid[1,\dots,2n]$ and, more sharply, $n\binom{2n}n\mid[1,\dots,2n]$ (**Nair's lemma**) —
  no counterexample for $n\le200$; this is a known theorem.
* Hence the Pareto frontier of $F_j$ is $(L\text{ layers},\ e)=(L,\max(0,j-L))$ for $L\ge1$ — verified
  exactly for $j\le10$, $L\le3$, $n\le200$. **In particular every $F_j$ is single-layer**, with sharp
  exponent $n^{j-1}$.
* The rate $2$ is sharp: $\binom{2n}n\nmid[1,\dots,cn]$ for $c<2$ (primes in $(n,2n]$), first failure
  at $n=2,3,6,12$ for $c=1,\tfrac32,1.8,1.9$.

### 2.2 The pure module: what collapses  [PROVED identities, verified numerically]

With $w=sx/(x-s)$ one has $(1-x/s)(1-w/s)=1$, whence

| statement | consequence |
|---|---|
| $\mathrm{Sym}^+\mathrm{Li}_1\equiv0$ | the **only** lcm-free pure candidate is the zero polynomial — an independent confirmation of $u_1=1$ |
| $\mathrm{Sym}^+\mathrm{Li}_2=F_2$ (Landen: $\mathrm{Li}_2(x)+\mathrm{Li}_2(\tfrac{-x}{1-x})=-\tfrac12\log^2(1-x)$) | no new weight-2 direction |
| $\mathrm{Sym}^+\mathrm{Li}_{1,1}=-2F_2$ | idem |
| $(x-w)\,\mathrm{Sym}^-\mathrm{Li}_1\in\mathbf Q(y)\,F_1$ (because $\log(1-x)=\pm2i\arcsin(\sqrt y/2)$) | the "$\mathrm{Sym}^-\log(1-x)/(x-w)$" candidate of the brief is the $F_1$ direction |
| $(4-y)F_0=2F_1+y$ **(exact)** | $F_0$ is $\mathbf Q(y)$-dependent on $\{1,F_1\}$ — this is why CDT never use $F_0$ |
| CDT's $B_3=F_2$, $B_6=F_3$ **(exact coefficient identities)** | their "factorial series" *are* the arcsine tower |
| CDT's $B_2=2\sum y^n/(n(n-1)\binom{2n}n)$ | single layer, $e=0$; $\mathbf Q(y)$-equivalent to $F_1$ |

Rank certificates (mod $p=2^{61}-1$, series to $y^{92}$, ranks of $\{y^jf_i\}_{j\le D}$):

| set | $D=0$ | $1$ | $2$ | $3$ | $\mathbf Q(y)$-dim |
|---|---|---|---|---|---|
| $\{1,F_0,F_1,B_2,(x{-}w)\mathrm{Sym}^-\mathrm{Li}_1\}$ | 5 | 7 | 9 | 11 | **2** |
| $+\{F_2,\mathrm{Sym}^+\mathrm{Li}_2\}$ | 5 | 8 | 11 | 14 | **3** |
| $+\{F_3,F_4,F_5\}$ | 8 | 14 | 20 | 26 | **6** |
| $+\{B_4,B_5,\mathrm{Sym}^+\mathrm{Li}_3\}$ | 12 | 21 | 30 | 39 | **9** |

### 2.3 The supply function $N(L,e)$  [VERIFIED]

$N(L,e):=\dim_{\mathbf Q(y)}\mathrm{span}\{1\}\cup\{f:\ f$ admits a type with $\le L$ layers of
$[1,\dots,2n]$ and exponent $\le e\}$, computed from **all** iterated integrals of weight $\le5$ on
$\mathbf P^1\setminus\{0,s,\infty\}$ (all $31$ admissible words in $\omega_0=\tfrac{dx}x$,
$\omega_1=\tfrac{dx}{x-s}$, both symmetrisations $\mathrm{Sym}^+$ and $(x-w)\mathrm{Sym}^-$: $63$ nonzero
series) together with $F_1,\dots,F_6$ and $B_5$:

| | $e\le0$ | $e\le1$ | $e\le2$ | $e\le3$ | $e\le4$ |
|---|---|---|---|---|---|
| $L=0$ | **1** | 1 | 1 | 1 | 1 |
| $L=1$ | 2 | 3 | 4 | 5 | 6 |
| $L=2$ | 4 | 6 | 7 | 8 | 9 |
| $L=3$ | 8 | — | — | — | — |

* the $L=0$ row is River's theorem ($u_1=1$) — **[PROVED]**;
* $N(1,e)=e+2$: **exactly one new single-layer function per exponent** — the arcsine tower, and nothing
  else, up to weight 4 **[VERIFIED]**;
* $N(2,0)=4$: $\{1,F_1,F_2,B_4\}$ from the words, plus $B_5$, which is an independent fifth direction
  (ranks $[12,16,20]\to$ dim 4 without it, $[15,20,25]\to$ dim 5 with it) — so **two** non-arcsine
  atoms at $(2,0)$;
* $N(L,e)\approx(2^L-1)+e+1$: the "$+e$" is the arcsine tower, the $2^L-1$ is the weight tower.

**Reading.** The only lever on $u_2$ is the arcsine tower, and it is a *ladder*: the $a$-th rung costs
exponent $a-1$.

---

## 3. The conditional orbit (`t3_orbit.py`, `t3_ode.py`, `t3_int.py`, `t3_verify.py`, `t3_x15.py`)

$\mathcal G=\mathrm{Sym}^+H$, $H=aH_A+bH_B+cH_C$ the hypothesis-fixed solution of CDT's Prop. 11.1.4
(rebuilt from their ODE; $H_A=1,3,15,93$, $H_C$ has $\tfrac{343}9$ — matches their printed data).
Series to $y^{160}$, ranks mod $2^{61}-1$, relations searched to polynomial degree $40$.

| question | answer | evidence |
|---|---|---|
| how many $\theta$-derivatives are $\mathbf Q(y)$-independent? | **4** ($\mathcal G,\mathcal G',\mathcal G'',\mathcal G'''$) | $\{1,\mathcal G,\dots,\mathcal G'''\}$ independent; $\{1,\mathcal G,\dots,\mathcal G''''\}$ has a relation of degree $\le9$ |
| how many integrations $\int\mathcal G\,dy/y^i$ ($e=1$ each)? | **3** ($i=0,1,2$) | $i\le2$ independent; adding $i=3$ gives a relation of degree $\le9$ |
| beyond that? | iterated integrals $\int\!\!\int$ are independent, at $e=2$ | verified |

So the **minimal inhomogeneous $y$-ODE of $\mathcal G$ has order 4**, and the orbit with $\max e_i\le1$
has exactly $7$ members — *precisely CDT's choice*. The structural reading:

$$\#\{\text{derivatives}\}=2\cdot\mathrm{ord}(\text{homogeneous }x\text{-ODE}),\qquad
\#\{\text{$e=1$ integrations}\}=\dim H^1_{\mathrm{dR}}\ (=3\text{ here}),$$

the factor $2$ being the degree of the descent $x\to y$ (the $y$-module is the pushforward of the
$x$-module, of rank $2\cdot2=4$; the "$+1$" for the rational inhomogeneity gives the order-4 *inhomogeneous*
equation on $\{1,\mathcal G,\dots,\mathcal G'''\}$).

**Weight-2 rows.** For the $X_1(5)$ $\mathrm{Sym}^2$ host: the weight-1 (Zagier D) row satisfies an
order-2 ODE and its Cauchy square $A_n=[t^n]F(t)^2=1,6,47,408,3745,\dots$ satisfies an **order-3** ODE
(coefficients of degree $\le6$) — the symmetric square, verified. Hence $D=6$ derivatives on the
$y$-line, as `CDT_FINDER.md` §6 guessed. The integration count is not computed (the conditional
function itself was not built); it is scanned as a parameter below and **makes no difference**.

---

## 4. The optimisation (`optimise.py`, `optimise2.py`, `final.py`, `frontier.py`, `thresh.py`, `tradeoff.py`)

An inventory is a multiset of $(L,e)$; $\tau$ depends on it only through $(m,u_1,\dots,u_k,\sum e_i,\max e_i)$,
$$\tau^\flat=2k-\frac2{m^2}\sum_{j=1}^k u_j^2,\qquad u_j=\#\{i:L_i<j\},$$
all $b_j=2$ (the descent doubles every layer). Supply as in §2.3 and §3; exponents assigned greedily.
Hosts: $\mathrm{ceiling}=\log(256/\lambda_2^{\rm norm})$, contour $=\mathrm{ceiling}+\log0.62922$,
$\mathrm{BC}=11.845+\log s$.

**Control.** The optimiser, given the *measured* supply on CDT's host with one conditional generator,
returns **exactly CDT's inventory** as the maximiser of the margin: $m=14$, $a=2$, top block
$6$ at $e=0$ and $5$ at $e=1$, $\tau=16603/3920$, margin $+0.0053$.

### 4.1 Best entry / best margin, by host and scenario

Scenarios: **B** = measured supply ($\nu=1$ single-layer per exponent, $q_2=2$ extra pure atoms at
$(2,0)$, $4$ pure at $(3,0)$, $D$ derivatives $+\,I$ integrations per generator);
**C** = measured $\nu=1$ but unbounded supply of $e=0$ atoms at the top level;
**D** = hypothetical $\nu=2,3$; **E** = $\nu=\infty$ (i.e. $u_2=m-1$: forbidden by §2.3).

*Scenario C is listed only to locate the binding constraint.* On H1 and H3 it returns margins of
$+107$ and $+103$ at $m=110$ — i.e. it would prove vastly more than CDT did — because it posits an
unbounded supply of top-level $e=0$ functions where §2.3 measures $N(2,0)=4$ and §3 caps the
conditional block at $4g$. Scenario B is the honest row, and on H1 it *is* CDT. On H2 scenario C
changes nothing at all, because there entry, not supply, is what stops $m$ from growing.

| host | scenario | best entry @ ceiling | best entry @ CDT contour | **best margin** | inventory at the margin optimum |
|---|---|---|---|---|---|
| **H1** CDT, $k=2$, $s=1$ | B, $g=1$ | $+2.6563$ | $+2.1930$ | $\mathbf{+0.0053}$ | $m=14$ = CDT's own |
| | B, $g=3$ | $+2.6563$ | $+2.1930$ | $+21.28$ | $m=50$ (3 generators) |
| | E (absurd) | $+3.4522$ | $+2.9889$ | $+137.7$ | — |
| **H2** Catalan lvl 8, $k=2$, $s=\tfrac14$ | **B, any $g$** | $\mathbf{+1.2700}$ | $\mathbf{+0.8067}$ | $\mathbf{-8.039}$ | $m=3$: $\{1,\ F_1(4y),\ \mathcal G\}$ |
| | C | $+1.2700$ | $+0.8067$ | $-8.039$ | unchanged |
| | D, $\nu=2$ | $+1.4089$ | $+0.9456$ | $-6.368$ | $m=6$ |
| | D, $\nu=3$ | $+1.5189$ | $+1.0556$ | $-3.789$ | $m=11$ |
| | E (absurd) | $+2.0659$ | $+1.6026$ | $+56.85$ | $m=42$, $a=40$ |
| **H3** $X_1(5)\,\mathrm{Sym}^2$, $k=3$ | **B, any $g$** | $\mathbf{+1.5452}$ | $\mathbf{+1.0819}$ | $\mathbf{-6.986}$ | $m=7$: $1+2$ single-layer $+3$ two-layer $+\mathcal G$ |
| | D, $\nu=2$ | $+1.9202$ | $+1.4569$ | $-6.017$ | $m=12$ |
| | D, $\nu=3$ | $+2.2092$ | $+1.7460$ | $+2.12$ | $m=8$ |
| | E (absurd) | $+3.3581$ | $+2.8948$ | $+133.8$ | — |

(For reference, CDT's $m=14$ inventory *transported*: H2 entry $-0.0766$/$-0.5398$, margin $-18.02$;
H3 entry $-0.4352$/$-0.8984$, margin $-24.42$ — reproducing `CDT_FINDER.md` §4A exactly.)

### 4.2 The frontier: best entry and margin at each $m$ (measured supply, $g=1$)

| $m$ | H1 entry@ceil | H1 margin | H2 entry@ceil | H2 entry@contour | H2 margin | H3 entry@ceil | H3 entry@contour | H3 margin |
|---|---|---|---|---|---|---|---|---|
| 3 | $+2.656$ | $-5.27$ | $\mathbf{+1.270}$ | $\mathbf{+0.807}$ | $\mathbf{-8.04}$ | $+1.545$ | $+1.082$ | $-8.60$ |
| 5 | $+1.995$ | $-4.19$ | $+0.609$ | $+0.146$ | $-9.73$ | $+1.275$ | $+0.812$ | $-7.79$ |
| 7 | $+1.749$ | $-2.84$ | $+0.363$ | $-0.100$ | $-11.16$ | $+1.157$ | $+0.694$ | $\mathbf{-6.99}$ |
| 10 | $+1.460$ | $-1.88$ | $+0.074$ | $-0.389$ | $-14.35$ | $+0.870$ | $+0.407$ | $-7.78$ |
| **14** | $+1.310$ | $\mathbf{+0.01}$ | $-0.077$ | $-0.540$ | $-18.02$ | $+0.471$ | $+0.008$ | $-11.74$ |
| 19 | $+1.082$ | $-0.08$ | $-0.304$ | $-0.767$ | $-25.04$ | $+0.093$ | $-0.371$ | $-18.89$ |
| 24 | $+0.911$ | $-1.10$ | $-0.475$ | $-0.939$ | $-32.99$ | $-0.279$ | $-0.742$ | $-29.65$ |
| 40 | $+0.002$ | $-30.29$ | $-1.384$ | $-1.847$ | $-84.36$ | $-0.838$ | $-1.302$ | $-63.91$ |

**$m=14$ is the unique $m$ with a positive margin on CDT's host**, and it is positive by $0.01$.

### 4.3 Why the single-layer ladder cannot be climbed

Leanest inventory: $1+F_1+\dots+F_a+\mathcal G$, $m=a+2$, forced exponents $(0,1,\dots,a-1)$:

| $a$ | $m$ | $\tau^\flat$ | $\tau^\sharp$ | $\tau$ | H1 entry@ceil | H2 entry@ceil | H2 margin |
|---|---|---|---|---|---|---|---|
| 0 | 2 | 3.0000 | 0 | 3.0000 | $+2.545$ | $+1.159$ | $-9.07$ |
| **1** | **3** | **2.8889** | **0** | **2.8889** | $\mathbf{+2.656}$ | $\mathbf{+1.270}$ | $\mathbf{-8.04}$ |
| 2 | 4 | 2.7500 | 0.4375 | 3.1875 | $+2.358$ | $+0.971$ | $-8.43$ |
| 4 | 6 | 2.5556 | 1.3750 | 3.9306 | $+1.615$ | $+0.228$ | $-11.87$ |
| 8 | 10 | 2.3600 | 2.7708 | 5.1308 | $+0.414$ | $-0.972$ | $-24.81$ |
| 16 | 18 | 2.2099 | 4.4039 | 6.6138 | $-1.069$ | $-2.455$ | $-62.99$ |

$\tau^\flat$ falls by $\le0.83$ in total (from $3$ to its infimum $2$) while $\tau^\sharp$ grows without
bound like $\log m$: **the second rung of the ladder already costs more than it saves at small $m$**,
and every rung after the second costs more than it saves at any $m$. This is the quantitative form of
River's expectation, sharper than "a few tenths": the optimum is $a=1$ for a lean inventory, $a=2$ once
the top block is large enough to amortise a single $e=1$.

### 4.4 What would have to be true

| host | lever | measured value | value needed for margin $>0$ |
|---|---|---|---|
| H2 Catalan | $\nu$ = single-layer functions per exponent | $\mathbf 1$ | $\mathbf 4$ (margin $+0.43$ at $m=18$, $a=16$) |
| H3 $X_1(5)$ | $\nu$ | $\mathbf 1$ | $\mathbf 2$ (margin $+0.13$ at $m=12$) |
| H3 $X_1(5)$ | $q_2$ = pure atoms with $\le2$ layers and $e=0$ | $\mathbf 4$ | $\mathbf{12}$ (margin $+1.29$ at $m=18$) |
| H2 Catalan | extra pure atoms at the *top* level | any | **never helps** (top level $=\sigma_m$: raises $m$, not $u_2$) |
| H2, H3 | conditional orbit size ($D$ derivatives, $I$ integrals, $g$ generators) | $D{=}4/6$, $I{=}3$, $g{\le}3$ | **doubling $D$ and $I$ changes nothing** (the optimum already sits below the supply cap) |

---

## 5. Verdict (Task 5)

**Is there any inventory of admissible functions that passes ENTRY?**

* **H2 (Catalan, level 8, $\lambda_2=4$): YES, comfortably.** Maximum entry $+1.2700$ at the hard
  ceiling $\log64$ and $+0.8067$ at CDT's realised contour, achieved by
  $$\{\,1,\quad F_1(4y)=\textstyle\sum_n\frac{4^ny^n}{n\binom{2n}n},\quad \mathcal G=\mathrm{Sym}^+H\,\}
  \qquad(m=3,\ u_1=1,\ u_2=2,\ \mathbf e=0,\ \tau=\tfrac{26}9).$$
  Entry stays positive at the ceiling up to $m=11$ and at CDT's contour up to $m=5$.
  **The finder's "$-0.0765$ at the ceiling" is an artefact of transporting CDT's $m=14$ array; the
  entry test is not the Catalan obstruction.**
* **H3 ($X_1(5)\,\mathrm{Sym}^2$, $k=3$): YES.** Maximum entry $+1.5452$ (ceiling), $+1.0819$ (contour),
  at $m=3$; the best margin configuration ($m=7$) still has entry $+1.1574$/$+0.6942$.
  This also corrects `CDT_FINDER.md` §6, whose $-0.435$/$-0.898$ were CDT-proportional-inventory numbers.

**Is there any inventory that passes the MARGIN? NO — on either host, by a wide and stable amount.**

| host | best margin over all inventories | in per-function units | what it would take |
|---|---|---|---|
| H2 Catalan level 8 | $\mathbf{-8.04}$ | the bound at that inventory is $m\le12.97$; the inventory has $m=3$ | $\nu=4$ single-layer functions per exponent |
| H3 $X_1(5)\,\mathrm{Sym}^2$ | $\mathbf{-6.99}$ | the bound at that inventory is $m\le17.1$; the inventory has $m=7$ | $\nu=2$, or $q_2=12$ |

The obstruction has moved but not weakened. It is now stated as a **dimension deficit in the pure
module**, and it is the same deficit on both hosts:

> On the descent orbifold the entire supply of admissible functions with $\max_ie_i\le1$ is $14$
> (one constant, two single-layer, $2+2$ pure at the top level, $4+3$ conditional per generator).
> CDT's host has enough archimedean room ($\log256$) to convert those $14$ into a contradiction with
> $0.005$ to spare. Catalan's host loses $\log4=1.386$ of ceiling and, with it, $14\cdot1.386\approx19$
> of margin — a loss no rearrangement of a $14$-function inventory can recover, and one that the
> arcsine ladder cannot pay for because each extra rung costs $\Delta\tau^\sharp>|\Delta\tau^\flat|$.

The one honest caveat is that the "$\nu=1$" statement (one new single-layer function per exponent) is
verified over the weight-$\le4$ symmetrised iterated-integral module plus the arcsine tower plus CDT's
factorial series, not proved. A *proof* would be the natural sequel to `INVENTORY_BOUND.md`: an
analogue of Pólya–Carlson for series with a single lcm layer.

---

## 6. Proved / verified / assumed

**[PROVED]**
1. $u_1=1$ (River, `INVENTORY_BOUND.md` Thm 2.1/Cor 2.2). Independently corroborated here:
   $\mathrm{Sym}^+\mathrm{Li}_1\equiv0$ because $(1-x/s)(1-w/s)=1$.
2. $n\binom{2n}n\mid[1,\dots,2n]$ (Nair) $\Rightarrow$ $F_j$ is single-layer with sharp exponent $n^{j-1}$.
3. $(4-y)F_0=2F_1+y$ exactly (at $s=1$; it rescales with $y\mapsto y/s$), so
   $F_0\in\mathbf Q(y)\{1,F_1\}$ — the arcsine tower loses its bottom rung.
4. $\mathrm{Sym}^+\mathrm{Li}_2=F_2$ (Landen); $(x-w)\mathrm{Sym}^-\mathrm{Li}_1\in\mathbf Q(y)F_1$;
   $B_3=F_2$, $B_6=F_3$ as coefficient identities.
5. The formulas (6.0.4)–(6.0.6) as implemented — reproduce CDT's $191/49$, $843/40$, $27/80$,
   $16603/3920$, $13.9938$.

**[VERIFIED — exact computation, no proof]**
6. Sharp denominator types of $F_e$ ($n\le200$, $\lambda_2\in\{1,2,4,8\}$, $L\le3$, $e\le10$), of
   $\mathrm{Sym}^\pm\mathrm{Li}_j$ ($j\le4$) and of $B_2,B_3,B_5,B_6$.
7. The supply table $N(L,e)$ over all weight-$\le5$ symmetrised iterated integrals ($63$ series),
   ranks mod $2^{61}-1$ to $y^{80}$; $B_5$ is a fifth $(2,0)$ direction.
8. The conditional orbit saturates at $4$ derivatives and $3$ integrations (series to $y^{160}$,
   relations to degree $40$); double integrations are independent at $e=2$.
9. The $X_1(5)$ $\mathrm{Sym}^2$ row satisfies an order-3 ODE (symmetric square); the weight-1 row order 2.

**[ASSUMED — flagged, and each is scanned as a parameter]**
10. **Admissibility.** That every $F_j$ and every symmetrised iterated integral is admissible on the
    orbifold, i.e. that $\varphi^*f$ is meromorphic on $\mathbf D$. They are holomorphic on
    $\mathbf P^1_y\setminus\{4s,\infty\}$ with $\mathbf Z/2$ monodromy at $4s$; $F_0,F_1$ have a pole
    in the local uniformiser $\sqrt{4s-y}$, exactly like CDT's own $\mathcal G',\mathcal G'',\mathcal G'''$,
    which is why this is taken to be allowed. *No independence proof in CDT's sense (their Lemma 12.1.1)
    is offered here for any inventory other than CDT's own.*
11. **$\nu=1$** — that there is exactly one new single-layer direction per exponent. Verified to weight 4;
    scanned at $\nu=2,3,\infty$ (§4.1, §4.4).
12. **Transport of the contour.** $\log|\varphi'(0)|=\mathrm{ceiling}+\log0.62922$ and
    $\mathrm{BC}=11.845+\log s$ on every host — `CDT_FINDER.md` §8's estimate 1, unchanged and untested.
    CDT's own improved numerators ($13.62$ vs $13.845$, i.e. $\mathrm{BC}\rightsquigarrow11.53$) would
    move the H2/H3 margins by $\approx+0.3$ and $+0.3$: no verdict changes.
13. **The number-field normalisation for H3** — the average over places, per `NUMBER_FIELD_HOLONOMY.md`.
14. **$D=6$ for weight-2 rows**, and the integration count for H3 (scanned; irrelevant, §4.4).

---

## 7. Scripts

| file | contents |
|---|---|
| `tau.py` | exact-rational $I_u^v(w)$, $\tau^\flat$ (both forms), $\tau^\sharp$ |
| `fasttau.py` | float version used inside the scans, cross-checked against `tau.py` |
| `calib.py` | Task 1: CDT's numbers reproduced |
| `qseries.py`, `dtypes.py` | exact $\mathbf Q$-series arithmetic, descent $x\to y$, denominator-type tests |
| `t2_arcsine.py`, `t2_pareto.py` | the arcsine tower's sharp types |
| `t2_polylog.py`, `t2_indep.py`, `t2_b5b.py` | symmetrised polylogs, CDT's $B_i$, the collapse identities |
| `t2_words.py`, `t2_supply.py` | all weight-$\le5$ iterated integrals; the supply table $N(L,e)$ |
| `t3_orbit.py`, `t3_ode.py`, `t3_int.py`, `t3_verify.py` | the conditional orbit and its two saturation points |
| `t3_x15.py` | the $X_1(5)$ $\mathrm{Sym}^2$ ODE order |
| `optimise.py`, `optimise2.py`, `final.py` | the inventory optimiser and the scenario tables |
| `frontier.py`, `thresh.py`, `tradeoff.py` | frontier by $m$, supply thresholds, the ladder trade-off |
