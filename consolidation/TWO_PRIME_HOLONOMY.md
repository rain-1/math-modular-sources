# Two-prime overconvergence and the adelic holonomy bound

*Fable (Opus 5), 2026-08-22.  Scripts and logs: `lattice/two_prime_holonomy/`.
Calibrated on and reusing `lattice/cdt_finder/` (`CDT_FINDER.md`),
`lattice/adelic_holonomy/` (`ADELIC_HOLONOMY.md`), `lattice/multi_prime/`
(`MULTI_PRIME_LATTICE.md`), `lattice/euler_criterion/` (`EULER_CRITERION.md`).
Literature in source: CDT, arXiv:2408.15403v2 §§6–8, 11–13, App. A.
Tags: **[proved]** = exact identity or argument; **[verified]** = exact
rational / $p$-adic computation over a stated range; **[measured]** = numerical
rate at finite index; **[open]**.
**No irrationality claim is made anywhere in this note.***

---

## 0. Verdict

| claim | verdict |
|---|---|
| **All fourteen of CDT's functions have $p$-adic slope $0$ at $p=2$ and at $p=3$.** $v_p$ of the $y$-coefficients is $O(\log n)$: every least-squares slope **halves** when the window doubles from $[100,400]$ to $[200,800]$ — the signature of a logarithm — and lies in $[-0.013,+0.003]$ | **[verified]** to $y$-order $830$ ($x$-order $1670$), §2.1 |
| Hence $\gamma_2=\gamma_3=0$ and $\tau_{\mathrm{ad}}=\tau$: **the adelic margin on CDT's own host is exactly CDT's, $+0.0053$**, unchanged to the last digit | **[verified]**, §2.3 |
| The host *does* carry a $3$-adic resource, $\sigma_3=v_3(9)=2$ — but it lives **only** in the direction $B-\xi_3A$, $\xi_3=\tfrac12\zeta_3(2)$, whose coefficients are $3$-adically irrational. A scan of all $865$ rational $(a,b,c)$ of height $\le6$ for the conditional function, and all rational $r$ of height $\le6$ for $B-rA$, finds **no** direction with positive slope; $B-\xi_3A$ has slope $1.99759\to2$ | **[verified]**, §2.2 |
| Were that resource realisable with rational coefficients on $u=7$ of the $14$ functions it would be worth $\gamma_3=+0.807$ and **$+11.3$ of margin** — an order of magnitude more than everything else in this file. It is not realisable | **[verified]** (the counterfactual), §2.4 |
| **Adding a well-poised decayer's generating function to the space contributes a *negative* $\gamma_2+\gamma_3$**, between $-0.25$ and $-5.69$, because two-prime-ness in the lattice sense is carried by $p$-power **denominators** ($\kappa_p>0$), which the holonomy bound reads as negative overconvergence slopes | **[proved]** + **[verified]**, §3 |
| It is also **not the same holonomic module and not the same singular set**. Exact minimal joint $(Q,P)$ recurrences were fitted over $\mathbf Q$ (nullity $1$, verified on held-out indices): at $\alpha=2$ the conductor-6 family has $\chi(x)=(x-1)(x+1024)$, **identical to conductor 3**, so $\Lambda=2^{10}$ exactly and the singular set is $\{0,1,-2^{-10},\infty\}$ against CDT's $\{0,\tfrac19,1,\infty\}$; the conductor-12 product form has $\Lambda=2^{20}$. Contour cost from the three-point bound $\lvert\varphi'(0)\rvert\le256\lvert q\rvert$: $\log\lvert\varphi'(0)\rvert\le-1.386$, entry $\le-5.62$ | **[verified]**, §3.1a, §3.3 |
| At $\alpha=1$, by contrast, $\deg c_i=(9,11,13)$: the Newton polygon has slope $2$, the coefficients decay factorially and **the generating function is entire** — no contour cost at all, and the entire cost is the denominators. This closes `MULTI_PRIME_LATTICE.md` §8 open item 1 at $\alpha\in\{1,2\}$, and gives the conductor-6 $Q_b$ in closed hypergeometric form | **[verified]**, §3.1a |
| The cheapest member of any family (conductor 6, $\alpha=1/4$) still costs $\mathbf{1.82}$ of entry against a headroom of $1.31$ at CDT's hard ceiling. **Every member of every family fails the entry test on every host tested** | **[verified]**, §3.4 |
| **A conductor-12 two-prime Catalan row exists and is new.** Poles at twelfths, four classes, one antisymmetry and one extra condition; cleanest in the *product* form $Q_b=160\,s^{(1)}_bs^{(5)}_b$. Exact identity $\sum_{t\ge0}R(t)=Q_bG-P_b$ | **[proved]** + **[verified]** to $10^{-75}$, §4.1 |
| Its two $p$-adic avatars are **exactly** the Theorem-F predictions, and they mirror the conductor-6 row for $\chi_{-3}$ term for term: $$\xi_2=\zeta_2(2),\qquad \xi_3=\tfrac9{10}L_3(2,\chi_{12}).$$ Verified to $v_2=463$ and $v_3=176$ at $b=60$ (Cauchy $451$/$168$), with every competing rational scalar giving $v\le1$ | **[verified]**, §4.4 |
| $\sigma_2=8$, $\sigma_3=3$ per unit of $b$ — **the same profile as the conductor-6 row**; $\kappa_2=8\alpha$ exactly | **[measured]**, §4.3 |
| **It has nothing to bridge against.** $\sigma_3=0$ on Zagier $\mathbf E$ and on Zudilin's row ($v_3$ of the increments is $0$ or $-2$ out to $n=400$), so `MULTI_PRIME_LATTICE.md` §8 open item 7 is now settled in the negative half | **[verified]**, §4.5 |
| **Putting it into the Catalan holonomy bounds makes them much worse:** level 8, $-7.97\to\mathbf{-67.4}$; level 16, $-7.55\to\mathbf{-67.2}$. In the per-function units of `CDT_NONCONGRUENCE.md`, deficit $-0.613\to-4.81$ and $-0.539\to-4.48$ | **[verified]**, §5 |
| **The pure-module multi-prime identity.** For the polylogarithm module $\mathrm{Li}_j(x/s)$, $\sum_p\varsigma_p\log p=\log\lvert\lambda_2\rvert$ **exactly**, so $$\underbrace{\log\tfrac{256}{\lvert\lambda_2\rvert}}_{\text{archimedean}}+\underbrace{\sum_p\varsigma_p\log p}_{\text{all primes}}=\log256$$ on *every* host. The two halves enter the entry condition at weights $1$ and $(u-1)^2/m^2=0.1837$, so $\lambda_2=\pm1$ is optimal and **a second prime in $\lambda_2$ is never worth anything** | **[proved]** + **[verified]**, §6 |

One sentence: *there are no two-prime overconvergent functions to carry into the
bound — $\sigma_p>0$ is a Casoratian invariant of a **pair** $(A,B)$, realised
only in a direction with $p$-adically irrational coefficients, whereas the
holonomy bound sees only individual $\mathbf Q$-rational functions; the one
genuinely new object produced here, a conductor-12 Catalan row with avatars at
both $2$ and $3$, confirms Theorem F to $463$ and $176$ digits and costs
$-60$ of margin if you try to use it.*

---

## 1. The distinction that decides everything

Two different numbers are both called "$p$-adic" in this project and they are not
the same thing.

* $\sigma_p$ — the **slope of a row**: $\sigma_p>0$ iff
  $v_p\bigl(b_n/a_n-b_{n-1}/a_{n-1}\bigr)$ grows linearly, i.e. iff the $p$-adic
  Apéry limit $\xi_p=\lim_pb_n/a_n$ exists.  It equals $w_p+2\kappa_p$
  (`ONE_CLASS_TWO_WORLDS.md` §6.1) and for an integral row $v_p(c)$
  (`CRYSTAL_THEOREM_F.md` Cor. 3.2).  It is an invariant of the **pair**
  $(A,B)$ — of the Casoratian.
* $\varsigma_p$ — the **overconvergence slope of one function**:
  $v_p(c_n)\ge\varsigma_pn-o(n)$ for the coefficients of a single
  $f\in\mathbf Q[\![x]\!]$.  This, and only this, is what
  $\gamma_p$ of `ADELIC_HOLONOMY.md` §2.2 is built from.

The bridge between them is exactly one function: $B-\xi_pA$, which has
$\varsigma_p=\sigma_p$.  Its coefficients are $\xi_p$-linear combinations of
rationals, and $\xi_p$ is (conjecturally, and for $p=2,3$ in the two cases at
issue by Calegari 2005) **irrational**.  So:

> **Principle.**  A row's $p$-adic slope is never directly usable in the
> holonomy bound.  The usable positive slopes are exactly two: the **pure**
> polylogarithm module ($\varsigma_p=v_p(\lambda_2)$, unconditional, rational)
> and a **doubly-small** class ($\varsigma_p>0$, rational, unconditional —
> `CATALAN_TWO_CLASSES.md` §3).  Everything else has $\varsigma_p\le0$.

This note is, in the end, that principle made quantitative at two primes.

---

## 2. Part (a): CDT's own function space at $p=2$ and $p=3$

### 2.1 The measurement

`01_lib.gp`, `01_main.gp`, `01_main_big.gp` (logs `out/01_*.log`).  CDT's
fourteen functions were rebuilt **exactly over $\mathbf Q$** — a port of
`lattice/cdt_finder/indep_check2.py`, which works mod $2^{61}-1$ and therefore
destroys all $p$-adic information.  Validation: $H_A=1,3,15,93,639,4653$;
$H_B$'s $x^2$-coefficient $23/4$; $H_C$'s $x^3$-coefficient $343/9$;
$B_4=-4y+\tfrac49y^2+\tfrac{31}{900}y^3$; and the exact $G,B_4,B_7$ and the two
$\int$'s reduce mod $2^{61}-1$ to `indep_check2.py`'s output digit for digit.

$y$-order $830$ ($x$-order $1670$).  "LS" = least-squares slope of $v_p(c_n)$
against $n$ on the window; "$\alpha$" = the fitted coefficient in
$v_p(c_n)\approx\alpha\log_pn$.

| $f$ | LS$_{[100,400]}(2)$ | LS$_{[200,800]}(2)$ | $\alpha_2$ | LS$_{[100,400]}(3)$ | LS$_{[200,800]}(3)$ | $\alpha_3$ |
|---|---|---|---|---|---|---|
| $B_1=1$ | — | — | — | — | — | — |
| $B_2$ | $-0.00284$ | $-0.00140$ | $-0.43$ | $-0.00115$ | $-0.00083$ | $-0.55$ |
| $B_3$ | $-0.00297$ | $-0.00145$ | $-0.44$ | $-0.00122$ | $-0.00081$ | $-0.54$ |
| $B_4$ | $-0.00960$ | $-0.00476$ | $-1.51$ | $-0.00460$ | $-0.00282$ | $-1.62$ |
| $B_5$ | $-0.00284$ | $-0.00142$ | $-0.44$ | $-0.00123$ | $-0.00080$ | $-0.53$ |
| $G$ | $-0.00979$ | $-0.00445$ | $-1.45$ | $-0.00459$ | $-0.00365$ | $-1.90$ |
| $G'$ | $-0.00955$ | $-0.00437$ | $-1.42$ | $-0.00472$ | $-0.00363$ | $-1.89$ |
| $G''$ | $-0.00931$ | $-0.00431$ | $-1.40$ | $-0.00497$ | $-0.00356$ | $-1.85$ |
| $G'''$ | $-0.00907$ | $-0.00420$ | $-1.36$ | $-0.00502$ | $-0.00349$ | $-1.80$ |
| $B_6$ | $-0.00310$ | $-0.00149$ | $-0.45$ | $-0.00122$ | $-0.00081$ | $-0.54$ |
| $B_7$ | $-0.00984$ | $-0.00479$ | $-1.52$ | $-0.00447$ | $-0.00285$ | $-1.64$ |
| $\int G$ | $-0.01009$ | $-0.00448$ | $-1.45$ | $-0.00453$ | $-0.00369$ | $-1.93$ |
| $\int(G-G_0)/y$ | $-0.01009$ | $-0.00446$ | $-1.45$ | $-0.00447$ | $-0.00370$ | $-1.94$ |
| $\int(G-G_0-G_1y)/y^2$ | $-0.01015$ | $-0.00448$ | $-1.45$ | $-0.00447$ | $-0.00369$ | $-1.93$ |

**Every slope halves when the window doubles.**  That is the exact signature of
$v_p(c_n)=\alpha\log n+O(1)$, i.e. **$\varsigma_p=0$** for all fourteen at both
primes.  Wide-window slopes $(v_p(c_{400})-v_p(c_{100}))/300$ lie in
$[-0.027,+0.003]$ and shrink to $[-0.013,+0.003]$ on $[200,800]$;
$\max_nv_p(c_n)/n\le0$ throughout.  Two clean clusters appear: the
hypergeometric $B_2,B_3,B_5,B_6$ at $\alpha\approx-0.44\log_2n$,
$-0.54\log_3n$; everything touching $\mathrm{Li}_2$ or $G$ at
$\approx-1.4\log_2n$, $-1.9\log_3n$.

In the $x$-coordinate, to $n\le800$: $v_2(a_n)=0$ for **every** $n$ ($H_A$'s
coefficients are all odd); $v_3(a_n)=3,6,8,8$ at $n=100,200,400,800$
($\alpha=+1.17\log_3n$); $H_B$ has $v_2=-12,-14,-16,-18$
($-2.15\log_2n$); $H_C$ has $v_3=-6,-4,-3,-5$.  All zero slope.

### 2.2 Where the $3$-adic resource actually is

CDT's host is Zagier's row $\mathbf C$, $(\alpha,\gamma,\delta)=(10,3,9)$,
$c=9$.  `01_host.gp`, `01_host_wide.gp`:

* $v_3(b_n/a_n-b_{n-1}/a_{n-1})/n=1.760,1.930,1.9667,1.9575,1.980,1.97625$ at
  $n=50,100,300,400,500,800$: $\sigma_3=2=v_3(9)$; exactly,
  $v_3(\mathrm{inc})-2n\in[-17,-7]$ for $n=100..110$.  **[verified]**
* $v_2$ of the same increments is $-2,-4,-6,-4,-8,-4,-6,-4,-10$ at
  $n=50,\dots,800$ — bounded: $\sigma_2=0=v_2(9)$.  **[verified]**
* $B-\xi_3A$ with $\xi_3=\lim_3b_n/a_n$: $v_3=195,392,596,792,995,1191,1389,1591$
  at $n=100,\dots,800$; LS slope $+1.99759$, wide window $+1.99833$,
  $\min v/n=+1.9535$.  **Slope exactly $2$.**  ($v_3(\xi_3)=-1$;
  $\xi_3=\tfrac12\zeta_3(2)$ by `EULER_CRITERION.md` §4.1.)
* **No rational direction works.**  Over all $r=p/q$ with $\lvert p\rvert\le6$,
  $1\le q\le6$, the best LS slope of $v_3$ of the coefficients of $B-rA$ is
  $+0.002890$ on $[100,400]$, falling to $+0.002375$ on $[200,800]$ — a pure
  $\log$ artefact (it is the $+1.17\log_3n$ of $v_3(a_n)$, and essentially every
  $r$ ties, since generically $v_3(b_n-ra_n)=\min(v_3(b_n),v_3(ra_n))$).
  At $p=2$ the best rational $r$ gives $-0.006671$, halving from $-0.013317$;
  and $B-\xi A$ with the archimedean $\xi$ gives $v_2\equiv-20$, constant.
* **Nor at the level of the conditional function.**  Over all $865$ primitive
  triples $(a,b,c)\in\{-6..6\}^3$ for $H=aH_A+bH_B+cH_C$: max LS slope
  $+0.000422\to+0.000214$ ($p=2$) and $+0.002890\to+0.002375$ ($p=3$) under
  window doubling.  **No triple beats slope $0$.**  CDT's own $(1,-3,5)$ is not
  among the near-winners; it inherits $-2.15\log_2n$ from $H_B$ and
  $-0.49\log_3n$ from $H_C$.

### 2.3 The adelic margin on CDT's host

`08_adelic.py` §1.  $\log\lvert\varphi'(0)\rvert=5.081908$ (CDT's realised
contour), $\mathrm{BC}=11.845$, $\tau=16603/3920=4.235459$.

| accounting | $m$ | $\tau$ | $\gamma_2$ | $\gamma_3$ | entry | margin |
|---|---|---|---|---|---|---|
| CDT, $\gamma\equiv0$ (their own numbers) | 14 | $4.2355$ | $0$ | $0$ | $+0.8464$ | $\mathbf{+0.0053}$ |
| **adelic, measured slopes (all $\varsigma=0$)** | 14 | $4.2355$ | $\mathbf 0$ | $\mathbf 0$ | $+0.8464$ | $\mathbf{+0.0053}$ |
| same, at the hard ceiling $\log256$ | 14 | $4.2355$ | $0$ | $0$ | $+1.3097$ | $+6.491$ |

> **The adelic refinement is worth exactly nothing on CDT's own host.**  This is
> the $\lambda_2=1$ case of the identity of §6: an Apéry-perfect host has all of
> its budget on the archimedean side, and CDT already spend all of it.

### 2.4 What the unreachable $3$-adic budget would be worth

Counterfactual: suppose $u$ of the $14$ functions were $\mathbf Q$-rational and
carried the row's resource, $\varsigma_x=2$, i.e. $\varsigma_y=4$ in the
symmetrised coordinate (`ADELIC_HOLONOMY.md` §3, $\varsigma_y\ge2\varsigma_x-2v_p(\lambda_2)$
with $v_3(\lambda_2)=v_3(1)=0$):

| $u$ | 1 | 3 | 5 | **7** | 10 | 14 |
|---|---|---|---|---|---|---|
| $\gamma_3$ | $0$ | $+0.0897$ | $+0.3587$ | $\mathbf{+0.8071}$ | $+1.816$ | $+4.081$ |
| margin | $+0.005$ | $+1.261$ | $+5.028$ | $\mathbf{+11.305}$ | $+25.43$ | $+57.13$ |

So the $3$-adic side of CDT's host is worth $+11.3$ of margin at CDT's own
$7{+}7$ split — twenty times their realised $+0.0053$ — and it is *entirely
inaccessible*, because the only vector realising it is $B-\tfrac12\zeta_3(2)A$.
This is the sharpest way to state what the irrationality of $\zeta_3(2)$ costs
the method, and it is the exact mirror of `ADELIC_HOLONOMY.md` §3's middle row.

---

## 3. Can a two-prime decayer be added to the space?

### 3.1 The three tests, and the answers

For $f=\sum_bc_bx^b$ to join CDT's inventory it must (i) be $\mathbf Q(x)$-linearly
independent of them and holonomic, (ii) satisfy the **denominator type** (6.0.9)
for the *same* $\mathbf b,\mathbf e$ array or an enlarged one, and (iii) be such
that $\varphi^*f$ is meromorphic on $\mathbf D$ for the *same* $\varphi$ — i.e.
the enlarged singular set must still be uniformised by CDT's $\varphi$.

| test | conductor 3 | conductor 6 | conductor 12 |
|---|---|---|---|
| same holonomic module? | **no** | **no** | **no** |
| same singular set? | **no** | **no** | **no** |
| denominator type? | yes, with $\varsigma_3=-3$ ($\mathrm{den}(Q_m)$ is a pure $3$-power) | yes, LCM layer of rate $\nu$ plus $\varsigma_2=-\kappa_2$, $\varsigma_3=-\kappa_3$ | yes (product form only — see §4.2) |

### 3.1a  The exact recurrences and singular sets

`02_gen_seqs.gp`, `02_gen_halves.gp`, `02_fit_recur.gp`, `02_run_all.gp`.  Joint
$(Q,P)$ recurrences $\sum_{i=0}^rc_i(n)u_{n+i}=0$ were fitted over $\mathbf Q$
(mod-$p$ nullity sweep, then exact rational reconstruction), with $\ge190$ terms
used in the fit and $8$ held out.  **Nullity $1$ in every case, so these are *the*
minimal recurrences; zero failures on the held-out indices.**  This closes
`MULTI_PRIME_LATTICE.md` §8 open item 1 for $\alpha\in\{1,2\}$.

| family | $r$ | $\deg c_i$ | $\chi(x)$ | finite nonzero singular points of $\sum u_nx^n$ | $\min\lvert q\rvert$ |
|---|---|---|---|---|---|
| conductor 3, $\alpha=2$ | $2$ | $9$ | $(x-1)(x+1024)$ | $1,\ -\tfrac1{1024}$ | $2^{-10}$ |
| **conductor 6, $\alpha=2$** | $2$ | $15$ | $\mathbf{(x-1)(x+1024)}$ | $1,\ -\tfrac1{1024}$ | $2^{-10}$ |
| conductor 6, $\alpha=1$ | $2$ | $(9,11,13)$ | — (degrees unequal) | **none: entire** | $\infty$ |
| conductor 12, $\alpha=2$ ($Q$) | $1$ | $18$ | $x-1048576$ | $\tfrac1{2^{20}}$ | $2^{-20}$ |
| conductor 12 halves, $\alpha=2$ | $2$ | $15$ | $(x-1)(x+1024)$ | $1,\ -\tfrac1{1024}$ | $2^{-10}$ |
| conductor 12 halves, $\alpha=1$ | $2$ | $(9,11,13)$ | — | **none: entire** | $\infty$ |

**[verified]**  Four facts worth recording separately.

1. **The conductor-6 $\alpha=2$ row has the *same* characteristic polynomial as
   the conductor-3 $\alpha=2$ row, $(x-1)(x+1024)$** — so $\Lambda=1024=2^{10}$
   exactly (the $6.775$ of `MULTI_PRIME_LATTICE.md` §2.2 is a finite-$b$ reading
   of $\log1024=6.9315$), and both generating functions are singular at exactly
   $x=1$ and $x=-2^{-10}$.  Neither is CDT's $\{0,\tfrac19,1,\infty\}$: they share
   $0,1,\infty$ and differ in the fourth point, $-2^{-10}$ against $\tfrac19$.
2. $Q_b$ alone satisfies a **first-order** (hypergeometric) recurrence at
   $\alpha=2$, degree $9$, exactly
   $$\frac{Q_{b+1}}{Q_b}=-\,\frac{(2b+1)(12b+1)^2(12b+5)^2(12b+7)^2(12b+11)^2}
   {72\,(b+1)^3(3b+1)(3b+5)(6b+1)^2(6b+5)^2}\ \longrightarrow\ -1024 ,$$
   the conductor-6 analogue of the conductor-3 closed form of
   `ONE_CLASS_TWO_WORLDS.md` §3 (check: $Q_1/Q_0=(-17787/16)/(135/2)=-148225/9000$,
   the formula at $b=0$).  **[verified]** on $b\le199$.
3. **At $\alpha=1$ there is no finite singularity at all.**  The degrees
   $(\deg c_0,\deg c_1,\deg c_2)=(9,11,13)$ are unequal; the Newton polygon of
   $(i,\deg c_i)$ has both edges of slope $2$, so $u_{b+1}/u_b\sim\rho b^{-2}$
   — factorial decay — and $\sum u_bx^b$ is **entire**.  The same holds for both
   conductor-12 half-rows at $\alpha=1$.  For those members the contour cost of
   §3.3 is therefore exactly **zero**, and their entire cost is the denominators.
   (The $\rho$'s are the roots of $185752092672x^2-6315571150848x+185752092672$,
   $\lvert\rho\rvert=0.02944,\ 33.97$.)
4. The conductor-12 product form inherits its root from its halves:
   $\Lambda=1024^2=2^{20}$ exactly at $\alpha=2$, confirming
   $Q_b=160\,s^{(1)}_bs^{(5)}_b$ structurally.

The **geometric** part of the denominator cannot be absorbed into an LCM layer:
$v_p([1,\dots,Bn])=O(\log n)$, so a $p^{\kappa n}$ denominator is a genuine
negative slope $\varsigma_p=-\kappa$ and the theorem of `ADELIC_HOLONOMY.md` §2.3
charges it as such.  This is the whole story:

> **$\kappa_p>0$ — the very thing that makes a decayer a decayer, and that makes
> it a *two-prime* decayer — is negative overconvergence.**

### 3.2 The scaling invariant

The identification of the decayer's own variable with the host's coordinate is
free up to a scalar $x\rightsquigarrow x/N$.  Under it the singularity moves out
by $N$ and every $\varsigma_p$ shifts by $-v_p(N)$, so for integer $N$ the
quantity

$$\rho_{\mathrm{int}}\ :=\ \log\Lambda\ +\ \nu\ +\ \sum_p\kappa_p\log p$$

($\Lambda=\lim\lvert Q_b\rvert^{1/b}$, $\nu$ the prime-to-$\{2,3\}$ denominator
rate) is **invariant**: it is the growth rate of the integrally normalised
coefficient sequence.  It is the single number that measures how expensive the
function is.  **[proved]** (this is the scale-covariance check of
`ADELIC_HOLONOMY.md` §2.4(2) applied to one function, and it is
`CDT_NONCONGRUENCE.md` §2's "$\lambda^n$ costs exactly $\log\lambda$" for a
single member of the inventory.)

| family | $\alpha$ | $\log\Lambda$ | $\kappa_2$ | $\kappa_3$ | $\nu$ | $\rho_{\mathrm{int}}$ |
|---|---|---|---|---|---|---|
| conductor 3 | $2$ | $+6.932$ | $0$ | $3.00$ | $0$ | $10.23$ |
| conductor 6 | $1/4$ | $-10.198$ | $4.41$ | $0$ | $9.84$ | $\mathbf{2.70}$ |
| conductor 6 | $1/2$ | $-7.371$ | $4.95$ | $0$ | $8.34$ | $4.40$ |
| conductor 6 | $1$ | $-2.062$ | $5.83$ | $0.83$ | $5.48$ | $8.36$ |
| conductor 6 | $2$ | $+6.9315$ (exact, $=\log2^{10}$) | $7.82$ | $2.85$ | $0.66$ | $16.14$ |
| conductor 12 | $1/4$ | $-23.541$ | $2.00$ | $0$ | $29.50$ | $7.35$ |
| conductor 12 | $1/2$ | $-18.711$ | $4.00$ | $0$ | $27.28$ | $11.34$ |
| conductor 12 | $1$ | $-8.214$ | $8.00$ | $1.73$ | $20.10$ | $19.34$ |
| conductor 12 | $2$ | $+13.8629$ (exact, $=\log2^{20}$) | $16.00$ | $5.73$ | $1.90$ | $33.15$ |

(conductor-6 rates from `MULTI_PRIME_LATTICE.md` §2.2; conductor-12 rates
measured here, §4.2.  $\rho_{\mathrm{int}}>0$ in every case: **there is no free
function anywhere in these families.**)

*$\rho_{\mathrm{int}}$ is the stable quantity.*  `13_rho_int.gp` measures
$\tfrac1b\log\lvert Q_b\cdot\mathrm{den}(Q_b)\rvert$ directly:

| $\alpha$ | $1/4$ | $1/2$ | $1$ | $2$ |
|---|---|---|---|---|
| conductor 6, $b=20/40/60/80$ | $2.95,2.80,2.85,2.70$ | $4.31,4.24,4.65,4.46$ | $7.82,8.36,8.45,8.47$ | $15.10,16.03,15.98,16.04$ |
| conductor 12, $b=20/40/60/80$ | $6.59,6.58,6.93,7.16$ | $10.07,10.71,10.75,10.75$ | $18.41,19.79,19.39,19.80$ | $32.17,33.04,32.82,33.29$ |

$\rho_{\mathrm{int}}$ is flat in $b$ while $\log\Lambda$ and $\nu$ separately
drift by several units (at $\alpha=1$ the coefficients decay factorially, so
$\log\Lambda\to-\infty$ and $\nu\to+\infty$ at the same $\log b$ rate and the
sum is finite).  **[verified]**  This is why the accounting of §§3.3–3.4 is
robust even though its inputs are not.

*Caveat on $\nu$.*  `12_c6_primes.gp` re-measures the conductor-6 family and finds
$\nu=5.63,7.66,9.08$ at $b=20,40,60$ ($\alpha=1/4$): $\nu$ drifts upward at
roughly a $\log b$ rate in both families, so all $\nu$'s here (and in
`MULTI_PRIME_LATTICE.md` §2.2) are finite-$b$ readings, not limits.  The prime
*support* is bounded and clean in both: $\le3b-1$ for $\mathrm{den}(Q_b)$ and
$\le6b-1$ for $\mathrm{den}(P_b)$ at conductor 6, $\le6b+1$ and $\le12b-1$ at
conductor 12.  **[verified]**

### 3.3 Contour cost of an extra singular point

If the added function branches at $y=q$ then $\varphi$ must factor through the
universal cover of $\mathbf P^1$ minus the enlarged singular set, which is
contained in $\mathbf P^1\setminus\{0,q,\infty\}$; the $\lambda$-uniformisation
of that gives, unconditionally,

$$\lvert\varphi'(0)\rvert\ \le\ 16\lvert q\rvert\qquad(\text{or } 256\lvert q\rvert
\ \text{if the } \mathbf Z/2 \text{ descent is available with } q \text{ outer}).$$

The bound is not vacuous: at CDT's own configuration $q=s=1$ it reads
$\log256=5.5452$ and CDT realise $5.0819$.

| family | $\lvert q\rvert=1/\Lambda$ | $\log(16\lvert q\rvert)$ | $\log(256\lvert q\rvert)$ | shortfall vs $\tau=4.2355$ |
|---|---|---|---|---|
| conductor 3, $\alpha=2$ | $9.77\cdot10^{-4}$ | $-4.159$ | $-1.386$ | $\mathbf{-5.62}$ |
| conductor 6, $\alpha=2$ | $2^{-10}=9.766\cdot10^{-4}$ | $-4.159$ | $-1.386$ | $-5.62$ |
| conductor 12, $\alpha=2$ | $2^{-20}=9.537\cdot10^{-7}$ | $-11.090$ | $-8.318$ | $-12.55$ |
| conductor 6, $\alpha=1$ | $\infty$ (entire) | — | — | none |
| conductor 6, $\alpha=1/4$ | $\ge2.68\cdot10^{4}$ | $+12.971$ | $+15.743$ | none |
| conductor 12, $\alpha=1/4$ | $\ge1.67\cdot10^{10}$ | $+26.314$ | $+29.086$ | none |

So the families split cleanly, and the split is governed by $\rho_{\mathrm{int}}$:
**large $\alpha$ pays in the contour** (a singularity within $10^{-3}$ of the
origin), **small $\alpha$ pays in the denominators** ($\nu$ up to $29.5$).
Neither end is affordable.  In the accounting of §3.4 the cheaper of the two is
always taken (an integral rescaling by $N=\lceil\Lambda s\rceil$ converts contour
cost into slope cost at the invariant rate of §3.2).

### 3.4 The numbers

`08_adelic.py` §3.  "$\mathrm{entryC}$" is at the host's hard ceiling,
"$\mathrm{marginC}$" the corresponding margin; $\gamma_2,\gamma_3$ are the
*total* $\gamma$'s of the enlarged $m=15$ (resp. $16$) inventory.

**On CDT's level-6 host** (baseline $m=14$: $\tau=4.2355$, $\gamma=0$,
entryC $=+1.3097$, marginC $=+6.491$):

| added function | $\tau$ | $\gamma_2$ | $\gamma_3$ | entryC | marginC |
|---|---|---|---|---|---|
| conductor 3, $\alpha=2$ | $4.2269$ | $-1.2631$ | $-0.6006$ | $-0.5454$ | $-20.03$ |
| **conductor 6, $\alpha=1/4$** | $5.4952$ | $-0.5570$ | $0$ | $\mathbf{-0.5070}$ | $-19.45$ |
| conductor 6, $\alpha=1/2$ | $5.3018$ | $-0.6252$ | $0$ | $-0.3819$ | $-17.57$ |
| conductor 6, $\alpha=1$ | $4.9332$ | $-0.7364$ | $-0.1662$ | $-0.2906$ | $-16.20$ |
| conductor 6, $\alpha=2$ | $4.3120$ | $-2.2508$ | $-0.5705$ | $-1.5881$ | $-35.67$ |
| conductor 12, $\alpha=1/4$ | $8.0292$ | $-0.2526$ | $0$ | $-2.7367$ | $-52.90$ |
| conductor 12, $\alpha=1$ | $6.8177$ | $-1.0105$ | $-0.3469$ | $-2.6299$ | $-51.29$ |
| conductor 12, $\alpha=2$ | $4.4713$ | $-4.5470$ | $-1.1477$ | $-4.6208$ | $-81.16$ |

Every row fails the **entry condition** — the theorem does not even apply.  The
best case, conductor 6 at $\alpha=1/4$, costs $1.82$ of entry against a headroom
of $1.31$.

**Break-even** (`08_adelic.py` §5).  On CDT's host a fifteenth function of
$p$-adic slope $0$ pays for itself iff its LCM denominator rate is
$\nu\le0.744$; with $\nu=0$ it pays for itself iff $\varsigma_2\ge-0.759$, i.e.
iff its geometric denominator rate is $\kappa_2\le0.759$.  The cheapest decayer
in the corpus has $\kappa_2=4.41$.

$$\boxed{\ \gamma_2+\gamma_3\ <\ 0\ \text{ for every member of every family}\ }$$

On CDT's host (where the baseline $\gamma$ is $0$) the contribution is between
$-0.253$ (conductor 12, $\alpha=1/4$) and $-5.695$ (conductor 12, $\alpha=2$).
On the level-16 Catalan host, where the baseline is $\gamma_2=+0.2619$, the
*change* in $\gamma_2+\gamma_3$ runs from $-0.270$ to $-5.145$.

---

## 4. Part (b): the conductor-12 Catalan row

### 4.1 Construction

> **Theorem $1''$ (conductor 12).**  **[proved]**  Fix $b\ge0$; let $R=N/D_b$ with
> $$D_b(t)=\prod_{j=0}^{b}\ \prod_{r\in\{1,5,7,11\}}\Bigl(t+j+\tfrac r{12}\Bigr)^{2},
> \qquad\deg D_b=8b+8,$$
> $\deg N\le8b+6$ and $N(-t-b-1)=-N(t)$.  The involution $\sigma(t)=-t-b-1$ sends
> the class $r$ to the class $12-r$, so $\sum A^{(12-r)}=-\sum A^{(r)}$ and
> $\sum B^{(12-r)}=+\sum B^{(r)}$; $\deg N\le\deg D-2$ forces $\sum_iB_i=0$.
> Impose the two further linear conditions
> $$\textstyle(\mathrm{C}1)\ \sum B^{(1/12)}=0,\qquad
> (\mathrm{C}2)\ \sum A^{(1/12)}=\sum A^{(5/12)} .$$
> Then
> $$\sum_{t\ge0}R(t)=Q\,G-P,\qquad Q=160\sum_jA_j^{(1/12)},\qquad P\in\mathbf Q,$$
> $G=L(2,\chi_{-4})$ Catalan's constant.  No $\zeta(2)$, no $L(2,\chi_{-3})$, no
> $\pi$, no logarithm.

*Proof.*  Exactly the proof of Theorem $1'$ of `MULTI_PRIME_LATTICE.md` §2.1 with
the conductor-12 evaluations
$$\zeta(2,\tfrac1{12})-\zeta(2,\tfrac{11}{12})=72\bigl(L_{12}(2,\chi_{-4})+L_{12}(2,\chi_{-3})\bigr),\qquad
\zeta(2,\tfrac5{12})-\zeta(2,\tfrac7{12})=72\bigl(L_{12}(2,\chi_{-4})-L_{12}(2,\chi_{-3})\bigr),$$
where $L_{12}(2,\chi_{-4})=\bigl(1-\chi_{-4}(3)3^{-2}\bigr)G=\tfrac{10}9G$ and
$L_{12}(2,\chi_{-3})=\tfrac54L(2,\chi_{-3})$; $(\mathrm C2)$ kills the
$\chi_{-3}$ component and $144\cdot\tfrac{10}9=160$.  $(\mathrm C1)$ is needed
because, unlike conductor 6, the $\psi$-terms do **not** cancel automatically:
the surviving combination is
$$\psi(\tfrac1{12})+\psi(\tfrac{11}{12})-\psi(\tfrac5{12})-\psi(\tfrac7{12})
=4\sqrt3\,\log(2-\sqrt3)\neq0 . \qquad\square$$

**Note the arithmetic significance of $160=144\cdot\tfrac{10}9$**: the Euler
factor $\mathcal E_3(2)=\tfrac{10}9$ of `EULER_CRITERION.md` is built into the
normalisation, which is exactly why this row has a $3$-adic avatar while
Zudilin's conductor-4 row does not — the precise mirror of the conductor-6
constant $45=36\cdot\tfrac54$.

**Two implementations, and only one of them is any good.**

* **(I) four classes, free even quartic.**  $N=(2t+b+1)\prod_{j=1}^a(t-j+1)(t+b+j)\cdot(u^4+c_1u^2+c_0)$
  with $u=t+\tfrac{b+1}2$; $(\mathrm C1),(\mathrm C2)$ are affine-linear in
  $(c_1,c_0)$ and are solved exactly (`c12row` in `lib12.gp`).  The identity
  $\sum_{t\ge0}R=QG-P$ checks against $\sum_{t\le40000}R(t)$ to
  $2.9\cdot10^{-75}$ at $(a,b)=(1,2)$ and $9.1\cdot10^{-82}$ at $(1,3)$.
  **[verified]**  This form has $\sigma_2=16$, $\sigma_3=6$ and the same
  $p$-adic limits (verified to $v_2=753$, $v_3=282$ at $b=48$) — but the linear
  solve destroys the arithmetic: at $\alpha=1$, $b=24$, $\mathrm{den}(Q)$ has
  $35$ distinct prime factors, the largest a **179-digit (pseudo)prime**.  It is
  therefore **not of CDT denominator type at all** (no $[1,\dots,Bn]$ layer can
  contain it) and useless for lattices.  Recorded as a negative.
* **(II) product form** (`c12rowB`).  Split the four classes into the two
  $\sigma$-pairs $\{1,11\}$ and $\{5,7\}$ and run a conductor-6-style well-poised
  row on each: $(\mathrm C1)$ then holds on each half automatically.  Tie them by
  cross-scaling — no division:
  $$\boxed{\ Q_b=160\,s^{(1)}_b\,s^{(5)}_b,\qquad P_b=s^{(5)}_bP^{(1)}_b+s^{(1)}_bP^{(5)}_b\ }$$
  with $s^{(r)}_b=\sum_jA_j^{(r/12)}$ of the $r$-half.  All denominators are
  products of two partial-fraction denominators, hence genuinely of LCM type:
  `11_chi12_primes.gp` finds the largest prime of $\mathrm{den}(Q_b)$ to be
  $113,241,359$ at $b=20,40,60$ (i.e. $\le6b+1$) and of $\mathrm{den}(P_b)$
  $239,479,719$ ($\le12b-1$), against the partial-fraction ceiling $12b+11$.
  **[verified]**  **This is the row.**

### 4.2 Integrality structure (product form, $\alpha=a/b$ equal on both halves)

`06_chi12_full.gp`, `10_chi12_small_alpha.gp`.  Wide-window rates per unit of $b$
($b:30\to60$, resp. $32\to64$ for $\alpha\le3/4$).

| $\alpha$ | $\log\Lambda$ | $\log\lambda$ | $\kappa_2$ | $\kappa_3$ | $\nu$ | $\eta$ |
|---|---|---|---|---|---|---|
| $1/4$ | $-23.541$ | $-26.077$ | $2.000$ | $0.000$ | $29.50$ | $55.98$ |
| $1/2$ | $-18.711$ | $-19.926$ | $4.000$ | $0.000$ | $27.28$ | $54.37$ |
| $3/4$ | $-13.677$ | $-14.007$ | $6.000$ | $1.063$ | $23.79$ | $52.72$ |
| $1$ | $-8.214$ | $-8.463$ | $8.000$ | $1.733$ | $20.10$ | $51.47$ |
| $3/2$ | $+2.526$ | $+2.561$ | $12.000$ | $3.733$ | $11.22$ | $46.55$ |
| $2$ | $+13.748$ | $+13.494$ | $16.000$ | $5.733$ | $1.90$ | $38.45$ |

$$\boxed{\ \kappa_2=8\alpha\ \text{exactly (all six values)},\qquad
\kappa_3\approx4\alpha-2.27\ \text{for}\ \alpha\ge3/4,\ =0\ \text{below}\ }$$
**[measured]**; $\eta=\tfrac1b\log\operatorname{lcm}(\mathrm{den}Q_b,\mathrm{den}P_b)$
and $\nu$ the prime-to-$\{2,3\}$ part of $\mathrm{den}Q_b$.  $\nu$ and $\eta$
still drift upward with $b$; the values above are the $b\le64$ readings and
should be read as lower estimates.  The *support* is settled — every prime of
$\mathrm{den}(Q_b)$ is $\le6b+1$ and of $\mathrm{den}(P_b)$ is $\le12b-1$
(`11_chi12_primes.gp`, $b=20,40,60$) — so the denominators are of LCM type and
the function does satisfy CDT's (6.0.9) with $\sum_jb_{i,j}\approx\nu$; what is
**[open]** is the sharp multiplicity, i.e. the smallest $k$ with
$\mathrm{den}(Q_b)\mid 12^\ast D_{6b+1}^{\,k}$ ($k=2$ and $k=4$ both fail
already at $b=20$; $\nu\approx27$ against $\log D_{6b+1}/b\to6$ suggests
$k\approx5$).  This is the conductor-12 form of
`MULTI_PRIME_LATTICE.md` §8 open item 2.

The row **decays** raw ($\log\lambda<0$ for $\alpha\le1$) but its integralised
linear forms grow at rate $\eta+\log\lambda\approx30$–$43$: there is no
Diophantine content whatever, and none is claimed.

### 4.3 The two slopes

From the alignment valuations of §4.4 (wide window $b:30\to60$):

$$\sigma_2=\frac{463-223}{30}=\mathbf{8.00},\qquad
\sigma_3=\frac{176-89}{30}=\mathbf{2.90}\to3 ,$$

uniformly over $\alpha\in\{1/2,1,3/2,2\}$.  **[measured]**

> **The conductor-12 Catalan row has exactly the conductor-6 profile,
> $\sigma_2=8$, $\sigma_3=3$ per unit of $b$** — the two constructions are
> arithmetic mirrors of one another across $\chi_{-3}\leftrightarrow\chi_{-4}$.

(The four-class form (I) has $\sigma_2=16$, $\sigma_3=6$: its index $b$ is worth
twice as much because $Q_b$ is not a product of two halves.)

### 4.4 The $p$-adic limits against Theorem F

**Prediction.**  $\psi=\chi_{-4}$, $\varphi=\mathbf 1$, $w=1$, $r_\infty=1$
(the row's archimedean limit is $G$ itself).  By
`EULER_CRITERION.md`'s value formula $\xi_p=r_p\Lambda_p$,
$\Lambda_p=L_p(w+1,\psi\omega^{-w})$, $r_p=r_\infty/\mathcal E_p(w+1)$:

| $p$ | $\omega$ | $\psi\omega^{-1}$ | $\mathcal E_p(2)=1-\chi_{-4}(p)p^{-2}$ | prediction |
|---|---|---|---|---|
| $2$ | $\chi_{-4}$ (Washington's mod-4 normalisation) | $\mathbf 1$ | $1$ ($\chi_{-4}(2)=0$) | $\xi_2=\zeta_2(2)=L_2(2,\mathbf1)$ |
| $3$ | $\chi_{-3}$ | $\chi_{-4}\chi_{-3}=\chi_{12}$ | $1+\tfrac19=\tfrac{10}9$ | $\xi_3=\tfrac9{10}L_3(2,\chi_{12})$ |

$\kappa_p=\tfrac12L_p(2,\psi\omega^{-1})$ in the notation of the criterion; both
values are computed by `lattice/euler_criterion/lp.gp` (Washington Thm 5.11,
exact integer binomials at integer $s$).

**Measured** (`06_chi12_full.gp`, product form, $\xi_b=P_b/Q_b$).  "Cauchy" is
$v_p(\xi_b-\xi_{b-1})$, the precision the row itself carries at that $b$:

| $\alpha$ | $b$ | $v_2\bigl(\xi_b-\zeta_2(2)\bigr)$ | Cauchy$_2$ | $v_3\bigl(\xi_b-\tfrac9{10}L_3(2,\chi_{12})\bigr)$ | Cauchy$_3$ |
|---|---|---|---|---|---|
| $1/2$ | $30$ | $223$ | $211$ | $89$ | $82$ |
| $1/2$ | $60$ | $463$ | $446$ | $176$ | $171$ |
| $1$ | $30$ | $223$ | $215$ | $90$ | $84$ |
| $1$ | $45$ | $343$ | $339$ | $130$ | $120$ |
| $1$ | $60$ | $\mathbf{463}$ | $451$ | $\mathbf{174}$ | $168$ |
| $3/2$ | $60$ | $463$ | $446$ | $176$ | $171$ |
| $2$ | $60$ | $467$ | $456$ | $176$ | $171$ |

Agreement **exceeds the row's own Cauchy precision at every entry**: the
prediction is confirmed to the full available accuracy, at both primes, for four
independent members of the family.  **[verified]**

**Controls** ($\alpha=1$, $b=60$; twelve rational scalars): $v_2(\xi-c\,\zeta_2(2))\le-1$
for every $c\ne1$ and $v_3(\xi-c\,L_3(2,\chi_{12}))\le1$ for every
$c\ne\tfrac9{10}$.  The four-class form (I) reproduces the same two limits
independently, to $v_2=753$, $v_3=282$ at $b=48$, with the same controls.

**Consequences.**

1. $\xi_2=\zeta_2(2)$ is *exactly* Zudilin's Catalan row's $2$-adic limit
   (`EULER_CRITERION.md` §4.1: $\xi_2^{\mathrm{Zud}}=\zeta_2(2)$) and twice
   Zagier $\mathbf E$'s ($\tfrac12\zeta_2(2)$): the conductor-12 row aligns
   $1{:}1$ with Zudilin and $2{:}1$ with $\mathbf E$ at $p=2$.
2. $\xi_3=\tfrac9{10}L_3(2,\chi_{12})$ is a **new constant in the corpus**.  The
   $\chi_{12}$ $p$-adic $L$-value had appeared only at $p=2$
   (row $\mathbf F$: $\tfrac12L_2(2,\chi_{12})$; conductor 6:
   $\tfrac45L_2(2,\chi_{12})$).  Its appearance at $p=3$ here is the exact mirror
   of its appearance at $p=2$ there, and it is forced by the same mechanism:
   $\chi_{-3}\omega_2^{-1}=\chi_{-4}\omega_3^{-1}=\chi_{12}$.
3. $L_3(2,\chi_{12})$ is a $3$-adic **unit**
   ($=1+2\cdot3+3^3+2\cdot3^5+\cdots$), so $\xi_3\neq0$: the $3$-adic Catalan
   slot is genuinely occupied, not vacuously so.

### 4.5 There is no partner

`09_catalan_3adic.gp`.  On Zagier $\mathbf E$ $(12,4,32)$ and on Zudilin's row,
$v_3$ of the Apéry increments is $0$ or $-2$ out to $n=400$ (resp. $m=300$),
while $v_2/n\to4.92$ resp. $7.90$.  **$\sigma_3=0$ on every Catalan-class
modular row**, as `MULTI_PRIME_LATTICE.md` §7 predicted from $c=2^\ast$.  So the
conductor-12 row has a $3$-adic avatar and nothing to bridge it against; open
item 7 of `MULTI_PRIME_LATTICE.md` §8 is answered — *the construction exists,
the partner does not*.

---

## 5. The conductor-12 function in the level-8 and level-16 adelic bounds

`08_adelic.py` §3.  Level-8 host: $s=\tfrac14$, ceiling $\log64=4.158883$,
$\mathrm{BC}=10.458706$, baseline $m=14$ ($7$ pure of $\varsigma_y=2$, $7$
conditional of $\varsigma_y=0$).  Level-16 host: same geometry
(`CATALAN_TWO_CLASSES.md` §5), baseline $m=15$ with the doubly-small function of
$\varsigma_y=1$.

| host | inventory | $m$ | $\tau$ | $\gamma_2$ | $\gamma_3$ | entryC | **marginC** | deficit $=\dfrac{\text{margin}}{m-1}$ |
|---|---|---|---|---|---|---|---|---|
| level 8 | CDT's $14$ (baseline) | $14$ | $4.2355$ | $+0.2546$ | $0$ | $+0.1780$ | $\mathbf{-7.97}$ | $-0.6128$ |
| level 8 | $+$ conductor 12, $\alpha=1$ | $15$ | $6.8177$ | $-0.7886$ | $-0.3469$ | $-3.7944$ | $\mathbf{-67.37}$ | $-4.812$ |
| level 8 | $+$ conductor 12, $\alpha=1/2$ | $15$ | $7.7434$ | $-0.2834$ | $0$ | $-3.8679$ | $-68.48$ | $-4.891$ |
| level 8 | $+$ conductor 12, $\alpha=2$ | $15$ | $4.4713$ | $-4.0726$ | $-1.1477$ | $-5.5327$ | $-93.45$ | $-6.675$ |
| level 16 | $14+$ doubly-small (baseline) | $15$ | $4.2269$ | $+0.2619$ | $0$ | $+0.1938$ | $\mathbf{-7.55}$ | $-0.5394$ |
| level 16 | $+$ conductor 12, $\alpha=1$ | $16$ | $6.6527$ | $-0.7229$ | $-0.3272$ | $-3.5439$ | $\mathbf{-67.16}$ | $-4.477$ |
| level 16 | $+$ conductor 12, $\alpha=1/2$ | $16$ | $7.5224$ | $-0.2464$ | $0$ | $-3.6099$ | $-68.22$ | $-4.548$ |

At CDT's realised contour every entry falls a further $0.4633$ and every margin a
further $\approx m\cdot0.4633$.

> **Answer to the task's question.**  Level 8: $-7.97\to-67.4$.  Level 16:
> $-7.55\to-67.2$.  In `CDT_NONCONGRUENCE.md`'s per-function units,
> $-0.613\to-4.81$ and $-0.539\to-4.48$.  The conductor-12 function is worth
> $-4.2$ per function, i.e. it costs about seven times the entire remaining
> deficit.  **It must not be added.**

For calibration, the same script reproduces `CATALAN_TWO_CLASSES.md` §6 exactly
($m=15$, $\tau=4.2269$, $\gamma_2=+0.2619$, entry $+0.1938$, margin $-7.551$)
and its "$\approx14$ doubly-small functions would be needed": with
$\varsigma_2=1$ the answer is **exactly $14$**; with $\varsigma_2=2$, $7$; with
$\varsigma_2=4$, $6$; with $\varsigma_2=8$, $5$.  The break-even for a fifteenth
function on the level-8 host is $\varsigma_2\ge+0.309$ at $\nu=0$: **an extra
function of CDT's own denominator type and no overconvergence makes the Catalan
margin worse**, because the $(u-1)^2/m^2$ weight on the pure orbit's $\gamma_2$
falls faster than $m$ grows.

---

## 6. What a second prime is actually worth: the pure-module identity

The one place where a genuinely multi-prime *admissible* gain exists is the pure
polylogarithm orbit $\mathrm{Li}_j(x/s)=\sum(\lambda_2x)^n/n^j$, whose slope is
$v_p(\lambda_2)$ at every $p$.  Hence, for $\lambda_2\in\mathbf Z\setminus\{0\}$,

$$\boxed{\ \underbrace{\log\frac{256}{\lvert\lambda_2\rvert}}_{\text{archimedean ceiling}}
\ +\ \underbrace{\sum_p\varsigma_p\log p}_{\text{all primes at once}}
\ =\ \log 256\quad\text{on every host.}\ }\qquad\textbf{[proved]}$$

This is `ADELIC_HOLONOMY.md` §5's ceiling made exact and prime-complete: there
the right-hand side was $\log(256\lambda_1)$ because $\sigma_p=v_p(c)=v_p(\lambda_1\lambda_2)$
was the *row's* slope; for the *pure module* the slope is $v_p(\lambda_2)$ and
the identity closes with no slack at all.

But the two halves are not weighted equally in the entry condition.  The
archimedean half enters at weight $1$; the $p$-adic half at
$(u-1)^2/m^2=36/196=0.1837$ for CDT's $7{+}7$ split (at best $1-1/m=0.93$, for a
uniform profile).  Therefore:

$$\text{entry}=\log256-\tau-\Bigl(1-\tfrac{(u-1)^2}{m^2}\Bigr)\log\lvert\lambda_2\rvert
=\log256-\tau-0.8163\log\lvert\lambda_2\rvert .$$

| host | $\lambda_2$ | ceiling | $\sum_pv_p(\lambda_2)\log p$ | total | $\gamma$ ($u{=}7,m{=}14$) | entry | margin |
|---|---|---|---|---|---|---|---|
| Apéry / CDT / $X_1(5)$ | $1$ | $5.5452$ | $0$ | $5.5452$ | $0$ | $+1.3097$ | $+6.491$ |
| Catalan level 8/16, Domb | $4$ | $4.1589$ | $1.3863$ | $5.5452$ | $+0.2546$ | $+0.1780$ | $-7.966$ |
| Zagier $\mathbf F$, $\sqrt{s_{10}}$ | $8$ | $3.4657$ | $2.0794$ | $5.5452$ | $+0.3819$ | $-0.3878$ | $-15.195$ |
| Cooper $s_{18}$ | $12$ | $3.0603$ | $2.4849$ ($=2\log2+\log3$) | $5.5452$ | $+0.4564$ | $-0.7188$ | $-19.423$ |
| $\sqrt{s_{18}}$ | $24$ | $2.3671$ | $3.1781$ ($=3\log2+\log3$) | $5.5452$ | $+0.5837$ | $-1.2846$ | $-26.652$ |

The last two rows are the only genuinely **two-prime** $\gamma$'s in the whole
corpus, and they are the *worst* hosts, not the best.  This is the clean answer
to "what is the second prime worth in the arithmetic holonomy bound":

> **Nothing.  Every prime factor of $\lambda_2$ buys $v_p(\lambda_2)\log p$ of
> $p$-adic gain at weight $0.18$ and costs the same $v_p(\lambda_2)\log p$ of
> archimedean ceiling at weight $1$; the net is $-0.82\,v_p(\lambda_2)\log p$ per
> prime.**  The Apéry-perfect hosts $\lambda_2=\pm1$ are optimal, and CDT are
> already sitting on one.

(Compare `MULTI_PRIME_LATTICE.md` §4.1's cancellation identity for the *lattice*
side: there the second prime on an integral engine is exactly self-cancelling,
`net`$(r)=r\log\Lambda_{\rm eng}-\sum_p\max(0,\sigma_p^{\rm eng}r-\sigma_p^{\rm dec})\log p$.
The holonomy analogue is worse than self-cancelling — it is strictly negative —
because of the $(u-1)^2/m^2$ discount.)

---

## 7. Ledger

**[proved]**
* Theorem $1''$ (§4.1) and its two implementations; the identity
  $\sum_{t\ge0}R=QG-P$ with $Q=160\sum A^{(1/12)}$; the necessity of
  $(\mathrm C1)$ via $\psi(\tfrac1{12})+\psi(\tfrac{11}{12})-\psi(\tfrac5{12})-\psi(\tfrac7{12})=4\sqrt3\log(2-\sqrt3)$.
* The scaling invariant $\rho_{\mathrm{int}}=\log\Lambda+\nu+\sum_p\kappa_p\log p$ (§3.2).
* The pure-module multi-prime identity $\log(256/\lvert\lambda_2\rvert)+\sum_pv_p(\lambda_2)\log p=\log256$ (§6).
* The principle of §1: geometric denominators are negative overconvergence, so
  $\kappa_p>0\Rightarrow\gamma_p<0$.

**[verified]** (exact arithmetic, stated ranges)
* All fourteen CDT functions rebuilt exactly over $\mathbf Q$ to $y$-order $830$,
  cross-checked mod $2^{61}-1$ against `indep_check2.py`; slopes zero at $p=2,3$.
* $\sigma_3=2$, $\sigma_2=0$ on row $\mathbf C$ to $n=800$; slope of $B-\xi_3A$
  equal to $2$ to five digits; no rational direction of height $\le6$ with
  positive slope, at either prime, at either level ($B-rA$ and $aH_A+bH_B+cH_C$).
* The conductor-12 identity to $10^{-75}$; its two $p$-adic limits to
  $v_2=463,753$ and $v_3=176,282$ with twelve-scalar controls.
* $\sigma_3=0$ on Zagier $\mathbf E$ ($n\le400$) and Zudilin's row ($m\le300$).
* The 179-digit prime in $\mathrm{den}(Q_{24})$ of implementation (I).
* The minimal joint $(Q,P)$ recurrences of §3.1a: nullity $1$ in the mod-$p$
  sweep, exact rational reconstruction, and zero failures on $8$ held-out indices
  in every case; $\chi=(x-1)(x+1024)$ for conductor 6 at $\alpha=2$ and for both
  conductor-12 halves at $\alpha=2$; $\chi=x-2^{20}$ for the conductor-12 product
  form; the first-order ratio formula for the conductor-6 $Q_b$; the entire-ness
  at $\alpha=1$ (Newton polygon of slope $2$).
* The stability of $\rho_{\mathrm{int}}$ in $b$ (§3.2).
* Every entry of every margin table (all of `08_adelic.py`, which reproduces
  CDT's $+0.0053$, `ADELIC_HOLONOMY.md`'s $-7.966$ and
  `CATALAN_TWO_CLASSES.md`'s $-7.551$ exactly).

**[measured]**
* $\kappa_2=8\alpha$, $\kappa_3$, $\nu$, $\eta$, $\log\Lambda$, $\log\lambda$ for
  the conductor-12 family ($b\le64$); $\sigma_2=8$, $\sigma_3=3$.  These are
  finite-$b$ readings; $\nu$ and $\eta$ still drift upward at a $\log b$ rate.
* The conductor-6 rates quoted in §3.2–3.4 are `MULTI_PRIME_LATTICE.md` §2.2's,
  not re-measured here.

**[open]**
1. **Recurrences at $\alpha\notin\{1,2\}$.**  §3.1a settles $\alpha=1$ and
   $\alpha=2$ exactly for conductor 6 and for the conductor-12 halves
   (`MULTI_PRIME_LATTICE.md` §8 open item 1, now closed at those two values).
   The intermediate members ($\alpha=1/4,1/2,3/2,\dots$) were not fitted; §3.3
   treats them by the measured $\Lambda$, which for $\alpha<1$ is a lower bound
   on the distance to the nearest singularity (and at $\alpha=1$ the truth is
   "entire", so the treatment is conservative).
2. **The sharp denominator exponent** for the conductor-12 product form.  The
   support is settled ($\le6b+1$ for $Q$, $\le12b-1$ for $P$), so the row *is*
   of CDT denominator type; the smallest $k$ with
   $\mathrm{den}(Q_b)\mid12^\ast D_{6b+1}^{\,k}$ is not known ($k=2,4$ fail at
   $b=20$).  Until it is, the $\nu$ used in §§3–5 is a measurement, not a bound,
   and §5's margins could move by a few units — in either direction, but they are
   $-60$ deep.
3. The exact values of $\sigma_2,\sigma_3$ as theorems ($8$ and $3$), and of
   $\kappa_2=8\alpha$.
4. Whether a *doubly-small* conductor-12-like object exists — a
   $\mathbf Q$-rational function of positive slope at **two** primes.  Nothing in
   this note exhibits one, and §6 shows that the pure module can never supply one
   with positive net value.  The only known positive-slope rational unconditional
   function on any Catalan host remains `CATALAN_TWO_CLASSES.md`'s single
   $\varsigma_y=1$ class at level 16.
5. `ADELIC_HOLONOMY.md` §6's items 1, 2, 4, 5 (softness of the theorem, number
   fields, functional independence, the contour) are untouched here.

**No irrationality claim.**  Nothing above proves anything about Catalan's
constant, $L(2,\chi_{-3})$, $\zeta_2(2)$, $\zeta_3(2)$ or $L_3(2,\chi_{12})$.
The conductor-12 row's integralised linear forms **grow** at rate $\approx30$–$43$
per unit of $b$; it is a decayer only in the raw normalisation and carries no
Diophantine content.

---

## 8. Reproduction

```sh
cd /home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy
mkdir -p out

# (a) CDT's 14 functions, exact over Q; slopes at p=2 and p=3
timeout 3000 gp -q 01_check.gp     </dev/null > out/01_check.log 2>&1
timeout 3000 gp -q 01_main.gp      </dev/null > out/01_main.log 2>&1      # y-order 412
timeout 6000 gp -q 01_main_big.gp  </dev/null > out/01_main_big.log 2>&1  # y-order 830
timeout 3000 gp -q 01_host.gp      </dev/null > out/01_host.log 2>&1
timeout 3000 gp -q 01_host_wide.gp </dev/null > out/01_host_wide.log 2>&1
timeout 3000 gp -q 01_scan.gp      </dev/null > out/01_scan.log 2>&1
timeout 3000 gp -q 01_scan_wide.gp </dev/null > out/01_scan_wide.log 2>&1

# exact recurrences and singular sets of the decayer generating functions
timeout 3000 gp -q 02_gen_seqs.gp     </dev/null > out/02_gen_seqs.log 2>&1
timeout 3000 gp -q 02_gen_halves.gp   </dev/null > out/02_gen_halves.log 2>&1
timeout 6000 gp -q 02_run_all.gp      </dev/null > out/02_run_all.log 2>&1

# (b) the conductor-12 Catalan row
timeout 3000 gp -q 03_chi12_build.gp        </dev/null > out/03_chi12_build.log 2>&1
timeout  900 gp -q 03b_chi12_productform.gp </dev/null > out/03b_chi12_productform.log 2>&1
timeout 3000 gp -q 04_chi12_rates.gp        </dev/null > out/04_chi12_rates.log 2>&1
timeout 3000 gp -q 05_chi12_padic.gp        </dev/null > out/05_chi12_padic.log 2>&1
timeout 3000 gp -q 06_chi12_full.gp         </dev/null > out/06_chi12_full.log 2>&1
timeout 3000 gp -q 07_chi12_denbound.gp     </dev/null > out/07_chi12_denbound.log 2>&1
timeout 3000 gp -q 09_catalan_3adic.gp      </dev/null > out/09_catalan_3adic.log 2>&1
timeout 3000 gp -q 10_chi12_small_alpha.gp  </dev/null > out/10_chi12_small_alpha.log 2>&1
timeout 1200 gp -q 11_chi12_primes.gp       </dev/null > out/11_chi12_primes.log 2>&1
timeout 1500 gp -q 12_c6_primes.gp          </dev/null > out/12_c6_primes.log 2>&1
timeout 2400 gp -q 13_rho_int.gp            </dev/null > out/13_rho_int.log 2>&1

# margins
python3 08_adelic.py > out/08_adelic.log 2>&1
```

Key single-line checks:

```sh
# the two-prime Catalan alignment certificate (463 / 176 at b=60)
grep -A6 'alpha  b    v2' out/06_chi12_full.log
# CDT's margin is unchanged by the adelic term
grep 'MEASURED slopes' out/08_adelic.log
# the level-8 / level-16 collapse
grep 'conductor 12, alpha=1  ' out/08_adelic.log
```

Nothing here needs **snake**; the longest run is about 25 minutes on the dev box.

### File map

| file | what |
|---|---|
| `01_lib.gp` | exact-$\mathbf Q$ rebuild of CDT's $H_A,H_B,H_C$, the symmetrisation, `to_y`, and the seven pure functions (port of `cdt_finder/indep_check2.py` off the finite field) |
| `01_check.gp`, `01_main.gp`, `01_main_big.gp` | validation against CDT's published coefficients and against `indep_check2.py` mod $2^{61}-1$; the slope tables of §2.1 |
| `01_host.gp`, `01_host_wide.gp` | row $\mathbf C$: $\sigma_3=2$, $\sigma_2=0$, and the slope of $B-\xi_3A$ |
| `01_scan.gp`, `01_scan_wide.gp` | the rational-direction scans of §2.2 |
| `02_gen_seqs.gp`, `02_gen_halves.gp`, `02_gen_c12_long.gp`, `02_fit_recur.gp`, `02_run_all.gp` | exact joint $(Q,P)$ recurrence fits over $\mathbf Q$ (mod-$p$ nullity sweep $+$ rational reconstruction $+$ held-out verification), characteristic polynomials and singular sets |
| `lib12.gp` | the conductor-12 library: `c12raw`, `c12cond`, `c12solve`, `c12row` (form I), `c12half`, `c12rowB` (form II, the product form) |
| `03_chi12_build.gp`, `03b_chi12_productform.gp` | the identity $\sum R=QG-P$ against direct summation, both forms |
| `04_chi12_rates.gp`, `06_chi12_full.gp`, `10_chi12_small_alpha.gp` | rate profiles $(\Lambda,\lambda,\kappa_2,\kappa_3,\nu,\eta)$ |
| `05_chi12_padic.gp`, `06_chi12_full.gp` | the $p$-adic limits against $\zeta_2(2)$ and $\tfrac9{10}L_3(2,\chi_{12})$, with controls |
| `13_rho_int.gp` | stability of the invariant $\rho_{\mathrm{int}}$ |
| `07_chi12_denbound.gp`, `11_chi12_primes.gp`, `12_c6_primes.gp` | denominator-type tests: rates, and the prime support of $\mathrm{den}(Q_b),\mathrm{den}(P_b)$ |
| `09_catalan_3adic.gp` | $\sigma_3=0$ on Zagier $\mathbf E$ and Zudilin's row |
| `08_adelic.py` | every margin in §§2–6; reuses `cdt_finder/cdt_bound.py` and `adelic_holonomy/adelic_bound.py` |

---

## 9. Corrections and additions this note forces on earlier documents

1. `MULTI_PRIME_LATTICE.md` §8 open item 7 ("whether a conductor-12 analogue
   exists for $L(2,\chi_{-4})$ … it would be a two-prime Catalan decayer — but
   there is no modular Catalan row with $\sigma_3>0$ to align it with").  **Both
   halves are now settled**: the construction exists (§4.1), has
   $\sigma_2=8,\sigma_3=3$ and the predicted avatars $\zeta_2(2)$ and
   $\tfrac9{10}L_3(2,\chi_{12})$ (§4.4); and the partner provably does not exist
   in the modular corpus (§4.5).  The "interesting half" — building a $3$-adic
   Catalan *modular* row — remains open and is now the only route.
2. `ADELIC_HOLONOMY.md` §5's structural ceiling
   $\log(256/\lambda_2)+\sigma_p\log p\le\log(256\lambda_1)$ should be read
   alongside the sharper §6 identity here: for the *pure module*, which is the
   only unconditional source of $\gamma$, the multi-prime total is exactly
   $\log256$, with **no** dependence on $\lambda_1$ and no slack.
3. `CDT_NONCONGRUENCE.md` §6's table gives $\sum_p\gamma_p=0.5837$ for
   $\sqrt{s_{18}}$ ($v_2=3$, $v_3=1$) — the corpus's only genuine two-prime
   $\gamma$.  §6 here explains why it does not help: the same $\log24$ is
   subtracted from the archimedean ceiling at weight $1$.
4. `CDT_FINDER.md` §2's calibration table can now record one more line: the
   adelic term on CDT's own host is $0$, measured and not merely argued.
5. `MULTI_PRIME_LATTICE.md` §8 open item 1 ("the recurrence of the conductor-6
   family … no analogue was fitted here, hence exact $\Lambda(\alpha)$,
   $\lambda(\alpha)$") is **closed at $\alpha=1$ and $\alpha=2$** by §3.1a:
   at $\alpha=2$ the joint $(Q,P)$ recurrence has order $2$, degree $15$,
   $\chi=(x-1)(x+1024)$, and $Q_b$ alone satisfies a first-order recurrence with
   the explicit ratio given there; at $\alpha=1$ the sequence is factorially
   decaying and there is no characteristic polynomial.  The §2.2 readings
   $\log\Lambda=6.775$ ($\alpha=2$) and $-2.062$ ($\alpha=1$) should be
   corrected to $\log1024=6.9315$ and to $-\infty$ respectively.
