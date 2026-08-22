# Multi-prime two-row lattices: a census of alignment primes, one new two-prime decayer, and what the second prime is worth

*Claude (Opus 5), 2026-08-22. Scripts and logs: `lattice/multi_prime/`.
Tags: **[proved]** = exact identity; **[verified]** = exact rational / $p$-adic computation over a
stated range, no floating point in the arithmetic; **[measured]** = numerical rate at finite index;
**[open]** as usual.
**No irrationality statement is claimed or implied anywhere below.** Every $\delta>1$ reading in
this file is reproduced verbatim by a rational surrogate for the period (§6) and therefore carries
no information about that period, exactly as `CATALAN_AUDIT.md` §4(a) and
`paper/sections/05_two_row.tex` Remark "the hypothesis $F>0$ is not decorative" require.*

---

## 0. Verdict

1. **Every decayer previously known to this project has exactly one alignment prime.**
   Zudilin's Catalan row: $\{2\}$, $\kappa_2=4$, $\sigma_2=8$. Nesterenko's $(4,7)$ row: $\{2\}$,
   $\kappa_2=14$, $\sigma_2=28$. The $\chi_{-3}$ conductor-3 hypergeometric family of
   `ONE_CLASS_TWO_WORLDS.md`: $\{3\}$ only. $\mathbf T$, Domb, Brown–Zudilin, $\eta$, AESZ 184:
   integral, one prime or none. **[verified]** §1. The unique row *among the twelve sporadic
   families* with two slope primes is $\mathbf F$ $(17,6,72)$, $c=72=2^3\cdot3^2$ — and it is
   integral, hence not a decayer. (Outside that list there is a second two-prime family, the
   cusp-move rows D.2, D.3, D.4, D.6 of §6.5, and it has no decayer either.)
   This confirms `paper/sections/05_two_row.tex` Remark (multi-prime bridging) *within* the
   canonical classification and extends it to the two hypergeometric decayers.

2. **$\kappa_p>0$ at two primes does not give slopes at two primes.** The conductor-3
   family at $\alpha\ne2$ carries genuine $2$-power denominators ($\kappa_2\approx1.4$–$5.7$ per
   unit of $k$, **[verified]** $k\le18$) — a fact not recorded in `ONE_CLASS_TWO_WORLDS.md`, which
   measured only the aggregate prime-to-3 rate $\nu$ — yet $\sigma_2=0$ exactly: the increments
   $v_2(P_k/Q_k-P_{k-1}/Q_{k-1})$ stay in $[-17,-8]$ for every member tested. The Casoratian eats
   the denominators, $w_2=-2\kappa_2$. **[verified]** §1.3. *The invariant that decides is
   $\sigma_p=w_p+2\kappa_p$, never $\kappa_p$ alone* — `ONE_CLASS_TWO_WORLDS.md` §6.1's corrected
   criterion, here confirmed at a second prime.

3. **A genuine two-prime decayer exists, and it is new: the conductor-6 well-poised row.**
   Theorem 1 of `ONE_CLASS_TWO_WORLDS.md` (double poles at thirds) generalises to **double poles at
   sixths**, giving exact two-term forms $Q_bL(2,\chi_{-3})-P_b$ whose denominators are supported at
   *both* $2$ and $3$. **[proved]** (§2, and **[verified]** against the direct sum $\sum_{t\ge0}R(t)$
   to $10^{-44}$). It has
   $$\boxed{\ \sigma_2=8,\qquad \sigma_3=3\quad\text{per unit of }b\ }$$
   uniformly across the family **[measured]**, and its two $p$-adic limits are exactly the values
   Theorem F predicts:
   $$\xi_3=\zeta_3(2)=2\,\xi_3^{\mathbf C}=2\,\xi_3^{\mathbf B}=\tfrac85\,\xi_3^{\mathbf F},
   \qquad
   \xi_2=\tfrac45L_2(2,\chi_{12})=\tfrac85\,\xi_2^{\mathbf F} .$$
   **[verified]** to $v_3=352$ and $v_2=951$ digits at $b=120$, with every competing rational
   scalar giving $v\le1$ (§3). This is the first row in the project that aligns with the modular
   world at two primes at once.

4. **The multi-prime lattice is therefore realisable, and it is worth a large, measured amount.**
   On the same pair of rows, the empirical box-constrained quality goes
   $$\delta_{\rm emp}: \quad 0.51\ (\text{no bridge})\ \longrightarrow\ 0.65\ (p=3)\ \longrightarrow\
   1.05\ (p=2\text{ and }3),$$
   with the divisor rate $\tfrac1n\log T_n$ rising from $0$ to $4.27$ to $9.84$ nats. **[measured]**
   §5. That is the largest single bridge gain measured in this project after the Catalan $G_2$
   bridge.

5. **But nothing follows.** The $\delta>1$ readings all sit at $F<0$, and the rational-surrogate
   control reproduces them **to every printed digit** (§6). In the regime where
   Theorem (two-minimum selection) actually applies, $F>0$, the honest best is
   $$\delta\approx0.94\ \ (\text{two primes, }\mathbf{c6}_{1/2}\times\mathbf{c6}_{1/3},\ n\ge60),$$
   against $0.9191$ (formula) / $0.893$ (measured) for the best single-prime construction in the
   class. So the second prime buys roughly $+0.02$ to $+0.05$ of *honest* worthiness on top of the
   existing record for $L(2,\chi_{-3})$, and the class still does not reach $\delta>1$ legitimately.

6. **A second two-slope-prime modular family exists (the $\sqrt{\text{Domb}}$ cusp-move rows
   D.2, D.3, D.4, D.6), and it too has no decayer.** $\sigma_2=4$ or $2$, $\sigma_3=1$, both
   primes produced by the *non-Zagier* linear term in the trailing coefficient rather than by
   $v_p(c)$ **[verified]** §6.5.3. The second prime is worth $+0.098$ there ($0.450\to0.548$) —
   the only place in the project where both rows of a pair are modular and both bridges fire — but
   $|\lambda_2|\in\{4,12\}$, so the linear forms grow and there is no content. Separately, the
   *mixed-period* rows D.5, D.7 have only **one** slope prime, and their $2$-adic limits are
   $\tfrac{15}{16}\Lambda_2$ and $\tfrac{15}{32}\Lambda_2$ — a new confirmation of the
   **linearity** of Theorem F, since the $\zeta(2)$ component of a mixed period has $p$-adic
   avatar $0$ at every prime **[verified]**, $947$ digits, §6.5.2.

7. **A cancellation identity explains why an integral modular engine can never gain from a second
   prime.** For an integral order-2 engine with roots $\Lambda\ge|\rho_2|$, $c=\Lambda\rho_2$,
   $\sigma_p=v_p(c)$,
   $$\textbf{net}(r)\;=\;\sum_p\min(\sigma^{\rm eng}_pr,\ \sigma^{\rm dec}_p)\log p\;-\;r\log\rho^{\rm eng}_2
   \;=\;r\log\Lambda_{\rm eng}-\sum_p\max\bigl(0,\ \sigma^{\rm eng}_pr-\sigma^{\rm dec}_p\bigr)\log p .$$
   **[proved]** §4.1. Every prime the engine contributes to the bridge is *exactly* repaid by the
   growth of its own linear form; the only net gain is $r\log\Lambda_{\rm eng}$, and it is capped by
   the **smallest** ratio $\sigma_p^{\rm dec}/\sigma_p^{\rm eng}$. Numerically: engine $\mathbf F$
   using both of its primes and engine $\mathbf C$ using its one prime give the *same* $\delta$ to
   four decimals, on every decayer tested. Real multi-prime gain needs a **non-integral engine**,
   where $\log\rho_2$ is not tied to $\sum_p\sigma_p\log p$ — which is precisely what the
   conductor-6 $\times$ conductor-6 pair supplies.

---

## 1. Census: alignment primes of every decayer

Definitions as in `paper/sections/05_two_row.tex` §5.3(ii) and `ONE_CLASS_TWO_WORLDS.md` §6.1.
For a row with linear form $Q_n\Theta-P_n$:
$$\kappa_p=-\lim_n\frac{v_p(Q_n)}{n}\ (\text{rate of }p\text{-power denominators}),\qquad
w_p=\lim_n\frac{v_p(\mathrm{Cas}_n)}{n},\qquad \sigma_p=w_p+2\kappa_p ,$$
and $\sigma_p>0$ iff $v_p(P_n/Q_n-P_{n-1}/Q_{n-1})$ grows linearly, i.e. iff a $p$-adic Apéry
limit exists.

### 1.1 The two Catalan decayers — one prime each

`lattice/multi_prime/census_1_zud.gp`, `census_2_nest.gp`.

| row | $\operatorname{den}(Q_n)$ | primes with $\kappa_p>0$ | $\sigma_p$ | $w_p$ |
|---|---|---|---|---|
| Zudilin $(Q_m,P_m)$ | pure power of $2$; $v_2=4m-2s_2(m)$ **exactly**, $0$ exceptions for $m\le300$ | $\{2\}$, $\kappa_2=4$ | $\sigma_2=8.000$ | $w_2=0$ |
| Nesterenko $(4,7)$ | pure power of $2$; $v_2=14n-O(\log n)$ | $\{2\}$, $\kappa_2=14$ | $\sigma_2=28.000$ | $w_2=0$ |

For $p=3,5,\dots,47$: $\kappa_p=0$ and the increment valuations are bounded
($[-6,1]$ resp. $[-9,0]$). $\operatorname{den}(P_n)$ does carry all small primes, but only at the
logarithmic $D_{2m-1}^2$ / $D_{6n}^2$ rate, so $\kappa_p=0$ there too. **[verified]**

### 1.2 The integral rows — the slope set is $\{p:p\mid c\}$

`lattice/multi_prime/census_3_introws.gp`, `census_6_bz.gp`. Measured slopes (wide-window
$(v(399)-v(100))/299$):

| row | $c$ | $\sigma_2$ | $\sigma_3$ | $\sigma_5$ | $\sigma_7$ |
|---|---|---|---|---|---|
| $\mathbf A(7,2,-8)$ | $-8$ | $2.973\to3$ | — | — | — |
| $\mathbf B(9,3,27)$ | $27$ | — | $2.967\to3$ | — | — |
| $\mathbf C(10,3,9)$ | $9$ | — | $1.967\to2$ | — | — |
| $\mathbf D(11,3,-1)$ | $-1$ | — | — | — | — |
| $\mathbf E(12,4,32)$ | $32$ | $4.987\to5$ | — | — | — |
| $\mathbf F(17,6,72)$ | $72$ | $2.977\to3$ | $1.967\to2$ | — | — |
| Cooper $s_{18}$ | $192$ | — | $0.957\to1$ | — | — |
| Domb $(10,4,64)$ | $64$ | $5.990\to6$ | — | — | — |
| $\mathbf T(12,4,16)$ | $16$ | $3.980\to4$ | — | — | — |
| Apéry $(17,5,1)$ | $1$ | — | — | — | — |
| $\delta(7,3,81)$ | $81$ | — | $3.940\to4$ | — | — |
| $\zeta(9,3,-27)$ | $-27$ | — | $2.940\to3$ | — | — |
| $\eta(11,5,125)$ | $125$ | — | — | $2.997\to3$ | — |
| Brown–Zudilin $\zeta(5)$ | — | none | none | none | none |

$\mathbf F$ is the **unique** row in the whole census with two slope primes, and it is integral
($\kappa\equiv0$), so it is an engine, never a decayer ($|\lambda_2|=8>1$). **[verified]**

### 1.3 $\kappa_p>0$ at $p=2$ without a slope: the conductor-3 family

`lattice/multi_prime/01_kappa_scan.gp`, `02_decayer_padic.gp`. Writing the family member as
$a=pk$, $b=qk$ (notation of `ONE_CLASS_TWO_WORLDS.md` §2.4), the $2$-part of
$\operatorname{den}(Q_k)$ is *not* trivial:

| $\alpha=a/b$ | $2$ | $9/5$ | $7/4$ | $5/3$ | $3/2$ | $4/3$ | $5/4$ | $1$ | $1/2$ |
|---|---|---|---|---|---|---|---|---|---|
| $\kappa_2$ (per $k$, at $k=18$) | $0$ | $1.67$ | $1.50$ | $1.61$ | $1.44$ | $3.67$ | $5.61$ | — | — |
| $\kappa_3$ (per $k$, at $k=18$) | $2.67$ | $12.72$ | $9.50$ | $6.67$ | $3.72$ | $4.67$ | $5.50$ | — | $0$ |

but the $2$-adic increments $v_2(\xi_k-\xi_{k-1})$ sit in $[-17,-8]$ for **every** member and every
$k$ tested (up to $b=72$), and $v_2\bigl(\xi_k-\tfrac85\xi_2^{\mathbf F}\bigr)$ likewise. So
$$\sigma_2=0\quad\text{exactly, i.e. } w_2=-2\kappa_2 . \qquad\textbf{[verified]}$$
The conductor-3 family is a **one-prime** decayer, and the $2$-power denominators it carries are
pure cost. (The corresponding $3$-adic column is the known one: $\sigma_3=3$ per unit of $b$,
uniformly, with $v_3(\xi_k-2\xi_3^{\mathbf C})=207$ at $b=72$ for $\alpha=5/3$ **[verified]**.)

---

## 2. The conductor-6 well-poised row

### 2.1 Statement

> **Theorem $1'$ (conductor 6).** **[proved]** Fix $b\ge0$ and set
> $$D_b(t)=\prod_{j=0}^{b}\Bigl(t+j+\tfrac16\Bigr)^{2}\Bigl(t+j+\tfrac56\Bigr)^{2},\qquad
> \deg D_b=4b+4 .$$
> Let $N\in\mathbf Q[t]$ satisfy $\deg N\le4b+2$ and $N(-t-b-1)=-N(t)$, and put $R=N/D_b$. Then
> $\sum_{t\ge0}R(t)$ converges and
> $$\sum_{t\ge0}R(t)=Q\,L(2,\chi_{-3})-P,\qquad Q,P\in\mathbf Q,\qquad
> Q=45\sum_{j=0}^{b}A_j^{(1/6)} ,$$
> $A^{(1/6)}_j$ the double-pole coefficient at $t=-(j+\tfrac16)$. No $\zeta(2)$, no $\pi\sqrt3$,
> no $\log$.

*Proof.* Verbatim the proof of Theorem 1 of `ONE_CLASS_TWO_WORLDS.md` §2.2 with the three
conductor-6 evaluations
$$\zeta(2,\tfrac16)+\zeta(2,\tfrac56)=24\,\zeta(2),\qquad
\zeta(2,\tfrac16)-\zeta(2,\tfrac56)=36\,L(2,\chi^{(6)}_{-3})=45\,L(2,\chi_{-3}),$$
$$\psi(\tfrac16)-\psi(\tfrac56)=-\pi\cot(\pi/6)=-\pi\sqrt3 ,$$
the second because the imprimitive character mod $6$ contributes the Euler factor
$1-\chi_{-3}(2)2^{-2}=\tfrac54$. The involution $\sigma(t)=-t-b-1$ maps $-(j+\tfrac16)$ to
$-\bigl((b-j)+\tfrac56\bigr)$ because $\tfrac16+\tfrac56=1$, so it permutes the pole set exchanging
the two classes, whence $A_{\sigma(i)}=-A_i$ and $B_{\sigma(i)}=+B_i$; the $\zeta(2)$ and
$\pi\sqrt3$ coefficients vanish and $45\,\Sigma A^{(1/6)}\cdot L$ survives. $\square$

**Note the arithmetic significance of the constant $45=36\cdot\tfrac54$**: the Euler factor
$\mathcal E_2(2)=\tfrac54$ of `EULER_CRITERION.md` is *built into the normalisation*, which is
exactly why this row has a $2$-adic avatar while the conductor-3 row does not.

Implementation: `chi6row(a,b)` in `lattice/multi_prime/lib.gp`, exact rational arithmetic, with the
two vanishings asserted at runtime (they hold for every instance computed). **[verified]** against
the direct sum $\sum_{t=0}^{20000}R(t)$: agreement $1.9\times10^{-44}$ at $(a,b)=(4,4)$, and
$3.1\times10^{-18}$ / $2.5\times10^{-9}$ at the smaller instances where the analytic tail
$\asymp T^{2a+2-4b-4}$ is itself the limiting error.

First values, $N_{a,b}(t)=(2t+b+1)\prod_{j=1}^{a}(t-j+1)(t+b+j)$:
$$Q_{0,0}=\tfrac{135}2,\quad Q_{1,1}=\tfrac{1107}4,\quad Q_{2,1}=-\tfrac{17787}{16},\quad
Q_{3,2}=-\tfrac{112969233}{8192},\quad Q_{4,4}=\tfrac{581850843843}{161480704}.$$

### 2.2 The rate profile

`lattice/multi_prime/04_chi6_rates.gp`, all rates per unit of $b$, at the largest $b$ reached.

| $\alpha$ | $b$ | $\log\Lambda$ | $\log\lambda$ | $\kappa_2$ | $\kappa_3$ | $\nu$ | $\eta$ | $\sigma_2$ | $\sigma_3$ | $w_2$ | $w_3$ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| $2$ | $60$ | $6.775$ | $-0.232$ | $7.82$ | $2.85$ | $0.66$ | $17.76$ | $7.85$ | $2.88$ | $-7.78$ | $-2.82$ |
| $3/2$ | $80$ | $1.599$ | $-4.487$ | $6.91$ | $1.81$ | $3.03$ | $20.85$ | $7.89$ | $2.84$ | $-5.94$ | $-0.79$ |
| $5/3$ | $96$ | $3.205$ | $-1.877$ | $7.28$ | $2.24$ | $2.22$ | $20.11$ | $7.96$ | $2.89$ | $-6.60$ | $-1.59$ |
| $4/3$ | $96$ | $-0.312$ | $-5.396$ | $6.62$ | $1.57$ | $4.50$ | $21.92$ | $7.97$ | $2.90$ | $-5.26$ | $-0.25$ |
| $1$ | $40$ | $-2.062$ | $-9.146$ | $5.83$ | $0.83$ | $5.48$ | $21.96$ | $7.83$ | $2.78$ | $-3.83$ | $+1.13$ |
| $2/3$ | $96$ | $-6.977$ | $-12.053$ | $5.28$ | $0.24$ | $8.77$ | $24.32$ | $7.97$ | $2.90$ | $-2.59$ | $+2.42$ |
| $1/2$ | $64$ | $-7.371$ | $-13.021$ | $4.95$ | $0.00$ | $8.34$ | $23.42$ | $7.95$ | $2.88$ | $-1.95$ | $+2.88$ |
| $1/4$ | $80$ | $-10.198$ | $-14.495$ | $4.41$ | $0.00$ | $9.84$ | $24.67$ | $7.90$ | $2.83$ | $-0.93$ | $+2.83$ |

Here $\eta=\tfrac1b\log\operatorname{lcm}(\operatorname{den}Q_b,\operatorname{den}P_b)$ is the full
integralising cost and $\nu$ the prime-to-$\{2,3\}$ part of $\operatorname{den}Q_b$.
$\sigma_2$ and $\sigma_3$ are read off the alignment valuations of §3, not fitted.

**The uniformity is the structural fact:** $\sigma_2=8$ and $\sigma_3=3$ per unit of $b$ for the
whole family, exactly as $\sigma_3=3$ per unit of $b$ is uniform for the conductor-3 family
(`ONE_CLASS_TWO_WORLDS.md` §5). Consequently $w_p=\sigma_p-2\kappa_p$ and the *net* resource at $p$
is $(\sigma_p-\kappa_p)\log p$, which for $\alpha=2$ is $\approx0$ at both primes (the row pays for
its own denominators, the pathology of `ONE_CLASS_TWO_WORLDS.md` §6.1) and rises to
$3.5\log2+2.8\log3=5.5$ at $\alpha=1/4$.

---

## 3. Two-prime alignment, verified

`lattice/multi_prime/08_best.gp`. Reference rows $\mathbf B,\mathbf C,\mathbf F$ generated exactly
to $N=700$, giving $v_3$ Cauchy precision $1378$ and $v_2$ Cauchy precision $2064$.

**Prediction (Theorem F, `EULER_CRITERION.md`).** The conductor-6 row has $r_\infty=1$ for
$\Theta=L(2,\chi_{-3})$, $\psi=\chi_{-3}$, $w=1$. Since $\chi_{-3}(3)=0$ and $\chi_{-3}(2)=-1$,
$$\mathcal E_3(2)=1,\quad r_3=1,\quad \xi_3=\Lambda_3=\zeta_3(2);\qquad
\mathcal E_2(2)=\tfrac54,\quad r_2=\tfrac45,\quad \xi_2=\tfrac45\Lambda_2=\tfrac45L_2(2,\chi_{12}).$$
With $r_\infty^{\mathbf C}=\tfrac12$, $r_\infty^{\mathbf F}=\tfrac58$ this says
$\xi_3=2\xi_3^{\mathbf C}$ and $\xi_2=\tfrac85\xi_2^{\mathbf F}$ — *the same rational ratio $8:5$
against $\mathbf F$ at both primes.*

**Measured** ($b=120$, four different family members):

| $\alpha$ | $v_3(\xi-2\xi_3^{\mathbf C})$ | $v_3(\xi-2\xi_3^{\mathbf B})$ | $v_2(5\xi-8\xi_2^{\mathbf F})$ |
|---|---|---|---|
| $2$ | $352$ | $352$ | $951$ |
| $3/2$ | $352$ | $352$ | $947$ |
| $1/2$ | $352$ | $352$ | $947$ |
| $1$ | $352$ | $352$ | $947$ |

**Controls** (same $\xi$, wrong scalar): $v_3(\xi-\xi_3^{\mathbf C})=-1$,
$v_3(\xi-4\xi_3^{\mathbf C})=-1$, $v_2(\xi-\xi_2^{\mathbf F})=0$, $v_2(\xi-2\xi_2^{\mathbf F})=1$;
a ten-scalar sweep at $b=60$ gives $v_2\le4$ for every $c\ne\tfrac85$ and $v_3\le0$ for every
$c\ne2$. **[verified]**

So the requested $\ge200$ $p$-adic digits are met at **both** primes ($352$ and $951$), and the
alignment is with $\mathbf B$, $\mathbf C$, $\mathbf F$, $s_{18}$ at $p=3$ and with $\mathbf F$
alone at $p=2$ ($\mathbf B$, $\mathbf C$, $s_{18}$ have $\sigma_2=0$; $\mathbf A$ carries a
different period and $\xi_2^{\mathbf A}=0$; $\mathbf E$ carries $G$).

---

## 4. The master formula with $G=\sum_p\min(\cdot,\cdot)\log p$

Notation of `paper/sections/05_two_row.tex` eq. (master): decayer at index $\gamma n$, engine at
$\alpha n$, $r=\alpha/\gamma$,
$$F=\tfrac12\Bigl[K+r\log\rho^{\rm eng}_2+\log\lambda_{\rm dec}-G\Bigr],\qquad
H=F+\log\frac{\Lambda_{\rm dec}}{\lambda_{\rm dec}},\qquad \delta=1-\frac FH,$$
$$G=\sum_{p}\min\bigl(\sigma^{\rm eng}_pr,\ \sigma^{\rm dec}_p\bigr)\log p ,$$
$K$ the rate of the common denominator clearing (measured as an exact `lcm`, not estimated).

### 4.1 The cancellation identity

> **Proposition (multi-prime cancellation).** **[proved]** Let the engine be an *integral* order-2
> row with characteristic roots $\Lambda\ge|\rho_2|$, $c=\Lambda\rho_2$, $\sigma^{\rm eng}_p=v_p(c)$.
> Then $\sum_p\sigma^{\rm eng}_p\log p=\log|c|=\log\Lambda+\log|\rho_2|$, hence
> $$\textbf{net}(r):=G-r\log\rho^{\rm eng}_2
> = r\log\Lambda_{\rm eng}\;-\;\sum_p\max\bigl(0,\ \sigma^{\rm eng}_pr-\sigma^{\rm dec}_p\bigr)\log p
> \;\le\; r\log\Lambda_{\rm eng},$$
> with equality iff $r\le\min_p\sigma^{\rm dec}_p/\sigma^{\rm eng}_p$.

Two consequences, both confirmed numerically in `out/05_design.log`:

* **Adding a prime to an integral engine is exactly self-cancelling.** Row $\mathbf F$ has
  $\rho_2=8$ and $\sigma_2=3$, so its $p=2$ bridge contributes $3r\log2=r\log8$ — precisely the
  cost of its own growing linear form. Engine $\mathbf F$ bridging at *both* primes and engine
  $\mathbf C$ bridging at its *one* prime give the identical $\delta$: e.g. against the conductor-6
  $\alpha=3/2$ decayer, $0.5510$ vs $0.5510$; against $\alpha=4/3$, $0.5680$ vs $0.5680$; against
  $\alpha=1$, $0.5954$ vs $0.5954$. To four decimals, on every member.
* **The binding constraint is the *smallest* ratio.** For $\mathbf F$,
  $\min(\sigma_2^{\rm dec}/3,\ \sigma_3^{\rm dec}/2)=\min(2.65,1.45)=1.45$, the same value as for
  $\mathbf C$. The second prime adds a constraint, not headroom.

**Therefore a real multi-prime gain requires a non-integral engine**, for which $\log\rho_2$ is
unrelated to $\sum_p\sigma_p\log p$. In the class $L(2,\chi_{-3})$ that means a *second*
conductor-6 row.

### 4.2 Formula values, exact `lcm` denominators

`lattice/multi_prime/06_design_full.gp` (grid $r\in\{0.25,\dots,3.5\}$). Selected rows; full table
in `out/06_design_full.log`.

| decayer | engine | $\delta$(no bridge) | $\delta(p=3)$ | $\delta(p=2)$ | $\delta$(both) |
|---|---|---|---|---|---|
| $\mathbf{c6}_{2/1}$ | $\mathbf C$ | $0.4451$ | $0.4719$ | $0.4451$ | $0.4719$ |
| $\mathbf{c6}_{2/1}$ | $\mathbf F$ | $0.4374$ | $0.4450$ | $0.4446$ | $\mathbf{0.4718}$ |
| $\mathbf{c6}_{2/1}$ | $\mathbf{c6}_{1/4}$ | $0.4752$ | $0.5328$ | $0.5844$ | $\mathbf{0.6740}$ |
| $\mathbf{c6}_{3/2}$ | $\mathbf{c6}_{1/4}$ | $0.4664$ | $0.5182$ | — | $\mathbf{0.6}$–$0.75$ |
| $\mathbf{c3}_{5/3}$ | $\mathbf C$ | $0.7845$ | $0.8049$ | $0.7845$ | $0.8049$ |
| $\mathbf{c3}_{1/2}$ | $\mathbf C$ | $0.8628$ | $0.8865$ | $0.8628$ | $0.8865$ |

(The $\mathbf{c3}\times\mathbf{c3}$ and some $\mathbf{c6}\times\mathbf{c6}$ cells of the full log
report $\delta>1$; every one of them has $F<0$ and is void — see §6.)

---

## 5. The actual lattice

`lattice/multi_prime/07_lattice.gp`, `08_best.gp`, `09_control.gp`. Exact integers throughout,
in the architecture of `paper/sections/05_two_row.tex` §5.1:

$$S_Q=\operatorname{lcm}_i\operatorname{den}(Q_i),\quad \alpha_i=S_QQ_i\in\mathbf Z;\qquad
S_P=\operatorname{lcm}_i\operatorname{den}(S_QP_i),\quad Y_i=S_QS_PP_i,\ X_i=S_P\alpha_i ;$$
$$T=\prod_{p\in\text{bridge}}p^{\min_{i<j}v_p(\alpha_iY_j-\alpha_jY_i)},\qquad M=S_PT,$$
$$K_n=\Bigl\{c\in\mathbf Z^k:\ \textstyle\sum c_i\alpha_i\equiv0\ (T),\ \sum c_iY_i\equiv0\ (M)\Bigr\},
\qquad q=\frac{c\cdot X}{M},\ \ p=\frac{c\cdot Y}{M}\in\mathbf Z .$$
Selection is inside the prescribed anisotropic box $|c_1|\le e^{xn}$, $|c_2|\le e^{(\sigma-x)n}$,
$x=\tfrac12(\sigma+E_2-E_1)$ — **never** unconstrained LLL (`ZETA3_TWO_LATTICE.md` §8).

### 5.1 Pipeline control: $\zeta(3)$, Domb $\times$ $\mathbf T$ at $2{:}3$

Known answer $\delta=0.9009532$, $F=+1.1627$ (`ZETA3_TWO_LATTICE.md` §2).

| $n$ | $\tfrac1n\log T$ | $\sigma$ | $F$ | $\delta_{\rm formula}$ | $\tfrac1n\log|q|$ | $\tfrac1n\log|q\zeta(3)-p|$ | $\delta_{\rm emp}$ |
|---|---|---|---|---|---|---|---|
| $20$ | $8.075$ | $16.130$ | $0.604$ | $0.9160$ | $11.082$ | $+0.621$ | $0.9440$ |
| $40$ | $8.196$ | $16.640$ | $0.816$ | $0.8905$ | $11.354$ | $+0.746$ | $0.9343$ |
| $60$ | $8.399$ | $17.064$ | $0.859$ | $0.8856$ | $11.416$ | $+0.855$ | $0.9251$ |
| $80$ | $8.257$ | $17.012$ | $0.993$ | $0.8701$ | $11.520$ | $+0.984$ | $0.9146$ |
| theory | $8.3178$ | $17.3178$ | $1.1627$ | $0.90095$ | $11.7392$ | $+1.1627$ | $0.90095$ |

Both columns converge to the published values, and $\tfrac1n\log|q\zeta(3)-p|>0$ throughout — the
linear forms grow, exactly as they must. **The pipeline is calibrated.** **[measured]**

### 5.2 The two-prime ladder on one pair of conductor-6 rows

Decayer $\mathbf{c6}_{1/4}$, engine $\mathbf{c6}_{1/3}$, both at index $n$ (sampling $1{:}1$).
Measured $\tfrac1n v_2(h_n)=8.04$ and $\tfrac1n v_3(h_n)=3.89$ at $n=108$ — the two-prime bridge is
present at the predicted rates simultaneously.

| bridge | $\tfrac1n\log T$ | $\sigma$ | $F$ | $\tfrac1n\log|q|$ | $\tfrac1n\log|q\Theta-p|$ | $\delta_{\rm emp}$ |
|---|---|---|---|---|---|---|
| none | $0.000$ | $11.926$ | $+4.706$ | $9.512$ | $+4.688$ | $0.5072$ |
| $p=3$ | $4.272$ | $16.198$ | $+2.570$ | $7.390$ | $+2.570$ | $0.6522$ |
| $p=2$ and $p=3$ | $9.844$ | $21.769$ | $-0.216$ | $4.589$ | $-0.220$ | $1.0480$ |

(at $n=108$; the same ladder at $n=24,36,\dots,96$ is in `out/09_control.log`.)

**This is the first end-to-end multi-prime two-row lattice in the project.** The divisor rate
$9.844=8.04\log2+3.89\log3$ is the sum of the two primes' contributions, the Smith-normal-form
index bound goes through verbatim for $T=2^{a}3^{b}$, and the outputs $q,p$ are integral at every
$n$ tested.

### 5.3 Three rows, two primes

`07_lattice.gp` §(D): decayer $\mathbf{c6}_{1/2}$ together with $\mathbf C$ (which bridges only at
$p=3$) and $\mathbf F$ (which bridges at $p=2$ and $p=3$), rank-3 correlated lattice, box
$x_i=\sigma/3+(\bar E-E_i)$.

| $n$ | bridge | $\tfrac1n\log T$ | $\sigma$ | $\delta_{\rm emp}$ |
|---|---|---|---|---|
| $48$ | $p=3$ | $2.220$ | $13.848$ | $0.3749$ |
| $48$ | $p=2,3$ | $5.484$ | $17.112$ | $0.2258$ |

**The three-row configuration is worse than the two-row one at every $n$ tested, and adding the
second prime makes it worse still.** Two reasons, both visible in the numbers: (i) the bridge
modulus is the *minimum* over all three pairs of mixed minors, and the $\mathbf C\times\mathbf F$
minor has only $\min(\sigma_2^{\mathbf C},\sigma_2^{\mathbf F})=0$ at $p=2$, so the third row caps
the $p=2$ divisor at what $\mathbf C$ can supply; (ii) the rank-3 box is a guess — Theorem
(two-minimum selection) is a rank-2 theorem and its rank-3 analogue is **[open]**, exactly as
`ZETA3_TWO_LATTICE.md` §6.3 warned. The numbers above should be read as "no evidence that rank 3
helps", not as a measurement of a rank-3 theory.

---

## 6. The $F>0$ hypothesis and the rational-$\Theta^*$ control

Every configuration in §5.2 with $\delta>1$ has $F<0$. `paper/sections/05_two_row.tex`
Remark (`rem:Fpos`) and `CATALAN_AUDIT.md` §4(a) say that in that regime the construction
degenerates to Dirichlet and proves nothing. The mandatory test: the congruence lattice never uses
$\Theta$, so replacing $\Theta$ by a rational surrogate must reproduce the numerics.

`lattice/multi_prime/09_control.gp`, $\Theta^*=\operatorname{bestappr}(\Theta,10^{1200})$:

| $n$ | $F$ | $\tfrac1n\log|q|$ | $\tfrac1n\log|q\Theta-p|$ | $\delta_{\rm emp}$ | $\tfrac1n\log|q\Theta^*-p|$ | $\delta^*$ |
|---|---|---|---|---|---|---|
| $24$ | $-0.7016$ | $4.134$ | $-0.8793$ | $1.2127$ | $-0.8793$ | $1.2127$ |
| $48$ | $-0.3729$ | $4.430$ | $-0.3730$ | $1.0842$ | $-0.3730$ | $1.0842$ |
| $72$ | $-0.2744$ | $4.539$ | $-0.2827$ | $1.0623$ | $-0.2827$ | $1.0623$ |
| $96$ | $-0.1965$ | $4.608$ | $-0.2217$ | $1.0481$ | $-0.2217$ | $1.0481$ |
| $108$ | $-0.2157$ | $4.589$ | $-0.2203$ | $1.0480$ | $-0.2203$ | $1.0480$ |

**Identical to every printed digit.** $1/\operatorname{den}(\Theta^*)\approx10^{-1200}$ while
$|q\Theta^*-p|\approx e^{-24}\approx5\cdot10^{-11}$, so there is no contradiction and no
information: the $\delta>1$ readings say nothing whatever about $L(2,\chi_{-3})$.

Note also the shape of $F$: $-0.70,-0.49,-0.37,-0.27,-0.27,-0.40,-0.20,-0.22$ at
$n=24,\dots,108$ — drifting toward $0$ from below. The honest reading is that this configuration
sits at, or just below, the **knife-edge $F=0$, worthiness exactly $1$** — the same place the
Catalan Zudilin $\times$ Nesterenko bridge lands (`ZETA3_TWO_LATTICE.md` §15.2, §16), and with the
same void status.

**The honest $F>0$ numbers.** Restricting to configurations where $F>0$ at every $n$ tested (so
that Theorem (selection) applies):

| configuration (engine $\times$ decayer) | $n$ | bridge | $F$ | $\delta_{\rm emp}$ | log |
|---|---|---|---|---|---|
| $\mathbf{c3}_{5/3}\times\mathbf C$ at $r=2$ (single prime, previous record) | $36$ | $p=3$ | $+0.93$ | $0.8928$ | `08_best` |
| $\mathbf{c6}_{1/3}\times\mathbf{c6}_{1/2}$ | $48$ | $p=3$ | $+2.81$ | $0.6677$ | `08_best` |
| $\mathbf{c6}_{1/2}\times\mathbf{c6}_{1/3}$ | $72$ | $p=2,3$ | $+0.34$ | $\mathbf{0.9444}$ | `09_control` |
| $\mathbf{c6}_{1/2}\times\mathbf{c6}_{1/3}$ | $96$ | $p=2,3$ | $+0.48$ | $\mathbf{0.9229}$ | `09_control` |
| $\mathbf{c6}_{1/2}\times\mathbf{c6}_{1/3}$ | $108$ | $p=2,3$ | $+0.37$ | $\mathbf{0.9401}$ | `09_control` |

(For $n\le48$ the same pair has $F<0$ and $\delta_{\rm emp}=1.020$–$1.045$; it crosses into the
$F>0$ regime at $n\approx60$ and stays there.)

So in the $F>0$ regime the two-prime bridge takes this family of rows from $0.668$ to $0.92$–$0.94$, and
$0.94$ is the best honest quality now known for $L(2,\chi_{-3})$ — marginally above the
conductor-3 single-prime record of $0.9191$ (formula) / $0.893$ (measured). **Still $<1$.**

---

## 6.5 The mixed $\zeta(2)/L(2,\chi_{-3})$ class of the $\sqrt{\text{Domb}}$ cusp-move orbit

*Added at the coordinator's request, following `CUSP_MOVE_PROGRAM.md` §4.1, §6.
Scripts `lattice/multi_prime/10_sqrtdomb.gp`, `11_mixed_class.gp`, `12_mixed_delta.gp`.*

The eight-member $\sqrt{\text{Domb}}$ orbit splits into two arithmetically different halves, and
the distinction matters for this file. **The two-prime members and the mixed-period members are
not the same rows.**

### 6.5.1 The mixed-period rows D.5, D.7 have only ONE slope prime

$$\mathrm{D.5}=(20,10,64):\ (n{+}1)^2u_{n+1}=(20n^2{+}26n{+}10)u_n-64n^2u_{n-1},\quad
\xi=\tfrac1{16}\bigl(15L(2,\chi_{-3})-6\zeta(2)\bigr),$$
$$\mathrm{D.7}=(20,4,64):\ (n{+}1)^2u_{n+1}=(20n^2{+}14n{+}4)u_n-64n^2u_{n-1},\quad
\xi=\tfrac1{32}\bigl(6\zeta(2)+15L(2,\chi_{-3})\bigr).$$

Both archimedean limits reproduced to $30$ digits at $n=500$ **[verified]**. Their trailing
coefficient is the pure $64n^2$, so the Casoratian is $-64^{\,n}/(n+1)^2$ and $v_3=0$:

| row | $\sigma_2$ (measured) | $\sigma_3$ (measured) |
|---|---|---|
| D.5 | $2966/500=5.93\to6$ | $v_3(\text{increment})\in\{0,-2\}$ — **no slope** |
| D.7 | $2963/500=5.93\to6$ | **no slope** |

So the class carrying the $(-6,15)$ and $(6,15)$ directions in the $(\zeta(2),L)$-plane is a
**one-prime** class. **[verified]** $n\le500$.

### 6.5.2 But they *do* align $2$-adically with the conductor-6 decayer, and the reason is Theorem F

Theorem F assigns the $\zeta(2)$ period the **outer** placement $(\psi,\varphi)=(\mathbf1,\chi_{-3})$
with $\varphi$ odd, hence $\kappa_p=0$ and $\xi_p=0$ at *every* prime — this is the mechanism behind
$\xi_2^{\mathbf A}=0$ (`ACF_ONE_SURFACE.md` Thm 3(c); re-measured here, $v_2(\xi^{\mathbf A})=1474$
at $N=500$, i.e. $3N$). By linearity of $\xi_p$ in the source, the $\zeta(2)$ component of a mixed
period is **$p$-adically invisible**, so
$$\xi_2^{\mathrm{D.5}}=\tfrac{15}{16}\Lambda_2,\qquad
\xi_2^{\mathrm{D.7}}=\tfrac{15}{32}\Lambda_2,\qquad
\Lambda_2=\tfrac45L_2(2,\chi_{12})=\xi_2^{\text{c6}} .$$

**Measured** ($b=120$ conductor-6 row against the $N=500$ modular rows):
$$v_2\bigl(16\,\xi^{\mathrm{D.5}}-15\,\xi^{\text{c6}}\bigr)=947,\qquad
v_2\bigl(32\,\xi^{\mathrm{D.7}}-15\,\xi^{\text{c6}}\bigr)=947,$$
against controls $v_2(16\xi^{\mathrm{D.5}}-16\xi^{\text{c6}})=3$,
$v_2(32\xi^{\mathrm{D.7}}-16\xi^{\text{c6}})=3$, $v_2(\xi^{\mathrm{D.5}}-\xi^{\text{c6}})=-1$;
at $p=3$ everything is flat ($v_3=-10$). **[verified]**, and the same three $\alpha$'s of the
conductor-6 family give the identical valuation. This is a clean new confirmation of Theorem F's
*linearity* on a mixed-period row — the first in the project — but it is a **one-prime** alignment.

It also does **not** give a rank-2 two-row lattice: D.5 and the conductor-6 row carry *different*
periods, so their linear forms live in $\langle1,\zeta(2),L\rangle$, not in $\langle1,\Theta\rangle$.
The only $\zeta(2)$-free rank-2 combination available is
$$a^{(7)}_n\!\cdot\!16\bigl(a^{(5)}_n\Theta_5-b^{(5)}_n\bigr)+a^{(5)}_n\!\cdot\!32\bigl(a^{(7)}_n\Theta_7-b^{(7)}_n\bigr)
=30L\,a^{(5)}_na^{(7)}_n-16a^{(7)}_nb^{(5)}_n-32a^{(5)}_nb^{(7)}_n,$$
an $L(2,\chi_{-3})$-row with $\Lambda=256$, $\rho_2=64=2^{\sigma_2}$, $\sigma_2=6$. By the
cancellation identity of §4.1 its $p=2$ bridge contributes $\min(6r,\sigma_2^{\rm dec})\log2\le6r\log2
=r\log\rho_2$ — **exactly its own cost** — so the net resource is $\le0$ and the pairing is worse than
no bridge at all. Dead end, but for a reason that is now an identity.

### 6.5.3 The genuinely two-prime rows are D.2, D.3, D.4, D.6 — and they carry a different period

These four have a *linear term* in the trailing coefficient, so the Casoratian is not $c^{\,n}$ and
both primes survive:
$$\mathrm{D.2}=(-28,-6,192),\ Q(n)=192n^2-96n;\qquad \mathrm{D.3}=(8,0,-48),\ Q(n)=-48n^2+24n;$$
$$\mathrm{D.4}=(8,6,-48),\ Q(n)=-48n^2-24n;\qquad \mathrm{D.6}=(-28,-12,192),\ Q(n)=192n^2+96n.$$

| row | $\sigma_2$ | $\sigma_3$ | $\Lambda$ | $\lambda_2$ | $k$ | $a_n\in\mathbf Z$ |
|---|---|---|---|---|---|---|
| D.2 | $1975/500=3.95\to4$ | $489/500=0.98\to1$ | $-16$ | $-12$ | $2$ | yes |
| D.3 | $976/500=1.95\to2$ | $492/500=0.98\to1$ | $16$ | $-4$ | $2$ | yes |
| D.4 | $977/500=1.95\to2$ | $495/500=0.99\to1$ | $16$ | $-4$ | $2$ | yes |
| D.6 | $1971/500=3.94\to4$ | $495/500=0.99\to1$ | $-16$ | $-12$ | $2$ | yes |

**[verified]** $n\le500$. **This is the second family in the project with two slope primes, and the
first that is not $\mathbf F$.** Unlike $\mathbf F$ the two primes do *not* come from $v_p(c)$ — the
"$c$" here is $64$ resp. $-48$ and the $3$-part is produced by the non-Zagier linear term, exactly
the mechanism of `paper/sections/04_padic.tex` Remark (non-Zagier normalisations). The internal
relations of `CUSP_MOVE_PROGRAM.md` §6 are reproduced:
$$v_2\bigl(\xi^{(2)}+2\xi^{(3)}\bigr)=981,\quad v_3(\cdot)=496;\qquad
v_2\bigl(\xi^{(4)}-\xi^{(6)}\bigr)=979,\quad v_3(\cdot)=496\qquad(N=500).$$

**But there is no decayer to pair them with, and none of them decays.**

* $|\lambda_2|=12$ resp. $4$: the linear forms **grow**. All four are engines.
* Their period is not in $\langle1,\zeta(2),L(2,\chi_{-3})\rangle$: `lindep` on
  $[\xi,1,\zeta(2),L]$ at $60$ digits returns coefficients of height $\sim3\cdot10^{14}$ for every
  one of the four — the random level, i.e. **no relation**. (Consistent with
  `CUSP_MOVE_PROGRAM.md` §5.2's finding that their archimedean periods are $\mathbf Q$-independent.)
* Consequently no factorial-world decayer aligns: a scan of $14$ rationals $c$ in
  $\xi^{\mathrm{D.2}}-c\,\xi^{\text{c6}}$ and $\xi^{\mathrm{D.3}}-c\,\xi^{\text{c6}}$ gives
  $v_2\le0$ and $v_3\le1$ at every $c$ (including $\tfrac{15}{16},\tfrac{15}{32},\tfrac58,\tfrac85$),
  and the conductor-3 row likewise. **[verified]** The direction $(-6,15)$ / $(6,15)$ that the
  coordinator suggested belongs to D.5/D.7, which are the one-prime members.

Master formula over the four two-prime rows (`12_mixed_delta.gp`), best over $r$:

| pair | none | $p=3$ | $p=2$ | **both** |
|---|---|---|---|---|
| D.3 $\times$ D.4 | $0.4492$ | $0.4500$ | $0.4502$ | $\mathbf{0.5479}$ |
| D.3 $\times$ D.2 | $0.4484$ | $0.4492$ | $0.4609$ | $0.5073$ |
| D.2 $\times$ D.6 | $0.1131$ | $0.1134$ | $0.1206$ | $0.1566$ |

The second prime is worth $+0.098$ here — a real multi-prime gain, and it is the *only* place in
the project where **both** rows of a pair are modular and both bridges fire. But $\delta\le0.548$,
$F>0$ throughout, and the "decayer" grows like $4^n$: there is no irrationality content whatever.

### 6.5.4 What would close it

A factorial-world decayer carrying the mixed period $\Theta_5=\tfrac1{16}(15L-6\zeta(2))$ is
constructible in principle and is a well-posed target. In the conductor-6 residue framework of
§2, dropping the antisymmetry and imposing instead the **two** linear conditions
$$\Sigma B^{(1/6)}=0\qquad\text{(kills }\pi\sqrt3\text{)},\qquad
\Sigma A^{(5/6)}+7\,\Sigma A^{(1/6)}=0\qquad\text{(fixes }\zeta(2){:}L=-6{:}15\text{)}$$
gives, from
$\zeta(2,\tfrac16)=12\zeta(2)+\tfrac{45}2L$, $\zeta(2,\tfrac56)=12\zeta(2)-\tfrac{45}2L$,
$$\sum_{t\ge0}R(t)=192\,\Sigma A^{(1/6)}\cdot\Theta_5-P,\qquad P\in\mathbf Q .$$
**[proved]**, and the numerator space has dimension $4b+3$ against two conditions, so the family is
large. What is **[open]** is whether any member *decays* — the antisymmetry $R\circ\sigma=-R$ is
what made the conductor-3 and conductor-6 forms small, and it is exactly the condition being traded
away. This was not attempted here.

The same construction with the ratio condition $\Sigma A^{(5/6)}=\tfrac{?}{}\,\Sigma A^{(1/6)}$
tuned to D.2's period cannot even be posed until that period is identified; `lindep` says it is not
in $\langle1,\zeta(2),L\rangle$ and `CUSP_MOVE_PROGRAM.md` §8 lists identifying the weight-three
sources of these rows as open.

---

## 7. The census table

$\Theta$ = the period; "decayer" = a row of the class with $|\lambda_2|<1$ (raw linear form
decays); $\kappa$-primes = primes carrying $p$-power denominators at linear rate;
$\sigma$-primes = primes at which the row has a $p$-adic Apéry limit.

| class $\Theta$ | decayer (factorial world) | $\kappa$-primes | $\sigma$-primes of the decayer | aligned modular rows, per prime | single-prime $\delta$ | multi-prime $\delta$ | verdict |
|---|---|---|---|---|---|---|---|
| $\zeta(2)$ | none with $\kappa_p>0$ (Beukers/Rhin–Viola $\zeta(2)$ rows are integral; $\mathbf D$ has $c=-1$) | — | — | $p=2$: $\mathbf A$, but $\xi_2^{\mathbf A}=0$ (Theorem F, $\varphi=\chi_{-3}$ odd); $s_7,s_{10}$ no slope | n/a | n/a | **no aligned decayer**; the class's only slope row has a vanishing $p$-adic avatar |
| $\zeta(3)$ | $\mathbf T\,(12,4,16)$, integral, $\lambda_2=12-8\sqrt2$ | none | $\{2\}$, $\sigma_2=4$ | $p=2$: Domb ($\sigma_2=6$), ratio $4{:}3$, $\xi_2=\tfrac13\zeta_2(3)$ vs $\tfrac14\zeta_2(3)$ | $\mathbf{0.9010}$ | **impossible**: $c^{\rm Domb}=64$, $c^{\mathbf T}=16$ are $2$-powers | deficit: need $k\le k^*=2.2248$ (have $3$), or an extra $2^{3.355n}$ |
| $L(2,\chi_{-4})$ (Catalan) | Zudilin ($\kappa_2=4$); Nesterenko ($\kappa_2=14$) | $\{2\}$ | $\{2\}$, $\sigma_2=8$ resp. $28$ | $p=2$: $\mathbf E$, ratio $2{:}1$, $\xi_2=\zeta_2(2)$ vs $\tfrac12\zeta_2(2)$; level-16 rows | $\mathbf{0.9025}$ (Zud $\times\mathbf E$, $5{:}8$); knife-edge $1$ (Zud $\times$ Nest) | **impossible**: every Catalan-class row has $c$ a $2$-power; no modular Catalan row has $\sigma_3>0$ | deficit at $p=3$: no partner exists, though $\mathcal E_3(2)=\tfrac{10}9$ and $L_3(2,\chi_{12})\ne0$ leave the slot open |
| $L(2,\chi_{-3})$ | conductor-3 family ($\kappa_3\le12.7$, $\kappa_2\le5.6$) | $\{2,3,5,7,\dots\}$ | $\{3\}$ only ($\sigma_2=0$: $w_2=-2\kappa_2$) | $p=3$: $\mathbf B,\mathbf C,\mathbf F,s_{18}$ | $0.9191$ (formula) / $0.893$ (measured) | — | $\kappa_2>0$ is pure cost |
| $L(2,\chi_{-3})$ | **conductor-6 family (new, §2)** | $\{2,3,5,7,\dots\}$ | $\{2,3\}$: $\sigma_2=8$, $\sigma_3=3$ per unit $b$ | $p=3$: $\mathbf B,\mathbf C,\mathbf F,s_{18}$ ($352$ digits); $p=2$: $\mathbf F$ ($951$ digits) | $0.668$ (same rows, $p=3$ only) | $\mathbf{0.94}$ ($F>0$); $1.05$–$1.21$ at $F<0$, **void** | **the one realisable multi-prime instance**; still $<1$ honestly |
| $\zeta(5)$ | Brown–Zudilin cellular row: $Q_n\in\mathbf Z$, $\kappa_p=0$, no slope at any $p\le47$ | none | none | $p=2$: level-16 rows, $\xi_2=\tfrac7{32}\zeta_2(5)$ — nothing to align with | n/a (counterfactual $\le0.875$) | n/a | deficit: needs $\log\tfrac1\lambda+\sigma_2\log2>k=5$; best available $2.47$ |
| $L(3,\chi_5)$ | none ($\eta$ and AESZ 184 are integral complex-fold rows, $|\lambda_2|=|\lambda_1|$) | none | $\{5\}$, $\sigma_5=3$ | $p=5$: $\eta\leftrightarrow$ AESZ 184, ratio $2{:}1$, exact | $\delta=0$ identically (no decay) | n/a | AESZ 184 has $c=2000=2^4\cdot5^3$ but $\sigma_2=0$ (non-Zagier normalisation eats the $2$-part) |
| $\tfrac1{16}(15L-6\zeta(2))$, $\tfrac1{32}(6\zeta(2)+15L)$ — $\sqrt{\text{Domb}}$ D.5, D.7 | none known | — | modular side: $\{2\}$ only, $\sigma_2=6$; no $3$-adic slope | $p=2$: conductor-6 row, $\xi_2^{\rm D.5}=\tfrac{15}{16}\Lambda_2$, $947$ digits **[verified]** | n/a (periods differ $\Rightarrow$ rank 3) | **impossible**: one prime; the $\zeta(2)$-free combination has $\rho_2=2^{\sigma_2}$, net $\le0$ | §6.5.1–2 |
| $\sqrt{\text{Domb}}$ D.2, D.3, D.4, D.6 (period **unidentified**, not in $\langle1,\zeta(2),L\rangle$) | none exists | — | $\{2,3\}$: $\sigma_2=4$ or $2$, $\sigma_3=1$ | each other: $\xi_p^{(2)}=-2\xi_p^{(3)}$ and $\xi_p^{(4)}=\xi_p^{(6)}$ at both primes ($981/496$, $979/496$ digits) | $0.450$ | $\mathbf{0.548}$ (D.3 $\times$ D.4) | two slope primes, **no decayer** ($|\lambda_2|=4,12$); $F>0$ always, no content | §6.5.3 |
| $\zeta(7)$ | none | — | level-24 row is Apéry-perfect; level-60 has no genus-zero host | — | n/a | n/a | see `ZETA7_LEVEL60.md` |
| $L(f,2)$, $f=\eta_2^3\eta_6^3$ | none known | — | cuspidal clause: $\xi_p=0$ | — | n/a | n/a | `ONE_CLASS_TWO_WORLDS.md` §8: structurally out of reach for the residue mechanism |
| $L(3,\chi_{-3})$ | none | — | $\zeta$-row $(9,3,-27)$: $\sigma_3=3$ but $\xi_3=0$ | — | n/a | n/a | single row only |

**Answer to the headline question: no class gets $\delta>1$ legitimately.** The two classes where a
multi-prime lattice is even conceivable are $L(2,\chi_{-3})$ (where it is now realised, and reaches
$0.94$) and, in principle, $L(2,\chi_{-4})$ (where the $p=3$ partner does not exist). The deficits
are:

| class | deficit |
|---|---|
| $\zeta(3)$ | $\delta=0.9010$; needs $k\le2.2248$ against the sharp $k=3$, i.e. a $26\%$ arithmetic saving, or a Domb $2$-adic slope of $7.68$ instead of $6$ |
| $L(2,\chi_{-4})$ | $\delta=0.9025$; the $G_2$ bridge already takes the Nesterenko pair to the knife-edge $F=0$, and everything past it is void by the $\Theta^*$ control |
| $L(2,\chi_{-3})$ | $\delta\approx0.94$ with two primes; the binding cost is $\eta\approx21$–$25$ per unit $b$ (the conductor-6 rows' own denominators), against a total net $p$-adic resource of at most $(\sigma_2-\kappa_2)\log2+(\sigma_3-\kappa_3)\log3\approx5.5$ |
| $\zeta(5)$, $\zeta(7)$ | no decayer with any $p$-adic slope exists; the $d_n^5$/$d_n^7$ cost is a factor $2$–$3$ beyond anything on offer |
| $L(3,\chi_5)$ | rich $p$-adic side ($3\log5=4.83$ of net resource), empty archimedean side: every known realisation is a complex fold |

---

## 8. What is proved, verified, measured, open

* **[proved]** Theorem $1'$ (§2.1): the conductor-6 well-poised construction and the vanishing of
  the $\zeta(2)$ and $\pi\sqrt3$ coefficients. The multi-prime cancellation identity (§4.1). The
  multi-prime form of the master formula, $G=\sum_p\min(\cdot,\cdot)\log p$
  (`paper/sections/05_two_row.tex`, proof of Theorem E — the Smith-normal-form bound needs only
  $T_n\mid h_n$, nothing about $T_n$ being a prime power).
* **[verified]** (exact arithmetic, stated ranges) the whole census of §1 ($m\le300$ for Zudilin,
  $n\le400$ for the integral rows, $n\le300$ for Brown–Zudilin, $k\le18$ for the conductor-3
  $\kappa_2$); the identity $\sum_{t\ge0}R(t)=QL-P$ for conductor 6; the two-prime alignment of §3
  to $v_3=352$, $v_2=951$ with ten- and eight-scalar controls; integrality of every lattice output
  $q,p$ at every $n$ tested; the surrogate control of §6.
* **[measured]** every rate in §2.2 and every $\delta$ in §§4–5. These are finite-index readings
  with visible $O(\log n/n)$ drift; $\delta$ should be read to two decimals at best. In particular
  $\sigma_2=8$ and $\sigma_3=3$ for the conductor-6 family are conjectural exact values fitted to
  $2.83\le\sigma_3\le2.94$, $7.83\le\sigma_2\le7.97$.
* **[open]**
  1. The recurrence of the conductor-6 family (the conductor-3 $\alpha=2$ member has the exact
     order-2, degree-9 recurrence with roots $-1024,1$; no analogue was fitted here), hence exact
     $\Lambda(\alpha)$, $\lambda(\alpha)$.
  2. A *proved* denominator bound for $P_b$ (the measured $k_d$ is roughly $\tfrac12$ of the
     provable partial-fraction bound $D_{6b+5}^2$, exactly as in the conductor-3 case; closing that
     gap is the single most valuable next step, and Rhin–Viola group actions are the technology).
  3. $\sigma_2=8$, $\sigma_3=3$ as theorems.
  4. The rank-$\ge3$ selection theorem. Everything in §5.3 is heuristic.
  5. Whether the mixed-period conductor-6 family of §6.5.4 (antisymmetry replaced by the two
     conditions $\Sigma B^{(1/6)}=0$, $\Sigma A^{(5/6)}+7\Sigma A^{(1/6)}=0$) contains a
     *decaying* member. The forms exist and are exact; whether any of them is small is the whole
     question, and the antisymmetry that made the §2 family small is precisely what has been
     traded away.
  6. The period of the cusp-move rows D.2, D.3, D.4, D.6 (`lindep` excludes
     $\langle1,\zeta(2),L(2,\chi_{-3})\rangle$ at height $10^{15}$); until it is identified no
     decayer can be built for the one class in the corpus that has two slope primes *and* a
     complete internal alignment.
  7. Whether a conductor-$12$ analogue exists for $L(2,\chi_{-4})$ (poles at
     $\mathbf Z+\tfrac j{12}$, $j\in\{1,5,7,11\}$, giving $\tfrac{10}9G$ with denominators at $2$
     *and* $3$). It would be a two-prime Catalan decayer — but §7 shows there is no modular Catalan
     row with $\sigma_3>0$ to align it with, so it would have nothing to bridge against. Building
     the missing $3$-adic Catalan row is the interesting half of that question.
* **No irrationality is claimed.** $\delta<1$ in every $F>0$ configuration; and even $\delta>1$
  would be void without an unconditional non-vanishing statement for the *selected* linear form
  (`CATALAN_AUDIT.md` §6).

---

## 9. Corrections this note forces on earlier documents

1. `paper/sections/05_two_row.tex`, Remark (multi-prime bridging): "*Multi-prime bridging is thus
   available in principle and realizable nowhere in the canonical classification — a clean
   negative, not an oversight*". The first clause stands and is confirmed (§1.2); the conclusion
   must be qualified: **it is realisable outside the canonical classification**, by the conductor-6
   hypergeometric decayer of §2, which is not a modular row and therefore not in the list of twelve
   $c$-values. The Remark's ceiling argument ("the resource is still at most
   $\gamma\log\Lambda_{\rm dec}$") is also an integral-row statement and is superseded by §4.1 for
   non-integral engines.
2. `ONE_CLASS_TWO_WORLDS.md` §3 records "*no other prime has a slope*" for the conductor-3
   $\alpha=2$ member and reports the aggregate prime-to-3 rate $\nu$ for the others. Both are
   correct, but the decomposition matters: the $\alpha\ne2$ members carry **$2$-power denominators
   at linear rate** ($\kappa_2$ up to $5.6$), and it is worth recording that this buys nothing
   ($\sigma_2=0$). §1.3.
3. `CUSP_MOVE_PROGRAM.md` §6, the row recording
   $\xi_p^{(2)}=-2\,\xi_p^{(3)}$ with $792$ digits at $p=2$ **and** $386$ at $p=3$: reproduced
   exactly here ($981$ and $496$ at $N=500$). Two readings should be attached to it. (i) The rows
   in question are D.2 and D.3, whose trailing coefficients carry linear terms; the two primes are
   *not* $v_p(c)$, and the same statement is **false** for the mixed-period members D.5, D.7 of the
   same orbit, which have no $3$-adic slope at all ($v_3$ of the increments is $0$ or $-2$ out to
   $n=500$). (ii) D.2/D.3's period is not in $\langle1,\zeta(2),L(2,\chi_{-3})\rangle$, so the
   $(-6,15)$/$(6,15)$ coefficient direction belongs to the *other* half of the orbit. §6.5.
4. `ZETA3_TWO_LATTICE.md` §13.2 ("**F** is the unique row with two slope primes … multi-prime
   bridging is realisable nowhere") — same qualification as (1). The observation that the per-prime
   term is a minimum, so a pair containing $\mathbf F$ loses its $p=2$ contribution, is exactly
   right and is now explained by the cancellation identity of §4.1: even when the $p=2$ term
   survives, it is worth zero. The clause "$\mathbf F$ is the unique row with two slope primes"
   is true of the twelve *sporadic* families but not of the corpus: the cusp-move rows
   D.2, D.3, D.4, D.6 also have two (§6.5.3).

---

## 10. Reproduction

All scripts are PARI/GP, exact rational and $p$-adic arithmetic, absolute paths, and each begins
with `default(parisizemax, 8000000000)`. **`lib.gp` must be `read()` and must *not* itself contain
a `default(parisizemax,...)` call** — a stack reallocation during `read` silently aborts the rest
of the parse (this cost an hour; it is why the `default` sits in every driver instead).

```sh
cd /home/ubuntu/code/math-modular-sources/lattice/multi_prime
mkdir -p out

# 1. census: which decayers have how many alignment primes
timeout 2400 gp -q 01_kappa_scan.gp     </dev/null > out/01_kappa_scan.log 2>&1     #  ~2 min
timeout 2400 gp -q 02_decayer_padic.gp  </dev/null > out/02_decayer_padic.log 2>&1  #  ~4 min
timeout 2400 gp -q census_1_zud.gp      </dev/null > out/census_1_zud.log 2>&1
timeout 2400 gp -q census_2_nest.gp     </dev/null > out/census_2_nest.log 2>&1
timeout 2400 gp -q census_3_introws.gp  </dev/null > out/census_3_introws.log 2>&1
timeout 2400 gp -q census_6_bz.gp       </dev/null > out/census_6_bz.log 2>&1
timeout 2400 gp -q census_7_summary.gp  </dev/null > out/census_7_summary.log 2>&1

# 2. the conductor-6 row: slopes at both primes, full rate profile
timeout 2400 gp -q 03_chi6_scan.gp      </dev/null > out/03_chi6_scan.log 2>&1      #  ~3 min
timeout 2400 gp -q 04_chi6_rates.gp     </dev/null > out/04_chi6_rates.log 2>&1     #  ~5 min

# 3. master formula, single vs multi prime
timeout 2400 gp -q 05_design.gp         </dev/null > out/05_design.log 2>&1         #  ~4 min
timeout 3000 gp -q 06_design_full.gp    </dev/null > out/06_design_full.log 2>&1    # ~12 min

# 4. the actual lattice: zeta(3) control, two-row multi-prime, three-row
timeout 2400 gp -q 07_lattice.gp        </dev/null > out/07_lattice.log 2>&1        #  ~8 min
timeout 2400 gp -q 08_best.gp           </dev/null > out/08_best.log 2>&1           # ~20 min
timeout 2400 gp -q 09_control.gp        </dev/null > out/09_control.log 2>&1        #  ~6 min

# 5. the mixed zeta(2)/L(2,chi_-3) class of the sqrt(Domb) cusp-move orbit
timeout 2400 gp -q 10_sqrtdomb.gp       </dev/null > out/10_sqrtdomb.log 2>&1       #  ~2 min
timeout 2400 gp -q 11_mixed_class.gp    </dev/null > out/11_mixed_class.log 2>&1    #  ~3 min
timeout  600 gp -q 12_mixed_delta.gp    </dev/null > out/12_mixed_delta.log 2>&1    #  seconds
```

Nothing here needs **snake**; the longest single run is about 20 minutes on the dev box. If it is
sent there anyway (`tools/SNAKE.md`), run each script detached with one log per job and pull back
`lattice/multi_prime/out/`.

Key single-line checks:

```sh
# the two-prime alignment certificate (should print 352 / 352 / 951)
grep 'alpha=2/1  b=120' out/08_best.log

# the surrogate control (the two delta columns must be identical)
grep 'c6_1/3 x c6_1/4 \[BOTH\]' out/09_control.log

# the zeta(3) pipeline calibration (d_emp must fall towards 0.9010)
grep 'Domb(2n) x T(3n)' out/07_lattice.log
```

### File map

| file | what |
|---|---|
| `lib.gp` | shared exact row library: `row2`, `row3`, `cooper18`, `zudrow`, `chi3row`, and the new `chi6row`, `cfrow` |
| `01_kappa_scan.gp` | factors $\operatorname{den}(Q_n)$ for every decayer; finds all primes with $\kappa_p>0$ |
| `02_decayer_padic.gp` | conductor-3 family: $p$-adic limits at $2$ and $3$ against $\mathbf C$, $\mathbf F$ — the $\sigma_2=0$ negative |
| `03_chi6_scan.gp` | conductor-6 family: first scan, both slopes |
| `04_chi6_rates.gp` | full rate profile $(\Lambda,\lambda,\kappa_2,\kappa_3,\nu,\eta,\sigma_2,\sigma_3)$ + scalar controls |
| `05_design.gp`, `06_design_full.gp` | master formula, no bridge / $p=3$ / $p=2$ / both, over all pairs, exact `lcm` denominators |
| `07_lattice.gp` | the real congruence lattice: $\zeta(3)$ control, two-row multi-prime, $\times\mathbf F$, three-row |
| `08_best.gp` | the $\ge200$-digit two-prime alignment certificate; empirical sweep over all conductor-6 pairs |
| `09_control.gp` | the rational-$\Theta^*$ control on every $F<0$ reading |
| `10_sqrtdomb.gp` | mixed-period rows D.5, D.7: limits, slopes, and their $2$-adic alignment with the conductor-6 decayer |
| `11_mixed_class.gp` | two-slope-prime rows D.2, D.3, D.4, D.6: slopes at both primes, `lindep` on the period, alignment scan |
| `12_mixed_delta.gp` | master formula on the D.2/D.3/D.4/D.6 pairs, none / $p{=}3$ / $p{=}2$ / both |
| `census_*.gp` | independent rebuild and slope scan of Zudilin, Nesterenko, the twelve sporadics, Brown–Zudilin |
