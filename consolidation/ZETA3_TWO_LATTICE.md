# The two-lattice / two-row 2-adic method applied to $\zeta(3)$ (Domb + T)

*Claude (Fable), 2026-08-21. Scripts: `lattice/zeta3_lattice/*.gp` (PARI, exact integer model).*

**Headline.** Running the Catalan two-row construction verbatim on the aligned
pair Domb $(10,4,64)$ + T $(12,4,16)$ gives construction quality
$$
\boxed{\ \delta \;=\; \frac{12\log(12+8\sqrt2)-24\log 2}{9+9\log(12+8\sqrt2)-20\log 2}
\;=\;0.9009531686558563\ldots\ }
$$
which is **$<1$**. So this is **not** a second proof of Apéry's theorem.
It is a clean calibration: the method behaves on $\zeta(3)$ almost exactly as it does
on Catalan ($0.90253$), landing $\approx 10\%$ short in both cases.
The near-coincidence $0.9010$ vs $0.9025$ is, as far as I can tell, an accident of
the numbers, but it is a useful sanity check that the machine was transcribed correctly.

---

## 0. Setup and notation

Zagier/Almkvist–Zudilin normalisation
$$(n+1)^3u_{n+1}=(2n+1)(an^2+an+b)u_n-c\,n^3u_{n-1},\qquad
a_0=1,\ a_1=b;\qquad b_0=0,\ b_1=1 .$$

| row | $(a,b,c)$ | $a_n\sim$ | limit | linear form | $\sigma_2=v_2(c)$ |
|---|---|---|---|---|---|
| **D** (Domb) | $(10,4,64)$ | $16^n$ | $b_n/a_n\to\tfrac7{24}\zeta(3)$ | $\bigl|\tfrac7{24}\zeta(3)a_n-b_n\bigr|\sim4^n$ (**grows**) | $6$ |
| **T** | $(12,4,16)$ | $\Lambda^n$, $\Lambda=12+8\sqrt2$ | $b_n/a_n\to\tfrac7{32}\zeta(3)$ | $\sim\lambda_2^n$, $\lambda_2=12-8\sqrt2=0.6863$ (**decays**) | $4$ |

$\Lambda\lambda_2=16$, $\log\Lambda=3.1490415\ldots$, $\log\lambda_2=-0.3764528\ldots$

**Integer model used throughout.** $a_n\in\mathbf Z$ and $b_n=c_n/(n!)^3$ with
$$c_{n+1}=(2n+1)(an^2+an+b)c_n-c\,n^6c_{n-1},\qquad c_0=0,\ c_1=1,$$
so all computations are pure integer arithmetic ( `lattice/zeta3_lattice/rows.gp` ).
Rows to $n=800$ take milliseconds.

**New identification (verified).** The T row is
$$a^T_n=\sum_k\binom nk^2\binom{2k}n^2 \qquad(1,4,40,544,8536,145504,\dots),$$
checked to $n=7$; Domb is $a^D_n=\sum_k\binom nk^2\binom{2k}k\binom{2n-2k}{n-k}$.

---

## 1. Task A — exact 2-adic arithmetic

All statements below are **verified** by exhaustive computation in the stated range,
not proved, unless marked *proved*.

### 1.1 Denominators

*Verified* ($n\le400$): $a^D_n,a^T_n\in\mathbf Z$; $d_n^3b_n\in\mathbf Z$ for both rows
($d_n=\operatorname{lcm}(1,\dots,n)$), and $k=3$ is **sharp**: $d_n^2b_n\notin\mathbf Z$.
Moreover the *true* denominator is essentially $d_n^3$ — the measured
$k_{\rm eff}(n)=\log\operatorname{den}(b_n)/n$ is
$2.94,2.94$ at $n=400$ against $\log d_n^3/n=2.98$, i.e.
$k_{\rm eff}\to3$. **There is no free denominator saving** (`dens.gp`).
No extra 2-power denominators beyond those inside $d_n^3$.

### 1.2 $2$-adic valuations of the $a_n$

*Verified for all $1\le n\le700$, no exceptions:*
$$\boxed{\,v_2(a^T_n)=3s_2(n)-[n\ \mathrm{odd}]\,}$$
(the exact analogue of the Catalan lemma $v_2(A_n)=2s_2(n)$).
For Domb only a bound: $v_2(a^D_n)-3s_2(n)\in[-5,9]$ for $n\le700$, so
$v_2(a^D_n)\le 3s_2(n)+9=O(\log n)$; $\max_{n\le700}v_2(a^D_n)=27$.

### 1.3 Exact tail valuations

*Proved* from the exact Casoratian $a_nb_{n-1}-a_{n-1}b_n=-c^{\,n-1}/n^3$:
$$v_2\Bigl(\frac{b_n}{a_n}-\frac{b_{n-1}}{a_{n-1}}\Bigr)
=(n-1)v_2(c)-3v_2(n)-v_2(a_n)-v_2(a_{n-1})=:\iota(n),$$
hence $b_n/a_n$ is $2$-adically Cauchy and
$v_2(\xi_2-b_m/a_m)=\min_{j>m}\iota(j)=:\tau(m)$.

For **T** everything collapses to a closed form (*verified $2\le n\le699$, no exceptions*):
$$\boxed{\ \iota_T(n)=4n-6-6\,s_2(n-1),\qquad
\tau_T(m)=\min_{j>m}\bigl(4j-6-6s_2(j-1)\bigr)\ }$$
so $\tau_T(m)=4m-O(\log m)$ with measured deficit $4m-\tau_T(m)\in[8,50]$ for $m\le400$.
For **D**: $\iota_D(n)=6(n-1)-3v_2(n)-v_2(a^D_n)-v_2(a^D_{n-1})$,
$6m-\tau_D(m)\in[6,52]$ for $m\le400$.
Note $\iota$ is **not** monotone (unlike the Catalan case), so the "first increment"
shortcut is wrong; one must take the min.

### 1.4 The cross determinant

$$\Delta_{m,n}:=3\,a^T_n\,b^D_m-4\,a^D_m\,b^T_n .$$
Writing $\xi^D_m=b^D_m/a^D_m$ etc. and $\varepsilon:=3\xi^D_2-4\xi^T_2\in\mathbf Q_2$,
$$\Delta_{m,n}=a^T_na^D_m\Bigl[3(\xi^D_m-\xi^D_2)-4(\xi^T_n-\xi^T_2)+\varepsilon\Bigr].$$
**If $\varepsilon=0$** the ultrametric inequality gives (*proved conditionally on $\varepsilon=0$*)
$$\boxed{\ v_2(\Delta_{m,n})\ \ge\ v_2(a^T_n)+v_2(a^D_m)+\min\bigl\{\tau_D(m),\,2+\tau_T(n)\bigr\}\ }$$
with equality whenever the two entries differ.
*Verified exhaustively for $1\le m,n\le60$*: **0 violations out of 3600**, and exact
equality in 3443/3600 cases (the 157 exceptions are ties, where $v_2$ is larger).

The two affine slopes $6m$ and $4n$ balance at $6m=4n$, i.e.
$$\boxed{\,m:n=2:3\,}\qquad\text{(Domb at }2n,\ \text{T at }3n).$$
(Confirmed by scanning $v_2(\Delta_{an,cn})/n$ over $a\le4$, $c\le6$: the values are
$\min(6a,4c)+O(1/n)$; $2{:}3\to11.75$, $3{:}5\to17.875$, $4{:}6\to23.75$ at $n=40$.)
Then
$$v_2(\Delta_{2n,3n})=12n-O(\log n),\qquad
12n-v_2(\Delta_{2n,3n})\in[-5,17]\ \text{for }n\le110,$$
with the deficit bounded by $5\log_2 n$ in that range.
Diagonal check: $v_2(\Delta_{n,n})-(4n-3)\in\{-1,0,1\}$ (so the previously reported
"exactly $4n-3$" is only true on the sampled multiples of $50$).

**Comparison with Catalan.** There $v_2(\Delta_{5n,8n})\ge20n$ *and* an extra
$2^{e_{5n}}=2^{20n}$ came from $v_2(Q_m)=-4m+2s_2(m)<0$, total $40n\log2$.
Here both $a_n$ are integers with $v_2=O(\log n)$, so there is **no second factor**:
the total 2-adic divisor is only $12n\log2$. That is the arithmetic reason the
$\zeta(3)$ configuration is weaker than it looks.

---

## 2. Task B — the lattice argument and the quality number

### 2.1 The construction (verbatim Catalan §3, with $k=3$)

Sample Domb at $m=2n$, T at $k=3n$. Put $S_n=d_{3n}^3$ ($\log S_n\sim9n$) and
$$(X_{1,n},Y_{1,n})=S_n\bigl(7a^D_{2n},\,24\,b^D_{2n}\bigr),\qquad
(X_{2,n},Y_{2,n})=S_n\bigl(7a^T_{3n},\,32\,b^T_{3n}\bigr)\in\mathbf Z^2 .$$
(The rational factors are chosen so that $X_r\zeta(3)-Y_r$ *is* the linear form:
$\tfrac7{24}$ for Domb, $\tfrac7{32}$ for T.) Removing the common $S_n$ from the
first coordinates, the mixed minor is
$$h_n=\alpha_{1,n}Y_{2,n}-\alpha_{2,n}Y_{1,n}=S_n\cdot 56\,\Delta_{2n,3n},$$
since $7\cdot24\cdot a^Tb^D-7\cdot32\cdot a^Db^T=56\,\Delta$. Hence
$T_n=2^{12n-O(\log n)}\mid h_n$, $\ \tfrac1n\log T_n\to12\log2$.
Correlated congruence lattice $K_n$, index $\le M_n=S_nT_n$, sublattice $L_n$ with
$M_n\le\operatorname{covol}(L_n)<2M_n$:
$$\sigma=9+12\log2=17.3177661667\ldots$$

### 2.2 Exponential rates

$$\begin{aligned}
A_1&=9+8\log2, & E_1&=9+4\log2 &&\text{(Domb, }16^{2n}\text{ and }4^{2n}),\\
A_2&=9+3\log\Lambda, & E_2&=9+12\log2-3\log\Lambda &&\text{(T)} .
\end{aligned}$$
Crossed-gap condition: $(A_2+E_1)-(A_1+E_2)=6\log\Lambda-16\log2=7.8039>0$ ✓
(Catalan: $50\log\varphi-8\log2$).

Two-successive-minima selection theorem with anisotropy
$x=\tfrac12(\sigma+E_2-E_1)=\tfrac12(9+20\log2-3\log\Lambda)=6.70690$:
$$H=A_2-x=\tfrac12\bigl(9+9\log\Lambda-20\log2\bigr)=11.7392151026,$$
$$F=x+E_1-\sigma=\tfrac12(9+20\log2-3\log\Lambda)-8\log2=1.1627320584,$$
$$H-F=A_2-E_2=6\log\Lambda-12\log2=10.5764830442 .$$

### 2.3 Quality

$$\delta=\frac{H-F}{H}=\frac{12\log\Lambda-24\log2}{9+9\log\Lambda-20\log2}
=0.9009531686558563\ldots<1 .$$

$F>0$: the linear forms **grow**, so there is no irrationality statement at all
(not even a weak one). No irrationality measure follows. For reference, had
$\delta>1$ held one would get $\mu(\zeta(3))\le\delta/(\delta-1)$ — e.g. $\delta=1.033$
would give $\mu\le31.4$, already worse than Apéry's $13.42$ and far from
Rhin–Viola's $5.5139$. **Even in the best case this architecture would not be
competitive as a measure; its only possible value would be as a new proof.**

### 2.4 Optimality of $2:3$ and the archimedean-only comparison

General sampling Domb at $an$, T at $cn$, denominator exponent $k$:
$$\sigma=k\max(a,c)+\min(6a,4c)\log2,\quad
A_1=k\max+4a\log2,\ E_1=k\max+2a\log2,$$
$$A_2=k\max+c\log\Lambda,\ E_2=k\max+4c\log2-c\log\Lambda,\quad
\delta(a,c)=\frac{2c(\log\Lambda-2\log2)}{A_2-\tfrac12(\sigma+E_2-E_1)} .$$
Numerical optimisation over $r=c/a\in(0.4,100)$ (`quality.gp`) gives a **unique
maximum at $r=3/2$** — the same ratio that balances the two 2-adic tails, exactly
as $5{:}8$ did for Catalan.

**Archimedean-only version** (drop the hidden 2-adic divisor, $T_n=1$):
$\delta_{\rm arch}(2{:}3)=0.6652671861$, and $\sup_r\delta_{\rm arch}=0.7278$
(approached as $r\to\infty$, not attained).
So the hidden 2-adic determinant is worth
$$0.9009532-0.6652672=+0.2356860$$
in quality. The mechanism is real and large; it is simply not large enough.

### 2.5 How far short is it? (two exact thresholds)

*Threshold in the denominator exponent.* With everything else fixed,
$\delta\ge1\iff$
$$k\ \le\ k^*=\log(12+8\sqrt2)-\tfrac43\log2=2.2248452944\ldots$$
Actual $k=3$, and §1.1 shows $k=3$ is sharp for these rows. A Rhin–Viola/Zudilin
"arithmetic method" saving would have to remove $26\%$ of the denominator to close
the gap — much more than such savings usually give for a modular row.

*Threshold in the 2-adic divisor.* With $k=3$, one needs
$\tfrac1n\log T_n\ge 9+16\log2-3\log\Lambda=10.6432$ instead of $12\log2=8.3178$,
i.e. an extra factor $2^{3.3549\,n}$ — equivalently a Domb 2-adic slope of $7.68$
instead of $6$. Nothing in sight supplies this.

*Robustness check performed.* I also ran the asymmetric variant in which each row is
cleared by its own denominator ($d_{2n}^3$ and $d_{3n}^3$ separately). All four rates
and $\sigma$ shift, but $H$ and $F$ come out **identical**; $\delta$ is unchanged.
The LCM-variant architecture of `CATALAN_TWO_ROW_FULL_LCM_PAPER.tex` gave a *lower*
quality for Catalan ($0.6588$ vs $0.9025$) and there is no reason to expect otherwise here.

---

## 3. Task C — what is proved, what is verified, what is open

### Proved (unconditionally)
1. The Casoratian identity and hence the exact increment valuation $\iota(n)$;
   $2$-adic convergence of $b_n/a_n$ for both rows with slopes $6$ and $4$,
   *given* $v_2(a_n)=O(\log n)$.
2. The reduction $\Delta_{m,n}=a^T_na^D_m[3(\xi^D_m-\xi^D_2)-4(\xi^T_n-\xi^T_2)+\varepsilon]$.
3. The lattice/selection machinery of the Catalan paper, once the divisor $T_n\mid h_n$
   is granted, yields $H,F$ as computed.

### Verified numerically (large range, no exceptions) — but not proved
4. $v_2(a^T_n)=3s_2(n)-[n\text{ odd}]$ and $\iota_T(n)=4n-6-6s_2(n-1)$, $n\le700$.
5. $v_2(a^D_n)\le3s_2(n)+9$, $n\le700$ (any $O(\log n)$ bound suffices).
6. $d_n^3b_n\in\mathbf Z$, $n\le400$, with $k=3$ sharp.
7. $\varepsilon=0$, i.e. $3\xi^D_2=4\xi^T_2$: from $v_2(\Delta_{330,330})=1314$ and
   $v_2(a^T_{330})+v_2(a^D_{330})\le20$ one gets $v_2(\varepsilon)\ge1294$.

### Open — and the key honest point

> **The task brief asked whether the determinant lower bound can be obtained
> *without* identifying the common 2-adic limit. It cannot.**
>
> If $\varepsilon\ne0$ then for all large $m,n$
> $v_2(\Delta_{m,n})=v_2(a^T_n)+v_2(a^D_m)+v_2(\varepsilon)=O(\log n)$ — bounded —
> and the whole construction collapses. The Casoratians control the two *tails*
> unconditionally, but they say nothing about the constant $\varepsilon$: two
> $2$-adically convergent sequences can have limits differing by a $2$-adic number
> of huge valuation. No finite computation can distinguish $\varepsilon=0$ from
> $v_2(\varepsilon)=10^6$. So the equality of the two $2$-adic limits **is** the gap,
> exactly as in the Catalan paper. (This corrects the optimistic reading in
> `THEORY_NOTES_03` §3 item 3.)

**Cleanest route to closing it.** Both rows are modular of the
$F\cdot D_q^{-3}\Phi$ type with $\Phi$ Eisenstein of weight 4
(book `v8/08_applications.tex`, Application B; `EISENSTEIN_SOURCE_THEOREM.tex`):
$$B(t(q))=F(q)\,D_q^{-3}\Phi(q),\qquad D_q=q\frac{d}{dq},\qquad \Phi=F\cdot D_qt\ \text{(wt 4)} .$$
Domb is the level-$6$ system ($t=(\eta_1\eta_6/\eta_2\eta_3)^{12}$-type Hauptmodul).
For T, $c=16=2^4$ forces the level to be $2$-power; **this is now RESOLVED in §9: T is the level-8 system
with $t=(\eta_1\eta_8/\eta_2\eta_4)^8$ and $F$ a weight-2 Eisenstein series on $\Gamma_0(8)$.**
The original failed attempt: an exhaustive search over
$t=q\prod_{d\mid8}(\eta(d\tau)q^{-d/24})^{r_d}$ with $\operatorname{ord}_\infty t=1$,
$|r_d|\le16$, found no $F=\sum a^T_nt^n$ that is an eta quotient on $\{1,2,4,8\}$
(`modular_id.gp`); the level-$16$ search timed out. **Open sub-task.**

Given the parametrisations, the proof route is:
$D_q^{-3}$ acting on a weight-4 Eisenstein series is a $2$-adic Eichler integral;
its constant term after Teichmüller normalisation should be a fixed multiple of
Calegari's $\zeta_2(3)$ in both cases, giving $\xi^D_2=\tfrac7{24}\zeta_2(3)$ and
$\xi^T_2=\tfrac7{32}\zeta_2(3)$, whence $3\xi^D_2=4\xi^T_2$. This is the exact
$p=2$, weight-4 analogue of the sketch in the Catalan paper and is, I believe,
genuinely reachable — but it is **not done here**, and note the irony: since
$\delta<1$, proving it would gain nothing.

### Ancillary open items
- Prove $v_2(a^T_n)=3s_2(n)-[n\text{ odd}]$ (a Kummer/Lucas argument on
  $\sum_k\binom nk^2\binom{2k}n^2$ should do it).
- Prove any $O(\log n)$ bound for $v_2(a^D_n)$.
- Prove $d_n^3b_n\in\mathbf Z$ (routine from the $D_q^{-3}$ representation once the
  parametrisation is fixed — the Catalan paper's "two integrations" argument with
  three integrations).

---

## 4. Verdict

The Domb+T configuration is a genuine Catalan-type two-lattice configuration for
$\zeta(3)$: aligned $2$-adic limits, complementary slopes $6$ and $4$, decaying
partner, optimal sampling $2{:}3$ forced by the arithmetic, and a real $+0.236$
quality gain from the hidden determinant. It nevertheless **misses irrationality
by $10\%$**, for a structural reason that is now explicit: the third-order rows pay
$d_n^3$ ($k=3$) while supplying only $\min(6a,4c)\log2=12\log2$ of $2$-adic divisor,
against Catalan's $k=2$ with $40\log2$. The deficit is not a near miss that better
bookkeeping will close; the two thresholds in §2.5 quantify exactly how far.

The useful conclusion for the programme: **the two-lattice method's ceiling is set by
$k\cdot\max(a,c)$ versus $\min(\sigma_2 a,\sigma'_2c)\log2$.** Weight-2 (second-order,
$k=2$) rows are where it can work; weight-3 (third-order, $k=3$) rows are structurally
handicapped. Searching for a *second-order* pair aligned at some prime, with a decaying
partner, is a better use of effort than pushing this one.

---

## 5. Addendum — the master formula for the whole architecture

Doing the $\zeta(3)$ case forced the bookkeeping into closed form, and the result is
worth more than the $\zeta(3)$ number itself. Let a two-row Catalan-type construction
have: **decayer** sampled at $\gamma n$ with $a_n\sim\Lambda^{\gamma n}$ and linear form
$\sim\lambda^{\gamma n}$ ($\lambda<1$), $p$-adic slope $\sigma_{\rm dec}$; **engine**
sampled at $\alpha n$ with linear form $\sim\rho_2^{\alpha n}$, slope $\sigma_{\rm eng}$;
common denominator $S=k\cdot(\text{max denominator index})$; any extra integralising
multiplier $\eta$ (rate per $n$); hidden divisor rate
$G=\min(\sigma_{\rm eng}\alpha,\ \sigma_{\rm dec}\gamma)\log p$. Then

$$\boxed{\ F=\tfrac12\Bigl[\underbrace{S+\eta+\alpha\log\rho_2}_{\text{COST}}
-\underbrace{\bigl(\gamma\log\tfrac1\lambda+G\bigr)}_{\text{RESOURCE}}\Bigr],\qquad
H=F+\gamma\log\tfrac\Lambda\lambda,\qquad \delta=1-\frac FH\ }$$

**$F$ is exactly one half of the cost–resource deficit.** Verified against both
instances to full precision:

| | $F$ | $H$ | $\delta$ |
|---|---|---|---|
| Catalan ($k=2$, $\eta=20\log2$, $\alpha{:}\gamma=8{:}5$) | $2.5985578256$ | $26.6591490786$ | $0.9025266029$ ✓ paper |
| $\zeta(3)$ ($k=3$, $\eta=0$, $\alpha{:}\gamma=2{:}3$) | $1.1627320584$ | $11.7392151026$ | $0.9009531687$ |

### 5.1 The design rule

For **integral** rows ($\eta=0$) at the balance point $\alpha\sigma_{\rm eng}=\gamma\sigma_{\rm dec}$,
the product formula $\Lambda\lambda=c$, $\log c=\sum_p\sigma_p\log p$ makes
$\gamma\log\frac1\lambda+G=\gamma\log\Lambda$ **exactly** (when $p$ is the only prime in $c$).
So the entire criterion collapses to

$$\boxed{\ \delta>1\iff \log\Lambda_{\rm dec}\;>\;k\max(r,1)+r\log\rho_2^{\rm eng},
\qquad r=\frac{\alpha}{\gamma}=\frac{\sigma_{\rm dec}}{\sigma_{\rm eng}}\ }$$

This is the "budget" invariant of `THEORY_NOTES_03` §1 vindicated in exact form:
**what the decaying row contributes is its numerator growth $\log\Lambda$, not its
archimedean decay** — the $p$-adic slope makes up the whole difference, automatically.
The method is a machine for converting the product formula into Diophantine decay.

Shopping list it produces:

| $k$ | $\sigma_{\rm eng}$ | $\sigma_{\rm dec}$ | $\log\rho_2^{\rm eng}$ | need $\Lambda_{\rm dec}>$ |
|---|---|---|---|---|
| 3 | 6 | 4 | $\log4$ | **50.6** (we have $23.31$ — the $\zeta(3)$ failure) |
| 3 | 6 | 6 | $\log4$ | 80.3 |
| 2 | 5 | 8 | $\log4$ | 225.4 (Catalan; non-integral, see $\eta$) |
| 2 | 5 | 5 | $\log4$ | 29.6 |
| 2 | 5 | 4 | $\log4$ | 22.4 |
| 2 | 3 | 3 | $\log2$ | **14.8** |

The last line is the interesting one: a $k=2$ pair with modest slopes and a *slowly
growing engine* ($\rho_2=2$) only needs $\Lambda_{\rm dec}>14.8$ — Zagier's row
D $(11,3,-1)$ has $\Lambda=11.09$, so we are within a factor $1.34$. **Lowering
$\log\rho_2^{\rm eng}$ is the cheapest lever in the whole formula** and it was never
optimised: every construction so far grabbed whatever engine was aligned, not the one
with the slowest-growing linear form.

---

## 6. What I think we should do next (ranked)

1. **Search on the design rule, not on alignment.** The criterion above is three
   numbers per candidate pair $(\log\Lambda_{\rm dec},\ \sigma$'s$,\ \log\rho_2^{\rm eng})$.
   Every second-order Apéry-like row in the literature (Zagier's six, Cooper's
   $s_7,s_{10},s_{18}$, the Beukers/Rivoal Padé rows, Rhin–Viola integrals, the
   Almkvist–Zudilin–Straub tables) can be scored in seconds. Do the whole census and
   sort by $\log\Lambda_{\rm dec}-k\max(r,1)-r\log\rho_2$. This is a few hours of work
   and it either finds a live candidate or closes the architecture for good. **Do this first.**
2. **Minimise the engine's $\rho_2$.** The engine is only there to supply $p$-adic slope
   and a second lattice direction; its linear form growing at $\rho_2^{\alpha n}$ is pure
   cost. Look for aligned rows with $\rho_2$ close to $1$ — i.e. rows whose two
   characteristic roots are nearly equal in the *other* direction, or non-integral rows
   where $\eta$ can be made to cancel.
3. **Rank $\ge3$ lattices.** With $r$ aligned rows the covolume is still $\approx S\cdot T$
   (two congruences, $r$ variables) while Minkowski gives $r$ successive minima of size
   $\mathrm{covol}^{1/r}$ instead of $\mathrm{covol}^{1/2}$. That is a genuine structural
   change, not a constant-factor one, and nobody has done the bookkeeping. The
   $L(2,\chi_{-3})$ triple $(B,C,F)$ aligned at $p=3$ is the natural test bed even though
   no member decays. **Cheapest way to settle it: skip the theory and run LLL on the actual
   integer rows for $n=20,\dots,60$**, measure $\log|q|$ and $\log|q\Theta-p|$ empirically,
   fit the slopes. That measures the true $\delta$ of any architecture in an afternoon and
   would have caught my $0.9010$ without any of §2.
4. **Finish the T-row modular identification** (level, Hauptmodul, weight-4 $\Phi$). It is
   needed for §7 below and it is the one loose computational end here.

---

## 7. The integral variation of mixed Hodge structure

I think this is the right home for the whole "alignment law", and it explains every
empirical fact in `THEORY_NOTES_03` §2 including the negative ones. Here is the picture
stated so that it can be attacked, plus what it predicts.

### 7.1 The object

For a modular Apéry-like row parametrised by $(t,F)$ on a genus-zero
$X_0(N)$, let $\mathcal H$ be the rank-2 VHS of the universal elliptic curve and
$\mathcal V=\mathrm{Sym}^{w}\mathcal H$ ($w=1$ for the second-order rows, $w=2$ for the
third-order ones). The companion construction $B=F\cdot D_q^{-(w+1)}\Phi$ with $\Phi$
Eisenstein of weight $w+2$ is an **Eichler integral**, and an Eichler integral is
precisely a section of an extension of variations of mixed Hodge structure

$$0\longrightarrow \mathcal V\longrightarrow \mathcal E\longrightarrow \mathbf Z(0)\longrightarrow 0,$$

i.e. a class $[\mathcal E]\in\mathrm{Ext}^1_{\mathrm{VMHS}/\mathbf Z}(\mathbf Z(0),\mathcal V)$.
The Eisenstein-ness of $\Phi$ is what makes the class **Beilinson–Kings**: it is the
motivic Eisenstein class, defined over $\mathbf Z$ (this is exactly the content of
`EISENSTEIN_SOURCE_THEOREM.tex` read backwards).

Now the two limits are *two realisations of one class*:

- **Archimedean (Betti/de Rham) realisation** at the point $t\to$ the singular fibre:
  the period of the extension is $\lim b_n/a_n$. That is the Apéry limit, $r\Theta$.
- **$p$-adic (syntomic / rigid) realisation**: the same class in
  $\mathrm{Ext}^1_{\text{syn}}$, whose period is $\xi_p$. That is the $p$-adic Apéry
  limit, and the existence of a slope $\sigma_p=v_p(c)>0$ is exactly the statement that
  the singular fibre is $p$-adically non-unit, i.e. the class is *overconvergent* there
  (Calegari's mechanism).

**Integrality is the whole point.** $\mathrm{Ext}^1$ in the category of *integral*
VMHS is a finitely generated $\mathbf Z$-module. Two rows carrying the same period
$\Theta$ give two classes in the same rank-one $\mathbf Z$-module, hence
$$q\,[\mathcal E_1]=q'\,[\mathcal E_2]\quad\text{in }\mathrm{Ext}^1\text{ for some }q,q'\in\mathbf Z .$$
A relation in $\mathrm{Ext}^1$ is realisation-independent. Therefore

> **Conjecture (realisation rigidity).** If two Apéry-like rows have archimedean limits
> $r\Theta$ and $r'\Theta$ for the same period $\Theta$, then for **every** prime $p$ at
> which both have positive slope, $\xi_p=r\Theta_p$ and $\xi'_p=r'\Theta_p$ for one
> $\Theta_p\in\mathbf Q_p$. In particular $r'\xi_p=r\xi'_p$ — *the same rational
> relation as in the archimedean realisation.*

This is not a re-labelling of the empirical law: it says the *specific rational number*
$r/r'$ transfers, and it says *why* rows with different periods never align (their
classes lie in different $\mathrm{Ext}^1$'s, so no relation exists to transfer) and why
row B aligns 3-adically with C despite having **no archimedean limit at all** (complex
conjugate roots kill the Betti period, but the class — and hence the syntomic period —
survives; the archimedean realisation is the degenerate one, not the $p$-adic one).

### 7.2 It is confirmed, sharply, by this project's data

- $\zeta(3)$: $r_D=\tfrac7{24}$, $r_T=\tfrac7{32}$, so $r_D/r_T=4/3$. The predicted
  and verified $2$-adic relation is $3\xi^D_2=4\xi^T_2$ — **the same $4:3$**, to
  $v_2\ge1294$ (§3.7). The conjecture got the exact rational number right.
- $L(2,\chi_{-3})$: $r_C=\tfrac12$, $r_F=\tfrac58$; the observed relation is
  $4\xi^F_3=5\xi^C_3$ — again the same ratio.
- Negative instances (A–E, A–F, E–F at $p=2$): different periods, no relation, none found.

Six for six, including three that could have falsified it.

### 7.3 A theorem one could actually prove from this

The conjecture has a form that is provable in the Beilinson–Kings framework, because
both sides are *specialisations of the motivic Eisenstein class at torsion sections*.
Concretely, for the $\zeta(3)$ pair the statement to prove is

> **Proposition (target).** Let $\Phi_D,\Phi_T$ be the weight-$4$ Eisenstein series
> attached to the Domb and $(12,4,16)$ systems. Then the $2$-adic Eichler integrals
> $D_q^{-3}\Phi$ have constant terms $\tfrac7{24}\zeta_2(3)$ and $\tfrac7{32}\zeta_2(3)$
> after Teichmüller normalisation.

and the mechanism is that $\zeta_2(3)$ is the value of the $2$-adic $L$-function of the
*same* Eisenstein class, so the two normalisations differ only by the ratio of the
Eisenstein constant terms — which is a computation with $q$-expansions, not a transcendence
input. **This is the piece of the Catalan paper that was only sketched, and it is a
strictly easier target for $\zeta(3)$**, because $\zeta_2(3)$ is Calegari's theorem
whereas the Catalan case needed $\zeta_2(2)$ via Rivoal's moving Padé.

### 7.4 What the framework buys beyond alignment

1. **It predicts $\sigma_p$ without computing anything**: slopes occur exactly at primes
   of bad reduction of the singular fibre, i.e. $\sigma_p=v_p(c)$, and $\sum_p\sigma_p\log p
   =\log|\lambda_1\lambda_2|$ is the product formula for the extension class. That is
   §5.1's design rule in disguise: *the resource a row can supply is a birational
   invariant of its extension class, $\log\Lambda$, and the two-lattice method just
   redistributes it between the archimedean and $p$-adic realisations.*
   The corollary is stark: **no amount of cleverness in the lattice can exceed
   $\gamma\log\Lambda_{\rm dec}$.** The $\zeta(3)$ failure is therefore not a failure of
   bookkeeping — it is a statement about the extension class of the $(12,4,16)$ motive.
2. **It says where to look for a better row**: maximise $\log\Lambda$, i.e. find an
   Eisenstein extension class whose singular fibre is as degenerate as possible, which
   for modular curves means high level with a cusp of large width. Concretely: scan
   $X_0(N)$ Hauptmoduls for large $|1/t_{\rm sing}|$. That is a finite, mechanical search
   over genus-zero levels and, as far as I can tell, **nobody in this project has run it.**
3. **It predicts new recurrences.** Every torsion section of the universal elliptic curve
   over a genus-zero $X_0(N)$, paired with an Eisenstein class of weight $w+2$, yields an
   Apéry-like row. The sporadic lists (Zagier's six, Almkvist–Zudilin's six, Cooper's
   three) are the ones that happen to have *integral* $a_n$; but the framework says the
   natural object is the full family and integrality is a separate, checkable condition
   (the class being integral at all finite places). Systematically enumerating
   $(N,\text{cusp},w)$ and testing integrality is a well-defined search for **new
   sporadic sequences** — the outcome that would be genuinely new mathematics
   independent of any irrationality application.

### 7.5 Honest status

Everything in §7 is a **conjectural framework with six confirming instances**. Nothing in
it is proved here. Its immediate cash value is §7.4(1) — the ceiling
$\gamma\log\Lambda_{\rm dec}$ — which is *not* conjectural: it follows from the product
formula and the master formula of §5, both of which are elementary. That ceiling is the
single most useful thing this investigation produced.

---

## 8. Empirical confirmation by LLL (`lattice/lll_harness/lll_zeta3.gp`)

Rather than trust §2's asymptotics, build the **actual** congruence lattice $K_n$ for each
$n$ (exact integers: $S_n=d_{3n}^3$, $T_n=2^{v_2(56\Delta_{2n,3n})}$, $M_n=S_nT_n$,
$K_n$ via `matsolvemod`), reduce it in the theory's anisotropic box
$|c_1|\le e^{xn}$, $|c_2|\le e^{(\sigma-x)n}$ with $x=6.70690$, $\sigma=17.31777$,
and measure $\log|q_n|/n$, $\log|q_n\zeta(3)-p_n|/n$ directly.

| $n$ | $\log q/n$ | $\log|q\zeta(3)-p|/n$ | $\delta_{\rm emp}$ | $v_2$ deficit | $\log M_n/n$ |
|---|---|---|---|---|---|
| 20 | 11.1007 | 0.5536 | 0.9501 | $-7$ | 16.71 |
| 40 | 11.3535 | 0.7465 | 0.9343 | $-7$ | 17.00 |
| 60 | 11.4090 | 0.8470 | 0.9258 | $+7$ | 17.39 |
| 80 | 11.5655 | 0.9955 | 0.9139 | $-7$ | 17.21 |
| 95 | 11.6796 | 1.1208 | **0.9040** | $-5$ | 17.46 |
| 100 | 11.5707 | 1.0132 | 0.9124 | $-3$ | 17.27 |
| **theory** | **11.7392** | **1.1627** | **0.90095** | $O(\log n)$ | **17.3178** |

All three rates converge to the predicted values, from above, with the expected
$O(\log n/n)$ corrections. **Every constant in §2 is confirmed independently.** In
particular $\log|q_n\zeta(3)-p_n|/n\to+1.16>0$: the linear forms provably grow, and
the failure is real, not an artefact of lossy estimates.

*Methodological warning worth keeping.* An earlier version of this harness LLL-reduced
with a huge weight on the error coordinate and reported $\delta_{\rm emp}\approx1.5$–$1.9$.
That is **meaningless**: a rank-2 lattice of covolume $M$ in the $(q,\ q\theta-p)$ plane
contains Dirichlet-quality vectors ($\delta\ge2$) for trivial Minkowski reasons, so
unconstrained LLL measures the *lattice's* density, not the *construction's* worthiness.
**Any empirical $\delta$ must be measured inside the prescribed anisotropic box.**
This trap will bite anyone who runs the "measure it with LLL in an afternoon" programme
of §6.3 on rank $\ge3$; state it prominently there.

---

## 9. RESOLVED: the modular parametrisation of the T $=(12,4,16)$ row

The open sub-task of §3 is closed. **T is a level-8 system.**

$$\boxed{\ t(\tau)=\left(\frac{\eta(\tau)\,\eta(8\tau)}{\eta(2\tau)\,\eta(4\tau)}\right)^{\!8}
=q-8q^2+\cdots,\qquad F(\tau)=\sum_{n\ge0}a^T_n\,t(\tau)^n\ }$$

with

$$F=-\tfrac16E_2(\tau)+\tfrac16E_2(2\tau)-\tfrac13E_2(4\tau)+\tfrac43E_2(8\tau),$$
$$F=1+4q+8q^2+16q^3+24q^4+24q^5+32q^6+32q^7+24q^8+52q^9+\cdots$$

*Verified* to $O(q^{26})$. Two independent confirmations that this is right:
1. **Modularity is exact, not fitted.** $E_2$ is quasi-modular with anomaly $\propto1/d$ for
   $E_2(d\tau)$; the combination satisfies $\sum_d c_d/d=-\tfrac16+\tfrac1{12}-\tfrac1{12}+\tfrac16=0$
   **exactly**, so $F$ is a genuine weight-$2$ modular form on $\Gamma_0(8)$ — as the
   $F\cdot D_q^{-3}\Phi$ theory demands. Constant term $-\tfrac16+\tfrac16-\tfrac13+\tfrac43=1$ ✓.
2. $a^T_n$ for odd $n$ satisfies $[q^n]F=4\sigma_1(n)$ ($n=1,3,5,7,9,25$: $4,16,24,32,52,124$) —
   pure Eisenstein behaviour.

Note the exact parallel with Apéry's row, which the same pipeline recovers as its
acceptance test:
$$\text{Apéry }(17,5,1):\quad t=\Bigl(\tfrac{\eta(\tau)\eta(6\tau)}{\eta(2\tau)\eta(3\tau)}\Bigr)^{12},
\quad F=1+5q+13q^2+23q^3+29q^4+30q^5+31q^6+\cdots\ \ (\Gamma_0(6)),$$
$$\text{T }(12,4,16):\quad t=\Bigl(\tfrac{\eta(\tau)\eta(8\tau)}{\eta(2\tau)\eta(4\tau)}\Bigr)^{8},
\quad F=1+4q+8q^2+16q^3+\cdots\ \ (\Gamma_0(8)).$$

**Why the §3 search failed — a methodological correction.** It assumed $F$ was itself an
*eta quotient* on the same divisor set and rejected everything else. $F$ is an *Eisenstein
combination*, not an eta quotient, so the test produced a **false negative** across all of
levels $6,8,12,16,24,32,48,64$ (8 exhaustive runs, 0 hits — all worthless). The correct
filter is far simpler and far more robust: accept $t$ iff the coefficients of
$F=\sum a_nt^n$ **grow slowly** (bounded by a small polynomial), which is what a modular
form of low weight must do and what a generic $t$ never does. Anyone running the §6.4
programme should use the growth filter, not a shape hypothesis.

**Consequence for §3's open problem.** Both rows now have explicit weight-$2$ $F$ and hence
explicit weight-$4$ $\Phi=F\cdot D_qt$: Domb on $\Gamma_0(6)$, T on $\Gamma_0(8)$. The route to
$\varepsilon=3\xi^D_2-4\xi^T_2=0$ is now fully concrete — compare the constant terms of the
two $2$-adic Eichler integrals $D_q^{-3}\Phi$ against $\zeta_2(3)$. Both levels are
$2$-power-divisible ($6=2\cdot3$, $8=2^3$), which is exactly the condition for the
$2$-adic slope, consistent with §7.4(1). (Domb's own eta-quotient Hauptmodul was **not**
recovered by this scan on $\{1,2,3,6\}$ or $\{1,2,3,4,6,12\}$ with $|r|\le24$; its level-6
parametrisation in the literature must use a different normalisation. Still open.)

---

## 10. Correction to §5.1 — the design rule for order $>2$ rows

§5.1 uses $\gamma\log\frac1\lambda+G=\gamma\log\Lambda$, which relies on
$\Lambda\lambda=c=\prod_i\lambda_i$. That is an **order-2 identity**. For a scalar
recurrence of order $m$ with roots $\lambda_1,\dots,\lambda_m$ the product formula reads
$$\sum_p\sigma_p\log p=\log\Bigl|\prod_{i=1}^m\lambda_i\Bigr|,$$
so the correct row score, valid at any order, is
$$\textbf{score}=\log\tfrac1{|\lambda_2|}-k,\qquad
\textbf{harvestable budget}=\textbf{score}+\log\Bigl|\prod_i\lambda_i\Bigr|,$$
which reduces to $\log|\lambda_1|-k$ only when $m=2$ (our Domb, T, Apéry, Zagier rows —
so §§1–9 are unaffected). Using $\log\lambda_1-k$ on a higher-order row **overstates**
the budget, sometimes wildly. Two instances found in `ROW_LEDGER.md` (see its CORRECTION
section): the level-12 $\zeta(5)$ row has harvestable budget $-5$, not $+5.54$; the
level-24 $\zeta(7)$ row has $\sigma_2=5$, not $0$ as measured. Both remain hopeless.

---

## 11. External cross-check: `certificates/padic/padic-apery-limits.tex`

River supplied the canonical source set on 2026-08-21. One document bears directly on
§§1, 3 and 7. Two independent confirmations and one sharp boundary.

### 11.1 Confirmations

Its Table `tab:families` (the canonical fifteen sporadic pairs, after Malik–Straub /
Straub 2023) lists

| # | label | $A(n)$ | Apéry limit | $w$ | $\chi$ |
|---|---|---|---|---|---|
| 7 | $(\alpha)$ | $\sum_k\binom nk^2\binom{2k}k\binom{2(n-k)}{n-k}$ (Domb) | $7\zeta(3)/24$ | 3 | $1$ |
| 8 | $(\gamma)$ | $\sum_k\binom nk^2\binom{n+k}n^2$ (Apéry) | $\zeta(3)/6$ | 3 | $1$ |
| 10 | $(\varepsilon)$ | $\sum_k\binom nk^2\binom{2k}n^2$ | $7\zeta(3)/32$ | 3 | $1$ |

So **T $=(12,4,16)$ is the family $(\varepsilon)$**, independently confirming the binomial
form $a^T_n=\sum_k\binom nk^2\binom{2k}n^2$ found in §0, the limit $\tfrac7{32}\zeta(3)$,
and the denominator exponent $k=w=3$ of §1.1 — all three derived here from scratch.

More usefully: $(\alpha),(\gamma),(\varepsilon)$ are the **only** $\zeta(3)$ families among
the canonical fifteen. Since $(\gamma)$ (Apéry) has $c=1$ and hence no slope anywhere,
**Domb + T is the only pair that can exist**, not merely the only one our census found.
The negative result of §4 is therefore sharper than stated: it is not "we looked and found
nothing better", it is "the canonical classification contains nothing better".

### 11.2 The boundary on §7 — the strong form of realisation rigidity is false

The paper's `Finding \ref{find:cross}`:

> \[verified: $p=7$, all pairs of the fifteen families, 12–21 digits\] Pairwise
> integer-relation searches for $(1,\Lambda^{\mathrm{fam}_1}_1,\Lambda^{\mathrm{fam}_2}_1)$
> find nothing below the noise floor.

This must be confronted head-on, because families $C$, $F$, $s_{18}$ all carry the period
$L_{-3}(2)$ (with archimedean ratios $\tfrac12,\tfrac58,\tfrac12$) and their limits are
*still* unrelated. **The strong form of §7.1 — "a relation in $\mathrm{Ext}^1$ transfers to
every realisation" — is therefore false as stated.**

It is not a refutation of the conjecture as actually written in §7.1, because the two
papers study **different $p$-adic objects**:

| | this file (§1.3) | `padic-apery-limits.tex` |
|---|---|---|
| object | $\xi_p=\lim_{n\to\infty}b_n/a_n$ in $\mathbf Q_p$ | $\Lambda_a=\lim_{s}\chi(p)^sp^{ws}B(ap^s)/A(ap^s)$ |
| index set | all $n$ | the tower $n=ap^s$ |
| primes | $p\mid c$ (slope $\sigma_p=v_p(c)>0$) — **ramified** | $p\ge5$, $p\nmid c$ — **unramified** |
| exists because | Casoratian $\Rightarrow$ Cauchy | descent congruence $(\mathrm{LB}_w^\chi)$ |

For Domb ($c=64$) and T ($c=16$) the only slope prime is $p=2$; at $p=7$ both have
$\sigma_7=0$ and $b_n/a_n$ has no naive limit at all. So the two findings are about
disjoint sets of primes and cannot contradict. §7.1 was already stated only for primes of
positive slope, and survives verbatim.

But the honest consequence is a **restriction, not a re-labelling**: the extension-class
relation is visible in the ramified/slope realisation and demonstrably invisible in the
unramified tower realisation. That is itself informative for §7 — it says $\Lambda_a$ is
*not* the syntomic realisation of the same $\mathrm{Ext}^1$ class in a way that respects
the $\mathbf Q$-structure. The paper's own reading ("an object about Frobenius rather than
about archimedean asymptotics") is consistent with $\Lambda_a$ being a Frobenius-eigenvalue
invariant living in the $\Gamma_p$-value ring, not a period of the extension.
**§7.4's claim that the framework predicts relations "at every prime" must be struck;
it predicts them at slope primes only.**

### 11.3 It does not close the $\varepsilon=0$ gap — but it shows its shape

The paper works at unramified $p\ge5$; §3's gap is at the ramified prime $p=2$. It offers
no theorem that applies. What it does offer is a strong structural hint: its engine is the
twisted descent law
$$(\mathrm{LB}_w^\chi)\qquad v_p\bigl(p^wB_nA_q-\chi(p)B_qA_n\bigr)\ \ge\ w,\qquad q=\lfloor n/p\rfloor,$$
verified for all 15 pairs, 28 primes $5\le p\le103$, with the floor *exactly* $w$ in every
cell. That is the same shape as our $\Delta_{m,n}=3a^T_nb^D_m-4a^D_mb^T_n$: a two-index
cross-determinant with a valuation floor. The difference is that theirs relates one family
to itself across a digit-shift, ours relates two families at a common prime.
**Proposed route to $\varepsilon=0$: find the $p=2$, cross-family analogue of
$(\mathrm{LB}_w^\chi)$.** A descent congruence $v_2(\lambda a^T_nb^D_q-\mu a^D_qb^T_n)\ge$
const with $q=\lfloor n/2\rfloor$ would give $\varepsilon=0$ by iteration, without ever
identifying either limit — which is exactly the unconditional route §3 argued was
unavailable by Casoratian methods alone. This is the most promising concrete attack on the
one open gap, and it did not exist before today.

---

## 12. A reduction of the $\varepsilon=0$ gap to a finite congruential statement

§3 argued that the determinant bound *cannot* be obtained without identifying the two
$2$-adic limits, and §11.3 suggested looking for a $p=2$ descent congruence. This section
carries that out. The gap is now reduced to **one inequality about integer sequences mod
powers of $2$**, with no transcendental input and no modular forms.

### 12.1 The single-family law does extend to $p=2$

`padic-apery-limits.tex` proves/verifies $(\mathrm{LB}_w^\chi)$ only for $p\ge5$. At the
ramified prime it still holds, with a smaller floor. *Verified* $2\le n\le200$, $q=\lfloor n/2\rfloor$:

$$v_2\bigl(2^3b_na_q-b_qa_n\bigr)\ \ge\ 2\qquad\text{for both Domb and T,}$$

with the minimum $2$ attained (so the floor is $2$, not $w=3$). This is a new data point for
that programme — the ramified prime was outside its stated range — but the floor does not
grow with $n$, so iterating it gives only the tower limit, not $\varepsilon=0$.

### 12.2 The doubling law

The useful structure is in the **cross-family** determinant. Write
$v(m,n)=v_2(\Delta_{m,n})$, $\Delta_{m,n}=3a^T_nb^D_m-4a^D_mb^T_n$. Along doubling chains
$(2^km_0,2^kn_0)$ the valuation *doubles*, exactly:

| base $(m_0,n_0)$ | $v$ along the chain $(m_0,n_0),(2m_0,2n_0),\dots$ |
|---|---|
| $(4,6)$ | $18,\ 42,\ 90,\ 186,\ 378,\ 762$ — exactly $v_{k+1}=2v_k+6$ |
| $(1,1)$ | $2,\ 5,\ 13,\ 29,\ 61,\ 125,\ 253,\ 509,\ 1021$ — exactly $v_{k+1}=2v_k+3$ from $k\ge1$ |
| $(3,2)$ | $11,\ 17,\ 33,\ 65,\ 129,\ 257,\ 513$ — exactly $v_{k+1}=2v_k-1$ |
| $(1,2)$ | $4,\ 9,\ 22,\ 47,\ 95,\ 191,\ 383,\ 767$ — exactly $v_{k+1}=2v_k+1$ |
| $(5,7)$ | $22,\ 49,\ 105,\ 217,\ 441,\ 889$ — exactly $v_{k+1}=2v_k+7$ |

**Control.** The same test on the *non-aligned* pair Apéry $(17,5,1)$ + Domb
(no common slope prime) gives $-4,-7,-10,-13,-16,-19$ along the chain from $(4,6)$:
linear decrease, no doubling. The test discriminates.

Grid data for the inequality $v(2m,2n)\ge2v(m,n)-C$:
min defect $-20$ over 16900 pairs $m,n\le130$; $-11$ over the 559 pairs in the balanced
band $|4n-6m|\le12$ (the only regime the argument below uses); band-by-band minima
$-9,-10,-8,-13,-12,-8$ with no downward drift.

### 12.3 The reduction

> **Lemma (reduction).** Suppose there exist $(m_0,n_0)$ and a constant $C$ such that
> $$v_2(\Delta_{2m,2n})\ \ge\ 2\,v_2(\Delta_{m,n})-C \tag{$\ast$}$$
> holds for every $(m,n)=(2^km_0,2^kn_0)$, $k\ge0$, and $v_2(\Delta_{m_0,n_0})>C$.
> Then $\varepsilon:=3\xi^D_2-4\xi^T_2=0$.

*Proof.* Put $v_k=v_2(\Delta_{2^km_0,2^kn_0})$ and $u_k=v_k-C$. By $(\ast)$,
$u_{k+1}\ge2u_k$, and $u_0>0$, so $v_k\ge C+2^k(v_0-C)\to\infty$.

Now suppose $\varepsilon\ne0$. From
$\Delta_{m,n}=a^T_na^D_m\bigl[3(\xi^D_m-\xi^D_2)-4(\xi^T_n-\xi^T_2)+\varepsilon\bigr]$
and the fact that the two tails $\to0$ $2$-adically, for all large $k$
$$v_k=v_2\bigl(a^T_{2^kn_0}\bigr)+v_2\bigl(a^D_{2^km_0}\bigr)+v_2(\varepsilon).$$
Here is the point: $s_2(2^kn)=s_2(n)$, so by §1.2
$v_2(a^T_{2^kn_0})=3s_2(n_0)$ is **constant** along the chain, and
$v_2(a^D_{2^km_0})\le3s_2(m_0)+9$ is bounded. Hence $v_k$ is *bounded* — contradicting
$v_k\to\infty$. Therefore $\varepsilon=0$. $\square$

With $(m_0,n_0)=(4,6)$: $v_0=18$, and $C=12$ suffices on all tested data. Note the
hypothesis needs $(\ast)$ **for every $k$**, which no finite computation supplies — that,
and only that, is what remains.

### 12.4 Why this is progress, and what is still missing

**Progress.** The open problem was "prove that two $2$-adic limits coincide", whose only
known route ran through $2$-adic Eichler integrals of weight-4 Eisenstein series and
Calegari-style overconvergence — analytic, heavy, and still not done for $p=2$. It is now
"prove one inequality between $v_2$ of two explicit integer expressions". No limits, no
modular forms, no transcendence. $(\ast)$ is exactly the shape of a Dwork/Frobenius
congruence, and congruences of that shape are known for Apéry-like sequences
(Beukers; Malik–Straub; and $(\mathrm{LB}_w^\chi)$ itself). The natural mechanism is a
matrix congruence $M_{2n}\equiv U_nM_n$ on the companion matrices
$M_n=\bigl(\begin{smallmatrix}a_n&b_n\\a_{n-1}&b_{n-1}\end{smallmatrix}\bigr)$, whose
determinant would square the valuation.

**Still missing — and I want to be exact about this.**
1. $(\ast)$ is **verified, not proved**, and verification can never settle it: if
   $\varepsilon\ne0$ with $v_2(\varepsilon)$ huge, every chain plateaus at $v_2(\varepsilon)$
   and $(\ast)$ fails only beyond computational reach. **The numerics in 12.2 are not
   evidence for $\varepsilon=0$** — they are consistent with any sufficiently large
   $v_2(\varepsilon)$. Their whole value is that they make $(\ast)$ a plausible target.
2. There is a real risk that proving $(\ast)$ for all $k$ is *equivalent* to the original
   problem rather than easier than it. I do not have an argument that it is strictly easier;
   I have an argument that it is better-posed.
3. Even if $(\ast)$ were proved, $\delta=0.9010<1$ and **no irrationality result follows**.
   Closing this gap would make §2's theorem unconditional; it would not make it useful.
   The honest value is methodological: it would validate the two-lattice machine on a case
   with a known answer, which was the original point of the exercise.

---

## 13. Multi-prime bridging, and what the complete classification says

### 13.1 The bridge generalises to several primes — for free

The $2$-adic limit is never an object of study in §2; it is a *device*. It forces a
divisibility on an integer determinant, that divisor enters a rational lattice, and
$\mathbf Q_2$ vanishes from the conclusion. The correlated-lattice step needs only
$T_n\mid h_n$ — nothing about $T_n$ being a prime power. So if $T_n=\prod_p p^{g_pn}$ the
Smith-normal-form index bound $[\mathbf Z^2:K_n]\le S_nT_n$ goes through verbatim, and the
master formula holds with

$$\boxed{\ G=\sum_p\min\bigl(\sigma_p^{\rm eng}\alpha,\ \sigma_p^{\rm dec}\gamma\bigr)\log p\ }$$

This is a real generalisation: neither the Catalan paper nor §2 has it. But note the
ceiling is unmoved. Since $\sum_p\sigma_p^{\rm dec}\log p=\log|c^{\rm dec}|$ and
$\log\frac1\lambda+\log|c|=\log\Lambda$, the resource is still at most
$\gamma\log\Lambda_{\rm dec}$. **Multi-prime bridging does not raise the ceiling; it only
lets you reach it when the engine's slopes dominate at *every* prime.** The balance ratio
$\alpha:\gamma$ is a single number, so at most one prime can be balanced; the others
contribute their unbalanced minimum.

### 13.2 It has no instance in the canonical fifteen

Computing $c$ and factoring, across all twelve sporadic families:

| | $c$ | slope primes |
|---|---|---|
| A | $-8$ | $2^3$ |
| B | $27$ | $3^3$ |
| C | $9$ | $3^2$ |
| D | $-1$ | none |
| E | $32$ | $2^5$ |
| **F** | $\mathbf{72}$ | $\mathbf{2^3\cdot3^2}$ |
| $(\delta)$ | $81$ | $3^4$ | 
| $(\zeta)$ | $-27$ | $3^3$ |
| $(\alpha)$ Domb | $64$ | $2^6$ |
| $(\eta)$ | $125$ | $5^3$ |
| $(\varepsilon)$ T | $16$ | $2^4$ |
| $(\gamma)$ Apéry | $1$ | none |

**F is the unique row with two slope primes.** Its period $(5/8)L(2,\chi_{-3})$ is shared
with $B$, $C$ and $s_{18}$ — all of which have $\sigma_2=0$. Since the per-prime term is a
*minimum* over the two rows, the $p=2$ contribution vanishes for every pair containing F.
So multi-prime bridging is available in principle and **realisable nowhere in the canonical
classification.** That is why neither paper needed it, and it is a clean negative rather
than an oversight.

### 13.3 The synthesis: budget is not the binding constraint

Ranking everything now known by budget $=\log\Lambda-k$:

| row | budget | has a same-period decayer? |
|---|---|---|
| **cusp-form row** ($L(f,2)$, $f=\eta_2^3\eta_6^3$) | $\mathbf{+0.773}$ | **no — none known** |
| Apéry $\zeta(3)$ | $+0.525$ | is its own (it *is* the 1979 proof) |
| D, Apéry $\zeta(2)$ | $+0.406$ | is its own |
| C, F ($L(2,\chi_{-3})$) | $+0.197$ | **no — no row in the class has $|\lambda_2|<1$** |
| T ($7\zeta(3)/32$) | $+0.149$ | yes — Domb is the engine $\Rightarrow$ $\delta=0.901$ |
| A, E | $+0.079$ | E: yes, Zudilin's row $\Rightarrow$ $\delta=0.9025$ |
| everything else | $<0$ | — |

The pattern across the complete classification is now unambiguous, and it is not about
budget:

> **The two-lattice method has produced a non-trivial construction exactly twice, and both
> times the decaying partner came from outside the modular world** — Zudilin's
> hypergeometric row for Catalan, and T (which decays only because
> $\lambda_2=12-8\sqrt2<1$ is an accident of its quadratic) for $\zeta(3)$. Every class
> where the method dies — $L(2,\chi_{-3})$ with four aligned rows, the new cusp-form row
> with the largest budget of all — dies for the *same* reason: **no decayer**, not
> insufficient budget.

This sharpens `THEORY_NOTES_03` §3 item 2 ("cross-world pairs are the template") from a
heuristic into a conclusion drawn over the complete list. It also redirects the search:
for the cusp-form row and for $L(2,\chi_{-3})$, **stop looking for modular partners and
look for hypergeometric ones** — Beukers/Rivoal Padé rows, Rhin–Viola integrals — which is
precisely where Catalan's decayer came from. A hypergeometric row with period $L(f,2)$ for
the weight-3 CM newform of level 12 would be the single most valuable object to find; the
engine on that side is already the best one known ($\sigma_2=6$, $k=2$, $\log\rho_2=\log4$),
and it needs only $\Lambda_{\rm dec}>14.8$.

---

## 14. CORRECTION: the product-formula ceiling is false for non-integral rows

River pointed out that the $G_2$ bridge got Catalan to $1-\varepsilon$ worthiness, which
cannot be reconciled with the ceiling asserted in §5.1, §7.4(1) and §13.1. He is right and
those claims are wrong. This section states the corrected law.

### 14.1 The measurement

Zudilin's Catalan row, computed directly (`/tmp/zud.gp`, $M=160$):

* $v_2(Q_m)=-4m+2s_2(m)$ — **exact, zero exceptions for $1\le m\le160$**. The row is *not
  integral*: it carries $2$-power denominators at rate $\kappa_2=4$ per index.
* measured $2$-adic convergence slope $\to8$ ($d/m=7.72,\,7.80,\,7.73$ at $m\approx100$–$120$,
  approaching $8$ with the expected $O(\log m/m)$ correction from the $2s_2(m)$ term).
* characteristic roots $\varphi^{5},\varphi^{-5}$, so $c=\lambda_1\lambda_2=1$ and $v_2(c)=0$.

A slope of $8$ at a row whose $\log|c|=0$. The product formula
$\sum_p\sigma_p\log p=\log|c|$ — the basis of the whole "budget" invariant — **fails**.

### 14.2 The corrected law

Let $\kappa_p\ge0$ be the rate of $p$-power denominators, $v_p(a_n)=-\kappa_pn+O(\log n)$.
The Casoratian slope law $(1)$ of `THEORY_NOTES_03` reads
$v_p(\text{increment})=(n-1)v_p(c)-\ldots-v_p(a_n)-v_p(a_{n-1})$, and the two $-v_p(a)$
terms each contribute $+\kappa_pn$. Hence

$$\boxed{\ \sigma_p=v_p(c)+2\kappa_p\ }$$

Zudilin: $0+2\cdot4=8$ — matches the measurement exactly. Domb: $6+0=6$ ✓. T: $4+0=4$ ✓.
Catalan's E row: $5+0=5$ ✓.

Clearing the denominators costs $\eta=\kappa_p\log p$ per index, so the **net** resource is
$(\sigma_p-\kappa_p)\log p=(v_p(c)+\kappa_p)\log p$, and the corrected budget is

$$\boxed{\ \textbf{budget}=\log\Lambda-k+\sum_p\kappa_p\log p\ }$$

The $\kappa$ term is a **second, independent source of $p$-adic resource, invisible to the
product formula**. An integral row ($\kappa\equiv0$) is capped at $\log\Lambda-k$ as §5.1
said; a non-integral row is not capped there at all.

### 14.3 The corrected ranking — and why Catalan wins

Per index, Zudilin's row: $\log\Lambda=5\log\varphi=2.40605$, $\kappa_2\log2=4\log2=2.77259$,
and $k=4$ (its $P_m$ needs $D_{2m-1}^2$, i.e. $2\cdot2$ per index $m$):

$$\textbf{budget}_{\rm Zudilin}=2.40605+2.77259-4=+1.17864 .$$

| row | $\kappa_2$ | budget | 
|---|---|---|
| **Zudilin (Catalan decayer)** | **4** | $\mathbf{+1.1786}$ |
| cusp-form row $L(f,2)$ | 0 | $+0.7726$ |
| Apéry $\zeta(3)$ | 0 | $+0.5255$ |
| D, Apéry $\zeta(2)$ | 0 | $+0.4061$ |
| C, F | 0 | $+0.1972$ |
| T | 0 | $+0.1490$ |
| A, E | 0 | $+0.0794$ |

**Zudilin's row has by far the largest budget of anything known, and every unit of its lead
comes from $\kappa_2$.** That is the mechanism behind $1-\varepsilon$, and it is exactly what
§13.3's table omitted, because that table tacitly assumed integrality.

### 14.4 What this changes

1. **§7.4(1) is retracted.** I billed "no amount of cleverness can exceed
   $\gamma\log\Lambda_{\rm dec}$ — the failure is a statement about the extension class" as
   the most useful output of the investigation. It is true only for integral rows. The
   $\zeta(3)$ failure is *not* a fact about the $(12,4,16)$ motive; it is a fact about the
   two rows we happened to have, both of which are integral.
2. **The search criterion inverts.** Do not look for rows with $p\mid c$. Look for rows with
   large $p$-power **denominators**. This also explains §13.3's "cross-world" observation at
   a deeper level: hypergeometric rows (Zudilin, Nesterenko, Rhin–Viola) acquire $\kappa_p>0$
   from their factorial/binomial normalisations, while modular rows have integral $a_n$ and
   $\kappa\equiv0$. "You need a cross-world partner" is really "**you need a partner with
   $\kappa_p>0$, and only the hypergeometric world supplies it**".
3. **A live target for $\zeta(3)$ that did not exist before.** Both Domb and T are integral —
   worse, $v_2(a^T_n)=3s_2(n)-[n\text{ odd}]\ge0$ is *positive*, which slightly *reduces* the
   slope. The $\zeta(3)$ construction has no $\kappa$ anywhere. So: **find a hypergeometric
   $\zeta(3)$ row with $2$-power denominators** — Rhin–Viola's $\zeta(3)$ integrals are the
   obvious place — and pair it with Domb as engine. That is a genuinely different attack on
   §2, and by §14.2 its budget could exceed anything in the table above.
4. §12's reduction is untouched: it concerns $\varepsilon=0$, not the quality number.

---

## 15. River's Catalan/Nesterenko CSVs — the bridge quantified, and one open discrepancy

`catalan-2-row-denominators/` contains computed Zudilin, Nesterenko and combined rows to
$n=98$, plus `catalan_arithmetic_diagnostics.csv` tracking the cross-determinant
$H_n=X_nU_n-Y_nV_n$. This is the best empirical data in the project and it settles several
things.

### 15.1 A fourth independent confirmation of the master formula

Measured worthiness $\delta=\text{digits}\cdot\log_210/\log_2q$:

| row | $\delta$ at $n=98$ |
|---|---|
| Zudilin alone | $0.5265$ |
| Nesterenko alone | $0.5128$ |
| **combined two-row lattice** | $\mathbf{0.6460}$ |

The Nesterenko-alone figure matches `NearCritical.lean`'s
$1-E_N/A_N=0.50968$. And the master formula (§5) applied to the Zudilin $\times$ Nesterenko
pair with $\sigma=12$ — i.e. **pure denominators, no $2$-adic bridge** — gives

$$F=\tfrac{E_Z+E_N-12}{2}=7.746367,\quad H=F+(A_N-E_N)=22.707996,\quad
\delta=\boxed{0.6588705}$$

against the LCM paper's independently derived $\mathbf{0.6588}$. Fourth instance, after
$0.857914$, $0.9025266$, $0.9009532$. (River's measured $0.6460$ sits $2\%$ low and is
drifting down, so the agreement is "same construction", not a precision match — finite-$n$
effects, or `best_q_n` is a best-vector rather than the prescribed box.)

### 15.2 What the bridge is worth, exactly

Same pair, same formula, with the $2$-adic divisor at the Lean threshold
$\sigma=12+k^*\log2$, $k^*=22.3512905953$:

$$F=0\ \text{exactly},\qquad H=14.961629,\qquad \delta=1 .$$

$$\boxed{\ \text{the }G_2\text{ bridge is worth } 0.6589\ \longrightarrow\ 1.0000\ \text{on this pair}\ }$$

This is the sharpest possible answer to River's point. The bridge is not a marginal
improvement on 0.9; on the Nesterenko pair it is the entire distance from 0.659 to
worthiness 1. And $\delta=1$ is a knife-edge: $F=0$ exactly, so **any** positive excess over
$k^*$ gives $F<0$ and irrationality.

### 15.3 The measured divisor — and why the obvious reading is wrong

`catalan_arithmetic_diagnostics.csv` gives $v_2(H_n)$. Fitting $v_2(H_n)=an+b\log_2n+c$
over $n\ge40$: $a=24.518$ (pure linear fit $24.114$; marginal increments $23.875$).

$$v_2(H_n)/n\ \approx\ 24.5\ \log 2\qquad\text{versus}\qquad k^*=22.351 .$$

**This is not a proof that Catalan is irrational, and the reason matters.** Roughly
$18\log2$ per $n$ of that valuation is a *common* $2$-power factor of the two rows — the
integrality clearing, $2^{e_m}$ with $e_m=4m$ on Zudilin (§14) and $2^{14n}$ on Nesterenko
(both visible in `NesterenkoUnconditional.lean` as `rate_pow_two_14`). A common factor of
$X$ and $U$ divides $H_n$ but **cancels**: dividing both rows by it changes nothing. Only
the excess is usable:

$$24.5-18\ \approx\ 6.5\log2\ \text{ per }n\ \text{net},$$

far short of $k^*=22.35$. This is exactly the $\kappa$/$\eta$ bookkeeping of §14 — a row with
$p$-power denominators gets slope $\sigma_p=v_p(c)+2\kappa_p$ but pays $\kappa_p$ to clear,
and the raw determinant valuation double-counts the part that is paid for.

### 15.4 RESOLVED — the bridge measured end-to-end

The apparent discrepancy is possibility (2), settled by direct computation
(`lattice/catalan_bridge/measure_bridge.py`).

**Cause.** `catalan_rows.py` already samples Zudilin at $3n$ (`zudilin_row_at_3n`), so
sampling was never the issue. Its `combine_rows` builds `congruence_lattice_basis(Y,U,S)` —
**one** congruence, modulus $S$ alone — and divides outputs by $S$. There is no $T_n$
anywhere. It implements the *LCM-square* lattice (as its own banner says), not the
correlated $2$-adic lattice. Hence $\sigma=\log S=12$ nats exactly and $\delta\to0.6589$.
**The CSVs measure the no-bridge construction; there is no conflict with the formalisation,
and my "factor of 2.3" alarm was an error of comparison.**

**Adding the missing congruence.** Implementing the architecture properly — $\alpha_r=X_r/S$,
mixed minor $h_n=\alpha_1U_n-\alpha_2Y_n$, $T_n=2^{v_2(h_n)}$, $M_n=S_nT_n$, and the *two*
congruences $\alpha\cdot c\equiv0\ (T_n)$, $(Y,U)\cdot c\equiv0\ (M_n)$ — yields integral
$q,p$ at every $n$ tested:

| $n$ | $v_2(h_n)/n$ | $\log_2M_n/n$ | $\log_2\mathrm{covol}/n$ | $\log q/n$ | $\log\lvert qG-p\rvert$ | $\delta$ | baseline |
|---|---|---|---|---|---|---|---|
| 20 | 25.050 | 41.985 | 41.985 | 14.268 | $-13.1$ | $1.04596$ | $0.65624$ |
| 40 | 24.600 | 41.818 | 41.530 | 14.558 | $-15.4$ | $1.02640$ | $0.65202$ |
| 60 | 24.483 | 41.758 | 41.601 | 14.693 | $-15.4$ | $1.01742$ | $0.64730$ |
| 80 | 24.337 | 41.760 | 41.652 | 14.781 | $-13.7$ | $1.01160$ | $0.64409$ |
| 98 | 24.316 | 41.426 | 41.385 | 14.770 | $-18.0$ | $1.01244$ | $0.64635$ |

**The bridge takes $\delta$ from $0.646$ to $\approx1.01$, converging to $1$ from above** —
the first end-to-end empirical measurement of the $G_2$ trick in this project.

**The mechanism is in the right-hand column.** $\log|q_nG-p_n|$ does *not* decrease: it sits
between $-13$ and $-23$ for all $n$. The linear form is $O(1)$, about $3\times10^{-7}$, and
does **not** tend to $0$. So $\delta\to1$ purely because $\log q$ grows while the numerator
stays bounded. This is exactly worthiness $1$ — the knife-edge $F=0$ — precisely what
`catalan_worthiness_one_sub_eps` asserts, and it is **not** irrationality, which needs
$\log|qG-p|\to-\infty$. Empirics and formalisation agree completely.

The $\approx1.2$ nats by which $\sigma$ appears to exceed $E_Z+E_N=27.49$ is a finite-$n$
artefact: $\log S_n/n=11.86$ at $n=98$ against its limit $12$, and $v_2(h_n)/n$ is still
drifting down ($25.05,24.60,24.48,24.34,24.32$), plausibly toward $k^*=22.351$.

**Caveat.** $\delta$ here is the best of four reduced-basis vectors, not the prescribed
anisotropic box (§8). The bounded-error signature is characteristic and hard to fake, but a
box-constrained rerun is the rigorous version.

---

## 16. UNRESOLVED CONTRADICTION — do not build on this until checked

Following §15.4 I ran the box-constrained selection (§8's methodology) on the
Zudilin $\times$ Nesterenko pair from River's CSVs. The result appears to produce integer
pairs with $|q_nG-p_n|\to0$ exponentially, which would prove Catalan's constant irrational.
**Catalan's irrationality is open, so I am near-certain this is an error of mine. I have not
found it.** Recording the full state so it can be checked.

### 16.1 What was measured

Inputs verified against `NearCritical.lean` (`lattice/catalan_bridge/`):
$\log|X_n|/n\to27.4$ ($A_Z=27.539$), $\log|X_nG-Y_n|/n\to12.98$ ($E_Z=13.100$),
$\log|V_n|/n\to29.16$ ($A_N=29.355$), $\log|V_nG-U_n|/n\to14.21$ ($E_N=14.393$).
Orientation correct; $q,p$ integral at every $n$ (checked $nq\bmod M=0$).

$v_2(h_n)=24n+O(\log n)$ — and $24=8\times3$ is exactly Zudilin's $2$-adic slope $\sigma_2=8$
(§14) times the sampling $3n$. So $T_n=2^{24n}$, $\sigma=\log S+24\log2\to12+16.636=28.636$.

The true index, from the Catalan paper's own gcd formula
$[\mathbf Z^2:K]=ST^2/\gcd(h,TY_1,TY_2,STa_1,STa_2,ST^2)$, confirms
$\log[\mathbf Z^2:K]/n=28.686$ at $n=98$ — i.e. index $=ST$, as the paper claims.

Hence $F=\tfrac12(E_Z+E_N-\sigma)=\tfrac12(27.493-28.636)=-0.571<0$.

### 16.2 The box-constrained measurement agrees with that prediction

| $n$ | $\log\lvert q\rvert/n$ | $\log\lvert qG-p\rvert/n$ | $\delta$ | predicted $F/n$ |
|---|---|---|---|---|
| 20 | 13.751 | $-1.160$ | 1.0844 | $-1.1192$ |
| 40 | 13.910 | $-0.959$ | 1.0689 | $-0.8473$ |
| 60 | 14.123 | $-0.838$ | 1.0593 | $-0.7763$ |
| 80 | 14.245 | $-0.730$ | 1.0513 | $-0.6721$ |
| 98 | 14.162 | $-0.796$ | 1.0562 | $-0.7676$ |

Measured $\log|qG-p|/n$ tracks the predicted $F/n$ closely, and both are **negative**.

### 16.3 Why this is almost certainly wrong

$\delta>1$ with $q_n\to\infty$ and $q_nG-p_n\ne0$ implies $G\notin\mathbf Q$. That is open.
Also `catalan_worthiness_one_sub_eps` is stated as $1-\varepsilon$ precisely because $\delta=1$
is the knife-edge; if $F<0$ were available the formalisation would say so.

Checks already done that did **not** find the error: rows match the Lean rates; orientation;
integrality of outputs; index by the paper's own formula; lattice is balanced
($\lambda_1\approx\lambda_2\approx\sqrt{\mathrm{covol}}$); the $q=0$ kernel vector
$\propto(V_n,-X_n)$ has size $e^{29n}$, far outside the box $e^{14.9n}$, so it cannot be the
Minkowski point; precision is ample (2500 dps against $q\sim2^{2003}$, error $\sim2^{-113}$).

Partial control: the *same* methodology in an independent PARI implementation (§8) on the
$\zeta(3)$ Domb+T pair returns $F_{\rm emp}=+1.12$ against theory $+1.16$ — correct sign and
value. So the machinery is not obviously broken.

### 16.4 Ranked suspects

1. **The Nesterenko CSV rows.** They are generated by `nesterenko_row()` via *numerical*
   hypergeometric evaluation (`J_n_hyper`, `dps_guard=50`), not by an exact recurrence.
   If $U_n$ is derived from a rounded quantity the pair may not be a genuine arithmetic row,
   and a tiny $|qG-p|$ would be an artefact. Against this: $|V_nG-U_n|\sim e^{14.2n}$ is
   enormous, not a rounding residue, and the CSV is internally consistent. **This is still
   the most likely culprit and the first thing to check** — regenerate the row from an exact
   recurrence or certified integer construction.
2. **$T_n=2^{v_2(h_n)}$ taken at its actual value rather than a uniform lower bound.** Using
   the exact per-$n$ valuation may smuggle in $n$-dependent structure the asymptotic
   bookkeeping does not license.
3. **My identification of $\sigma$ with both the covolume and the division modulus** — the
   step that lets one copy of $T_n$ pay for two congruences. Verified numerically here, but it
   is the conceptual heart of the construction and deserves re-derivation.

### 16.5 Status

**Do not treat this as a result.** It is a discrepancy between a computation and a known-open
problem, which in my experience means the computation is wrong. It is recorded because the
error is not where I looked, and because if suspect 1 is the cause it invalidates §15.4's
measured $\delta\approx1.01$ as well.

### 16.6 Suspect 1 eliminated — and the likely real resolution

**The Nesterenko row is sound.** Although `nesterenko_row()` defines
$U_n=\mathrm{round}(V_nG-\lambda_n)$ using the numerical value of $G$, the rounding lands on
the true arithmetic integer. Decisive test (`lattice/catalan_bridge/perturbation_test.py`):

| $n$ | $v_2(h_n)$ | with $U_n+1$ | with $U_n-1$ | with $V_n+S$ | $24n$ |
|---|---|---|---|---|---|
| 20 | 501 | 17 | 17 | 28 | 480 |
| 60 | 1469 | 19 | 19 | 34 | 1440 |
| 98 | 2383 | 20 | 20 | 37 | 2352 |

A $\pm1$ perturbation destroys the valuation completely ($2383\to20$). A rounding artefact
could not produce a $2383$-bit alignment, and a wrong integer could not survive it. So
$U_n$ is exact and $v_2(h_n)=24n+O(\log n)$ is a genuine arithmetic fact.

**Revised reading — this is probably not a contradiction at all.** The Lean development
defines $k^*:=(E_Z+E_N-12)/\log2$, i.e. *the break-even point*, and proves
worthiness $\ge1-\varepsilon$ from a $2$-adic package supplying rate $\ge k^*$. The
*observed* rate is $24$, and $24>k^*=22.351$. The natural explanation is that

> $k^*$ is what is **provable**; $24$ is what is **true**.

On that reading there is no inconsistency: the formalisation proves the strongest statement
its proved $2$-adic lower bound supports (worthiness $1$, hence $1-\varepsilon$), while the
actual valuation is larger, and the excess is precisely where Catalan's irrationality sits.

**This makes the target identical in shape to §12.** The observed rate has a structural
explanation — $24=8\times3$ is Zudilin's $2$-adic slope $\sigma_2=8$ (§14) times the
sampling $3n$ — so the target is:

> **Prove $v_2(h_n)\ge(22.352)\,n$ for all large $n$** (anything strictly above
> $k^*=22.3512905953$ suffices; the observed rate $24n$ leaves margin $1.65n$).
> This would prove Catalan's constant irrational.

As in §12, no limits, no modular forms — a valuation lower bound on an explicit integer
sequence. Unlike §12, here the *conclusion* is a famous open problem rather than a known one,
so the prior against it succeeding is correspondingly heavy.

**Residual caution.** I could not find an error, but "I could not find an error" is weak
evidence when the conclusion is this strong. The specific thing to check with the Lean
project: **what $2$-adic lower bound do the `nestFormInputs` fields actually prove?** If they
prove a rate $>k^*$, then either irrationality follows and has been missed, or my
$\sigma$/covolume identification (§16.4 suspect 3) is wrong. That single question decides it.

### 16.7 Control: the same pipeline on $\zeta(3)$, where the answer is known

`lattice/catalan_bridge/zeta3_control.py` runs the *identical* Python box-constrained
pipeline on the $\zeta(3)$ Domb+T pair, whose correct answer is known from §2 and §8
($F=+1.163$, $\delta=0.9010<1$):

| $n$ | $\log\lvert q\rvert/n$ | $\log\lvert q\zeta(3)-p\rvert/n$ | $\delta$ | predicted $F/n$ | $v_2(h)/n$ |
|---|---|---|---|---|---|
| 20 | 9.170 | $+2.510$ | 0.7263 | $+0.632$ | 12.400 |
| 40 | 9.423 | $+2.731$ | 0.7102 | $+0.840$ | 12.275 |
| 50 | 9.398 | $+2.767$ | 0.7055 | $+0.828$ | 12.360 |

**Correct sign, correct side of 1**, and $v_2(h)/n\to12$ exactly as §1.4 predicts for the
$2{:}3$ sampling. The predicted $F/n$ climbs toward $+1.163$ ($0.448,0.632,0.758,0.840$).

So the same code returns $F>0,\ \delta<1$ on $\zeta(3)$ and $F<0,\ \delta>1$ on Catalan.
**The pipeline discriminates correctly; the Catalan sign is not a coding artefact.** Combined
with §16.6 (the Nesterenko row is exact), suspects 1 and the generic-bug hypothesis are both
eliminated. What remains is either the "provable versus true $2$-adic rate" reading of §16.6,
or suspect 3 (the $\sigma$ = covolume = division-modulus identification).

### 16.8 Lean audit — §16.6 was wrong; the question sharpens

An audit of `lean-catalan-worthiness` (sonnet subagent, sources read not built) returns:

* **`dvd_reduced_cross_N`** (`RequestProject/NesterenkoCross.lean:249–283`):
  $$2^{24n}\ \big|\ \Delta_n,$$
  exponent exactly $24n$, linear, no correction term. **Proved**, via
  `nest_four_B_v2`, `nestJform_ne_zero`, `nestJform_v2` down to explicit hypergeometric
  valuation computations; no `sorry`, no `axiom`, no threaded hypothesis in the chain.
  So the measurement $v_2(h_n)=24n$ of §16.1 reproduces a *formalised theorem*.
* **`kstarNC` is not a $2$-adic quantity.** It is the zero of the *archimedean*
  $F_{\rm NC}(k)$ (`NearCritical.lean:211`). §16.6's reading — "$k^*$ provable, $24$ true" —
  is therefore **wrong and is retracted**.
* `TexpNC_le` (`NearCriticalAssembly.lean:80–86`) licenses $2^{\lfloor kn\rfloor}\mid\Delta_n$
  for **any** $k\le24$. Yet `exists_k_deltaNC_gt_one_sub` (`NearCritical.lean:303–330`)
  always selects $k=k^*-t<k^*$.
* Sorry/axiom inventory: one `sorry`, `Challenge.lean:20`, the comparator placeholder
  discharged by `Solution.lean`. No axioms, no `admit`. The `#guard_msgs` assertions record
  `catalan_worthiness_one_sub_eps_nesterenko`, `catalan_worthiness_857914`,
  `catalan_worthiness58_9025` as depending only on `propext, Classical.choice, Quot.sound`
  (as checked in, not re-executed).

Independently, the box membership check (`lattice/catalan_bridge/`) confirms the vector my
selection finds lies **inside** the prescribed box at every $n$ tested
($\log_2|c_1|=2114.1$ vs bound $2116.8$ at $n=98$; comparable margins at $40,60,80$).

**So the state is:** $k=24$ is proved and licensed; $F_{\rm NC}(k)>0$ iff $k<k^*$; the
development declines the licensed $k$ and takes $k<k^*$. Everything now turns on one
question, which decides between my suspect 3 and a real discrepancy:

> **Is the cap $k<k^*$ forced by a hypothesis of the selection/Minkowski theorem that fails
> at $k=24$, or is it merely the strongest conclusion the authors chose to state?**

Follow-up audit in flight. Until it returns, §16 remains "do not build on this".
