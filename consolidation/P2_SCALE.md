# P2$'$ at scale: the parity, the ratio $\rho$, and the continued fraction to $n=10\,000$

*Fable, 2026-08-24.  Scripts: `lattice/p2_scale/`, data in `lattice/p2_scale/data/`.
Extends `P2_STRUCTURE.md` ($n\le120$, 351 instances) and `P2_HOLONOMIC.md` ($n\le200$, 591
instances) to $4\le n\le10\,000$ at four values of $k$ — **39 988 instances**, every one of them
carrying an exactness certificate.  Tags: **[proved]**, **[verified]** = exact integer
computation over a stated range, **[measured]**, **[open]**.  All lattice arithmetic is exact in
$\mathbb Z^2$; Catalan's constant enters only through a dyadic bracket on one real number, and
every instance is computed twice, once at each end of the bracket, and accepted only if the two
runs return identical integer vectors.*

**No claim of irrationality is made anywhere below.**

---

## 0. Verdict

1. **The headline number.  $\log\rho_n$ has no detectable exponential growth: the slope of
   $\log\rho$ against $n$ is $-2.2\cdot10^{-6}\pm1.4\cdot10^{-6}$ pooled over the four $k$**, i.e.
   $$\bigl|\text{exponential growth rate of }\rho\bigr|\ <\ 5\cdot10^{-6}\quad(2\sigma),$$
   against $|F(23.9)|=0.5367$ and $|F(22.4)|=0.0169$: factors of $1.1\cdot10^{5}$ and
   $3.4\cdot10^{3}$ of headroom.  The $n\le120$ bound of `P2_STRUCTURE.md` §4.2 was $\pm0.0023$;
   this is $790$ times sharper.  All four per-$k$ OLS slopes, all four Theil–Sen slopes and all
   four winsorised slopes lie in $[-3.7\cdot10^{-6},+1.1\cdot10^{-6}]$.  **P2$'$ is not
   falsified; the Catalan ledger is unchanged.**  **[measured]**
2. **The statement P2$'$ actually needs is verified outright, for *every* $n$, not just along a
   subsequence.**  Since $\rho\ge1$ always, the content of P2$'$ is that $\tfrac1n\log\rho_n$
   returns to $0$.  Measured maxima of $\tfrac1n\log\rho_n$ **over all $n$ in a range**:
   $$\max_{8001\le n\le10\,000}\tfrac1n\log\rho_n=0.00090,\ 0.00067,\ 0.00099,\ 0.00071
     \qquad(k=22.35,\,22.4,\,23.0,\,23.9),$$
   a factor $17$–$25$ below $|F(22.4)|=0.0169$ at *every single $n$* in that range, and decaying
   like $C/n$ across the whole scan.  **[verified]**
3. **The run record does not survive, and its failure is exactly what a fair coin predicts.**
   The longest run of consecutive $n$ with $v_1$ outside the cone is now $12,10,13,11$ (it was
   $5,4,4$ at $n\le120$ and $6$ at $n\le200$).  A $10^5$-replicate fair-coin null of the same
   length and frequency gives mean longest runs $12.63,12.58,12.63,12.85$ with
   $\Pr[\text{null}\ge\text{observed}]=0.71,0.99,0.46,0.94$: the observed records are, if
   anything, *shorter* than chance, and the whole run-length spectrum matches $N2^{-(L+2)}$ term
   by term.  Longest runs grow like $\log_2N$, so "no bad run longer than $5$" was an artefact of
   $N=117$.  The right replacement is **in every window of twelve consecutive $n\le10\,000$
   there is an $n$ with $\rho\le\sqrt2$** — exact: the largest window-minimum of $\log\rho$ over
   windows of $12$ is $0.3465,0.3449,0.3381,0.3454$ against $\log\sqrt2=0.346574$, and over
   windows of $8$ it is $0.82,0.65,0.55,0.76$, so twelve is sharp.  (The longest run of $n$ with
   $\rho>\sqrt2$ is $11,10,11,11$ — shorter than the odd-parity run, because Corollary D lets an
   odd-parity $n$ still have $\rho\le\sqrt2$.)  **[verified + measured]**
4. **The parity is still a fair coin, at fifty times the range.**  Even-frequency
   $0.4998,0.5017,0.4996,0.4937$ ($9\,997$ values of $n$ each; pooled
   $0.49872\pm0.00250$); binomial $p=0.98,0.75,0.95,0.22$; runs tests pass; $96$
   autocorrelations at lags $1..24$ have $\max|z|=2.63$; $236$ tests of $n\bmod m$ for $m\le60$
   produce no Bonferroni survivor; the four parities are mutually independent; and there is no
   drift of the even-density with $n$ ($|{\rm slope}|\le1.9\cdot10^{-6}$ per unit).  The
   criterion $v_1\in\pm\mathcal P\iff i(n)$ even holds in **39 988 of 39 988** instances.
   **[measured + verified]**
5. **Gauss–Kuzmin holds to four decimals on $4.3\cdot10^{9}$ partial quotients**
   ($\chi^2=4.29$ on $8$ df, $p=0.83$), and the lengths obey Lévy to four figures — but
   **the partial quotient *at* the balance index is not Gauss–Kuzmin distributed, and this is
   new.**  The selected rung is size-biased: $a_{i(n)}$ has $\Pr[a=1]=0.4992$ against $0.4150$
   and $\Pr[a\ge9]=0.1118$ against $0.1520$, while $a_{i(n)+1}$ has $\Pr[a=1]=0.0758$ and
   $\Pr[a\ge9]=0.4188$, mean $149.9$, maximum $826\,227$ ($\chi^2=1\,379$ and $32\,360$ on $8$
   df).  **The balance falls just before a large partial quotient** — the mechanism behind every
   large $\rho$ — and $\log a_{i(n)+1}$ has no trend in $n$ (slopes
   $\le8.7\cdot10^{-6}\pm5.3\cdot10^{-6}$), with $\max_{n>5000}\tfrac1n\log a_{i(n)+1}\le0.0021$:
   the escape clause of P2$'_{\rm par}$ is quantified, not assumed.  Moreover **the parity and
   the size of $a_{i(n)+1}$ are independent** (point-biserial $|r|\le0.013$, $p\ge0.18$), which
   is what closes the one route by which P2$'$ could plausibly fail.  **[measured]**
6. **The midpoint law was mis-stated, and the corrected form fits.**  The measured
   $i(n)/L(n)=0.52146,0.52142,0.52108,0.52058$ ($n>1000$, s.e. $1\cdot10^{-5}$) is $200\sigma$
   from the $0.5238$ of `P2_STRUCTURE.md` §2.3.  The prediction must carry the Absorption
   content $m_Z$:
   $$\frac{i(n)}{L(n)}=\frac12+\frac{\log m_Z-\log h_{22}-(E_2-E_1)n}{2\bigl(\log h_{11}-\log m_Z\bigr)}
     \;=\;0.52141,\ 0.52138,\ 0.52103,\ 0.52053,$$
   which matches to $\le5\cdot10^{-5}$.  **[measured]**
7. **The $2$-adic law is exactly $24-k$, not $24.02-k$.**  Fitted slopes of
   $v_2(q_n\xi_2-p_n)$ against $n$: $1.649946,\,1.599949,\,0.999945,\,0.099943$ (s.e.
   $1.3\cdot10^{-5}$) against $24-k$; the residual $v_2-(24-k)n$ is genuinely $O(1)$ — mean
   $\approx+2.3$, sd $3.7$, and it never leaves $[-12,+21]$ over 39 988 instances.
   **[verified + measured]**
8. **$c=0$ is now decisive.**  $\log\gcd(X_n,Y_n,V_n,U_n)=49.20$ at $n=10\,000$ (maximum
   $72.82$ over the range), against the $0.196n+11.7$ fit that `POSITIVITY_PROGRAM.md` §3.4
   could not exclude at $n\le120$ and which predicts $1\,972$ — out by a factor $40$.  Likewise
   $h_{22}$ is polynomial ($e^{4.38}n^{1.39}$, never exceeding $e^{49.4}$) and
   $\tfrac1n\log m_Z=2.37225$ is identical to five decimals at all four $k$.
   **[verified $n\le10\,000$]**
9. **Every structural claim of the two predecessor documents survives verbatim at 39 988
   instances**: the closed forms for $h_{11},h_{22},h_{12}$ (Proposition 1 of `P2_HOLONOMIC.md`),
   $39\,988/39\,988$; $v_1$ a continued-fraction convergent at the balance index (Theorem 1 of
   `P2_STRUCTURE.md`), $39\,988/39\,988$ — the reduction loop takes zero steps every time;
   Proposition B's bound $|j|\le3$, attained values $\{2,3\}$ only; Proposition C and
   Corollary D, no violation.  **[verified]**
10. **The near-critical value $k=22.35$ behaves like the others.**  It is *below*
   $k_*=22.3512905953$, so asymptotically $F=+0.00045>0$ and it can never prove anything; yet at
   $n\le10\,000$ the finite-$n$ gcd deficit still makes the measured rate negative
   ($-0.00316$), and every parity, ratio and continued-fraction statistic is indistinguishable
   from those at $k=22.4,23.0,23.9$.  **The threshold is invisible to the lattice**; it lives
   entirely in the rate constant.  **[measured]**

---

## 1. What was computed

Rows, lattice and cone exactly as in `P2_STRUCTURE.md` §1 and `P2_HOLONOMIC.md` §1.1:
$$X_n=2^{e_{3n}}D_{6n}^2Q_{3n},\quad Y_n=2^{e_{3n}}D_{6n}^2P_{3n},\quad
  V_n=4^{7n+1}D_{6n}^2B_n,\quad U_n=4^{7n}D_{6n}^2C_n,$$
$S_n=D_{6n}^2$, $T_n=2^{\lfloor kn\rfloor}$, $M_n=S_nT_n$,
$\mathcal K_n=\{c\in\mathbb Z^2: M_n\mid c_ZX_n+c_NV_n,\ M_n\mid c_ZY_n+c_NU_n\}$,
oriented coordinates $(u,v)=(s_Zc_Z,c_N)$ with $s_Z=(-1)^n$, metric
$\mathrm{diag}(\lambda_Z,\lambda_N)$, $\lambda_Z=|X_nG-Y_n|/M_n$, $\lambda_N=(V_nG-U_n)/M_n$,
cone $=$ closed first quadrant, $\rho=(\text{cone minimum})/\lambda_1$.

**Range.**  $4\le n\le10\,000$, $k\in\{22.35,\,22.4,\,23.0,\,23.9\}$ — $39\,988$ instances.  The
first $k$ is the near-critical value ($k_*=22.3512905953$, so $k=k_*-0.0013$); the other three
are those of the two predecessor documents.  At $n=10\,000$ the integers involved have
$\log_{10}X_n=119\,677$ digits and the continued fraction of $h_{12}/h_{11}$ reaches
$L=221\,373$ partial quotients.

### 1.1 The exactness certificate

$G$ is the only real number in the problem, and it enters only through
$r=\lambda_N/\lambda_Z$.  Instead of computing in floating point we **bracket $r$ between two
dyadic rationals** and run every instance twice:

* Catalan's constant is truncated to $P=\lceil7.10\,n\rceil+200$ decimal digits, giving exact
  integers $A,W$ with $|G-A/W|\le2/W$.  ($7.10\,n$ dominates both horizons of
  `P2_HOLONOMIC.md` §2.3 — $6.267n$ for $s_Z$, $6.492n$ for the partial quotients — plus the
  $0.562n$ growth of $r$ itself.  At $n=10\,000$ that is $71\,200$ digits of $G$.)
* From $A,W$ one gets exact integer bounds $Z_{\rm lo}\le W|X_nG-Y_n|\le Z_{\rm hi}$ and
  $N_{\rm lo}\le W(V_nG-U_n)\le N_{\rm hi}$, hence integers $R_-\le r\,2^{256}\le R_+$.
  **In all $39\,988$ instances $R_+-R_-=1$**: $r$ is pinned to one unit in the last of $256$
  bits.
* The balance index, the reduced basis, the cone minimiser and the norm minimiser are computed
  **twice**, with the integer quadratic form $2^{512}x_1^2+R^2x_2^2$ and linear form
  $2^{256}x_1+Rx_2$ at $R=R_-$ and at $R=R_+$.  An instance is certified `ok` only if the two
  runs return **identical integer vectors**.  Since the true $r$ lies between the two, the
  returned vectors are then correct for the true metric.
* Additionally certified per instance: $s_Z=(-1)^n$; the two congruences that pin $h_{12}$;
  $[\mathbb Z^2:\mathcal K_n]=h_{11}h_{22}$ against an independent Hermite computation; and the
  integrality of the output pair $(q_n,p_n)$.

$$\boxed{\ \texttt{ok}=1\ \text{in}\ 39\,988\ \text{of}\ 39\,988\ \text{instances.}\ }$$

No floating-point comparison decides anything anywhere.  Reals are used only to print
logarithms.

### 1.2 What made $n=10\,000$ affordable

The $n\le200$ pipeline costs $13.5$ s per instance at $n=1000$ in `matkerint` alone.  Three
replacements:

| step | $n\le200$ pipeline | here |
|---|---|---|
| $h_{11},h_{12},h_{22}$ | `matkerint` + `mathnf` | the closed forms of `P2_HOLONOMIC.md` §1.2 and one modular inverse $h_{12}\equiv-s_Zh_{22}U_nY_n^{-1}$ |
| the convergent ladder | store all $L(n)\approx21n$ convergents | the three-term recurrence $e_i=a_ie_{i-1}+e_{i-2}$ for $e_i=h_{12}q_i-h_{11}p_i$, streamed twice: bit-length proxies first, exact squares only on the window where the proxy is within $12$ bits of its minimum |
| Gauss reduction | Lagrange from the Hermite basis, $O(L)$ steps in $3\,000$-digit reals | start from the ladder pair $(w_{i},w_{i-1})$, which Theorem 1 says is already reduced — the loop is still run and its verdict recorded |

Both shortcuts are *tested*, not assumed.

* The column `redst` counts reduction steps; it is $0$ in all $39\,988$ instances, i.e.
  **Theorem 1 of `P2_STRUCTURE.md` is re-verified at every instance.**
* The windowed balance search is controlled by brute force (`checkwin.gp`): re-evaluating the
  objective at *every* rung, with no window, at
  $n=7,37,113,250,400,555,700,850,999,1200,1500,1777,2000,2345,2700,3000,3500,4200,5000,5700,
  6000,7000,8500,10\,000$ and all four $k$ — ladders of up to $L=221\,373$ rungs — reproduces
  the windowed argmin in **96 / 96** cases.

Cost: $1.0$ s of wall time for all four $k$ at $n=1000$; $40$ min for $n\le3000$, $1$ h $45$ for
$3001\le n\le6000$, $3$ h $40$ for $6001\le n\le10\,000$, on 12 cores.  The last stage used an
*interleaved* split (worker $j$ takes $n\equiv j$ mod 12) so that the completed range is
contiguous at every moment rather than only at the end.

### 1.3 Validation against the published tables

The new pipeline reproduces the old ones exactly on their whole range
($4\le n\le120$ or $200$, $k\in\{22.4,23.0,23.9\}$, $351$ common instances):

| against | columns compared | mismatches |
|---|---|---|
| `p2_structure/data/struct_n120.csv` | `incone`, $\rho$ (all six printed digits), $\lambda_2/\lambda_1$ | **0 / 351** |
| `p2_holonomic/data/h2_balance200.csv` | $i(n)$, $L(n)$ | **0 / 351** |
| `p2_holonomic/data/h4_adelic.csv` | $v_2(q_n\xi_2-p_n)$ | **0 / 351** |

### 1.4 Validation of the extended rows

The rows past $n=200$ come from the fitted order-3, degree-14 P-recursion of
`p2_holonomic/rowrec.gp`, which was *fitted* on $n\le80$ and *checked* on $n\le120$.  The
extension is tested against three exact $2$-adic identities of `P2_HOLONOMIC.md` §5.1 that the
fit cannot have encoded (`verify.gp`), with $\xi_2=\zeta_2(2)$ built to $8\cdot(6n)$ binary
digits in a single non-storing pass:

$$v_2(X_n\xi_2-Y_n)=v_2(X_n)+24n-1-4s_2(3n),\qquad
  v_2\bigl(V_n\xi_2-U_n\bigr)=v_2(V_n)+28n-2s_2(n)-2s_2(3n)-1,$$
$$v_2(X_nU_n-V_nY_n)=v_2(X_n)+v_2(V_n)+24n-1-4s_2(3n).$$

**0 failures**, on 60 samples $n=50,100,\dots,3000$, 60 samples $n=100,200,\dots,6000$ and
20 samples $n=500,1000,\dots,10\,000$ — $420$ identity checks in all.  The middle identity is a
congruence to $28n+O(\log n)$ binary digits involving *only* the Nesterenko entries $B_n,C_n$ —
exactly the part of the cache the fitted operator produces — so it cannot hold by accident.
Also checked: $X_n,Y_n,V_n,U_n\in\mathbb Z$ at **every** $n\le10\,000$ (the recurrence runs over
$\mathbb Q$ throughout), and agreement with the cached $n\le200$ rows digit for digit.  Rates at $n=10\,000$:
$\tfrac1n\log X_n=27.5565$, $\tfrac1n\log V_n=29.3743$,
$\tfrac1n\log|X_nG-Y_n|=13.1202$ (against $E_1=13.0995887908$),
$\tfrac1n\log(V_nG-U_n)=14.4127$ (against $E_2=14.3931452672$).  The residual $\pm0.02$ is
the Chebyshev fluctuation of $\tfrac1n\log D_{6n}$ and does not decay monotonically; over the
$140$ samples it stays in $[-0.11,+0.10]$ for $E_1$ and $[-0.17,+0.08]$ for $E_2$, with no
trend.
**[verified]**

---

## 2. The headline: the growth rate of $\rho$

### 2.1 The regression

`POSITIVITY_PROGRAM.md` §3.3 regressed $\tfrac1n\log\rho$ on $n$; `P2_STRUCTURE.md` §4.2
corrected the specification to $\log\rho$ on $n$.  The latter is the one with content.

| $k$ | $|F(k)|$ | OLS slope, $4\le n\le10\,000$ | Theil–Sen (95% CI) | winsorised at $99\%$ |
|---|---|---|---|---|
| 22.35 | 0.00045 | $+0.55\cdot10^{-6}\pm2.91\cdot10^{-6}$ | $-0.04\cdot10^{-6}\ [-1.31,+1.18]\cdot10^{-6}$ | $+1.04\cdot10^{-6}$ |
| 22.4 | 0.01688 | $-3.65\cdot10^{-6}\pm2.81\cdot10^{-6}$ | $-0.03\cdot10^{-6}\ [-1.30,+1.18]\cdot10^{-6}$ | $-2.81\cdot10^{-6}$ |
| 23.0 | 0.22483 | $-2.76\cdot10^{-6}\pm2.85\cdot10^{-6}$ | $-1.02\cdot10^{-6}\ [-2.54,+0.15]\cdot10^{-6}$ | $-3.20\cdot10^{-6}$ |
| 23.9 | 0.53674 | $-2.67\cdot10^{-6}\pm2.84\cdot10^{-6}$ | $-0.23\cdot10^{-6}\ [-1.56,+0.93]\cdot10^{-6}$ | $-3.06\cdot10^{-6}$ |
| **pooled** | | $\mathbf{-2.17\cdot10^{-6}\pm1.43\cdot10^{-6}}$ | $2\sigma$ band $[-5.0,+0.7]\cdot10^{-6}$ | |

How the bound tightened as the range grew (at $k=23.9$; s.e. of the slope):

| range | $\le120$ | $\le200$ | $\le500$ | $\le1000$ | $\le2000$ | $\le3000$ | $\le6000$ | $\le10\,000$ |
|---|---|---|---|---|---|---|---|---|
| s.e. | $2.2{\cdot}10^{-3}$ | $1.0{\cdot}10^{-3}$ | $2.5{\cdot}10^{-4}$ | $8.6{\cdot}10^{-5}$ | $3.1{\cdot}10^{-5}$ | $1.7{\cdot}10^{-5}$ | $6.1{\cdot}10^{-6}$ | $2.8{\cdot}10^{-6}$ |

Split-half differences of the mean of $\log\rho$ (each over $4\,998$ units of $n$):
$-0.0002,\,-0.0284,\,-0.0172,\,-0.0101$, each $\pm0.0165$ — all four negative, at most
$1.8\sigma$.  The mis-specified model $\tfrac1n\log\rho$ against $n$ returns $-2\cdot10^{-7}$ at
every $k$, which is exactly the $c/n$ artefact that `P2_STRUCTURE.md` §4.2 identified.

**A cautionary reading of the intermediate ranges.**  At $n\le3000$ the pooled slope was
$+2.0\cdot10^{-5}\pm0.9\cdot10^{-5}$, a nominal $2.2\sigma$ *positive* value; doubling the range
to $6000$ moved it to $-2\cdot10^{-6}\pm3\cdot10^{-6}$ and tripling it to $10\,000$ to
$-2.2\cdot10^{-6}\pm1.4\cdot10^{-6}$.  That is what a fluctuation looks like, and it is a useful
calibration of how little a $2\sigma$ excursion means in this data.  **Any exponential growth of
$\rho$ is below $5\cdot10^{-6}$ per unit $n$ at $2\sigma$.**  **[measured]**

### 2.2 The statistic P2$'$ actually asks for

$\rho\ge1$ always (the cone minimum is a lattice vector), so $\tfrac1n\log\rho_n\ge0$ and P2$'$
— $\liminf\tfrac1n\log\rho_n=0$ — asks only that the sequence come back to $0$.  The data say
far more: it never leaves.

| range of $n$ | $\max_n\tfrac1n\log\rho_n$, $k=22.35$ | $22.4$ | $23.0$ | $23.9$ |
|---|---|---|---|---|
| $4$–$100$ | 0.24604 | 0.30905 | 0.37821 | 0.18253 |
| $101$–$300$ | 0.03837 | 0.04061 | 0.03516 | 0.03858 |
| $301$–$1000$ | 0.01637 | 0.01669 | 0.00907 | 0.01063 |
| $1001$–$2000$ | 0.00513 | 0.00952 | 0.00414 | 0.00465 |
| $2001$–$3000$ | 0.00377 | 0.00290 | 0.00294 | 0.00258 |
| $3001$–$4500$ | 0.00185 | 0.00200 | 0.00233 | 0.00240 |
| $4501$–$6000$ | 0.00114 | 0.00129 | 0.00164 | 0.00147 |
| $6001$–$8000$ | 0.00102 | 0.00112 | 0.00115 | 0.00102 |
| $8001$–$10\,000$ | **0.00090** | **0.00067** | **0.00099** | **0.00071** |

The maximum decays like $C/n$: $\log\rho_n=O(1)$ — the *stronger* form conjectured in
`POSITIVITY_PROGRAM.md` §3.3 — holds along **every** $n\le10\,000$, with
$\max_n\log\rho_n=10.363$.

Window form, exact:
$$\max_{\text{windows of }12\text{ consecutive }n\le10\,000}\ \min_{n\in\text{window}}\log\rho_n
\ =\ 0.3465,\ 0.3449,\ 0.3381,\ 0.3454\ \le\ \log\sqrt2=0.346574,$$
while windows of $8$ give $0.82,0.65,0.55,0.76$.  So **in every window of twelve consecutive
$n\le10\,000$ there is an $n$ at which the cone minimum is within $\sqrt2$ of the first
minimum**, and twelve is sharp.  This replaces the "window of six" of `P2_STRUCTURE.md` §4.2 and
the "window of seven" of `P2_HOLONOMIC.md` §6.  **[verified]**

### 2.3 Distribution of $\rho$, and the random-lattice comparison

| $k$ | median $\rho$ | $\Pr[\rho\le\sqrt2]$ | $\Pr[\rho>5]$ | $\Pr[\rho>10]$ | $\Pr[\rho>20]$ | $\Pr[\rho>50]$ | $\Pr[\rho>100]$ | median $\lambda_2/\lambda_1$ | $\max\rho$ |
|---|---|---|---|---|---|---|---|---|---|
| 22.35 | 1.4118 | 0.5411 | 0.1139 | 0.0570 | 0.0292 | 0.0110 | 0.0055 | 1.9242 | 7 051 ($n=2635$) |
| 22.4 | 1.4115 | 0.5417 | 0.1116 | 0.0553 | 0.0257 | 0.0095 | 0.0045 | 1.9152 | 31 670 ($n=1088$) |
| 23.0 | 1.4116 | 0.5423 | 0.1101 | 0.0556 | 0.0269 | 0.0102 | 0.0050 | 1.9184 | 7 198 ($n=8957$) |
| 23.9 | 1.4119 | 0.5360 | 0.1117 | 0.0550 | 0.0274 | 0.0109 | 0.0053 | 1.9357 | 4 655 ($n=5730$) |
| Haar $0.54/R$ | $\sqrt2$ | — | 0.1080 | 0.0540 | 0.0270 | 0.0108 | 0.0054 | $1.9336$ | — |

The Haar tail law of `P2_STRUCTURE.md` §4.3, fitted there on $117$ points per $k$, is now
confirmed on $9\,997$ per $k$ across more than two decades of $R$, and the median
$\lambda_2/\lambda_1$ agrees with the Haar value $1.9336$ to two decimals.  The median of $\rho$
is $\sqrt2$ to three decimals at every $k$, as Proposition C predicts.

**Proposition C and Corollary D.**  Four instances print $\rho=1.414214$ with
$v_1\in\pm\mathcal P$ ($k=22.35$ at $n=5938,6544$; $k=22.4$ at $n=2296$; $k=23.9$ at $n=6046$);
in each, $\rho_2=1$ exactly and $b_1$ lies on the diagonal, so $\rho=\sqrt2$ exactly and the
excess is the rounding of the printed sixth digit.  **No genuine violation.**  Corollary D:
$414,401,426,423$ instances have $\rho\le\sqrt2$ with $v_1$ outside the cone, and **all** of
them have $\lambda_2/\lambda_1\le2\sqrt2/\sqrt3=1.633$, as required.  **[verified]**

The record values of $\rho$ are instructive: each occurs at odd parity with
$\lambda_2/\lambda_1\approx\rho$ and a huge $a_{i(n)+1}$ — e.g. $k=22.4$, $n=1088$:
$a_{i+1}=222\,270$, $\lambda_2/\lambda_1=31\,347$, $\rho=31\,670$; $k=23.0$, $n=8957$:
$a_{i+1}=71\,183$, $\rho=7\,198$.  That is Proposition C's second branch and §4.2 below, in one
line.

---

## 3. The parity

### 3.1 No law, at fifty times the range of `P2_HOLONOMIC.md` §3.1

| $k$ | even $i(n)$ | binomial $p$ | runs $z$ ($p$) | longest odd run (starts at) | fair-coin null: mean, $\Pr[\ge\text{obs}]$ | longest run $\rho>\sqrt2$ |
|---|---|---|---|---|---|---|
| 22.35 | $4997/9997=0.4998$ | 0.984 | $-0.19$ (0.849) | 12 ($n=4035$) | 12.63, 0.708 | 11 |
| 22.4 | $5015/9997=0.5017$ | 0.749 | $-0.09$ (0.929) | 10 ($n=4125$) | 12.58, 0.992 | 10 |
| 23.0 | $4995/9997=0.4996$ | 0.952 | $-1.09$ (0.276) | 13 ($n=8827$) | 12.63, 0.458 | 11 |
| 23.9 | $4936/9997=0.4937$ | 0.215 | $+1.15$ (0.252) | 11 ($n=2836$) | 12.85, 0.936 | 11 |

Pooled even-density $19\,943/39\,988=0.49872\pm0.00250$.

* **The parity criterion $v_1\in\pm\mathcal P\iff i(n)$ even holds in $39\,988/39\,988$
  instances** — no exception over three and a half orders of magnitude of $n$.  **[verified]**
* **Autocorrelation**, lags $1..24$ at four $k$ ($96$ statistics): $\max|z|=2.63$.  For $96$
  independent standard normals the expected maximum is $2.7$.  Nothing.
* **$n\bmod m$, $m=2..60$** ($236$ contingency tests): smallest $p$-value $0.0033$
  ($k=22.35$, $m=53$), i.e. $p_{\rm Bonf}=0.19$; **no survivor at any $k$**, no replication
  across $k$.  The nominal hits at $m=3,9,11$ reported at $n\le200$ are gone.
* **Cross-$k$ agreement**: $0.4986,0.4990,0.4969,0.4974,0.5026,0.4929$ over the six pairs
  ($9\,997$ values each), $\max|z|=1.41$.  The four parities are independent.
* **No drift**: even-density by block of $1000$ stays in $[0.470,0.526]$ with a linear trend
  $\le1.9\cdot10^{-6}$ per unit $n$ at every $k$.
* The near-critical $k=22.35$ is statistically indistinguishable from the other three at every
  one of these tests.

### 3.2 The run record, and why it had to fail

`P2_STRUCTURE.md` §4.2 recorded "longest run of bad $n$ $=5,4,4$" at $n\le120$ and
`P2_HOLONOMIC.md` §6 upgraded it to $6$ at $n\le200$, framed as "in every window of seven
consecutive $n$ there is a good one".  **That claim does not survive**: the record is now $13$.

The fair-coin null makes the failure uninformative.  The longest head-run in $N$ fair tosses
concentrates at $\log_2N$: $6.9$ at $N=117$, $7.6$ at $N=197$, $13.3$ at $N=9997$.  The whole
run-length spectrum matches, term by term (observed against $N2^{-(L+2)}$):

| $L$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| $k=22.35$ | 1250 | 618 | 300 | 182 | 70 | 41 | 10 | 11 | 5 | 2 | 5 | 1 | — |
| $k=22.4$ | 1255 | 619 | 325 | 134 | 77 | 43 | 28 | 8 | 5 | 3 | — | — | — |
| $k=23.0$ | 1208 | 616 | 335 | 163 | 77 | 33 | 22 | 9 | 7 | 2 | — | — | 1 |
| $k=23.9$ | 1264 | 628 | 322 | 161 | 67 | 47 | 16 | 12 | 5 | 5 | 1 | — | — |
| fair coin | 1249.8 | 624.9 | 312.4 | 156.2 | 78.1 | 39.1 | 19.5 | 9.8 | 4.9 | 2.4 | 1.2 | 0.6 | 0.3 |

The honest form of the statement is therefore a rate, not a constant:

> **[measured]**  The longest run of consecutive $n\le N$ with $v_1(n)\notin\pm\mathcal P$ grows
> like $\log_2N$, as for a fair coin.  P2$'$ needs only that this stays $o(n)$, which the
> fair-coin model gives with enormous margin and which nothing observed contradicts.

A run of length $\ell$ costs P2$'$ nothing as long as $\ell=o(n)$: it merely puts gaps of that
size in the good subsequence.  And the run of $n$ with $\rho>\sqrt2$ — the quantity that
actually matters — is shorter still ($11,10,11,11$), because Corollary D lets an odd-parity $n$
have $\rho\le\sqrt2$ when the lattice is not skew.

---

## 4. The continued fraction at scale

### 4.1 Gauss–Kuzmin on $4.3\cdot10^{9}$ partial quotients

Pooled over all $39\,988$ instances (the continued fraction of $h_{12}/h_{11}$, length
$L(n)\approx21n$):

| $a$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | $\ge9$ |
|---|---|---|---|---|---|---|---|---|---|
| observed | .4150 | .1699 | .0931 | .0589 | .0406 | .0297 | .0227 | .0179 | .1520 |
| Gauss–Kuzmin | .4150 | .1699 | .0931 | .0589 | .0406 | .0297 | .0227 | .0179 | .1520 |

$4\,298\,332\,092$ partial quotients, $\chi^2=4.29$ on $8$ df, $p=0.830$.  (At $n\le120$:
$460\,365$ quotients, $\chi^2=5.09$.)  The lengths obey Lévy to four figures:

| $k$ | $L(n)/n$ measured ($n>6000$) | $0.8427\bigl(\tfrac1n\log h_{11}-\tfrac1n\log m_Z\bigr)$ |
|---|---|---|
| 22.35 | 21.164 | 21.162 |
| 22.4 | 21.194 | 21.192 |
| 23.0 | 21.544 | 21.542 |
| 23.9 | 22.070 | 22.068 |

**[measured]**

### 4.2 The partial quotient *at* the balance index is size-biased — new

This is the empirical content of P2$'$'s escape clause.  Over all $39\,988$ instances:

| $a$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | $\ge9$ | mean | median | max |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| $a_{i(n)}$ | .4992 | .1674 | .0812 | .0495 | .0327 | .0246 | .0184 | .0152 | .1118 | 10.50 | 2 | 17 547 |
| $a_{i(n)+1}$ | .0758 | .1367 | .1060 | .0787 | .0627 | .0504 | .0387 | .0322 | **.4188** | 149.87 | 6 | 826 227 |
| Gauss–Kuzmin | .4150 | .1699 | .0931 | .0589 | .0406 | .0297 | .0227 | .0179 | .1520 | $\infty$ | 2 | — |

$\chi^2=1\,379$ and $32\,360$ on $8$ df: **both are decisively non-Gauss–Kuzmin, in opposite
directions.**  The balance index prefers a rung preceded by a *small* quotient and followed by a
*large* one — exactly what the geometry dictates, since $A_i\approx h_{11}/q_{i+1}$ drops by the
factor $a_{i+1}$ at the next rung, and a large drop is what makes $i$ rather than $i+1$ the
argmin of $A_i^2+tB_i^2$.

Consequences for P2$'_{\rm par}$ ("$i(n)$ even, **or** $a_{i(n)+1}=e^{o(n)}$"):

* $\lambda_2/\lambda_1\in[0,1.52]\cdot a_{i(n)+1}$ at every $k$, correlation of logs
  $0.728$–$0.746$, confirming the $\lambda_2/\lambda_1\asymp a_{i+1}$ of `P2_STRUCTURE.md` §4.1
  on $9\,997$ points per $k$ instead of $117$.
* $\log a_{i(n)+1}$ has **no trend in $n$**: slopes $+3.3\cdot10^{-6},\,-2.3\cdot10^{-6},
  +8.7\cdot10^{-6},\,-6.2\cdot10^{-6}$, all with s.e. $5.3\cdot10^{-6}$.
* $\max_{n>5000}\tfrac1n\log a_{i(n)+1}=0.0016,\,0.0019,\,0.0021,\,0.0017$: even the single
  worst rung in the top half of the range is $e^{0.002n}$.  **The escape clause is satisfied
  with an order of magnitude to spare wherever it is needed.**
* Largest partial quotient *near* the balance: median $29$ within $i\pm5$, median $93$ within
  $i\pm20$, against a median of $194\,725$ over the whole continued fraction (global maximum
  $3\,303\,464\,313$).  The neighbourhood of the balance is *not* where the giant quotients
  live; the selected one is a moderate local record.
* **The parity and the size of $a_{i(n)+1}$ are independent.**  Point-biserial correlations
  $-0.0032,\,-0.0133,\,+0.0098,\,+0.0008$ ($p=0.75,0.18,0.33,0.93$; Mann–Whitney
  $p=0.89,0.25,0.38,0.77$).  This matters: the only way P2$'$ could fail is a run of odd
  parities coinciding with a run of giant $a_{i+1}$, and the two are uncorrelated.
  **[measured]**

### 4.3 The corrected midpoint law

`P2_STRUCTURE.md` §2.3 predicted $i/L=\tfrac12+\tfrac{E_2-E_1}{2\kappa}=0.5238$ with
$\kappa=\tfrac1n\log h_{11}$.  With $9\,997$ points per $k$ the measurement is $200\sigma$ from
that.  The derivation must use the *reduced* denominator $h_{11}/m_Z$ (the continued fraction is
that of $(h_{12}/m_Z)/(h_{11}/m_Z)$) and the balance condition
$\lambda_Z|e_i|=\lambda_Nq_ih_{22}$ with $|e_i|\approx m_Zh_{11}'/q_{i+1}$, giving
$q_iq_{i+1}\approx m_Zh_{11}'/(h_{22}r)$ and hence, with $\log q_i\approx(i/L)\log h_{11}'$,

$$\boxed{\ \frac{i(n)}{L(n)}=\frac12+\frac{\log m_Z-\log h_{22}-(E_2-E_1)n}
{2\bigl(\log h_{11}-\log m_Z\bigr)}\ }$$

| $k$ | measured ($n>1000$) | s.e. | old prediction | corrected prediction | residual |
|---|---|---|---|---|---|
| 22.35 | 0.52146 | 0.00001 | 0.52354 | **0.52141** | $+5\cdot10^{-5}$ |
| 22.4 | 0.52142 | 0.00001 | 0.52351 | **0.52138** | $+4\cdot10^{-5}$ |
| 23.0 | 0.52108 | 0.00001 | 0.52316 | **0.52103** | $+5\cdot10^{-5}$ |
| 23.9 | 0.52058 | 0.00001 | 0.52265 | **0.52053** | $+5\cdot10^{-5}$ |

The old form is out by $2.1\cdot10^{-3}$, the corrected one by $5\cdot10^{-5}$ — a uniform
residual at all four $k$, i.e. what is left is a genuine $O(1/\log h_{11})$ correction to the
Lévy step, not a mis-specification.  The sd of $i/L$ across $n$ is $0.00141$ at $n>1000$ (it was
$0.0137$ at $n\le200$, shrinking like $n^{-1/2}$ as it must).  **There is no drift**: the slope
of $i/L$ in $n$ is $|{\cdot}|\le2.2\cdot10^{-8}\pm0.9\cdot10^{-8}$.  **[measured]**

---

## 5. The $2$-adic side at scale

### 5.1 $v_2(q_n\xi_2-p_n)=(24-k)\,n+O(1)$, exactly

`P2_HOLONOMIC.md` §5.2 measured the slope as $24.02-k$ on $n\le200$.  At $n\le10\,000$:

| $k$ | fitted slope | s.e. | $24-k$ | difference | $v_2-(24-k)n$: mean, sd, min, max |
|---|---|---|---|---|---|
| 22.35 | 1.649946 | 0.000013 | 1.65 | $-5.4\cdot10^{-5}$ | $+2.41$, $3.70$, $-12$, $+21$ |
| 22.4 | 1.599949 | 0.000013 | 1.60 | $-5.1\cdot10^{-5}$ | $+2.34$, $3.73$, $-12$, $+18$ |
| 23.0 | 0.999945 | 0.000013 | 1.00 | $-5.5\cdot10^{-5}$ | $+1.95$, $3.75$, $-11$, $+19$ |
| 23.9 | 0.099943 | 0.000013 | 0.10 | $-5.7\cdot10^{-5}$ | $+2.38$, $3.74$, $-10$, $+19$ |

The "$24.02$" was a finite-$n$ artefact: **the surplus of the bridge over the modulus is exactly
$24-k$, and the remainder is bounded** ($|v_2-(24-k)n|\le21$ over $39\,988$ instances, sd
$3.7$).  This is consistent with the exact valuation formulas of `P2_HOLONOMIC.md` §5.1 ($24n$
from the Zudilin row, $\lfloor kn\rfloor$ absorbed by $T_n$), and it is the sharpest available
form of "the $2$-adic quality of the output is precisely what the construction did not spend".
**[verified + measured]**

### 5.2 The adelic quality

Measured on the cone minimisers at $n>4000$:

| $k$ | $\tfrac1n\log q_n$ | $\tfrac1n\log|q_nG-p_n|$ | $F(k)$ | $\delta_\infty$ | $\delta_2$ | $\delta_\infty+\delta_2$ |
|---|---|---|---|---|---|---|
| 22.35 | 14.9582 | $-0.00316$ | $+0.00045$ | 0.0002 | 0.0765 | **0.0767** |
| 22.4 | 14.9409 | $-0.02049$ | $-0.01688$ | 0.0014 | 0.0742 | **0.0756** |
| 23.0 | 14.7329 | $-0.22846$ | $-0.22483$ | 0.0155 | 0.0471 | **0.0626** |
| 23.9 | 14.4210 | $-0.54035$ | $-0.53674$ | 0.0375 | 0.0048 | **0.0423** |

Unchanged in substance from `P2_HOLONOMIC.md` §5.2 ($0.084,0.071,0.050$ at $n\le200$): the
construction is a factor $13$–$24$ below the trivial adelic Dirichlet bound
$\delta_\infty+\delta_2=1$ at every $k$, and the deficit worsens with $k$.  The $2$-adic place
is not a hidden resource.  **[measured]**

The measured archimedean rate sits $0.0036$ *below* $F(k)$ at all four $k$ (it was $0.0088$ at
$n\le3000$ and $0.0068$ at $n\le6000$).  That gap is $-c_n/2$ with
$c_n=\tfrac1n\log\gcd(X_n,Y_n,V_n,U_n)$, and its steady shrinkage as the range grows is the
statement $c=0$ of §6, seen from the other side.

---

## 6. Divisor data: $c=0$, decisively

`POSITIVITY_PROGRAM.md` §3.4 could not distinguish $\tfrac1n\log\gcd=0.196+11.7/n$ from
$-0.196+7.55\log n/n$ at $n\le120$; `P2_STRUCTURE.md` §2.1 argued for the second.  At
$n\le10\,000$ the question is settled by nearly two orders of magnitude:

| model | fit | residual sd | $R^2$ | prediction at $n=10\,000$ |
|---|---|---|---|---|
| logarithmic | $11.051+3.367\log n$ | 6.49 | 0.208 | $42.1$ |
| linear | $34.171+0.00091\,n$ | 6.81 | 0.128 | $43.3$ |
| the $n\le120$ linear fit | $11.7+0.196\,n$ | — | — | $1\,972$ |

**Observed: $\log\gcd(X_n,Y_n,V_n,U_n)=49.20$ at $n=10\,000$**, maximum $72.82$ over the whole
range.  The $n\le120$ linear reading is out by a factor $40$.  Hence
$$c=\lim_n\tfrac1n\log\gcd(X_n,Y_n,V_n,U_n)=0,\qquad \kappa_n=\sigma(k)-O(\log n/n),$$
confirming `P2_STRUCTURE.md` §0.6(b) and closing ranked next step #1 of
`POSITIVITY_PROGRAM.md`.  **[verified $n\le10\,000$]**

Two further divisor facts, both new at this range:

* $h_{22}=e^{4.38}\,n^{1.39}$ — **polynomial** in $n$, with $\log h_{22}\in[3.47,49.39]$ over
  $4\le n\le10\,000$, and $k$-independent.
* $\tfrac1n\log m_Z=2.37225$ at $n>6000$, **identical to five decimals at all four $k$**,
  confirming that the Absorption directional mass is an $S_n$-effect and carries no dependence
  on the $2$-adic modulus $T_n$.  (At $n\le120$ the estimate was $2.3911$, at $n\le3000$
  $2.37513$; the limit is approached from above like $\log n/n$.)

---

## 7. What changes in the record

1. **`P2_STRUCTURE.md` §4.2 / `P2_HOLONOMIC.md` §6, "longest run of bad $n$ is $5$ (resp. $6$)"
   — withdrawn.**  It is $13$ at $n\le10\,000$ and grows like $\log_2N$.  Replacement: *in every
   window of twelve consecutive $n\le10\,000$ there is an $n$ with $\rho\le\sqrt2$*, and the
   whole run-length spectrum matches a fair coin term by term.
2. **`P2_STRUCTURE.md` §2.3, $i(n)/L(n)=0.5238$ — corrected** to the $m_Z$-carrying form of
   §4.3, which fits to $5\cdot10^{-5}$.
3. **`P2_HOLONOMIC.md` §5.2, "$v_2(q_n\xi_2-p_n)/n=24.02-k$" — sharpened** to exactly $24-k$
   with a bounded remainder.
4. **`POSITIVITY_PROGRAM.md` §3.4, the gcd-loss constant — settled**: $c=0$, the linear
   alternative excluded by a factor $40$ at $n=10\,000$.
5. **`POSITIVITY_PROGRAM.md` §3.3 / `P2_STRUCTURE.md` §4.2, the slope of $\log\rho$** — the
   $2\sigma$ bound improves from $\pm0.0046$ ($n\le120$) to $\pm5.6\cdot10^{-6}$ per $k$
   ($\pm2.9\cdot10^{-6}$ pooled), and the "$\rho=e^{o(n)}$" reading is now supported at every
   single $n$, not merely on average.
6. **New: the balance index is size-biased** (§4.2) — $a_{i(n)}$ biased small and $a_{i(n)+1}$
   biased large, both decisively — and the parity is independent of that bias.
7. **Nothing in the Catalan ledger moves.**  P2$'$ was, and remains, the one unproved input of
   the supercritical construction; this computation removes the possibility that it fails for a
   visible reason below $n=10\,000$, and quantifies by how much.

---

## 8. Status of every claim

**Verified (exact, $4\le n\le10\,000$, $k\in\{22.35,22.4,23.0,23.9\}$, 39 988 instances).**
The exactness certificate of §1.1 (`ok` in $39\,988/39\,988$; bracket width $1$ in every
instance); Proposition 1 of `P2_HOLONOMIC.md` (both congruences and the index identity),
$39\,988/39\,988$; Theorem 1 of `P2_STRUCTURE.md` ($v_1$ is the ladder vector at the balance
index — the reduction loop takes zero steps in every instance); the parity criterion
$v_1\in\pm\mathcal P\iff i(n)$ even, $39\,988/39\,988$; Proposition B ($|j|\le3$; the attained
bound is $2$ or $3$), Proposition C (no violation) and Corollary D; integrality of every
$(q_n,p_n)$; the window statement of §2.2; the three $2$-adic row identities of §1.4 on 140
samples; the brute-force control of the balance window, 96/96; agreement with the three
published tables on all 351 common instances.

**Measured.**  The slope of $\log\rho$ and all its robustness variants; the parity statistics
and their fair-coin null; Gauss–Kuzmin and the two size-biased distributions at the balance
index; the corrected midpoint law; the $(24-k)n$ law; the divisor fits.

**Open, unchanged.**  P2$'$ itself; the density of even $i(n)$ (measured $0.49872\pm0.00250$
pooled, conjecturally $\tfrac12$); a proof that $\log\gcd(X_n,Y_n,V_n,U_n)=O(\log n)$; a proof
that $\lambda_2/\lambda_1\le Ca_{i(n)+1}$; the Nesterenko half of Proposition 3 of
`P2_HOLONOMIC.md`.

---

## 9. Scripts and data

`lattice/p2_scale/` — PARI/GP unless stated.  Run pattern (concatenate, do not `read()`):

    cat lattice/positivity/rows_pos.gp lattice/p2_structure/p2core.gp \
        lattice/p2_scale/scorel.gp lattice/p2_scale/srun.gp \
        PARAMS lattice/p2_scale/run_scan.gp > run.gp && gp -q run.gp

| file | prepend | what |
|---|---|---|
| `build_rows.gp` | `rows_pos.gp`, `p2core.gp`, `p2_holonomic/rowrec.gp` | exact rows to $n=N$; re-fits and re-checks the order-3 operator, re-verifies against the cached $n\le200$ |
| `scorel.gp` | `rows_pos.gp`, `p2core.gp` | `hermd` (closed-form Hermite data + three exactness checks), `gredx` (exact integer Gauss reduction), `conex`/`eightx` (exact cone scan), `balx` (two-pass streamed balance index) |
| `srun.gp` | + `scorel.gp` | `anal2`: one certified CSV line per $(k,n)$, 41 columns |
| `run_scan.gp` | + `srun.gp` | driver (`NL, NH, KLIST, OUTF, ROWF`, optional `NSTEP`) |
| `verify.gp` | `rows_pos.gp`, `p2core.gp` | the three exact $2$-adic row identities of §1.4; `zudsnap` builds $\xi_2$ in one non-storing pass |
| `checkwin.gp` | + `scorel.gp` | brute-force control for the windowed balance search |
| `split_rows.py` | — | per-worker row split, contiguous (equal cost) or interleaved (contiguous completed prefix) |
| `launch.sh`, `front.sh` | — | parallel launch; the contiguous front of an interleaved run |
| `stats.py`, `stats2.py`, `stats3.py` | — | §§2–6 |

| data file | content |
|---|---|
| `data/rows_scale.txt` | exact $n,X_n,Y_n,V_n,U_n$, $4\le n\le10\,000$ (2.5 GB; not kept in git, regenerable in $474$ s) |
| `data/scan_c*.csv`, `scan2_c*.csv`, `scan3_c*.csv` | the scan, 39 988 lines, 41 columns, header in each file |
| `data/verify.csv`, `verify6.csv`, `verify10.csv` | the $2$-adic row validation |
| `data/summary.json`, `summary2.json`, `summary3.json` | every number quoted above |

**Cost.**  Row cache to $n=10\,000$: $474$ s, $2.5$ GB, peak $4.9$ GB.  Scan: $40$ min
($n\le3000$) $+$ $1$ h $45$ ($3001$–$6000$) $+$ $3$ h $40$ ($6001$–$10\,000$), on 12 cores.
Row validation and window controls: about $1$ h more.  Precision: Catalan's constant to
$7.10\,n+200$ decimal digits, used only to produce the exact dyadic bracket; no floating-point
comparison decides anything.

---

## 10. Honest assessment

The reduction of `P2_STRUCTURE.md` — P2$'$ is the parity of a continued-fraction index — and
the negative result of `P2_HOLONOMIC.md` — that parity has no arithmetic memory surviving one
step in $n$ — are both confirmed at fifty times the range and sixty-eight times the number of
instances.  Nothing here brings a proof closer.  What has changed is the quality of the
evidence: the exponential-growth rate of $\rho$ is bounded by $5\cdot10^{-6}$ where it was
bounded by $0.0046$; the escape clause of P2$'_{\rm par}$ is quantified rather than assumed; and
one published constant ($0.5238$) and one published record ("no bad run longer than $5$") are
corrected — the first to a formula that fits, the second to a growth rate that a fair coin
explains.

The single genuinely new phenomenon is §4.2: **the balance index is not a uniformly chosen rung
of the ladder.**  It sits systematically just before a large partial quotient —
$\Pr[a_{i+1}=1]=0.076$ against Gauss–Kuzmin's $0.415$, and $\Pr[a_{i+1}\ge9]=0.419$ against
$0.152$.  That is what makes $\lambda_2/\lambda_1$ occasionally enormous, and $\rho$ with it.  It
is also the exact mechanism by which P2$'$ could conceivably fail: not through the parity, which
is a fair coin, but through a run of odd parities coinciding with a run of giant $a_{i+1}$.
Over $4\le n\le10\,000$ those two are independent to within $p=0.18$ at worst, the density of
even indices is $\tfrac12$ to $\pm0.0025$, and the size-bias of $a_{i+1}$ is $n$-independent.  A
proof of P2$'$ would need a lower bound on the density of even $i(n)$ — as `P2_STRUCTURE.md` §8
already concluded — and nothing at this scale suggests where such a bound would come from.
**[open]**
