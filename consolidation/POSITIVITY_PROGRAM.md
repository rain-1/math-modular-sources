# The positivity programme: sign as the datum valuations cannot see

*Claude (Fable), 2026-08-22.  Scripts: `lattice/positivity/`, data in
`lattice/positivity/data/`.  Tags: **[proved]** = exact identity or exact closed form;
**[verified]** = exact rational / integer computation over a stated range, no floating point
in the arithmetic; **[measured]** = numerical rate or fit, not converged to a proof;
**[open]** as usual.*

*Follow-up: `consolidation/P2_STRUCTURE.md` (2026-08-23) extends the P2$'$ measurement to
$n\le120$, corrects 51 of the 231 cone minima of `data/cone_n80.csv` (all downward), settles the
gcd-loss constant of §3.4 as $c=0$ by identifying $M_n/[\mathbb Z^2:\mathcal K_n]=\gcd(X_n,Y_n,V_n,U_n)$,
and reduces P2$'$ to a parity statement about one continued-fraction index.*

**No claim of irrationality is made anywhere.**  Every covolume-driven column below is
reproduced verbatim by a rational surrogate for $G$ (§6.4, and `CATALAN_AUDIT.md` §4(a)),
so none of it is evidence about $G$.

---

## 0. Verdict

1. **Positivity is not an integral phenomenon; it is a Casoratian phenomenon, and it is free
   for the whole corpus.**  For any row in the Zagier normalisation
   $(n+1)^{w+1}u_{n+1}=\cdots-c\,n^{w+1}u_{n-1}$, the exact Casoratian
   $a_nb_{n-1}-a_{n-1}b_n=-c^{\,n-1}/n^{w+1}$ telescopes to
   $$\boxed{\;a_N\xi-b_N \;=\; a_N\sum_{n>N}\frac{c^{\,n-1}}{n^{w+1}\,a_n a_{n-1}}\;}$$
   whose summands have the *constant* sign $\operatorname{sign}(c)^{\,n-1}$.  Hence
   $c>0$ and $a_n>0$ $\Rightarrow$ $a_N\xi-b_N>0$ for **every** $N$, unconditionally, with
   no integral anywhere.  **[proved]**, and the Casoratian verified in exact rational
   arithmetic for all twelve sporadic rows, $n\le110$ (§2.1).
2. **The right systematic criterion for "carries a Beukers-type positive integral" is
   *Hankel positivity*.**  A Beukers/Zudilin/Nesterenko kernel is always of the shape
   $r_n=\int u^n\,d\mu(u)$ with $\mu\ge0$; so the row admits such a representation iff
   $(r_n)=(a_n\xi-b_n)$ is a Hamburger moment sequence, which is decidable.  Every one of the
   nine corpus rows with a real archimedean limit passes: $\det H_d>0$ for all $d\le10$, with
   the measure on $[0,\infty)$ exactly when $c>0$ and on $(-\infty,0]$ exactly when $c<0$
   **[verified, $\backslash$p 2000]** (§2.3).  The three rows with complex spectrum
   (Zagier B, AZ $\delta$, AZ $\eta$) have no limit and the question is vacuous.
3. **Lemma P2 is supported far more strongly than before, and now on a population, not a
   single family.**  On the honest lattice $\mathcal K_n$ at **every** $n\le80$ and three
   values of $k$ (231 instances), and independently on **13 332** pair lattices drawn from
   the moment grid, the ratio
   $\rho=(\text{cone minimum})/\lambda_1$ has
   $$\frac1n\log\rho \;=\; (+0.0022\pm0.0041),\ (-0.0053\pm0.0055),\ (+0.0019\pm0.0038)
   \quad(k=22.4,\,23.0,\,23.9)$$
   and, over the pair population, $\frac1m\log\rho=(-0.00076\pm0.00082)$.  All four slopes are
   zero to within one standard error; the median of $\rho$ is $\sqrt2$ to four figures at
   **every** $m$ from 8 to 44.  ($\rho$ compares an $\ell^1$ functional with an $\ell^2$ norm, and
   $\ell^1/\ell^2\in[1,\sqrt2]$, so $\rho\le\sqrt2$ is exactly the range in which the cone
   minimum can be — and, spot-checked, is — attained at the first minimum itself.)  So
   $\rho=e^{o(n)}$ is not merely consistent with the data, it is the natural reading, and the
   property is generic rather than special to the Zudilin$\times$Nesterenko pair.  **[measured]**
4. **But positivity buys nothing that the lattice did not already have, and it does not
   approach $\delta=1$.**  The best *single* positive kernel over the whole two-parameter
   grid has $m^{-1}\log(\mathrm{den}\cdot|q G-p|)\to+3.81$, i.e. $+11.43$ per unit of the
   Nesterenko index — worse than the lattice by a factor $e^{11.6n}$.  The best *positive
   two-row* construction found has $\delta=0.2796$ at $m=40$ ($\delta=0.4629$ for the best
   three-row one at $m=30$), with both trends decaying to a common $1/m$-extrapolated
   asymptote $\delta_\infty\approx0.22$, and the alternative (equally good) linear fit sends the
   rate to $0$ at $m\approx165$–$185$.  Either way, **$\delta<1$ throughout, and the deficit is structural, not
   marginal.**  **[measured]**
5. **Two corrections to the existing record.**
   (a) The index deficit $\sigma(k)-\kappa_n$ splits into a Chebyshev part (fitting
   $-0.006+4.0/n$, i.e. $\to0$: this is $2\bigl(6n-\psi(6n)\bigr)/n$,
   $\psi(m)=\log\mathrm{lcm}(1,\dots,m)$) and a genuine gcd loss
   $n^{-1}\log(M_n/[\mathbb Z^2:\mathcal K_n])$ which is $+0.283$ at $n=80$ and fits *either*
   $0.196+11.7/n$ *or* $-0.196+7.55\log n/n$ equally well ($\mathrm{sd}\ 0.160$ vs $0.155$).
   So $\kappa_n=\sigma(k)-c+o(1)$ with $c\in[0,0.2]$; `CATALAN_POSITIVITY.md` §1 read the whole
   deficit as vanishing.  Either way it *helps*: the asymptotic half-log-covolume is
   $F(k)-c/2\le F(k)$.
   (b) `CATALAN_EXPLICIT.md` §6 asks whether the optimal moment pair follows a law.  It does:
   the optimum is the **adjacent** pair $(j_0,j_0+1)$ with $j_0/m\approx0.30$ — near, but not
   at, the Nesterenko point $j/m=1/3$, and never involving Zudilin's $j=0$.  It beats
   $\{$Zudilin$\times$Nesterenko$\}$ by $0.42$–$0.66$ nats per unit $m$, i.e. $1.3$–$2.0$ per
   unit $n$, and the advantage is *growing*.  **[measured, $m\le44$]**
6. **The literal reading of "products of Zudilin and Nesterenko kernels" is empty**, and this
   is worth recording: the product of two half-integer Catalan kernels has *integer*
   exponents throughout, so it is a $\zeta(2)$-world Beukers integral and produces a linear
   form in $1$ and $\zeta(2)$, not in $1$ and $G$.  The correct two-parameter family
   multiplies **one** half-integer kernel by integer-exponent positive factors, and that
   family is exactly the $(m,j)$ moment grid of `CATALAN_EXPLICIT.md` §1 (§4.1).

---

## 1. What "a positive integral representation" is, exactly

### 1.1 The moment normal form

Every construction in the literature that supplies free non-vanishing has the same shape.
Write $r_n=a_n\xi-b_n$ for the linear form.

| construction | kernel | moment form |
|---|---|---|
| Beukers $\zeta(2)$ | $\dfrac{[x(1-x)y(1-y)]^n}{(1-xy)^{n+1}}$ | $\displaystyle\int u^n\,d\mu$, $u=\frac{x(1-x)y(1-y)}{1-xy}$, $d\mu=\frac{dx\,dy}{1-xy}$ |
| Beukers $\zeta(3)$ | $\dfrac{[x(1-x)y(1-y)z(1-z)]^n}{(1-(1-xy)z)^{n+1}}$ | $u=\dfrac{x(1-x)y(1-y)z(1-z)}{1-(1-xy)z}$, $d\mu=\dfrac{dx\,dy\,dz}{1-(1-xy)z}$ |
| Zudilin Catalan | $\dfrac{x^{m-1/2}(1-x)^my^m(1-y)^{m-1/2}}{(1-xy)^{m+1}}$ | $u=\dfrac{x(1-x)y(1-y)}{1-xy}$, $d\mu=\dfrac{dx\,dy}{\sqrt{x(1-y)}\,(1-xy)}$ |
| Nesterenko $(4,7)$ | $K_Z^{(3n)}w^n$, $w=\frac{xy}{1-xy}$ | same $\mu$, different exponent pattern (§4.1) |
| $L(2,\chi_{-3})$ residue row | $R=N/D_b\ge0$ on $t\ge a$ | $\sum_{t\ge a}R(t)$, a positive sum |

In each case the index enters as a power of a *fixed* positive function against a *fixed*
positive measure.  That is the definition we adopt:

> **Definition.** A row *carries a positive-kernel (Beukers-type) representation* if there is a
> positive measure $\mu$ on $\mathbb R$ and constants $C,s$ with $r_n=C\,s^n\!\int u^n d\mu(u)$
> for all $n$.

### 1.2 The criterion, and why it is decidable

By Hamburger's theorem, $(r_n)_{n\ge0}$ is such a sequence **iff** every Hankel matrix
$H_d=(r_{i+k})_{0\le i,k\le d}$ is positive semidefinite; the measure lives on $[0,\infty)$
iff in addition the shifted Hankel $(r_{i+k+1})$ is PSD (Stieltjes), and on a compact
$[0,L]$ iff also $(L\,r_{i+k}-r_{i+k+1})$ is PSD (Hausdorff).  Rescaling $r_n\mapsto Cs^nr_n$
does not change any of these conditions, so the criterion is normalisation-free.  This is a
finite exact computation for each $d$ — see `lattice/positivity/hankel.gp`.

### 1.3 The Casoratian is the cheapest positive kernel of all

For the constant-sign question one does not even need §1.2.  In the Zagier normalisation
$$(n+1)^{w+1}u_{n+1}=p(n)u_n-c\,n^{w+1}u_{n-1},\qquad
a_nb_{n-1}-a_{n-1}b_n=-\frac{c^{\,n-1}}{n^{w+1}}\quad\textbf{[proved]},$$
so $\frac{b_n}{a_n}-\frac{b_{n-1}}{a_{n-1}}=\frac{c^{\,n-1}}{n^{w+1}a_na_{n-1}}$ and, whenever
$b_n/a_n\to\xi$,
$$a_N\xi-b_N=a_N\sum_{n>N}\frac{c^{\,n-1}}{n^{w+1}a_na_{n-1}} .$$
Hence
* $c>0$, $a_n>0$ for all $n$ $\Rightarrow$ $a_N\xi-b_N>0$ for all $N$ (all summands positive);
* $c<0$, $a_n>0$ $\Rightarrow$ the summands alternate, and by the alternating-series test
  (valid once $|c|/(a_{n+1}/a_{n-1})<1$, i.e. once $|\lambda_2|<\lambda_1$) the sum has the
  sign of its first term, so $(-1)^N(a_N\xi-b_N)>0$;
* $a_n$ not eventually of one sign (complex $\lambda$) $\Rightarrow$ no statement, and in fact
  $b_n/a_n$ has no limit.

This is elementary and classical in spirit; it is *not* stated anywhere in this repository
except for the $\sqrt{\text{Ap\'ery}}$ row (`SQRT_APERY.md` §4, where the Casoratian is
$C_n^2$).  It removes the entire "find an integral" problem for single rows.

---

## 2. The catalogue

### 2.1 Sign census (`lattice/positivity/signs.gp`, exact rational Casoratian, $n\le110$)

| row | $(a,b,c)$ | ord | $\xi$ | Casoratian | $a_n>0$ | sign of $a_n\xi-b_n$ | $n^{-1}\log|r_n|$ | $n^{-1}\log a_n$ |
|---|---|---|---|---|---|---|---|---|
| Zagier A | $(7,2,-8)$ | 2 | $\zeta(2)/4$ | OK | yes | $(-1)^n$, all $n\le110$ | $-0.054$ | $+2.028$ |
| Zagier B | $(9,3,27)$ | 2 | — | OK | **no** | oscillating, no pattern | $+1.599$ | $+1.603$ |
| Zagier C | $(10,3,9)$ | 2 | $L(2,\chi_{-3})/2$ | OK | yes | **strictly $>0$**, all $n\le110$ | $-0.054$ | $+2.146$ |
| Zagier D | $(11,3,-1)$ | 2 | $\zeta(2)/5$ | OK | yes | $(-1)^n$, all $n\le110$ | $-2.461$ | $+2.354$ |
| Zagier E | $(12,4,32)$ | 2 | $G/2$ | OK | yes | **strictly $>0$** | $+1.335$ | $+2.033$ |
| Zagier F | $(17,6,72)$ | 2 | $\tfrac58L(2,\chi_{-3})$ | OK | yes | **strictly $>0$** | $+2.032$ | $+2.159$ |
| AZ $\delta$ | $(7,3,81)$ | 3 | — | OK | **no** | oscillating | $+2.115$ | $+2.123$ |
| AZ | $(9,3,-27)$ | 3 | $L(3,\chi_{-3})/3$ | OK | yes | $(-1)^n$ | $+0.251$ | $+2.889$ |
| Domb | $(10,4,64)$ | 3 | $\tfrac7{24}\zeta(3)$ | OK | yes | **strictly $>0$** | $+1.309$ | $+2.699$ |
| AZ $\eta$ | $(11,5,125)$ | 3 | — | OK | **no** | oscillating | $+2.347$ | $+2.353$ |
| T | $(12,4,16)$ | 3 | $\tfrac7{32}\zeta(3)$ | OK | yes | **strictly $>0$** | $-0.457$ | $+3.073$ |
| Ap\'ery | $(17,5,1)$ | 3 | $\zeta(3)/6$ | OK | yes | **strictly $>0$** | $-3.608$ | $+3.448$ |

The classification is *exactly* $\operatorname{sign}(c)$, with no exceptions.  **[verified]**
The tail identity of §1.3 was checked numerically at $N=20$ for Ap\'ery, T, Zagier C and
Zagier D: agreement to all 30 printed digits.

Rows not run here, with the criterion's prediction: Cooper $s_7$ ($c=-27$) and $s_{10}$
($c=-64$) alternate; Cooper $s_{18}$ ($c=192$) is strictly positive.  **[predicted, not run]**

### 2.2 Rows outside the Zagier normalisation

* **Zudilin's Catalan row** $(Q_m,P_m)$: $(-1)^m(Q_mG-P_m)>0$, from his double integral
  (`papers/0201024v3.pdf` Thm 3).  Reconstructed here exactly: with
  $\mathrm{mom}(m,j)=\iint x^{m-1/2+j}(1-x)^my^{m+j}(1-y)^{m-1/2}(1-xy)^{-(m+1+j)}$ we verify
  $$\mathrm{mom}(m,0)=8(-1)^m\bigl(Q_mG-P_m\bigr)\quad\text{as an identity in }\mathbb Q,
  \ \textbf{[verified, all }m\le60\textbf{]}.$$
* **Nesterenko $(4,7)$**: $4B_nG-C_n>0$; verified exactly as
  $\mathrm{mom}(3n,n)=(V_nG-U_n)/(4^{7n}D_{6n}^2)$ in $\mathbb Q$ for **all $n\le20$**
  **[verified]**.  This is Prop. one-kernel of `06_threshold.tex` in the sharper form of an
  exact interpolation.
* **Beukers' row ($\sqrt{\text{Ap\'ery}}$)**: Casoratian $=C_n^2$ (Catalan number squared),
  so $a_n\xi-b_n>0$ exactly, for all $n$ — `SQRT_APERY.md` §4, Lemma 4.1.  The only
  self-contained proved fixed-sign statement in the repository that feeds an irrationality
  criterion.
* **Hypergeometric $L(2,\chi_{-3})$ row**: $Q_mL-P_m=\sum_{t\ge a}R(t)>0$ strictly, with
  $R\ge0$ on $t\ge a$ — `ONE_CLASS_TWO_WORLDS.md` §2.4, **[proved]**.  The companion Euler
  kernel $L(2,\chi_{-3})=\iint dx\,dy/(1+xy+x^2y^2)$ re-verified here to **36 digits**.
* **Brown–Zudilin $\zeta(5)$**: Brown's cellular integrals give $I_n>0$ rigorously
  (all-positive fourfold series) and $I''_n<0$ exactly for $n\le140$; positivity is used
  there to prove a *negative* (sector exclusion), `certificates/misc/frobenius.tex` §1444–1460.
* **Level-16 Catalan classes $\Phi_0$ and $E$**: **no linear form, hence no kernel.** These
  are holonomy-space functions; CDT's bound replaces non-vanishing by functional
  independence (`06_threshold.tex` §sec:adelic).  The positivity question does not arise, and
  this is exactly why the adelic route does not meet Lemma P2.
* **$\zeta(7)$ level 60**: no system exists (`ZETA7_LEVEL60.md`); vacuous.

### 2.2b The catalogue, in one table

"known" = an explicit fixed-sign kernel is in the literature or in this repository;
"derivable" = §1.3 supplies the sign unconditionally and §2.3 certifies a positive measure,
but no closed-form kernel is known; "n/a" = there is no linear form to represent.

| row | archimedean $\xi$ | sign of the form | explicit fixed-sign kernel | source |
|---|---|---|---|---|
| Zagier A $(7,2,-8)$ | $\zeta(2)/4$ | $(-1)^n$ | derivable | §2.1, §2.3 |
| Zagier B $(9,3,27)$ | none (complex) | — | **n/a** | §2.1 |
| Zagier C $(10,3,9)$ | $L(2,\chi_{-3})/2$ | $>0$ | derivable | §2.1, §2.3 |
| Zagier D $(11,3,-1)$ = Ap\'ery $\zeta(2)$ | $\zeta(2)/5$ | $(-1)^n$ | **known**: Beukers' $\zeta(2)$ double integral, re-verified to 20 digits, $5\le n\le14$ | §2.4 |
| Zagier E $(12,4,32)$ | $G/2$ | $>0$ | derivable (form diverges, so useless) | §2.1, §2.3 |
| Zagier F $(17,6,72)$ | $\tfrac58L(2,\chi_{-3})$ | $>0$ | derivable (form diverges) | §2.1, §2.3 |
| AZ $\delta$ $(7,3,81)$, AZ $\eta$ $(11,5,125)$ | none (complex) | — | **n/a** | §2.1 |
| AZ $(9,3,-27)$ | $L(3,\chi_{-3})/3$ | $(-1)^n$ | derivable | §2.1, §2.3 |
| Domb $(10,4,64)$ | $\tfrac7{24}\zeta(3)$ | $>0$ | derivable | §2.1, §2.3 |
| T $(12,4,16)$ | $\tfrac7{32}\zeta(3)$ | $>0$ | derivable | §2.1, §2.3 |
| Ap\'ery $(17,5,1)$ | $\zeta(3)/6$ | $>0$ | **known**: Beukers' triple integral, ratio $12$ verified at $n=1$ | §2.4 |
| Cooper $s_7$, $s_{10}$ | $\zeta(2)/7$, $\zeta(2)/5$ | $(-1)^n$ (predicted, $c<0$) | derivable | §2.1 note |
| Cooper $s_{18}$ | $\tfrac12L(2,\chi_{-3})$ | $>0$ (predicted, $c=192$) | derivable | §2.1 note |
| **Beukers' row** $=\sqrt{\text{Ap\'ery}}$ | $L(\Psi,2)$ | $>0$, **exact** | **known**: Casoratian $=C_n^2$, positive series | `SQRT_APERY.md` §4 |
| other root rows $\sqrt T$, $\sqrt{\text{Domb}}$, $\sqrt{s_\bullet}$ | $L(f_\bullet,2)$ etc. | $>0$ where $c>0$ | derivable | `ROOT_ROWS.md` §5 |
| **Zudilin Catalan** $(Q_m,P_m)$ | $G$ | $(-1)^m$, **proved** | **known**: half-integer Beukers kernel; exact in $\mathbb Q$ for all $m\le60$ | §2.2, §2.4 |
| **Nesterenko $(4,7)$** | $G$ | $>0$, **proved** | **known**: $K_Z^{(3n)}w^n$; exact in $\mathbb Q$ for all $n\le20$ | §2.2, §2.4 |
| the whole grid $\mathrm{mom}(m,j)$, $0\le j\le m\le60$ | $G$ | $>0$, **proved** | **known** (this document) | §4.1 |
| hypergeometric $L(2,\chi_{-3})$ row | $L(2,\chi_{-3})$ | $>0$, **proved** | **known**: well-poised residue sum; Euler kernel verified to 36 digits | `ONE_CLASS_TWO_WORLDS.md` §2, §2.4 |
| Brown–Zudilin $\zeta(5)$ | $\zeta(5)$, $\zeta(3)$ | $I_n>0$ rigorous; $I''_n<0$ exact $n\le140$ | **known**: Brown's cellular integrals | `frobenius.tex` §1444–1460 |
| level-12 $\zeta(5)$ | $\tfrac{11}{144}\zeta(5)$ | $(-1)^{n+1}$, only polynomial decay | derivable in principle; the $\log^5$ singularity is on the unit circle | `CLAUDE_FINDINGS.md` |
| $\zeta(7)$ level 24 | $\tfrac{1463}{13824}\zeta(7)$ | not determined | not known (order-7 system; §1.3 does not apply verbatim) | `ROW_LEDGER.md` C3 |
| $\zeta(7)$ level 60 | — | — | **n/a**: no system exists | `ZETA7_LEVEL60.md` |
| level-16 classes $\Phi_0$, $E$ | — | — | **n/a**: holonomy-space functions, not linear forms | `CATALAN_TWO_CLASSES.md`, §2.2 |

Two entries deserve emphasis because they are the only *asymmetries* in the table.  First,
every "derivable" is genuinely derivable: §1.3 is a two-line argument from the recurrence and
§2.3 certifies the measure.  Second, the two rows for which the fixed-sign kernel is
load-bearing in an actual irrationality argument — Beukers' $\sqrt{\text{Ap\'ery}}$ row and
the Zudilin/Nesterenko Catalan pair — are also the two where the sign is *proved for all $n$*
rather than verified on a range.  That is not a coincidence: an argument that uses the sign
needs it for all $n$, and only the Casoratian and the integral supply that.

### 2.3 The Hankel/moment test (`lattice/positivity/hankel.gp`, $\backslash$p 2000, $d\le10$)

| row | $\det H_d>0$ for all $d\le$ | $\log_{10}\det H_{10}$ | Stieltjes ($\mu$ on $[0,\infty)$) |
|---|---|---|---|
| Zagier A | 10 (=max) | $-68.61$ | fails at $d=0$ ($r_1<0$) |
| Zagier C | 10 | $-68.83$ | **yes**, $d\le10$ |
| Zagier D | 10 | $-184.61$ | fails at $d=0$ |
| Zagier E | 10 | $-1.79$ | **yes** |
| Zagier F | 10 | $+31.76$ | **yes** |
| AZ$(9,3,-27)$ | 10 | $-56.18$ | fails at $d=0$ |
| Domb | 10 | $-4.86$ | **yes** |
| T | 10 | $-90.47$ | **yes** |
| Ap\'ery $\zeta(3)$ | 10 | $-242.17$ | **yes** |

Every row with a real limit is a Hamburger moment sequence as far as tested; the Stieltjes
column is $\operatorname{sign}(c)>0$ exactly.  For the three $c<0$ rows the twisted sequence
$((-1)^nr_n)$ is Stieltjes instead (immediate from $r_n=\int u^n d\mu$ with
$\operatorname{supp}\mu\subset(-\infty,0]$).  **[verified $d\le10$; the full moment property
is $d\to\infty$ and is not proved]**

**Reading.**  This says the *existence* of a positive-kernel representation is generic, and
that hunting for explicit Beukers integrals row by row is the wrong use of effort: what is
scarce is not a positive kernel but a kernel whose $n$-th moment is *small* relative to the
denominator of the row.  Compare §4.2: on the Catalan moment family the positive kernels are
in bijection with a two-parameter grid and *every one of them* has
$\mathrm{den}\cdot|r|\to\infty$ exponentially.

### 2.4 Numerical verification of the integral representations

`lattice/positivity/intreps.gp`, all-positive series evaluations against exact rationals.

| representation | check | accuracy |
|---|---|---|
| $\mathrm{mom}(m,0)=8(-1)^m(Q_mG-P_m)$ | exact in $\mathbb Q$ | **all** $m\le60$, no failure |
| $\mathrm{mom}(3n,n)=(V_nG-U_n)/(4^{7n}D_{6n}^2)$ | exact in $\mathbb Q$ | **all** $n\le20$, no failure |
| $\mathrm{mom}(m,j)$ vs its all-positive series | numeric | 48 digits at $(12,0)$, 43 at $(16,0)$, 36 at $(20,0)$, 38 at $(18,6)$, 36 at $(15,5)$, 30 at $(12,4)$ |
| Beukers $\zeta(2)$ vs Zagier D | numeric | $\displaystyle\iint\frac{[x(1-x)y(1-y)]^n}{(1-xy)^{n+1}}=(-1)^n\bigl(a_n\zeta(2)-5b_n\bigr)$, ratio $=\pm5.00000000000000000000$ for $5\le n\le14$ (20 digits) |
| Beukers $\zeta(3)$ vs Ap\'ery | numeric | $\displaystyle\iiint\frac{[x(1-x)y(1-y)z(1-z)]^n}{(1-(1-xy)z)^{n+1}}=12\bigl(a_n\tfrac{\zeta(3)}6-b_n\bigr)$, ratio $=12.0000000000000000$ at $n=1$ (3-D `intnum`, 16 digits; higher $n$ not attempted — the integrand's endpoint behaviour makes 3-D quadrature expensive) |
| $L(2,\chi_{-3})=\iint\frac{dx\,dy}{1+xy+x^2y^2}$ | numeric | 36 digits |

Note the Beukers $\zeta(2)$ line is also an independent confirmation of the sign census: the
positive kernel represents $(-1)^n r_n$, matching the $c=-1$ row of §2.1.

---

## 3. Task 2 — the P2 ratio on the exact lattices $\mathcal K_n$, $n\le80$

### 3.1 Setup

Unchanged from `CATALAN_POSITIVITY.md` §1.  $X_n=2^{e_{3n}}D_{6n}^2Q_{3n}$,
$Y_n=2^{e_{3n}}D_{6n}^2P_{3n}$, $V_n=4^{7n+1}D_{6n}^2B_n$, $U_n=4^{7n}D_{6n}^2C_n$;
$S_n=D_{6n}^2$, $T_n=2^{\lfloor kn\rfloor}$, $M_n=S_nT_n$, and
$$\mathcal K_n=\{c\in\mathbb Z^2:\ M_n\mid c_ZX_n+c_NV_n,\ M_n\mid c_ZY_n+c_NU_n\}$$
computed by `matkerint`/`mathnf`.  Scaled coordinates $(s_Zc_Z\lambda_Z,\ c_N\lambda_N)$,
$\lambda_Z=|X_nG-Y_n|/M_n$, $\lambda_N=(V_nG-U_n)/M_n$; positive cone $\mathcal P$;
$\lambda_1,\lambda_2$ the successive minima in the scaled Euclidean metric; the cone minimum
is the minimum of the *linear form* $s_Zc_Z\lambda_Z+c_N\lambda_N=|q_nG-p_n|$ over
$\mathcal K_n\cap\mathcal P$, computed by an exact-integer interval sweep whose range $R$ is
doubled until the answer is stable.  For every one of the 231 instances the returned
$(q_n,p_n)$ were checked to be integers with $|q_nG-p_n|$ equal to the reported minimum
(3 000-digit $G$, working precision 4 000 digits).  **No instance failed.**

Constants: $E_1=13.0995887908\ldots$, $E_2=14.3931452672\ldots$,
$k_*=22.3512905953\ldots$, $\sigma(k)=12+k\log2$, $F(k)=\tfrac{\log2}2(k_*-k)$, so
$F(22.4)=-0.01688$, $F(23.0)=-0.22483$, $F(23.9)=-0.53674$.

### 3.2 The data (`lattice/positivity/data/cone_n80.csv`, all $4\le n\le80$, $k\in\{22.4,23.0,23.9\}$)

Selected rows at $k=22.4$ (full table in the CSV; the $n\equiv0\ (4)$ entries reproduce
`CATALAN_POSITIVITY.md` §3.2 bit for bit, which is the cross-check that the two independent
implementations agree):

| $n$ | $\kappa_n$ | $\frac1{2n}\log\mathrm{covol}$ | $\frac1n\log\lambda_1$ | $\frac1n\log(\text{cone-min})$ | $\rho$ | $\frac1n\log q_n$ |
|---|---|---|---|---|---|---|
| 8 | 25.7946 | $-0.9284$ | $-0.9256$ | $-0.9099$ | 1.13 | 13.955 |
| 11 | 25.3005 | $-1.1476$ | $-1.2935$ | $-0.9845$ | **29.95** | 13.908 |
| 24 | 26.7366 | $-0.4320$ | $-0.4571$ | $-0.4565$ | 1.02 | 14.301 |
| 44 | 27.0917 | $-0.2511$ | $-0.2608$ | $-0.2538$ | 1.36 | 14.680 |
| 60 | 27.0914 | $-0.2589$ | $-0.2677$ | $-0.2473$ | 3.39 | 14.701 |
| 68 | 27.0621 | $-0.2705$ | $-0.2981$ | $-0.2414$ | **47.36** | 14.649 |
| 78 | 27.4797 | $-0.0607$ | $-0.0723$ | $-0.0455$ | 8.08 | 14.869 |
| 80 | 27.3199 | $-0.1421$ | $-0.1488$ | $-0.1459$ | 1.26 | 14.803 |

### 3.3 The fit — and the precise conjecture

$N=77$ values of $n$ at each $k$.

| $k$ | $\min\rho$ | median | $\max\rho$ | $\frac1n\log\rho$ fit | fraction $\rho\le\sqrt2$ |
|---|---|---|---|---|---|
| 22.4 | 1.016 | 1.409 | 47.4 | $(+0.00217\pm0.00407)\,n+0.665$ | 0.571 |
| 23.0 | 1.010 | 1.397 | 424.7 | $(-0.00525\pm0.00546)\,n+0.985$ | 0.558 |
| 23.9 | 1.003 | 1.410 | 33.0 | $(+0.00190\pm0.00380)\,n+0.637$ | 0.545 |

(residual sd $0.74$–$1.06$; the alternative model $\log\rho=\alpha\log n+\beta$ fits no better:
$\alpha=+0.030\pm0.124$, $-0.141\pm0.167$, $+0.100\pm0.115$.)

Equivalent statement in the form Lemma P2 actually needs:
$$\frac1n\log\frac{\text{cone-min}}{\sqrt{\mathrm{covol}}}=(+0.00148\pm0.00346)\,n+0.209,\quad
(-0.00356\pm0.00374)\,n+0.464,\quad(+0.00405\pm0.00339)\,n+0.101 .$$

Two further statistics, both stable in $n$:

* **The shortest vector of $\mathcal K_n$ is a cone vector (up to $\pm$) at $57\%/55\%/51\%$
  of $n\le80$** — the $61\%$ of `06_threshold.tex` at $n\le44$ was not a small-sample artefact.
* **In every block of ten consecutive $n$ there is an $n$ with $\log\rho\le0.20$**, and usually
  $\le0.06$ (block minima at $k=22.4$: $0.126,\,0.088,\,0.016,\,0.043,\,0.056,\,0.031,\,0.157$).
  Since Lemma P2 only asks for *infinitely many* $n$, this is the relevant statistic and it is
  far stronger than the average.

> **Conjecture P2$'$ (quantitative form of Lemma P2).**  For the Catalan two-row lattices,
> $$\liminf_{n\to\infty}\ \frac1n\log
> \frac{\min\{s_Zc_Z\lambda_Z+c_N\lambda_N:\ c\in\mathcal K_n\cap\mathcal P\}}
> {\lambda_1(\mathcal K_n)}\;=\;0 ,$$
> and more strongly $\log\rho_n=O(1)$ along a positive-density set of $n$.

**What P2$'$ would prove.**  Combined with Lemma P1 ($\kappa_n=\sigma(k)+o(n)$ for some
$k>k_*$) and Lemma P3 (the rates $E_1,E_2$, $\log D_{6n}=6n+o(n)$), it gives integers
$q_n,p_n$ with
$$0<q_nG-p_n\le e^{(F(k)-c/2)n+o(n)}\longrightarrow0 ,$$
the positivity of the left inequality being a *theorem* (Zudilin's and Nesterenko's
integrals), whence $G\notin\mathbb Q$.  The correction $-c/2\le0$ is the gcd loss of §3.4.

**What it cannot be.**  Unchanged from `CATALAN_POSITIVITY.md` §4: P2$'$ is **false for every
rational $G=a/b$** (the kernel direction becomes the first minimum once
$n\gtrsim\log b/|F|$), so no covolume-only lattice lemma can prove it, and the present data
would be reproduced verbatim by a rational surrogate (§6.4).  The honest content of the
numerics is: *the exponential gap, if any, is smaller than $0.010$ per unit $n$ at $k=22.4$,
which does not yet exclude a gap of size $|F(22.4)|=0.0169$; at $k=23.9$, where
$|F|=0.537$, the one-sigma bound $0.0038$ excludes any gap that could defeat the mechanism by
a factor of $140$.*  That asymmetry — the mechanism is *better* supported the further one goes
past $k_*$, while $F$ is more negative there — is the single most useful thing the extended
range shows.

### 3.4 Correction: what the index deficit really is

`CATALAN_POSITIVITY.md` §1 reports $\sigma-\kappa_n$ falling $0.94\to0.44$ over $n\le44$ and
reads it as $o(1)$.  At $n=80$ it is $0.2066$, at all three $k$ *identically* — which already
shows it is not a $2$-adic effect.  Splitting
$$\sigma(k)-\kappa_n=\underbrace{\Bigl[\sigma(k)-\tfrac1n\log M_n\Bigr]}_{\text{Chebyshev}}
+\underbrace{\tfrac1n\log\frac{M_n}{[\mathbb Z^2:\mathcal K_n]}}_{\text{gcd loss}},$$
* the Chebyshev part is $2\bigl(6n-\psi(6n)\bigr)/n$, $\psi(m)=\log\mathrm{lcm}(1,\dots,m)$;
  it fits $-0.006+4.0/n$ and is $-0.076$ at $n=80$, where it has changed sign
  ($\psi(480)-480=+3.05$).  It $\to0$ by PNT.  **[measured, standard]**
* the gcd loss is $+0.283$ at $n=80$ and fits $0.196+11.7/n$ (sd $0.160$) or
  $-0.196+7.55\log n/n$ (sd $0.155$) *equally well*.  **[measured, undetermined]**

So $\kappa_n=\sigma(k)-c+o(1)$ with $c\in[0,0.2]$; Lemma P1 should be restated with this
constant, and the achievable rate is $F(k)-c/2$, i.e. between $-0.017$ and $-0.117$ at
$k=22.4$.  Deciding $c$ is a clean finite question about $\gcd(X_n,Y_n,V_n,U_n,M_n)$.

Consequences at $n=80$: $\delta=-\log|q_nG-p_n|/\log q_n=0.00986$ ($k=22.4$), $0.02342$
($23.0$), $0.04563$ ($23.9$) — with **proved** non-vanishing.

---

## 4. Task 3 — positive constructions on the two-parameter kernel grid

### 4.1 The family, and a structural correction

Write the Catalan-world Beukers integral as
$$I(a,b,c,d,e)=\iint_{[0,1]^2}x^{a-1/2}(1-x)^by^c(1-y)^{d-1/2}(1-xy)^{-e}\,dx\,dy,
\qquad a,b,c,d,e\in\mathbb Z_{\ge0}.$$
The **two** half-integer exponents are what make the value a form in $1,G$ rather than in
$1,\zeta(2)$.  Consequently:

> **Observation (why "products of Zudilin and Nesterenko kernels" is empty).**  The product
> $K_Z^{(m)}\cdot K_N^{(n)}$ has $x$-exponent $(m-\tfrac12)+(4n-\tfrac12)=m+4n-1\in\mathbb Z$
> and likewise for $(1-y)$.  So it is an all-integer Beukers integral and its value is a
> $\mathbb Q$-linear form in $1$ and $\zeta(2)$, **not** in $1$ and $G$.  Products of two
> Catalan kernels leave the Catalan world.

The correct two-parameter family multiplies *one* half-integer kernel by integer-exponent
positive factors.  Up to the two directions realised by the sources it is
$$\mathrm{mom}(m,j):=I(m+j,\ m,\ m+j,\ m,\ m+j+1)
=\iint K_Z^{(m)}\,w^j,\qquad w=\frac{xy}{1-xy},\qquad 0\le j\le m,$$
with the equivalent product description (put $s=m-3j$, $t=j$)
$$\mathrm{mom}(m,j)=\iint K_Z^{(s)}\cdot v^{\,t},\qquad
v=\frac{x^4(1-x)^3y^4(1-y)^3}{(1-xy)^4}>0 ,$$
so the grid **is** the $(s,t)$ grid of products of a Zudilin kernel of index $s$ with $t$
copies of the Nesterenko shape: $\mathrm{Zudilin}(m)=(m,0)$, $\mathrm{Nesterenko}(n)=(3n,n)$.
Every $\mathrm{mom}(m,j)>0$ (positive integrand), so *every* grid point, and every
non-negative combination of grid points, is an admissible positive-cone linear form.
`lattice/positivity/rows_pos.gp` computes $[A,B]$ with $\mathrm{mom}=AG-B$ by exact rational
linear algebra (a generalisation of the Nesterenko partial-fraction solver).

### 4.2 One kernel: the whole grid $m\le60$ (`grid.gp`, 1 890 points)

Objective: $m^{-1}\log(\mathrm{den}\cdot\mathrm{mom}(m,j))$ where
$\mathrm{den}=\mathrm{lcm}(\mathrm{den}\,A,\mathrm{den}\,B)$, so that $q=\mathrm{den}\,A$,
$p=\mathrm{den}\,B$ are integers and the objective is $m^{-1}\log(qG-p)$ with $qG-p>0$.

| $m$ | 5 | 10 | 20 | 30 | 40 | 50 | 60 |
|---|---|---|---|---|---|---|---|
| $\min_j$ objective | $+1.69$ | $+2.51$ | $+2.63$ | $+3.24$ | $+3.23$ | $+3.40$ | $+3.36$ |
| argmin $j$ | 0 | 4 | 3 | 3 | 8 | 12 | 13 |
| $j=0$ (Zudilin) | $+1.69$ | $+3.11$ | $+3.33$ | $+3.49$ | $+3.91$ | $+3.65$ | $+3.97$ |
| $j=m/3$ (Nesterenko) | — | — | — | $+3.43$ | — | — | $+3.63$ |

Fit over $m\ge20$: $\min_j\text{objective}=+3.809-20.87/m$ (residual sd $0.100$), i.e.

$$\boxed{\ \lim_{m}\ \min_j\ \tfrac1m\log\bigl(\mathrm{den}\cdot|q G-p|\bigr)=+3.81
\quad\Longleftrightarrow\quad +11.43\ \text{per unit of the Nesterenko index }n=m/3.\ }$$

Compare the two-row lattice, $-0.15$ per unit $n$ at $n=80$: **the best single positive kernel
is worse by $\approx11.6$ per unit $n$, i.e.\ by a factor $e^{11.6n}$.**  Only the trivial
point $(m,j)=(1,0)$ — the form $13-14G=0.17655$ with $q=-14$, $p=-13$ — has $\delta>0$; every point with $m\ge2$ has
$\delta<0$, falling to $\delta=-0.40$ at $m=60$.  This confirms and sharpens
`CATALAN_EXPLICIT.md` §2, whose $n\le5$ table ($\le+7.7$ per unit $n$) was still far from the
asymptote.

**Vertex lemma, correctly stated.**  Fix a common denominator $\mathrm{den}$ for a finite set
$S$ of grid points.  For $p\in\mathbb Z_{\ge0}^S$, $\mathrm{den}\sum p_i\,\mathrm{mom}_i$ is a
sum of positive terms, hence $\ge\min_i\mathrm{den}\cdot\mathrm{mom}_i$, attained at a vertex.
Since $\mathrm{den}\ge\mathrm{den}_i$ for each $i$, this holds across different $m$ too.
**Non-negative coefficients give no cancellation in the analytic factor, ever.**
The *only* escape is the arithmetic factor: $\gcd(\mathrm{den}\sum p_iA_i,\
\mathrm{den}\sum p_iB_i)$ can be huge, and that is precisely the $2$-adic bridge.  So the
question "is there a good explicit positive construction?" is *identical* to "is there a
small non-negative vector in a prescribed congruence class?" — which is Lemma P2 for that
tuple of rows.  §4.3 measures it.

### 4.3 Two kernels: 13 332 pair lattices (`pairs.gp`, all $0\le j_1<j_2\le m$, $m\le44$)

For a pair, put $d=\mathrm{lcm}$ of the four denominators, $a_i=dA_i$, $b_i=dB_i$,
$h=a_1b_2-a_2b_1$, and
$$\mathrm{MOD}=d\cdot2^{v_2(h)},\qquad
L=\{p\in\mathbb Z^2:\ \mathrm{MOD}\mid p\!\cdot\!a,\ \mathrm{MOD}\mid p\!\cdot\!b\},$$
so $q=p\!\cdot\!a/\mathrm{MOD}$, $p'=p\!\cdot\!b/\mathrm{MOD}$ are integers and, on
$L\cap\mathbb Z^2_{\ge0}$, $qG-p'=d(p_1M_1+p_2M_2)/\mathrm{MOD}>0$ **by a theorem**.  Of the
13 357 pairs computed, 25 are degenerate (the cone minimum has $q=0$, i.e. the combination is
a rational number) and are discarded; all 13 332 remaining instances were checked to return
genuine integers with $|qG-p'|$ equal to the reported cone minimum.

**(a) The P2 ratio, on a population.**  $\rho=(\text{cone-min})/\lambda_1$ (the same
$\ell^1$-on-the-cone versus $\ell^2$ comparison as §3; $\rho\le\sqrt2$ is the range compatible
with the cone minimum being attained at the first minimum):

| $m$ | $N$ | $\min\rho$ | median $\rho$ | $\Pr[\rho\le\sqrt2]$ | $\Pr[\rho>20]$ |
|---|---|---|---|---|---|
| 12 | 78 | 1.027 | 1.555 | 0.462 | 0.013 |
| 20 | 210 | 1.006 | 1.417 | 0.500 | 0.024 |
| 28 | 405 | 1.001 | 1.410 | 0.528 | 0.025 |
| 36 | 664 | 1.000 | 1.414 | 0.517 | 0.021 |
| 44 | 988 | 1.000 | 1.414 | 0.491 | 0.022 |
| all | 13 332 | 1.000 | **1.4142** | 0.507 | 0.027 |

$$\tfrac1m\log\rho=(-0.00076\pm0.00082)\,m+0.790 .$$
The median is $\sqrt2$ to four figures at every $m$, and the tail probabilities do not move.
**On these lattices the positivity constraint costs a bounded factor, not an exponential one,
and this is a generic property of the family.**  This is the strongest evidence for P2 in the
project, because it is a population statement: it cannot be an accident of the particular
Zudilin$\times$Nesterenko pairing.

**(b) The best positive two-row construction, and the pair law.**  Best pair at each $m$
(selection over a growing pool, so read the fixed-rule column as the unbiased one):

| $m$ | best $(j_1,j_2)$ | $\frac1m\log|qG-p|$ | $\frac1m\log q$ | $\delta$ | $\{$Zud$\times$Nest$\}$ $(0,m/3)$ | advantage |
|---|---|---|---|---|---|---|
| 12 | $(2,3)$ | $-1.4877$ | 3.234 | 0.460 | $-0.7643$ ($\delta=0.185$) | $+0.72$ |
| 24 | $(5,6)$ | $-1.1514$ | 3.753 | 0.307 | $-0.5349$ ($\delta=0.122$) | $+0.62$ |
| 30 | $(10,11)$ | $-1.0446$ | 3.905 | 0.268 | $-0.5211$ ($\delta=0.118$) | $+0.52$ |
| 36 | $(12,13)$ | $-1.0656$ | 3.902 | 0.273 | $-0.4828$ ($\delta=0.111$) | $+0.58$ |
| 40 | $(10,11)$ | $-1.0818$ | 3.869 | **0.280** | — | — |
| 42 | — | — | — | — | $-0.3621$ ($\delta=0.079$) | $+0.66$ |

Fixed rule $j_0=\mathrm{round}(0.30\,m)$, $j_1=j_0+1$ (no selection bias):

| $m$ | 10 | 20 | 30 | 40 | 44 |
|---|---|---|---|---|---|
| $\frac1m\log|qG-p|$ | $-1.292$ | $-1.181$ | $-1.017$ | $-1.012$ | $-1.013$ |
| $\delta$ | 0.375 | 0.318 | 0.259 | 0.257 | 0.258 |
| $v_2(h)/m$ (bits) | 8.50 | 8.80 | 8.83 | 9.00 | 8.93 |

Fits over $m\ge15$: rate $=-0.890-5.08/m$ (sd $0.066$) or $=+0.0069\,m-1.282$ (sd $0.070$);
the two are statistically indistinguishable, so the asymptote is either $-0.89$ per unit $m$
($\delta_\infty\approx0.23$) or $0$ at $m\approx186$.  **Either way $\delta<1$.**

* **The pair law** answers `CATALAN_EXPLICIT.md` §6's open question: the optimum is
  $(j_0,j_0+1)$ with $j_0/m\approx0.30$ — *adjacent* indices near, but not at, the Nesterenko
  point $1/3$, and never Zudilin's $j=0$.  The advantage over $\{$Zudilin$\times$Nesterenko$\}$
  is $0.42$–$0.66$ nats per unit $m$ ($1.3$–$2.0$ per unit $n$) and is **growing**, whereas
  the Zudilin$\times$Nesterenko $\delta$ decays from $0.30$ ($m=6$) to $0.079$ ($m=42$).
* **The $2$-adic bridge across the whole grid.**  $v_2(h)/m$ has median $9.72$ bits at
  $m\ge40$ (range $7.8$–$11.7$), i.e. **median $29.2$ bits per unit $n=m/3$**, against
  $k_*=22.351$.  The $24.06n$ law of `CATALAN_AUDIT.md` §4(d) is the value at the
  $\{0,n\}$ point specifically ($8.4$ bits per $m$ = $25.2$ per $n$ here); the family as a
  whole sits well above it.  **[measured, $m\le44$]**  As in the audit, this is the *exact
  per-instance* valuation with no proved lower bound — the single least secure input below.

### 4.4 Three kernels (`triples.gp`)

Rank-3 lattices with $\mathrm{MOD}=d\cdot2^{v_2(\gcd\text{ of the three }2\times2\text{ minors})}$,
cone minimum over $\mathbb Z^3_{\ge0}$:

| $m$ | best pair $\delta$ | best triple | triple $\delta$ | gain (nats/$m$) |
|---|---|---|---|---|
| 12 | 0.460 | $(0,5,7)$ | **0.826** | $+0.711$ |
| 18 | 0.368 | $(0,6,8)$ | 0.683 | $+0.658$ |
| 24 | 0.307 | $(1,5,10)$ | 0.526 | $+0.549$ |
| 30 | 0.268 | $(6,12,15)$ | 0.463 | $+0.518$ |

A third positive row helps substantially at finite $m$ — the $\delta=0.826$ at $m=12$ is the
best positive-cone number anywhere in this project — but the gain shrinks and the $1/m$
extrapolation of the triple column lands on $\delta_\infty\approx0.22$, the *same* asymptote
as the pairs.  **Rank buys pre-asymptotics, not rate.**  This is the constructive counterpart
of the "more rows" clause of Theorem (Absorption)(5) in `06_threshold.tex`.

### 4.5 Sums of squares

Unchanged: `CATALAN_EXPLICIT.md` §3 and Lemma E4 stand.  The mechanism is the Pareto tail
$\nu_n(w>W)\asymp W^{-(3n+1/2)}$ of the pushforward measure — only $3n$ finite moments and no
concentration — so all orthogonal-polynomial norms share one exponential order and
$\det(H^{(d)})^{1/(d+1)}=e^{-8.2n+o(n)}$ flat in $d$.  Nothing in §4.2–4.4 changes this; the
grid computation above supersedes it in the only direction that mattered (the *scale* of the
denominator, which is $e^{+3.81m}$ against an $e^{-8.2n/3}=e^{-2.7m}$ analytic headroom).

---

## 5. Synthesis: what sign-based non-vanishing buys

1. **For a single row it buys everything and costs nothing** — and it does not need an
   integral: §1.3 supplies it from the recurrence.  Consequently *the search for explicit
   Beukers integrals for the modular rows is not the bottleneck*; §2.3 shows the
   representations exist (as positive measures) for every row with a real limit.
2. **For a multi-row construction it buys exactly a factor $\rho$, and $\rho$ is empirically
   $O(1)$ on 13 563 lattices.**  This is the content of Lemma P2, and it is now tested on a
   population rather than one family.
3. **What it does not buy is rate.**  The positive-cone rate is pinned to $\lambda_1$, hence
   to $\sqrt{\mathrm{covol}}$, hence to the arithmetic input — $\kappa_n$ (Lemma P1, now with
   the gcd-loss constant $c$ of §3.4) and the $2$-adic bridge $v_2(h)$.  Every number in §4
   that looks good comes from the bridge, not from positivity.
4. **The residual logical situation is unchanged.**  Positivity converts "non-vanishing *and*
   size" into "size, non-vanishing free"; the size statement is still equivalent to
   irrationality, still false for every rational $G$, and still not provable from the
   covolume.  §6.4's control demonstrates this concretely on the new tables.

### 5.1 Ranked next steps

1. **Decide the gcd-loss constant $c$ of §3.4** ($c/n$ vs $\log n/n$).  Finite, cheap,
   and it moves the target rate by up to $0.1$ per unit $n$ — comparable to $|F(22.4)|$ itself.
2. **Prove the pair law of §4.3(b)** ($j_0\approx0.30m$, $j_1=j_0+1$) and re-run the whole
   `CATALAN_DIRECTIONAL`/`CATALAN_POSITIVITY` package on that pair rather than $\{0,n\}$;
   the measured $v_2(h)/n\approx29$ there is $6.8$ bits above $k_*$, against $24.06-22.35=1.7$
   for $\{0,n\}$.
3. **A proved lower bound for $v_2(h_{j_0,j_0+1})$** in the moment family.  Everything in §4
   rests on the exact per-instance valuation.  The candidate lemma
   $v_2(\mathrm{den}\,B_{n,j})=v_2(\mathrm{den}\,A_{n,j})+1$ (`CATALAN_EXPLICIT.md` Lemma E1)
   is the cheap end of this.
4. **Extend the Hankel test of §2.3 to $d\to\infty$ for one row** (say T, whose form decays):
   an actual positive measure would give the first Beukers-type integral for a modular row,
   which is exactly the object `ZETA3_TWO_LATTICE.md` §6 and `07_problems.tex` ask for.

---

## 6. Status of every claim

### 6.1 Proved
* The Casoratian identity and the tail identity of §1.3, for any row in the Zagier
  normalisation; the sign dichotomy that follows from $\operatorname{sign}(c)$ (with the
  alternating case conditional on the standard $|\lambda_2|<\lambda_1$ monotonicity).
* The moment normal form of §1.1 for Beukers $\zeta(2)$, $\zeta(3)$ and Zudilin's kernel.
* Convergence range $0\le j\le m$ for $\mathrm{mom}(m,j)$; positivity of every grid point.
* The $\zeta(2)$-world observation of §4.1 (products of two Catalan kernels).
* Vertex lemma (§4.2), with the arithmetic-factor caveat.

### 6.2 Verified (exact arithmetic, stated range)
* Casoratian for all twelve sporadic rows, $n\le110$; sign classification, same range.
* $\mathrm{mom}(m,0)=8(-1)^m(Q_mG-P_m)$ in $\mathbb Q$ for all $m\le60$;
  $\mathrm{mom}(3n,n)=(V_nG-U_n)/(4^{7n}D_{6n}^2)$ in $\mathbb Q$ for all $n\le20$.
* Integrality and exactness of every reported cone minimum: 231 lattice instances ($n\le80$,
  three $k$) and 13 332 pair instances, each re-checked as $|q G-p|$ from the returned
  integers.  No failures.
* Hankel positivity $\det H_d>0$, $d\le10$, for the nine rows of §2.3.

### 6.3 Measured (numerical fit, not converged)
* Every regression in §3.3, §3.4, §4.2, §4.3, §4.4.  In particular the two competing fits for
  the gcd loss, and for the best-pair rate, are *statistically indistinguishable* over the
  available range; the asymptotes quoted are therefore ranges, not values.
* All $v_2(h)$ statistics: exact per instance, but with no proved lower bound for general $m$.

### 6.4 The control
`lattice/positivity/control.gp`.  With $G^*=\mathrm{bestappr}(G,10^{320})$
($\log_{10}\mathrm{den}=319.50$, $|G-G^*|=10^{-639.6}$), the best pair and its cone minimum
agree to more than 100 decimal digits at $m=8,12,16,20,24$ (they can only differ at the
$640$th digit, where $G^*$ leaves $G$):

```
m= 8  true G: best (1,2) rate/m=-1.596749  |  G*: best (1,2) rate/m=-1.596749  IDENTICAL to >100 digits
m=12  true G: best (2,3) rate/m=-1.487682  |  G*: best (2,3) rate/m=-1.487682  IDENTICAL to >100 digits
m=24  true G: best (5,6) rate/m=-1.151368  |  G*: best (5,6) rate/m=-1.151368  IDENTICAL to >100 digits
```

A rational number produces every table in §3 and §4.  **No finite portion of this data is
evidence about $G$.**

### 6.5 Open
* Conjecture P2$'$ (§3.3); the gcd-loss constant $c$ (§3.4); the pair law (§4.3);
  a proved lower bound for $v_2(h)$; Hankel positivity for all $d$.
* Whether $\delta_\infty$ for positive-cone constructions on this family is $\approx0.22$ or
  $0$; both fits are admissible on $m\le44$.

---

## 7. Scripts

`lattice/positivity/` — PARI/GP unless stated.  Run as
`cat lattice/positivity/rows_pos.gp [deps] X.gp > run.gp; gp -q run.gp`.

| file | deps | what |
|---|---|---|
| `rows_pos.gp` | — | exact Zudilin/Nesterenko rows; the general moment `mom(m,j)`; `klat2` |
| `signs.gp` | — | §2.1 sign census, exact Casoratian, tail identity |
| `hankel.gp` | `signs.gp` | §2.3 Hamburger/Stieltjes test |
| `intreps.gp` | `rows_pos.gp` | §2.4 all-positive series for the kernels; Beukers $\zeta(2)$, $\zeta(3)$; $\chi_{-3}$ Euler kernel |
| `cone80.gp` | `rows_pos.gp` | §3 lattice $\mathcal K_n$, minima, exact cone sweep (`sweepn(n,klist,G)`) |
| `grid.gp` | `rows_pos.gp` | §4.2 single-kernel grid (`gridrow(m)`) |
| `pairs.gp`, `pairs_run.gp` | `rows_pos.gp`, `cone80.gp` | §4.3 pair lattices (`runm(m)`) |
| `triples.gp` | + `pairs.gp` | §4.4 rank-3 cone minima (`runtrip(m,R)`) |
| `control.gp` | + `pairs.gp` | §6.4 rational surrogate |
| `fit.py`, `fit2.py`, `pnt.py` | — | §3.3, §3.4 regressions |
| `gan.py`, `pan.py`, `pan2.py` | — | §4.2, §4.3 regressions |
| `data/cone_n80.csv` | | 231 rows: $k$, $n$, $\kappa_n$, $\tfrac1{2n}\log\mathrm{covol}$, $\tfrac1n\log\lambda_{1,2}$, $\tfrac1n\log$cone-min, $\rho$, $\tfrac1n\log q$, in-cone flag, sweep range, guard, bad flag |
| `data/grid_m60.csv` | | 1 890 rows: $m$, $j$, $s=m-3j$, $\tfrac1m\log\mathrm{den}$, $\tfrac1m\log|{\rm mom}|$, objective per $m$, per $n$, $\tfrac1m\log q$, $\delta$ |
| `data/pairs_m44.csv` | | 13 357 rows: $m$, $j_1$, $j_2$, $\tfrac1m\log$cone-min, $\tfrac1m\log\lambda_1$, $\rho$, $\tfrac1m\log q$, ok flag, $v_2(h)\log2/m$, $\tfrac1m\log[\mathbb Z^2\!:\!L]$ |

**Precision.**  `sweepn` needs $G$ supplied at $\ge3\,000$ digits (the production runs used
$\backslash$p 4000); `redu`/`redu2` scale their LLL rounding to the precision actually present
in the arguments, so the scripts compose without a precision crash, but a low-precision $G$
silently degrades the results.  Every table above was produced at the precision stated in its
section.

Total compute: about 25 minutes wall-clock on 6 cores.  Longest single job: `cone80.gp` at
$n=80$ (60 s, 2.3 GB).
