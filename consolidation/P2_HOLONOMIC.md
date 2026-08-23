# Is $\mathcal K_n$ holonomic?  The Hermite data of the Catalan congruence lattices, and the adelic form of P2$'$

*Fable, 2026-08-23.  Scripts: `lattice/p2_holonomic/`, data in `lattice/p2_holonomic/data/`.
Companion to `P2_STRUCTURE.md` (which reduced P2$'$ to the parity of a continued-fraction
index); this document asks whether the arithmetic that builds $\mathcal K_n$ controls that
continued fraction.  Tags: **[proved]**, **[verified]** = exact integer computation over a
stated range, **[measured]**, **[open]**.  All lattice arithmetic is exact in $\mathbb Z^2$;
$G$ enters only through two real weights.*

**No claim of irrationality is made anywhere below.**

---

## 0. Verdict

1. **The lattice is completely $G$-free, and the continued fraction is an exact modular-inverse
   object.**  In oriented coordinates $\mathcal K_n=\mathbb Z(h_{11},0)+\mathbb Z(h_{12},h_{22})$ with
   $$h_{11}=\frac{M_n}{\gcd(M_n,X_n,Y_n)},\qquad
     h_{22}=\frac{\gcd(M_n,X_n,Y_n)}{\gcd(X_n,Y_n,V_n,U_n)},$$
   and $h_{12}$ the unique residue mod $h_{11}$ with
   $$s_Zh_{12}Y_n+h_{22}U_n\equiv0,\qquad s_Zh_{12}X_n+h_{22}V_n\equiv0 \pmod{M_n},$$
   the first congruence alone already having modulus $M_n/\gcd(M_n,Y_n)=h_{11}$.  Equivalently
   $$\boxed{\;h_{12}\equiv-s_Z\,h_{22}\,U_nY_n^{-1}\ \ \bigl(\mathrm{mod}\ h_{11}\bigr).\;}$$
   All four identities verified with **0 failures in 591 instances** ($4\le n\le200$,
   $k\in\{22.4,23.0,23.9\}$).  **[verified]**
2. **$h_{11},h_{22}$ are divisor data; $h_{12}$ is a modular inverse.  Nothing here is
   P-recursive.**  $X,Y$ and $V,U$ *are* holonomic up to elementary factors — new here: the
   Nesterenko $(4,7)$ entries $B_n,C_n$ satisfy a **common order-$3$ recurrence with polynomial
   coefficients of degree $14$** (fitted on $n\le80$, exact on $n=81\ldots120$), which replaces the
   $O(n^{4.1})$ partial-fraction solve and takes the row cache to $n=200$ in $0.2$ s instead of
   $\approx15$ h.  But $h_{12}$ is a modular inverse of a holonomic quantity, and that operation
   destroys holonomy: **[verified]** §1.4, sliding the $2$-adic exponent $E=v_2(T_n)$ by **one bit**
   replaces the $\approx2500$ partial quotients of $h_{12}/h_{11}$ by an independent string (mean
   common prefix $1.2$ partial quotients).
3. **$G$ enters at exactly one place and its price is now exact.**  The integers
   $h_{11},h_{12},h_{22}$, hence the entire continued fraction, do not involve $G$ at all; $G$
   enters only through the single real number $t=(\lambda_N/\lambda_Z)^2$, which selects the
   balance index $i(n)$.  Computing the exact interval $(t_{\rm lo},t_{\rm hi})$ of $t$ on which
   the index is constant gives, for each $(k,n)$, the exact number $d^*(n)$ of decimal digits of
   $G$ that decide $i(n)$:
   $$d^*(n)=6.2674\,n+0.01\pm0.53\ \text{digits}\qquad(=14.430\,n\ \text{nats}),$$
   the same at all three $k$, with **intercept zero**.  This upgrades the $\xi-E_1=14.43$ of
   `P2_STRUCTURE.md` §5.3 from a four-point horizon grid to a per-$n$ exact margin over $591$
   instances.  **[verified + measured]**
4. **No law for the parity, with more power than before.**  Over $4\le n\le200$
   ($591$ instances) the even-frequency is $0.548,0.533,0.503$ (binomial $p=0.20,0.39,1.00$);
   runs tests pass; autocorrelations at lags $1..12$ are noise; $n\bmod m$ for $m\le12$ shows
   nothing surviving multiple-testing; the parities at the three $k$ are mutually independent;
   and the parity is uncorrelated with every $2$-adic invariant tried, including
   $v_2(X_nU_n-V_nY_n)$, $v_2(X_n\xi_2-Y_n)$, $s_2(n)$, $s_2(3n)$.  **No candidate exact law
   matches better than $60\%$.**  **[measured]**
5. **But the parity is *not* memoryless in the $2$-adic direction — and that explains the
   absence of a law across $n$.**  Sliding $E=v_2(T_n)$ by one bit preserves the parity in
   $0.811$ of steps ($1331$ instances), with the agreement decaying to the per-$n$ baseline by lag
   $\approx8$–$16$ bits.
   A fair-coin control that reproduces the exact structure of the slide (§3.3) gives $0.745$
   decaying to $0.50$ by lag $7$: the memory is a generic index-$2$-sublattice-chain effect, not
   Catalan arithmetic.  **One unit of $n$ costs $\lfloor k\rfloor\approx22$ fresh bits of
   modulus, and the memory is $\approx8$ bits — hence no correlation between consecutive $n$.**
   **[measured + control]**
6. **The adelic picture, and the answer to "what is $\rho_2$".**  There is no fixed pair
   $(\rho_\infty,\rho_2)$: the coefficient-ratio targets $-U_n/Y_n\in\mathbb Q_2$ and
   $-(V_nG-U_n)/(X_nG-Y_n)\in\mathbb R$ both move with $n$ (the $2$-adic one has non-convergent
   unit part; the real one grows like $e^{(E_2-E_1)n}$).  What *does* converge, and is the
   right adelic object, is the pair of **targets of the output**:
   $$q_n\to G\ \text{in }\mathbb R,\qquad \frac{p_n}{q_n}\to\xi_2=\zeta_2(2)\ \text{in }\mathbb Q_2 .$$
   $\xi_2=\zeta_2(2)$ is the *common $2$-adic limit of both rows*; this is exactly the content of
   the $24n$ cross-divisibility.  Exact valuation formulas, verified $4\le n\le200$:
   $$v_2\bigl(X_n\xi_2-Y_n\bigr)=v_2(X_n)+24n-1-4s_2(3n)\quad\text{[proved: \texttt{ZUDILIN\_2ADIC.md}]},$$
   $$\boxed{\,v_2\bigl(V_n\xi_2-U_n\bigr)=v_2(V_n)+28n-2s_2(n)-2s_2(3n)-1\,}\quad\text{[verified, new]},$$
   $$\boxed{\,v_2\bigl(X_nU_n-V_nY_n\bigr)=v_2(X_n)+v_2(V_n)+24n-1-4s_2(3n)\,}\quad\text{[verified, new]}.$$
   The last one turns the empirical "$v_2(h_n)\approx24.06\,n$" of `CATALAN_AUDIT.md` §4(d) into an
   **exact closed form**, and reduces the audit's open gap to a single clean $2$-adic statement
   about the Nesterenko row alone.  **[verified $4\le n\le200$]**
7. **The adelic quality is far below every threshold, and the two places do not add up to
   anything.**  With $\delta_\infty=-\log|q_nG-p_n|/\log q_n$ and
   $\delta_2=v_2(q_n\xi_2-p_n)\log2/\log q_n$, measured on the cone minimisers,
   $$v_2(q_n\xi_2-p_n)=(24-k)\,n+O(1)\quad\text{exactly the surplus of the bridge over the modulus},$$
   $$\delta_\infty+\delta_2=0.084,\ 0.071,\ 0.050\ \ (k=22.4,23.0,23.9),$$
   decreasing in $k$: **raising $k$ buys archimedean quality at half the rate at which it spends
   $2$-adic quality**.  Adelic Dirichlet already gives $\delta_\infty+\delta_2=1$ for *any* target
   pair, so the construction is a factor $12$–$20$ below the trivial bound at both places jointly.
   A $2$-adic irrationality proof would need $\delta_2>1$; Ridout / $p$-adic Roth need an
   algebraic target and $\xi_2$ shows no algebraic relation of degree $\le8$ at $600$ bits.
   **The adelic place is not a hidden resource.**  **[measured + proved]**

---

## 1. The objects, explicitly

### 1.1 Setup and notation

Rows exactly as in `P2_STRUCTURE.md` §1 and `POSITIVITY_PROGRAM.md` §3.1:
$$X_n=2^{e_{3n}}D_{6n}^2Q_{3n},\quad Y_n=2^{e_{3n}}D_{6n}^2P_{3n},\quad
  V_n=4^{7n+1}D_{6n}^2B_n,\quad U_n=4^{7n}D_{6n}^2C_n,$$
$e_m=\min(6m,4m+3+\lfloor\log_2(2m-1)\rfloor)$, $S_n=D_{6n}^2$ (the *denominator* modulus, **not**
the $2$-adic one), $T_n=2^{\lfloor kn\rfloor}$ (the $2$-adic modulus), $M_n=S_nT_n$,
$$\mathcal K_n=\{c\in\mathbb Z^2:\ M_n\mid c_ZX_n+c_NV_n,\quad M_n\mid c_ZY_n+c_NU_n\}.$$
Oriented coordinates $(u,v)=(s_Zc_Z,s_Nc_N)$ with $s_Z=\operatorname{sign}(X_nG-Y_n)=(-1)^n$
(Zudilin positivity, a theorem) and $s_N=+1$; weights $\lambda_Z=|X_nG-Y_n|/M_n$,
$\lambda_N=(V_nG-U_n)/M_n$.  Rate constants, now in closed form:
$$E_1=\lim\tfrac1n\log|X_nG-Y_n|=12\log2+12-15\log\varphi=13.0995887908,$$
$$E_2=\lim\tfrac1n\log(V_nG-U_n)=14\log2+12+2\log T^-=14.3931452672,\qquad E_2-E_1=1.2935564764,$$
$\varphi$ the golden ratio, $2\log T^-=-7.3109152606$ Nesterenko's $(4,7)$ rate.  **[verified to
10 digits]**

### 1.2 The Hermite data

> **Proposition 1 [verified, $591/591$ instances, $4\le n\le200$, $k\in\{22.4,23.0,23.9\}$].**
> $$h_{11}=\frac{M_n}{\gcd(M_n,X_n,Y_n)},\qquad
>   h_{22}=\frac{\gcd(M_n,X_n,Y_n)}{\gcd(X_n,Y_n,V_n,U_n)},\qquad
>   [\mathbb Z^2:\mathcal K_n]=h_{11}h_{22}=\frac{M_n}{\gcd(X_n,Y_n,V_n,U_n)},$$
> and $h_{12}$ is the unique residue $0\le h_{12}<h_{11}$ satisfying **both**
> $$s_Zh_{12}Y_n+h_{22}U_n\equiv0\pmod{M_n},\qquad s_Zh_{12}X_n+h_{22}V_n\equiv0\pmod{M_n}.$$
> Moreover $\gcd(M_n,Y_n)=\gcd(M_n,X_n,Y_n)$ in every instance, so the **first** congruence alone
> has modulus exactly $h_{11}$ and pins $h_{12}$:
> $$h_{12}\equiv-s_Z\,h_{22}\,U_n\,Y_n^{-1}\pmod{h_{11}} .$$
> (The $X$-congruence has the smaller modulus $M_n/\gcd(M_n,X_n)$ — since $S_n\mid X_n$ it is
> essentially purely $2$-adic — and is implied.)

*Reading.*  $\mathcal K_n$ is the single congruence $c_ZY_n+c_NU_n\equiv0\pmod{M_n}$ after the
$X$-congruence has been absorbed; $h_{12}/h_{11}$ is literally the residue $-U_nY_n^{-1}$ scaled by
$h_{22}$.  This is the $M_n$-analogue of the statement recorded in `P2_STRUCTURE.md` §6 for the
single-congruence lattice $\mathcal K^S_n$.

**Sizes** ($k=22.4$, $n=200$; `data/h1_hermite200.csv`):
$$\tfrac1n\log h_{11}=27.224,\quad v_2(h_{11})=4460=\lfloor kn\rfloor+2v_2(D_{6n})-v_2(\gcd(M_n,Y_n)),
\quad \tfrac1n\log(\text{odd part of }h_{11})=11.767,$$
$$h_{22}=e^{13.4}\ \text{(sub-exponential, $k$-independent)},\qquad
  m_Z:=\gcd(h_{11},h_{12})=e^{2.382n},\ \ v_2(m_Z)=0 .$$
The last line sharpens `P2_STRUCTURE.md` §2.1: **the Absorption directional mass $m_Z$ is entirely
odd** — it comes from $S_n=D_{6n}^2$, never from the $2$-adic factor $T_n$.  So the $2$-adic bridge
buys size, and the odd modulus buys the forced divisibility; they are disjoint.  **[verified]**

### 1.3 What is holonomic and what is not

| object | status |
|---|---|
| $Q_m,P_m$ (Zudilin) | P-recursive, order $2$, Zudilin's own recurrence (1.3) — **[proved]** |
| $B_n,C_n$ (Nesterenko $(4,7)$) | **P-recursive, order $3$, coefficient degree $14$; the *same* operator annihilates both** — fitted on $n=4..80$, exact on $n=81..120$ (`rowrec.gp`), then used to build $n\le200$ — **[verified]** |
| $D_{6n}^2$, $2^{e_{3n}}$, $4^{7n}$ | elementary, not P-recursive but explicit |
| $\gcd(M_n,X_n,Y_n)$, $\gcd(X_n,Y_n,V_n,U_n)$, hence $h_{11},h_{22}$ | divisor data: $e^{O(\log n)}$ deviations from the elementary factors, no recurrence |
| $h_{12}$ | a **modular inverse** $-s_Zh_{22}U_nY_n^{-1}\bmod h_{11}$: not holonomic, and §1.4 shows it is maximally unstable |
| $i(n)$, $v_1(n)$ | depend on $G$; see §2 |

The order-$3$ operator for $B_n,C_n$ is worth isolating: it makes the whole row cache an $O(n)$
computation, and it is the reason $n\le200$ is available here at all.  Building
$n=4..200$ from the cache at $n\le120$ took **$0.2$ s**, and reproduced the cached rows
$n\le120$ digit-for-digit.  **[verified]**

### 1.4 One bit of modulus destroys the continued fraction

Slide $E:=v_2(T_n)$ by one at fixed $n$ (`run_h7.gp`, `data/h7_slide.csv`).  Then $h_{22}$ is
unchanged, $h_{11}$ doubles exactly ($\Delta\log h_{11}=\log2$ in all $360$ steps), and
$h_{12}$ either stays or gains $h_{11}$ — **one fresh bit**.  Result:

| $n$ | $L(n)$ (partial quotients) | mean common prefix of the CF with $E-1$ |
|---|---|---|
| 20 | 401 | 179.8 |
| 40 | 844 | 1.35 |
| 60 | 1270 | 1.28 |
| 80 | 1679 | 1.15 |
| 100 | 2122 | 1.18 |
| 120 | 2526 | 1.12 |

(the $n=20$ line is the degenerate case where $2^E$ already exceeds the odd part).  Two
independent random continued fractions share $\approx1.35$ initial partial quotients; the
observed value **is** that.  So the arithmetic that builds $\mathcal K_n$ does *not* control the
continued fraction of $h_{12}/h_{11}$ in any stable sense: one bit of $T_n$ resets it.  This is the
structural reason no recursion $v_1(n-1)\to v_1(n)$ can exist (`P2_STRUCTURE.md` §2.4), since
$n\to n+1$ changes $E$ by $\approx22$ bits.  **[verified]**

## 2. Where $G$ enters, exactly

### 2.1 The balance index and its exact stability interval

By Theorem 1 of `P2_STRUCTURE.md` the first minimum is a convergent, so it lies on the *ladder*
$$w_i=\bigl(A_i,\;B_i\bigr),\qquad A_i=|h_{12}q_i-h_{11}p_i|,\qquad B_i=q_ih_{22},$$
$p_i/q_i$ the convergents of $\theta=h_{12}/h_{11}$.  The weighted length is
$\lambda_Z^2A_i^2+\lambda_N^2B_i^2$, so with $t=(\lambda_N/\lambda_Z)^2$
$$i(n)=\arg\min_i\ \bigl(A_i^2+t\,B_i^2\bigr),$$
a purely integer problem in $t$.  The set of $t$ giving a fixed argmin is an open interval
$(t_{\rm lo},t_{\rm hi})$ with **exactly computable rational endpoints**
$t=(A_i^2-A_j^2)/(B_j^2-B_i^2)$.  `run_h2.gp` computes $i(n)$ this way for
$4\le n\le200$, $k\in\{22.4,23.0,23.9\}$; it agrees with the independent weighted Gauss
reduction of `p2_structure/p2run.gp` in **$351/351$** overlapping instances.  **[verified]**

### 2.2 The exact price of $G$

Only $t$ involves $G$.  A perturbation $\delta$ of $G$ moves $\log\lambda_Z$ by
$\delta X_n/|X_nG-Y_n|=\delta e^{(\xi-E_1)n}$, so the index changes as soon as
$\delta>{\rm marg}\cdot e^{-(\xi-E_1)n}$ where ${\rm marg}=\tfrac12\min(\log t-\log t_{\rm lo},
\log t_{\rm hi}-\log t)$.  Writing $d^*(n)$ for the number of decimal digits of $G$ needed:

| $k$ | least-squares fit of $d^*(n)$, $4\le n\le200$ | residual sd |
|---|---|---|
| 22.4 | $6.2668\,n+0.004$ | $0.548$ digits |
| 23.0 | $6.2681\,n+0.002$ | $0.534$ digits |
| 23.9 | $6.2675\,n+0.014$ | $0.506$ digits |

i.e. $d^*(n)=6.267\,n$ **with intercept $0$ and a spread of half a digit** — the balance-margin
term contributes only $O(1)$.  In nats: $14.430\,n$, $14.433\,n$, $14.431\,n$, against
$\xi-E_1=14.43$.  **[verified: exact rational margins; measured: the fit]**

### 2.3 The perturbation experiment, and what changes below the horizon

`run_h3.gp` recomputes *everything* — orientation, $\mathcal K_n$, the full partial-quotient
string, the index — from truncations $G_d$ of $G$ at $d$ decimal digits.  Three regimes,
verified at $n=8,16,24,40,60,80,100,120$:

| regime | $d$ | orientation $s_Z$ | partial quotients | index |
|---|---|---|---|---|
| (i) below the $Z$-horizon | $d<(\xi-E_1)n/\log10=6.267n$ | **wrong** | identical (both signs flip together) | wrong, off by $\sim30$ |
| (ii) between the horizons | $6.267n<d<6.492n$ | right | **different from the 2nd quotient on** | wrong |
| (iii) above both | $d>6.492n$ | right | identical | **correct, for every larger $d$** |

The second horizon $6.492\,n$ digits is $(\tfrac1n\log V_n-E_2)=29.246-14.297=14.949$ nats$/n$, the
Nesterenko analogue.  Sample ($n=120$, $k=23.9$; true index $1366$): at $d=722,\dots,746$ the index
is $1337$ with the orientation wrong; at $d=752,\dots,776$ the CF differs and the index runs
$1336,1341,1352,1357,1362$; from $d=782$ on everything is correct.

**So: no partial quotient of $h_{12}/h_{11}$ depends on any digit of $G$.**  The premise
"the first $O(1)$ partial quotients are set by the archimedean sizes and the later ones by the
$2$-adic residue structure" is false in both halves: the whole string is a $G$-free integer
object, and it is set by the residue $-U_nY_n^{-1}\bmod h_{11}$, which §1.4 shows is scrambled
by a single bit of the modulus.  $G$'s only job is to point at one of the $L(n)\approx21n$ rungs
of the ladder.  **[verified]**

### 2.4 The information budget

Specifying $v_1(n)$ given the rows costs $\log_2L(n)\approx12.1$ bits at $n=200$
($L(200)\approx4272$); specifying its parity costs $1$ bit.  Producing either costs
$14.43\,n$ nats of $G$ — at $n=200$, $1255$ decimal digits for $1$ bit of output.  This is the
sharp form of the obstruction of `P2_STRUCTURE.md` §8.2.  Regime (ii) above shows the parity is
no cheaper than the index: within $30$ digits of the horizon the index sweeps a range of $\approx30$
and its parity flips $\approx15$ times.  **[measured]**

## 3. The parity across $n$

### 3.1 No law ($4\le n\le200$, $591$ instances)

`parity.py`, on `data/h2_balance200.csv`:

| $k$ | even $i(n)$ | binomial $p$ | runs $z$ | longest run of odd $n$ |
|---|---|---|---|---|
| 22.4 | $108/197=0.548$ | $0.20$ | $-1.09$ | 6 |
| 23.0 | $105/197=0.533$ | $0.39$ | $-0.30$ | 6 |
| 23.9 | $99/197=0.503$ | $1.00$ | $+0.79$ | 6 |

* **Autocorrelation** of the parity sequence, lags $1..12$, three $k$: $36$ statistics, largest
  $|z|=1.81$.  Nothing.
* **$n\bmod m$**, $m=2..12$: $33$ contingency tests, three nominal hits
  ($k=22.4$: $m=3$, $p=0.011$; $m=9$, $p=0.005$; $k=23.0$: $m=11$, $p=0.019$), none surviving
  Bonferroni ($0.005\times33=0.17$) and none replicated at another $k$.
* **Cross-$k$**: pairwise agreement $0.477,0.528,0.462$ ($|z|\le1.07$) — the three $k$ are
  *independent* parities, i.e. changing $\lfloor kn\rfloor$ by $\approx1.5n$ bits fully decorrelates.
* **$2$-adic and arithmetic covariates**, point-biserial over $39$ tests:
  $s_2(n)$, $s_2(3n)$, $v_2(n)$, $v_2(X_nU_n-V_nY_n)$, $v_2(X_n\xi_2-Y_n)$, $v_2(V_n\xi_2-U_n)$,
  $v_2(X_n)$, $v_2(D_{6n}^2)$, $L(n)$, $\tfrac1n\log h_{11}$, $\tfrac1n\log h_{22}$, $a_{i+1}$,
  the balance margin.  Two nominal hits at $k=23.9$ ($s_2(3n)$, $p=0.002$; $v_2(n)$, $p=0.010$),
  neither replicated at the other two $k$, neither surviving Bonferroni.
* **Exact-law candidates**: parity of $i(n)$ against the parity of $L(n)$, $n$, $\lfloor kn\rfloor$,
  $s_2(3n)$, $v_2(X_nU_n-V_nY_n)$, $L(n)-i(n)$ — **no candidate matches more than $60\%$**.

The sign of the Casoratian-type determinant $X_nU_n-Y_nV_n$ is constant (both rows are positive
forms and $s_Zs_N=(-1)^n$ is a theorem), so it carries no information; it is the *valuation* of
that determinant that matters, and §5 gives it exactly — and it too is uncorrelated with the
parity.  **[measured]**

### 3.2 But the parity has a $2$-adic memory of $\approx8$ bits

`run_h8.gp` slides $E=v_2(T_n)$ over $120$ consecutive bits at $n=20,30,\ldots,120$.  Under
$E\to E+1$: $h_{22}$ fixed, $h_{11}$ doubles, $h_{12}\in\{h_{12},h_{12}+h_{11}\}$ (one bit), the
metric ratio is unchanged and $\lambda_Z,\lambda_N$ both halve.  Measured parity agreement:

| $n$ | even fraction over $E$ | agreement $E\to E+1$ | $v_1$ survives |
|---|---|---|---|
| 20 | 0.917 | 0.925 | 0.100 |
| 30 | 0.157 | 0.892 | 0.142 |
| 40 | 0.388 | 0.900 | 0.208 |
| 50 | 0.636 | 0.808 | 0.225 |
| 60 | 0.512 | 0.775 | 0.250 |
| 70 | 0.537 | 0.717 | 0.292 |
| 80 | 0.430 | 0.775 | 0.325 |
| 90 | 0.504 | 0.700 | 0.300 |
| 100 | 0.521 | 0.842 | 0.400 |
| 110 | 0.463 | 0.792 | 0.333 |
| 120 | 0.471 | 0.792 | 0.350 |
| **pooled** | **0.503** | **0.811** | **0.266** |

and pooled over $n$, the agreement at lag $j$ in $E$:

| lag $j$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 12 | 16 | 32 | 64 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Catalan | .811 | .730 | .676 | .644 | .621 | .621 | .600 | .593 | .616 | .570 | .549 | .493 |

(the plateau near $0.57$ is the per-$n$ scatter of the even fraction, $\mathrm{sd}=0.19$ over $11$
values of $n$ with $\approx24$ effectively independent $E$ each — it is not a residual correlation.)


The parent's first minimum survives into the index-$2$ sublattice in $0.27$ of steps, and even
conditioned on it *not* surviving the parity is preserved in $0.74$ of steps.  So the parity is a
*shape* statistic that index-$2$ descent tends to preserve, not a fresh coin.

### 3.3 The control: this memory is generic

`control_slide.py` runs the exact same model with the fresh bit replaced by a fair coin — i.e.
$\theta_{j+1}\in\{\theta_j/2,\ \theta_j/2+\tfrac12\}$, the reversed binary expansion of a random
real, with the balance target doubling each step.  $400$ chains $\times\,120$ steps, $\theta$ of
bit-length $400$:

| lag | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 12 |
|---|---|---|---|---|---|---|---|---|---|
| control agreement | .745 | .657 | .575 | .550 | .517 | .514 | .506 | .508 | .501 |

even fraction $0.4991$.  **The Catalan slide is the control plus per-$n$ noise.**  The memory is
index-$2$-sublattice geometry, not arithmetic; and since one unit of $n$ costs $\lfloor k\rfloor\approx22$
fresh bits of modulus while the memory is $\approx8$ bits, the parities at consecutive $n$ are
necessarily uncorrelated — which is what §3.1 measures.  **This is the quantitative reason there
is no law, and it is a negative result about the whole approach: the parity carries no arithmetic
memory that survives one step in $n$.**  **[measured + control]**

## 4. $v_1(n)$, explicitly

> **Proposition 2 [proved, given Theorem 1 of `P2_STRUCTURE.md`].**  In oriented coordinates
> $$v_1(n)=\bigl(\pm A_{i(n)},\ q_{i(n)}h_{22}\bigr),\qquad
>   A_i=|h_{12}q_i-h_{11}p_i|,$$
> with $p_i/q_i$ the convergents of $h_{12}/h_{11}$ and the sign $(-1)^{i(n)}$; in the original
> coordinates $c_Z=s_ZA_{i(n)}\cdot(-1)^{i(n)}$, $c_N=q_{i(n)}h_{22}$.  So $v_1$ needs **no
> separate description**: it is the pair (continued fraction, index), and the continued fraction is
> the modular inverse of §1.2.

Consequences, all measured on $4\le n\le200$:

* $i(n)/L(n)=0.5208\pm0.0137$, $0.5199\pm0.0131$, $0.5200\pm0.0139$ at the three $k$;
  $L(n)/n=21.04,21.29,21.89$; $i(n)/n=10.99,11.10,11.42$.  The midpoint law of
  `P2_STRUCTURE.md` §2.3 holds unchanged to $n=200$.
* $h_{22}=e^{O(\log n)}$ and $k$-independent, so $c_N(v_1)$ is $q_{i(n)}$ up to a sub-exponential
  factor; $\tfrac1n\log|c_Z|\approx14.3$–$14.8$.
* **Erratic, and quantifiably so.**  Given the rows, $v_1(n)$ carries $\log_2L(n)\approx12$ bits.
  Given $v_1(n-1)$ it carries all $12$: the continued fraction at step $n$ is independent of the
  one at step $n-1$ (§1.4 with $\Delta E\approx22\gg$ the $\approx1.3$-quotient overlap of two
  independent CFs), and `P2_STRUCTURE.md` §2.4 already showed that the only surviving arithmetic
  correlation between consecutive $n$ is the forced content $m_Z$.  There is no constant, no
  linear-in-$n$, no holonomic description of $c_Z(n),c_N(n)$.  **[measured]**

## 5. The adelic reformulation

### 5.1 The $2$-adic period is $\zeta_2(2)$, and it is what the bridge *is*

Let $\xi_2:=\lim_{m}P_m/Q_m\in\mathbb Q_2$.  By `ZUDILIN_2ADIC.md` (**proved**) $\xi_2=\zeta_2(2)$
and $v_2(\xi_2-P_m/Q_m)=8m-1-4s_2(m)$ exactly; reverified here at $m=40,80,\dots,240$
($6/6$), and $v_2(\xi_2)=-1$.  New here, by exact computation at $2$-adic precision $7011$ bits:

> **Proposition 3 [verified, $197/197$, $4\le n\le200$].**
> $$v_2\Bigl(\xi_2-\frac{C_n}{4B_n}\Bigr)=28n-2s_2(n)-2s_2(3n)-1 ,$$
> i.e. $v_2(V_n\xi_2-U_n)=v_2(V_n)+28n-2s_2(n)-2s_2(3n)-1$.  Consequently, since the Nesterenko
> rate $28$ exceeds the Zudilin rate $24$,
> $$v_2\bigl(X_nU_n-V_nY_n\bigr)=v_2(X_n)+v_2(V_n)+24n-1-4s_2(3n),$$
> and with $h_n=(X_nU_n-V_nY_n)/S_n$,
> $$v_2(h_n)=v_2(X_n)+v_2(V_n)-2v_2(D_{6n})+24n-1-4s_2(3n)=24n+O(\log n).$$

**Both rows converge $2$-adically to the same number $\zeta_2(2)$; the $24n$ cross-divisibility on
which the whole supercritical construction rests is precisely the *Zudilin* convergence rate, and
the Nesterenko row is $2$-adically better by $4n$.**  This is the "common $2$-adic Catalan period"
lemma of the $5{:}8$ paper, made exact.

*What it does to the audit.*  `CATALAN_AUDIT.md` §4(d) lists as an independent gap that
$v_2(h_n)/n$ is "measurably decreasing ($30.5$ at $n=2\to24.32$ at $n=98$) and nothing proves it
stays above the threshold $22.35$".  Proposition 3 explains the decrease exactly — it is
$v_2(X_n)+v_2(V_n)-2v_2(D_{6n})=O(\log n)$ divided by $n$, plus $-4s_2(3n)/n$ — and shows the limit
is exactly $24$.  The gap now reduces to the single statement
$$v_2\bigl(\zeta_2(2)-C_n/(4B_n)\bigr)\ \ge\ 24n+O(\log n),$$
a $2$-adic convergence statement about the Nesterenko row alone, with the Zudilin half already a
theorem.  **[verified $n\le200$; the reduction is proved, the remaining half is open]**

### 5.2 The output is a simultaneous $(\mathbb R,\mathbb Q_2)$ approximation

With $q_n=(c_ZX_n+c_NV_n)/M_n$, $p_n=(c_ZY_n+c_NU_n)/M_n$,
$$q_n\xi_2-p_n=\frac{c_Z(X_n\xi_2-Y_n)+c_N(V_n\xi_2-U_n)}{M_n},$$
so the *same* coefficient vector produces an archimedean form against $G$ and a $2$-adic form
against $\zeta_2(2)$.  Measured on the cone minimisers, $4\le n\le200$ (`data/h4_adelic.csv`):

| $k$ | $\tfrac1n\log q_n$ | $\tfrac1n\log|q_nG-p_n|$ | $v_2(q_n\xi_2-p_n)/n$ | $\delta_\infty$ | $\delta_2$ | $\delta_\infty+\delta_2$ |
|---|---|---|---|---|---|---|
| 22.4 | 14.827 | $-0.122$ | $1.622$ | 0.0083 | 0.0756 | **0.0839** |
| 23.0 | 14.617 | $-0.331$ | $1.024$ | 0.0227 | 0.0481 | **0.0708** |
| 23.9 | 14.308 | $-0.641$ | $0.117$ | 0.0448 | 0.0056 | **0.0505** |

with $v_2(q_n\xi_2-p_n)/n=1.6219,\,1.0236,\,0.1170$, i.e. $24.02-k$ at all three $k$: **the $2$-adic quality of
the output is exactly the surplus of the bridge ($24$) over the modulus ($k$) that the construction
did not spend.**  So the exact statement equivalent to the whole construction is:

> **The two-row device is a converter, at an exchange rate of $\tfrac12$.**  Measured over
> $100\le n\le200$,
> $$\frac{\mathrm d}{\mathrm dk}\Bigl(-\tfrac1n\log|q_nG-p_n|\Bigr)=0.3461
>   \quad\text{against}\quad \tfrac{\log2}{2}=0.34657,$$
> $$\frac{\mathrm d}{\mathrm dk}\Bigl(\tfrac1n v_2(q_n\xi_2-p_n)\log2\Bigr)=-0.6954
>   \quad\text{against}\quad -\log2=-0.69315 .$$
> **The archimedean gain is exactly half the $2$-adic loss** (the $\tfrac{\log2}{2}$ is
> $|F'(k)|$, the $\log2$ is one bit of $T_n$ per unit $k$).  Hence the total adelic quality
> decreases monotonically through the whole supercritical window and is maximised at the smallest
> admissible $k$.  **[measured, $n\le200$]**

### 5.3 Why this is not a Diophantine resource

* **There is no fixed target pair.**  The task's $(\rho_\infty,\rho_2)$ do not exist: the
  archimedean coefficient target $-(V_nG-U_n)/(X_nG-Y_n)$ grows like $e^{(E_2-E_1)n}$
  ($E_2-E_1=1.2935564764=2\log2+15\log\varphi+2\log T^-$, closed form), and the $2$-adic
  coefficient target $-U_n/Y_n$ has $v_2\in[-11,-5]$ but a unit part that does **not** converge
  (verified: its first $40$ $2$-adic digits jump at every $n$; and $-U_n/Y_n\equiv-V_n/X_n$ to
  $24n$ bits — the two targets agree, which *is* the bridge, but neither converges).
  Ridout, Mahler and $p$-adic Roth are theorems about a **fixed** algebraic target; there is no
  fixed target here.  **[verified]**
* **The one fixed target, $\zeta_2(2)$, shows no algebraicity.**  `algdep` at $2$-adic precision
  $600$ bits returns, for degrees $1..8$, minimal polynomials of $\log_{10}$-height
  $90,60,45,36,31,27,24,22$ — exactly the LLL noise floor $180.6/(d+1)$ for a generic $2$-adic
  number.  So no algebraic relation of degree $\le8$ below the detection threshold; Ridout /
  $p$-adic Roth are inapplicable even in principle.  ($\zeta_2(2)\notin\mathbb Q$ is Calegari 2005;
  its transcendence is open.)  **[measured]**
* **The qualities are an order of magnitude short.**  A $2$-adic irrationality proof for
  $\zeta_2(2)$ from this device would need $\delta_2>1$ (if $\xi=a/b$ then a nonzero $q\xi-p$ has
  $|q\xi-p|_2\ge1/(2b\max(|p|,|q|))$); it delivers $0.076$.  An archimedean irrationality proof
  for $G$ needs $\delta_\infty>0$ *plus* non-vanishing, and non-vanishing is exactly P2$'$.
  Jointly, adelic Dirichlet/Minkowski already gives $\delta_\infty+\delta_2=1$ for **any** pair
  $(\alpha_\infty,\alpha_2)\in\mathbb R\times\mathbb Q_2$, and the device gives $0.05$–$0.084$.
  **[proved + measured]**
* **The exact trade-off identity.**  For any $c\in\mathbb Z^2$, writing
  $\ell_\infty(c)=c_Z(X_nG-Y_n)+c_N(V_nG-U_n)$ and $\ell_2(c)=c_Z(X_n\xi_2-Y_n)+c_N(V_n\xi_2-U_n)$,
  $$\ell_\infty(c)\,(V_n\xi_2-U_n)-\ell_2(c)\,(V_nG-U_n)=c_Z\,(X_nU_n-V_nY_n)\,(\xi_2-G),$$
  an identity in $\mathbb Q[G,\xi_2]$.  In particular the rationality-kernel vector
  $c^{\rm ker}=(aV_n-bU_n,-(aX_n-bY_n))$ of a hypothetical $G=a/b$, which annihilates
  $\ell_\infty$ identically, has
  $$\ell_2(c^{\rm ker})=(X_nU_n-V_nY_n)\,(a-b\xi_2)\ \ne0\quad\text{because }\zeta_2(2)\notin\mathbb Q,$$
  of valuation $v_2(X_nU_n-V_nY_n)+v_2(a-b\xi_2)$ — the **largest possible**.  The vector that is
  archimedean-worst is $2$-adically best, exactly.  This is the same wall as
  `CATALAN_OBSTRUCTION.md` (2.3)/(E2): the archimedean period $G$ and the $2$-adic period
  $\zeta_2(2)$ are different numbers, and Calegari's theorem, which one might hope to import,
  is precisely what makes the kernel vector $2$-adically invisible.  **[proved]**

### 5.4 The adelic statement equivalent to P2$'$

> **P2$'$ (adelic form).**  For infinitely many $n$ there is $(c_Z,c_N)\in\mathbb Z^2$ with
> $(-1)^nc_Z\ge0$, $c_N\ge0$, satisfying
> $$\bigl|\,c_Z\tfrac{X_n}{M_n}+c_N\tfrac{V_n}{M_n}\,\bigr|_2\le1,\qquad
>   \bigl|\,c_Z\tfrac{Y_n}{M_n}+c_N\tfrac{U_n}{M_n}\,\bigr|_2\le1\qquad(\text{$2$-adic and odd parts}),$$
> $$\bigl|c_Z(X_nG-Y_n)+c_N(V_nG-U_n)\bigr|_\infty\le e^{o(n)}M_n\lambda_1(\mathcal K_n).$$
> Equivalently (§5.2): $p_n/q_n$ must approximate $G$ archimedeanly *and* $\zeta_2(2)$
> $2$-adically at the rates of the table, **with both coefficients of one sign**.  The
> congruence side is automatic — it is Proposition 3 — and the sign side is the parity of §3.
> No place-theoretic device changes the sign side: that is Absorption, and §5.3 shows the
> $2$-adic place is already saturated by the construction itself.

## 6. What is new here, and what is corrected

**New.**
1. The closed forms for $h_{11},h_{22}$ and the single-congruence characterisation of $h_{12}$
   (Proposition 1), $591/591$.
2. The order-$3$ P-recursion for the Nesterenko $(4,7)$ entries $B_n,C_n$; row cache to $n=200$
   in $0.2$ s.
3. $d^*(n)=6.267\,n$ exactly, with intercept $0$ and spread $0.5$ digit — the price of $G$ per
   unit $n$, from exact rational stability intervals rather than a horizon grid.
4. The two-horizon structure of the perturbation experiment ($6.267n$ for the orientation,
   $6.492n$ for the partial quotients).
5. Proposition 3: the exact $2$-adic tail of the Nesterenko row, the exact valuation of the
   Catalan cross determinant, and hence an exact closed form for $v_2(h_n)$ — reducing
   `CATALAN_AUDIT.md` gap (d) to one clean $2$-adic statement.
6. The identification of the $2$-adic target of the output as $\zeta_2(2)$, with
   $v_2(q_n\xi_2-p_n)=(24-k)n+O(1)$, and the converter law of §5.2.
7. The $\approx8$-bit $2$-adic memory of the parity, and its fair-coin control (§3.2–3.3).

**Corrections / clarifications to the record.**
* The task framing "$S_n$ = the $2$-adic modulus" is wrong: $S_n=D_{6n}^2$ is the *denominator*
  modulus (mostly odd), $T_n=2^{\lfloor kn\rfloor}$ is the $2$-adic one.  The forced directional
  mass $m_Z$ has $v_2(m_Z)=0$ in all $591$ instances — it is an $S_n$-effect, not a $T_n$-effect.
* `P2_STRUCTURE.md` §2.4's "no recursion" is now explained rather than observed: one bit of
  modulus resets the continued fraction (§1.4), and $n\to n+1$ costs $22$ bits.
* The longest run of consecutive $n$ with $v_1$ outside the cone is $6$ at all three $k$ over
  $n\le200$ (it was $5,4,4$ over $n\le120$).  In every window of seven consecutive $n\le200$
  there is an $n$ with the parity even.

**Still open.**  P2$'$; a proof that $v_2(\zeta_2(2)-C_n/(4B_n))\ge24n+O(\log n)$ (Proposition 3's
Nesterenko half); a proof that $\log\gcd(X_n,Y_n,V_n,U_n)=O(\log n)$; the density of even $i(n)$.

## 7. Scripts and data

`lattice/p2_holonomic/` — PARI/GP unless stated.  Run pattern (concatenate, do not `read()`):

    cat lattice/positivity/rows_pos.gp lattice/p2_structure/p2core.gp \
        lattice/p2_holonomic/hcore.gp lattice/p2_holonomic/run_h2.gp > run.gp && gp -q run.gp

| file | prepend | what |
|---|---|---|
| `hcore.gp` | `rows_pos.gp`, `p2core.gp` | `hdat` (oriented Hermite data), `hpred` (the closed forms), `hchk12`, `ladder`, `balance` (exact stability interval) |
| `rowrec.gp` | — | `fitrec2`, `extend`: the order-3 recurrence for $B_n,C_n$ |
| `adelic.gp` | + `hcore.gp` | `xi2at` ($\xi_2=\zeta_2(2)$ from the Zudilin row), `vform` |
| `build_rows200.gp` | + `rowrec.gp` | rebuilds and extends the row cache to $n=200$ |
| `run_h1.gp` / `run_h1b.gp` | + `hcore.gp` | Proposition 1, $n\le120$ / $n\le200$ |
| `run_h2.gp` / `run_h2b.gp` | + `hcore.gp` | balance index, stability interval, $d^*(n)$ |
| `run_h3.gp` | + `hcore.gp` | the truncation experiment of §2.3 |
| `run_h4.gp` | + `adelic.gp` | the adelic table of §5.2 |
| `run_h5.gp` | + `adelic.gp` | target (non-)convergence, `algdep` on $\xi_2$ |
| `run_h7.gp`, `run_h8.gp` | + `hcore.gp` | the $E$-slide of §1.4, §3.2 |
| `parity.py` | — | all statistics of §3.1 and §4 |
| `control_slide.py` | — | the fair-coin control of §3.3 |

| data file | rows | content |
|---|---|---|
| `data/rows_n200.txt` | 197 | exact $n,X_n,Y_n,V_n,U_n$ for $4\le n\le200$ (identical to `p2_structure/data/rows_all.txt` on $n\le120$) |
| `data/h1_hermite200.csv` | 591 | Proposition 1 checks and the sizes of §1.2 |
| `data/h2_balance200.csv` | 591 | $i(n)$, $(t_{\rm lo},t,t_{\rm hi})$, margins, $d^*(n)$ |
| `data/h3_perturb.csv` | 104 | §2.3 |
| `data/h4_adelic.csv` | 591 | §5.2 (its $v_2(V\xi_2-U)$ column saturates for $n\gtrsim175$; use `h6_val.csv`) |
| `data/h5_targets.txt` | 197 | §5.3 |
| `data/h6_val.csv` | 197 | the exact $2$-adic valuations of Proposition 3 |
| `data/h7_slide.csv`, `data/h8_slide.csv` | 366 / 1331 | the $E$-slide |
| `data/nest_recurrence.txt` | 4 | the order-3 operator annihilating $B_n$ and $C_n$ |

**Cost.**  Everything above runs in under four minutes total on one core, the row build included
(the $\approx1.8$ h partial-fraction solve of `P2_STRUCTURE.md` is no longer needed).  Precision:
$G$ at `\p 3000`; $\xi_2$ at $7011$ $2$-adic bits (from $m=880$ of the Zudilin row).
