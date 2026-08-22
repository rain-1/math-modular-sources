# A $p$-adic irrationality census

*Claude (Opus 5), 2026-08-22.  Scripts and logs: `lattice/padic_irrationality/`.
Companion to `SLOPE_CENSUS.md`, `ROW_LEDGER.md`, `EULER_CRITERION.md`,
`ZUDILIN_2ADIC.md`, `ZETA5_TWO_ROW.md`, `AESZ207.md`, `ACF_ONE_SURFACE.md`.
Source paper: F. Calegari, "Irrationality of certain $p$-adic periods for small $p$",
`papers/calegari_0408214.txt` (= math/0408214).*

---

## 0. Verdict

Score every census cell $(\text{row},p)$ that has a proved or verified $p$-adic Apéry
limit $\xi_p$ by Calegari's criterion, in the form

$$\boxed{\ S_p\;=\;\sigma_p\log p\;-\;k\;-\;\kappa_p\log p\;-\;\log\lambda_1\;
=\;\bigl(v_p(c)+\kappa_p\bigr)\log p-k-\log\lambda_1\ }$$

($\xi_p\notin\mathbb Q$ follows if $S_p>0$; §1 derives this and proves
$S_p>0\iff\theta_p>1$ in Calegari's normalisation).  Then:

* **The criterion is calibrated exactly.**  Calegari's own rows were rebuilt from
  their modular definitions (`calegari_rows.gp`) — the printed $a_n,b_n$ of his
  §3 and §4 are reproduced digit for digit — and the three inputs
  $\sigma_p,\log\lambda_1,k$ were *measured* from those coefficients.  They
  reproduce his $\theta=1.1618804316$ ($\zeta_2(3)$), $1.0469892839$ ($\zeta_3(3)$),
  $1.1618804316$ ($L_2(2,\chi_{-4})=\zeta_2(2)$), $0.9081638111$ ($\zeta_2(5)$)
  and $0.8917942081$ ($\zeta_5(3)$) to all printed digits.
* **No cell in the census has $S_p>0$ except Calegari's own three.**  The census
  maximum is $S_3(\mathbf B)=-0.3521$ ($\theta_3=0.9035$).  **The list of candidate
  new theorems is empty.**
* **The best near miss is `AZ` $\eta$ at $p=5$**, $\xi_5=\tfrac12\zeta_5(3)$,
  $S_5=-0.5858$, $\theta_5=0.8917942081$ — *numerically identical* to Calegari's
  own failed $X_0(5)$ cell, because $\eta$'s level-lowering twin **is** a
  $\Gamma_0(5)$ row (`THEOREM_F_HYPOTHESES.md` §3.1).  Its target
  $\zeta_5(3)$ is, as far as this repo's reading of the literature goes, still
  open: Beukers 2008 reaches only $\zeta_5(3)\pm L_5(3,\chi_5)$.
* **Neither lever helps.**  (i) The cusp-move orbit of `ACF_ONE_SURFACE.md`
  Theorem 1 leaves $v_p(c)\log p-\log\lambda_1$ *invariant* on every placement
  that keeps the slope prime (`scores.gp`, placement scan), so the "smallest
  $\lambda_1$ among realisations" is already the one in the table.  (ii) The
  measured coefficient slope $\sigma_p$ **equals** the overconvergence-radius
  bound in every case, including Calegari's own — there is no slack to harvest.
* **The obstruction is exact and structural** (§7): for an order-2 row
  $S_p\le\frac12 v_p(c)\log p+\kappa_p\log p-k$ with equality iff the row is
  *Fricke-balanced* ($|\lambda_1|=|\lambda_2|$, $|c|=p^{v_p(c)}$).  Calegari's
  rows are balanced with $v_p(c)=\frac{12}{p-1}$ resp. $\frac{24}{p^2-1}$; the
  census's balanced rows have $v_p(c)=3$ or $4$, which is simply too small.

---

## 1. The criterion in our normalisation

### 1.1 The data attached to a cell

A **cell** is a pair (row, prime $p$) together with two sequences of *rationals*
$A_n,B_n$ ($n\ge0$) such that $B_n/A_n\to\xi_p$ in $\mathbb Q_p$.  Attach to it:

| symbol | definition | meaning |
|---|---|---|
| $\kappa_p$ | $\displaystyle\lim_n\frac{-v_p(A_n)}{n}\ \ (\ge0)$ | rate of $p$-power **denominators** of the row |
| $\delta_n$ | the prime-to-$p$ part of a common denominator of $A_n,B_n$ | integrality clearing |
| $k$ | $\displaystyle\lim_n\frac1n\log\delta_n$ (nats per index) | denominator exponent |
| $\lambda_1$ | $\displaystyle\exp\lim_n\frac1n\log\max(|A_n|,|B_n|)$ | archimedean growth |
| $\sigma_p$ | $\displaystyle\lim_n\frac1n v_p\!\Bigl(\xi_p-\frac{B_n}{A_n}\Bigr)$ | $p$-adic quality slope |

For every *modular* census row $A_n\in\mathbb Z$ with $v_p(A_n)=O(\log n)$
(hypothesis **H1** of `CONJ_D_PROOF.md`), so $\kappa_p=0$; and
$\delta_n=d_n^{\,k_{\rm exp}}$ with $d_n=\operatorname{lcm}(1,\dots,n)$, so by the
prime number theorem $k=k_{\rm exp}$ and the two meanings of "$k$" agree.  For the
*hypergeometric* rows they do **not**: Zudilin's $P_m$ needs $D_{2m-1}^2$, i.e.
$k=4$ nats per index but exponent $2$; Nesterenko's needs $D_{6n}^2$, i.e. $k=12$
but exponent $2$.  **$k$ in this document always means the nats-per-index rate.**

### 1.2 What $q_n$ is

Calegari's Lemma 2.1/2.2 needs *integers*.  Put

$$\boxed{\ q_n:=p^{\lceil\kappa_pn\rceil}\,\delta_n\,A_n\in\mathbb Z,\qquad
p_n:=p^{\lceil\kappa_pn\rceil}\,\delta_n\,B_n\in\mathbb Z.\ }$$

For a modular row ($\kappa_p=0$) this is Calegari's own choice with the roles of
$a$ and $b$ interchanged: he has $b_n\in\mathbb Z$ and $[1,\dots,n]^{2k+1}a_n\in\mathbb Z$
and takes $2a_n/b_n=p_n/q_n$; we have $A_n\in\mathbb Z$ and $d_n^{\,k}B_n\in\mathbb Z$
and take $B_n/A_n=p_n/q_n$.  The extra factor $p^{\lceil\kappa_pn\rceil}$ is the
only new ingredient, and it is needed exactly for the non-integral (hypergeometric)
rows.

**The two estimates.**

1. *($p$-adic — the $p$-power clearing is free.)*  $p_n/q_n=B_n/A_n$ as a rational
   number, so
   $$\Bigl|\xi_p-\frac{p_n}{q_n}\Bigr|_p=\Bigl|\xi_p-\frac{B_n}{A_n}\Bigr|_p
   =p^{-\sigma_pn+o(n)} .$$
   Equivalently, since $v_p(q_n)=\lceil\kappa_pn\rceil+v_p(A_n)=o(n)$,
   $$v_p\bigl(q_n\xi_p-p_n\bigr)=\sigma_p\,n+o(n),$$
   which is the sense in which "$\sigma_p=\lim v_p(a_n\xi_p-b_n)/n$" of the brief
   is correct: it is true for the **cleared integer pair**, not for the raw
   $(A_n,B_n)$, for which $v_p(A_n\xi_p-B_n)=(\sigma_p-\kappa_p)n+o(n)$.
   *(Zudilin: $\sigma_2=8$ but $v_2(Q_m\xi_2-P_m)=4m+O(\log m)$ — `zudilin_check.gp`.)*
2. *($\infty$-adic — the $p$-power clearing is the whole cost.)*
   $$\max\bigl(|p_n|,|q_n|\bigr)=\exp\Bigl(\bigl(\log\lambda_1+k+\kappa_p\log p\bigr)n+o(n)\Bigr).$$

### 1.3 The criterion

> **Lemma 1 (Calegari Lemma 2.2 in our normalisation).**  Suppose $\lambda_1>1$
> and $A_n\xi_p-B_n\ne0$ for all large $n$.  Put
> $$\theta_p:=\frac{\sigma_p\log p}{\log\lambda_1+k+\kappa_p\log p},\qquad
> S_p:=\sigma_p\log p-\log\lambda_1-k-\kappa_p\log p .$$
> If $\theta_p>1$ — equivalently $S_p>0$ — then $\xi_p\notin\mathbb Q$.

*Proof.*  With $q_n,p_n$ as above, $q_n\to\infty$ and for every $\varepsilon>0$ and
$n\gg0$
$$0<\Bigl|\xi_p-\frac{p_n}{q_n}\Bigr|_p\le p^{-(\sigma_p-\varepsilon)n}
=\Bigl(e^{(\log\lambda_1+k+\kappa_p\log p)n}\Bigr)^{-\theta_p+\varepsilon'}
\le\bigl(\max(|p_n|,|q_n|)\bigr)^{-\theta_p+\varepsilon''}.$$
Apply Calegari's Lemma 2.2 with $\delta=\theta_p-1-\varepsilon''>0$. $\square$

Since the denominator $\log\lambda_1+k+\kappa_p\log p$ is positive,
$\theta_p>1\iff S_p>0$: **the two formulations are literally equivalent**, and
$S_p$ is $\theta_p$'s numerator minus its denominator.  When $\kappa_p=0$ this is
the brief's $S_p=\sigma_p\log p-k-\log\lambda_1$ verbatim.

*(Nonvanishing.  Calegari proves $a_n-\eta b_n\ne0$ by a modularity argument
(his Lemma 3.2: otherwise $H$ is a polynomial in $f$ and $E_{-2k}$ classical).
For our rows the same statement is available whenever $\xi_p\ne0$ and the
$\mathbb Q$-linear form $B(t)-\xi_pA(t)$ is not a polynomial; it is *not* needed
for any cell below, since no cell has $S_p>0$.  It is recorded here as a
hypothesis that any future positive cell must discharge.)*

### 1.4 The slope law, and the $\kappa$ correction

For an order-2 row in Zagier normalisation
$(n+1)^2u_{n+1}=(an^2+an+b)u_n-cn^2u_{n-1}$ the Casoratian is exact,
$A_nB_{n+1}-A_{n+1}B_n=c^{\,n}/(n+1)^2$, whence
$$v_p\Bigl(\frac{B_{n+1}}{A_{n+1}}-\frac{B_n}{A_n}\Bigr)
=n\,v_p(c)-v_p(A_n)-v_p(A_{n+1})+O(\log n)
=\bigl(v_p(c)+2\kappa_p\bigr)n+O(\log n),$$
and summing the (ultrametrically dominant) tail,

$$\boxed{\ \sigma_p=v_p(c)+2\kappa_p\ }$$

— the law of `ZETA3_TWO_LATTICE.md` §14.2.  Substituting into Lemma 1,

$$\boxed{\ S_p=\bigl(v_p(c)+\kappa_p\bigr)\log p-k-\log\lambda_1\ }$$

so a row with $p$-power denominators gets $2\kappa_p$ of slope but pays $\kappa_p$
back in the archimedean clearing: **the net resource is $\kappa_p\log p$, not
$2\kappa_p\log p$.**  For order $>2$ the two-term Casoratian argument fails and
$\sigma_p$ must be measured (or read off the $p$-adic Newton polygon of the
recurrence); in particular $\sigma_p\ne v_p(\prod_i\lambda_i)$ in general — see
AESZ 207, where $\sigma_2=12$ while $v_2(\prod\lambda_i)=48$.

**Verification.**  `census_scores.gp` measures $\sigma_p$ as
$\frac{1}{200}\bigl[v_p(\xi_{500}-\xi_{400})-v_p(\xi_{500}-\xi_{200})\bigr]$ with
$\xi_N=B_N/A_N$ (exact rational arithmetic, $N=500$ for order 2, $420$ for order 3),
and $\kappa_p$ as $-v_p(A_{400})/400$:

| row | $p$ | $\sigma_p$ measured at $n=400$ | $v_p(c)$ | $-v_p(A_{400})/400$ | $\log\max(|A|,|B|)/n$ at $n=420$ | $\log\lambda_1$ |
|---|---|---|---|---|---|---|
| $\mathbf A$ | 2 | $3.010$ | 3 | $-0.022$ | $2.0650$ | $2.07944$ |
| $\mathbf B$ | 3 | $2.980$ | 3 | $-0.012$ | $1.6351$ | $1.64792$ |
| $\mathbf C$ | 3 | $1.990$ | 2 | $-0.008$ | $2.1830$ | $2.19722$ |
| $\mathbf E$ | 2 | $5.000$ | 5 | $-0.024$ | $2.0661$ | $2.07944$ |
| $\mathbf F$ | 2 | $3.005$ | 3 | $-0.022$ | $2.1858$ | $2.19722$ |
| $\mathbf F$ | 3 | $1.990$ | 2 | $-0.008$ | $2.1858$ | $2.19722$ |
| $\delta$ | 3 | $3.960$ | 4 | $-0.019$ | $2.1737$ | $2.19722$ |
| $\zeta$ | 3 | $2.950$ | 3 | $-0.024$ | $2.9401$ | $2.96488$ |
| $\alpha$ Domb | 2 | $6.000$ | 6 | $-0.024$ | $2.7486$ | $2.77259$ |
| $\eta$ | 5 | $3.000$ | 3 | $-0.012$ | $2.3932$ | $2.41416$ |
| $\varepsilon$=T | 2 | $4.000$ | 4 | $-0.029$ | $3.1242$ | $3.14904$ |
| Cooper $s_{18}$ | 3 | $0.965$ | 1 | $-0.030$ | $2.7504$ | $2.77259$ |

Every entry confirms $\sigma_p=v_p(c)$, $\kappa_p=0$ ($v_p(A_n)>0$ but $O(\log n)$
— the small negative numbers are $O(\log n)/n$), and the archimedean rate
converging to $\log\lambda_1$ from below at rate $O(\log n/n)$.  The two
non-integral rows are checked separately in `zudilin_check.gp`, which reproduces
the two exact laws of `ZUDILIN_2ADIC.md` with zero exceptions at
$m=64,100,128,200,300,400$:
$$v_2(Q_m)=-4m+2s_2(m),\qquad v_2\bigl(\zeta_2(2)-P_m/Q_m\bigr)=8m-1-4s_2(m),$$
i.e. $\kappa_2=4$, $\sigma_2=8=v_2(c)+2\kappa_2$ with $v_2(c)=0$; and confirms
$D_{2m-1}^2P_m\in2^{-e_m}\mathbb Z$ ($k=4$) and $\frac1m\log|Q_m|\to5\log\varphi$.

---

## 2. Calibration against Calegari 2005

`calegari_rows.gp` rebuilds Calegari's rows from their definitions, with no input
from his printed tables:

* **§3, $\zeta_p(3)$ on $X_0(p)$** ($p=2,3,5,7,13$, all genus zero):
  $f=\bigl(\Delta(p\tau)/\Delta(\tau)\bigr)^{1/(p-1)}$,
  $E_2^*=E_2-pV_pE_2$, $E'_{-2}=\theta^{-3}(E_4-V_pE_4)$,
  $H=E_2^*(E'_{-2}+\eta)=\sum(a_n+\eta b_n)f^n$ after $\times24$.
* **§4, the $2$-adic Catalan row on $X_1(4)$**:
  $z=(\Delta(4\tau)/\Delta(\tau))^{1/3}$, $F_1$, $F'_{-1}$ as in his §4.

**Reproduction of his printed sequences** (exact rationals, $q$-series to $O(q^{201})$):

$$b_n\ (p=2):\ 1,\ 24,\ -552,\ 19392,\ -810024,\ 37210944,\ -1815620160$$
$$a_n\ (p=2):\ 0,\ 1,\ 1,\ -\tfrac{8072}{27},\ \tfrac{160841}{9},\ -\tfrac{1088512616}{1125},\ \tfrac{175310024408}{3375}$$
$$b_n\ (X_1(4)):\ \mp\bigl(1,\ 4,\ -28,\ 272,\ -3036,\ 36624,\ -464368\bigr)$$
$$a_n\ (X_1(4)):\ 0,\ 1,\ -3,\ \tfrac{116}{9},\ -\tfrac{331}{9},\ -\tfrac{99116}{225},\ \tfrac{3133076}{225}$$

All agree with the paper (the Catalan $b_n$ up to his global sign; his sixth
$\zeta_2(3)$ entry "$-1088512616/1125$" is printed twice in the arXiv text —
the true $a_6$ is $175310024408/3375$).

**Measured invariants versus his asserted ones.**

| row | $\sigma_p$ measured at $n=50,100,150$ | asserted $\tfrac{12}{p-1}$ | $\log\lambda_1$ measured | asserted $\tfrac{6}{p-1}\log p$ | $k$ measured | asserted $2k{+}1$ |
|---|---|---|---|---|---|---|
| $X_0(2)$ | $11.64,\ 11.79,\ 11.86$ | $12$ | $4.0228,\ 4.0805,\ 4.1025$ | $4.15888$ | $3$ | $3$ |
| $X_0(3)$ | $5.74,\ 5.92,\ 5.89$ | $6$ | $3.1664,\ 3.2196,\ 3.2405$ | $3.29584$ | $3$ | $3$ |
| $X_0(5)$ | $2.82,\ 2.90,\ 2.93$ | $3$ | $2.3218,\ 2.3562,\ 2.3687$ | $2.41416$ | $3$ | $3$ |
| $X_0(7)$ | $1.88,\ 1.95,\ 1.96$ | $2$ | $1.8589,\ 1.8902,\ 1.9036$ | $1.94591$ | $3$ | $3$ |
| $X_0(13)$ | $0.90,\ 0.96,\ 0.97$ | $1$ | $1.2285,\ 1.2177,\ 1.2509$ | $1.28247$ | $3$ | $3$ |
| $X_1(4)$ | $7.76,\ 7.86,\ 7.91$ | $8$ | $2.6449,\ 2.6993,\ 2.7201$ | $2.77259$ | $2$ | $2$ |

**This settles the "overconvergence radius versus coefficient rate" question in
the negative:** Calegari's $p$-adic rate is proved by Buzzard's continuation
theorem to be *at least* $p^{12/(p-1)}$ (resp. $2^8$), and the coefficients say it
is *exactly* that.  There is no gap between the two, in his rows or in ours; the
overconvergence radius is pinned by the position of the second cusp, where
$\Theta$ genuinely has a pole.

**Reproduction of his $\theta$'s** (`scores.gp`):

| Calegari statement | $\sigma_p$ | $k$ | $\log\lambda_1$ | $S_p$ | $\theta_p$ | his printed $\theta$ |
|---|---|---|---|---|---|---|
| Thm 3.3, $\zeta_2(3)$ | 12 | 3 | $6\log2$ | $\mathbf{+1.15888}$ | $1.161880$ | $1.1618804316$ |
| Thm 3.4, $\zeta_3(3)$ | 6 | 3 | $3\log3$ | $\mathbf{+0.29584}$ | $1.046989$ | $1.0469892839$ |
| Thm 4.2, $\zeta_2(2)=L_2(2,\chi_{-4})$ | 8 | 2 | $4\log2$ | $\mathbf{+0.77259}$ | $1.161880$ | $1.1618804316$ |
| §3.1 remark, $\zeta_2(5)$ | 12 | 5 | $6\log2$ | $-0.84112$ | $0.908164$ | $0.9081638111$ |
| §3.2 remark, $\zeta_5(3)$ | 3 | 3 | $\tfrac32\log5$ | $-0.58584$ | $0.891794$ | $0.8917942081$ |
| ($\zeta_7(3)$, $\zeta_{13}(3)$, $\zeta_2(7)$) | 2,1,12 | 3,3,7 | — | $-1.05409,\ -1.71753,\ -2.84112$ | $0.786876,\ 0.598941,\ 0.745394$ | — |

Every printed digit is reproduced.  **$S_p>0$ is exactly Calegari's $\theta>1$.**

---

## 3. Normalisation of $\xi_p$: which $p$-adic $L$-value is which

Two conventions collide and must be kept apart (`DWORK_CRYSTALS_PRIMER.md` §4.3,
`SOURCES_S18_ZUDILIN.md` §6.1, both confirmed against Beukers' own introduction):

* $\zeta_p(s):=L_p(s,\mathbf 1)$, the Kubota–Leopoldt $p$-adic zeta function.
  Its even values are **nonzero**.  At $p=2$ the Teichmüller character is
  $\omega=\chi_{-4}$, so $\zeta_2(2)=L_2(2,\mathbf 1)$ **is** the number Calegari
  writes $L_2(2,\chi)$ with $\chi$ of conductor $4$; at $p=3$, $\omega=\chi_{-3}$
  and $\zeta_3(2)=L_3(2,\mathbf 1)$.
* The Beilinson–Vologodsky/CDvS "$\zeta_p(m)=L_p(m,\omega^{1-m})$" vanishes for
  even $m$.  Different object.

Theorem F (`EULER_CRITERION.md` §0) gives every census $\xi_p$ in closed form:
$\xi_p=-Q(w+1)\kappa_p^{\rm KL}$, $\kappa_p^{\rm KL}=\tfrac12L_p(w+1,\psi\omega^{-w})$
if $\varphi=\mathbf 1$ and $0$ otherwise.  The identifications used below are those
verified there to $\ge p^{2991}$ digits.

---

## 4. The table

$\log\lambda_1$ is the archimedean growth rate of $\max(|A_n|,|B_n|)$;
$k$ is in nats per index; $\kappa_p$ is the $p$-power denominator rate.
Cells with $\xi_p=0$ are marked **deg** (degenerate: $\xi_p$ is rational, so no
irrationality statement is available at any $S_p$; listed for completeness).

| # | cell | $p$ | $\xi_p$ | $\sigma_p$ | $\kappa_p$ | $k$ | $\log\lambda_1$ | $S_p$ | $\theta_p$ |
|---|---|---|---|---|---|---|---|---|---|
| C1 | Calegari $X_0(2)$, $\zeta_p(3)$ | 2 | $\tfrac12\zeta_2(3)$ | 12 | 0 | 3 | $4.15888$ | $\mathbf{+1.15888}$ | $1.161880$ |
| C2 | Calegari $X_1(4)$, Catalan | 2 | $\tfrac12\zeta_2(2)$ | 8 | 0 | 2 | $2.77259$ | $\mathbf{+0.77259}$ | $1.161880$ |
| C3 | Calegari $X_0(3)$, $\zeta_p(3)$ | 3 | $\tfrac12\zeta_3(3)$ | 6 | 0 | 3 | $3.29584$ | $\mathbf{+0.29584}$ | $1.046989$ |
| 1 | Zagier $\mathbf B$ $(9,3,27)$ | 3 | $\tfrac12\zeta_3(2)$ | 3 | 0 | 2 | $1.64792$ | $-0.35208$ | $0.903484$ |
| 2 | AZ $\eta$ $(11,5,125)$ | 5 | $\tfrac12\zeta_5(3)$ | 3 | 0 | 3 | $2.41416$ | $-0.58584$ | $0.891794$ |
| — | Calegari $X_0(5)$ (fails) | 5 | $\tfrac12\zeta_5(3)$ | 3 | 0 | 3 | $2.41416$ | $-0.58584$ | $0.891794$ |
| 3 | Zagier $\mathbf E$ $(12,4,32)$ | 2 | $\tfrac12\zeta_2(2)$ | 5 | 0 | 2 | $2.07944$ | $-0.61371$ | $0.849561$ |
| 4 | AZ $\delta$ $(7,3,81)$ | 3 | $\tfrac14\zeta_3(3)$ | 4 | 0 | 3 | $2.19722$ | $-0.80278$ | $0.845538$ |
| — | Calegari $X_0(2)$, $\zeta_p(5)$ (fails) | 2 | $\tfrac12\zeta_2(5)$ | 12 | 0 | 5 | $4.15888$ | $-0.84112$ | $0.908164$ |
| 5 | $\alpha$ Domb $(10,4,64)$ | 2 | $\tfrac13\zeta_2(3)$ | 6 | 0 | 3 | $2.77259$ | $-1.61371$ | $0.720454$ |
| 6 | Zagier $\mathbf C$ $(10,3,9)$ | 3 | $\tfrac12\zeta_3(2)$ | 2 | 0 | 2 | $2.19722$ | $-2.00000$ | $0.523495$ |
| 7 | Zagier $\mathbf F$ $(17,6,72)$ | 3 | $\tfrac58\zeta_3(2)$ | 2 | 0 | 2 | $2.19722$ | $-2.00000$ | $0.523495$ |
| 8 | Zagier $\mathbf A$ $(7,2,-8)$ **deg** | 2 | $0$ | 3 | 0 | 2 | $2.07944$ | $-2.00000$ | $0.509737$ |
| 9 | Zagier $\mathbf F$ $(17,6,72)$ | 2 | $\tfrac12L_2(2,\chi_{12})$ | 3 | 0 | 2 | $2.19722$ | $-2.11778$ | $0.495433$ |
| 10 | AZ $\zeta$ $(9,3,-27)$ **deg** | 3 | $0$ | 3 | 0 | 3 | $2.96488$ | $-2.66904$ | $0.552541$ |
| 11 | $\varepsilon=$T $(12,4,16)$ | 2 | $\tfrac14\zeta_2(3)$ | 4 | 0 | 3 | $3.14904$ | $-3.37645$ | $0.450898$ |
| 12 | cusp $L(f,2)$, level 12 **deg** | 2 | $0$ | 2 | 0 | 2 | $2.77259$ | $-3.38629$ | $0.290470$ |
| 13 | Zudilin Catalan row | 2 | $\zeta_2(2)$ | 8 | **4** | 4 | $2.40606$ | $-3.63347$ | $0.604139$ |
| 14 | Cooper $s_{18}$ | 3 | $\tfrac12\zeta_3(2)$ | 1 | 0 | 2 | $2.77259$ | $-3.67398$ | $0.230192$ |
| 15 | $\zeta(5)$ level 16, row II | 2 | $\tfrac7{32}\zeta_2(5)$ | 1 | 0 | 5 | $1.38629$ | $-5.69315$ | $0.108537$ |
| 16 | $\zeta(5)$ level 16, row I | 2 | $\tfrac7{32}\zeta_2(5)$ | 1 | 0 | 5 | $2.03560$ | $-6.34245$ | $0.098520$ |
| 17 | AESZ 207 (rank 4) | 2 | *unidentified* | 12 | 0 | 4 | $11.40234$ | $-7.08458$ | $0.540032$ |
| 18 | Nesterenko $(4,7)$ row | 2 | $\zeta_2(2)$ | 28 | **14** | 12 | $7.65071$ | $-9.94665$ | $0.661157$ |

**Cells with no $p$-adic Apéry limit, hence not scored** (recorded so the census
is closed): Zagier $\mathbf D$ and $\gamma$ Apéry ($c=\pm1$, $\sigma_p=0$ at every
prime); Cooper $s_7,s_{10}$ ($\sigma_p=0$ at $p=2,3,5,7$, `cooper_scores.gp`
confirms slopes $\le0.03$ at $n=400$); the level-12 $\zeta(5)$ rows
($\prod\lambda_i=-1$, all $\sigma_p=0$, *proved*); the Brown–Zudilin cellular
$\zeta(5)$ row ($\sigma_p=0$ everywhere, extension splits); the level-24
$\zeta(7)$ row (measured $\sigma_p=0$ at $p\le7$ to $n=218$ — note that
`ROW_LEDGER.md` correction **C3**'s "$\sigma_2=5$ is forced by
$\prod\lambda_i=-32$" is an over-reading: $\sigma_p$ is $v_p$ of the *relevant
root*, and here $x^2-8x+8$ has both roots of valuation $3/2$ and $x^2+4x-4$ both
of valuation $1$, so $\sigma_2\le\frac32$ at best, and is measured to be $0$);
the Domb cuspidal apparatus and rows $\mathbf B,\mathbf C$ at $p=2$,
$\gamma,\delta,\eta$ at $p=2$, $\alpha$ at $p=3$ (`EULER_CRITERION.md` §4.3).

---

## 5. Cells with $S_p>0$

$$\textbf{Three, and all three are Calegari's own theorems.}$$

* **C1** $S_2=6\log2-3=+1.158883$: this *is* Calegari Theorem 3.3, $\zeta_2(3)\notin\mathbb Q$.
* **C2** $S_2=4\log2-2=+0.772589$: this *is* Calegari Theorem 4.2,
  $\zeta_2(2)=L_2(2,\chi_{-4})\notin\mathbb Q$ (also Beukers 2008 Cor. 7.3).
* **C3** $S_3=3\log3-3=+0.295837$: this *is* Calegari Theorem 3.4, $\zeta_3(3)\notin\mathbb Q$.

**The list of candidate new theorems is empty.**  No row constructed in this
project — no Zagier sporadic, no Almkvist–Zudilin row, no Cooper row, no Domb or
$\varepsilon$ row, no level-16 $\zeta(5)$ row, no AESZ rank-4 row, and neither of
the two hypergeometric Catalan rows (Zudilin, Nesterenko) — attains $S_p>0$ at
any prime.  Consequently there is nothing here for a literature audit to check
for novelty against Calegari 2005 or Beukers 2008: the only positive cells are
literally those two papers' theorems, recovered.

---

## 6. Near misses, ranked

"Near" is taken as $S_p>-1$ (the cut $S_p>-0.5$ of the brief isolates only #1).
For each, (i) is the target already known irrational, and (ii) can the
archimedean companion be improved?

| rank | cell | $S_p$ | $\theta_p$ | $\xi_p$ | (i) already proved irrational? |
|---|---|---|---|---|---|
| 1 | Zagier $\mathbf B$, $p=3$ | $-0.35208$ | $0.903484$ | $\tfrac12\zeta_3(2)$ | **Yes** — Beukers 2008 Cor. 7.3 ($\zeta_3(2)$). Not Calegari. |
| 2 | AZ $\eta$, $p=5$ | $-0.58584$ | $0.891794$ | $\tfrac12\zeta_5(3)$ | **Apparently not** — Calegari fails at $p=5$; Beukers 2008 Thm. 11.2/Cor. 11.4 gives only $\zeta_5(3)\pm L_5(3,\chi_5)$. **Open target.** |
| 3 | Zagier $\mathbf E$, $p=2$ | $-0.61371$ | $0.849561$ | $\tfrac12\zeta_2(2)$ | Yes — Calegari Thm. 4.2 (= C2, the *same number*, with $\theta=1.1619$). |
| 4 | AZ $\delta$, $p=3$ | $-0.80278$ | $0.845538$ | $\tfrac14\zeta_3(3)$ | Yes — Calegari Thm. 3.4 (= C3, the *same number*, with $\theta=1.0470$). |

(Then a gap: #5 Domb at $-1.61371$, $\tfrac13\zeta_2(3)$ = C1's number.)

### 6.1 (i) Is it the same number Calegari already did?

Sorting the whole table by $\xi_p$ up to $\mathbb Q^\times$ — which is all that
matters, since $\xi\notin\mathbb Q\iff r\xi\notin\mathbb Q$ for $r\in\mathbb Q^\times$ —
gives six families, and in **every** family Calegari's own row is the best
realisation:

| $\xi_p$ up to $\mathbb Q^\times$ | realisations, best $S_p$ first | status of the number |
|---|---|---|
| $\zeta_2(3)$ | Calegari $X_0(2)$ $+1.1589$; Domb $-1.6137$; $\varepsilon$ $-3.3765$ | proved (Cal. Thm 3.3; Beukers Thm 11.2) |
| $\zeta_2(2)$ | Calegari $X_1(4)$ $+0.7726$; $\mathbf E$ $-0.6137$; Zudilin $-3.6335$; Nesterenko $-9.9467$ | proved (Cal. Thm 4.2; Beukers Cor 7.3) |
| $\zeta_3(3)$ | Calegari $X_0(3)$ $+0.2958$; $\delta$ $-0.8028$ | proved (Cal. Thm 3.4; Beukers Thm 11.2) |
| $\zeta_3(2)$ | $\mathbf B$ $-0.3521$; $\mathbf C=\mathbf F$ $-2.0000$; $s_{18}$ $-3.6740$ | proved (Beukers Cor 7.3), **not** by Calegari |
| $\zeta_5(3)$ | $\eta=$ Calegari $X_0(5)$ $-0.5858$ | **open** (Beukers has only $\zeta_5(3)\pm L_5(3,\chi_5)$) |
| $\zeta_2(5)$ | Calegari $X_0(2)$ $-0.8411$; level-16 II $-5.6931$; level-16 I $-6.3425$ | **open** |
| $L_2(2,\chi_{12})$ | $\mathbf F$ $-2.1178$ | **open** as far as this repo can tell (Beukers explicitly cannot reach $L_3(2,\chi_{12})$ by his condition (A); $\chi_{12}$ at $p=2$ is not in his list either) |
| AESZ 207's $\xi_2$ | AESZ 207 $-7.0846$ | *the number itself is unidentified* (`AESZ207.md` §4.2 rules out 188 Dirichlet targets, the $\zeta_p(m)$, Iwasawa logs, and the $\widehat\Gamma$-span) |

Two entries deserve a note.

* **$\eta$ *is* Calegari's $X_0(5)$ cell.**  Its $(\sigma_5,k,\log\lambda_1)=(3,3,\frac32\log5)$
  coincide exactly with Calegari's, so $\theta_5=0.8917942081$ to all digits.
  This is not a coincidence: `THEOREM_F_HYPOTHESES.md` §3.1 shows $\eta$'s
  level-lowering twin ($q\mapsto-q$) is a $\Gamma_0(5)$ row of degree $1$ in $t$,
  i.e. the *same* genus-zero curve Calegari uses.  Row $\eta$ therefore adds no
  information about $\zeta_5(3)$ beyond the paper.
* **$\mathbf B$ is the $p=3$, level-$9$ analogue of Calegari's $X_1(4)$ Catalan
  row.**  Its twin is a $\Gamma_0(9)$ row (`THEOREM_F_HYPOTHESES.md` §3.1,
  $\Gamma_{\mathbf B}=\gamma\Gamma_0(9)\gamma^{-1}$).  On $X_0(p^2)$ with
  $f=(\Delta(p^2\tau)/\Delta(\tau))^{1/(p^2-1)}$ the Fricke involution sends
  $p^{24/(p^2-1)}f$ to $1/f$, so the Fricke exponent is $w_F=24/(p^2-1)$ and the
  balanced score is
  $$S_p=\tfrac12 w_F\log p-2=\frac{12\log p}{p^2-1}-2:\qquad
  p=2:\ +0.77259,\quad p=3:\ -0.35208,\quad p=5:\ -1.19528 .$$
  The $p=2$ value is C2 and the $p=3$ value is $\mathbf B$ **exactly**.  So $\mathbf B$
  is the natural continuation of Calegari's §4 to $p=3$, and it fails for the same
  reason his §3 fails at $p=5$: the Fricke scale shrinks like $p^{1/p^2}$ while
  $k$ stays fixed.

### 6.2 (ii) Can the archimedean companion be improved?

Two levers were checked.

**Lever A — cusp moves / other placements** (`ACF_ONE_SURFACE.md` Theorem 1).
Every second-order Apéry system admits two cusp moves
$$(a,b,c)\longmapsto(\mu-2\lambda,\ b-\lambda,\ \lambda^2-\lambda\mu),\qquad
\text{roots }-\lambda,\ \mu-\lambda ,$$
one for each characteristic root; the three placements of a system put,
respectively, the three non-MUM fibres at $t=\infty$.  Since the gauge factor is
a modular *unit*, all placements carry the same $\xi_p$ up to a rational scalar
(when they stay in the same Eisenstein plane), so the smallest $\lambda_1$ among
them wins.  `scores.gp` runs the scan; the result is:

| system | placements $(c,\lambda_1)$ | $S_p$ per placement |
|---|---|---|
| $\{\mathbf A,\mathbf C,\mathbf F'\}$ | $(-8,8)$, $(9,9)$, $(72,9)$ | $p=2$: $-2.0000$, —, $-2.1178$;  $p=3$: —, $-2.0000$, $-2.0000$ |
| $\mathbf E$ | $(32,8)$, $(-16,4)$, $(32,8)$ | $p=2$: $-0.6137$, $-0.6137$, $-0.6137$ |
| $\mathbf D$ | $(-1,11.090)$, $(1,11.180)$, $(124,11.180)$ | no $p$-power $c$ except $p=2$: $-3.0279$ |
| $\mathbf B$ | complex roots: the move is **not defined over $\mathbb Q$** (`ACF_ONE_SURFACE.md` Cor. 1.1) | only $(27,\sqrt{27})$: $-0.3521$ |

So the cusp move **never improves the score**: on the $\mathbf E$ orbit
$v_2(c)\log2-\log\lambda_1$ is exactly invariant ($5\log2-\log8=4\log2-\log4=2\log2$),
and on the $\mathbf A/\mathbf C/\mathbf F$ orbit it is invariant on every placement
that retains the slope prime ($2\log3-\log9=0$ for both $\mathbf C$ and
$\mathbf F'$), while the placement that changes the slope prime also changes the
Eisenstein family and hence the $p$-adic number itself ($\mathbf A$ sits in the
$\mathcal T$-family with $\xi_2^{\mathbf A}=0$, `ACF_ONE_SURFACE.md` Theorem 3).
**The invariance is structural**: writing the two singular points of the
Picard–Fuchs equation as $t=1/\lambda,1/\mu$, the cusp move is the Möbius map
that sends one of them to $\infty$, which rescales *both* the archimedean radius
and the $p$-adic radius by the same modular unit.

The one genuinely different realisation of a shared $\xi_p$ is a different *row*,
not a different placement — and the table of §6.1 shows that the best row in each
family is already listed.

**Lever B — the overconvergence radius.**  Calegari's $p$-adic rate is not read
off coefficients: it is Buzzard's theorem (analytic continuation of finite-slope
overconvergent eigenforms) applied to the ordinary locus of $X_0(p)$/$X_1(4)$.
The question is whether that bound can exceed $\sigma_p=v_p(c)$.  It cannot:

* On **his** rows the two agree to the measurement precision (§2 table:
  $\sigma_2\to12$, $\sigma_3\to6$, $\sigma_5\to3$, $\sigma_7\to2$, $\sigma_{13}\to1$,
  $\sigma_2^{\rm Cat}\to8$), so his use of Buzzard is *sharp*, not lossy.
* On **our** rows $\sigma_p=v_p(c)$ is measured (§1.4 table) and equals the
  distance from $t=0$ to the second cusp: $R^{\mathbf C}_\xi$ converges exactly on
  $|t_{\mathbf C}|_3<3^{\sigma_3}=9$ (`PADIC_PERIOD.md` §3), and $t_{\mathbf C}=1/9$
  is the value at the $I_6$ cusp $\tau=0$ (`ACF_ONE_SURFACE.md` §3: $1-9t_{\mathbf C}$
  has a simple zero there).  Past that point $\Theta$ has a pole, so no
  continuation theorem can help.
* The only place where the naive product bound *over*-states the truth is order
  $>2$: AESZ 207 has $\sigma_2=12=v_2(53248)$ while $v_2(\prod_i\lambda_i)=48$;
  the slope is $v_p$ of the *one* root that carries it (there, the rational
  singular point $z=1/53248$ with exponents $(0,1,3,4)$), not of the product.

**Conclusion for §6:** neither lever moves any cell, and in particular neither
moves #2, the only near miss whose target is open.

---

## 7. Why every cell fails, exactly

For an order-2 row, $c=\lambda_1\lambda_2$, hence
$\log\lambda_1\ge\frac12\log|c|\ge\frac12v_p(c)\log p$, with equality throughout
iff $|\lambda_1|=|\lambda_2|$ **and** $|c|=p^{v_p(c)}$.  Call such a row
*Fricke-balanced at $p$*.  Therefore

$$\boxed{\ S_p\ \le\ \tfrac12\,v_p(c)\log p\;+\;\kappa_p\log p\;-\;k,
\qquad\text{equality iff Fricke-balanced.}\ }$$

Geometrically: the archimedean radius of convergence is the distance to the
nearest fold, the $p$-adic radius is the distance to the second ordinary
component, and balance means the fold sits at the *geometric mean* of the two —
which is precisely Calegari's picture (the elliptic point of $X_0(2)$ at
$f=-2^{-6}$, exactly halfway between $\|f\|=1$ and $\|f\|=2^{12}$).

**What a positive cell needs.**  With $\kappa_p=0$ the condition $S_p>0$ is
$v_p(c)\log p>2k$, i.e.

| $k$ (weight/Hodge depth) | $p=2$ | $p=3$ | $p=5$ | $p=7$ |
|---|---|---|---|---|
| $k=2$ (weight-3 source) | $v_2(c)\ge6$ ($c\ge64$) | $v_3(c)\ge4$ ($\ge81$) | $v_5(c)\ge3$ ($\ge125$) | $v_7(c)\ge3$ |
| $k=3$ (weight-4 source) | $v_2(c)\ge9$ ($\ge512$) | $v_3(c)\ge6$ ($\ge729$) | $v_5(c)\ge4$ | $v_7(c)\ge4$ |
| $k=5$ (weight-6, $\zeta(5)$) | $v_2(c)\ge15$ | $v_3(c)\ge10$ | $v_5(c)\ge7$ | — |

Calegari's rows: $X_0(2)$ has $v_2(c)=12\ge9$ ✓; $X_0(3)$ has $v_3(c)=6\ge6$ ✓
(marginally — this is why $\theta_3=1.047$ only); $X_0(5)$ has $v_5(c)=3<4$ ✗;
$X_1(4)$ has $v_2(c)=8\ge6$ ✓.

The census's balanced rows are $\mathbf B$ ($c=3^3$), $\delta$ ($c=3^4$),
$\eta$ ($c=5^3$) — all short.  The census's large-$v_p(c)$ rows
($\mathbf E$: $2^5$; Domb: $2^6$; $\varepsilon$: $2^4$; $\mathbf F$: $2^3\!\cdot\!3^2$)
are *unbalanced*, and the unbalancing costs exactly
$\log\lambda_1-\frac12\log|c|$ (Domb: $\log16-\log8=\log2$; $\varepsilon$:
$\log 23.31-\log4=1.76$).  **Every failure in the census is one of these two
defects, and usually both.**

**Three concrete targets this analysis points at** (none realised in the
census; recorded as the open problem, not as results):

1. **An order-2, $k=2$, Fricke-balanced row with $c=2^6,2^7$ or $3^4,3^5$ or $5^3$.**
   Scores would be $+0.07944$, $+0.42602$, $+0.19722$, $+0.74653$, $+0.41416$.
   Zagier's classification of *integral* order-2 sporadics in the normalisation
   $(n+1)^2u_{n+1}=(an^2+an+b)u_n-cn^2u_{n-1}$ contains no such row — but note
   that **Calegari's own Catalan row is not in that normalisation**:
   $$(n+1)^2u_{n+1}=(4-32n^2)u_n-256(n-1)^2u_{n-1}$$
   has the shift $(n-1)^2$ in place of $n^2$ and $a$-coefficient $-32$ with
   $n$-coefficient $0$, i.e. $a^2=4c$ (double root $-16$) with $c=2^8$.  A search
   over this *shifted* normalisation — which the census never performed — is the
   obvious next step, and is the only place in this analysis where a genuinely new
   positive cell could live.
2. **A $\kappa_p$ row with small $k$.**  Since the net resource is $\kappa_p\log p$
   and the cost is $k$, the hypergeometric world could in principle win — but both
   instances here have $k$ inflated by the *range* of the lcm ($D_{2m-1}^2$,
   $D_{6n}^2$), so $\kappa_2\log2=2.77$ against $k=4$ (Zudilin) and
   $9.70$ against $k=12$ (Nesterenko).  Wanted: a $p$-adic Padé row whose
   denominators run only to $\operatorname{lcm}(1..n)$.
3. **$\zeta_5(3)$ and $\zeta_2(5)$** are the two named open targets that the
   census actually touches.  For $\zeta_5(3)$ one needs a $p=5$ realisation with
   $v_5(c)\ge4$ balanced (Calegari and $\eta$ both give $3$); for $\zeta_2(5)$ one
   needs $v_2(c)\ge15$ balanced at $k=5$ (Calegari's $X_0(2)$ gives $12$, the
   level-16 rows give $1$).

---

## 8. Honest status of each input

| input | status |
|---|---|
| Lemma 1 (the criterion) | **proved** here, from Calegari Lemma 2.2, modulo the nonvanishing clause, which no cell needs |
| $\sigma_p=v_p(c)+2\kappa_p$ | **proved** for order-2/3 rows from the exact Casoratian (`CRYSTAL_THEOREM_F.md` Cor. 3.2); **measured** to $n=400$ here for all twelve |
| Calegari's $a_n,b_n,\theta$ | **reproduced exactly** from his modular definitions (`calegari_rows.gp`) |
| $\xi_p$ identifications | **verified** to $\ge p^{2991}$ in `EULER_CRITERION.md` §4.1; Theorem F's hypothesis (b) is open for $\mathbf E,\mathbf F,\alpha,\varepsilon$ at $p=2$ |
| $\xi_2^{\rm Zud}=\zeta_2(2)$ | **proved** (`ZUDILIN_2ADIC.md`), depends on Beukers Props 5.1, 6.1 |
| Nesterenko's $v_2$ laws | **hypotheses of a Lean structure**, not theorems (`CATALAN_DIRECTIONAL.md` §162–166) — flagged; the cell is far from positive either way |
| AESZ 207 $\sigma_2=12$, $\lambda_1$, $k=4$ | **measured/verified** (`AESZ207.md`); $\xi_2$ **unidentified** |
| level-16 $\zeta(5)$ $\sigma_2=1$ | **measured**, not derived (`ZETA5_TWO_ROW.md` caveat 4) |
| "$\zeta_5(3)$, $\zeta_2(5)$, $L_2(2,\chi_{12})$ are open" | **read off this repo's summary of Beukers 2008** (`SOURCES_S18_ZUDILIN.md` §6.1) — **must be confirmed by the literature audit**, which this document does not perform |

**Novelty flag.**  §5 is empty, so there is nothing to audit for novelty.  The
statements that *do* need the literature audit are the negative ones in the last
row of the table above: that $\zeta_5(3)$, $\zeta_2(5)$ and $L_2(2,\chi_{12})$ are
not already known irrational.  If any of them is known, near miss #2 loses even
its aspirational value.

---

## Scripts

All under `lattice/padic_irrationality/`, PARI/GP, absolute paths, `quit;` at the end.

* `calegari_rows.gp` — rebuilds Calegari's §3 rows on $X_0(p)$, $p=2,3,5,7,13$, and
  his §4 row on $X_1(4)$, from the modular definitions; prints his $a_n,b_n$;
  measures $\sigma_p$, $\log\lambda_1$, $k$ at $n=50,100,150$ ($q$-series to $O(q^{201})$).
* `census_scores.gp` — exact integer models for the six Zagier and six AZ
  sporadics to $n=500$/$420$; denominator-exponent scan, $\kappa_p$, measured
  $\sigma_p$ at $n=400$, archimedean rate; prints $S_p,\theta_p$ per cell.
* `cooper_scores.gp` — Cooper $s_7,s_{10},s_{18}$ to $n=420$: $k$, integrality of
  $A_n$, archimedean rate, slopes at $p=2,3,5,7$.
* `zudilin_check.gp` — Zudilin's Catalan row to $m=420$: reproduces
  $v_2(Q_m)=-4m+2s_2(m)$ and $v_2(\zeta_2(2)-P_m/Q_m)=8m-1-4s_2(m)$ exactly,
  confirms $k=4$, $\kappa_2=4$, $\log\lambda_1=5\log\varphi$; the working example
  of the $\kappa$-corrected criterion.
* `scores.gp` — the assembly: the table of §4 and §5, sorted by $S_p$; plus the
  cusp-move placement scan of §6.2.

*(GP note for whoever edits these: outside a `{ }` block gp terminates a statement
at end of line, and a function body written as `f(x) = a; b;` swallows `b` into the
closure.  Both bit this session.)*
