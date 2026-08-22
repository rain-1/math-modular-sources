# Explicit moment constructions for the Catalan two-row problem

*Follow-up: `consolidation/POSITIVITY_PROGRAM.md` §4 extends the moment grid to $m\le60$,
answers the open question of §6 (the optimal pair is $(j_0,j_0+1)$ with $j_0/m\approx0.30$),
and gives the asymptote $+3.81$ per unit $m$ for the best single moment.*


**No claim of irrationality is made anywhere below.** All numbers are finite-$n$ measurements;
by `CATALAN_AUDIT.md` §4(a) every covolume-driven column here is reproduced verbatim by a
rational surrogate for $G$, so none of it is evidence about $G$.

---

## Verdict

1. **The "free non-vanishing" idea works, and it is not new: it *is* Mechanism B.** The
   positive-cone vector of `CATALAN_POSITIVITY.md` §3 already *is* an $I_n(P)$ with
   $P=\alpha+\beta w^n$, $\alpha,\beta>0$ — verified: the cone-minimum coefficient vector has
   $(s_Zc_Z,s_Nc_N)>0$ at every $n\le44$. So positivity is already free in the existing
   construction, and the open point is only that $\alpha,\beta$ are *congruence-selected*.
2. **Making $P$ explicit in the naive senses fails, quantitatively and for a structural reason.**
   With a **fixed common denominator** the objective $\log\bigl(\mathrm{den}\cdot I_n(P)\bigr)$ is a
   *positive sum* over $j$ and is therefore minimised at a vertex: a non-negative-coefficient $P$
   is never better than the best single moment, which measures
   $n^{-1}\log(\mathrm{den}_j|M_{n,j}|)\ \ge\ +5.2\,(n{=}2),\ +6.4\,(3),\ +7.1\,(4),\ +7.7\,(5)$ —
   positive and **increasing**. Sums of squares $P=Q^2$ (and $P=\sigma_0+w\sigma_1$, the full
   cone of polynomials $\ge0$ on $(0,\infty)$) do allow cancellation, but they hit a hard floor:
   $\min_{Q\in\mathbb Z^{d+1}}Q^{\mathsf T}H_nQ = e^{-8.25n+o(n)}$ **independently of $\deg Q$**,
   because $n^{-1}\log\det H_n^{(d)}\approx-8.2\,(d+1)$, i.e. Minkowski's
   $\det^{1/(d+1)}$ is flat in $d$. Higher degree buys nothing in size and costs denominator.
   Best found over $P=Q^2,\,wQ^2$ with true (post-cancellation) denominators:
   $n^{-1}\log(\mathrm{den}\cdot I_n(P)) = +6.12\,(n{=}3),\ +6.52\,(4),\ +7.31\,(5)$ — never negative.
3. **The one genuinely new and useful finding is §5: the 2-adic bridge belongs to the whole
   moment family, and the Zudilin/Nesterenko pair $\{0,n\}$ is *not* the best pair in it.**
   The mixed-minor valuation $v_2(h_{j_0,j_1})/n$ lies in $[22,33]$ for *every* pair, clustering at
   $24$–$25$ — the same $24.06n$ law as $\{0,n\}$ (`CATALAN_AUDIT.md` §4(d)). But the
   half-log-covolume predictor $F$ is smaller for other pairs by $\approx2$ per $n$:

   | $n$ | best pair | $F_{\text{best}}$ | $v_2(h)/n$ | $F(0,n)$ | $v_2(h_{0,n})/n$ |
   |---|---|---|---|---|---|
   | 3 | $(2,5)$ | **−0.411** | 24.67 | +1.458 | 23.67 |
   | 4 | $(2,3)$ | **−0.362** | 24.00 | +2.507 | 24.50 |
   | 5 | $(4,5)$ | **−0.082** | 25.20 | +2.245 | 23.40 |

   This is the only actionable output: **if the two-row programme is worth pushing, it should be
   pushed on a better-chosen pair of moments than $\{0,n\}$.** (It changes nothing logically:
   Lemma P2 is unaffected, and §4(a) of the audit applies unchanged.)
4. **Lemma P2 is untouched.** Nothing here attacks it, and nothing here is a lattice/valuation
   substitute for it.

---

## 1. The moment family, exactly

$$K_n(x,y):=K_Z^{(3n)}(x,y)=\frac{x^{m-1/2}(1-x)^m y^m(1-y)^{m-1/2}}{(1-xy)^{m+1}},\quad m=3n,\qquad
w=\frac{xy}{1-xy},$$
$$M(m,j):=\iint_{[0,1]^2}K_Z^{(m)}w^j\,dx\,dy = I\!\left(m-\tfrac12+j,\;m,\;m+j,\;m-\tfrac12,\;m+1+j\right).$$

**Range.** Near $(1,1)$ the integrand is $\asymp u^m v^{m-1/2}(u+v)^{-(m+1+j)}$, so
$$M(m,j)<\infty\iff j<m+\tfrac12\iff 0\le j\le m=3n .$$
The Nesterenko exponent $j=n$ sits comfortably inside; $j=3n$ is the exact ceiling
(confirmed: `nestgen(2,3)` returns a *negative* value, i.e. the analytic continuation, and
`nestgen(2,2)` still agrees with the series).

**Exact evaluation.** `lattice/catalan_explicit/moments.gp` generalises the Nesterenko
partial-fraction solver `nest()` of `lattice/catalan_positivity/rows_common.gp` from the special
case $(m,j)=(3n,n)$ to arbitrary $(m,j)$: replace $N=3n\to m$, $4n\to m+j$, and the constant by
$c=m!\,\Gamma(m+\frac12)/\bigl(\Gamma(\frac32)(m+j)!\bigr)$. It returns $A,B\in\mathbb Q$ with
$M(m,j)=AG-B$, by exact rational linear algebra (no numerics anywhere).

**Endpoint checks (exact, `moments.gp` + `t14`).**
* $j=0$: $M(m,0)=8(-1)^m\bigl(Q_mG-P_m\bigr)$ for $m\le8$, ratio $=8.000000000000000$.
  Zudilin's Theorem 3 (`papers/0201024v3.pdf`, eq. (6)) reads
  $u_nG-v_n=\frac{(-1)^n}{4}\iint x^{n-1/2}(1-x)^ny^n(1-y)^{n-1/2}(1-xy)^{-n-1}$; the repo's
  $(Q,P)$ with $Q_0=1,Q_1=7/4,P_0=0,P_1=13/8$ therefore differ from the paper's $(u,v)$ by a
  factor $2$. Irrelevant — it cancels in every ratio used below. (`integrals.gp` recorded the
  same kind of normalisation slack.)
* $j=n$: $M(3n,n)=\bigl(V_nG-U_n\bigr)/\bigl(4^{7n}D_{6n}^2\bigr)$ for $n\le4$, ratio $=1.000000000000000$.

So Mechanism A's identity is now an exact *interpolation*: the Zudilin row and the Nesterenko
$(4,7)$ row are the $j=0$ and $j=n$ members of one $(3n+1)$-term family of linear forms in $1,G$.

## 2. Denominators and 2-adic structure ( `momtab.gp`, $n\le5$ )

For every tabulated $(n,j)$:
$$\boxed{\,v_2\bigl(\mathrm{den}\,B_{n,j}\bigr)=v_2\bigl(\mathrm{den}\,A_{n,j}\bigr)+1\,}$$
without exception (a clean, cheap, provable-looking lemma candidate; it is the moment-family
form of the extra factor $2$ between $Q_m$ and $P_m$ denominators).

$n^{-1}\log(\mathrm{den}_{n,j}\cdot|M(3n,j)|)$, $\mathrm{den}=\mathrm{lcm}(\mathrm{den}A,\mathrm{den}B)$:

| $n$ | min over $j$ | at $j=$ | value at $j=0$ | value at $j=n$ | $\min_j n^{-1}\log|M|$ |
|---|---|---|---|---|---|
| 1 | +1.207 | 2 | 2.112 | 4.097 | −7.46 |
| 2 | +5.235 | 2 | 6.726 | 5.235 | −7.81 |
| 3 | +6.402 | 2 | 8.320 | 8.166 | −7.79 |
| 4 | +7.105 | 2 | 9.485 | 8.137 | −7.75 |
| 5 | +7.721 | 5 | 8.738 | 7.721 | −7.72 |

All positive and rising. The size term is pinned: $n^{-1}\log|M(3n,j)|\in[-8.3,-4.6]$, with the
asymptotic $-3\cdot2.40606=-7.218$ for fixed $j$ (Zudilin's rate at index $3n$); the true
denominators are $e^{(12.8\ \text{to}\ 23)n}$, far larger than the $e^{-8}$ headroom.

## 3. Positivity is free — and, at fixed denominator, worthless

$P\ge0$ on $(0,\infty)$ $\Rightarrow$ $I_n(P)>0$: true, unconditional, and exactly the content of
`CATALAN_POSITIVITY.md` §2.1 extended from $\{1,w^n\}$ to $\mathbb R[w]$.

**Vertex lemma.** Fix any common denominator $\mathrm{den}$ with
$\mathrm{den}\cdot A_{n,j},\mathrm{den}\cdot B_{n,j}\in\mathbb Z$ for all $j\le J$. For $p\in\mathbb Z^{J+1}_{\ge0}$,
$\mathrm{den}\cdot I_n(P)=\sum_jp_j\,\mathrm{den}\,M_{n,j}$ is a sum of positive terms, hence
$\ \ge\min_j \mathrm{den}\,M_{n,j}$, attained at a vertex $p=e_j$. **Non-negative coefficients give
no cancellation at all.** (Table above: floor $\ge+5.2$ for $n\ge2$.)

**Sums of squares.** $P\ge0$ on $[0,\infty)$ iff $P=\sigma_0+w\sigma_1$ with $\sigma_i$ SOS, so
$I_n(P)=Q^{\mathsf T}H\,Q+Q'^{\mathsf T}H'Q'$ with $H=(M_{n,i+k})$, $H'=(M_{n,i+k+1})$, both
positive definite (Hankel matrices of a positive measure) — block diagonal, so the minimum is the
minimum over the two blocks separately. The objective $\log(\mathrm{den}\cdot|I_n(P)|)$ is
invariant under $P\mapsto cP$, $c\in\mathbb Q^\times$, so integer $Q$ is WLOG.

`sos.gp`: with $\mathrm{lcm}$-denominators,

| $n$ | $n^{-1}\log\det H^{(d)}$ | $n^{-1}\log\min_{Q\in\mathbb Z^{d+1}}Q^{\mathsf T}HQ$ | best objective |
|---|---|---|---|
| 3 | $-7.64,\,-15.93,\,-24.28,\,-32.08,\,-38.45$ | $-7.64,-8.13,-8.22,-8.22,-8.22$ | $+8.53$ ($d{=}1$) |
| 4 | $-7.61,\,-15.79,\,-24.20,\,-32.49,\,-40.31,\dots$ | $-7.61,-8.01,-8.23,-8.26,-8.26$ | $+9.60$ ($d{=}1$) |
| 5 | $-7.57,\,-15.67,\,-24.02,\,-32.43,\,-40.66,\dots$ | $-7.57,-7.91,-8.13,-8.25,-8.25$ | $+9.62$ ($d{=}1$) |

The minimising $Q$ are $(-1,1)$, $(1,-2,1)$, $(1,-4,4,-1)$, i.e. $P=(1-w)^2,(1-w)^4,\dots$ —
exactly the "Chebyshev-like, positive" candidates the idea proposed. **They saturate:**
$\det(H^{(d)})^{1/(d+1)}=e^{-8.2n+o(n)}$ for every $d$, so no degree helps. Structurally: pushing
forward to $u=xy$, the measure $\nu_n$ on $w\in(0,\infty)$ has a Pareto tail
$\nu_n(w>W)\asymp W^{-(3n+1/2)}$ — only $3n$ finite moments and no concentration — so its
orthogonal-polynomial norms $\|\pi_d\|^2$ are all of the same exponential order. Gaussian
quadrature gives no exponential gain here, unlike the Legendre/Beukers situation.

`sossearch.gp` (true post-cancellation denominators, $|q_i|\le5$–$6$, $d\le3$, both shifts) gives
best objectives $+5.23\,(n{=}2)$, $+6.12\,(3)$, $+6.52\,(4)$, $+7.31\,(5)$.

**Important correction to the naive reading.** The vertex lemma is *only* about a fixed common
denominator. Once one allows the *numerators* $A_P=\sum p_jA_{n,j}$ and $B_P$ to acquire a large
common factor, non-negative $p$ **can** produce arbitrarily large cancellation — and that is
precisely the 2-adic bridge. So "positive $\Rightarrow$ no cancellation" is false; the correct
statement is "positive $\Rightarrow$ no cancellation *in the analytic factor*, all the gain must
come from the arithmetic factor", which is `CATALAN_POSITIVITY.md`'s Mechanism B verbatim.

## 4. The Minkowski/cone coefficients: sign pattern and law ( `signs.gp`, $4\le n\le44$ )

* The cone-minimum vector has $(s_Zc_Z,\,s_Nc_N)$ **both positive at every $n$** (by construction);
  the *unrestricted* first minimum $\lambda_1$ lies in the cone (up to $\pm$) at $25$ of $41$
  values of $n$ ($61\%$) — i.e. it is not a law, just the near-orthogonality of §3.2.
* $n^{-1}\log|c_Z|\to\approx14.0$, $n^{-1}\log|c_N|\to\approx12.9$, so
  $\log|c_N/c_Z|\approx-1.15\,n$ (from $-10.4$ at $n=10$ to $-50.8$ at $n=44$).
* The invariant $\log\bigl|c_N\lambda_N/(c_Z\lambda_Z)\bigr|$ is $O(1)$ throughout — range
  $[-4.83,+4.00]$, no drift. **So the only law obeyed by the coefficients is "balance the two
  terms", $c_N/c_Z\approx\lambda_Z/\lambda_N$.** That ratio is a *provable* quantity ($E_1,E_2$,
  Lemma P3), so the *scale* of an explicit $P$ is prescribed; what is not prescribed is the
  congruence class, which is the whole of Lemma P2.
* 2-adic structure of the coefficients: $v_2(c_Z)\in\{0,1,2,3,4\}$ (mostly $0$),
  $v_2(c_N)\in[5,12]$ and **$\ge5$ at every $n\in[4,44]$**, drifting up like $\log n$. Nothing
  matched a recognisable sequence (`bestappr`, OEIS-shaped ratios, binomial transforms all
  negative). The coefficients look like generic congruence solutions.

## 5. The 2-adic bridge across the family, and a better pair ( `bridge.gp`, `pairs.gp` )

For a pair $S=\{j_0,j_1\}$, clear denominators ($d=\mathrm{lcm}$, $a_i=dA_i$, $b_i=dB_i\in\mathbb Z$) and
form the mixed minor $h=a_{j_0}b_{j_1}-a_{j_1}b_{j_0}$ — the exact analogue of the $h_n$ of
`CATALAN_AUDIT.md` §4(d). Predictor (as in `CATALAN_POSITIVITY.md` §3.2, cone-min
$\asymp\sqrt{\mathrm{covol}}$):
$$F(j_0,j_1)=\frac1{2n}\Bigl[\log|M_{j_0}|+\log|M_{j_1}|+2\log d-v_2(h)\log2\Bigr].$$

Findings over all $\binom{3n+1}{2}$ pairs, $n=2,\dots,5$:

* $v_2(h)/n\in[22,33]$ for **every** pair, monotone increasing in $j_0+j_1$; the pairs that
  minimise $F$ have $v_2(h)/n\approx24$–$25$. The $24.06n$ bridge is a property of the whole
  Beukers family at level $3n$, not of the Zudilin$\times$Nesterenko pairing.
* The $\{0,n\}$ pairing is **not** optimal. Best pairs and their $F$ (table in the Verdict):
  $(2,5),(2,3),(4,5)$ with $F=-0.411,-0.362,-0.082$ against $F(0,n)=+1.458,+2.507,+2.245$.
  The gain is $\approx2$ per $n$ and comes from the denominator term $2\log d$, not from $v_2(h)$:
  interior $j$ have markedly smaller $\mathrm{den}_{n,j}$ than $j=0$.
* Without any extra 2-adic division (`pairs.gp`, canonical lattice
  $L_S=\{p:\sum p_jA_j\in\mathbb Z,\ \sum p_jB_j\in\mathbb Z\}$, cone minimum over $p\ge0$) the best
  pairs give $n^{-1}\log(\text{cone-min}) = -0.074\,(n{=}2,\ (0,1))$, $+0.108\,(3,\ (2,5))$,
  $+0.157\,(4,\ (2,3))$, $+0.325\,(5,\ (3,4))$, versus $+0.777,+1.163,+1.457,+1.526$ for $\{0,n\}$.
  Same ordering, same $\approx1$–$1.2$ per $n$ advantage.

**Caveat, stated plainly.** $F$ is a covolume heuristic, and $F$ is *rising* with $n$ in the table
(as $\kappa_n\uparrow\sigma$ does in `CATALAN_POSITIVITY.md` §3.2), so the negative values at
$n=3,4$ are pre-asymptotic. Extending the family beyond $n=5$ needs a faster `nestgen`
(the solver is $O((6n)^3)$ in big rationals and already dominates the runtime at $m=15$).

## 6. Lemma candidates

* **Lemma E1 (cheap, likely provable).** $v_2(\mathrm{den}B_{n,j})=v_2(\mathrm{den}A_{n,j})+1$ for all
  $0\le j\le3n$. Verified for all $(n,j)$, $n\le5$.
* **Lemma E2 (convergence).** $\iint K_Z^{(m)}w^j<\infty\iff j\le m$. Proved above by the
  homogeneity count at $(1,1)$; only worth recording because it fixes $J_{\max}=3n$ exactly.
* **Lemma E3 (negative, proved).** For fixed common denominator, $\min\{\mathrm{den}\cdot I_n(P):
  P\in\mathbb Z_{\ge0}[w],P\ne0\}=\min_j\mathrm{den}\cdot M_{n,j}$ (vertex of a positive sum). Hence
  non-negative-coefficient constructions cannot beat the best single moment.
* **Lemma E4 (negative, measured).** $\det(H_n^{(d)})^{1/(d+1)}=e^{-8.2n+o(n)}$ uniformly in
  $d\le3n/2$, so the SOS route has an $n$-rate floor of $-8.25$ against denominators
  $\ge e^{12.8n}$. A proof would follow from the Pareto tail exponent $3n+\frac12$ of $\nu_n$.
* **Open, and the useful one.** Does $\min_{j_0<j_1}F(j_0,j_1)$ stay below $F(0,n)$ by $\Theta(1)$
  as $n\to\infty$, and does the minimising pair follow a law ($j_0\approx m/4$, $j_1=j_0+1$)?
  If yes, the whole `CATALAN_DIRECTIONAL`/`CATALAN_POSITIVITY` programme should be re-run on that
  pair; Lemma P1's threshold $k_*=22.35129$ would be tested against a $v_2(h)/n$ that is
  *measurably larger* than $24.06$ for the interior pairs.

**None of this weakens or strengthens Lemma P2**, and none of it is evidence about $G$
(the `control.gp` argument of `CATALAN_AUDIT.md` §4(a) applies verbatim to every table here).

Scripts: `lattice/catalan_explicit/{moments.gp, momtab.gp, sos.gp, sossearch.gp, signs.gp,
control.gp, pairs.gp, bridge.gp}`; each run as
`cat lattice/catalan_positivity/rows_common.gp lattice/catalan_explicit/moments.gp X.gp > run.gp; gp -q run.gp`
(`signs.gp` needs only `rows_common.gp`). Longest single run $\approx110$ s (`pairs.gp`).
