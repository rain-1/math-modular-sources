# Positivity and direction in the Catalan two-row lattice

**Verdict.**

* **Mechanism A works as an identity and fails as a criterion.** The two source integrands are
  literally the *same* Beukers kernel: with
  $K_Z^{(m)}(x,y)=\frac{x^{m-1/2}(1-x)^m y^m(1-y)^{m-1/2}}{(1-xy)^{m+1}}$ one has the exact
  identity (verified to 57 digits at 200 random points, $n\le4$, `integrals.gp`)
  $$\boxed{\;K_N^{(n)}(x,y)=K_Z^{(3n)}(x,y)\cdot w(x,y)^{\,n},\qquad w=\frac{xy}{1-xy}\;}$$
  so the combined form is a **single** integral
  $c_Z\Lambda_Z+c_N\Lambda_N=\int_0^1\!\!\int_0^1 K_Z^{(3n)}(x,y)\bigl[\alpha+\beta\,w(x,y)^n\bigr]dx\,dy$
  with $\alpha=\kappa_Z(-1)^{3n}c_Z$, $\beta=\kappa_N c_N$, $\kappa_\bullet>0$ explicit.
  Since $w$ maps the open square **onto $(0,\infty)$**, the bracket is of one sign iff
  $\alpha,\beta$ are of one sign: for mixed signs the integrand changes sign on the hyperbola
  $xy=\text{const}$ and *both* pieces have positive measure. So the integral representation
  delivers exactly the trivial cone statement (10) and **nothing more**; for mixed-sign
  coefficients, deciding the sign of the integral is literally deciding whether
  $\alpha/\beta\lessgtr-\Lambda_N/\Lambda_Z$, i.e. the question itself. No Beukers-style free
  non-vanishing is available.
* **Mechanism B is the real content, and it is *not* refuted.** Modelled exactly on the honest
  lattice $\mathcal K_n=\{c:\ q_n,p_n\in\mathbb Z\}$ for $n\le44$: the shortest vector of
  $\mathcal K_n$ **already lies in the positive cone** at essentially every $n$
  (cone-min $/\lambda_1\in[1.0,1.7]$ at 7 of 10 tabulated $n$, $\le5.0$ always), and
  $\lambda_2/\lambda_1\le2.4$ throughout — the Catalan lattice is *near-orthogonal* in the
  scaled metric, so Sol's degenerate model does **not** occur here. The headline number is
  $$n^{-1}\log|q_nG-p_n|\ \text{ for the best positive-cone vector}\ =\ -0.25\ \text{at }n=44,
  \ \text{tracking } \tfrac1{2n}\log\operatorname{covol}\ \text{to within }0.03 .$$
  It is **negative at every $n$** and its predicted limit is $F(k)<0$. But the *pigeonhole /
  Dilworth* version of the idea gives only the second-minimum scale (§3.3) — the observed
  $\sqrt{\operatorname{covol}}$ behaviour is empirical, not a theorem, and §4 shows no theorem
  of that shape can be proved from covolume alone.
* **Mechanism C is dead.** For every modulus tested ($3,5,7,16,64$) and every realised class,
  **100 %** of the classes contain vectors of *both* signs (`character.gp`). This is not an
  empirical accident: sign is a half-plane condition and a congruence class is a coset of a
  finite-index sublattice, which meets both open half-planes. No character-valued invariant of
  $(c_Z,c_N)$ can predict $\operatorname{sign}(c_Z\Lambda_Z+c_N\Lambda_N)$.

**No claim of irrationality is made.** The net effect is a reduction of the supercritical
obstruction to the two explicit lemmas of §5.

---

## 1. Setup (exact, reusing `lattice/catalan_audit`)

Rows as in `CATALAN_DIRECTIONAL.md` §1: $X_n=2^{e_{3n}}D_{6n}^2Q_{3n}$, $Y_n=2^{e_{3n}}D_{6n}^2P_{3n}$,
$V_n=4^{7n+1}D_{6n}^2B_n$, $U_n=4^{7n}D_{6n}^2C_n$; $S_n=D_{6n}^2$, $T_n=2^{\lfloor kn\rfloor}$,
$M_n=S_nT_n$, $\sigma(k)=12+k\log2$.

I deliberately **do not** use the surrogate congruence pair
($c\cdot(Y,U)\equiv0\,(S)$, $c\cdot(a^Z,a^N)\equiv0\,(T)$): it does *not* by itself make
$p_n=(c_ZY+c_NU)/M_n$ integral (one only gets $T_n\mid a^Z_n(c_ZY+c_NU)$, and
$v_2(a^Z_n)=11\text{–}17$ for $n\le24$, `diag.gp`). Instead
$$\mathcal K_n=\{c\in\mathbb Z^2:\ M_n\mid c_ZX_n+c_NV_n\ \text{ and }\ M_n\mid c_ZY_n+c_NU_n\},$$
computed by `matkerint`/`mathnf`, so $q_n,p_n\in\mathbb Z$ by construction (checked at every
tabulated $n$). Its index satisfies $n^{-1}\log[\mathbb Z^2:\mathcal K_n]=\kappa_n$ with
$\kappa_n\uparrow$: $25.79\,(n{=}8)$, $26.57\,(20)$, $26.90\,(28)$, $27.09\,(44)$, against
$\sigma(22.4)=27.5265$ — deficit $0.94\to0.44$, consistent with $\kappa_n=\sigma+o(1)$ but
**not proved**.

**Signs.** $\operatorname{sign}(X_nG-Y_n)=(-1)^{3n}$ and $V_nG-U_n>0$, verified numerically for
$n\le5$ and guaranteed by (2.14)–(2.15) and Zudilin (8). Write $s_Z=(-1)^{3n}$, $s_N=+1$,
$\lambda_Z=|X_nG-Y_n|/M_n$, $\lambda_N=(V_nG-U_n)/M_n$. The **positive cone** is
$\mathcal P=\{c:\ s_Zc_Z\ge0,\ c_N\ge0\}\setminus\{0\}$, on which
$q_nG-p_n=s_Zc_Z\lambda_Z+c_N\lambda_N>0$ **unconditionally**.

## 2. Mechanism A in detail

Zudilin (his eq. (8), index $m$) and Nesterenko ((2.14), index $n$) are both
$I(a,b,c,d,e)=\iint x^a(1-x)^by^c(1-y)^d(1-xy)^{-e}$:
$$Z_m=I(m-\tfrac12,m,m,m-\tfrac12,m+1),\qquad
J_n=I(4n-\tfrac12,3n,4n,3n-\tfrac12,4n+1).$$
At $m=3n$ the parameter difference is $(n,0,n,0,n)$, i.e. exactly the factor $w^n$ above.
Numerically `i2`-style checks give $4^{7n}D_{6n}^2J_n=V_nG-U_n$ **exactly** (ratio
$1.000000000000$), and the Zudilin row matches $4(-1)^m2^{e_m}D_{6n}^2Z_m$ up to the constant
$2^5$ coming from the $Q_0,Q_1$ normalisation — the *sign* $(-1)^m$ is exact.

Consequences.
1. Positivity of the combination is available **iff** $(s_Zc_Z,c_N)\ge0$. (Mechanism B's cone.)
2. Push forward along $u=xy$: the form becomes $\int_0^1(\alpha+\beta(u/(1-u))^n)\,d\mu_n(u)$ for
   a positive measure $\mu_n$. $\{1,w^n\}$ is a Chebyshev system on $(0,1)$, so a mixed-sign
   bracket has **exactly one** sign change — the most favourable possible situation for a sign
   argument, and it still yields nothing: one sign change places no constraint on the sign of
   the integral. The two moments $\mu_n(1)=\Lambda_Z$, $\mu_n(w^n)=\Lambda_N$ are precisely the
   two numbers whose ratio is unknown.
3. A common-kernel Padé/Rivoal reading changes nothing: any third row in the same family adds
   $w^{n'}$, so all rows are moments of one measure, and the sign of a signed combination of
   moments is again equivalent to knowing the moment ratios.

## 3. Mechanism B, quantified

### 3.1 Method
Work in the scaled coordinates $(u,v)=(s_Zc_Z\lambda_Z,\ c_N\lambda_N)$; $\mathcal K_n$ becomes a
planar lattice $L_n$ of covolume $\operatorname{covol}=[\mathbb Z^2:\mathcal K_n]\lambda_Z\lambda_N$,
and the task is to minimise the linear functional $u+v$ over $L_n\cap(\text{closed first quadrant})$.
LLL-reduce in that metric, then for each $i\in[-600,600]$ solve the two cone inequalities for an
interval in $j$ and evaluate at its endpoints (the functional is linear, so the minimum over the
segment is at an endpoint). Exact integers throughout, $G$ to 3000 digits;
$|q_nG-p_n|$ recomputed from the returned $(q_n,p_n)$ as a consistency check.

### 3.2 Data (`cone_sweep.gp`, $k=22.4$, 8 s total)

| $n$ | $\kappa_n$ | $\frac{1}{2n}\log\mathrm{covol}$ | $\frac1n\log\lambda_1$ | $\frac1n\log\lambda_2$ | $\frac1n\log(\text{cone-min})$ | ratio to $\lambda_1$ | $\frac1n\log q_n$ |
|---|---|---|---|---|---|---|---|
| 8 | 25.7946 | −0.9284 | −0.9256 | −0.9181 | **−0.9099** | 1.13 | 13.955 |
| 12 | 26.0792 | −0.7408 | −0.7773 | −0.7027 | **−0.6793** | 3.24 | 14.190 |
| 16 | 25.9565 | −0.8372 | −0.8380 | −0.8313 | **−0.8066** | 1.65 | 14.064 |
| 20 | 26.5720 | −0.5474 | −0.5771 | −0.5172 | **−0.4968** | 4.98 | 14.402 |
| 24 | 26.7366 | −0.4320 | −0.4571 | −0.4067 | **−0.4565** | 1.02 | 14.301 |
| 28 | 26.9019 | −0.3618 | −0.3979 | −0.3257 | **−0.3870** | 1.36 | 14.510 |
| 32 | 26.9168 | −0.3415 | −0.3463 | −0.3366 | **−0.3373** | 1.33 | 14.565 |
| 36 | 26.7910 | −0.4131 | −0.4238 | −0.4025 | **−0.4141** | 1.41 | 14.506 |
| 40 | 26.8806 | −0.3751 | −0.3855 | −0.3643 | **−0.3559** | 3.27 | 14.576 |
| 44 | 27.0917 | −0.2511 | −0.2608 | −0.2411 | **−0.2538** | 1.37 | 14.681 |

Reading of the table.

* **The cone costs nothing.** cone-min $\asymp\lambda_1$ (factor $\le5$, usually $\le1.7$), and
  $\lambda_2/\lambda_1\le e^{0.87}=2.4$. The shortest vector of $\mathcal K_n$ *is* (up to a
  bounded multiple) a positive-cone vector. Compare with Sol's model $x+\alpha y\equiv0\ (d)$,
  where $\lambda_2/\lambda_1$ is unbounded and the cone forces the second-minimum scale: that
  degeneracy is absent at every $n\le44$.
* **The trend.** $n^{-1}\log(\text{cone-min})$ agrees with $\frac1{2n}\log\mathrm{covol}
  =\frac12(\kappa_n+E_1+E_2-2\sigma)$ to within $0.03$–$0.05$; it is negative at every $n$ but
  **drifts upward**, because $\kappa_n\uparrow\sigma$. If $\kappa_n=\sigma+o(1)$ the limit is
  exactly $F(k)=\frac{\log2}{2}(k_*-k)$: $-0.0169$ at $k=22.4$, $-0.537$ at $k=23.9$
  (the $k=23.9$ sweep is in the script history; at $n=38$ it reads $-0.970$, again $\approx$
  $\frac1{2n}\log\mathrm{covol}=-0.930$). So the *observed* $-0.25$ at $n=44$ is **not** the
  asymptotic rate; it is the pre-asymptotic bonus from $\sigma-\kappa_n\approx0.44$.
* Note $\delta=-\log|q_nG-p_n|/\log|q_n|\approx0.017$ is *tiny* — and irrelevant. For
  irrationality one needs only $0<|q_nG-p_n|\to0$ with $q_n,p_n\in\mathbb Z$; the quality
  exponent plays no role. This is why the cone route is worth taking seriously even at
  $F$ infinitesimally negative.

### 3.3 What Dilworth/pigeonhole actually gives (Sol's version)
Take the box $[0,L]^2$ in the scaled coordinates; it holds $N\approx L^2/\mathrm{covol}$ points of
$L_n$. A comparable pair gives a difference vector in the closed positive quadrant of size $\le L$.
By Dilworth one needs $N>|\text{max antichain}|$. A maximal antichain is a staircase, and its size
is bounded by the number of distinct first coordinates realised, i.e. $\approx L/\ell_1$ where
$\ell_1$ is the minimal positive first coordinate — which is $\asymp\lambda_1$ in the worst case.
So the criterion reads $L^2/\mathrm{covol}>L/\lambda_1$, i.e.
$$L>\frac{\mathrm{covol}}{\lambda_1}\ \asymp\ \lambda_2 .$$
**The pigeonhole bound is exactly the second minimum** — it does not beat ordinary Minkowski, and
in the degenerate regime ($\lambda_1\ll\lambda_2$) it is as bad as the naive bound. The
combinatorics adds nothing; the entire gain in §3.2 comes from the *measured* near-orthogonality
$\lambda_2\asymp\lambda_1\asymp\sqrt{\mathrm{covol}}$, which is a property of these particular rows.

## 4. The control, and why a decaying positive-cone form *would* be legitimate

`control.gp`: replace $G$ by $G^*=\mathrm{bestappr}(G,10^{320})=a/b$, $\log b=735.68$,
$|G-G^*|=e^{-1472.8}$. **Every column of the table is bit-identical** for $n\le38$. So the finite-$n$
data is *not* evidence of anything — exactly as in `CATALAN_AUDIT.md` §4(a).

But the *logical* status is different from the box-selected vector, and it is worth being precise
about why.

* For the box/Minkowski vector, the two needed facts — $\ell_n\ne0$ and $|\ell_n|\le e^{F n}$ —
  are never simultaneously available (audit §4(c)). One observes $\ell_n\ne0$; observation is not
  proof, and the observation is reproduced by $G^*$.
* For a positive-cone vector, $\ell_n>0$ is a **theorem** (Zudilin's and Nesterenko's integrals).
  Hence a *proved* upper bound $\ell_n\le e^{-cn}$, $c>0$, for infinitely many $n$ would give
  $0<|q_nG-p_n|\to0$ with $q_n,p_n\in\mathbb Z$, hence $G\notin\mathbb Q$ — with no observational
  step anywhere. The control does not refute this: if $G=a/b$ then every cone vector obeys
  $q_nG-p_n\ge1/b$, so the decay must stop; with $c=|F(22.4)|=0.0169$ and $\log b=735.68$ that
  happens only past $n\approx43\,500$, far beyond any computation. A proof of the decay for all
  large $n$ therefore *proves* $\operatorname{den}(G)\ge e^{cn}$ for all $n$ — which is the
  irrationality statement, obtained legitimately rather than assumed.
* The flip side, and the reason for scepticism: because the conclusion is equivalent to
  irrationality, **the required lattice-geometry lemma cannot be a generic one.** If $G=a/b$, the
  kernel direction $(r^N,-r^Z)$, $r^Z=aX_n-bY_n$, $r^N=aV_n-bU_n$, has scaled length
  $\asymp b\,e^{2F(k)n}$, which drops below $\sqrt{\mathrm{covol}}=e^{F(k)n}$ once
  $n>\log b/|F(k)|$; from then on $\lambda_1\ll\lambda_2$ and the cone minimum is pushed up to
  $\mathrm{covol}/\lambda_1\ge1/b$. So "$\text{cone-min}\asymp\sqrt{\mathrm{covol}}$ for all large
  $n$" is **false for every rational $G$**, and no proof of it can be covolume-only. This is the
  honest cost of the mechanism: it converts one impossible requirement (nonvanishing *and* size)
  into one hard requirement (size, with nonvanishing free), and the hard requirement still has to
  see something about $G$ that a rational cannot fake. It is a genuine reduction, not a proof, and
  it is not circular in form — an *explicit construction* of a short cone vector would discharge it.

## 5. What this reduces to

Non-vanishing is removed from the supercritical Catalan problem. What remains is:

* **Lemma P1 (arithmetic, already isolated).** $\kappa_n=\sigma(k)+o(n)$ for some $k>k_*=22.35129$,
  i.e. the exact 2-adic slopes assumed in `NesterenkoCross.lean`'s `NesterenkoPadicInputs`
  (`val_four_BN`, `val_tail`) hold for all $n$. Unchanged from `CATALAN_DIRECTIONAL.md` §4(b);
  measured deficit here $\sigma-\kappa_n=0.44$ at $n=44$, still decreasing.
* **Lemma P2 (geometric, new, and the whole remaining content).** For infinitely many $n$,
  $$\min\{\,s_Zc_Z\lambda_Z+c_N\lambda_N:\ c\in\mathcal K_n\cap\mathcal P\,\}\ \le\
  e^{o(n)}\sqrt{[\mathbb Z^2:\mathcal K_n]\,\lambda_Z\lambda_N}.$$
  Verified for $n\le44$ with constant $\le5.0$. Equivalently: the first minimum of the scaled
  Catalan lattice is attained (up to $e^{o(n)}$) in the positive quadrant. By §4 this cannot follow
  from the covolume; it must come from the structure of $\mathcal K_n$ — the 2-adic reduction
  modulo $T_n$, the shape of $X_n,Y_n,V_n,U_n$ modulo $D_{6n}^2$, or an explicit construction.
* **Lemma P3 (bookkeeping).** $\lambda_Z\lambda_N\,[\mathbb Z^2:\mathcal K_n]=e^{2F(k)n+o(n)}$ with
  $F(k)<0$, i.e. the source-form rates $E_1,E_2$ and $\log D_{6n}=6n+o(n)$. Standard, already in the
  papers.

P1+P2+P3 $\Rightarrow$ $0<q_nG-p_n\le e^{F(k)n+o(n)}\to0$ with integer $q_n,p_n$.
Lemma P2 is the only genuinely new object, and §4 shows it is of irrationality strength; it is
stated here so that the remaining difficulty is explicit, not because it is believed to be near.

## 6. Mechanism C: refuted

`character.gp`, $n=12$ and $n=20$, $\sim14\,600$ lattice vectors each.

| modulus | realised classes | classes containing both signs |
|---|---|---|
| 3 | 3 | 3 (100 %) |
| 5 | 5 | 5 (100 %) |
| 7 | 7 | 7 (100 %) |
| 16 | 16 | 16 (100 %) |
| 64 | 64 | 64 (100 %) |

A battery of invariants (Legendre symbols $\left(\frac{c_Z}{\ell}\right)$,
$\left(\frac{c_N}{\ell}\right)$, $\left(\frac{c_Zc_N}{\ell}\right)$ for $3\le\ell\le29$; parities;
$\chi_{-4}(c_Z+c_N)$) achieves chance accuracy; the only "predictors" at accuracy $1$ are
$\operatorname{sign}(c_Z)$ and $\operatorname{sign}(c_N)$, i.e. the cone itself — which is
Mechanism B, not an arithmetic invariant.

The reason is structural and needs no experiment: $\operatorname{sign}(c_Z\Lambda_Z+c_N\Lambda_N)$
is the indicator of an open half-plane through the origin with irrational-slope boundary, while a
congruence condition mod $m$ defines a coset of $m\mathcal K_n$, a full-rank sublattice; every such
coset meets both open half-planes in infinitely many points. Hence **no** congruence- or
character-valued invariant of $(c_Z,c_N)$ can determine the sign. The window law
$\ell^2P_m/Q_m\equiv\chi_{-4}(\ell)$ is a statement about the *rows*, not about the coefficient
vector, and is invariant under $c\mapsto-c$; it can carry no orientation. (Side observation: only
$\ell$ of the $\ell^2$ classes mod $\ell$ are realised — one more visible face of the absorption
theorem of `CATALAN_DIRECTIONAL.md` §3.)

## 7. Summary of the single most informative number

$$\frac1n\log\bigl(\text{best positive-cone }|q_nG-p_n|\bigr)\Big|_{k=22.4}:\quad
-0.91\,(8),\ -0.50\,(20),\ -0.39\,(28),\ -0.34\,(32),\ -0.25\,(44),$$
negative throughout, provably a *lower* bound on nothing (it is an achieved value with proved
positivity), and drifting toward the asymptotic $F(22.4)=-0.0169$ from below as $\kappa_n\to\sigma$.
The achievable $\delta$ with **guaranteed non-vanishing** is therefore $\approx0.0012$ asymptotically
— useless as a quality exponent, sufficient in principle as an irrationality mechanism, and
contingent entirely on Lemmas P1 and P2.

Scripts: `lattice/catalan_positivity/{rows_common.gp, diag.gp, integrals.gp, cone_honest.gp,
cone_sweep.gp, control.gp, character.gp}`; each run as
`cat rows_common.gp X.gp > run.gp; gp -q run.gp`, total runtime $<30$ s.
