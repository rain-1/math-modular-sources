# The two unidentified sources: Cooper's $s_{18}$ and Zudilin's Catalan row

*Answers to the two "verified only because their SOURCES are unidentified" entries
of `EULER_CRITERION.md` §4.1.  Scripts: `lattice/sources_s18_zud/`
(`s18_setup.gp`, `s18_order.gp`, `s18_order2test.gp`, `s18_struct.gp`, `s18_mf2.gp`,
`s18_pole.gp`, `s18_padic.gp`, `zud_row.gp`, `zud_check.gp`, `zud_ode.gp`,
`zud_ode2.gp`, `zud_padic.gp`, `zud_beukers.gp`, `zud_moving.gp`, `zud_mirror.gp`,
`ode_fit.gp`).  All PARI computations are exact over $\Q$ unless stated.*

---

## 0. Verdicts

1. **Cooper $s_{18}$ has no Eisenstein source, and Theorem F does not apply to it.**
   Its source $\Phi=F\,\thq t$ is *genuinely meromorphic*: it has a triple pole at the
   CM point $\tau_0=\tfrac{3+i}{6}$ (disc $-36$), the fixed point of the Atkin–Lehner
   involution $W_9$ on $\Gamma_0(18)$, where $x_{18}=-3$ and $t_{18}=\infty$.  §1–§3.
2. **The census entry for $s_{18}$ is numerically right and structurally wrong.**
   $\xi_3=\tfrac12\zeta_3(2)$ is re-verified here independently to $3^{219}$
   (§4), but the recorded explanation "$(\psi,\varphi)=(\chi_{-3},\mathbf 1)$,
   $Q(2)=-1$" **cannot** be correct: $s_{18}$ is a **weight $w=2$** row (minimal
   operator of order $3$, $F$ of weight $2$, $\Phi$ of weight $4$, $\Theta=\thq^{-3}\Phi$
   — all proved in §1), for which Theorem F would force
   $\kappa_3=\tfrac12L_3(3,\psi)$, i.e. $\tfrac12\zeta_3(3)$; and
   $v_3(\xi_3-\tfrac12\zeta_3(3))=-1$.  So $s_{18}$ realises a **weight drop**:
   a $w=2$ row whose $3$-adic (and archimedean) period is a $w=1$ value.
   The drop is caused by the principal part, and is not explained by any theorem
   we have. **Open.**
3. **The good news for the meromorphic extension:** the pole sits at $t=\infty$, hence
   outside every affinoid $\{|t|_3\le\rho\}$; equivalently $\tau_0$ has CM by an order
   in $\Q(i)$ and $3$ is inert in $\Q(i)$, so $\tau_0$ is **supersingular at $3$**.
   The principal part therefore lives in the supersingular locus and does **not**
   obstruct rigid-analyticity of $H_\xi=F(\Theta-\xi)$ on a strict neighbourhood of
   the ordinary locus.  This is why a $3$-adic limit exists at all.  §3.3.
4. **The archimedean Eisenstein-source proof of $s_{18}\to\tfrac12L(2,\chi_{-3})$ is
   not available**: the decomposition (§2) exhibits the *holomorphic* weight-$4$ factor
   of $\Phi$ as
   $\tfrac13\bigl(4E_4(6\tau)-E_4(3\tau)\bigr)+4\eta(3\tau)^8-16\eta(6\tau)^8$
   — a **trivial-character** $\zeta(s)\zeta(s-3)$ Eisenstein part plus the level-$9$
   weight-$4$ **CM newform**.  Neither factor carries $L(s,\chi_{-3})$, so no
   Eisenstein argument produces $\tfrac12L(2,\chi_{-3})$.  **Open.**
5. **Zudilin's Catalan row is not a modular Apéry row at all.**  The minimal linear
   ODE annihilating $\sum_m Q_mt^m$ has **order 4** (no annihilator of order $\le3$
   with polynomial coefficients of degree $\le10$), exponents $\{0,0,\tfrac12,\tfrac12\}$
   at $t=0$, and leading symbol $(1-11T-T^2)(1+10T+T^2)$, $T=16t$.  Its canonical
   coordinate $t(q)$ is violently non-integral (denominators $3,15,315,\dots$, and
   primes $3,5,7,11,13,\dots$, so no rescaling repairs it).  There is therefore **no**
   $(\Gamma,t,F)$ with $F$ of weight $1$, and the census entry
   "$(\psi,\varphi)=(\chi_{-4},\mathbf 1)$, $Q(2)=-2$" is a fit with no source. §5.
6. **The Beukers link is real, but it does not prove the alignment.**
   $Q_m=q_m(\tfrac12-m)$ **exactly**, where $q_n(x)$ are Beukers' Padé denominators for
   $\Theta(x)$: the 5:8 paper's "Rivoal's moving Padé specialisation is exactly
   Zudilin's sequence" is **confirmed for the denominators** (§6.2).  But the
   numerators are *not* the specialised Padé numerators, and the moving-point passage
   claimed in that paper is **invalid as written**: we exhibit the failure explicitly
   ($p_m(\tfrac12-m)/Q_m$ has $v_2$ of successive differences $\equiv4$, i.e. it does
   **not** converge $2$-adically, while $P_m/Q_m$ converges at slope $8$).  §6.4.
   $\xi_2^{\rm Zud}=\zeta_2(2)=2\,\xi_2^{\mathbf E}$ is re-verified to $2^{398}$.
   ~~remains **numerical**. **Open.**~~ **CLOSED — proved in `ZUDILIN_2ADIC.md`** (§6.3–§6.5 below are superseded).
7. **Normalisation, settled by Beukers himself.**  There is no $p$-adic Catalan
   constant: $L_2(2,\chi_{-4})\equiv0$ (odd character), and Calegari's "$2$-adic
   Catalan constant" *is* $\zeta_2(2)$.  This confirms the convention already used in
   `EULER_CRITERION.md` §4.1. §6.1.

---

## 1. Cooper $s_{18}$: the row is $w=2$, $r=3$

Cooper's parametrisation (Ramanujan J. **29** (2012); transcribed in
`packages/phase 1 zeta math package/sporadic_eisenstein_cooper_research_notes.txt`
§8):
$$x=x_{18}=\frac{\eta(3\tau)^4\eta(6\tau)^4}{\eta(\tau)^2\eta(2\tau)^2\eta(9\tau)^2\eta(18\tau)^2}
=q^{-1}+2+7q+\cdots,\qquad
t=t_{18}=\frac{x}{(x+3)^2},$$
$$F=F_{18}=\tfrac14\bigl(18E_2(18\tau)-9E_2(9\tau)-12E_2(6\tau)+6E_2(3\tau)+2E_2(2\tau)-E_2(\tau)\bigr)\in M_2(\Gamma_0(18)).$$
($\sum c_d=0$ in $\sum c_d\,dE_2(d\tau)$, so $F$ is a genuine holomorphic weight-$2$ form;
`s18_mf2.gp` writes it in $M_2(\Gamma_0(18))$, $\dim=7$.)

**Verified exactly** (`s18_order.gp`, $q$-precision $46$):

* $F_{18}-\sum_{n\le20}A_nt_{18}^n=O(q^{21})$ with $A_0=1,A_1=6$ and Cooper's
  recurrence $(n+1)^3A_{n+1}=2(2n+1)(7n^2+7n+3)A_n-12n(16n^2-1)A_{n-1}$;
  $A=1,6,54,564,6390,76356,948276,\dots$
* $\thq\log t_{18}=F_{18}\sqrt{(1-16t)(1-12t)}$ to $O(q^{26})$.
* **Order of integration.**  $[t^n]\bigl(F\thq^{-k}\Phi\bigr)$ equals
  $B_n$ ($B_0=0,B_1=1$) **iff $k=3$**: for $k=3$ the two agree
  ($0,1,\tfrac{51}4,\tfrac{313}2,\tfrac{31245}{16},\dots$), for $k=1,2$ they do not.
  Hence $r=3$, $w=2$, $\Theta=\thq^{-3}\Phi$, $\wt\Phi=4$.
* **Minimal operator has order 3, not 2** (`s18_order2test.gp`, $50$–$90$ terms):
  no annihilator of order $2$ with $\deg p_i\le10$ exists; the unique order-$3$,
  $t$-degree-$2$ annihilator is
  $$L_{18}=\theta^3-2t(2\theta+1)(7\theta^2+7\theta+3)+12t^2(\theta+1)(4\theta+3)(4\theta+5),$$
  leading symbol $1-28t+192t^2=(1-16t)(1-12t)$.
  **Correction to `paper/sections/02_sources.tex`, Def. 2.1**: the claim that Cooper's
  $d\neq0$ rows "are of operator order two despite the cubic-looking recurrence"
  is false for $s_{18}$.

---

## 2. The exact structure of the source

Write $Z:=\thq\log x_{18}\in M_2(\Gamma_0(18))$ (an $E_2$-combination with
$\sum c_d=0$; constant term $-1$).  All of the following are verified to $O(q^{30})$
or better (`s18_struct.gp`, `s18_mf2.gp`):

$$\boxed{\ Z=-F\sqrt{1-16t},\qquad \frac{x-3}{x+3}=\sqrt{1-12t},\qquad
t=\frac1{\,x+9/x+6\,},\ }$$

the last identity showing that $t$ is invariant under the Atkin–Lehner involution
$x\mapsto 9/x$ (i.e. $t$ is a coordinate on the $W_9$-quotient of $X_0(18)$, while
$x$ is the Hauptmodul of $\Gamma_0(18)$ itself).  Consequently

$$\boxed{\ \Phi=F\,\thq t=W\cdot t\sqrt{1-12t}=W\cdot\frac{x(x-3)}{(x+3)^3}
=W\Bigl(\frac1{x+3}-\frac9{(x+3)^2}+\frac{18}{(x+3)^3}\Bigr),\ }$$
$$W:=F^2\sqrt{1-16t}=-F\,\thq\log x_{18}\ \in M_4(\Gamma_0(18)).$$

$W=1+4q-16q^2-80q^3-32q^4-400q^6+80q^7+128q^8-2240q^9-\cdots$, and `mftobasis`
(exact, residual $0$ through $q^{40}$) gives the **complete decomposition**

$$\boxed{\ W=\underbrace{\tfrac13\bigl(4E_4(6\tau)-E_4(3\tau)\bigr)}_{\text{Eisenstein}}
\;+\;\underbrace{4\,\eta(3\tau)^8-16\,\eta(6\tau)^8}_{\text{cuspidal}}\ }$$

with $E_4=1+240\sum\sigma_3(n)q^n$ and $\eta(3\tau)^8$ the CM newform of
$S_4(\Gamma_0(9))$ (CM by $\Q(\sqrt{-3})$, $a_2=a_3=0$).

* The Eisenstein part has **trivial characters**: $L=\tfrac{240}3(4\cdot6^{-s}-3^{-s})\zeta(s)\zeta(s-3)$,
  i.e. $P(V)=80\,V_3(4V_2-1)$, entirely inside the image of $V_3$.  In particular
  $\mathcal E_3(s)=1-3^{-s}$ does **not** divide $P(s)$: criterion (a) of Theorem F
  fails for this factor.
* The cuspidal part has $a_3=0$, so the cuspidal clause would give $\kappa_3=0$.
* Neither factor involves $L(s,\chi_{-3})$.  This is the concrete reason why the
  Eisenstein-source route to $\xi_\infty=\tfrac12L(2,\chi_{-3})$ (task 1, second half)
  **does not close**.

---

## 3. The principal part

### 3.1 Where the pole is
$\Phi$ is *not* holomorphic: its coefficients grow exponentially,
$$\Phi=q-10q^2+54q^3-236q^4+870q^5-3132q^6+\cdots,\qquad
|c(m)|^{1/m}\to 2.8497\ldots$$
(measured $2.858$ at $m=43$).  The pole is at $x_{18}=-3$, i.e. $t_{18}=\infty$, i.e.
(`s18_pole.gp`, $50$ digits)
$$\boxed{\ \tau_0=\tfrac12+\tfrac i6=\tfrac{3+i}{6},\qquad x_{18}(\tau_0)=-3.000\ldots,\qquad
q_0=e^{2\pi i\tau_0}=-e^{-\pi/3}=-0.3509198\ldots,\ }$$
and $1/|q_0|=e^{\pi/3}=2.84965\ldots$ — exactly the measured growth rate, and the sign
alternation is $\operatorname{Re}\tau_0=\tfrac12$.  $\tau_0$ is the fixed point of $W_9$
($9/x=x$ at $x=-3$); it is a CM point of discriminant $-36$ (order $\Z[3i]\subset\Q(i)$).
Since $t$ has a double zero denominator $(x+3)^2/x$ there, $t$ has a double pole and
$\Phi=F\thq t$ a **triple** pole at $\tau_0$ and its $\Gamma_0(18)$-translates.
The other $W_9$ fixed point $x=3$ is $t=1/12$, the *sub-dominant* singularity of $L_{18}$;
the dominant one $t=1/16$ is $x\in\{1,9\}$.

### 3.2 Consequence for Theorem A
Theorem A ("$\Phi$ is a holomorphic Eisenstein series, cuspidal projection zero")
is **false** for $s_{18}$ — not by a small margin but qualitatively.  The correct
statement is the boxed factorisation of §2: *holomorphic weight-4 form times a
degree-3 rational function of the Hauptmodul with a triple pole at a CM point.*
This is the precise form of the "meromorphic source with controlled principal part"
proposed in the Cooper notes §12; the notes' hope that the principal part is a
"total derivative plus boundary data" is not what happens.

### 3.3 Why a $3$-adic limit nevertheless exists (the affinoid question)
The task's question — *do the poles lie outside the relevant affinoid?* — has a clean
affirmative answer, twice over:

* **In the $t$-line.**  The pole is at $t=\infty$.  Hypothesis (c) of Theorem F puts
  the ordinary locus inside $\{|t|_3\le1\}$, and $\infty\notin\{|t|_3\le\rho\}$ for every
  $\rho<\infty$.  So $H_\xi=F(\Theta-\xi)$ can be rigid-analytic on a strict
  neighbourhood of $D_\infty$ despite the pole.
* **On the modular curve.**  $\tau_0$ has CM by an order in $K=\Q(i)$, and
  $3\equiv3\pmod 4$ is **inert** in $K$; a CM point whose field is non-split at $p$
  reduces to a **supersingular** point.  Hence the principal part is supported in the
  supersingular locus, which is disjoint from every strict neighbourhood of the
  ordinary locus.  (Caveat: $3$ divides the conductor $f=3$ of $\Z[3i]$; the
  supersingularity conclusion uses only $K$, not the order, and the $t$-line argument
  above is independent of it.)

This is a genuinely positive structural finding: **meromorphy at supersingular CM points
is harmless for the $p$-adic Apéry mechanism.**  It also predicts that $s_{18}$ should
have *no* obstruction at $p=3$ and none of the usual protection at $p=2$ — and indeed
`SLOPE_CENSUS.md` §1 records slope $1$ at $p=3$ and no slope at $2,5,7$.

### 3.4 What is *not* explained: the weight drop
$\Theta=\thq^{-3}\Phi$ has weight $-w=-2$ (Bol/Coleman bookkeeping: $F$ has weight $2$,
$H$ must have weight $0$).  The constant term of a *Coleman–Mazur Eisenstein* family at
weight $-2$ is $\tfrac12L_p(3,\psi\omega^{-2})$, i.e. a $\zeta_3(3)$-type number; the
measured $\xi_3$ is $\tfrac12\zeta_3(2)$.  Both statements are exact (§4), so the
overconvergent avatar of $\Theta_{s_{18}}$ is **not** an oldform combination of an
Eisenstein family — it is a general overconvergent form of weight $-2$ whose constant
term happens to be $\tfrac12\zeta_3(2)$.  Nothing in the present theory predicts it.

*The only principled derivation currently available* is by rigidity: `SLOPE_CENSUS.md`
§1 measures exact $3$-adic alignment of $s_{18}$ with Zagier $\mathbf B$ and
$\mathbf C$ ($v_3\approx n$ out to $n=300$); given the *proved* value
$\xi_3^{\mathbf C}=\tfrac12\zeta_3(2)$ (Theorem F, $w=1$) the rigidity lemma of
`RIGIDITY_PROOF.md` transfers it to $s_{18}$.  That reduces one numerical statement to
another, and is recorded as the honest status.

One small exact curiosity, possibly a handle: $\thq^{-1}\Phi=\sum_m\frac{c(m)}m q^m$
is **integral** — $1,-5,18,-59,174,-522,1520,-4421,\dots$ — verified for $m\le43$
(`s18_struct.gp`); equivalently $\sum_n A_nt^{n+1}/(n+1)\in\Z[[q]]$.

---

## 4. Independent re-verification of $\xi_3(s_{18})$

`s18_padic.gp`, $N=339$, $L_3$ from `lattice/euler_criterion/lp.gp` at precision $3^{220}$:

| quantity | value |
|---|---|
| Cauchy $v_3(\xi_N-\xi_{N-1})$ | $328$ (slope $\approx1$) |
| $v_3\bigl(\xi_N-\tfrac12\zeta_3(2)\bigr)$ | $\mathbf{219}$ (= full available precision) |
| $v_3\bigl(\xi_N-\tfrac12\zeta_3(3)\bigr)$ | $-1$ |
| $v_3\bigl(\xi_N+\tfrac12\zeta_3(2)\bigr)$ | $-1$ |

So $\xi_3=\tfrac12\zeta_3(2)$, confirming `EULER_CRITERION.md` §4.1 row 12 **as a value**
and refuting the weight-$2$ (i.e. $\zeta_3(3)$) prediction that Theorem F would make for
a $w=2$ row.

---

## 5. Zudilin's Catalan row: no modular structure

Zudilin, *An Apéry-like difference equation for Catalan's constant*
(arXiv:math/0201024).  With $p(n)=20n^2-8n+1$,
$q(n)=3520n^6+5632n^5+2064n^4-384n^3-156n^2+16n+7$,
$$(2n{+}1)^2(2n{+}2)^2p(n)\,u_{n+1}-q(n)\,u_n-(2n{-}1)^2(2n)^2p(n{+}1)\,u_{n-1}=0,$$
$Q_0=1,Q_1=\tfrac74$, $P_0=0,P_1=\tfrac{13}8$; $P_m/Q_m\to G$ (checked: $43$ digits at $m=20$).

**Exact $2$-adic denominator law** (`zud_check.gp`, $m\le14$, and used by the 5:8 paper):
$$v_2(Q_m)=-4m+2s_2(m),\qquad s_2=\text{binary digit sum};$$
equivalently $16^mQ_m/\binom{2m}m^2$ is a $2$-adic unit.  Also, exactly,
$v_2(P_m)=v_2(Q_m)-1$ for $1\le m\le14$.  This is sharper than Zudilin's
$2^{4m+3}Q_m\in\Z$ and is the source of $\kappa_2=4$, $\sigma_2=2\kappa_2=8$ ($c=-1$).

**Holonomic rank** (`zud_ode.gp`, `zud_ode2.gp`; $U_m:=16^mQ_m\in\Z$, $70$ terms):

* No annihilator of order $2$ or $3$ with $\deg p_i\le10$. (**Order $\ge4$.**)
* A unique order-$4$, $t$-degree-$4$ annihilator exists:
  $p_4=262144t^4+344064t^3+112640t^2+64t-4$, …, with
  $$p_4(t)=(256t^2+160t+1)(256t^2+176t-1)\ \Longrightarrow\
  \text{symbol}\;=\;(1-11T-T^2)(1+10T+T^2),\quad T=16t.$$
  The first factor is **exactly Apéry's $\zeta(2)$ / Zagier $\mathbf D$ symbol**
  $\lambda^2-11\lambda-1$ (whence Zudilin's characteristic roots
  $((1\pm\sqrt5)/2)^5$); the second, $\lambda^2+10\lambda+1$ (roots $-5\pm2\sqrt6$),
  is a new strand.
* Indicial polynomial at $t=0$: $-\theta^2(2\theta-1)^2$, exponents
  $\{0,0,\tfrac12,\tfrac12\}$ — half-integer exponents, consistent with a tensor
  product $L_{\rm MUM}\otimes L_{1/2}$ of two rank-$2$ systems (exponents $\{0,0\}$ and
  $\{0,\tfrac12\}$).  This is the precise sense in which "the half-integer
  hypergeometric parameters" of the task statement are a quadratic-transformation
  phenomenon.
* **Mirror map** (`zud_mirror.gp`): solving $L(h)=-(\partial L/\partial\theta)y_0$ for
  the exponent-$0$ log solution and forming $q=t\exp(h/y_0)$, $t(q)=\operatorname{rev}(q)$:
  $$t(q)=q-120q^2+\tfrac{38276}3q^3-\tfrac{19763392}{15}q^4+\tfrac{43158218618}{315}q^5-\cdots$$
  The denominators contain $3,5,7,11,13,\dots$, so **no rescaling $t\mapsto\lambda t$,
  $q\mapsto\mu q$ makes $t(q)$ integral**.  There is no canonical coordinate, hence no
  Hauptmodul, hence no $(\Gamma,t,F)$.

**Verdict.** $\sum_mQ_mt^m$ is *not* a modular form in a Hauptmodul of any level;
in particular it is not a level-$8/16/32$ pullback of $\mathbf E$, and no
binomial/Hadamard transform of $\mathbf E$'s $a_n$ can produce it (those preserve
order $2$ / order $3$, not order $4$ with $\{0,0,\frac12,\frac12\}$).
The Theorem-F entry for this row in `EULER_CRITERION.md` §4.1 must be withdrawn as a
*derivation*; the value stands as a measurement.

**$2$-adic value, re-verified** (`zud_padic.gp`, $N=59$, precision $2^{400}$):
Cauchy $v_2(\xi_N-\xi_{N-1})=451$ (slope $\approx8$),
$v_2(\xi-\zeta_2(2))=\mathbf{398}$ (full precision), $v_2(\xi-\tfrac12\zeta_2(2))=-2$.
So $\xi_2^{\rm Zud}=\zeta_2(2)=2\xi_2^{\mathbf E}$, matching $G:\tfrac12G$.

---

## 6. Beukers' $p$-adic Padé construction and the alignment

### 6.1 What Beukers proves
F. Beukers, *Irrationality of some $p$-adic $L$-values*, Acta Math. Sinica **24** (2008)
663–686 (= arXiv:math/0603277; author version
`webspace.science.uu.nl/~beuke106/padicL.pdf`).  In the field of formal Laurent series
in $1/x$ he studies Stieltjes' three continued fractions for
$$R(x)=\sum_{n\ge0}B_n(-1/x)^{n+1},\quad
\Theta(x)=\sum_{n\ge0}t_n(-1/x)^{n+1}\ (t_n=(2^{n+1}-2)B_n),\quad T(x)=R'(x),$$
with $\Theta(x)=R(x/2)-2R(x)$ (so the 5:8 paper's $\mathcal T=2R(x)-R(x/2)=-\Theta$;
its normalisation is consistent).  For $\Theta$:
$$\Theta(x)=\cfrac1{x^2-x+a_1-\cfrac{b_1}{x^2-x+a_2-\cdots}},\qquad
a_n=2n^2-2n+1,\ b_n=n^4,$$
$$\boxed{\ (n+1)^2w_{n+1}=(2n(n+1)+1-x+x^2)w_n-n^2w_{n-1}\ }$$
with $q_0=1,q_1=x^2-x+1$ and $p_0=0,p_1=1$, and
$$q_n(x)=\sum_{k=0}^n\binom nk\binom{-x}k\binom{k-x}k,\qquad
p_n(x)-q_n(x)\Theta(x)=O(x^{-2n-2}).$$
Substituting $x=a/F$ with $p\mid F$, $p\nmid a$ makes everything converge $p$-adically;
Prop. 7.1 gives $|p_n(a/F)-\Theta_p(a/F)q_n(a/F)|_p\le p^2n^2p^{-2n(r+1/(p-1))}$,
$|F|_p=p^{-r}$ — **uniformly in $a$**.  Dictionary (Prop. 5.1):
$$\Theta_2(1/2)=-8\zeta_2(2),\quad \Theta_2(1/6)=-40\zeta_2(2),\quad
\Theta_2(1/4)=-16L_2(2,\chi_8),\quad \Theta_3(1/3)=-\tfrac{27}2\zeta_3(2).$$
Cor. 7.3: **$\zeta_2(2)$, $\zeta_3(2)$, $L_2(2,\chi_8)$ are irrational**; and (Thm. 11.2,
Cor. 11.4) $\zeta_2(3)$, $\zeta_3(3)$, $\zeta_5(3)\pm L_5(3,\chi_5)$, $\zeta_2(3)\pm L_2(3,\chi_8)$.
He notes explicitly that the convergents at $x=\tfrac12$ **are** Calegari's approximations,
and that $\Theta_3(1/6)=-36L_3(2,\chi_{12})$ is *not* reachable by his condition (A)
(so $L_3(2,\chi_{12})$ is not among his irrationality results; note our row $\mathbf F$
carries $\tfrac12L_2(2,\chi_{12})$ at $p=2$, a different number).

**Normalisation, decisively.**  Beukers' introduction: Calegari's "irrationality of the
$2$-adic Catalan constant" is *not* about $L_2(2,\chi_4)$, which "vanishes identically"
(odd character); it is $\zeta_2(2)$.  This is exactly the convention of
`EULER_CRITERION.md` §4.1 note ($\Lambda_2=\zeta_2(2)$ because $\omega=\chi_{-4}$ at $p=2$)
— **confirmed against the literature.**

### 6.2 The moving specialisation: confirmed for denominators
Beukers remarks that $x=-n$ recovers Apéry's $\zeta(2)$ numbers and that
$x=-n+\tfrac12$ "gives numbers which play a role in approximations of Catalan's
constant" (Rivoal, *Nombres d'Euler, approximants de Padé et constante de Catalan*).
We verify this **exactly** (`zud_check.gp`, `zud_beukers.gp`):
$$q_n(-n)=1,3,19,147,1251,11253,\dots=\text{Apéry }\zeta(2),\qquad
\boxed{\ q_n(-n+\tfrac12)=Q_n\ }\ (n=0,\dots,10,\ \text{exact}),$$
$Q_n=1,\tfrac74,\tfrac{649}{64},\tfrac{19471}{256},\dots$.  So the 5:8 paper's
identification is **true for the denominator sequence**.


> **SUPERSEDED (see `ZUDILIN_2ADIC.md`).** $\xi_2^{\rm Zud}=\zeta_2(2)$ is now **proved**.
> The correction missing below is the exact identity
> $P_m=\tfrac{(-1)^m}8p_m(\tfrac12-m)+G_mq_m(\tfrac12-m)$, $G_m=\sum_{j\le m}\frac{(-1)^{j-1}}{(2j-1)^2}$;
> the moving target drifts by exactly $8(-1)^{m+1}G_m$, which $P_m$ already carries.
> Also proved there: $v_2(Q_m)=-4m+2s_2(m)$ and $v_2(\zeta_2(2)-P_m/Q_m)=8m-1-4s_2(m)$.
> The verdicts of §1–§5 (no modular source for the Zudilin row, the $s_{18}$ analysis) are unaffected.

### 6.3 …but not for the numerators
$P_n\ne p_n(-n+\tfrac12)$: e.g. $p_1=1$ vs $P_1=\tfrac{13}8$; $p_2(-\tfrac32)=\tfrac{35}{16}$ vs
$P_2=\tfrac{10699}{1152}$.  The difference satisfies $v_2\bigl((P_m-p_m(x_m))/Q_m\bigr)=-1$
for **all** $2\le m\le32$ — a bounded, never-vanishing correction.

### 6.4 The moving-point passage is invalid as written
The 5:8 paper argues: *"At the moving point $x_m=\tfrac12-m$, Rivoal's exact Padé identity
and the shift equation show that $P_m/Q_m$ converges to the single constant $\zeta_2(2)$.
Beukers's Padé estimate is uniform for the odd numerator $1-2m$, so the moving-point
passage is legitimate."*  Uniformity is real (§6.1) but insufficient, because the
**target moves**: $1/x_m-2=4m/(1-2m)$, so
$$|\Theta_2(x_m)-\Theta_2(\tfrac12)|_2\ \asymp\ 2^{-2-v_2(m)},$$
which is bounded away from $0$ along odd $m$.  Uniformity therefore only gives
$p_m(x_m)/Q_m\approx\Theta_2(x_m)$, and that has **no limit**.  Numerically
(`zud_moving.gp`, $m=2,\dots,32$):

| | $v_2$ of successive differences |
|---|---|
| $P_m/Q_m$ (Zudilin) | $3,11,15,27,31,39,43,59,\dots,227$ — slope $\approx8$, converges |
| $p_m(x_m)/Q_m$ (Beukers, moving) | $4,4,4,4,4,4,\dots,4$ — **does not converge** |

So the two divergences (moving target, and the correction $P_m-p_m(x_m)$) must cancel
exactly; the proof of that cancellation is missing, and it is the whole content of the
claimed corollary.  The shift equation $\Theta(x+1)+\Theta(x)=-2/x^2$ makes
$\Theta_2(x_m)=(-1)^m\Theta_2(\tfrac12)+(\text{explicit elementary rational})$, so the
repair is plausible — but it must be carried out with the elementary terms, and it must
explain the sign $(-1)^m$ and the exact factor $2$ (not $\tfrac12$).

### 6.5 Status of the alignment $\xi_2^{\rm Zud}=2\,\xi_2^{\mathbf E}$
* $\xi_2^{\mathbf E}=\tfrac12\zeta_2(2)$: **proved** by Theorem F modulo hypothesis (b)
  (`EULER_CRITERION.md` §4.1, row $\mathbf E$).
* $\xi_2^{\rm Zud}=\zeta_2(2)$: **measured** to $2^{398}$ here, $2^{3015}$ in the census.
  Not proved: Zudilin's row has no modular source (§5) and the Beukers route has the gap
  of §6.4.
* Hence the "common-period lemma" of the 5:8 paper — and with it the $5{:}8$ diagonal —
  currently rests on a numerical identification on the Zudilin side.
* Beukers *does* supply the irrationality of $\zeta_2(2)$ itself, so the constant is not
  in doubt; only the identification of Zudilin's $2$-adic limit with it is.

---

## 7. What to change in the paper / ledger

1. `paper/sections/02_sources.tex`, Def. 2.1: drop "of operator order two" for Cooper's
   $d\neq0$ rows; $s_{18}$ is order $3$ (§1).
2. `EULER_CRITERION.md` §4.1: mark the $s_{18}$ and Zudilin rows **value-only**; delete
   the conjectural $(\psi,\varphi)$/$Q(2)$ columns for both.  For $s_{18}$ note that the
   Theorem-F prediction for its true weight is $\tfrac12\zeta_3(3)$ and is **refuted**.
3. `EULER_CRITERION.md` §5 "Open": replace "the identification of the Eisenstein sources
   of $s_{18}$ and of Zudilin's row" by "**there are none**" (§2, §5), and add the two
   new open problems: the weight drop (§3.4) and the moving-point passage (§6.4).
4. Add the positive lemma of §3.3 (poles at supersingular CM points are harmless) to the
   Theorem-F section as the first step of a meromorphic extension.
