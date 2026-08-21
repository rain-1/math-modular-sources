# The period identity behind Conjecture D for (C, F) at p = 3

**Status.** This document does *not* prove anything new. It (a) makes the
archimedean side of the reduction in `RIGIDITY_PROOF.md` completely explicit and
verifies it to 37 digits, (b) formulates the 3-adic twin as a statement about the
single curve C, and (c) tests that statement pointwise at rational points to
3-adic precision $3^{1600}$. Scripts: `lattice/padic_period/`.

Recall the reduction (review note of `RIGIDITY_PROOF.md`): the transfer identity
$$B^{\mathbf F}-\tfrac54\xi A^{\mathbf F}=G\,R^{\mathbf C}_\xi(u)+\tfrac14G^2R^{\mathbf C}_\xi(v)$$
is linear in $\xi$ and holds for *every* $\xi$, so the input **H3**
($w$-invariance of the right-hand side at $\xi=\xi_3^{\mathbf C}$) is *equivalent*
to the theorem $\xi_3^{\mathbf F}=\tfrac54\xi_3^{\mathbf C}$, not an ingredient of
an independent proof.

---

## 1. The involution, made explicit  (`involution.gp`)

$$\boxed{\ \gamma=\begin{pmatrix}8&-3\\12&-4\end{pmatrix},\qquad
\gamma\tau=\frac{8\tau-3}{12\tau-4},\qquad \det\gamma=4,\qquad
\gamma^2=4\begin{pmatrix}7&-3\\12&-5\end{pmatrix}\in 4\,\Gamma_0(12).\ }$$

So $\gamma$ is the Atkin–Lehner involution $W_4$ of $X_0(12)$. Normalised to
$\det=1$ it is $\tfrac12\gamma=\begin{psmallmatrix}4&-3/2\\6&-2\end{psmallmatrix}$,
with automorphy factor
$$j(\tau):=c\tau+d=6\tau-2 .$$
It was found by solving the two coset conditions
$\gamma\in\Gamma_0(6)W_2^{(6)}\delta_2$ and $\delta_2\gamma\in\Gamma_0(6)W_2^{(6)}$,
where $W_2^{(6)}=\begin{psmallmatrix}2&-1\\6&-2\end{psmallmatrix}$ is the
Atkin–Lehner involution of $X_0(6)$ acting on the hauptmodul by
$t_{\mathbf C}\mapsto\sigma(t_{\mathbf C})$, $\sigma(t)=\dfrac{t-1}{9t-1}$
(this last fact is itself checked numerically), and $\delta_2:\tau\mapsto2\tau$.

**VERIFIED to $10^{-37}$** at three values of $\tau$:
$$t_{\mathbf C}(\gamma\tau)=\bar u,\qquad t_{\mathbf C}(2\gamma\tau)=\bar v,\qquad
t_{\mathbf F}(\gamma\tau)=t_{\mathbf F}(\tau),$$
$$F_{\mathbf F}(\gamma\tau)=j(\tau)\,F_{\mathbf F}(\tau),\qquad
\Phi_{\mathbf F}(\gamma\tau)=j(\tau)^3\,\Phi_{\mathbf F}(\tau).$$
That is: $F_{\mathbf F}\in M_1(\Gamma_0(12),\chi_{-3})$ and
$\Phi_{\mathbf F}\in M_3(\Gamma_0(12),\chi_{-3})$ both have $W_4$-eigenvalue $+1$.
(NB: $\gamma/2$ has trace $2$, i.e. it is parabolic in $SL_2(\mathbb R)$ with the
double fixed point $\tau=1/2$; it is an involution of $X_0(12)$, not of $\mathbb H$.)

## 2. Proposition (archimedean $W$-period identity)  (`arch_period.gp`)

Let $\Theta_{\mathbf F}=D^{-2}\Phi_{\mathbf F}=\Theta_{\mathbf C}+\tfrac14V_2\Theta_{\mathbf C}$,
a weight $-1$ Eichler integral, and let
$$P_\gamma(\tau):=\bigl(\Theta_{\mathbf F}|_{-1}\gamma\bigr)(\tau)-\Theta_{\mathbf F}(\tau)
=(6\tau-2)\,\Theta_{\mathbf F}(\gamma\tau)-\Theta_{\mathbf F}(\tau).$$
Then $P_\gamma$ is a polynomial of degree $1$ (as Eichler predicts for weight $3$),
and
$$\boxed{\ P_\gamma(\tau)=\tfrac54\,\xi_\infty\,\bigl(j(\tau)-1\bigr)
\;-\;i\,\frac{\pi^2}{6\sqrt3}\,\bigl(j(\tau)+1\bigr),\qquad
j(\tau)=6\tau-2,\quad \xi_\infty=\tfrac12L(\chi_{-3},2). \ }$$
Equivalently $P_\gamma(\tau)=\tfrac54\xi_\infty(6\tau-3)-i\frac{\pi^2}{6\sqrt3}(6\tau-1)$.

**VERIFIED**: $|P_\gamma(\tau)-\text{RHS}|<3\cdot10^{-37}$ at
$\tau=0.37+0.31i,\ 0.30+0.26i,\ 0.42+0.35i,\ 0.36+0.40i,\ 0.55+0.29i,\ 0.25+0.5i,\ -0.1+0.45i$
(38-digit arithmetic, $q$-series to $O(q^{900})$; the coefficients were fitted from
two $\tau$'s and the other five are held out).

The two coefficients are the two critical values of $L(\Phi_{\mathbf F},s)=P_{\mathbf F}(s)L(\chi_{-3},s)\zeta(s-2)$,
$P_{\mathbf F}(s)=1-7\cdot2^{-s}-8\cdot4^{-s}$:
* $s=2$: $P_{\mathbf F}(2)=-\tfrac54$, $\zeta(0)=-\tfrac12$, so
  $L(\Phi_{\mathbf F},2)=\tfrac58L(\chi_{-3},2)=\tfrac54\xi_\infty$ — the **first** coefficient, exactly the term predicted in the task brief;
* $s=1$: $P_{\mathbf F}(1)=-\tfrac92$, $\zeta(-1)=-\tfrac1{12}$,
  $L(\chi_{-3},1)=\pi/(3\sqrt3)$, so $L(\Phi_{\mathbf F},1)=\pi/(8\sqrt3)\neq0$ and
  $\tfrac{\pi^2}{6\sqrt3}=\tfrac{4\pi}{3}L(\Phi_{\mathbf F},1)=\tfrac\pi2 L(\chi_{-3},1)$ — the **second** coefficient.

**Correction to the working hypothesis.** The brief conjectured
$L(\Phi_{\mathbf F},1)=0$ "by purification" and $P_\gamma(\tau)=\tfrac54\xi_\infty(j-1)$
outright. That is false: $L(\Phi_{\mathbf F},1)=\pi/(8\sqrt3)\ne0$ and the second
term is present. The $\xi_\infty$-term is exactly as predicted.

**Correction 2 (important).** The naive statement "$H=w(H)$ over $\mathbb C$", if
$w(H)$ is read as substitution $\tau\mapsto\gamma\tau$ with modular automorphy, is
**false**: a direct test of
$\bigl(F_{\mathbf F}\Theta_{\mathbf F}(\tau)-F_{\mathbf F}\Theta_{\mathbf F}(\gamma\tau)\bigr)/\bigl(F_{\mathbf F}(\tau)-F_{\mathbf F}(\gamma\tau)\bigr)$
gives a non-constant function of $\tau$, and with the correct weights the deficiency
is exactly $F_{\mathbf F}(\tau)P_\gamma(\tau)\neq0$. The reason is that
over $\mathbb C$ the branch values $\bar u,\bar v\to1$ as $x\to0$, far outside the
disc of convergence $|t|<1/9$ of $A^{\mathbf C},B^{\mathbf C}$: the second sheet is
reached only by analytic continuation *past the fold* $t=1/9$, and that continuation
is not given by $F_{\mathbf C}(\gamma\tau)$ (indeed
$F_{\mathbf C}(\gamma\tau)/(jF_{\mathbf C}(\tau))\notin\{\pm1\}$). So §4/Prop. 3 of
`RIGIDITY_PROOF.md` ("*Reason.* Over $\mathbb C$ …") is **not** a proof of H3;
its archimedean shadow is the period identity above, not a vanishing statement.

## 3. The 3-adic twin

3-adically the situation is genuinely different and genuinely better: $|\bar u|_3,|\bar v|_3\le1$,
while $R^{\mathbf C}_{\xi_3}$ converges on $|t|_3<3^{\sigma_3}=9$. So the second
sheet is *inside* the domain of convergence and no continuation is needed.

> **3-adic $W$-period identity (statement, OPEN).** Let $\xi\in\mathbb Q_3$ and
> $R^{\mathbf C}_\xi=B^{\mathbf C}-\xi A^{\mathbf C}$. For $x_0\in\mathbb Q_3$ with
> $|x_0|_3\le3^\delta$, $0<\delta<2$, let $(u,v)$, $(\bar u,\bar v)$ be the two points
> of the cover over $x_0$ and $G=\frac{1+3v}{1-3u}$, $\bar G=\frac{1+3\bar v}{1-3\bar u}$. Then
> $$\Xi_{x_0}(\xi):=\Bigl[G R^{\mathbf C}_\xi(u)+\tfrac14G^2R^{\mathbf C}_\xi(v)\Bigr]
> -\Bigl[\bar G R^{\mathbf C}_\xi(\bar u)+\tfrac14\bar G^2R^{\mathbf C}_\xi(\bar v)\Bigr]=0
> \quad\text{for }\xi=\xi_3^{\mathbf C}.$$

Two remarks make this the right "single-curve" object.

* **It only mentions C.** $\mathbf F$ has disappeared; the statement is a functional
  equation for the overconvergent function $R^{\mathbf C}_{\xi}$ under the
  correspondence $w$ on the $X_0(12)$-cover of the C-curve. This is exactly the shape
  Calegari's method wants: one curve, one overconvergent form, one $U_p$/AL relation.
* **It determines $\xi_3$.** $\Xi_{x_0}$ is affine in $\xi$; the truncations
  $\Xi_{x_0}^{(K)}(\xi)=c_0^{(K)}-c_1^{(K)}\xi$ have $v_3(c_1^{(K)})$ *bounded*
  (measured: $0$–$9$), so the single scalar equation $\Xi_{x_0}=0$ pins $\xi$ down.
  Hence **"$\xi_3^{\mathbf C}$ is the unique constant for which $R^{\mathbf C}_\xi$
  satisfies the $W_4$-functional equation"** is a correct reformulation, and it is
  equivalent to Conjecture D for the pair $(\mathbf C,\mathbf F)$ at $p=3$.

## 4. Numerical evidence  (`padic_point.gp`, `padic_point_hi.gp`)

Rational points with $(1-8x_0)$ a square, so that all four branch values are rational:

| $x_0$ | $u$ | $v$ | $\bar u$ | $\bar v$ | $G$ | $\bar G$ |
|---|---|---|---|---|---|---|
| $-1$ | $3/2$ | $1/4$ | $-3/5$ | $1/25$ | $-1/2$ | $2/5$ |
| $-3$ | $15/7$ | $9/49$ | $-5/4$ | $1/16$ | $-2/7$ | $1/4$ |
| $-6$ | $14/5$ | $4/25$ | $-21/11$ | $9/121$ | $-1/5$ | $2/11$ |

All five cover identities of §3 of `RIGIDITY_PROOF.md` hold exactly at these points,
and $|u|_3,|v|_3,|\bar u|_3,|\bar v|_3\le1<9$.

With $\xi=\xi_3^{\mathbf C}$ taken as $b_N/a_N$, $N=1800$ (accurate to $\approx3^{3500}$),
series truncated at $K=900$, working precision $3^{1600}$:

$$v_3\bigl(\Xi_{x_0}(\xi_3)\bigr)\ \ge\ 1599\quad\text{at all three points}$$
i.e. **zero to full working precision $3^{1600}$**. (At $M=700,K=650$: $v_3\ge699$.
`RIGIDITY_PROOF.md` reported $3^{44}$ coefficientwise in $x$; this is the same
statement evaluated pointwise, ~36× deeper.)

Conversely, solving the *truncated* equation for $\xi$ recovers $\xi_3^{\mathbf C}$:

| $K$ | 50 | 100 | 200 | 300 | 400 | 600 | 900 |
|---|---|---|---|---|---|---|---|
| $v_3\bigl(c_0^{(K)}/c_1^{(K)}-\xi_3\bigr)$ | 93 | 192 | 388 | 594 | 690\* | 1185 | 1594\* |

(\* capped by the working precision $3^{700}$ resp. $3^{1600}$.) The recovery rate is
$\approx1.98K$, i.e. the full slope $\sigma_3=2$: the $W_4$-functional equation at a
*single* rational point determines $\xi_3^{\mathbf C}$ at the maximal possible rate.
Same table at $x_0=-3,-6$. **No failure of any kind was observed.**

For reference, $v_3(\xi_3^{\mathbf C})=-1$ and $\xi_\infty^{\mathbf C}=\tfrac12L(\chi_{-3},2)=0.390651206448243148\ldots$
(this identification of the archimedean constant is itself confirmed to 37 digits by §2).

*Limitation.* No rational $x_0$ with $|x_0|_3=3$ and split cover exists
($1-8x_0=s^2$ forces $v_3(x_0)\in\{0,-2\}$ over $\mathbb Q$), so the tests all sit at
$|x_0|_3\le1$; probing $1<|x|_3<9$ needs quadratic (or ramified) points. That is the
obvious next numerical step and is cheap.

## 5. What Calegari/Coleman would provide, and what remains

*Not consulted directly* (no network access was used here); from the known structure:

* **Coleman**: $p$-adic Eichler integrals / the $p$-adic Eichler–Shimura map give a
  cocycle $\gamma\mapsto P_\gamma^{(p)}$ on $\Gamma_0(12)$ (and its AL-extension)
  valued in $\mathbb Q_p$-polynomials of degree $\le k-2=1$, whose coefficients are
  $p$-adic $L$-values of $\Phi_{\mathbf F}$ obtained from the *unit-root* splitting of
  the $U_p$-operator. What we need is exactly: the $W_4$-component of that cocycle,
  evaluated on the overconvergent avatar of $\Theta_{\mathbf C}$, is
  $\tfrac54\xi_3(j-1)$ with the second (quasi-period, $L(\cdot,1)$) coefficient
  *absent* — its absence is the 3-adic input with no archimedean analogue, and it is
  what our numerics are seeing.
* **Calegari 2005** (`Irrationality of certain p-adic periods for small p`, IMRN):
  the relevant mechanism is that an Apéry-type limit $\lim b_n/a_n$ is the constant
  term / value of an overconvergent modular form of small weight on the ordinary
  locus, and irrationality/identification follows from the $U_p$-eigenvalue
  decomposition. Reading the actual paper is the single highest-value next step and
  was not possible here.

**What would remain even with the above.** One still needs (i) the identification of
$\xi_3$ with a value of the $p$-adic Eichler integral — this is *not* proved
anywhere in the project; Proposition C only gives the limit's existence and slope;
(ii) the vanishing of the $L(\cdot,1)$-coefficient 3-adically; (iii) H1 ($\kappa_3=0$)
and H2 (Sturm certification of the cover identities), both routine.

## 6. Verdict (candid)

The reduction is now clean and the target is a single, sharply-stated,
*single-curve* fact: $\xi_3^{\mathbf C}$ is the unique $\xi\in\mathbb Q_3$ for which
$R^{\mathbf C}_\xi$ satisfies the $W_4$-functional equation on the double cover, and
this is equivalent to Conjecture D for $(\mathbf C,\mathbf F)$ at $p=3$. The
archimedean analogue is now an exact, verified proposition — but it is an analogue
with a *nonzero* period, so it does **not** transfer, and the previous claim that the
complex-analytic argument "proves H3 over $\mathbb C$" is withdrawn. The 3-adic
statement holds numerically to $3^{1600}$ at three rational points, which is strong
evidence and nothing more. **Conjecture D for (C,F) at $p=3$ is open.**
