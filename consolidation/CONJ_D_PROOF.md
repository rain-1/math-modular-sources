# Conjecture D at $p=3$ for Zagier's rows B, C, F, and the identification
# $\xi_3=\tfrac12\zeta_3(2)$

**Summary.** The $3$-adic Apéry limits of the three Zagier rows B, C, F are
*computed in closed form*:
$$\boxed{\ \xi_3^{\mathbf B}=\xi_3^{\mathbf C}=\tfrac12\zeta_3(2),\qquad
\xi_3^{\mathbf F}=\tfrac58\zeta_3(2)\ }$$
where $\zeta_3=L_3(\,\cdot\,,\mathbf 1)$ is the Kubota–Leopoldt $3$-adic zeta
function and $\zeta_3(2)$ is its value on the branch interpolating
$L(1-n,\chi_{-3})$ at odd $n$ (see §3). The rational factors $\tfrac12,\tfrac12,\tfrac58$
are *exactly the archimedean factors* $r_{\mathbf B},r_{\mathbf C},r_{\mathbf F}$
of $\Lambda=L(2,\chi_{-3})$. So Conjecture D holds for all three pairs at $p=3$,
with $\Lambda_3=\zeta_3(2)$ the $3$-adic avatar of $L(2,\chi_{-3})$.

**Logical status.**

* For **C** and **F** the result is *proved*, modulo (i) the two standard
  citations recorded in §4 (Coleman's Eisenstein family; Katz's theorem that the
  weight of a $p$-adic modular form is determined by its $q$-expansion) and
  (ii) the input **H1** ($v_3(a_n)=O(\log n)$), verified for $n\le400$.
  In particular the $3$-adic period identity **H3** of `PADIC_PERIOD.md`,
  previously listed as OPEN, is now proved (§6): it is the statement that the
  weight-$(-1)$ *overconvergent* form $\Theta_{\mathbf F}-\xi^{\mathbf F}_3$ is an
  Atkin–Lehner eigenvector for $W_4$ on $X_0(12)$, which is a two-line
  computation on the oldform basis $\{\mathcal E,V_2\mathcal E,V_4\mathcal E\}$
  once the Eichler integral has been replaced by the overconvergent Eisenstein
  series.
* For **B** everything except the *last* geometric step goes through verbatim;
  the gap is stated precisely in §8. Numerically $\xi_3^{\mathbf B}=\tfrac12\zeta_3(2)$
  holds to $3^{328}$ (full working precision).

The method is Calegari's \[Cal05\]; the novelty is only that the source
$\Phi_X=P_X(V_2)S$ is an *oldform combination* of an Eisenstein series of
weight $3$, so the Eichler integral $\Theta_X=D^{-2}\Phi_X$ is $P_X(V_2/4)$
applied to a *single* overconvergent Eisenstein series of weight $-1$, and the
Mellin polynomial $P_X$ evaluated at the $V_2$-eigenvalue $2^{-2}=\tfrac14$ is
what produces the rational factor.

---

## 0. Setting and notation

$q=e^{2\pi i\tau}$, $D=q\,d/dq$, $(V_df)(\tau)=f(d\tau)$, $p=3$,
$\chi=\chi_{-3}$ the quadratic character of conductor $3$. Note $\chi_{-3}=\omega$,
the Teichmüller character at $p=3$.

$$S=E_{3,\chi,1}=\sum_{m\ge1}c(m)q^m,\qquad c(m)=\sum_{d\mid m}\chi(m/d)d^2,
\qquad L(S,s)=L(\chi,s)\zeta(s-2).$$

| $X$ | recurrence $(a,b,c)$ | level $N_X$ | $t_X,F_X$ | $P_X(Y)$, $\Phi_X=P_X(V_2)S$ | $P_X(\tfrac14)$ |
|---|---|---|---|---|---|
| **C** | $(10,3,9)$ | $6$ | $t_{\mathbf C}=\frac{\eta_1^4\eta_6^8}{\eta_2^8\eta_3^4}$, $F_{\mathbf C}=\frac{\eta_2^6\eta_3}{\eta_1^3\eta_6^2}$ | $1-8Y$ | $-1$ |
| **F** | $(17,6,72)$ | $12$ | $t_{\mathbf F}=\frac{\eta_1^5\eta_3\eta_4^5\eta_6^2\eta_{12}}{\eta_2^{14}}$, $F_{\mathbf F}=\frac{\eta_2^{15}\eta_3^2\eta_{12}^2}{\eta_1^6\eta_4^6\eta_6^5}$ | $(1+Y)(1-8Y)$ | $-\tfrac54$ |
| **B** | $(9,3,27)$ | $36$ | $t_{\mathbf B}=\frac{\eta_1^3\eta_4^3\eta_{18}^9}{\eta_2^9\eta_9^3\eta_{36}^3}$, $F_{\mathbf B}=\frac{\eta_2^9\eta_3\eta_{12}}{\eta_1^3\eta_4^3\eta_6^3}$ | $1-6Y-8Y^2$ | $-1$ |

$F_X\in M_1(\Gamma_0(N_X),\chi)$; $\Theta_X:=D^{-2}\Phi_X$; the rows are
$A^X(t_X)=F_X=\sum a_nt_X^n$, $B^X(t_X)=F_X\Theta_X=\sum b_nt_X^n$;
$R_\xi^X:=B^X-\xi A^X$; $\xi_3^X=\lim b_n/a_n$ (Proposition C, $\sigma_3>0$).

$P_C(\tfrac14)=1-2=-1$; $P_F(\tfrac14)=\tfrac54\cdot(-1)=-\tfrac54$;
$P_B(\tfrac14)=1-\tfrac32-\tfrac12=-1$. **VERIFIED.**

---

## 1. Lemma 1 (overconvergence lemma; `RIGIDITY_PROOF.md` §2)

Let $a_n\in\mathbb Z_p$ with $v_p(a_n)=O(\log n)$ (**H1**) and slope $\sigma_p>0$. Then
$R_{\xi_p}$ converges on $|t|_p<p^{\sigma_p}$; $A$ has radius of convergence
exactly $1$ and diverges for every $|t|_p>1$; and **$\xi_p$ is the unique
$\xi\in\mathbb Q_p$ for which $R_\xi$ converges on some disc $|t|_p\le\rho$ with
$\rho>1$.**

*Status of* **H1**: $a_n\in\mathbb Z$ for all three rows, and
$\max_{n\le400}v_3(a_n)=10,13,11$ for B, C, F respectively (`radius.gp`, rerun
here). So $v_3(a_n)$ is *bounded* over the tested range, far stronger than
$O(\log n)$. A proof would follow from Katz's weight-rigidity as in §4(c); we do
not write it out.

For C: $c=9,\ \sigma_3=2$; F: $c=72,\ \sigma_3=2$; B: $c=27,\ \sigma_3=3$.

---

## 2. The Eichler integral is $P_X(V_2/4)$ applied to one Eisenstein series

Since $DV_d=dV_dD$, on the zero-constant-term sector $D^{-2}V_d=d^{-2}V_dD^{-2}$
(Prop. `prop:mellinshift`). Hence with
$$\mathcal E:=D^{-2}S=\sum_{n\ge1}\Bigl(\sum_{e\mid n}\chi(e)e^{-2}\Bigr)q^n
\qquad\Bigl(\text{because }\tfrac{c(n)}{n^2}=\sum_{d\mid n}\chi(n/d)(d/n)^2
=\sum_{e\mid n}\chi(e)e^{-2}\Bigr)$$
one has $\Theta_X=P_X(V_2/4)\,\mathcal E$, explicitly

$$\Theta_{\mathbf C}=(1-2V_2)\mathcal E,\qquad
\Theta_{\mathbf F}=(1-\tfrac74V_2-\tfrac12V_4)\mathcal E,\qquad
\Theta_{\mathbf B}=(1-\tfrac32V_2-\tfrac12V_4)\mathcal E .$$

**VERIFIED exactly** to $O(q^{120})$, together with
$\Phi_{\mathbf C}=(1-8V_2)S$, $\Phi_{\mathbf F}=(1-7V_2-8V_4)S$,
$\Phi_{\mathbf B}=(1-6V_2-8V_4)S$ and the divisor-sum identity for $\mathcal E$.

Because $\chi(3)=0$, the divisor sum defining $\mathcal E$ is automatically
$3$-*depleted*: $\sum_{e\mid n}\chi(e)e^{-2}=\sum_{e\mid n,\,3\nmid e}\chi(e)e^{-2}$.
This is exactly the shape of Calegari's $\sigma^*_\kappa$.

---

## 3. $\mathcal E$ is the $q$-part of an overconvergent Eisenstein series

Let $\mathcal W=\operatorname{Hom}_{\rm cont}(\mathbb Z_3^\times,\mathbb C_3^\times)$
be weight space and $\mathcal W^+=\{\kappa:\kappa(-1)=1\}$. Write $a=\omega(a)\langle a\rangle$.

**Definition.** $\kappa_0(a):=\omega(a)a^{-1}=\langle a\rangle^{-1}$.

$\kappa_0(-1)=\langle-1\rangle^{-1}=1$, so $\kappa_0\in\mathcal W^+$; in the
parametrisation of $\mathcal W^+$ by classical weights $a\mapsto a^{2k}=\langle a\rangle^{2k}$
(note $\omega^{2k}=1$ for $p=3$), $\kappa_0$ is the point "$2k=-1$", i.e. **weight
$-1$ with nebentypus $\omega=\chi_{-3}$**. With Calegari's
$\sigma^*_\kappa(n)=\sum_{d\mid n,\,(d,3)=1}\kappa(d)d^{-1}$ one gets
$$\sigma^*_{\kappa_0}(n)=\sum_{d\mid n,\,3\nmid d}\omega(d)d^{-2}
=\sum_{e\mid n}\chi(e)e^{-2},$$
i.e. $\mathcal E$ is precisely the $q$-part of $E_{\kappa_0}$.

> **Theorem (Coleman; \[Cal05, Thm. 2.3\]).** There is a function $\zeta_3(\kappa)$ on
> $\mathcal W^+$, rigid analytic away from $\kappa=1$, such that
> $E_\kappa=\tfrac12\zeta_3(\kappa)+\sum_{n\ge1}\sigma^*_\kappa(n)q^n$ varies rigid
> analytically over $\mathcal W^+$ and specialises at every $\kappa\in\mathcal W^+$ to an
> **overconvergent** eigenform of weight $\kappa$ and level $\Gamma_0(3)$; at
> $\kappa:a\mapsto a^{2k}$, $2k\ge2$, it is the classical $E^*_{2k}$.

*(References: R. Coleman, "$p$-adic Banach spaces and families of modular forms",
Invent. Math. **127** (1997) 417–479; Coleman–Mazur, "The eigencurve", LMS Lect.
Notes **254** (1998). Calegari's Thm. 2.3 is the statement in the form we need.
Note that **no character of $p$-power conductor has to be invoked**: because
$\chi_{-3}=\omega$, the weight-character $\kappa_0=\langle\cdot\rangle^{-1}$ lies on the
*same* connected component of $\mathcal W^+$ as the trivial character, so the
classical Eisenstein family with trivial nebentypus already passes through it.)*

**Corollary (constant term).** Put
$$\kappa:=\tfrac12\zeta_3(\kappa_0)=\tfrac12\zeta_3(2),\qquad\text{so}\qquad
\mathcal E+\kappa=E_{\kappa_0}\ \ \text{is overconvergent of weight }\kappa_0 .$$
By \[Cal05, Lem. 2.4\] (continuity of $\zeta_3$ on $\mathcal W^+$),
$\zeta_3(\kappa_0)=\lim\zeta(1-2k)$ over even $2k\to\infty$ with $2k\to-1$ in
$\mathbb Z_3$; equivalently, in Kubota–Leopoldt normalisation
(Washington, *Introduction to Cyclotomic Fields*, Thm. 5.11, valid for **all**
$n\ge1$),
$$\zeta_3(1-n)=L_3(1-n,\mathbf 1)=
\begin{cases}-(1-3^{n-1})B_n/n,& n\ \text{even},\\[2pt]
-B_{n,\omega}/n=L(1-n,\chi_{-3}),& n\ \text{odd},\end{cases}$$
so that $\zeta_3(2)$ is the value at $n=-1$ of the odd branch, i.e. the
$3$-adic interpolation of $L(1-n,\chi_{-3})$ — the honest $3$-adic avatar of
$L(2,\chi_{-3})$. *(This is **not** the "$\zeta_p(2)=0$" of \[Cal05, §4\]: that
one is $L_3(2,\omega)$, attached to an odd character, hence identically zero.
The nonzero object is $L_3(2,\mathbf 1)$.)*

Both branches were implemented and **cross-validated exactly** (see §7):
$v_3(\kappa)=-1$ and
$$\kappa\equiv 3^{-1}\cdot 8386265965554334030 \pmod{3^{39}} .$$

---

## 4. Overconvergence of $H_\xi$

Set $H_\xi:=F_X\,(\Theta_X-\xi)$, a weight-$0$ object of level $N_X$.

**(a) $V_d$ preserves overconvergence for $p\nmid d$.** The degeneracy map
$\tau\mapsto d\tau$ is induced by a finite flat map of rigid spaces
$X_1(Nd)\to X_1(N)$ which commutes with the canonical-subgroup section and
therefore carries (strict neighbourhoods of) the ordinary component containing
$\infty$ to itself. (Coleman, *loc. cit.*, §B; Coleman–Mazur §2. This is the
step that fails for $p\mid d$, cf. §9.) Hence
$P_X(V_2/4)E_{\kappa_0}=\Theta_X+P_X(\tfrac14)\kappa$
is overconvergent of weight $\kappa_0$ and level $\operatorname{lcm}(3,4\cdot3)=12$
(resp. $6$ for **C**, where only $V_1,V_2$ occur).

**(b) Product with a classical form.** $F_X$ is classical of weight $1$, hence
overconvergent; the product of an overconvergent form of weight $1$ (nebentypus
$\chi$) and one of weight $\kappa_0$ (nebentypus $\omega=\chi$) is overconvergent of
weight $1\cdot\kappa_0=$ trivial and trivial nebentypus, i.e. a **rigid analytic
function on a strict neighbourhood of the ordinary component containing $\infty$**.
Therefore
$$\xi=\xi^*_X:=-P_X(\tfrac14)\,\kappa
\quad\Longrightarrow\quad H_{\xi^*_X}\ \text{is overconvergent of weight }0 .$$

**(c) Conversely,** for $\xi\ne\xi^*_X$ the difference
$H_\xi-H_{\xi^*_X}=(\xi^*_X-\xi)F_X$ would then also be overconvergent of weight
$0$; but $F_X\ne0$ is classical of weight $1$, and by Katz's theorem that the
weight of a $p$-adic modular form is determined by its $q$-expansion (N. Katz,
*p-adic properties of modular schemes and modular forms*, LNM 350, Cor. 4.4.2 —
the reference Calegari cites as "\[8\], 4.4"), $1\ne0$ in $\mathcal W$ forces a
contradiction. So

$$\boxed{\ H_\xi\ \text{is overconvergent}\iff \xi=\xi^*_X=-P_X(\tfrac14)\,\kappa. \ }$$

Numerically: $\xi^*_{\mathbf C}=\xi^*_{\mathbf B}=\kappa=\tfrac12\zeta_3(2)$ and
$\xi^*_{\mathbf F}=\tfrac54\kappa=\tfrac58\zeta_3(2)$.

*Sign check against the archimedean side.* The archimedean statement is
$\xi_\infty^X=L(\Phi_X,2)=P_X(2)\,L(\chi,2)\zeta(0)$ with $\zeta(0)=-\tfrac12$,
giving $\tfrac12,\tfrac12,\tfrac58$ times $L(2,\chi)$. Since $P_X(2)$ means
$V_2\mapsto2^{-2}=\tfrac14$, this is literally $-P_X(\tfrac14)\cdot\tfrac12L(2,\chi)$:
*the same formula*, with $\tfrac12L(2,\chi)$ (the constant term of the classical
weight-$(-1)$ Eisenstein series) replaced by $\kappa=\tfrac12\zeta_3(2)$.
The minus sign is not decorative: it comes from $H_\xi=F_X(\Theta_X-\xi)$, so
$-\xi$ plays the role of the constant term.

---

## 5. The radius jump (Calegari §3.1)

**Lemma 2 (ordinary disc).** For each $X$, every point of $X_0(N_X)^{\rm an}$ with
$0<|t_X|_3<1$ satisfies $|j|_3>1$; consequently the region $\{|t_X|_3<1\}$ lies in
the residue disc of the cusp $\infty$, is connected, and is contained in the
connected component $D_\infty$ of the ordinary locus containing $\infty$.

*Proof.* $j=q^{-1}+744+\dots\in q^{-1}\mathbb Z[[q]]$ and $t_X=q+O(q^2)\in q\mathbb Z[[q]]$,
and the cusp $\infty$ of $\Gamma_0(N)$ has width $1$; hence $t_X\,j\in\mathbb Z[[t_X]]$
with constant term $1$. **VERIFIED** for B, C, F (coefficients
$1,741,196884,\dots$ / $1,740,196878,\dots$ / $1,739,196878,\dots$). For
$0<|t_X|<1$ all higher terms have $|\cdot|<1$, so $|t_Xj|=1$ and $|j|=|t_X|^{-1}>1$.
A point with $|j|>1$ has potentially multiplicative reduction, hence is ordinary.
$\square$

Since $D_\infty$ is a *closed* affinoid containing the open disc $\{|t_X|<1\}$
(in the branch of $t_X^{-1}$ through $\infty$), it contains $\{|t_X|\le1\}$, and a
strict neighbourhood of it contains $\{|t_X|\le\rho\}$ for some $\rho>1$.

**Lemma 3 (descent to the $t_X$-line).** Suppose $t_X$ is a hauptmodul for a
genus-zero group $\Gamma_X\supseteq\Gamma_0(N_X)$ and $H_{\xi^*_X}$ is
$\Gamma_X$-invariant. Then the $t_X$-expansion of $H_{\xi^*_X}$ converges on
$|t_X|_3\le\rho$ for some $\rho>1$.

*Proof.* $H_{\xi^*}$ is a rigid function on a strict neighbourhood $U\supseteq D_\infty$
and descends to $X_{\Gamma_X}\cong\mathbb P^1_{t_X}$, on which the image of $U$
contains $\{|t_X|\le\rho\}$, $\rho>1$, by Lemma 2. On a disc a rigid function is a
convergent power series. $\square$

**Theorem.** Under the hypotheses of Lemma 3 (and **H1**),
$$\xi_3^X=\xi^*_X=-P_X(\tfrac14)\cdot\tfrac12\zeta_3(2).$$

*Proof.* $H_{\xi^*_X}$ has $t_X$-expansion $R^X_{\xi^*_X}=B^X-\xi^*_XA^X$
(the same formal $q\to t_X$ substitution). By Lemma 3 it converges on a disc of
radius $>1$. By Lemma 1(3) that forces $\xi^*_X=\xi_3^X$. $\square$

---

## 6. Verification of the hypothesis of Lemma 3 for C and F

**C.** $t_{\mathbf C}$ has degree $1$ on $X_0(6)$ (Ligozat: order $+1$ at the cusp
$\infty$, $-1$ at the cusp $d=2$, $0$ elsewhere), i.e. it is a hauptmodul of
$\Gamma_0(6)$ itself. $\Theta_{\mathbf C}-\xi^*=(1-2V_2)E_{\kappa_0}$ has level $6$,
$F_{\mathbf C}\in M_1(\Gamma_0(6),\chi)$, so $H_{\xi^*}$ is a $\Gamma_0(6)$-invariant
weight-$0$ form. Take $\Gamma_{\mathbf C}=\Gamma_0(6)$. $\checkmark$

**F.** $t_{\mathbf F}$ has degree $2$ on $X_0(12)$ (Ligozat: $+1$ at $\infty$, $+1$ at
the cusp $d=3$, $-2$ at the cusp $d=2$), and it is a hauptmodul of
$\Gamma_{\mathbf F}=\langle\Gamma_0(12),W_4\rangle$ — this is the statement
$t_{\mathbf F}(\gamma\tau)=t_{\mathbf F}(\tau)$ for
$\gamma=\bigl(\begin{smallmatrix}8&-3\\12&-4\end{smallmatrix}\bigr)=W_4$
verified to $10^{-37}$ in `PADIC_PERIOD.md` §1. So we must show
$H_{\xi^*}|W_4=H_{\xi^*}$.

*Action of $W_4$ on the oldform basis.* Let $g$ have level $3$, weight $k$,
nebentypus $\chi$, and let $W_4=\bigl(\begin{smallmatrix}4&1\\12&4\end{smallmatrix}\bigr)$.
Writing $\delta_e=\bigl(\begin{smallmatrix}e&0\\0&1\end{smallmatrix}\bigr)$ and
$V_eg=e^{-k/2}g|_k\delta_e$, one computes
$\delta_1W_4=\bigl(\begin{smallmatrix}1&1\\3&4\end{smallmatrix}\bigr)\delta_4$,
$\delta_4W_4=4\bigl(\begin{smallmatrix}4&1\\3&1\end{smallmatrix}\bigr)$,
$\delta_2W_4=2\bigl(\begin{smallmatrix}2&1\\3&2\end{smallmatrix}\bigr)\delta_2$, all
inner matrices in $\Gamma_0(3)$, whence
$$V_1\mapsto 4^{k/2}V_4,\qquad V_2\mapsto\chi(2)V_2=-V_2,\qquad V_4\mapsto4^{-k/2}V_1 .$$
(The other choice of representative, e.g. $\gamma$ above, differs by an element of
$\Gamma_0(12)$ with $\chi(d)=-1$ and hence flips the global sign; this is the usual
$\pm$ ambiguity of $W_Q$ in the presence of a nebentypus, and it cancels in the
product below.)

* Weight $k=3$: $1-7V_2-8V_4\mapsto 8V_4+7V_2-1=-(1-7V_2-8V_4)$, i.e.
  $\Phi_{\mathbf F}|W_4=-\Phi_{\mathbf F}$ for this representative
  ($=+\Phi_{\mathbf F}$ for $\gamma$, matching the numerical check in `PADIC_PERIOD.md`).
* Weight $k=-1$: $V_1\mapsto\tfrac12V_4$, $V_2\mapsto-V_2$, $V_4\mapsto2V_1$, so
  $$1-\tfrac74V_2-\tfrac12V_4\ \longmapsto\ \tfrac12V_4+\tfrac74V_2-1
  =-\bigl(1-\tfrac74V_2-\tfrac12V_4\bigr).$$

So $\Theta_{\mathbf F}-\xi^*_{\mathbf F}=P_{\mathbf F}(V_2/4)E_{\kappa_0}$ is a
$W_4$-eigenform of weight $-1$ with the **same** eigenvalue as $\Phi_{\mathbf F}$
has in weight $3$; and $F_{\mathbf F}$ is a $W_4$-eigenform of weight $1$ with that
same eigenvalue (verified to $10^{-37}$, `PADIC_PERIOD.md` §1:
$F_{\mathbf F}(\gamma\tau)=(6\tau-2)F_{\mathbf F}(\tau)$, $\Phi_{\mathbf F}(\gamma\tau)=(6\tau-2)^3\Phi_{\mathbf F}(\tau)$).
Hence $H_{\xi^*}=F_{\mathbf F}\cdot(\Theta_{\mathbf F}-\xi^*)$ has $W_4$-eigenvalue
$(\pm1)^2=+1$: it is $\Gamma_{\mathbf F}$-invariant. $\checkmark$

**This is exactly the input H3 of `RIGIDITY_PROOF.md`/`PADIC_PERIOD.md`.** Over
$\mathbb C$ it is *false* for the Eichler integral $\Theta_{\mathbf F}$ — the period
polynomial $\Theta_{\mathbf F}|_{-1}\gamma-\Theta_{\mathbf F}=\tfrac54\xi_\infty(6\tau-3)
-i\tfrac{\pi^2}{6\sqrt3}(6\tau-1)$ is nonzero (`PADIC_PERIOD.md` §2). Three-adically
the obstruction disappears because $\Theta_{\mathbf F}-\xi^*_{\mathbf F}$ is not a formal
primitive at all: it is a *genuine* modular form of weight $-1$, on which $W_4$ acts
by the oldform formula. The single constant $\xi^*_{\mathbf F}$ is precisely what turns
the non-modular primitive into a modular object. This also explains the
archimedean/$3$-adic asymmetry noted in `PADIC_PERIOD.md` §6.

---

## 7. Numerical verification

Scripts run for this note (PARI/GP; kept in the scratch log, values reproducible
from the recipes below).

**(a) Independent computation of $\kappa$.** Washington's explicit formula
(*Introduction to Cyclotomic Fields*, Thm. 5.11)
$$L_3(s,\mathbf 1)=\frac1{3(s-1)}\sum_{a\in\{1,2\}}\langle a\rangle^{1-s}
\sum_{j\ge0}\binom{1-s}{j}\Bigl(\frac3a\Bigr)^{j}B_j ,$$
implemented to $3$-adic precision $3^{320}$. It was validated against the
interpolation property **exactly** (difference $=0$ as a rational number) at
$n=2,4,6$ (even branch, $-(1-3^{n-1})B_n/n$) and $n=3,5,7,9$ (odd branch,
$-B_{n,\omega}/n$ with $B_{n,\omega}=3^{n-1}\sum_{a=1}^{3}\chi(a)B_n(a/3)$).

**(b) Main identification.** With $\kappa=\tfrac12L_3(2,\mathbf 1)$ and rows to
$N=300$ (so $b_N/a_N$ is accurate to $\approx3^{600}$):

| | $v_3(\xi^X_3-\xi^*_X)$ |
|---|---|
| $\xi^{\mathbf C}_3-\kappa$ | $328$ |
| $\xi^{\mathbf B}_3-\kappa$ | $328$ |
| $\xi^{\mathbf F}_3-\tfrac54\kappa$ | $328$ |

$328$ is the working precision of the $L_3$ series (truncation $K=320$), i.e. the
agreement is **exact to full precision**; it is not a measurement of a discrepancy.
$v_3(\kappa)=-1$ and
$\kappa\equiv3^{-1}\cdot8386265965554334030\pmod{3^{39}}$, matching the value of
$\xi^{\mathbf C}_3$ recorded in `RIGIDITY_PROOF.md` §5 digit for digit.

**(c) Kummer-congruence cross-check (no $p$-adic $L$-function code).**
$\kappa_m:=-B_{n,\chi}/(2n)$ at $n=2\cdot3^m-1$ (exact rationals, `bernpol`):

| $m$ | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| $n$ | 5 | 17 | 53 | 161 | 485 | 1457 |
| $v_3(\kappa_m-\xi^{\mathbf C}_3)$ | 0 | 1 | 2 | 3 | 4 | 5 |

i.e. $\kappa_m\to\xi_3^{\mathbf C}$ at exactly the Kummer rate $3^{m-1}$. This is an
*independent* confirmation of $\xi_3^{\mathbf C}=\lim L(1-n,\chi_{-3})/2$.

**(d) Structural checks** (all **EXACT** to $O(q^{120})$): the three source
identities $\Phi_X=P_X(V_2)S$; $\mathcal E$'s divisor sum; the three identities
$\Theta_X=P_X(V_2/4)\mathcal E$; $t_X\,j\in\mathbb Z[[t_X]]$ with constant term $1$
for $X=$ B, C, F. Also $P_{\mathbf B}(\tfrac14)=P_{\mathbf C}(\tfrac14)=-1$,
$P_{\mathbf F}(\tfrac14)=-\tfrac54$.

**(e) Cusp divisors** (Ligozat, $\operatorname{ord}_{a/d}=\frac{N}{24\gcd(d^2,N)}\sum_\delta\frac{\gcd(d,\delta)^2r_\delta}{\delta}$):
$\deg t_{\mathbf C}=1$ on $X_0(6)$; $\deg t_{\mathbf F}=2$ on $X_0(12)$;
$\deg t_{\mathbf B}=6$ on $X_0(36)$ (zeros $+1$ at $\infty$, $+1$ at $d=9$, $+4$ at
$d=18$).

---

## 8. What remains for row B

Everything in §§1–5 applies verbatim to **B**:
$\Theta_{\mathbf B}-\xi^*_{\mathbf B}=(1-\tfrac32V_2-\tfrac12V_4)E_{\kappa_0}$ is
overconvergent of weight $-1$ and level $12$, and it is even a $W_4$-eigenform by
the same computation as in §6 (the condition is $c_4=-\tfrac12c_1$, satisfied by
both $\mathbf F$ and $\mathbf B$); $F_{\mathbf B}\in M_1(\Gamma_0(36),\chi)$; so
$H_{\xi^*_{\mathbf B}}$ is an overconvergent modular *function* on $X_0(36)$, and
Lemma 2 gives $\{|t_{\mathbf B}|_3<1\}\subseteq D_\infty$.

**The gap** is Lemma 3's hypothesis. $X_0(36)$ has genus $1$ and
$\deg t_{\mathbf B}=6$; thus $\mathbb Q(t_{\mathbf B})$ has index $6$ in the function
field, and one must identify the corresponding genus-zero group
$\Gamma_{\mathbf B}\supseteq\Gamma_0(36)$ (index $6$) and check that
$H_{\xi^*_{\mathbf B}}$ is $\Gamma_{\mathbf B}$-invariant. $W_4$ accounts for a factor
$2$; the remaining factor $3$ is not identified here. (Note this is the *same*
degree-$3$ phenomenon flagged in `RIGIDITY_PROOF.md` §7: exactly,
$(9x-1)t_{\mathbf B}(27t_{\mathbf B}^2-9t_{\mathbf B}+1)+x(8x-1)^2=0$ with
$x=t_{\mathbf F}$, re-verified here to $O(q^{80})$, so
$[\mathbb Q(t_{\mathbf B},t_{\mathbf F}):\mathbb Q(t_{\mathbf F})]=3$.)

Given the present proof, however, the B case no longer needs a new *arithmetic*
idea — only the group-theoretic identification above. And
$\xi^{\mathbf B}_3=\tfrac12\zeta_3(2)$ holds numerically to $3^{328}$.

---

## 9. The general statement, and what fails at $p\mid d$

**Theorem R$''$.** Let $(\Gamma,t,F)$ be a modular Apéry system of weight
$w$, level $N$, with source $\Phi=P(V_{d_1},\dots)E$, where

1. $E$ is an Eisenstein series of weight $w+2$ whose $p$-adic avatar
   $E+\kappa_E$ is the specialisation of the Coleman Eisenstein family at a point
   $\kappa_E\in\mathcal W^+$ (equivalently: the nebentypus of $E$, twisted by
   $\omega^{-(w+2)}$, has conductor prime to $p$ — for $p=3$, $\chi_{-3}=\omega$
   makes this automatic);
2. every $d_i$ occurring in $P$ is prime to $p$;
3. $t$ is a hauptmodul of a genus-zero $\Gamma'\supseteq\Gamma$, $H_{\xi^*}$ is
   $\Gamma'$-invariant, and $t\,j\in\mathbb Z_{(p)}[[t]]$ with unit constant term;
4. $v_p(a_n)=O(\log n)$.

Then $\xi_p=-P(p^{-(w+1)}\text{-eigenvalues})\cdot\kappa_E$, i.e.
$\xi_p=P(\text{Mellin at }s=w+1)\cdot\Lambda_p$ for one constant
$\Lambda_p=-\kappa_E$ independent of $P$. In particular Conjecture D holds
throughout the prime-to-$p$ Hecke orbit of a single Eisenstein source, and
$\xi_p=0$ whenever the Mellin polynomial vanishes at $s=w+1$ (the $\Psi_0$
phenomenon of `RIGIDITY_PROOF.md` §6).

**$p\mid d$.** Hypothesis (2) is essential: at $p\mid d$ the map $\tau\mapsto d\tau$
is the Frobenius/Verschiebung direction, does not commute with the canonical
subgroup, and does not preserve strict neighbourhoods of the ordinary locus. In
the $\zeta(3)$ pair at $p=2$ the connecting correspondence is
$(1-V_2)(1-4V_2)(1-16V_2)$ *at $p=2$*, with
$P_\alpha(3)/P_\varepsilon(3)=(-\tfrac7{12})/(-\tfrac7{16})=\tfrac43$, matching the
archimedean ratio and the measured $v_2(\Delta_n)=4n+O(1)$. Calegari's own §4
(Catalan at $p=2$, level $\Gamma_1(4)=\Gamma_1(p^2)$) is **not** a $p\mid d$ case:
there the *level* is divisible by $p$ but the *source* is a single Eisenstein
series with no $V_p$ in it, exactly as our C, F, B at $p=3$ have level divisible
by $3$ and only $V_2,V_4$ in the source. So Calegari's §3/§4 machinery applies to
our situation but not to the $\zeta(3)$ pair; there one must replace "both
branches" by the canonical-subgroup section and $P_C$ by its unit-root factor.
That remains open.

---

## 10. Summary of what is proved and what is cited

* **Proved here (new):** $\xi_3^{\mathbf C}=\tfrac12\zeta_3(2)$ and
  $\xi_3^{\mathbf F}=\tfrac58\zeta_3(2)$; hence
  $\xi_3^{\mathbf F}=\tfrac54\xi_3^{\mathbf C}$ (Conjecture D for $(\mathbf C,\mathbf F)$
  at $p=3$, previously only *reduced* to the open input H3), and the
  identification of the $3$-adic period as a Kubota–Leopoldt $L$-value.
  H3 itself is proved (§6).
* **Cited:** Coleman's Eisenstein family and its overconvergence
  \[Coleman 1997; Coleman–Mazur 1998; = Cal05, Thm. 2.3, Lem. 2.4\];
  Katz's weight rigidity \[Katz LNM 350, 4.4\]; $V_d$ preserves overconvergence
  for $p\nmid d$ \[Coleman 1997 §B\]; Washington Thm. 5.11 for the interpolation
  formula; Ligozat for eta-quotient divisors.
* **Assumed (measured):** **H1**, $v_3(a_n)=O(\log n)$ (bounded by $13$ for
  $n\le400$).
* **Open:** row **B**'s group-theoretic step (§8); the $p\mid d$ case (§9).

---

**Update (see `THEOREM_F_HYPOTHESES.md` §3).**  The row-$\mathbf B$ gap is
closed.  $\Gamma_{\mathbf B}=\gamma\,\Gamma_0(9)\,\gamma^{-1}$ with
$\gamma=\bigl(\begin{smallmatrix}1&1/2\\0&1\end{smallmatrix}\bigr)$
(index $72/12=6$ over $\Gamma_0(36)$, $=\deg t_{\mathbf B}$), equivalently:
row $\mathbf B$ has a second modular presentation at level $9$,
$t=\eta_9^3/\eta_1^3$ (a hauptmodul of $\Gamma_0(9)$),
$F=\eta_1^3/\eta_3$, source $E^{\chi_{-3},\mathbf 1}_3$ *primitive*, in which
$a_n\mapsto(-1)^na_n$, $b_n\mapsto-(-1)^nb_n$, so
$\xi^{\mathbf B}_3=-(-\tfrac12\zeta_3(2))=\tfrac12\zeta_3(2)$.  Hypothesis
(b) is then vacuous ($\deg t=1$).  The $W_4$-computation of §6 is the special
case $Q=4$ of a general lemma: the Mellin shift $V_d\mapsto d^{-(w+1)}V_d$
preserves Atkin–Lehner eigenvalues while moving the weight from $w+2$ to $-w$.
