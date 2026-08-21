# The Euler-factor criterion for $p$-adic Apéry limits

*Companion to `CONJ_D_PROOF.md`, which is the special case $\psi(p)=0$ (Euler
factor $\equiv1$) at $p=3$.  Scripts: `lattice/euler_criterion/`
(`lp.gp`, `rows.gp`, `criterion.gp`, `verify.gp`); logs `criterion.log`,
`euler3000.log` (run on **snake**, `tools/SNAKE.md`).*

---

## 0. Summary

Let $(\Gamma,t,F)$ be a modular Apéry row of weight $w$ whose companion source
is the Eisenstein series $\Phi=P(V)E$.  Write $E$'s two characters as
$(\psi,\varphi)$ (definitions in §1) and let

$$\mathcal E_p(s)\;=\;1-\psi(p)\,p^{-s}$$

be the $p$-Euler factor of the *first* $L$-factor of $L(\Phi,s)$.  Then:

> **Criterion.**  The row has a $p$-adic Apéry limit exactly when
> $\mathcal E_p(s)$ divides the Mellin polynomial $P(s)=\sum_d c_dd^{-s}$ in
> $\mathbb Q[p^{-s}]$ (a vacuous condition when $p\mid\operatorname{cond}\psi$).
>
> **Value.**  With $Q=P/\mathcal E_p$,
> $$\boxed{\ \xi_p=-\,Q(w+1)\,\kappa_p,\qquad
> \kappa_p=\begin{cases}\tfrac12L_p\bigl(w+1,\ \psi\omega^{-w}\bigr)&\varphi=\mathbf 1,\\[2pt]
> 0&\varphi\neq\mathbf 1,\end{cases}}$$
> $L_p$ the Kubota–Leopoldt $p$-adic $L$-function.  For a **cuspidal** source
> with $U_p$-eigenvalue $a_p$ the Euler factor is $1-a_pp^{-s}$, $\kappa_p=0$,
> and hence $\xi_p=0$ whenever it exists.
>
> **Consequence (Conjecture D).**  If $\xi_\infty=r_\infty\Lambda$ with
> $\Lambda=L(\psi,w+1)$, then $\xi_p=r_p\Lambda_p$ with
> $\Lambda_p=L_p(w+1,\psi\omega^{-w})$ and
> $$\boxed{\ r_p=\frac{r_\infty}{\mathcal E_p(w+1)}\ }$$
> The correction factor depends only on $(\psi,p,w)$, **not on the row**.
> Hence any two rows sharing the Eisenstein series $E$ have
> $r_p/r'_p=r_\infty/r'_\infty$: Conjecture D.

**Status.** The criterion and the value formula are *verified exactly* on the
whole census — every row with a positive slope, at every slope prime — to
$\ge p^{2991}$ (§4).  The criterion reproduces the slope set
$\{p:\sigma_p>0\}=\{p:p\mid c\}$ for all twelve Eisenstein families and all
$p\in\{2,3,5\}$, with no exception (§3).  The derivation (§2) is Calegari's
overconvergence argument of `CONJ_D_PROOF.md` with one new input — $V_p$
preserves overconvergence — and the same two standard citations; hypothesis
(b) (descent to $\mathbb P^1_t$) has to be re-checked at level $Np$ in the four
rows where $Q$ still contains $V_p$, and is recorded there as **Open**.

---

## 1. Notation: the two placements

Write $q=e^{2\pi i\tau}$, $D=q\,d/dq$, $(V_df)(\tau)=f(d\tau)$, $w=r-1$,
so the source has weight $w+2$.  A holomorphic Eisenstein series of weight
$k=w+2$ with characters $(\psi,\varphi)$ is

$$E^{\psi,\varphi}_{k}=\delta(\psi)\,L(1-k,\varphi)+\sum_{n\ge1}
\Bigl(\sum_{d\mid n}\psi(n/d)\varphi(d)\,d^{\,w+1}\Bigr)q^n,
\qquad L\bigl(E^{\psi,\varphi},s\bigr)=L(\psi,s)L(\varphi,s-w-1),$$

$\delta(\psi)=1$ iff $\psi=\mathbf 1$ (Diamond–Shurman, Thm. 4.5.1).  In the
notation of `paper/sections/02_sources.tex`, Table 1:

* **inner placement** $\sigma^{(w+1)}_\chi(m)=\sum_{d\mid m}\chi(m/d)d^{w+1}$
  gives $\psi=\chi$, $\varphi=\mathbf 1$;
* **outer placement** $\widetilde\sigma^{(w+1)}_\chi(m)=\sum_{d\mid m}\chi(d)d^{w+1}$
  gives $\psi=\mathbf 1$, $\varphi=\chi$.

Thus $\psi$ is the character of the **first** $L$-factor $L(\psi,s)$ and
$\varphi$ that of the **shifted** factor $L(\varphi,s-w-1)$.  Applying
$D^{-(w+1)}$ moves the exponent onto the *other* divisor variable:

$$\mathcal E:=D^{-(w+1)}E=\sum_{n\ge1}\Bigl(\sum_{e\mid n}\psi(e)\,e^{-(w+1)}\varphi(n/e)\Bigr)q^n,
\qquad L(\mathcal E,s)=L(\varphi,s)L(\psi,s+w+1).$$

So in $\mathcal E$ the character $\psi$ sits on the variable carrying the
negative power, and $\varphi$ on the co-divisor.  This is the "$\chi'$" of the
task statement: $\chi'=\psi$ ($=\chi$ inner, $=\mathbf 1$ outer).

**$p$-depletion.**  Deleting the terms with $p\mid e$ multiplies
$L(\psi,s+w+1)$ by its inverse $p$-Euler factor, i.e.

$$\mathcal E^{[p]}:=\bigl(1-\psi(p)p^{-(w+1)}V_p\bigr)\,\mathcal E ,\qquad
\mathcal E=\bigl(1-\psi(p)p^{-(w+1)}V_p\bigr)^{-1}\mathcal E^{[p]} .$$

**Mellin polynomial.**  $\Phi=P(V)E$ with $P(V)=\sum_dc_dV_d$; write
$P(s)=\sum_dc_dd^{-s}$, so that $L(\Phi,s)=P(s)L(\psi,s)L(\varphi,s-w-1)$ and,
since $D^{-(w+1)}V_d=d^{-(w+1)}V_dD^{-(w+1)}$ on $q\mathbb Q[[q]]$,

$$\Theta:=D^{-(w+1)}\Phi=P\bigl(V_\bullet/\bullet^{\,w+1}\bigr)\mathcal E,
\qquad P\bigl(V_\bullet/\bullet^{\,w+1}\bigr)\,c=P(w+1)\,c\ \ (c\text{ constant}).$$

---

## 2. The theorem

> **Theorem E (Euler-factor criterion).**  Let $(\Gamma,t,F)$ be a modular
> Apéry row of weight $w$ and level $N$ with Eisenstein source
> $\Phi=P(V)E^{\psi,\varphi}_{w+2}$, $A=F$, $B=F\Theta$, $\Theta=D^{-(w+1)}\Phi$.
> Fix a prime $p$ and assume
>
> **(a)** *(Euler-factor divisibility)* $\mathcal E_p(s)=1-\psi(p)p^{-s}$
> divides $P(s)$ in $\mathbb Q[p^{-s}]$; put $Q:=P/\mathcal E_p$;
>
> **(b)** *(descent)* $t$ is a Hauptmodul of a genus-zero group
> $\Gamma_t\supseteq\Gamma_0(N')$, $N'=\operatorname{lcm}(N,p\cdot\text{levels in }Q)$,
> and $H_{\xi^*}:=F\,(\Theta-\xi^*)$ is $\Gamma_t$-invariant;
>
> **(c)** *(ordinary disc)* $t\,j\in\mathbb Z_{(p)}[[t]]$ with unit constant
> term (so $\{|t|_p\le1\}\subseteq D_\infty$);
>
> **(d)** *(H1)* $v_p(a_n)=O(\log n)$ and the slope $\sigma_p>0$.
>
> Then $b_n/a_n$ converges in $\mathbb Q_p$ and
> $$\xi_p=-\,Q(w+1)\,\kappa_p,\qquad
> \kappa_p=\delta(\varphi)\cdot\tfrac12L_p\bigl(w+1,\psi\omega^{-w}\bigr).$$
> Conversely, if (a) fails then $\Theta$ is not a finite oldform combination of
> overconvergent forms and no $\xi$ makes $H_\xi$ overconvergent.
>
> **Cuspidal clause.**  If instead $\Phi=P(V)f$ with $f$ a normalised cuspidal
> eigenform of weight $w+2$ and $U_pf=a_pf$, the Euler factor is
> $1-a_pp^{-s}$ and, when it divides $P(s)$, $\Theta=Q(V_\bullet/\bullet^{w+1})\,
> D^{-(w+1)}f^{[p]}$ is overconvergent **with no constant term**; hence
> $\xi_p=0$.

### Proof

**Step 1 (the depleted Eisenstein series is overconvergent).**
$\mathcal E^{[p]}$ has $q$-coefficients
$\sum_{e\mid n,\,p\nmid e}\psi(e)e^{-(w+1)}\varphi(n/e)$.  Decompose
$\psi=\psi^{(p)}\psi_p$ into prime-to-$p$ and $p$-power parts.  Then on
$e\in\mathbb Z_p^\times$,
$$\psi(e)e^{-(w+1)}=\psi^{(p)}(e)\cdot\bigl[\psi_p(e)e^{-w}\bigr]\cdot e^{-1},$$
so $\mathcal E^{[p]}$ is the $q$-part of the Coleman–Mazur Eisenstein family
with tame characters $(\varphi,\psi^{(p)})$ at the weight-character
$$\kappa_0=\psi_p\,\omega^{-w}\langle\cdot\rangle^{-w}\in\mathcal W,$$
which is the point "classical weight $-w$".  By Coleman's theorem
[Coleman 1997; Coleman–Mazur 1998; = Cal05, Thm. 2.3] the family is rigid
analytic over $\mathcal W^+$ and every specialisation is an **overconvergent**
eigenform of that weight.  For $p=3$, $\psi=\chi_{-3}=\omega$, $w=1$ this is
exactly $\kappa_0=\langle\cdot\rangle^{-1}$ of `CONJ_D_PROOF.md` §3.

**Step 2 (the constant term).**  The constant term of the family is the
constant term of its classical specialisations, which for
$E^{\varphi,\psi}_{k}$ (co-divisor character $\varphi$) is
$\delta(\varphi)L(1-k,\psi)$.  Hence:

* if $\varphi\neq\mathbf 1$, **every** classical member has constant term $0$,
  so the family — and $G^{(p)}:=\mathcal E^{[p]}+\kappa_p$ — has
  $\kappa_p=0$;
* if $\varphi=\mathbf 1$, the $p$-stabilised constant term at weight $k$ is
  $-\tfrac12(1-\psi(p)p^{k-1})B_{k,\psi}/k=\tfrac12L_p(1-k,\psi\omega^{k})$
  (Washington, Thm. 5.11, valid for all $k\ge1$), and by continuity of $L_p$
  on $\mathcal W^+$ the value at $k=-w$ is
  $$\kappa_p=\tfrac12L_p\bigl(w+1,\psi\omega^{-w}\bigr).$$
  *Parity is automatic:* $E^{\psi,\varphi}_{w+2}$ requires
  $\psi\varphi(-1)=(-1)^{w}$, so $\varphi=\mathbf 1$ forces
  $\psi(-1)=(-1)^w$, i.e. $\psi\omega^{-w}$ is **even** and $L_p(\cdot,\psi\omega^{-w})
  \not\equiv0$.  Conversely on the census $\varphi\neq\mathbf 1$ always has
  $\varphi$ odd, so the two vanishing mechanisms ($\delta(\varphi)=0$ and
  "$L_p\equiv0$ for odd character") agree; we take $\delta(\varphi)$ as the
  definition since it is the one that is a statement about the family.

**Step 3 (the criterion).**  By §1,
$\Theta=P(V_\bullet/\bullet^{w+1})\mathcal E$ and
$(1-\psi(p)p^{-(w+1)}V_p)$ is exactly the operator
$\mathcal E_p(s)\big|_{p^{-s}\mapsto V_p/p^{w+1}}$.  So (a) holds iff
$$\Theta=Q\bigl(V_\bullet/\bullet^{w+1}\bigr)\,\mathcal E^{[p]}
=Q\bigl(V_\bullet/\bullet^{w+1}\bigr)\bigl(G^{(p)}-\kappa_p\bigr)
=Q\bigl(V_\bullet/\bullet^{w+1}\bigr)G^{(p)}\;-\;Q(w+1)\,\kappa_p ,$$
a **finite** oldform combination of an overconvergent form.  If (a) fails,
$P/\mathcal E_p$ is an infinite power series $\sum_k\psi(p)^kp^{-k(w+1)}V_p^k$
whose partial sums leave every fixed strict neighbourhood of the ordinary
locus (each $V_p$ shrinks the radius of overconvergence by $p$), so no finite
$\xi$ makes $H_\xi$ overconvergent.

**Step 4 ($V_p$ is harmless).**  This is the point at which the present
statement goes beyond `CONJ_D_PROOF.md` §9, whose hypothesis (2) demanded
$p\nmid d$ for every $V_d$ in $P$.  On overconvergent forms, $V_p$ is the
Frobenius: for $0<v<p/(p+1)$ it maps $M^\dagger_k(v)\to M^\dagger_k(v/p)$
(Coleman, *Classical and overconvergent modular forms*, Invent. Math. **124**
(1996), §§2–3; Coleman, Invent. Math. **127** (1997), §B).  It therefore
**preserves overconvergence**, merely with a smaller strict neighbourhood;
the level is multiplied by $p$.  Calegari's own $E^*_{2k}$ at $p=2$ already
contains $V_2$.  Since Lemma 2 of `CONJ_D_PROOF.md` only requires *some*
$\rho>1$, the shrinkage is immaterial.  Consequently the four census rows in
which $Q$ retains $V_p$ (**E**, **F** at $p=2$; $\alpha$, $\varepsilon$ at
$p=2$) are covered.

**Step 5 (conclusion, as in `CONJ_D_PROOF.md` §§4–5).**  $F$ is classical of
weight $w$ and nebentypus $\chi_F$; the product
$H_{\xi^*}=F\cdot\bigl(\Theta+Q(w+1)\kappa_p\bigr)$ is overconvergent of weight
$w+(-w)=0$ and trivial nebentypus, i.e. a rigid function on a strict
neighbourhood of $D_\infty$.  By (c) that neighbourhood contains
$\{|t|_p\le\rho\}$, $\rho>1$, and by (b) $H_{\xi^*}$ descends to
$\mathbb P^1_t$, so its $t$-expansion $B-\xi^*A$ converges there.  By the
overconvergence lemma (`RIGIDITY_PROOF.md` §2, restated as Lemma 1 of
`CONJ_D_PROOF.md`) and (d), $\xi_p$ is the **unique** constant with that
property, whence $\xi_p=\xi^*=-Q(w+1)\kappa_p$.  Uniqueness in the other
direction uses Katz's weight rigidity (LNM 350, Cor. 4.4.2): $H_\xi-H_{\xi^*}
=(\xi^*-\xi)F$ would be overconvergent of weight $0$ and of weight $w\neq0$.
$\square$

**Cuspidal clause.**  By Coleman's $\theta^{k-1}$-surjectivity onto
$U_p$-depleted overconvergent forms (Invent. Math. **124** (1996)),
$D^{-(k-1)}f^{[p]}$ is overconvergent of weight $2-k$ for
$f^{[p]}=f-a_pV_pf$; it has no constant term because a cusp form has none and
$\theta^{k-1}$ is a bijection on $q$-expansions without constant term.  Hence
$H_\xi=F(\Theta-\xi)$ is overconvergent iff $\xi=0$.  The divisibility
condition is $\bigl(1-a_pY\bigr)\mid P(Y)$, $Y=V_p/p^{k-1}$; it is automatic
when $p^2\mid N$ (then $a_p=0$).

### Corollary (Conjecture D, and the Euler correction of the rational factor)

Suppose $\varphi=\mathbf 1$.  Then
$$\xi_\infty=L(\Phi,w+1)=P(w+1)\,L(\psi,w+1)\,\zeta(0)=-\tfrac12P(w+1)\,\Lambda,
\qquad \Lambda:=L(\psi,w+1),$$
$$\xi_p=-\tfrac12\,Q(w+1)\,\Lambda_p,\qquad
\Lambda_p:=L_p\bigl(w+1,\psi\omega^{-w}\bigr).$$
Hence $r_\infty=-\tfrac12P(w+1)$, $r_p=-\tfrac12Q(w+1)$ and
$r_p=r_\infty/\mathcal E_p(w+1)$: the *same* correction for every row with the
same $E$ and the same $p$.  Two such rows therefore satisfy
$r\,\xi'_p=r'\,\xi_p$ — **Conjecture D** — and $\Lambda_p$ is the $p$-adic
avatar of $\Lambda$.  When $\varphi\neq\mathbf 1$, $\xi_p=0$: the archimedean
period survives but its $p$-adic avatar is $0$.

---

## 3. The criterion reproduces the slope set

`lattice/euler_criterion/criterion.gp`, exact, on the twelve-family source
table of `paper/sections/02_sources.tex`.  For each family and each
$p\in\{2,3,5\}$ it tests $\mathcal E_p\mid P$ and compares with $p\mid c$.

| fam | $w$ | $P(w+1)$ | $c$ | $p=2$ | $p=3$ | $p=5$ |
|---|---|---|---|---|---|---|
| **A** | 1 | $3/4$ | $-8$ | **YES**, $Q=1$ / $2\mid c$ | no / no | no / no |
| **B** | 1 | $-1$ | $27$ | no / no | **YES**, $Q=-1$ / $3\mid c$ | no / no |
| **C** | 1 | $-1$ | $9$ | no / no | **YES**, $Q=-1$ / $3\mid c$ | no / no |
| **D** | 1 | $1$ | $-1$ | no / no | no / no | no / no |
| **E** | 1 | $-1$ | $32$ | **YES**, $Q=-1$ / $2\mid c$ | no / no | no / no |
| **F** | 1 | $-5/4$ | $72$ | **YES**, $Q=-1$ / $2\mid c$ | **YES**, $Q=-5/4$ / $3\mid c$ | no / no |
| $\alpha$ (Domb) | 2 | $-7/12$ | $64$ | **YES**, $Q=-2/3$ / $2\mid c$ | no / no | no / no |
| $\gamma$ (Apéry) | 2 | $-1/3$ | $1$ | no / no | no / no | no / no |
| $\delta$ | 2 | $-13/27$ | $81$ | no / no | **YES**, $Q=-1/2$ / $3\mid c$ | no / no |
| $\varepsilon$ | 2 | $-7/16$ | $16$ | **YES**, $Q=-1/2$ / $2\mid c$ | no / no | no / no |
| $\zeta$ | 2 | $1$ | $-27$ | no / no | **YES**, $Q=1$ / $3\mid c$ | no / no |
| $\eta$ | 2 | $-1$ | $125$ | no / no | no / no | **YES**, $Q=-1$ / $5\mid c$ |

**Every entry agrees**: $\mathcal E_p\mid P\iff p\mid c$, $36$ tests, no
exception.  This is a new, purely algebraic characterisation of the slope set
of Proposition C ($\sigma_p=v_p(c)+2\kappa_p$) for integral rows.

Explicit factorisations (all exact, `criterion.gp`), with $X=2^{-s}$, $Y=3^{-s}$:
$$P_{\mathbf A}=1-X,\quad P_{\mathbf F}=(1+X)(1-8X),\quad
P_\alpha=(1-X)(1-16X)(1-9Y),$$
$$P_\varepsilon=(1-X)(1-4X)(1-16X),\quad P_\delta=(1-Y)(1-14X+16X^2),$$
and the non-divisibilities (each is $P$ evaluated at the root of
$\mathcal E_p$, the remaining variable being set to its $s=w+1$ value):
$P_{\mathbf C}|_{X=-1}=9$, $P_{\mathbf B}|_{X=-1}=-1$, $P_\eta|_{X=-1}=-1$,
$P_\gamma|_{X=1}=-27(1-Y)|_{Y=1/27}=-26$,
$P_\gamma|_{Y=1}=64(1-X)|_{X=1/8}=56$,
$P_\alpha|_{Y=1}=-8(1-X)(1-16X)|_{X=1/8}=7$,
$P_\delta|_{X=1}=3(1-Y)|_{Y=1/27}=26/9$.
These are exactly the checks in Fable's census list, **all confirmed**.

---

## 4. Census verification (exact, snake, `euler3000.log`)

Method.  (i) $\kappa_p$ from Washington's Thm. 5.11 series with $s\in\mathbb Z$
(then $\binom{1-s}{j}$ is an exact integer and no precision is lost),
implemented to $p$-adic precision $p^{3060}$ and **validated exactly** against
the interpolation property $L_p(1-n,\chi)=-(1-\chi\omega^{-n}(p)p^{n-1})
B_{n,\chi\omega^{-n}}/n$ for $n=1,\dots,9$ and all four characters used
($\mathbf 1$ at $p=2,3,5$ and $\chi_{12}$ at $p=2$): difference $=O(p^{\text{prec}})$
in every case.  (ii) Rows generated exactly over $\mathbb Q$ to $N$ with
$\sigma_pN\ge3040$; $\xi_N=b_N/a_N$.  (iii) Report
$v_p(\xi_N-\xi_{N-1})$ (the Cauchy precision actually achieved) and
$v_p(\xi_N-\text{prediction})$.

Independent cross-check of every $\kappa_p$ by the **Kummer sequence** of exact
generalised Bernoulli numbers $\kappa^{(m)}=-(1-\psi(p)p^{k-1})B_{k,\psi}/(2k)$
at $k\equiv-w \pmod{(p-1)p^m}$ (resp. $\bmod\ 2^{m+1}$ for $p=2$):
$v_p(\kappa^{(m)}-\kappa_p)$ increases by exactly $1$ per $m$ in all six cases
— the Kummer rate — confirming $\kappa_p$ independently of the $L_p$ code.

### 4.1 Main table

$\Lambda_p=L_p(w+1,\psi\omega^{-w})$; $\kappa_p=\tfrac12\Lambda_p$ or $0$;
$\xi_p=-Q(w+1)\kappa_p$.  "digits" $=v_p(\xi_N-\xi_p^{\rm pred})$, to be
compared with the Cauchy precision of the row at that $N$ (next column);
equality means **agreement to full available precision**, not a discrepancy.

| row | $p$ | $(\psi,\varphi)$ | $\mathcal E_p(w{+}1)$ | $Q(w{+}1)$ | $\kappa_p$ | $\xi_p$ | $N$ | digits | Cauchy |
|---|---|---|---|---|---|---|---|---|---|
| **A** $(7,2,-8)$ | 2 | $(\mathbf 1,\chi_{-3})$ | $3/4$ | $1$ | $0$ | $\mathbf 0$ | 1014 | **3003** | 3003 |
| **E** $(12,4,32)$ | 2 | $(\chi_{-4},\mathbf 1)$ | $1$ | $-1$ | $\tfrac12\zeta_2(2)$ | $\tfrac12\zeta_2(2)$ | 608 | **3026** | 3005 |
| **F** $(17,6,72)$ | 2 | $(\chi_{-3},\mathbf 1)$ | $5/4$ | $-1$ | $\tfrac12L_2(2,\chi_{12})$ | $\tfrac12L_2(2,\chi_{12})$ | 1014 | **3005** | 3005 |
| **F** $(17,6,72)$ | 3 | $(\chi_{-3},\mathbf 1)$ | $1$ | $-5/4$ | $\tfrac12\zeta_3(2)$ | $\tfrac58\zeta_3(2)$ | 1520 | **3025** | 3027 |
| **B** $(9,3,27)$ | 3 | $(\chi_{-3},\mathbf 1)$ | $1$ | $-1$ | $\tfrac12\zeta_3(2)$ | $\tfrac12\zeta_3(2)$ | 1014 | **3029** | 3024 |
| **C** $(10,3,9)$ | 3 | $(\chi_{-3},\mathbf 1)$ | $1$ | $-1$ | $\tfrac12\zeta_3(2)$ | $\tfrac12\zeta_3(2)$ | 1520 | **3025** | 3027 |
| $\alpha$ Domb $(10,4,64)$ | 2 | $(\mathbf 1,\mathbf 1)$ | $7/8$ | $-2/3$ | $\tfrac12\zeta_2(3)$ | $\tfrac13\zeta_2(3)$ | 507 | **2991** | 2992 |
| $\varepsilon$ $(12,4,16)$ | 2 | $(\mathbf 1,\mathbf 1)$ | $7/8$ | $-1/2$ | $\tfrac12\zeta_2(3)$ | $\tfrac14\zeta_2(3)$ | 760 | **3000** | 2986 |
| $\delta$ $(7,3,81)$ | 3 | $(\mathbf 1,\mathbf 1)$ | $26/27$ | $-1/2$ | $\tfrac12\zeta_3(3)$ | $\tfrac14\zeta_3(3)$ | 760 | **3030** | 3029 |
| $\zeta$ $(9,3,-27)$ | 3 | $(\chi_{-3},\chi_{-3})$ | $1$ | $1$ | $0$ | $\mathbf 0$ | 1014 | **3023** | 3019 |
| $\eta$ $(11,5,125)$ | 5 | $(\chi_5,\mathbf 1)$ | $1$ | $-1$ | $\tfrac12\zeta_5(3)$ | $\tfrac12\zeta_5(3)$ | 1014 | **3027** | 3025 |
| Cooper $s_{18}$ | 3 | (conj. $(\chi_{-3},\mathbf 1)$) | $1$ | $-1$ | $\tfrac12\zeta_3(2)$ | $\tfrac12\zeta_3(2)$ | 3040 | **3025** | 3026 |
| Zudilin Catalan | 2 | (conj. $(\chi_{-4},\mathbf 1)$) | $1$ | $-2$ | $\tfrac12\zeta_2(2)$ | $\zeta_2(2)$ | 380 | **3015** | 3003 |
| cusp $L(f,2)$, $L12$ | 2 | cuspidal, $a_2=0$ | $1$ | $1$ | $0$ | $\mathbf 0$ | 1520 | **3027** | 3016 |

Notation: $\zeta_p(s)=L_p(s,\mathbf 1)$.  In Washington's normalisation
$\omega$ is the character mod $4$ at $p=2$, so $\zeta_2(2)=L_2(2,\mathbf 1)$
**is** the number Calegari writes $L_2(2,\chi_{-4})$: both are
$\lim_{n\to-1,\,n\ \rm odd}\,-B_{n,\chi_{-4}}/n$.  Likewise
$\zeta_3(2)=L_3(2,\mathbf 1)=\lim_{n\to-1,\,n\ \rm odd}-B_{n,\chi_{-3}}/n$,
and $3\kappa_3\equiv8386265965554334030\pmod{3^{40}}$ — digit for digit the
value recorded in `CONJ_D_PROOF.md` §7.

### 4.2 Reading of the table

1. **$\alpha$ (Domb) and $\varepsilon$ at $p=2$** are both rational multiples
   of $\zeta_2(3)/2$ — Calegari's number — with $\xi_2^\alpha/\xi_2^\varepsilon
   =\tfrac13/\tfrac14=\tfrac43$, exactly the archimedean ratio
   $(7/24)/(7/32)$.  Conjecture D at $p=2$ for the $\zeta(3)$ pair, previously
   the open $p\mid d$ case of `CONJ_D_PROOF.md` §9, is now a computed identity
   to $2^{2991}$, and the $4/3$ is explained: both rational factors are divided
   by the *same* $\mathcal E_2(3)=7/8$.
2. **E at $p=2$** gives $\xi_2^{\mathbf E}=\tfrac12\zeta_2(2)$, i.e. exactly
   Calegari's $L_2(2,\chi_{-4})/2$ (§4 of \[Cal05\]) — the normalisation
   matches on the nose.
3. **A at $p=2$** (outer placement, $\xi_\infty=\zeta(2)/4$) has
   $\xi_2^{\mathbf A}=0$: $P_{\mathbf A}$ *is* the Euler factor, $Q=1$, and the
   co-divisor character $\varphi=\chi_{-3}$ is non-trivial, so the Eisenstein
   family through this point has identically zero constant term.  The row has
   $\sigma_2=3>0$ and $b_n/a_n\to0$ at that rate: $v_2(b_N/a_N)=3003$ at
   $N=1014$.  **The archimedean period $\zeta(2)$ has $2$-adic avatar $0$ here.**
   The same happens for $\zeta$ at $p=3$ ($\varphi=\chi_{-3}$): $\xi_3^\zeta=0$
   to $3^{3023}$, even though $L(\Phi_\zeta,3)=L(\chi_{-3},3)L(\chi_{-3},0)\ne0$.
   The second vanishing mechanism agrees: $\Lambda_2^{\mathbf A}=L_2(2,\mathbf 1\cdot\omega^{-1})
   =L_2(2,\omega)$ and $\Lambda_3^\zeta=L_3(3,\chi_{-3}\omega^{-2})=L_3(3,\chi_{-3})$
   are $L_p$'s of **odd** characters, hence identically zero.
4. **F at $p=2$** is the first genuine $p\mid d$ case handled: $\psi(2)=
   \chi_{-3}(2)=-1$, $\mathcal E_2=1+X$, $P_{\mathbf F}=(1+X)(1-8X)$, so the
   Euler factor is *cancelled* and $Q=1-8V_2$ still contains $V_2=V_p$.  The
   value $\xi_2^{\mathbf F}=\tfrac12L_2(2,\chi_{12})$ confirms both the
   cancellation and Step 4 ($V_p$ preserves overconvergence).  Note
   $\chi_{12}=\chi_{-3}\chi_{-4}=\chi_{-3}\omega$ is the even character of
   conductor $12$; the twist by $\omega^{-w}=\omega$ is essential and visible.
5. **$\eta$ at $p=5$**: $\chi_5$ here is the **quadratic** (Legendre) character
   mod $5$, even, matching the $L(\chi_5,s)$ of Table 1 (the quartic characters
   $\psi_1,\psi_2$ belong to row **D**).  Since $\omega^{-2}=\omega^2=\chi_5$ at
   $p=5$, $\psi\omega^{-w}=\chi_5^2=\mathbf 1$ and $\Lambda_5=\zeta_5(3)$.
6. **$\zeta$ family**, $L(\chi_{-3},s)L(\chi_{-3},s-3)$: a product of two
   $L$-functions with the *same* non-trivial character; the $3$-adic avatar is
   $0$ by item 3.  This is the cleanest instance of the $\varphi\neq\mathbf 1$
   clause.
7. **Zudilin's Catalan row** ($c=1$, $\kappa_2=4$ non-integral, $\sigma_2=8$)
   is *not* in the twelve-family table, but the value formula still applies
   with $Q(2)=-2$: $\xi_2=\zeta_2(2)=2\,\xi_2^{\mathbf E}$, matching the
   archimedean ratio $G:\tfrac12G$.  **This is a new Conjecture-D pair**
   (Zudilin, **E**) at $p=2$ with $\Lambda_2=\zeta_2(2)$, and it shows the
   criterion/value formula survives non-integrality ($\kappa_p>0$).
8. **Cooper $s_{18}$**: $\xi_3=\tfrac12\zeta_3(2)$ to $3^{3025}$, so the fourth
   row on the $L(2,\chi_{-3})$ line at $p=3$ is confirmed with the *same*
   $\Lambda_3$ as **B**, **C**, **F** — upgrading the "aligned" entry of
   `SLOPE_CENSUS.md` §1 to an identified value.  (Its Eisenstein source is
   still unidentified, Remark 2.x of the paper; the fit $Q(2)=-1$ is therefore
   an observation, not a derivation.)
9. **The cuspidal rows.**  The level-12 weight-3 cusp row
   ($f=\eta_2^3\eta_6^3$, $a_2=0$, recurrence
   $(n{+}1)^2a_{n+1}=(20n^2{+}10n{+}2)a_n-16(2n{-}1)^2a_{n-1}$, $a_1=2$) has
   $\xi_2=0$ to $2^{3027}$ at $N=1520$ (slope $2$, cf. Remark on non-Zagier
   normalisations) and no $3$-adic slope ($a_3=-3$, $P=1$ not divisible by
   $1+3Y$) — exactly the cuspidal clause.

### 4.3 Rows/primes where the criterion fails: no convergence

Increments $v_p\bigl(\tfrac{b_n}{a_n}-\tfrac{b_{n-1}}{a_{n-1}}\bigr)$ at
$n=120,\dots,150$ (step 6), all bounded and negative — no $p$-adic limit:

| row | $p$ | increments |
|---|---|---|
| **C** | 2 | $-6,-2,-4,-2,-8,-2$ |
| **B** | 2 | $-6,-2,-4,-2,-8,-2$ |
| $\gamma$ (Apéry) | 2 | $-9,-3,-6,-3,-12,-3$ |
| $\gamma$ (Apéry) | 3 | $-3,-6,-3,-3,-6,-3$ |
| $\delta$ | 2 | $-9,-3,-6,-3,-12,-3$ |
| $\alpha$ Domb | 3 | $-3,-6,-3,-3,-6,-3$ |
| $\eta$ | 2 | $-9,-3,-6,-3,-12,-3$ |
| cusp $L(f,2)$ | 3 | $-1,-2,-1,-1,-2,-1$ |
| Domb$^*$ cuspidal | 2 | $-15,-15,-8,-6,-11,-11$ |
| Domb$^*$ cuspidal | 3 | $-2,-1,3,0,-2,2$ |

The last two are the **Domb cuspidal apparatus** of
`packages/phase 0 math campaign/domb/main6_fully_completed_v2.tex`
($A_n$ = Domb, $B^*$ from the inhomogeneous recurrence
$m^3B^*_m=(2m{-}1)(10(m{-}1)^2{+}10(m{-}1){+}4)B^*_{m-1}-64(m{-}1)^3B^*_{m-2}
+\binom{2m-2}{m-1}$, $B^*_n/A_n\to L(f_6,3)/2$).  Source
$f_*=f_6-4V_2f_6$, $f_6=(\eta_1\eta_2\eta_3\eta_6)^2$, $a_2=-2$, $a_3=-3$:
$P=1-4V_2=1-32Y$ is divisible by neither $1+2Y$ nor (at $p=3$) $1+3Y$, so the
criterion predicts **no $p$-adic slope at $2$ or $3$** — confirmed above.
(The row is inhomogeneous, so Proposition C does not apply and the prediction
is a genuine test of the criterion rather than of $v_p(c)$.)

### 4.4 The P-strand (task 3): towers do not rescue the failed cases

For **C** at $p=2$ and Domb at $p=3$ we tested convergence of
$\xi_{ap^s}=b_{ap^s}/a_{ap^s}$ along $p$-power towers, $a\in\{1,3,5\}$ resp.
$\{1,2,4\}$, $n\le600$.  The increments $v_p(\xi_{ap^s}-\xi_{ap^{s-1}})$ are
**strictly decreasing**, e.g. $-2,-4,-6,\dots,-18$ (C, $p=2$, $a=1$) and
$-3,-6,-9,-12,-15$ (Domb, $p=3$, $a=1$).  So there is no tower limit either:
when the Euler factor fails to divide, the geometric series
$\sum_k\psi(p)^kp^{-k(w+1)}V_p^k$ diverges in every sense we can measure.
**Negative result, recorded.**

---

## 5. What is proved, cited, assumed, open

* **Proved here (modulo the citations):** the criterion, the value formula
  $\xi_p=-Q(w+1)\kappa_p$ with $\kappa_p=\delta(\varphi)\tfrac12L_p(w+1,\psi\omega^{-w})$,
  the cuspidal clause $\xi_p=0$, and Conjecture D in the strong form
  $r_p=r_\infty/\mathcal E_p(w+1)$ for any two rows sharing an Eisenstein source.
* **Cited:** Coleman's Eisenstein family and overconvergence
  [Coleman, Invent. Math. **127** (1997); Coleman–Mazur 1998; Cal05 Thm. 2.3];
  $V_p$ preserves overconvergence, $M^\dagger_k(v)\to M^\dagger_k(v/p)$
  [Coleman, Invent. Math. **124** (1996)]; $\theta^{k-1}$ surjectivity onto
  $U_p$-depleted forms [ibid.]; Katz weight rigidity [LNM 350, 4.4.2];
  Washington Thm. 5.11; Diamond–Shurman Thm. 4.5.1 for the constant term.
* **Assumed (measured):** H1, $v_p(a_n)=O(\log n)$.
* **Open:** hypothesis (b) at level $Np$ for the four rows in which $Q$ retains
  $V_p$ (**E**, **F** at $p=2$; $\alpha$, $\varepsilon$ at $p=2$) — the
  Atkin–Lehner/Hauptmodul check of `CONJ_D_PROOF.md` §6 has to be redone one
  level higher; also row **B**'s group-theoretic step (unchanged from
  `CONJ_D_PROOF.md` §8), and the identification of the Eisenstein sources of
  $s_{18}$ and of Zudilin's row.  In all of these the value is confirmed
  numerically to $\ge p^{2991}$, so a failure could only be a failure of (b),
  not of the formula.

---

**Update (see `THEOREM_F_HYPOTHESES.md`).**  The three items under **Open**
above are now settled or reduced: $N'=N$ always (the "level $Np$" caveat is
vacuous); hypothesis (b) is a theorem for eleven of the fourteen rows, via an
Atkin–Lehner/Mellin-shift lemma plus, for rows $\mathbf B$, $\eta$, $\delta$,
a level-lowering $q\mapsto-q$ twin in which $t$ becomes a hauptmodul of
$\Gamma_0(9)$, $\Gamma_0(5)$, $\Gamma_0(6)$ — which in particular identifies
row $\mathbf B$'s index-six group and completes Conjecture D for $\mathbf B$;
and (d) is reduced to $v_p(a_n)=o(n)$, needed only for *existence* of the limit,
the value being forced by (a)–(c) alone whenever a limit exists.
