> **Review note (Fable).** The transfer identity of §4 is linear in ξ and holds for every ξ (the A- and B-parts transfer separately), hence H3 ⇔ "R^F_{(5/4)ξ} has radius > 1" ⇔ the theorem. Theorem R is therefore a *reduction*, not a proof: Conjecture D for (C,F) at p=3 is equivalent to the 3-adic period identity for Θ_C under the Atkin–Lehner involution of X_0(12) (the 3-adic analogue of the fold lemma, which over ℂ supplies the w-invariance through the period polynomial). Everything else in this document stands.

# Hecke functoriality of $p$-adic Apéry limits, and Conjecture D for Zagier's C, F at $p=3$

**Status of this document.** Problem 7.1 of `paper/main.pdf` (Conjecture D of
`paper/sections/04_padic.tex`) asks for a proof of $p$-adic rigidity for the Zagier
rows B, C, F at $p=3$. This document

* **proves** (modulo two explicitly flagged finite/measured inputs) that
  $\xi_3^{\mathbf F} = \tfrac54\,\xi_3^{\mathbf C}$, i.e. Conjecture D for the pair
  $(\mathbf C,\mathbf F)$ at $p=3$ — this is the first case of Conjecture D
  established by an argument rather than a computation;
* isolates the general mechanism (**Theorem R**, prime-to-$p$ Hecke functoriality of
  $p$-adic Apéry limits) and shows it is *not* a formal consequence of the source
  identity: it needs an affinoid overconvergence input which is exactly what fails
  when $p\mid d$;
* **shows that the pair $(\mathbf B,\mathbf C)$ is not a prime-to-$3$ case at all**,
  even though $\Phi_{\mathbf B}-\Phi_{\mathbf C}$ is built from $V_1,V_2,V_4$ only:
  Zagier's B coordinate $t_{\mathbf B}$ satisfies an irreducible **cubic** over
  $\mathbb Q(t_{\mathbf F})$, so $\mathbf B$ is a *degree-3 reparametrisation* of a
  level-12 source. Conjecture D for $(\mathbf B,\mathbf C)$ therefore splits into a
  prime-to-3 part (covered by Theorem R) and a genuinely $p$-adic
  *coordinate-invariance* gap, of the same nature as the $p\mid d$ obstruction at
  $p=2$ for Domb vs. $\varepsilon$.

All computations are in `lattice/rigidity/`; every identity below marked **EXACT**
was checked in exact rational arithmetic.

---

## 0. Notation

$q=e^{2\pi i\tau}$, $D=\theta_q=q\,d/dq$, $(V_df)(\tau)=f(d\tau)$. A *modular Apéry
system* is a triple $(\Gamma,t,F)$ with $t=q+O(q^2)$ a modular function, $F=1+O(q)$
a form of weight $w$; the rows are
$$A(t)=\sum a_nt^n := F,\qquad B(t)=\sum b_nt^n := F\cdot\Theta,\quad \Theta=D^{-(w+1)}\Phi,\quad \Phi=F\,Dt$$
(Theorem A and Theorem 2.4 of `ProjAnchors`, recalled in `02_sources.tex`).
For a row with $\kappa_p=0$ and slope $\sigma_p>0$, $\xi_p=\lim_n b_n/a_n\in\mathbb Q_p$
(Proposition C).

Throughout $S=E_{3,\chi_{-3},1}$, $c_S(m)=\sum_{d\mid m}\chi_{-3}(m/d)d^2$,
$L(S,s)=L(\chi_{-3},s)\zeta(s-2)$.

---

## 1. The source picture (all EXACT to $O(q^{220})$)

Script: `lattice/rigidity/setup.gp`, `cover_checks.gp`.

With the Zagier/eta parametrisations of Table `tab:eta`
($t_{\mathbf C}=\eta_1^4\eta_6^8/(\eta_2^8\eta_3^4)$, $F_{\mathbf C}=\eta_2^6\eta_3/(\eta_1^3\eta_6^2)$;
$t_{\mathbf F}=\eta_1^5\eta_3\eta_4^5\eta_6^2\eta_{12}/\eta_2^{14}$,
$F_{\mathbf F}=\eta_2^{15}\eta_3^2\eta_{12}^2/(\eta_1^6\eta_4^6\eta_6^5)$;
$t_{\mathbf B},F_{\mathbf B}$ at level 36):

$$\Phi_{\mathbf C}=(1-8V_2)S,\qquad \Phi_{\mathbf F}=(1-7V_2-8V_4)S=(1+V_2)\Phi_{\mathbf C},$$
$$\Phi_{\mathbf B}=(1-6V_2-8V_4)S=\Phi_{\mathbf C}+2V_2\Psi_0,\qquad \Psi_0:=(1-4V_2)S .$$

Mellin factors at $s=w+1=2$:
$P_{\mathbf C}(2)=P_{\mathbf B}(2)=-1$, $P_{\mathbf F}(2)=-\tfrac54$, $P_{2V_2(1-4V_2)}(2)=0$;
with $\zeta(0)=-\tfrac12$ this reproduces the archimedean $r$-values
$\tfrac12,\tfrac12,\tfrac58$ of the verification record.

Since $D^{-2}V_d=d^{-2}V_dD^{-2}$ on the zero-constant-term sector
(Proposition `prop:mellinshift`), writing $\Theta_X=D^{-2}\Phi_X$:
$$\boxed{\ \Theta_{\mathbf F}=\Theta_{\mathbf C}+\tfrac14 V_2\Theta_{\mathbf C}\ }\qquad\textbf{EXACT}$$
$$\Theta_{\mathbf B}=\Theta_{\mathbf C}+\tfrac12 V_2 D^{-2}\Psi_0 .\qquad\textbf{EXACT}$$

---

## 2. Overconvergence lemma

**Lemma 1 (overconvergence and uniqueness).**
Let $A(t)=\sum a_nt^n$, $B(t)=\sum b_nt^n$ be the rows of a system with
$a_n\in\mathbb Z_p$, denominator rate $\kappa_p=0$ (i.e. $v_p(a_n)=O(\log n)$)
and slope $\sigma_p>0$. Then:

1. $R_\xi:=B-\xi A$ converges on $|t|_p<p^{\sigma_p}$ for $\xi=\xi_p$;
2. $A$ has $p$-adic radius of convergence exactly $1$ and $\sum a_nt^n$ diverges for every $|t|_p>1$;
3. $\xi_p$ is the **unique** $\xi\in\mathbb Q_p$ for which $R_\xi$ converges on some disc $|t|_p\le\rho$ with $\rho>1$.

*Proof.* (1) By Proposition C, $v_p(b_n/a_n-\xi_p)=\sigma_pn+O(\log n)$, hence
$v_p(b_n-\xi_pa_n)=v_p(a_n)+\sigma_pn+O(\log n)=\sigma_pn+O(\log n)$; so
$\limsup|b_n-\xi_pa_n|^{1/n}=p^{-\sigma_p}$.
(2) $v_p(a_n)=O(\log n)$ gives $|a_n|_p^{1/n}\to1$; and for $|t|_p=p^{\varepsilon}$, $\varepsilon>0$,
$|a_nt^n|_p=p^{\varepsilon n-O(\log n)}\to\infty$.
(3) If $R_\xi$ and $R_{\xi_p}$ both converge on $|t|\le\rho$, $\rho>1$, then so does
$(\xi_p-\xi)A$, contradicting (2) unless $\xi=\xi_p$. $\square$

For $\mathbf C$: $c=9$, $\sigma_3=2$; for $\mathbf F$: $c=72$, $\sigma_3=2$; both integral,
$\kappa_3=0$ **VERIFIED** ($\max_{n\le400}v_3(a_n)=13$ resp. $11$; `radius.gp`).
This is the one *measured* hypothesis in the argument below (input **H1**).

---

## 3. The common cover of the C- and F-systems

Put $u=t_{\mathbf C}(\tau)$, $v=t_{\mathbf C}(2\tau)$, $x=t_{\mathbf F}(\tau)$. All the
following are **EXACT** identities of $q$-expansions to $O(q^{220})$
(`cover_checks.gp`); each is an identity between eta quotients of level dividing 12 and
so is provable by a Sturm/Ligozat check in the sense of Theorem A (the corresponding
Sturm bounds are $\le 24$; the certification is not written out here — input **H2**).

**(3.1) Modular equation of the C-curve.**
$$v^2-(9u^2-8u+1)\,v+u^2=0 .$$

**(3.2) The F-coordinate is rational on the cover.**
$$x=\frac{v-u}{9v-1}=\frac{v-u^2}{8v},\qquad\text{equivalently}\qquad v=\frac{x-u}{9x-1}.$$

**(3.3) The C-coordinate is quadratic over the F-coordinate.**
$$(9x-1)u^2+(1-8x)u+x(8x-1)=0,\qquad \operatorname{disc}=(1-8x)(1-6x)^2 .$$
Hence $\mathbb Q(u,v)=\mathbb Q(x)(u)$ is a **degree-2** extension of $\mathbb Q(x)$, with
$$e_1:=u+\bar u=\frac{8x-1}{9x-1},\qquad e_2:=u\bar u=\frac{x(8x-1)}{9x-1},$$
$$e_1':=v+\bar v=\frac{18x^2-10x+1}{(9x-1)^2},\qquad e_2':=v\bar v=\frac{x^2}{(9x-1)^2}.$$
The nontrivial automorphism $w$ of the cover acts by
$$\boxed{\ w(u)=\bar u=\frac{v-1}{9v-1},\qquad w(v)=\bar v=\frac{u-1}{9u-1}.\ }$$
(In modular terms $w$ is the Atkin–Lehner-type involution of the level-12 curve fixing
$t_{\mathbf F}$; it interchanges the two 2-isogeny directions.)

**(3.4) The weight-one forms.**
$$\frac{F_{\mathbf C}(2\tau)}{F_{\mathbf C}(\tau)}=\frac{1-3u}{1+3v},\qquad
G:=\frac{F_{\mathbf C}(\tau)}{F_{\mathbf C}(2\tau)}=\frac{1+3v}{1-3u},\qquad
\boxed{\ F_{\mathbf F}=\frac{F_{\mathbf C}(\tau)^2}{F_{\mathbf C}(2\tau)}\ }$$

The last identity is the exact statement of "the F-system is the $(1+V_2)$-twist of the
C-system on the common cover": $F_{\mathbf F}=G\cdot F_{\mathbf C}$.

---

## 4. The transfer identity

**Proposition 2 (transfer).** As formal power series in $x$ (branch $u=x+O(x^2)$),
$$A^{\mathbf F}(x)=G\cdot A^{\mathbf C}(u),\qquad
B^{\mathbf F}(x)=G\cdot B^{\mathbf C}(u)+\tfrac14G^2\,B^{\mathbf C}(v),$$
and consequently, for every $\xi$,
$$\boxed{\ B^{\mathbf F}(x)-\tfrac54\xi A^{\mathbf F}(x)
= G\cdot R^{\mathbf C}_\xi(u)\;+\;\tfrac14G^2\cdot R^{\mathbf C}_\xi(v)\ }\tag{4.1}$$
where $R^{\mathbf C}_\xi=B^{\mathbf C}-\xi A^{\mathbf C}$.

*Proof.* $A^{\mathbf F}(x)=F_{\mathbf F}=G F_{\mathbf C}=G A^{\mathbf C}(u)$ by (3.4).
For $B$: $B^{\mathbf F}(x)=F_{\mathbf F}\Theta_{\mathbf F}
=F_{\mathbf F}\bigl(\Theta_{\mathbf C}(\tau)+\tfrac14\Theta_{\mathbf C}(2\tau)\bigr)$ by §1, and
$$F_{\mathbf F}\Theta_{\mathbf C}(\tau)=G\,F_{\mathbf C}(\tau)\Theta_{\mathbf C}(\tau)=G\,B^{\mathbf C}(u),\qquad
F_{\mathbf F}\Theta_{\mathbf C}(2\tau)=G^2F_{\mathbf C}(2\tau)\Theta_{\mathbf C}(2\tau)=G^2B^{\mathbf C}(v),$$
using $F_{\mathbf F}=G^2F_{\mathbf C}(2\tau)$. Finally $G\,A^{\mathbf C}(v)\cdot G=G A^{\mathbf C}(u)=A^{\mathbf F}$,
so the $\xi$-terms assemble to $(1+\tfrac14)\xi A^{\mathbf F}$. $\square$

The factor $\tfrac54=1+2^{-(w+1)}=P_{\mathbf F}(2)/P_{\mathbf C}(2)$ is exactly the Mellin
eigenvalue ratio, as Proposition `prop:mellinshift` predicts.

**Proposition 3 ($w$-invariance).** Identity (4.1) also holds on the conjugate branch:
$$B^{\mathbf F}(x)-\tfrac54\xi A^{\mathbf F}(x)
= \bar G\cdot R^{\mathbf C}_\xi(\bar u)+\tfrac14\bar G^2 R^{\mathbf C}_\xi(\bar v),
\qquad \bar G=\frac{1+3\bar v}{1-3\bar u},$$
where $R^{\mathbf C}_\xi(\bar u)=\sum_k r_k\bar u^{\,k}$ is understood $3$-adically
(it converges by Lemma 1 and $|\bar u|_3\le1$).

*Reason.* Over $\mathbb C$ the derivation of Proposition 2 is an identity of analytic
functions of $\tau$ on the upper half-plane; the two branches of the cover over a point
$x$ are $\tau$ and $\gamma\tau$ for the Atkin–Lehner element $\gamma$ realising $w$,
and $t_{\mathbf F}(\gamma\tau)=t_{\mathbf F}(\tau)$, so the right-hand side of (4.1)
takes the same value on both. The transfer of this to the $3$-adic setting is input
**H3**. It was **VERIFIED** directly $3$-adically: with $\xi=\xi_3^{\mathbf C}$, both
$G R^{\mathbf C}(u)+\tfrac14G^2R^{\mathbf C}(v)$ and its $w$-conjugate agree with
$\sum_n(b^{\mathbf F}_n-\tfrac54\xi a^{\mathbf F}_n)x^n$ **to $3$-adic precision
$3^{44}$ in every $x$-coefficient up to $x^{25}$** (`galois_test.gp`).

---

## 5. Affinoid estimate and the main theorem

Fix $0<\delta<1$ rational and let $\mathbb D_\delta=\{|x|_3\le3^{\delta}\}$, with Tate
algebra $T_\delta$ and sup-norm $\|\cdot\|$.

**Lemma 4.** On $\mathbb D_\delta$: $|9x-1|=1$, $\|e_1\|,\|e_2\|,\|e_1'\|\le3^{\delta}$ and $\|e_2'\|\le3^{2\delta}$, hence every branch value satisfies
$|u|,|\bar u|,|v|,|\bar v|\le3^{\delta}<9$, and $|3u|,|3v|\le3^{\delta-1}<1$, so
$\|G\|=\|\bar G\|=1$.

*Proof.* $|9x|\le3^{\delta-2}<1$, so $|9x-1|=1$. Then
$|e_1|\le\max(|8x|,1)=\max(3^{\delta},1)=3^{\delta}$, $|e_2|\le3^{\delta}$, similarly for
$e_1',e_2'$. A root of $T^2-e_1T+e_2$ has $|u|\le\max(|e_1|,|e_2|^{1/2})\le3^{\delta}$.
Finally $|3u|\le3^{-1+\delta}<1$ gives $|1\pm3u|=1$. $\square$

**Theorem R (prime-to-$p$ Hecke functoriality; the case at hand).**
Assume **H1** ($\kappa_3=0$ for the F-row), **H2** (the cover identities of §3), **H3**
(Proposition 3). Then
$$\boxed{\ \xi_3^{\mathbf F}=\tfrac54\,\xi_3^{\mathbf C}=\frac{P_{\mathbf F}(2)}{P_{\mathbf C}(2)}\,\xi_3^{\mathbf C}. \ }$$
In particular Conjecture D holds for the pair $(\mathbf C,\mathbf F)$ at $p=3$.

*Proof.* Let $\xi=\xi_3^{\mathbf C}$ and $r_k=b^{\mathbf C}_k-\xi a^{\mathbf C}_k$; by
Lemma 1(1), $v_3(r_k)=2k+O(\log k)$, so $|r_k|\,3^{\delta k}\to0$ for $\delta<2$.
Let $\mathcal O=T_\delta[T]/(T^2-e_1T+e_2)$, a free $T_\delta$-algebra of rank 2 in which
$u,\bar u$ are the two roots; by Lemma 4 all of $u,v,\bar u,\bar v,G,\bar G$ are
power-bounded elements of norm $\le3^\delta$ resp. $=1$. Hence
$$H:=\sum_{k\ge0}r_k\Bigl(Gu^k+\tfrac14G^2v^k\Bigr)$$
converges in $\mathcal O$ (the $k$-th term has norm $\le|r_k|3^{\delta k}\to0$). By
Propositions 2 and 3, $w(H)=H$, so $H\in T_\delta$; and by Proposition 2 its expansion
in the branch $u=x+O(x^2)$ is the formal series
$R^{\mathbf F}:=\sum_n(b^{\mathbf F}_n-\tfrac54\xi a^{\mathbf F}_n)x^n$. Therefore
$R^{\mathbf F}\in T_\delta$, i.e. $B^{\mathbf F}-\tfrac54\xi A^{\mathbf F}$ converges on
$|x|_3\le3^{\delta}$ with $3^\delta>1$. By Lemma 1(3) applied to the F-row,
$\tfrac54\xi=\xi_3^{\mathbf F}$. $\square$

**Remark (why this is not formal).** The proof uses *three* separate facts about $p=3$:
(i) the correspondence $V_2$ is prime to 3, so the cover $\mathbb Q(u,v)/\mathbb Q(x)$ is
degree 2 and its discriminant $(1-8x)(1-6x)^2$, though it vanishes inside
$\mathbb D_\delta$, is harmless because we only use *symmetric* (Galois-invariant)
expressions; (ii) the modular units $1-3u$, $1+3v$ are $3$-adically close to 1 on
$\mathbb D_\delta$, which is what makes $G$ a unit of the Tate algebra; (iii) the slope
$\sigma_3=2$ of the C-row exceeds $\delta$, giving room in the composition. None of
(i)–(iii) survives $p\mid d$.

**Numerical confirmation** (`padic_limits.gp`, $N=400$, exact rationals):
$$v_3\bigl(4\,\xi^{\mathbf F}_3-5\,\xi^{\mathbf C}_3\bigr)=786,\qquad
v_3\bigl(\xi^{\mathbf B}_3-\xi^{\mathbf C}_3\bigr)=784,$$
where $\xi^X_3$ is replaced by $b^X_N/a^X_N$ ($\approx 2N=800$, the predicted rate).
Also $v_3(\xi)=-1$ for all three rows, and
$\xi^{\mathbf C}_3\equiv 3^{-1}\cdot 8386265965554334030$,
$\xi^{\mathbf F}_3\equiv 3^{-1}\cdot 4403999727414453137 \pmod{3^{39}}$,
$\xi^{\mathbf B}_3\equiv\xi^{\mathbf C}_3$.

---

## 6. The general statement

The proof of Theorem R used nothing about $\{1,V_2\}$ beyond the shape of the cover.
In general:

**Theorem R$'$ (expected form; proved here only in the case above).**
Let $(\Gamma,t,F)$ be a modular Apéry system with source $\Phi$ of weight $w+2$, and let
$C=\sum_{d\in\mathcal D}c_dV_d$ with $\gcd(d,p)=1$ for all $d\in\mathcal D$. Let
$(\Gamma',t',F')$ be a modular Apéry system whose source is $C\Phi$ and which is defined
on a cover $X'\to X$ that is *prime to $p$* in the sense that

1. $t(d\tau)$, $d\in\mathcal D$, generate a finite extension $L/\mathbb Q(t')$ of degree
   prime to $p$ whose branch points $x$ lie in $|x|_p\le p^{\delta}$ with all branch
   *values* $|t(d\tau)|_p\le p^{\delta}<p^{\sigma_p}$;
2. the modular-unit ratios $F(d\tau)/F'$ are units of the Tate algebra of
   $|x|_p\le p^{\delta}$;

then, with $\Theta'=D^{-(w+1)}C\Phi=\sum_dc_dd^{-(w+1)}V_d\Theta$,
$$\xi_p(C\Phi)=P_C(w+1)\,\xi_p(\Phi)\quad\text{whenever }P_C(w+1)\neq0,\qquad
\xi_p(C\Phi)=0\text{ if }P_C(w+1)=0 .$$
Equivalently: there is a single $\Lambda_p\in\mathbb Q_p$ (a $p$-adic avatar of
$L(\chi_{-3},2)$, resp. of the relevant period) with $\xi_p(C\Phi)=P_C(w+1)\Lambda_p$
for every prime-to-$p$ correspondence $C$ — which is precisely Conjecture D restricted
to the Hecke orbit of a single Eisenstein source.

**The $P_C(w+1)=0$ clause is not vacuous.** For $\Psi_0=(1-4V_2)S$ we have
$P(2)=1-4\cdot2^{-2}=0$; building the companion of $\Psi_0$ *on the C-curve*
(i.e. the row $\psi_n:=[t_{\mathbf C}^n]\,F_{\mathbf C}\cdot D^{-2}\Psi_0$) gives
(`psi0_test.gp`)
$$v_3(\psi_n/a^{\mathbf C}_n)= 18,36,56,75,93,115,135,153,176 \ \text{ at }n=10,\dots,90,$$
i.e. $2n-O(\log n)$: the $3$-adic Apéry limit of the target-zero class is **exactly $0$**,
at the full slope. This is the $p$-adic shadow of $L(\Psi_0,2)=0$ and is the strongest
single piece of evidence that $\xi_p$ is the $p$-adic realisation of the same
extension class as $\Lambda$.

---

## 7. Why B is *not* a prime-to-3 case

$\Phi_{\mathbf B}=\Phi_{\mathbf C}+2V_2\Psi_0$ involves only $V_1,V_2,V_4$, so at the
level of *sources* B is a prime-to-3 modification of C. But the B **coordinate** is not.
Two exact findings (`btest.gp`, PARI kernel computations over $\mathbb Q$):

* $t_{\mathbf B}\notin\mathbb Q(u,v)$: no relation
  $\sum_{i\le6,\,j\le1}(\alpha_{ij}+\beta_{ij}t_{\mathbf B})u^iv^j=0$ exists.
* $t_{\mathbf B}$ satisfies an irreducible **cubic** over $\mathbb Q(t_{\mathbf F})$:
  $$\boxed{\ (9x-1)\,t_{\mathbf B}\bigl(27t_{\mathbf B}^2-9t_{\mathbf B}+1\bigr)+x(8x-1)^2=0,\qquad x=t_{\mathbf F}.\ }\quad\textbf{EXACT}$$

So Zagier's B row is a degree-3 (i.e. $p$-power) reparametrisation of a level-12 source.
Correspondingly, the naive C-curve expansion of $\Theta_{\mathbf B}$ does **not**
converge $3$-adically: with $\tilde b_n:=[t_{\mathbf C}^n]F_{\mathbf C}\Theta_{\mathbf B}$
one measures $v_3(\tilde b_n/a^{\mathbf C}_n-\xi^{\mathbf C}_3)=-2,-1,-2$ at $n=10,20,30$
(`bcase_test.gp`) — bounded, not growing. (The reason is structural, not numerical:
$F_{\mathbf C}(\tau)\Theta_{\mathbf C}(2\tau)$ is a function on the *cover*, not on the
C-curve.)

Consequently Conjecture D for $(\mathbf B,\mathbf C)$ at $p=3$ decomposes as

* **(B1)** a prime-to-3 statement on the level-12 cover, of the type proved in Theorem R,
  giving $\xi_3$ of the source $\Phi_{\mathbf B}$ *in any level-12 parametrisation*
  equal to $\xi_3^{\mathbf C}$ (using $P_{2V_2(1-4V_2)}(2)=0$);
* **(B2)** **OPEN:** invariance of $\xi_p$ under a change of Apéry parametrisation
  $(t,F)\to(t',F')$ of the *same* source along a degree-$p$ cover.

(B2) is the same obstruction as the $p\mid d$ case of §8, and (B2) is what the
$3^{784}$ agreement of §5 is really testing.

---

## 8. The $p\mid d$ case: $\zeta(3)$ at $p=2$

For $\Lambda=\zeta(3)$, $r=3$, $w=2$: Domb $\alpha$ (level 12) has
$P_\alpha(s)=1-17\cdot2^{-s}-9\cdot3^{-s}+16\cdot4^{-s}+153\cdot6^{-s}-144\cdot12^{-s}$
and $\varepsilon=(12,4,16)$ (level 8) has
$P_\varepsilon(s)=(1-2^{-s})(1-4\cdot2^{-s})(1-16\cdot2^{-s})$. At $s=3$:
$$P_\alpha(3)=-\tfrac7{12},\qquad P_\varepsilon(3)=-\tfrac7{16},\qquad
\frac{P_\alpha(3)}{P_\varepsilon(3)}=\frac43=\frac{7/24}{7/32},$$
matching both the archimedean ratio $r_\alpha/r_\varepsilon$ and the measured $2$-adic
alignment $v_2(\Delta_n)=4n+O(1)$ of the verification record.

Here the connecting correspondence is $(1-V_2)(1-4V_2)(1-16V_2)$ at the prime $p=2$
itself. Every one of the three ingredients of Theorem R fails:

* $V_2$ at $p=2$ is the Frobenius/Verschiebung direction: the cover
  $X_0(2N)\to X_0(N)$ is **not** étale at $p$, and the "two branches" over a point of the
  base are not both in the residue disc — the whole point of the canonical subgroup
  theory is to select **one** of them.
* the discriminant of the cover is divisible by $p$, so Galois descent from the double
  cover to the base loses a factor of $p$ per application; the trace argument of
  Theorem R degenerates.
* the modular units $1\pm pu$ are no longer $\equiv1$: $\|G\|>1$.

The extra input needed is exactly Calegari-type overconvergence
(`Cal05`, and Coleman's theory): one must replace the affinoid $\mathbb D_\delta$ by the
*ordinary locus* $X(v)$, replace "both branches" by the canonical-subgroup section
$s:X(v)\to X_0(pN)(v)$, and re-prove Proposition 3 in the form
"$U_p$-eigen-decomposition of the Eichler integral": the transfer identity then reads
$\xi_p(C\Phi)=P_C^{(\mathrm{ord})}(w+1)\,\xi_p(\Phi)$ where $P_C^{(\mathrm{ord})}$ is the
Mellin polynomial evaluated on the **unit-root** factor only. Since
$\sigma_p=v_p(c)$ and $c=\lambda_1\lambda_2$ (Corollary `cor:product`), the unit-root
splitting is visible already in the recurrence, which is why the *numerics* work at
$p\mid d$ even though the geometric argument does not. Making this precise is the
natural next problem, and it would simultaneously settle (B2) of §7.

---

## 9. Verification record

Scripts in `lattice/rigidity/` (PARI/GP; run e.g. `timeout 560 gp -q <script>`):

| script | what | result |
|---|---|---|
| `lib.gp`, `setup.gp` | eta quotients $t_X,F_X$, $\Phi_X=F_XDt_X$ | — |
| `cover_checks.gp` | source identities; §3 (3.1)–(3.4); $\bar u,\bar v$; discriminant | all **EXACT** to $O(q^{220})$ |
| `psi0_test.gp` | $3$-adic limit of the $\Psi_0$-companion on the C-curve | $v_3(\psi_n/a_n)=2n-O(\log n)$, $n\le90$ |
| `galois_test.gp` | $w$-invariance of the transfer expression; $H=\bar H=R^{\mathbf F}$ | agree mod $3^{44}$, $x$-order $25$ |
| `bcase_test.gp` | $\Theta_{\mathbf B}$ decomposition; C-curve rows $=$ Zagier rows; failure of the naive B expansion | as stated |
| `btest.gp` | $t_{\mathbf B}\notin\mathbb Q(u,v)$; cubic for $t_{\mathbf B}$ over $\mathbb Q(t_{\mathbf F})$ | **EXACT** |
| `padic_limits.gp` | $\xi^{B,C,F}_3$ to $N=400$ | $v_3(4\xi^F-5\xi^C)=786$, $v_3(\xi^B-\xi^C)=784$ |
| `radius.gp` | $\kappa_3=0$ evidence | $\max_{n\le400}v_3(a_n)\le13$ |

## 10. Summary of logical status

* **Proved** (given H1–H3): $\xi_3^{\mathbf F}=\tfrac54\xi_3^{\mathbf C}$ (Theorem R).
  H1 = $\kappa_3=0$ for the F-row (verified $n\le400$); H2 = the eta-quotient cover
  identities of §3 (verified to $O(q^{220})$; a Sturm certification in the style of
  Theorem A would upgrade these to proofs and is routine); H3 = $w$-invariance of the
  transfer expression (proved over $\mathbb C$, verified $3$-adically).
* **Proved unconditionally:** the algebraic transfer identity (4.1); Lemma 1;
  Lemma 4; the identification of the ratio $\tfrac54$ with $P_{\mathbf F}(2)/P_{\mathbf C}(2)$.
* **New structural fact:** $t_{\mathbf B}$ is cubic over $\mathbb Q(t_{\mathbf F})$, so
  $(\mathbf B,\mathbf C)$ at $p=3$ is a $p\mid d$-type problem in disguise.
* **Open:** (B2) coordinate invariance of $\xi_p$ along degree-$p$ covers; the $p\mid d$
  case in general (Domb vs $\varepsilon$ at $p=2$); Conjecture D in full.
