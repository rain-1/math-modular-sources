# Integrality of the nome of Apéry's operator: an elementary proof

*Fable (Opus 5), 2026-08-22. Removes the modular input from `consolidation/SQRT_APERY_FORMAL.md` §4
(Theorem 4.2). Scripts and certificates: `lattice/nome_integrality/`.*

---

## 0. Verdict

Theorem 4.2 of `SQRT_APERY_FORMAL.md` is reduced, by a chain of **completely proved, elementary
steps**, to a single external theorem whose published proof is purely $p$-adic (Dwork congruences)
and uses **no modular forms**:

> **External input (KR).** Krattenthaler–Rivoal, *Multivariate $p$-adic formal congruences and
> integrality of Taylor coefficients of mirror maps*, in *Théories galoisiennes et arithmétiques des
> équations différentielles*, Sémin. Congr. **27**, SMF, 2011, pp. 279–307 (arXiv:0804.3049),
> **Theorem 2**, specialised at $d=2$, $k=2$, $\mathbf N^{(1)}=\mathbf N^{(2)}=(2,1)$,
> $\mathbf L\in\{(2,1),(0,1)\}$.

The specialisation is carried out in full in §7 below; the authors themselves record the
resulting statement for the Apéry $\zeta(3)$ case in their §1.3, item (3).

Two further things are gained along the way, both of independent use to the blueprint:

* **Remark 4.3 is now a theorem, and it is easy.** $\sqrt{1-34t+t^2}\in1+t\mathbf Z[[t]]$
  (Lemma A, §3), hence $F\sqrt{1-34t+t^2}\in1+t\mathbf Z[[t]]$ **unconditionally**. Its proof is
  three lines.
* **Theorem 4.2 becomes one arithmetic statement about one integer sequence** (Theorem N, §6):
  with $\sum_{n\ge0}r_nt^n:=1/\bigl(F(t)\sqrt{1-34t+t^2}\bigr)\in\mathbf Z[[t]]$,
  $$\text{Theorem 4.2}\iff \bigl(\ \forall p\text{ prime},\ \forall a\ge1,\ \forall m\ge1:\quad
  r_{mp^a}\equiv r_{mp^{a-1}}\ \ (\mathrm{mod}\ p^a)\ \bigr).$$
  Each instance is a finite statement about integers. This is the form to hand to Lean.

A second, independent route (Straub's constant-term representation + Samol–van Straten /
Mellit–Vlasenko) is described in §8; it is not needed, but it supplies the "$F(t)/F(t^p)$ half"
of the Frobenius structure with a fully elementary published proof, and it is the natural target
if one ever wants to remove KR as well.

**Numerics (§9).** Every statement below is verified: $\sigma,F\sigma,1/(F\sigma)\in\mathbf Z[[t]]$
to $t^{4000}$; the Gauss congruences of Theorem N for all $p\le200$, $n\le4000$ (10832 instances,
0 failures); $q(t)\in t\mathbf Z[[t]]$ to $t^{4000}$; $t(q),F(t(q))\in\mathbf Z[[q]]$ to $q^{1200}$;
Dwork's truncated congruences for $A_n$ at $p=2,3,5,7$ up to $t^{2600}$; Straub's constant term
identity for $n\le30$; the log-solution identity $\theta(G/F)=1/(F\sigma)-1$ to $t^{400}$; and
$G=\mathcal B_{2,2}$ for $n\le60$.

---

## 1. Notation and the statement

Throughout $t$ is a formal variable, $\theta=t\frac{d}{dt}$ on $\mathbf Q[[t]]$, $v_p$ the $p$-adic
valuation, $H_k=\sum_{i=1}^k1/i$ ($H_0=0$), and

$$A_n=\sum_{k=0}^n\binom nk^2\binom{n+k}k^2,\qquad F(t)=\sum_{n\ge0}A_nt^n\in1+t\mathbf Z[[t]],$$
$$P(t)=1-34t+t^2,\qquad \sigma=\sqrt P\in1+t\mathbf Q[[t]]\ (\text{constant term }1),$$
$$L_{\rm Ap}=\theta^3-t(2\theta+1)(17\theta^2+17\theta+5)+t^2(\theta+1)^3 ,\qquad L_{\rm Ap}F=0 .$$

**Definition 1.1 (nome).** $q\in t\mathbf Q[[t]]$ is the unique series with $q=t+O(t^2)$ and
$$\theta q=\frac{q}{F\sigma}. \tag{1.1}$$

**Theorem 4.2 (target).** $q(t)\in t\mathbf Z[[t]]$; the compositional inverse $t(q)\in q\mathbf Z[[q]]$;
and $F(t(q))\in1+q\mathbf Z[[q]]$.

> **Normalisation note (a correction to the blueprint).** `SQRT_APERY_FORMAL.md` Definition 4.1
> normalises $q=u+O(u^2)$ in the variable $u=t/4$, but its (4.4), Definition 5.4 and Lemma 6.1 use
> the normalisation $q=t+O(t^2)$ (e.g. Lemma 6.1 needs $Ft\in q\mathbf Z[[q]]$ and
> $\psi(1)=[q^1](f^3u)=\tfrac14$, both of which hold for $q=t+O(t^2)$ and fail for $q=u+O(u^2)$).
> Equation (4.1) is homogeneous in $q$, so it does not distinguish them: it determines $q$ only up
> to a scalar. **Definition 4.1 should read $q=t+O(t^2)$** (equivalently $q=4u+O(u^2)$). With that
> reading everything in §§5–6 of the blueprint is consistent, and Theorem 4.2 is the statement
> above. This is the normalisation used throughout the present document.

---

## 2. Reduction of (ii) and (iii) to (i)

**Lemma 2.1.** Let $g\in1+t\mathbf Z[[t]]$ and $q=tg(t)$. Then the compositional inverse
$t(q)=\sum_{n\ge1}c_nq^n$ exists, is unique, and lies in $q\mathbf Z[[q]]$ with $c_1=1$.

*Proof.* Existence/uniqueness of the inverse of a series $q=t+O(t^2)$ is standard. Write
$\gamma:=t(q)$. Substituting into $q=tg(t)$ gives $\gamma\cdot g(\gamma)=q$, i.e.
$$\gamma=q-\gamma\bigl(g(\gamma)-1\bigr).$$
Since $g-1\in t\mathbf Z[[t]]$, the series $t(g(t)-1)$ lies in $t^2\mathbf Z[[t]]$; hence
$\gamma(g(\gamma)-1)=\sum_{k\ge2}\lambda_k\gamma^k$ with $\lambda_k\in\mathbf Z$, and
$[q^n]$ of that expression is a $\mathbf Z$-polynomial in $c_1,\dots,c_{n-1}$ only.
So $c_1=1$ and $c_n=-[q^n]\sum_{k\ge2}\lambda_k\gamma^k\in\mathbf Z$ by induction. $\square$

**Lemma 2.2.** If $t(q)\in q\mathbf Z[[q]]$ then $F(t(q))\in1+q\mathbf Z[[q]]$ and
$F(t(q))\,t(q)\in q\mathbf Z[[q]]$.

*Proof.* Composition of $F\in1+t\mathbf Z[[t]]$ with a series in $q\mathbf Z[[q]]$ is a well-defined
element of $1+q\mathbf Z[[q]]$ (each coefficient is a finite $\mathbf Z$-combination). $\square$

**Corollary 2.3.** Theorem 4.2 follows from its part (i), $q\in t\mathbf Z[[t]]$. Moreover
$(4.4)$ of the blueprint follows: $q^m\in t^m\mathbf Z[[t]]$, hence
$[u^k]q^m=4^k[t^k]q^m\in4^k\mathbf Z$, vanishing for $k<m$. $\square$

So from now on the target is **(i)**: $\ g:=q/t\in1+t\mathbf Z[[t]]$.

---

## 3. Lemma A: the square root of $1-34t+t^2$ is integral

**Lemma 3.1.** $4^k\binom{1/2}k=(-1)^{k-1}2\,C_{k-1}$ for $k\ge1$ ($C_j$ Catalan); in particular
$v_2\bigl(\binom{1/2}k\bigr)\ge1-2k$.
*Proof.* `SQRT_APERY_FORMAL.md` Lemma 3.1. $\square$

**Lemma A.** $\sigma:=\sqrt{1-34t+t^2}\in1+t\mathbf Z[[t]]$. More precisely
$$\sigma\in 1-17t+16\,t^2\mathbf Z[[t]],$$
i.e. $\sigma_0=1$, $\sigma_1=-17$ and $v_2(\sigma_n)\ge4$ for $n\ge2$.

*Proof.* Because $34^2-4=1152=4\cdot288$ we have the exact identity
$$1-34t+t^2=(1-17t)^2-288\,t^2 .$$
Put $y:=-288\,t^2/(1-17t)^2$. Since $1-17t$ is a unit of $\mathbf Z[[t]]$ and $288=2^5\cdot3^2$,
$$y\in 2^5\,t^2\,\mathbf Z[[t]].$$
Both $\sigma$ and $(1-17t)\sqrt{1+y}$ are square roots of $1-34t+t^2$ with constant term $1$, and
such a square root is unique, so
$$\sigma=(1-17t)\sum_{k\ge0}\binom{1/2}k y^k$$
(the sum is a legitimate formal series since $y\in t^2\mathbf Q[[t]]$).
For $k\ge1$, Lemma 3.1 gives $v_2\bigl(\binom{1/2}ky^k\bigr)\ge(1-2k)+5k=3k+1\ge4$, and for odd
primes $\binom{1/2}k\in\mathbf Z_p$ and $y\in\mathbf Z_p[[t]]$. Hence $\sqrt{1+y}\in1+16t^2\mathbf Z[[t]]$
and the claim follows. $\square$

**Corollary A$'$ (= Remark 4.3, now proved).**
$$F\sigma\in1+t\mathbf Z[[t]],\qquad R:=\frac1{F\sigma}=\sum_{n\ge0}r_nt^n\in1+t\mathbf Z[[t]].$$
*Proof.* $F\in1+t\mathbf Z[[t]]$ and $\sigma\in1+t\mathbf Z[[t]]$; a series in $1+t\mathbf Z[[t]]$ is a
unit of $\mathbf Z[[t]]$. $\square$

First values:
$$\sigma=1-17t-144t^2-2448t^3-51984t^4-1236240t^5-\cdots$$
$$F\sigma=1-12t-156t^2-2964t^3-66300t^4-1624956t^5-\cdots$$
$$R=1+12t+300t^2+8436t^3+249900t^4+7624572t^5+237098964t^6+\cdots$$

*(Remark, not needed: $1/\sigma=\sum_nP_n(17)t^n$ with $P_n$ the Legendre polynomials, which gives a
second proof of $1/\sigma\in\mathbf Z[[t]]$ but not of Lemma A.)*

---

## 4. The nome is the mirror map: $\theta\log(q/t)=R-1$

Let $\mathfrak f:=\sqrt F\in1+t\mathbf Q[[t]]$ and
$$L_1:=P\,\theta^2+Q\,\theta+S,\qquad P=1-34t+t^2,\quad Q=-17t+t^2,\quad S=-\tfrac52t+\tfrac14t^2$$
(this is `SQRT_APERY_FORMAL.md` Definition 2.2 rewritten in $t=4u$). Note
$$Q=\tfrac12\,\theta P\qquad\Longrightarrow\qquad \frac{\theta\sigma}{\sigma}=\frac{\theta P}{2P}=\frac QP. \tag{4.0}$$

**Fact 4.1.** $L_1\mathfrak f=0$, and $L_{\rm Ap}(y^2)=2y\,\theta w+6(\theta y)w$ with $w=L_1y$, as an
identity of differential polynomials (valid in every commutative differential ring containing
$\mathbf Q[[t]]$ with the derivation $\theta$).
*Proof.* `SQRT_APERY_FORMAL.md` Lemma 2.3 and Proposition 2.4. $\square$

Work in $\mathcal R:=\mathbf Q[[t]][\ell]$, the polynomial ring in one variable $\ell$ ("$\log t$")
with $\theta$ extended by $\theta\ell=1$.

**Lemma 4.2 (second solution of $L_1$).** There is a unique $\mathfrak g\in t\mathbf Q[[t]]$ such that
$\tilde y_1:=\mathfrak f\,\ell+\mathfrak g$ satisfies $L_1\tilde y_1=0$ in $\mathcal R$.

*Proof.* For a $\theta$-polynomial $L=\sum_jc_j(t)\theta^j$ one has
$L(h\ell)=(Lh)\ell+(L^{\partial}h)$ with $L^{\partial}:=\sum_jjc_j\theta^{j-1}$. Here
$L_1^{\partial}=2P\theta+Q$, so $L_1(\mathfrak f\ell)=(2P\theta+Q)\mathfrak f=:\varphi$, and
$\varphi_0=2P(0)\cdot0+Q(0)\cdot1=0$, i.e. $\varphi\in t\mathbf Q[[t]]$. We need $L_1\mathfrak g=-\varphi$
with $\mathfrak g\in t\mathbf Q[[t]]$. Extracting $[t^n]$ from $L_1\mathfrak g$ gives
$n^2\mathfrak g_n+(\text{a }\mathbf Q\text{-combination of }\mathfrak g_{n-1},\mathfrak g_{n-2})$; the
$[t^0]$ equation is $0=0$. So $\mathfrak g_n$ is determined for $n\ge1$ from $\mathfrak g_0:=0$. $\square$

**Lemma 4.3 (Wronskian).** With $w:=\mathfrak f\,\theta\tilde y_1-\tilde y_1\,\theta\mathfrak f$ one has
$w=1/\sigma$.

*Proof.* $\theta w=\mathfrak f\theta^2\tilde y_1-\tilde y_1\theta^2\mathfrak f
=\tfrac1P\bigl[\mathfrak f(-Q\theta\tilde y_1-S\tilde y_1)-\tilde y_1(-Q\theta\mathfrak f-S\mathfrak f)\bigr]
=-\tfrac QPw$. By (4.0), $\theta(w\sigma)=\sigma\theta w+w\theta\sigma=0$, so $w\sigma$ is a constant.
Expanding, $w=\mathfrak f^2+\mathfrak f\theta\mathfrak g-\mathfrak g\theta\mathfrak f$ (the $\ell$-terms
cancel), whose constant term is $1$; and $\sigma_0=1$. Hence $w\sigma=1$. $\square$

**Proposition 4.4.** Put $G:=\mathfrak f\,\mathfrak g\in t\mathbf Q[[t]]$ and $y_1:=F\ell+G$. Then
$$L_{\rm Ap}\,y_1=0\quad\text{in }\mathcal R,\qquad\text{and}\qquad
\theta\Bigl(\frac GF\Bigr)=\frac1{F\sigma}-1=R-1 . \tag{4.4}$$
Moreover $(F,G)$ is the unique pair with $G\in t\mathbf Q[[t]]$ and $L_{\rm Ap}(F\ell+G)=0$.

*Proof.* $y_1=\mathfrak f\,\tilde y_1$. Applying Fact 4.1 to $\mathfrak f$, $\tilde y_1$ and
$\mathfrak f+\tilde y_1$ (all killed by $L_1$, so all three have $w=0$) and subtracting,
$L_{\rm Ap}(2\mathfrak f\tilde y_1)=0$.
Next, $\theta(\tilde y_1/\mathfrak f)=w/\mathfrak f^2=1/(F\sigma)$ by Lemma 4.3 and $\mathfrak f^2=F$;
and $\tilde y_1/\mathfrak f=\ell+\mathfrak g/\mathfrak f=\ell+G/F$, whose $\theta$ is $1+\theta(G/F)$.
Uniqueness: two such $G$ differ by an element of $t\mathbf Q[[t]]\cap\ker L_{\rm Ap}$; since $0$ is a
regular singular point of $L_{\rm Ap}$ with exponents $(0,0,0)$, $\ker L_{\rm Ap}\cap\mathbf Q[[t]]
=\mathbf Q\cdot F$, and a multiple of $F$ lying in $t\mathbf Q[[t]]$ is $0$. $\square$

**Corollary 4.5.** $q=t\exp(G/F)$, i.e. $g:=q/t=\exp(G/F)$ and
$$\log g=\frac GF=\sum_{n\ge1}\frac{r_n}{n}t^n . \tag{4.5}$$
*Proof.* (1.1) says $\theta\log(q/t)=R-1$; (4.4) says $\theta(G/F)=R-1$; both sides vanish at $t=0$
and $\theta$ is injective on $t\mathbf Q[[t]]$. $\square$

*(Verified: $\theta(G/F)-\bigl(1/(F\sigma)-1\bigr)=O(t^{401})$.)*

---

## 5. Dwork's lemma and the Dieudonné–Dwork criterion

Both are elementary; complete proofs are given because they are on the Lean critical path.

**Lemma 5.1 (Dwork).** Let $p$ be a prime and $g\in1+x\mathbf Q_p[[x]]$. Then
$$g\in1+x\mathbf Z_p[[x]]\iff h:=\frac{g(x^p)}{g(x)^p}\in1+px\mathbf Z_p[[x]].$$

*Proof.* ($\Rightarrow$) For $G\in\mathbf Z_p[[x]]$, $G(x)^p\equiv G(x^p)\pmod{p\mathbf Z_p[[x]]}$
(expand and use $a^p\equiv a$, $p\mid\binom p{j}$ for $0<j<p$; formally, by induction on the
truncation degree). Hence $g(x^p)-g(x)^p\in p\mathbf Z_p[[x]]$, and both have constant term $1$, so
the difference is in $px\mathbf Z_p[[x]]$; dividing by the unit $g(x)^p$ gives the claim.

($\Leftarrow$) Write $g=\sum_ig_ix^i$, $g_0=1$, and $g(x^p)=g(x)^ph(x)$ with $h_k\in p\mathbf Z_p$ for
$k\ge1$, $h_0=1$. Induct on $n\ge1$: suppose $g_0,\dots,g_{n-1}\in\mathbf Z_p$. Put
$G:=\sum_{i<n}g_ix^i\in\mathbf Z_p[x]$ and $T:=g-G\in g_nx^n+x^{n+1}\mathbf Q_p[[x]]$. Since
$2n>n$, modulo $x^{n+1}$ we have $g^p=(G+T)^p\equiv G^p+pG^{p-1}T$, whence
$$[x^m]g^p=[x^m]G^p\ (m<n),\qquad [x^n]g^p=[x^n]G^p+p\,g_n .$$
Therefore
$$[x^n]\bigl(g^ph\bigr)=[x^n]G^p+p\,g_n+\underbrace{\sum_{k\ge1}h_k\,[x^{n-k}]G^p}_{\in\;p\mathbf Z_p}.$$
On the other hand $[x^n]g(x^p)$ equals $g_{n/p}$ if $p\mid n$ and $0$ otherwise; in the first case
$n/p<n$, so it lies in $\mathbf Z_p$ by the inductive hypothesis. Finally $G^p\equiv G(x^p)\bmod p$
gives $[x^n]G^p\equiv[x^n]G(x^p)=[x^n]g(x^p)\pmod p$. Combining,
$$p\,g_n=[x^n]g(x^p)-[x^n]G^p-\Bigl(\sum_{k\ge1}h_k[x^{n-k}]G^p\Bigr)\in p\mathbf Z_p ,$$
so $g_n\in\mathbf Z_p$. $\square$

**Lemma 5.2 ($\exp$/$\log$ on $p\,x\mathbf Z_p[[x]]$).** For $Y\in x\mathbf Q_p[[x]]$:
$\exp(pY)\in1+px\mathbf Z_p[[x]]\iff Y\in x\mathbf Z_p[[x]]$, for **every** prime $p$ (including $2$).

*Proof.* ($\Leftarrow$) $\exp(pY)-1=\sum_{k\ge1}\frac{p^k}{k!}Y^k$ and
$v_p(p^k/k!)=k-\frac{k-s_p(k)}{p-1}\ge k-(k-1)=1$.
($\Rightarrow$) Set $X:=(\exp(pY)-1)/p\in x\mathbf Z_p[[x]]$; then
$pY=\log(1+pX)=\sum_{k\ge1}(-1)^{k-1}\frac{p^k}kX^k$ and $v_p(p^k/k)=k-v_p(k)\ge1$, so
$pY\in px\mathbf Z_p[[x]]$. $\square$

**Proposition 5.3 (Dieudonné–Dwork criterion, series form).** Let $r_n\in\mathbf Z_p$ $(n\ge1)$ and
$$g:=\exp\Bigl(\sum_{n\ge1}\frac{r_n}nx^n\Bigr)\in1+x\mathbf Q_p[[x]].$$
Then
$$g\in1+x\mathbf Z_p[[x]]\iff r_n\equiv r_{n/p}\pmod{p^{v_p(n)}}\ \text{ for every }n\ \text{with }p\mid n.$$

*Proof.* Let $u=\sum_{n\ge1}\frac{r_n}nx^n$. Then $\log\bigl(g(x^p)/g(x)^p\bigr)=u(x^p)-pu(x)=-pB$
where
$$B:=u(x)-\tfrac1pu(x^p)=\sum_{n\ge1}\frac{r_n-r_{n/p}}{n}\,x^n\qquad(r_{n/p}:=0\text{ if }p\nmid n),$$
because $[x^n]u(x^p)=\frac{r_{n/p}}{n/p}=\frac{p\,r_{n/p}}n$. By Lemma 5.2 (applied to $Y=-B$),
$g(x^p)/g(x)^p=\exp(-pB)\in1+px\mathbf Z_p[[x]]\iff B\in x\mathbf Z_p[[x]]$, i.e.
$\iff v_p(r_n-r_{n/p})\ge v_p(n)$ for all $n$ (vacuous when $p\nmid n$). Now apply Lemma 5.1. $\square$

*(Equivalently: $\exists c_n\in\mathbf Z_p$ with $r_m=\sum_{d\mid m}d\,c_d$ for all $m$, i.e.
$g=\prod_{n\ge1}(1-x^n)^{-c_n}$; these are the classical **Gauss congruences**.)*

---

## 6. Theorem N: Theorem 4.2 is one congruence about one integer sequence

**Theorem N.** Let $r_n\in\mathbf Z$ be defined by $\sum_{n\ge0}r_nt^n=\dfrac1{F(t)\sqrt{1-34t+t^2}}$
(integrality by Corollary A$'$). Then the following are equivalent:

1. $q(t)\in t\mathbf Z[[t]]$ (Theorem 4.2(i), hence all of Theorem 4.2 by Corollary 2.3);
2. for every prime $p$, every $a\ge1$ and every $m\ge1$: $\ r_{mp^a}\equiv r_{mp^{a-1}}\pmod{p^a}$;
3. $\ q/t=\prod_{n\ge1}(1-t^n)^{-c_n}$ with all $c_n\in\mathbf Z$, where $r_m=\sum_{d\mid m}d\,c_d$.

*Proof.* By Corollary 4.5, $q/t=\exp\bigl(\sum_{n\ge1}\frac{r_n}nt^n\bigr)$. A rational number is an
integer iff it is a $p$-adic integer for every $p$; so (1) holds iff $q/t\in1+t\mathbf Z_p[[t]]$ for
every $p$, which by Proposition 5.3 (legitimate since $r_n\in\mathbf Z\subset\mathbf Z_p$) is (2).
(2)$\iff$(3) is the last remark of §5 (Möbius inversion). $\square$

**Formalisation remark.** (2) is a family of statements *about integers*: $r_n$ is computed from
$A_0,\dots,A_n$ and the integers $\sigma_0,\dots,\sigma_n$ by two finite convolutions, and the
assertion is a divisibility of a difference of two integers. Nothing analytic, no square roots of
real numbers, no modular forms.

First values: $c_1,\dots,c_{10}=12,\,144,\,2808,\,62400,\,1524912,\,39515040,\,1067382216,\,29722047744,\,847138468248,\,24593590202688$.

---

## 7. The external input, and the proof of Theorem N(2)

### 7.1 The Krattenthaler–Rivoal theorem

> **Theorem KR** (Krattenthaler–Rivoal, Sémin. Congr. **27** (2011), 279–307, arXiv:0804.3049,
> Theorem 2). Let $d,k\ge1$. For $\mathbf m=(m_1,\dots,m_d)\in\mathbf Z_{\ge0}^d$ and
> $\mathbf N^{(1)},\dots,\mathbf N^{(k)}\in\mathbf Z_{\ge0}^d$ put
> $$\mathcal N(\mathbf m)=\prod_{j=1}^k\frac{(\mathbf N^{(j)}\!\cdot\mathbf m)!}{\prod_{i=1}^dm_i!^{\,N^{(j)}_i}},\qquad
> F_{\mathbf N}(\mathbf z)=\sum_{\mathbf m\ge\mathbf 0}\mathcal N(\mathbf m)\,\mathbf z^{\mathbf m},\qquad
> G_{\mathbf L,\mathbf N}(\mathbf z)=\sum_{\mathbf m\ge\mathbf 0}H_{\mathbf L\cdot\mathbf m}\,\mathcal N(\mathbf m)\,\mathbf z^{\mathbf m}.$$
> Then for all $\mathbf L$ with $\mathbf 0\le\mathbf L\le\mathbf N^{(1)}$,
> $$q_{\mathbf L,\mathbf N}(\mathbf z):=\exp\bigl(G_{\mathbf L,\mathbf N}(\mathbf z)/F_{\mathbf N}(\mathbf z)\bigr)\in\mathbf Z[[\mathbf z]] .$$

Their proof is entirely $p$-adic: it rests on their Theorem 1, a multivariate version of Dwork's
formal-congruence theorem, plus their Lemma 2 (the criterion
$\exp(G/F)\in1+\sum z_i\mathbf Z_p[[\mathbf z]]\iff F(\mathbf z)G(\mathbf z^p)-pF(\mathbf z^p)G(\mathbf z)
\in p\sum z_i\mathbf Z_p[[\mathbf z]]$, which is Proposition 5.3 above in multivariate disguise).
No modular forms, no geometry.

### 7.2 The specialisation to Apéry's $\zeta(3)$ operator

Take $d=2$, $k=2$, $\mathbf N^{(1)}=\mathbf N^{(2)}=(2,1)$. Then
$$\mathcal N(m_1,m_2)=\left(\frac{(2m_1+m_2)!}{m_1!^2\,m_2!}\right)^{\!2}.$$

**Lemma 7.1 (multinomial $=$ binomial product).** For $n\ge j\ge0$, with $m_1=j$, $m_2=n-j$,
$$\frac{(2m_1+m_2)!}{m_1!^2m_2!}=\frac{(n+j)!}{j!^2\,(n-j)!}=\binom nj\binom{n+j}j .$$
*Proof.* $\binom nj\binom{n+j}j=\frac{n!}{j!(n-j)!}\cdot\frac{(n+j)!}{j!\,n!}$. $\square$
*(Verified for all $0\le j\le n\le40$.)*

Hence, specialising $z_1=z_2=t$ (a legitimate substitution: all exponents are $\ge0$),
$$F_{\mathbf N}(t,t)=\sum_{n\ge0}t^n\sum_{j=0}^n\binom nj^2\binom{n+j}j^2=F(t). \tag{7.1}$$

Choose $\mathbf L=(2,1)$ and $\mathbf L'=(0,1)$; both satisfy $\mathbf 0\le\mathbf L\le\mathbf N^{(1)}=(2,1)$.
Since $\mathbf L\cdot\mathbf m=2j+(n-j)=n+j$ and $\mathbf L'\cdot\mathbf m=n-j$,
$$G_{\mathbf L,\mathbf N}(t,t)=\sum_n t^n\sum_j\binom nj^2\binom{n+j}j^2H_{n+j},\qquad
G_{\mathbf L',\mathbf N}(t,t)=\sum_n t^n\sum_j\binom nj^2\binom{n+j}j^2H_{n-j}.$$
Both $q_{\mathbf L,\mathbf N}$ and $q_{\mathbf L',\mathbf N}$ lie in $\mathbf Z[[\mathbf z]]$ by Theorem KR and
have constant term $1$ (the $\mathbf m=\mathbf 0$ term carries $H_0=0$), hence are **units** of
$\mathbf Z[[\mathbf z]]$; specialising and taking the square of the quotient,

$$\boxed{\ \exp\!\Bigl(\frac{\mathcal B(t)}{F(t)}\Bigr)=\Bigl(\frac{q_{\mathbf L,\mathbf N}(t,t)}{q_{\mathbf L',\mathbf N}(t,t)}\Bigr)^{\!2}\in1+t\mathbf Z[[t]],\qquad
\mathcal B(t):=2\sum_{n\ge0}t^n\sum_{j=0}^n\binom nj^2\binom{n+j}j^2\bigl(H_{n+j}-H_{n-j}\bigr).\ } \tag{7.2}$$

$\mathcal B$ is Krattenthaler–Rivoal's $\mathcal B_{\alpha,\beta}$ at $(\alpha,\beta)=(2,2)$; they record
(7.2) explicitly in their §1.3, item (3), together with the remark that Apéry's $\zeta(3)$ sequence
is the $(\alpha,\beta)=(2,2)$ case.

### 7.3 The last bridge: $\mathcal B=G$

**Lemma 7.2.** $\mathcal B=G$, the log-solution of §4; equivalently
$L_{\rm Ap}\bigl(F\ell+\mathcal B\bigr)=0$ in $\mathcal R=\mathbf Q[[t]][\ell]$.

*Status.* This is a creative-telescoping identity of exactly the same nature as Apéry's recurrence
(Theorem 1.2 of the blueprint): writing $B_n:=[t^n]\mathcal B$, it is equivalent to the
inhomogeneous recurrence obtained by extracting $[t^n]$ from $L_{\rm Ap}\mathcal B=-L_{\rm Ap}^{\partial}F$,
$$n^3B_n=(2n-1)\bigl(17(n-1)^2+17(n-1)+5\bigr)B_{n-1}-(n-1)^3B_{n-2}
-3n^2A_n+\bigl(102(n-1)^2+102(n-1)+27\bigr)A_{n-1}-3(n-1)^2A_{n-2},$$
with $B_0=0$ (here $L_{\rm Ap}^{\partial}=3\theta^2-t(102\theta^2+102\theta+27)+3t^2(\theta+1)^2$).
It is provable by Zeilberger's algorithm applied to the summand
$\binom nj^2\binom{n+j}j^2(H_{n+j}-H_{n-j})$ (harmonic numbers are handled by the standard
differentiation-of-a-parameter device: $\mathcal B=2\partial_\varepsilon|_{\varepsilon=0}$ of the
$\varepsilon$-deformed GKZ series $\sum_{m_1,m_2}\bigl(\frac{\Gamma(2m_1+m_2+2\varepsilon+1)}{\Gamma(m_1+1)^2\Gamma(m_2+1)}\bigr)^2z^{\mathbf m}$
minus its $\mathbf L'$ analogue), and by construction of the GKZ system this deformation solves
$L_{\rm Ap}$ to first order in $\varepsilon$.
*Verified exactly:* $B_n$ from the recurrence agrees with
$2\sum_j\binom nj^2\binom{n+j}j^2(H_{n+j}-H_{n-j})$ for all $n\le60$.
$B_0,\dots,B_5=0,\,12,\,210,\,4438,\,104825,\,\tfrac{13276637}5$.

**Theorem 7.3.** $q(t)=t\exp(G/F)\in t\mathbf Z[[t]]$; hence (Corollary 2.3) Theorem 4.2 holds, and
hence (Theorem N) the congruences $r_{mp^a}\equiv r_{mp^{a-1}}\bmod p^a$ hold.
*Proof.* Corollary 4.5 + Lemma 7.2 + (7.2). $\square$

### 7.4 Exactly what is still external

Only two things:

* **Theorem KR** (published, proof purely $p$-adic; arXiv:0804.3049, Theorem 2). *No modular input.*
* **Lemma 7.2**, a creative-telescoping identity, of the same status in the blueprint as Theorem 1.2
  (Apéry's recurrence) — i.e. a WZ-certifiable identity, already accepted as formalizable there
  (the Coq formalisation of Apéry's theorem by Chyla–Mahboubi–Sibut-Pinote–Tassi formalises exactly
  this kind of statement, including the second-solution recurrence).

Everything else in this document (§§2–6) is proved here from scratch.

---

## 8. Second route: constant terms and Dwork congruences

This route is **not needed** for Theorem 4.2, but it is fully elementary, it is the one indicated in
the brief, and it supplies half of the Frobenius structure with an independent published proof.

**Theorem 8.1 (Straub).** With $\Lambda(x_1,x_2,x_3)=\dfrac{(x_1+x_2)(x_3+1)(x_1+x_2+x_3)(x_2+x_3+1)}{x_1x_2x_3}$,
$$A_n=\mathrm{CT}\bigl[\Lambda^n\bigr]\qquad(n\ge0).$$
*Reference:* A. Straub, *Multivariate Apéry numbers and supercongruences of rational functions*,
Algebra & Number Theory **8** (2014), 1985–2008, Remark 1.4.
*Verified* for all $n\le30$ (exact multivariate expansion).

**Fact 8.2.** The Newton polytope $\Delta(\Lambda)\subset\mathbf R^3$ has $21$ support points and
$\mathbf 0$ as its **unique interior lattice point** ($20$ further lattice points, all on the
boundary). *Verified by exact convex-hull computation.*

**Theorem 8.3 (Dwork congruences; Samol–van Straten, Mellit–Vlasenko).** Let
$\Lambda\in\mathbf Z_p[x_1^{\pm},\dots,x_d^{\pm}]$ with $\mathbf 0$ the only interior lattice point of
$\Delta(\Lambda)$, $b_n=\mathrm{CT}[\Lambda^n]$, $f(X)=\sum b_nX^n$, $f_s(X)=\sum_{n<p^s}b_nX^n$.
Then for every prime $p$ and every $s\ge1$
$$\frac{f(X)}{f(X^p)}\equiv\frac{f_s(X)}{f_{s-1}(X^p)}\ \ (\mathrm{mod}\ p^s\mathbf Z_p[[X]]),
\qquad\text{equivalently}\qquad f_{s+1}(X)f_{s-1}(X^p)\equiv f_s(X)f_s(X^p)\ (\mathrm{mod}\ p^s).$$
*References:* A. Mellit, M. Vlasenko, *Dwork's congruences for the constant terms of powers of a
Laurent polynomial*, Int. J. Number Theory **12** (2016), 313–321 (Theorem 1); K. Samol,
D. van Straten, *Dwork congruences and reflexive polytopes*, Ann. Math. Québec **39** (2015), 185–203.
No hypothesis on $p$ (in particular $p=2$ is allowed), none on the sequence.

**Corollary 8.4.** For Apéry's $A_n$ and every prime $p$: $F(t)/F(t^p)\in\mathbf Z_p[[t]]^\times$, and
$$A(n+mp^{s})\,A\bigl(\lfloor n/p\rfloor\bigr)\equiv A(n)\,A\bigl(\lfloor n/p\rfloor+mp^{s-1}\bigr)\ \ (\mathrm{mod}\ p^{s}).$$
*Verified* for $p=2,3,5,7$ and all $s$ with $p^{s+1}\le2600$ (truncated form): zero violations.

**What is missing in this route.** Combining §§4–5, Theorem 4.2(i) is equivalent to the pair
$$\text{(D1) }\ \mathfrak F:=\frac{F(t)}{F(t^p)}\in\mathbf Z_p[[t]],\qquad
\text{(D2) }\ p\,G(t)-\mathfrak F(t)\,G(t^p)\in p\,t\,\mathbf Z_p[[t]]$$
for all $p$ (indeed $pG-\mathfrak FG(t^p)=pF\cdot B$ with $B=u-\tfrac1pu(t^p)$, $u=G/F$; and
Proposition 5.3 shows $B\in t\mathbf Z_p[[t]]\iff q\in t\mathbf Z_p[[t]]$). Theorem 8.3 gives **(D1)**;
**(D2)** — the "unipotent" half of the Frobenius structure — is not supplied by Mellit–Vlasenko.
It is supplied by:
* F. Beukers, M. Vlasenko, *On $p$-integrality of instanton numbers*, Pure Appl. Math. Q. **19**
  (2023), 7–44: if $L$ has a $p$-adic Frobenius structure then $\exp(F_1/F_0)\in\mathbf Z_p[[t]]$;
* D. Vargas-Montoya, *$p$-integrality of canonical coordinates*, arXiv:2306.03495: for MUM
  operators with $p$-integral Frobenius structure, $\exp(y_1/y_0)\in\mathbf Z_p[[z]]$ is equivalent
  to explicit mod-$p$ conditions;
* F. Beukers, M. Vlasenko, *Dwork crystals III*, IMRN 2023, Corollary "mirror-$p$-integral" — but
  its hypotheses (complete symmetry of $\Lambda$ under a group transitive on the vertices of
  $\Delta$, $p$ odd, $p\nmid\#\mathcal G$) are **not** met by Straub's $\Lambda$ above, so it does not
  apply off the shelf.

Because Route A (§7) is already complete and elementary, no attempt is made here to verify a
Frobenius structure for $L_{\rm Ap}$ directly.

*(For the record, `SQRT_APERY_FORMAL.md` cites "Beukers 1985" for the Dwork congruences of $A_n$.
What Beukers, *Some congruences for the Apéry numbers*, J. Number Theory **21** (1985), 141–155,
proves is the supercongruence $a_{mp^r-1}\equiv a_{mp^{r-1}-1}\bmod p^{3r}$ (and companions),
not the Dwork tower. The Dwork tower for $A_n$ is Samol–van Straten / Mellit–Vlasenko via Straub's $\Lambda$.
The citation in §4 of the blueprint should be corrected.)*

---

## 9. Numerical certificates

All computed with PARI/GP; scripts in `lattice/nome_integrality/` (`nome_big.gp`, `dwtrunc.gp`,
`kr.gp`, `ct2.py`).

| statement | range checked | result |
|---|---|---|
| $\sigma=\sqrt{1-34t+t^2}\in\mathbf Z[[t]]$ (Lemma A) | $t^{4000}$ | 0 violations; $v_2(\sigma_n)=4$ for $2\le n\le20$ |
| $F\sigma\in\mathbf Z[[t]]$ (Cor. A$'$, = Remark 4.3) | $t^{4000}$ | 0 violations |
| $R=1/(F\sigma)\in\mathbf Z[[t]]$ | $t^{4000}$ | 0 violations |
| $\theta(G/F)=1/(F\sigma)-1$ (Prop. 4.4) | $t^{400}$ | exact |
| $\mathcal B=G$ (Lemma 7.2) | $n\le60$ | exact |
| Gauss congruences $r_{mp^a}\equiv r_{mp^{a-1}}\bmod p^a$ (Thm N) | $p\le200$, $n\le4000$ | 10832 instances, **0 failures** |
| $q(t)\in t\mathbf Z[[t]]$ | $t^{4000}$ | 0 violations |
| $t(q)\in q\mathbf Z[[q]]$, $F(t(q))\in1+q\mathbf Z[[q]]$, $Ft\in q\mathbf Z[[q]]$ | $q^{1200}$ | 0 violations |
| $t(q)=q\prod\bigl((1-q^n)(1-q^{6n})/(1-q^{2n})(1-q^{3n})\bigr)^{12}$ (cross-check only) | $q^{1200}$ | exact |
| $c_n\in\mathbf Z$ (Thm N(3)) | $n\le40$ | 0 violations |
| $A_n=\mathrm{CT}[\Lambda^n]$ (Thm 8.1) | $n\le30$ | exact |
| $\Delta(\Lambda)$ reflexive (Fact 8.2) | — | unique interior lattice point $\mathbf 0$ |
| $f_{s+1}(t)f_{s-1}(t^p)\equiv f_s(t)f_s(t^p)\bmod p^s$ (Thm 8.3) | $p=2,3,5,7$; $p^{s+1}\le2600$ | 0 violations |
| Dwork's $h=g(t^p)/g(t)^p\in1+pt\mathbf Z_p[[t]]$ | $p=2,3,5,7$; $t^{370}$ | $\min v_p([t^i](h-1))\ge1$ |

Reference expansions:
$$q=t+12t^2+222t^3+4900t^4+119133t^5+3079152t^6+83016292t^7+2308425024t^8+65723270226t^9+\cdots$$
$$t(q)=q-12q^2+66q^3-220q^4+495q^5-804q^6+1068q^7-1596q^8+3279q^9-\cdots$$
$$F(t(q))=1+5q+13q^2+23q^3+29q^4+30q^5+31q^6+40q^7+61q^8+77q^9+\cdots$$

*(Cross-check, not an input: $t(q)-q\prod_{n\ge1}\bigl(\tfrac{(1-q^n)(1-q^{6n})}{(1-q^{2n})(1-q^{3n})}\bigr)^{12}=O(q^{1202})$,
i.e. Beukers' eta quotient $t(q)=(\eta_1\eta_6/\eta_2\eta_3)^{12}$, and
$F(t(q))=(\eta_2\eta_3)^7/(\eta_1\eta_6)^5$. The modular identification is used nowhere in the proof.)*

---

## 10. What to formalise, in order

```
Lemma 3.1  (4^k C(1/2,k) = ±2·Catalan)            [already in SQRT_APERY_FORMAL §3]
   └─▶ Lemma A: sqrt(1-34t+t^2) ∈ 1+tZ[[t]]        (3 lines; the identity (1-17t)^2-288t^2)
          └─▶ Cor A': F·sigma ∈ 1+tZ[[t]],  R=1/(F·sigma) ∈ 1+tZ[[t]],  r_n ∈ Z
Fact 4.1  (L_1 f = 0; symmetric-square identity)   [already in SQRT_APERY_FORMAL §2]
   └─▶ Lemma 4.2 (second solution f·ℓ + g of L_1)
   └─▶ Lemma 4.3 (Wronskian w = 1/sigma)
          └─▶ Prop 4.4  (L_Ap(F ℓ + G) = 0,  θ(G/F) = R - 1)
                 └─▶ Cor 4.5  q = t·exp(G/F),  log(q/t) = Σ r_n t^n / n
Lemma 5.1 (Dwork)  +  Lemma 5.2 (exp/log)  ──▶ Prop 5.3 (Dieudonné–Dwork criterion)
Cor A' + Cor 4.5 + Prop 5.3  ──▶  **Theorem N**:  Thm 4.2 ⟺ r_{mp^a} ≡ r_{mp^{a-1}} mod p^a
Lemma 2.1 (Lagrange inversion) + Lemma 2.2 (composition)  ──▶ Thm 4.2 (ii),(iii) from (i)
                                          ┌── Theorem KR  [EXTERNAL, p-adic, arXiv:0804.3049 Thm 2]
Lemma 7.1 (multinomial = binomial) ───────┤
Lemma 7.2 (B = G; creative telescoping) ──┘──▶ **Theorem 7.3 = Theorem 4.2**
```

Lean-facing shapes of the two remaining non-self-contained items:

* **Theorem KR** — may be taken as an axiom exactly as stated in §7.1 (a statement about explicit
  multivariate factorial sums), or formalised; its proof is Dwork-congruence combinatorics on
  factorials and harmonic numbers.
* **Lemma 7.2** — a single inhomogeneous three-term recurrence for
  $B_n=2\sum_j\binom nj^2\binom{n+j}j^2(H_{n+j}-H_{n-j})$, i.e. an identity between explicit finite
  sums of rationals; WZ-certifiable, same class as Apéry's recurrence.

If one prefers to keep a *single* axiom, take **Theorem N(2)** — the Gauss congruences for $r_n$ —
as the axiom: it is one clean arithmetic statement, and §§2–6 then derive all of Theorem 4.2 from it
with no further input.
