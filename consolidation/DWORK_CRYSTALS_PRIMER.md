# Dwork crystals: a primer for this project

*Fable (Opus 5), 2026-08-22. Written for River. Assumes modular forms and Apéry
sequences; assumes **no** prior contact with Dwork's $p$-adic theory. Everything
in §§1–3 is literature, cited with theorem numbers from the arXiv/journal
versions I actually read (sources downloaded to the session scratchpad and read
in full: arXiv:1903.11155v4, 1907.10390v2, 2105.14841v3, 2302.09603v2,
1306.5811, 1401.0854). §4 is our project. §5 is a reading order and three
exercises. Tags: **[thm]** = published theorem; **[ours]** = proved in this
project; **[obs]** = verified numerically here, not proved; **[heur]** =
heuristic.*

*New computations made while writing this note (all reproduced below, scripts in
the session scratchpad): $\mathrm{CT}(\Lambda^n)=A_n$ for $n\le5$; reflexivity
and the full facet list of $\Delta(\Lambda)$; $\mathscr C_2(\Lambda)$ and
$\mathscr C_2(\Lambda^2)$ exactly; the Dwork limit $\beta_{2^s}\sigma(\beta_{2^{s-1}})^{-1}$
against $F(t)/F(t^2)$ to $s=4$; row $C$'s unit roots at $p=5$ to $5^4$.*

---

## 1. The objects, concretely

### 1.1 A Laurent polynomial with a reflexive Newton polytope

Fix a Laurent polynomial $g(\mathbf x)\in\mathbf Z[x_1^{\pm1},\dots,x_n^{\pm1}]$,
write $\mathrm{supp}(g)\subset\mathbf Z^n$ for its exponent set and
$\Delta=\mathrm{conv}(\mathrm{supp}\,g)\subset\mathbf R^n$ for its **Newton
polytope**. $\Delta$ is *reflexive* if it has maximal dimension, contains
$\mathbf 0$ in its interior, and every codimension-one face $\tau$ is cut out by
$\ell_\tau(\mathbf u)=1$ for an *integral* linear functional $\ell_\tau$
(Beukers–Vlasenko, arXiv:2302.09603 §2). The hypothesis that actually appears in
the congruence theorems is weaker: **$\mathbf 0$ is the only interior lattice
point of $\Delta$.**

**The running example.** Straub's Laurent model of the Apéry $\zeta(3)$ numbers
(*Multivariate Apéry numbers and supercongruences of rational functions*, ANT 8
(2014) 1985–2008, Remark 1.4):

$$\Lambda(x,y,z)=\frac{(x+y)(z+1)(x+y+z)(y+z+1)}{xyz},\qquad
A_n=\sum_{k=0}^n\binom nk^2\binom{n+k}k^2 .$$

**Verified here** (exact expansion): $\mathrm{CT}(\Lambda^n)=A_n$ for $n\le5$,

| $n$ | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| $\mathrm{CT}(\Lambda^n)$ | $1$ | $5$ | $73$ | $1445$ | $33001$ | $819005$ |
| $A_n$ | $1$ | $5$ | $73$ | $1445$ | $33001$ | $819005$ |

and $\Delta=\Delta(\Lambda)\subset\mathbf R^3$ **is reflexive**: $11$ vertices,
$8$ facets with
$$\ell_\tau\in\{\pm e_1,\ \pm e_2,\ \pm e_3\ (\text{those occurring}),\ \pm(1,1,0),\ \pm(1,1,1)\},$$
exactly $21$ lattice points, of which $\mathbf 0$ is the unique interior one and
the other $20$ lie on the boundary; moreover $\mathrm{supp}(\Lambda)=\Delta\cap\mathbf Z^3$
(all $21$ lattice points carry a nonzero coefficient, e.g. the coefficient at
$\mathbf 0$ is $5=A_1$). This sharpens `NOME_INTEGRALITY.md` Fact 8.2, which
recorded the unique-interior-point property but not reflexivity.

### 1.2 The module of "regular functions", and the formal expansion

Put $f(\mathbf x)=1-t\,g(\mathbf x)$, a Laurent polynomial over
$R\subseteq\mathbf Z_p[[t]]$. BV's basic object is the $R$-module of rational
functions of controlled pole order and support,
$$\Omega_f(\mu)=\Bigl\{(m-1)!\,\frac{A(\mathbf x)}{f(\mathbf x)^m}\ \Big|\
m\ge1,\ \mathrm{supp}(A)\subset m\mu\Bigr\},$$
where $\mu\subseteq\Delta$ is **open** in the "face topology" (closed sets =
unions of faces); the smallest nonempty open set is the interior $\Delta^\circ$,
and $\mu_{\mathbf Z}:=\mu\cap\mathbf Z^n$. For a reflexive $\Delta$,
$\Delta^\circ_{\mathbf Z}=\{\mathbf 0\}$, so $h:=\#\mu_{\mathbf Z}=1$: **the
Apéry crystal has rank one.** The factorials $(m-1)!$ are not decoration — they
are what makes Griffiths–Dwork pole reduction work over a ring rather than a
field (BV-I Remark 2.1).

Multiplying by $\tfrac{dx_1}{x_1}\wedge\cdots\wedge\tfrac{dx_n}{x_n}$ turns
$\Omega_f$ into $n$-forms on $\mathbf T^n\setminus Z_f$, and
$\mathscr W_f:=\Omega_f/d\Omega_f$ is BV's **Dwork module**; by Batyrev
(Duke 69 (1993), Thm 7.13) it is $H^n_{\mathrm{dR}}(\mathbf T^n\setminus Z_f)$
when $f$ is $\Delta$-regular. BV assume no regularity.

**Periods.** Two expansions are used, and it is worth keeping them apart.

* **$\mathbf x$-adic** (BV-I §2): fix a vertex $\mathbf b$ of $\Delta$ whose
  coefficient is a unit, expand $g/f^m$ as a Laurent series supported in the cone
  $C(\Delta-\mathbf b)$. This embeds $\Omega_f\hookrightarrow\Omega_{\rm formal}$.
* **$t$-adic** (BV-II §2, Def. 2.1): when one coefficient of $f$ dominates
  $t$-adically — here the constant $1$, sitting at the *interior* point
  $\mathbf 0$, **not** at a vertex — expand in powers of $t$. The period map
  $p_{\mathbf 0}:\Omega_f\to\mathbf Z_p[[t]]$ takes the constant term of the
  $t$-adic Laurent expansion, and
  $$\boxed{\ p_{\mathbf 0}\!\left(\tfrac1f\right)=\sum_{k\ge0}\mathrm{CT}\bigl(g^k\bigr)\,t^k\ }$$
  which for $g=\Lambda$ is $F(t)=\sum_nA_nt^n$. This is why the mirror-family
  set-up needs BV-II and not just BV-I: for $f=1-tg$ *no vertex coefficient is a
  unit* in $\mathbf Z_p[[t]]$.

### 1.3 Why these are periods of $g=1/t$

$\Delta$ reflexive means the toric variety $\mathbf P_\Delta$ is Gorenstein
Fano, and the closure in $\mathbf P_\Delta$ of the affine hypersurface
$\{g(\mathbf x)=1/t\}\subset\mathbf T^n$ is (Batyrev) a family of Calabi–Yau
$(n-1)$-folds degenerating maximally at $t=0$ — the mirror/MUM point.
The generating function is literally a period of that family:
$$y_0(t)=\frac1{(2\pi i)^n}\oint\!\cdots\!\oint\frac{1}{1-t\,g(\mathbf x)}\,
\frac{dx_1}{x_1}\cdots\frac{dx_n}{x_n}=\sum_{k\ge0}\mathrm{CT}(g^k)\,t^k,$$
the integral of the toric residue form over the invariant real $n$-torus
(arXiv:2302.09603, eq. after Thm 1.4). So **"Apéry numbers = constant terms" is
the statement "the Apéry generating function is the holomorphic period at a MUM
point"**, and everything below is the arithmetic of that period.

---

## 2. Dwork's Cartier operator, and the crystal

### 2.1 The operator

On formal Laurent series, the **Cartier operator** is embarrassingly simple:
$$\mathscr C_p\Bigl(\sum_{\mathbf k}a_{\mathbf k}\mathbf x^{\mathbf k}\Bigr)
=\sum_{\mathbf k}a_{p\mathbf k}\mathbf x^{\mathbf k}. \tag{2.1}$$
This is Dwork's $\psi$ (*On the congruence properties of the zeta function of
algebraic varieties*, Crelle 203 (1960), §2) and Reich's $\Psi$; BV-I §3 is the
observation that this two-line operator, applied to *rational functions* rather
than to Dwork's exponential/Banach spaces, already carries the whole theory.
Equivalently, on rational functions (BV-II §2)
$$\mathscr C_p(\omega)(\mathbf x)=\frac1{p^n}\sum_{\mathbf y^p=\mathbf x}\omega(\mathbf y),$$
so $\mathscr C_p$ is a *trace*, i.e. the transpose of Frobenius $\mathbf x\mapsto\mathbf x^p$.

Three facts do all the work:

* **[thm]** $\mathscr C_p\circ\theta_i=p\,\theta_i\circ\mathscr C_p$ for
  $\theta_i=x_i\partial_{x_i}$; hence $\mathscr C_p(d\Omega)\subseteq p\,d\Omega$.
* **[thm]** $\mathscr C_p(g(\mathbf x^p)h)=g(\mathbf x)\mathscr C_p(h)$ (BV-I Lemma 3.2).
* **[thm]** (BV-I Prop. 3.3) For $p>2$, $\mathscr C_p(\Omega_f)\subseteq\widehat\Omega_{f^\sigma}$,
  where $\sigma$ is a Frobenius lift on $R$ and $f^\sigma$ is $f$ with $\sigma$
  applied to coefficients. The proof is the expansion
  $f^p=f^\sigma(\mathbf x^p)-pG(\mathbf x)$ plus the estimate
  $\mathrm{ord}_p(F_{\mathbf u,\mathbf v})\ge\frac{p-2}{p-1}(v_0-1)$, which is
  where $p>2$ enters and *only* where it enters.

> **Caveat we must keep flagged.** BV-I Prop. 3.3, and hence most of BV-I/III,
> assumes $p>2$. Mellit–Vlasenko's congruence theorem (§2.3) has **no**
> hypothesis on $p$. A large part of our census lives at $p=2$ (rows $E,F,\alpha,
> \varepsilon$, Zudilin), so we lean on the Mellit–Vlasenko half there.
> (`CRYSTAL_THEOREM_F.md` §6 item 14.)

### 2.2 Hasse–Witt, the unit-root quotient, "crystal"

**Definition** (BV-I Def. 3.5). For $m\ge1$,
$$(\beta_m(\mu))_{\mathbf u,\mathbf v\in\mu_{\mathbf Z}}
=\text{coefficient of }\mathbf x^{m\mathbf v-\mathbf u}\text{ in }f(\mathbf x)^{m-1},$$
$\beta_1=\mathrm{Id}$. $\beta_p(\mu)$ is the **Hasse–Witt matrix**. For our
rank-one case $\mu=\Delta^\circ$, $\mu_{\mathbf Z}=\{\mathbf 0\}$,
$$\beta_m(t)=\mathrm{CT}\bigl((1-tg)^{m-1}\bigr)=\sum_{k<m}(-1)^k\binom{m-1}kA_kt^k,
\qquad \gamma_m(t)=\sum_{k<m}A_kt^k,$$
and $\beta_p\equiv\gamma_p\equiv F_p(t)\pmod p$ (BV-II Lemma 4.3): the
Hasse–Witt "matrix" is the mod-$p$ truncation of the period. For Apéry at $p=2$:
$\beta_2=1-5t$, $\beta_4=1-15t+219t^2-1445t^3$, $\gamma_4=1+5t+73t^2+1445t^3$.

**Definition** (BV-I Def. 4.1). $U_f(\mu)=\{\omega:\mathscr C_p^s\omega\equiv0
\bmod p^s\ \forall s\ge1\}$. **[thm]** (BV-I Prop. 4.2) this equals
$\widehat\Omega_f(\mu)\cap d\Omega_{\rm formal}$ — the forms that "die on formal
expansion" in Katz's phrase. Note $d\Omega_f\subsetneq U_f$ in general: $U_f$ is
*formally* exact, not exact.

> **[thm] (BV-I Theorem 4.3 — the unit-root crystal.)** If $\beta_p(\mu)$ is
> invertible over $R$, then
> $$Q_f(\mu):=\widehat\Omega_f(\mu)/U_f(\mu)$$
> is a **free $R$-module of rank $h=\#\mu_{\mathbf Z}$**, with basis the images of
> $\omega_{\mathbf u}=\mathbf x^{\mathbf u}/f$, $\mathbf u\in\mu_{\mathbf Z}$.

By Remark 4.4, $\mathscr C_p:Q_f(\mu)\to Q_{f^\sigma}(\mu)$ is **invertible**
(its matrix is $\equiv\beta_p^{\mathrm t}\bmod p$), while $\mathscr C_p$ is
divisible by $p$ on $U_f(\mu)$. That dichotomy is the name: $Q_f$ is the
**unit-root quotient**, $U_f$ the largest subcrystal on which $\mathscr C_p$ is
$\equiv0\bmod p$. A **crystal** over $R$ (BV-I Def. 4.9) is just: a differential
$R$-module $M_f$ for every $f$, plus for every Frobenius lift $\sigma$ an
$R$-linear $\mathscr C_p:M_f\to M_{f^\sigma}$ commuting with the connection.
"Frobenius" here means $\mathscr C_p$ — the *inverse* Cartier, i.e. Frobenius on
the dual $Q_f^\vee$ (BV-I Remark 5.4). Ordinarity is *exactly* invertibility of
$\beta_p\bmod p$; there is no other hypothesis anywhere.

Write $\Lambda_\sigma$ for the (transposed) matrix of $\mathscr C_p$ on $Q_f(\mu)$
and $N_\delta$ for that of a derivation $\delta$ of $R$ (defined inside BV-I Thm 5.3):
these are the **Frobenius** and the **Gauss–Manin connection** of the unit-root
crystal.

### 2.3 Dwork's congruences

> **[thm] (BV-I Theorem 5.3.)** With $\beta_p(\mu)$ invertible, for all $m,s\ge1$
> $$\beta_{mp^s}(\mu)\equiv\Lambda_\sigma\,\beta^\sigma_{mp^{s-1}}(\mu)\pmod{p^s},
> \qquad\text{so}\qquad
> \Lambda_\sigma\equiv\beta_{p^s}(\mu)\,\beta^\sigma_{p^{s-1}}(\mu)^{-1}\pmod{p^s},$$
> and $\delta(\beta_{mp^s})\equiv N_\delta\beta_{mp^s}\pmod{p^s}$.

The proof is three lines once the framework is set: apply the **period map mod
$m$** $\tau_{m\mathbf v}(\omega)=\mathrm{CT}(f^{v_0}\mathbf x^{-\mathbf v}\omega)$
(BV-I Lemma 5.1: it kills derivatives mod $m$ and commutes with $\delta$ mod $m$)
to $\mathscr C_p(\omega_{\mathbf u})\equiv\sum\lambda_{\mathbf u\mathbf w}
\omega^\sigma_{\mathbf w}\bmod p\,U_{f^\sigma}$, and use
$\tau_{m\mathbf v}\equiv\tau^\sigma_{m\mathbf v/p}\circ\mathscr C_p\bmod
p^{\mathrm{ord}_p m}$ (Prop. 5.2). BV-I Theorem 5.7 is the same statement read
on $\mathbf x$-expansion coefficients, and recovers Katz's *Internal
reconstruction of unit-root $F$-crystals* (Ann. ENS 18 (1985)) Thm 6.2:
$\mathbf a_{p^s\mathbf k}\equiv\Lambda_\sigma\mathbf a^\sigma_{p^{s-1}\mathbf k}\bmod p^s$.

For the mirror family the $t$-adic version is what we want:

> **[thm] (BV-II Cor. 3.1 + Cor. 4.4; = Mellit–Vlasenko, IJNT 12 (2016), Thm 1;
> Samol–van Straten, Ann. Math. Québec 39 (2015) 185–203.)**
> Let $g\in\mathbf Z_p[x^{\pm}]$ have $\mathbf 0$ as the only interior lattice
> point of its Newton polytope, $b_n=\mathrm{CT}(g^n)$, $f(X)=\sum b_nX^n$,
> $f_s(X)=\sum_{n<p^s}b_nX^n$. Then $\Lambda_\sigma=f(t)/f(t^\sigma)$ and for all
> $s\ge1$
> $$\frac{f(X)}{f(X^p)}\equiv\frac{f_s(X)}{f_{s-1}(X^p)}\pmod{p^s\mathbf Z_p[[X]]},
> \quad\Longleftrightarrow\quad
> f_{s+1}(X)f_{s-1}(X^p)\equiv f_s(X)f_s(X^p)\pmod{p^s}.$$
> **No hypothesis on $p$.** BV-II extends this to
> $\gamma_m\equiv\Lambda_\sigma\sigma(\gamma_{m/p})\bmod p^{\mathrm{ord}_p m}$ for
> every $m$ and every Frobenius lift.

Reading the $t^{mp^s}$ coefficient gives the **Dwork/Gauss tower**
$$a_{mp^s}\equiv a_{mp^{s-1}}\pmod{p^s},$$
and $s=1$ gives the **Lucas congruence** $a_n\equiv\prod_ia_{n_i}\bmod p$ over
the base-$p$ digits of $n$ (Straub, Remark 1.4).

*Citation correction for our files.* `CRYSTAL_THEOREM_F.md` §1.3 and §6 attribute
the Gauss/Dwork tower to "BV-I Cor. 5.9". BV-I Cor. 5.9 is a different statement
(it concerns $g/f$ where **every** lattice point of $\Delta$ is a vertex, and
gives $a_{p^s\mathbf k}\equiv a_{p^{s-1}\mathbf k}$ for its $\mathbf x$-expansion
coefficients). The correct citations for our tower are **Mellit–Vlasenko Thm 1 /
Samol–van Straten**, or **BV-II Cor. 3.1 + Cor. 4.4**.

### 2.4 Supercongruences

Dwork's congruences are mod $p^s$. For Apéry one has mod $p^{3s}$:

* **[thm]** (Beukers, *Some congruences for the Apéry numbers*, JNT 21 (1985)
  141–155): for $p\ge5$ and $m,r\ge1$,
  $A(mp^r-1)\equiv A(mp^{r-1}-1)\pmod{p^{3r}}$.
* **[thm]** (Coster, adapting Beukers' proof; case $r=1$ earlier by Gessel and by
  Mimura; $m=r=1$ mod $p^2$ conjectured by Chowla–Cowles–Cowles):
  $A(mp^r)\equiv A(mp^{r-1})\pmod{p^{3r}}$ for $p\ge5$.
  The two are the *same* statement, by $A(-n)=A(n-1)$ under the $\Gamma$-extension
  of the binomials (Straub, Remark 1.2).

I re-verified all of this while writing: for $n<2000$, the Dwork tower holds with
$0$ failures at $p=2,3,5,7,11$ ($1990,996,496,330,198$ instances) and the mod
$p^{3r}$ supercongruences hold with $0$ failures at $p=5,7,11,13$ in both forms.

**BV's mechanism for supercongruences** is Dwork crystals III. Replace $U_f=\mathrm{fil}_1$
by
$$\mathrm{fil}_k=\{\omega\in\widehat\Omega_f:\mathscr C_p^s(\omega)\in p^{ks}\widehat\Omega_{f^{\sigma^s}}\ \forall s\ge1\},$$
the module of **$k$-th formal derivatives** — concretely, $\omega=\sum a_{\mathbf u}\mathbf x^{\mathbf u}$
lies in $\mathrm{fil}_k$ iff $\gcd(u_1,\dots,u_n)^k\mid a_{\mathbf u}$ (BV-III §3).

> **[thm] (BV-III Theorem 4.2 + Cor. 5.9.)** For $1\le k<p$, under the *$k$-th
> Hasse–Witt condition* (the image of $\mathscr C_p$ mod $p^k$ is as large as
> possible; equivalently certain determinants $hw^{(j)}$, $j\le k$, are units),
> $$\widehat\Omega_f(\mu)\cong\Omega_f^{(k)}(\mu)\oplus\mathrm{fil}_k(\mu),\qquad
> \mathscr C_p(\widehat\Omega_f)\subseteq\Omega^{(k)}_{f^\sigma}+p^k\mathrm{fil}^\sigma_k,$$
> with $\Omega^{(k)}_f$ free on $\mathbf x^{\mathbf u}/f^m$, $m\le k$.
>
> **[thm] (BV-III Prop. 5.11, "supercongruences".)** Consequently the expansion
> coefficients satisfy
> $$\mathbf a_{p^s\mathbf m}\equiv\Lambda\,\sigma(\mathbf a_{p^{s-1}\mathbf m})\pmod{p^{sk}},
> \qquad \delta(\mathbf a_{p^s\mathbf m})\equiv N_\delta\,\mathbf a_{p^s\mathbf m}\pmod{p^{sk}} .$$

Honesty: BV-III **did not** obtain a mod-$p^{2s}$ version of the
Mellit–Vlasenko congruence itself; that is their Conjecture 7.5, verified
numerically for the hypercubic/hyperoctahedral ($m=2$) and simplicial ($m=n+1$)
families. So the general mechanism behind $p^{3r}$ for Apéry is still open, and
Beukers/Coster remain the only proofs.

---

## 3. What the theory computes

### 3.1 Integrality of $a_n$, and of the nome

Two different integralities, with two different mechanisms.

**(a) $a_n\in\mathbf Z$ and $F(t)/F(t^p)\in\mathbf Z_p[[t]]^\times$.** Immediate
from the constant-term model plus Mellit–Vlasenko. This is the "$F$-half" of the
Frobenius structure and is *unconditional* for every row with a Laurent model.

**(b) The nome $q=t\exp(G/F)$ is integral.** This is the "log-half" and needs
more. Here is **Dwork's lemma**, the engine, with its proof (also
`NOME_INTEGRALITY.md` Lemma 5.1, Prop. 5.3).

> **Lemma (Dwork).** For a prime $p$ and $g\in1+x\mathbf Q_p[[x]]$,
> $$g\in1+x\mathbf Z_p[[x]]\iff h:=\frac{g(x^p)}{g(x)^p}\in1+px\mathbf Z_p[[x]].$$

*Proof.* ($\Rightarrow$) For $G\in\mathbf Z_p[[x]]$ one has $G(x)^p\equiv G(x^p)
\bmod p$ (expand; $a^p\equiv a$ and $p\mid\binom pj$ for $0<j<p$). So
$g(x^p)-g(x)^p\in p\mathbf Z_p[[x]]$, both sides have constant term $1$, and
dividing by the unit $g^p$ gives the claim.
($\Leftarrow$) Induct on $n$: suppose $g_0,\dots,g_{n-1}\in\mathbf Z_p$, put
$G=\sum_{i<n}g_ix^i$ and $T=g-G\in g_nx^n+x^{n+1}\mathbf Q_p[[x]]$. Since $2n>n$,
mod $x^{n+1}$ we get $g^p\equiv G^p+pG^{p-1}T$, so $[x^m]g^p=[x^m]G^p$ for $m<n$
and $[x^n]g^p=[x^n]G^p+pg_n$. Hence
$[x^n](g^ph)=[x^n]G^p+pg_n+\sum_{k\ge1}h_k[x^{n-k}]G^p$, the last sum in
$p\mathbf Z_p$. But $[x^n](g^ph)=[x^n]g(x^p)$, which is $g_{n/p}\in\mathbf Z_p$
or $0$; and $[x^n]G^p\equiv[x^n]G(x^p)=[x^n]g(x^p)\bmod p$. Combining,
$pg_n\in p\mathbf Z_p$. $\square$

Applied to $g=q/t=\exp(\sum r_nt^n/n)$ this becomes the **Dieudonné–Dwork
criterion**: $q/t\in1+t\mathbf Z_p[[t]]$ iff $r_n\equiv r_{n/p}\bmod p^{v_p(n)}$
— the classical *Gauss congruences*. BV-III use the same lemma in the form
"$g(t)-\tfrac1pg(t^\sigma)\in\mathbf Z_p[[t]]\iff\exp g\in\mathbf Z_p[[t]]$"
(BV-III Lemma 7.10) to prove:

> **[thm] (BV-III Cor. 7.11.)** For a *completely symmetric* Calabi–Yau family
> ($\mathcal G$ transitive on the vertices of $\Delta$, all non-constant
> coefficients of $g$ equal to a unit $\gamma$, $p$ odd, $p\nmid\#\mathcal G$),
> the canonical coordinate $q(t)=t\exp(G(t)/F(t))$ lies in $t\mathbf Z_p[[t]]$;
> equivalently the mirror map is $p$-integral.

Note (as `NOME_INTEGRALITY.md` §8 does) that **Straub's $\Lambda$ does not
satisfy BV-III's complete-symmetry hypothesis**, so Cor. 7.11 does not apply to
Apéry off the shelf; our proof of nome integrality goes through
Krattenthaler–Rivoal instead.

### 3.2 The unit root and the Hasse invariant

Specialise $t\mapsto t_0\in\mathbf Z_p$ with $\gamma_p(t_0)\in\mathbf Z_p^\times$
(i.e. the Hasse invariant $\beta_p\equiv F_p$ does not vanish: the fibre is
*ordinary*). Then

> **[thm] (BV-II Remark 4.5; BV-I Appendix A.)**
> $\lim_{s\to\infty}\gamma_{p^s}(t_0)/\gamma_{p^{s-1}}(t_0)$ exists and equals the
> **unit root** of the zeta function of the fibre $\{g=1/t_0\}$ over $\mathbf F_p$.

(BV-I Thm A.1: $\mathrm{Tr}\bigl((q^s-1)^n\mathscr C_q^s\mid\widehat\Omega_f\bigr)
=\#(\mathbf T^n\setminus Z_f)(\mathbf F_{q^s})$; hence
$\det(1-T\mathscr C_q\mid Q_f)/(1-T)$ is, up to sign, the unit-root part of the
zeta function of $Z_f$.) When $\beta_p(t_0)\equiv0$ the residue disc is
**supersingular** and there is no unit root at all. Exercise 3 below does this
for row $C$ at $p=5$.

### 3.3 The excellent Frobenius lift

The Frobenius lift $\sigma$ on $R$ is a *choice*. BV-III's discovery:

> **[thm] (BV-III Thm 7.3, Def. 7.14, Thm 7.15.)** For a completely symmetric CY
> family with $hw^{(1)},hw^{(2)}$ invertible, there is a **unique** Frobenius lift
> with
> $$\mathscr C_p\Bigl(\frac1f\Bigr)\equiv\frac{F(t)}{F(t^\sigma)}\cdot\frac1{f^\sigma}
> \pmod{p^2\mathrm{fil}_2^\sigma}$$
> — i.e. the off-diagonal entry $\lambda_1$ vanishes, one level deeper than
> BV-I gives. It is characterised in the canonical coordinate by
> $$q^\sigma=\gamma^{p-1}q^p ,$$
> and $t^\sigma$ lies in the $p$-adic completion of
> $\mathbf Z_p[t,1/hw^{(1)},1/hw^{(2)}]$.

The identity that makes this transparent is BV-III Lemma 7.13,
$\lambda_1(t)=\frac{F(t)F(t^\sigma)}{W(t^\sigma)}\log\bigl(\gamma^{p-1}q^p/q^\sigma\bigr)$
($W$ = Wronskian). Dwork had this for the Legendre family (*p-adic cycles*,
Thm 8.1); BV's point is that it is not special to elliptic curves.

### 3.4 The Frobenius matrix at the MUM point, and why $\zeta_p(3)$

Now the object that our Theorem F is a rank-one extension of. Let
$L=\theta^n+a_1\theta^{n-1}+\dots+a_n$ with $\theta=t\frac{d}{dt}$ and $a_j(0)=0$
— a **MUM-type operator**. It has a unique *standard basis* of solutions
$$y_i=F_0\frac{\log^it}{i!}+F_1\frac{\log^{i-1}t}{(i-1)!}+\dots+F_i,\qquad i=0,\dots,n-1,$$
with $F_0(0)=1$, $F_i(0)=0$ for $i>0$. A **$p$-adic Frobenius structure over $R$**
(arXiv:2302.09603 Def. 1.1) is $\mathscr A=\sum_{j<n}A_j(t)\theta^j\in R[\theta]$
with $A_0(0)=1$ such that $\mathscr A(y_i(t^p))$ solves $L$ again.

In plain terms: **$\mathscr A$ is the matrix expressing "Frobenius-pulled-back
solutions" in terms of solutions**, in the one basis that the MUM degeneration
canonically supplies. It is the $p$-adic analogue of the monodromy matrix at
$t=0$, and it is *the* thing a crystal has that a plain differential equation
does not.

> **[thm] (arXiv:2302.09603 Prop. 1.3.)** There are $\alpha_0=1,\alpha_1,\dots,
> \alpha_{n-1}\in\mathbf Z_p$, namely $\alpha_j=A_j(0)$, with
> $$\mathscr A\bigl(y_i(t^p)\bigr)=p^i\sum_{j=0}^i\alpha_j\,y_{i-j}(t),\qquad i=0,\dots,n-1 .$$

So in the standard basis the Frobenius at the MUM point is **upper triangular with
diagonal $(1,p,p^2,\dots,p^{n-1})$**, and unipotent-times-diagonal; the whole
content sits in $(\alpha_1,\dots,\alpha_{n-1})$. The diagonal is forced by the
degeneration (a Tate/nodal fibre): each $\log t$ costs a factor $p$ because
$\log(t^p)=p\log t$.

> **[thm] (Thm 1.4, simplicial family $g=x_1+\dots+x_n+1/(x_1\cdots x_n)$, $p>n+1$.)**
> $\alpha_j=[x^j]\dfrac{\Gamma_p(x)}{\Gamma_p(x/(n+1))^{n+1}}$.
>
> **[thm] (Thm 1.5, hyperoctahedral family $g=\sum(x_i+1/x_i)$, $p>n$, $L$ irreducible.)**
> $\alpha_j=[x^j]\exp\Bigl(-\sum_{m\ge2}\frac{\zeta_p(m)}mx^m\Bigr)$, so
> $$\alpha_1=\alpha_2=0,\quad\alpha_3=-\tfrac13\zeta_p(3),\quad\alpha_4=0,\quad
> \alpha_5=-\tfrac15\zeta_p(5),\quad\alpha_6=\tfrac1{18}\zeta_p(3)^2,\ \dots$$

**Why $\zeta_p(3)$ appears, in one paragraph.** The $\alpha_j$ are computed from
the $p$-adic $\widehat\Gamma$-class of the family: they are the Taylor
coefficients of a ratio of Morita $\Gamma_p$'s built from the exponent data of
$g$. And $\log\Gamma_p$ has only **odd** zeta values,
$$\log\Gamma_p(x)=\Gamma_p'(0)\,x-\sum_{m\ge2}\frac{\zeta_p(m)}mx^m,
\qquad \zeta_p(m):=L_p\bigl(m,\omega^{1-m}\bigr),$$
and $\zeta_p(m)=0$ for even $m$ because $\omega^{1-m}$ is then an *odd*
character and Kubota–Leopoldt $L_p(\cdot,\chi)\equiv0$ for odd $\chi$
(arXiv:2302.09603 (1.6), citing Schikhof §§58,61; independently verified in
`MUM_SURVEY.md` §3 at $p=2,3,5,7,11,13,17$ to 33–83 digits). This is the exact
$p$-adic mirror of the complex monodromy conjecture, where the standard-basis
monodromy entries of a rank-$n$ CY operator are conjectured to lie in the ring
generated by $\zeta(k)/(2\pi i)^k$, $2\le k<n$ (van Straten). For $n=4$ this is
the Candelas–de la Ossa–van Straten conjecture $\alpha_1=\alpha_2=0$,
$\alpha_3\in\mathbf Q\,\zeta_p(3)$ (arXiv:2104.07816 §4.4), proved for the
quintic by Shapiro (IMRN 2009; JNT 2012).

---

## 4. Where our theory sits inside this

### 4.1 The companion row $b_n$ is the period of an *extension*

Every result above is about the **pure** crystal $\mathcal M$ — the Picard–Fuchs
module of the family. Our Apéry rows have a second sequence. For an (R2)/(R3)
row, $B(t)=\sum b_nt^n$ satisfies not $LB=0$ but
$$L\,B=\rho\ne0$$
(the defect at $n=0$). **[ours]** So $B$ is the period of an extension in the
category of $F$-isocrystals / mixed motives
$$0\longrightarrow\mathcal M\longrightarrow\widetilde{\mathcal E}\longrightarrow\mathbf Q_p(0)\longrightarrow0,$$
with class $\rho$ (`CRYSTAL_THEOREM_F.md` §1.5, `THEORY_NOTES_06` §1). Because
$B$ carries **no logarithm**, Frobenius in the basis $(y_0,\dots,y_{r-1},B)$ has
the shape
$$\widetilde\Phi=\begin{pmatrix}\Phi_{\rm pure}&\widetilde\xi\,e_1\\0&\varepsilon\end{pmatrix},
\qquad \Phi(B)=\varepsilon B+\widetilde\xi\,y_0,$$
the $y_1,\dots,y_{r-1}$ components vanishing by log-freeness. **Splitting
criterion [ours, elementary]:** $\widetilde\xi$ can be killed by
$B\mapsto B-cy_0$ iff $\varepsilon\ne1$ (solve $c(\varepsilon-1)=\widetilde\xi$).

> The extension splits $p$-adically **iff** the Frobenius eigenvalue on the
> trivial quotient differs from the unit root $1$ of $\mathcal M$. When they
> coincide the Frobenius is unipotent and $\widetilde\xi$ is a genuine invariant:
> the $p$-adic period of the extension.

This is the same shape, one step up, as BV-III's $\lambda_1$ (§3.3): there the
off-diagonal entry vanishes for the excellent lift; here it is an *obstruction*
that no choice of lift removes.

### 4.2 The slope law and the trichotomy as Frobenius statements

**[ours, proved]** Define
$$\mu_1=\lim_n\frac{v_p(a_n)}n,\quad
\mu_{\det}=\lim_n\frac{v_p(a_{n+1}b_n-a_nb_{n+1})}n,\quad
\sigma_p:=\mu_{\det}-2\mu_1 .$$
$\mu_1$ is the unit-root slope, $\mu_{\det}$ the slope of the determinant
(Wronskian) of the crystal, $\sigma_p$ the **gap between the two Frobenius
slopes**. From $b_{n+1}/a_{n+1}-b_n/a_n=-\mathrm{Cas}_n/(a_na_{n+1})$:
$$v_p\bigl(\xi_p-b_n/a_n\bigr)=n\sigma_p+O(\log n),$$
so $\xi_p=\lim b_n/a_n$ exists **iff $\sigma_p>0$**; and for the ten (R2)/(R3)
rows $\sigma_p=v_p(c)$ exactly (`CRYSTAL_THEOREM_F.md` §3, verified in $84$
cells).

**[obs]** The tower trichotomy, in Frobenius language. With
$\rho^A_s=a_{ap^{s+1}}/a_{ap^s}$, $\rho^B_s=b_{ap^{s+1}}/b_{ap^s}$:
$$\lim_s\rho^A_s=1\ \text{always}\quad[\text{{\bf thm}: this \emph{is} Dwork's congruence}],$$
$$\lim_s\rho^B_s=\begin{cases}
\varepsilon(p)\,p^{-(w+1)}, & \sigma_p=0\ \text{(good prime; distinct eigenvalues, extension splits)},\\
1, & \sigma_p>0,\ \xi_p\ne0\ \text{(unipotent; the obstruction is the invariant)},\\
0, & \sigma_p>0,\ \xi_p=0 .\end{cases}$$
Here $w+1$ is the **Hodge depth** — the denominator exponent
($d_n^{w+1}b_n\in\mathbf Z$), i.e. the Tate twist in
$0\to\mathcal V(w+1)\to\mathcal E\to\mathbf Q(0)\to0$ — **not** the rank of $L$
and not the modular weight. Census: $164$ cells, zero exceptions
(`GOOD_PRIME_TOWERS.md`, `CRYSTAL_THEOREM_F.md` §2.2, §5.5).

Note that $\lim_s\rho^A_s=1$ **refutes** the natural guess that the tower sees a
Dwork unit root of a fibre: the tower runs into the cusp $t=0$, where the fibre
is nodal and its unit root is $1$. The genuine unit roots (§3.2, Exercise 3) are
$t$-dependent and often $\not\equiv1$, and sometimes do not exist.

**[heur, well-supported]** BV-III Def. 7.14 normalises the excellent lift by
$q^\sigma=\gamma^{p-1}q^p$. This predicts $\varepsilon(p)=\gamma^{p-1}$ with
$\gamma$ the constant relating the row's coordinate to the canonical one — which
is consistent with everything measured ($\varepsilon\equiv1$ exactly for rows with
rational mirror map, $\varepsilon=\chi_{-4}$ for the Catalan row). Untested
directly.

### 4.3 Theorem F's Kubota–Leopoldt value, and the shape mismatch

**Theorem F** identifies the unipotent obstruction with a $p$-adic $L$-value:
$$\xi_p=-\tfrac12\,Q(w+1)\,L_p\bigl(w+1,\ \psi\omega^{-w}\bigr)$$
for a modular row with Eisenstein source $\Phi=P(V)E^{\psi,\varphi}_{w+2}$
(14 identifications to ~3000 digits; `THEORY_NOTES_04` §1). In crystal language
$\xi_p=\widetilde\xi$, the off-diagonal Frobenius entry of
$\widetilde{\mathcal E}$: **crystal theory gives existence and rate, the value
stays modular/syntomic.** No paper in the BV series computes an extension
period.

**The shape mismatch, and why it is not a contradiction.** BV's $\alpha_j$ involve
only *odd* $\zeta_p(m)=L_p(m,\omega^{1-m})$, while our $\xi_p$ can be an
even-index value — Zudilin's Catalan row converges $2$-adically to
$\zeta_2(2)$. Two things are going on, and they should not be confused.

1. **A naming collision.** "$\zeta_p(2)$" means two different numbers.
   BV/CDvS write $\zeta_p(m)=L_p(m,\omega^{1-m})$; for even $m$ the character
   $\omega^{1-m}$ is odd, so this vanishes identically. Calegari's and Beukers'
   "$\zeta_2(2)$" is the Kubota–Leopoldt $\zeta_p(s)=L_p(s,\chi_0)$ at $s=2$,
   whose character is *even*, so it does not vanish. Beukers spells this out
   himself (*Irrationality of some $p$-adic $L$-values*, Acta Math. Sinica 24
   (2008), Introduction): Calegari's result is about the Kubota–Leopoldt
   $\zeta_2(2)$, **not** about $L_2(2,\chi_{-4})$, which vanishes because
   $\chi_{-4}$ is odd. So "$\zeta_p$ has no even values" and "$\zeta_2(2)\ne0$"
   are both true, of different objects.
2. **A structural difference.** BV's $\alpha_j$ are entries of the Frobenius of
   the **pure** crystal in the standard basis; their character is forced to be the
   diagonal Teichmüller twist $\omega^{1-m}$, hence the odd-only parity, mirroring
   the odd $\zeta(k)$ of the complex monodromy conjecture. Our $\xi_p$ is a period
   of an **extension**, and the character it carries is the character $\psi$ of the
   Eisenstein class the extension realises; the non-vanishing condition is that
   $\psi\omega^{-w}$ be even, which happens for both parities of $w+1$.
   **Even-index $L_p$-values are exactly what the pure crystal cannot see.**
   (`CRYSTAL_THEOREM_F.md` §6 item 11.) This is a real gap in the literature, not
   a bookkeeping error on either side.

### 4.4 Zudilin's non-modular row has a crystal

`SOURCES_S18_ZUDILIN.md` §5 concluded, correctly, that Zudilin's Catalan row has
no *modular* parametrisation. **[obs]** That does not survive rescaling: with
$R_m:=16^mQ_m$,
$$R_0,\dots,R_6=1,\ 28,\ 2596,\ 311536,\ 41759524,\ 5953698928,\ 882922368784,$$
one finds $R_m\in\mathbf Z$ for $m\le1200$ and the **full Dwork tower**
$R_{mp^s}\equiv R_{mp^{s-1}}\bmod p^s$ at $p=2,3,5,7$ (2288 instances, 0
failures). By Mellit–Vlasenko this is exactly the observable signature of a
constant-term sequence of a Laurent polynomial with $\mathbf 0$ the unique
interior lattice point. Slopes: $\mu_1=-4$ exactly, $\mu_{\det}=0$, so
$\sigma_2=8>0$ and $\sigma_p=0$ for $p=3,5,7,11,13$ — precisely the observed
good/slope split, and the existence of $\xi_2^{\rm Zud}$ without any modular
input. **Open, and well posed: find $\Lambda$ with $\mathrm{CT}(\Lambda^m)=16^mQ_m$.**
Brown–Zudilin's $\zeta(5)$ row likewise is a Dwork row with $\sigma_p=0$ at every
tested prime — a *structural* obstruction to the two-row method through it.

### 4.5 The quintic's $p$-adic $\widehat\Gamma$-class is the same structure at rank 4

**[ours, verified + lit]** The quintic has no Apéry limit at all (its recurrence
has order 1); its $\zeta(3)$ lives in the connection coefficients between the MUM
Frobenius basis and the conifold, i.e. in the $\widehat\Gamma$-class:
$$V=\frac{\sqrt5}{4\pi^2}\bigl(\hat y_3-10\zeta(2)\hat y_1-40\zeta(3)\hat y_0\bigr),
\qquad -40=\chi/H^3 ,$$
verified to 116 decimal digits (`MUM_SURVEY.md` §2). The $p$-adic version is the
same expression with $\Gamma\to\Gamma_p$:
$$[\rho^3]\log\frac{\Gamma_p(1+5\rho)}{\Gamma_p(1+\rho)^5}=-40\,\zeta_p(3)
\qquad\text{at every }p,$$
which is exactly BV's Thm 1.4 at $n=4$
($\Gamma_p(x)/\Gamma_p(x/5)^5=1-\tfrac8{25}\zeta_p(3)x^3+\cdots$, and
$-8/25=-40/5^3$ after $\varphi=t^5$), proved by Shapiro. So: **the $\alpha_j$ of
§3.4 at rank 4 are the $p$-adic $\widehat\Gamma$-class, and $\chi/H^3$ is the same
rational number at every place.** Our rank-2/3 Apéry $\xi_p$ is the *extension*
analogue of that; the quintic shows what the pure part looks like when the row
degenerates to nothing.

---

## 5. Reading order, and three exercises

### 5.1 Read in this order

1. **Mellit–Vlasenko**, arXiv:1306.5811 (9 pages). Theorem 1 is the whole
   congruence, proof entirely elementary. Read it first: it is the statement our
   project actually uses, at every prime including $2$.
2. **BV, Dwork crystals I**, arXiv:1903.11155. §§2–3 (the module, the Cartier
   operator, Prop. 3.3), §4 (Thm 4.3, the crystal), §5 (Thm 5.3, Thm 5.7,
   Example 5.5 = the Legendre family). Skip §§6–7 and the appendix on first pass.
3. **BV, Dwork crystals II**, arXiv:1907.10390, §§2–4. This is where $f=1-tg$ and
   the $t$-adic period map live, and where the MV theorem is re-proved
   (Cor. 3.1, Cor. 4.4). Short.
4. **BV, arXiv:2302.09603**, §§1–2 only. Prop. 1.3 (the triangular Frobenius at
   the MUM point) and Thms 1.4–1.5 are what §4 of this primer is an extension of.
5. **BV, Dwork crystals III**, arXiv:2105.14841, §§1, 3, 7. The $\mathrm{fil}_k$
   filtration, the excellent lift (Def. 7.14), the mirror-map integrality
   (Cor. 7.11), and honest Conjecture 7.5. Heavy; skim §§4–5.
6. **Dwork, *p-adic cycles*** (Publ. IHÉS 37 (1969)): §1 Cor. 2, §2 Thm 2 (the
   hypergeometric congruences), (6.29) (the unit root as a limit), Thm 8.1
   (the excellent lift for the Legendre family). Read *after* BV — BV is the
   Rosetta stone.
7. **Beukers**, JNT 21 (1985) and **Coster** — for the supercongruences; and
   **Samol–van Straten**, arXiv:0911.0797, for the reflexive-polytope statement.
8. **Straub**, arXiv:1401.0854, Remark 1.4 and Remark 1.2 — one page, and it is
   where our $\Lambda$ comes from.

### 5.2 Three exercises

**Exercise 1 (the Cartier operator, by hand, at $p=2$).**
Compute $\mathscr C_2(\Lambda^n)$ for $n\le2$ directly from (2.1) — i.e. keep the
monomials of $\Lambda^n$ (as a Laurent polynomial in $x,y,z$) whose exponent
vectors are *all even*, and halve them.

*Answers* (verified exactly here). $\mathscr C_2(1)=1$. For $n=1$: the only
exponent vector of $\Lambda$ with all coordinates even is $\mathbf 0$, so
$$\mathscr C_2(\Lambda)=5\quad(\text{a constant!}).$$
For $n=2$: $\mathscr C_2(\Lambda^2)$ is a Laurent polynomial supported on
*exactly* the $21$ lattice points of $\Delta$, with constant term $73=A_2$:
coefficients $(1,6,1,1,19,6,1,1,1,6,1,6,\mathbf{73},19,6,6,1,6,1,1,1)$ against
$\Lambda$'s $(1,2,1,1,3,2,1,1,1,2,1,2,\mathbf 5,3,2,2,1,2,1,1,1)$ in the same
order. Now reduce mod $2$: **they agree.** That is the lesson:
$$\Lambda^p\equiv\Lambda^{(p)}(\mathbf x^p)\ (\mathrm{mod}\ p)
\ \Longrightarrow\ \mathscr C_p(\Lambda^p)\equiv\Lambda\ (\mathrm{mod}\ p),$$
i.e. $\mathscr C_p$ is a **lift of the inverse Cartier / a Frobenius descent**.
Everything in BV-I §3 is bookkeeping around this one congruence.
*(Caveat: BV-I Prop. 3.3 needs $p>2$; the operator (2.1) itself is defined at
every $p$, and it is the $p=2$ case that our census needs.)*

Follow-up: check
$\Lambda_\sigma\equiv\beta_{p^s}\sigma(\beta_{p^{s-1}})^{-1}\bmod p^s$ against
$F(t)/F(t^2)$. I did this to $s=4$: with $\beta_2=1-5t$, $\beta_4=1-15t+219t^2-1445t^3$,
$\beta_8,\beta_{16}$, the ratios agree with
$F(t)/F(t^2)=1+5t+4t^2+12t^3+12t^4+20t^5+\cdots$ modulo $2,4,8,16$ respectively,
on the nose.

**Exercise 2 (Dwork mod $p$ from Lucas).**
Show: (i) Lucas for constant terms. Over $\mathbf F_p$, $\Lambda(\mathbf x)^p
=\Lambda(\mathbf x^p)$, so writing $n=n_0+n_1p+\dots+n_\ell p^\ell$ with
$0\le n_i<p$,
$$\Lambda^n\equiv\Lambda(\mathbf x)^{n_0}\Lambda(\mathbf x^p)^{n_1}\cdots
\Lambda(\mathbf x^{p^\ell})^{n_\ell}\pmod p .$$
Because the exponent lattices $p^i\mathbf Z^d$ are nested and the supports are
bounded by $n_i\Delta$ with $n_i<p$, the constant term of the product factorises,
giving
$$A_n\equiv A_{n_0}A_{n_1}\cdots A_{n_\ell}\pmod p .$$
(ii) Deduce Dwork mod $p$: the base-$p$ digits of $mp$ are those of $m$ shifted,
with a new digit $0$ at the bottom, and $A_0=1$; hence
$A_{mp}\equiv A_0\cdot\prod A_{m_i}\equiv A_m\bmod p$. That is the $s=1$ case of
the tower. Then convince yourself why the $s\ge2$ case is genuinely deeper:
it needs the *quotient* $\gamma_{p^s}(t)/\gamma_{p^{s-1}}(t^p)$ to stabilise, i.e.
it needs the unit-root crystal to be free of rank $1$, not just a mod-$p$
factorisation. (Straub Remark 1.4 makes the same point.)

**Exercise 3 (the unit root of row $C$ at $p=5$).**
Row $C$ is Zagier's $(a,b,c)=(10,3,9)$ row, $a_n=\sum_k\binom nk^2\binom{2k}k$,
$$a_0,\dots,a_6=1,\ 3,\ 15,\ 93,\ 639,\ 4653,\ 35169 .$$
The Hasse invariant is $\beta_5\equiv\gamma_5(t)=\sum_{k<5}a_kt^k
\equiv1+3t+0\cdot t^2+3t^3+4t^4\pmod5$. Evaluate at $t_0=1,2,3,4$ and compute
$\lim_s\gamma_{5^s}(t_0)/\gamma_{5^{s-1}}(t_0)$.

*Answers* (computed here to $5^4$, i.e. $s\le4$, $5^4=625$ terms):

| $t_0$ | $\gamma_5(t_0)\bmod5$ | $\gamma_{5^s}/\gamma_{5^{s-1}}\bmod 5^s$, $s=1,2,3,4$ | unit root |
|---|---|---|---|
| $1$ | $1$ | $1,\ 1,\ 1,\ 1$ | $\lambda_0=1$ |
| $2$ | $0$ | — | **supersingular disc**, no unit root |
| $3$ | $0$ | — | **supersingular disc**, no unit root |
| $4$ | $4$ | $4,\ 24,\ 124,\ 624$ | $\lambda_0=-1$ |

So two of the four ordinary-candidate residue discs at $p=5$ are supersingular,
and on the other two the unit root is exactly $\pm1$. Contrast this with the
*tower* unit root $\lim_s a_{ap^{s+1}}/a_{ap^s}=1$, which is identically $1$ for
every row, prime and base: the tower is not sitting at a Teichmüller point of any
residue disc, it is running into the cusp $t=0$ where the fibre is nodal. That
contrast is, in my view, the single most instructive computation in this file.

---

## 6. Bibliography, with the numbers used

* F. Beukers, M. Vlasenko, **Dwork crystals I**, arXiv:1903.11155, IMRN 2021
  no. 12, 8807–8844. Rmk 2.1 (factorials), Lemma 3.1, Lemma 3.2, **Prop. 3.3**
  ($p>2$), Prop. 3.4 (face topology), **Def. 3.5** (Hasse–Witt), Prop. 3.6,
  Def. 4.1, Prop. 4.2, **Thm 4.3**, Rmk 4.4, Prop. 4.5, Def. 4.9 (crystal),
  Lemma 5.1, Prop. 5.2, **Thm 5.3**, Ex. 5.5 (Legendre), Cor. 5.6 (formal
  groups), **Thm 5.7** (= Katz Thm 6.2), Cor. 5.9, Thm 6.1 (weight filtration),
  Thm A.1 (point counting).
* — **Dwork crystals II**, arXiv:1907.10390, IMRN 2021, 4427–4444. Def. 2.1
  ($t$-adic period map), Prop. 2.2, **Thm 2.3**, **Cor. 3.1**
  ($\Lambda_\sigma=q(t)/q(t^p)$), **Thm 3.2** (= Mellit–Vlasenko), Rmk 3.3,
  Prop. 4.1, Thm 4.2, Lemma 4.3, **Cor. 4.4**, Rmk 4.5 (unit root as a limit),
  Thm 5.3 / Cor. 5.4 (A-hypergeometric).
* — **Dwork crystals III: from excellent Frobenius lifts towards
  supercongruences**, arXiv:2105.14841, IMRN 2023, 20433–20483. Def. 3.1
  ($\mathrm{fil}_k$), Lemma 3.3, **Thm 4.2**, Cor. 4.3, Def. 5.1, Prop. 5.7, Cor. 5.9,
  **Prop. 5.11** (supercongruences), Prop. 5.12, Def. 7.1, **Thm 7.3**
  (excellent lift), Conj. 7.5, Prop. 7.7, Lemma 7.10 (Dieudonné–Dwork),
  **Cor. 7.11** (mirror map $p$-integral), Lemma 7.13, **Def. 7.14**, Thm 7.15.
* — **Frobenius structure and $p$-adic zeta values**, arXiv:2302.09603,
  Adv. Math. 480 (2025). Def. 1.1, **Prop. 1.3**, **Thms 1.4–1.5**, Def. 2.1
  (cyclic MUM), Props 2.2–2.3.
* — *On $p$-integrality of instanton numbers*, Pure Appl. Math. Q. 19 (2023),
  7–44 (cited as [IN] by 2302.09603; Prop. 4.2, Thm 3.3, Props 5.1, 6.1).
* A. Mellit, M. Vlasenko, *Dwork's congruences for the constant terms of powers
  of a Laurent polynomial*, arXiv:1306.5811, IJNT 12 (2016) 313–321, **Thm 1**.
* K. Samol, D. van Straten, *Dwork congruences and reflexive polytopes*,
  arXiv:0911.0797, Ann. Math. Québec 39 (2015) 185–203.
* A. Straub, *Multivariate Apéry numbers and supercongruences of rational
  functions*, arXiv:1401.0854, ANT 8 (2014) 1985–2008. **Remark 1.4** (our
  $\Lambda$; Dwork and Lucas congruences), Remark 1.2 ($A(-n)=A(n-1)$),
  Thm 1.1, Thm 1.3.
* F. Beukers, *Some congruences for the Apéry numbers*, JNT 21 (1985) 141–155
  ($A(mp^r-1)\equiv A(mp^{r-1}-1)\bmod p^{3r}$, $p\ge5$); M. Coster (thesis)
  for $A(mp^r)\equiv A(mp^{r-1})\bmod p^{3r}$.
* F. Beukers, *Irrationality of some $p$-adic $L$-values*, Acta Math. Sinica 24
  (2008) 663–686 — Introduction (the $\zeta_2(2)$ vs $L_2(2,\chi_{-4})$ warning),
  Props 5.1, 6.1.
* B. Dwork, *p-adic cycles*, Publ. IHÉS 37 (1969) 27–115; *On the congruence
  properties of the zeta function of algebraic varieties*, Crelle 203 (1960)
  130–142 (the operator $\psi$); *A deformation theory for the zeta function of a
  hypersurface*, ICM Stockholm 1962 (1963), 247–259 (congruence (12)).
* N. Katz, *Internal reconstruction of unit-root $F$-crystals via expansion
  coefficients*, Ann. ENS 18 (1985) 245–285, Thm 6.2.
* V. Batyrev, Duke Math. J. 69 (1993), Thm 7.13.
* C. Krattenthaler, T. Rivoal, arXiv:0804.3049, Thm 2 (used in
  `NOME_INTEGRALITY.md` for the Apéry nome).
* Ph. Candelas, X. de la Ossa, D. van Straten, arXiv:2104.07816, §4.4 and
  App. B; I. Shapiro, IMRN 2009 no. 13, 2519–2545 and JNT 132 (2012) 1770–1779.

**Project files this primer connects to.** `CRYSTAL_THEOREM_F.md` (Theorem F$'$,
the slope law, the trichotomy, Zudilin's crystal), `NOME_INTEGRALITY.md`
(Dwork's lemma, Gauss congruences, Straub's $\Lambda$), `GOOD_PRIME_TOWERS.md`
(the eigenvalue census), `MUM_SURVEY.md` (the $\widehat\Gamma$-class at rank 4),
`THEORY_NOTES_04/06` (the extension-class picture).
