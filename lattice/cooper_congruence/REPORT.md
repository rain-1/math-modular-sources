# Cooper's magnetic congruence: reduction, mechanism, and two unconditional cells

*Working note, 2026-09-02.  All computations in PARI/GP 2.15.4, exact over $\mathbf Q$
unless a precision is quoted.  Scripts `lib.gp`, `01_beta.gp` … `10_exp.gp` with logs
`0*.log`, `10_exp.log`; data `beta_s*.txt`, `gamma_s*.txt`.  Background:
`lattice/cooper_sources/REPORT.md`, `consolidation/COMPANION_ARITHMETIC.md` §§4.5, 5,
`paper/companions/main.tex` Thm `cooper`.  Nothing outside this directory was modified.*

Claims are tagged **[proved]**, **[verified, range]** (exact rational computation over the
stated range), **[num, digits]**, or **[conjectural]**.

---

## 0. Verdict table

| # | statement | status |
|---|---|---|
| V1 | $\Phi\|U_{p^n}\equiv(\psi(p)p)^n\Phi\pmod{p^{n+2}}$ for all $n\ge1$, all $m$ $\iff$ $(S)$: $c'(pm)\equiv\psi(p)c'(m)\pmod{p^2}$ for all $m$ | **[proved]** §2.1 |
| V2 | $(S)$ at $p$ $\iff$ $p^2\mid\beta(n)$ for every $n$ with $p\mid n$, where $\beta:=c'\star(\mu\psi)$, i.e. $\beta(n)=\sum_{d\mid n}\mu(d)\psi(d)c'(n/d)$ | **[proved]** §2.2 |
| V3 | $(S)$ at **all** $p$ $\iff$ $\operatorname{rad}(n)^2\mid\beta(n)$ for all $n$ | **[proved]** §2.2 |
| V4 | **Master congruence** $n^2\mid\beta(n)$, equivalently $n^3\mid\sum_{d\mid n}\mu(d)\psi(d)\,d\,c(n/d)$. Implies V3, hence eq:magnetic, hence Cooper's free integration and the companion Lucas law | **[verified, $n\le1500$]** §2.3, §5 |
| V5 | $(S)$ itself, all three rows | **[verified, all $p\le199$, all $m\le1500/p$]** §5 |
| V6 | $(S)$ **fails** for Paşol–Zudilin's level-one magnetic forms $\Delta/E_4^2$, $E_4\Delta/E_6^2$, $E_6\Delta/E_4^3$ — so the congruence is *not* a consequence of magnetism, and no argument that only uses "meromorphic weight $4$, double pole at a CM point" can prove it | **[verified, $m\le200$]** §3 |
| V7 | Lagrange–Bürmann: $c(m)=[x^{m-1}]F(x)G(x)^m$, $G=x/q(x)$; equivalently $c'(m)=\operatorname{Res}_x\bigl(\eta\,q(x)^{-m}\bigr)$ with $\eta=\dfrac{l(x)\,dx}{x\sqrt{P(x)}\,F(x)}$, $P=1-2Bx+(B^2-4C)x^2$ | **[proved]**, **[verified, $m\le56$]** §4.1 |
| V8 | **The mechanism.** $(S)\bmod p$ $\iff$ $\mathcal C(\eta)\equiv\psi(p)\eta\pmod p$, $\mathcal C$ the Cartier operator; i.e. $a_{pj+p-1}\equiv\psi(p)a_j\pmod p$ for $\eta=\sum a_jx^j\,dx$ | **[proved]** (equivalence) + **[verified, all $p\le101$, sharp]** §4.2 |
| V9 | $a_{p-1}\equiv\psi(p)\pmod p$ — the $j=0$ slot; this *is* where the character $\psi$ comes from on the $x$-side | **[verified, $p\le199$]** §4.2 |
| V10 | The full Lucas congruence $a_{pj+r}\equiv a_ja_r\pmod p$ **fails**: only the $r=p-1$ slot holds. So V8 is a genuine Cartier-eigenvector statement, not a Dwork factorisation | **[verified, $p\le23$]** §4.3 |
| V11 | $\Phi\big|_4W_N=-\Phi$ for all three rows (from $u\|W_N=1/(Cu)$) | **[proved]** §6.1 |
| V12 | AL signs of $u$: $s_{10}$: $u\|W_2=u$, $u\|W_5=1/(Cu)$; $s_{18}$: $u\|W_2=1/(Cu)$, $u\|W_9=u$ | **[proved]** ($s_{10}$), **[num, 60 digits]** §6.1 |
| V13 | For $p\,\|\,N$: $\Phi\|U_p=\varepsilon_p\,p\,(\operatorname{Tr}^N_{N/p}\Phi-\Phi)$, so eq:magnetic at $p$ $\iff$ $\operatorname{Tr}^N_{N/p}\Phi\equiv(1+\varepsilon_p\psi(p))\Phi\pmod{p^2}$ | **[proved]** §6.2 |
| V14 | **THEOREM.** $\Phi_{s_7}\|U_7=7\,\Phi_{s_7}$ and $\Phi_{s_{10}}\|U_5=5\,\Phi_{s_{10}}$ *identically*. Hence eq:magnetic, $(S)$, V4 all hold **unconditionally** at $p=7$ for $s_7$ and at $p=5$ for $s_{10}$, for every $n$ and every $m$ | **[proved]** §6.3 |
| V15 | The exact cells are exactly the cells with $\varepsilon_p\psi(p)=-1$; equivalently $p\mid N$ and $p$ split in the CM order with class number one | **[proved]** ($\Leftarrow$), **[verified]** ($\Rightarrow$) §6.3–6.4 |
| V16 | $n\mid\beta(n)$, so $\Xi=D\log g$ with $g=\prod_n(1-q^n)^{-\beta(n)/n}\in\mathbf Z[[q]]$ ($\psi=\mathbf1$ rows), and $\Xi=\tfrac1{\sqrt{-3}}D\log G$, $G=\prod_n\bigl(\tfrac{1-\omega^2q^n}{1-\omega q^n}\bigr)^{\beta(n)/n}\in\mathbf Z[\omega][[q]]$ for $s_{18}$ | **[proved]** (the identity), **[verified, $n\le400$]** (integrality) §7 |
| V17 | $s_{18}$ at $p=3$ ($\psi(3)=0$) degenerates much further: $v_3(c(3^nm))\ge3n$, i.e. $v_3(c'(3^nm))\ge2n$, sharp | **[verified, $3^n m\le994$]** §6.5 |
| V18 | The exceptional cells $(s_7,p=2)$ and $(s_7,p=5)$ come from the accident $\gamma(p):=\beta(p)/p^2=-p$ (so $\beta(p)=-p^3$), not from a structural extra congruence | **[proved]** given the values, §6.6 |
| **GAP** | A proof of V8 (equivalently of $(S)\bmod p$) at a general prime, and its lift to $\bmod\,p^2$ | **open** §8 |

**Headline.** The task's target reduces exactly (V1–V3) to a divisibility of one integer
sequence, and the sharp form of that divisibility is $n^2\mid\beta(n)$ (V4) — one clean
statement that implies Cooper's free integration, the strong $p$-magnetic property and the
companion Lucas law simultaneously.  The mechanism is **not** modular-generic: it fails for
Paşol–Zudilin's level-one magnetic forms (V6).  It lives on the $x$-side, where it is a
Cartier-eigenvector congruence for one explicit differential $\eta$ built from the
Apéry-like series (V8), with the character $\psi$ appearing as $a_{p-1}\bmod p$ (V9).  Two
cells of the target are proved outright (V14), by an Atkin–Lehner trace plus rigidity.

---

## 1. Setting and conventions

Rows, from `lattice/cooper_sources/lib.gp` (re-copied into `lib.gp` here):

| row | $N$ | $B$ | $C$ | $u$ | $\psi$ |
|---|---|---|---|---|---|
| $s_7$ | $7$ | $13$ | $49$ | $\bigl(\eta(7\tau)/\eta(\tau)\bigr)^4$ | $\mathbf 1$ |
| $s_{10}$ | $10$ | $6$ | $25$ | $\bigl(\eta(5\tau)\eta(10\tau)/(\eta(\tau)\eta(2\tau))\bigr)^2$ | $\mathbf 1$ |
| $s_{18}$ | $18$ | $14$ | $1$ | $\bigl(\eta(2\tau)\eta(3\tau)^2\eta(18\tau)/(\eta(\tau)\eta(6\tau)^2\eta(9\tau))\bigr)^6$ | $\chi_{-3}$ |

$F=D\log u$, $g(u)=1+Bu+Cu^2$, $x=u/g(u)$, $\rho(u)=u\,dx/du=u(1-Cu^2)/g(u)^2$,
$$\Phi=F\,Dx=F^2\rho(u)=\sum_{m\ge1}c(m)q^m,\qquad
\Xi:=D^{-1}\Phi=\sum_{m\ge1}c'(m)q^m,\quad c'(m)=c(m)/m .$$
$A_n=[x^n]F$; $l(x)=\sum_{n\ge0}A_nx^{n+1}/(n+1)\in\mathbf Z[[x]]$ (Bogner), and
$\Xi=l(x(q))$.  $\psi$ is completely multiplicative; $\psi(p)=0$ only for $s_{18}$, $p=3$.
$U_p$ acts on $q$-expansions by $\sum a_mq^m\mapsto\sum a_{pm}q^m$.

---

## 2. The reduction: from $U_{p^n}$ to one divisibility

### 2.1 eq:magnetic $\iff$ $(S)$ **[proved]**

**Proposition 2.1.** Fix $p$.  The following are equivalent.

1. $\Phi|U_{p^n}\equiv(\psi(p)p)^n\Phi\pmod{p^{n+2}}$ for every $n\ge1$ (all coefficients);
2. $c(p^nm)\equiv\psi(p)^np^n c(m)\pmod{p^{n+2}}$ for every $n\ge1$ and every $m$ with $p\nmid m$;
3. $(S)$: $\;c'(pm)\equiv\psi(p)c'(m)\pmod{p^2}$ for every $m\ge1$;
4. $\Xi|U_{p^n}\equiv\psi(p)^n\,\Xi\pmod{p^2}$ for every $n\ge1$.

*Proof.*  (1)$\Leftrightarrow$(2) is the definition of $U_{p^n}$ on $q$-expansions together
with the observation that for $p\mid m$ the statement at $(n,m)$ follows from the statement
at $(n+v_p(m),m/p^{v_p(m)})$.  (2)$\Rightarrow$(3): write $m=p^vm_0$, $p\nmid m_0$.  Dividing
$c(p^jm_0)\equiv\psi(p)^jp^jc(m_0)\pmod{p^{j+2}}$ by $p^jm_0$ (a unit times $p^j$) gives
$c'(p^jm_0)\equiv\psi(p)^jc'(m_0)\pmod{p^2}$ for every $j\ge0$; subtracting the cases
$j=v+1$ and $j=v$ gives (3).  (3)$\Rightarrow$(4): $U_p$ preserves $p^2\mathbf Z_p[[q]]$, so
from $\Xi|U_p=\psi(p)\Xi+p^2R$, $R\in\mathbf Z_p[[q]]$, induction gives
$\Xi|U_{p^n}=\psi(p)^n\Xi+p^2R_n$.  (4)$\Rightarrow$(2): multiply the coefficient of $q^m$
by $p^nm$.  $\square$

So **the theorem to prove is $(S)$**, a $U_p$-eigenform-mod-$p^2$ statement for
$\Xi=D^{-1}\Phi$, which is not a modular form.

### 2.2 $(S)$ $\iff$ a divisibility of the Lambert coefficients **[proved]**

Since $\psi$ is completely multiplicative with $\psi(1)=1$, it is invertible for Dirichlet
convolution with inverse $\mu\psi$.  Define
$$\boxed{\;\beta(n):=\sum_{d\mid n}\mu(d)\,\psi(d)\,c'(n/d)\;}\qquad\text{so}\qquad
c'=\psi\star\beta,\quad \Xi=\sum_{e\ge1}\beta(e)\,\Lambda_\psi(q^e),\ \
\Lambda_\psi(t)=\sum_{d\ge1}\psi(d)t^d .$$
$\beta(n)\in\mathbf Z$ for all $n$ $\iff$ $c'(m)\in\mathbf Z$ for all $m$ $\iff$ magnetism
$m\mid c(m)$.

**Proposition 2.2 [proved].**  For every $n=p^vm$ with $p\nmid m$,
$$c'(pn)-\psi(p)c'(n)=\sum_{de=m}\psi(d)\,\beta\bigl(p^{v+1}e\bigr).$$
Consequently $(S)$ at $p$ holds $\iff$ $p^2\mid\beta(n)$ for every $n$ divisible by $p$;
and $(S)$ at every $p$ holds $\iff$ $\operatorname{rad}(n)^2\mid\beta(n)$ for every $n$.

*Proof.*  $c'(pn)=\sum_{de=pn}\psi(d)\beta(e)$ and
$\psi(p)c'(n)=\sum_{d'e=pn,\ p\mid d'}\psi(d')\beta(e)$ (substitute $d'=pd$), so the
difference is the sum over $de=pn$ with $p\nmid d$; those are exactly $e=p^{v+1}e'$,
$de'=m$.  For the second statement, fix $p$ and $v\ge0$: the family of congruences
$\sum_{de'=m}\psi(d)\beta(p^{v+1}e')\equiv0\pmod{p^2}$, $m$ running over integers prime to
$p$, is equivalent (invert the $\psi$-convolution over $\mathbf Z$, which is unitriangular)
to $\beta(p^{v+1}e')\equiv0\pmod{p^2}$ for all $e'$ prime to $p$.  Letting $v$ vary gives
exactly $p^2\mid\beta(n)$ for all $p\mid n$.  $\square$

This is the Paşol–Zudilin shape: if $\Phi=\Psi_D(f)$ is a Shimura–Borcherds lift of a
weight-$5/2$ form $f=\sum a(n)q^n$, then $c(n)=\sum_{d\mid n}\psi(d)d\,a(|D|n^2/d^2)$, so
$n\beta(n)=a(|D|n^2)$, and V3 reads $v_p\bigl(a(|D|n^2)\bigr)\ge v_p(n)+2$ for $p\mid n$.

### 2.3 The sharp form **[verified, $n\le1500$]**

Numerically the truth is stronger and cleaner (`02_gamma.log`, `05_verify.log`):
$$\boxed{\;n^2\mid\beta(n)\quad\text{for all }n\ \ \Bigl(\iff\ n^3\ \Bigm|\ \sum_{d\mid n}\mu(d)\psi(d)\,d\,c(n/d)\Bigr)\;}$$
verified for all $n\le1500$ in all three rows, and **sharp**: $\gamma(n):=\beta(n)/n^2$ is a
$p$-adic unit for some $n$ divisible by $p$, for every $p\le43$ and every row, except the
identically-zero cells $(s_7,7)$, $(s_{10},5)$ and the two cells $(s_7,p=2),(s_7,p=5)$ where
the minimum of $v_p(\gamma)$ is $1$.

$\gamma(n)$, $n=1,\dots,12$:

| row | $\gamma(n)$ |
|---|---|
| $s_7$ | $1,-2,1,2,-5,6,0,-18,39,-30,-57,234$ |
| $s_{10}$ | $1,-1,-1,2,0,-3,3,4,-11,0,27,-26$ |
| $s_{18}$ | $1,-1,2,-4,7,-14,31,-70,158,-367,877,-2120$ |

$\gamma$ is **not** multiplicative ($2286/2286$ coprime pairs fail for $s_{18}$).
$\gamma(p)\bmod p$ is the quantity $\lambda(p)$ of `cooper_sources/REPORT.md` §4.3.

**Hierarchy.**  $n^2\mid\beta(n)\ \Rightarrow\ \operatorname{rad}(n)^2\mid\beta(n)$ (=$(S)$,
V3) **and** $\Rightarrow\ n\mid\beta(n)$ (= the integral $q$-product, §7); the last two
statements are independent of each other.  So V4 is the single master statement.

---

## 3. The decisive negative test: it is not magnetism **[verified, $m\le200$]**

`03_pz.log`.  For Paşol–Zudilin's level-one magnetic weight-$4$ forms
$F_{4a}=\Delta/E_4^2$ (double pole at the disc $-3$ point $\rho$) and
$F_{4b}=E_4\Delta/E_6^2$ (disc $-4$), and for the doubly magnetic $F_6=E_6\Delta/E_4^3$:

* magnetism $m\mid c(m)$ holds (as it must);
* $\operatorname{rad}(n)^2\mid\beta(n)$ **fails already at $n=2$ or $n=3$**, for $\psi=\mathbf1$,
  $\chi_{-3}$ and $\chi_{-4}$ alike.  E.g. $F_{4a}$: $c'(1)=1$, $c'(2)=-252$, so
  $\beta(2)=-253$ and $4\nmid253$;
* the direct test gives $v_p\bigl(\gcd_m(c'(pm)-c'(m))\bigr)=0$ for $p=2,3,5,11$
  ($F_{4a}$) — the congruence is not even true modulo $p$.

Rescaling cannot repair this: $\beta$ is linear in $\Phi$, and one would need $p^2\mid\lambda$
for every $p$.

**Consequence.**  $(S)$ is *not* implied by "meromorphic of weight $4$ on $\Gamma_0(N)$ with
a double pole at a CM point and holomorphic at the cusps".  Any proof must use the extra
structure that Cooper's sources have and $\Delta/E_4^2$ lacks: that
$\Phi=F\,Dx$ with $F$ a weight-$2$ form which is *simultaneously* a power series in the
modular function $x$ with integral coefficients.  This kills route (D) of the task as a
stand-alone strategy and reorients everything to route (A).

---

## 4. The mechanism: a Cartier eigenvector on the $x$-line

### 4.1 The residue form of $c'$ **[proved; verified $m\le56$]**

Let $q(x)\in x+x^2\mathbf Z[[x]]$ be the compositional inverse of $x(q)$ and $G=x/q(x)$.

**Lemma 4.1 [proved].**  $c(m)=[x^{m-1}]\bigl(F(x)G(x)^m\bigr)$, i.e.
$c(m)=\operatorname{Res}_x\bigl(F(x)q(x)^{-m}dx\bigr)$.

*Proof.*  Lagrange inversion applied to $H=l$: $[q^m]H(x(q))=\tfrac1m[x^{m-1}]\bigl(H'(x)(x/q(x))^m\bigr)$,
and $H'=l'=F$, $[q^m]l(x(q))=c'(m)=c(m)/m$.  $\square$  (Verified for $m\le56$, `05_verify.log`.)

Since $Dx=x\sqrt{P(x)}\,F$ with $P(x)=1-2Bx+(B^2-4C)x^2$ (this is the identity
$\Phi=x\sqrt P\,F^2$ of `cooper_sources` §1, and $\sqrt P=(1-Bx)-2Cx\,u(x)\in\mathbf Z[[x]]$),
one has $dq/q=dx/(x\sqrt P\,F)$ and hence

$$\boxed{\;c'(m)=\operatorname{Res}_x\bigl(\eta\cdot q(x)^{-m}\bigr),\qquad
\eta:=\frac{l(x)}{x\sqrt{P(x)}\,F(x)}\,dx=\sum_{j\ge0}a_jx^j\,dx\in\mathbf Z[[x]]\,dx .}$$

First coefficients $a_0,\dots,a_9$ (`06_cartier.log`):

| row | $a_j$ |
|---|---|
| $s_7$ | $1, 11, 217, 4795, 112091, 2704835, 66613737, 1663943451, 41998320505, 1068520200601$ |
| $s_{10}$ | $1, 5, 70, 885, 12426, 177494, 2594740, 38428405, 575094490, 8673763130$ |
| $s_{18}$ | $1, 11, 138, 1829, 25054, 351318, 5016012, 72666075, 1065412062, 15778428914$ |

### 4.2 The Cartier congruence **[proved (equivalence) + verified, $p\le43$]**

For a $1$-form $\omega=\sum a_jt^j\,dt$ over $\mathbf F_p$ the Cartier operator is
$\mathcal C(\omega)=\sum_j a_{pj+p-1}^{1/p}t^j\,dt$; over $\mathbf F_p$ the $p$-th root is the
identity.  $\mathcal C$ is **intrinsic** (independent of the coordinate), and
$\operatorname{Res}_t\bigl(\omega\,Q(t^p)\bigr)=\operatorname{Res}_s\bigl(\mathcal C\omega\cdot Q(s)\bigr)$.

**Proposition 4.2 [proved].**  In the coordinate $q$ one has $\eta=\Xi\,dq/q$ and
$\mathcal C(\eta)=(\Xi|U_p)\,dq/q$.  Hence
$$(S)\bmod p\quad\Longleftrightarrow\quad \mathcal C(\eta)\equiv\psi(p)\,\eta\pmod p
\quad\Longleftrightarrow\quad a_{pj+p-1}\equiv\psi(p)\,a_j\pmod p\ \ (j\ge0).$$

**[verified]** `06_cartier.log`, `11_cartier_wide.log`: the congruence
$a_{pj+p-1}\equiv\psi(p)a_j\pmod p$ holds for every prime $p\le101$ and every row, over all
$j$ with $pj+p-1\le1200$ ($600$ tests at $p=2$ down to $11$ at $p=101$), and is **sharp** ($v_p$ of the gcd of the defects is exactly $1$,
except $2$ in the five cells $(s_7,3),(s_7,7),(s_{10},2),(s_{10},5),(s_{18},2)$).  The
control test with the opposite character gives $v_p=0$ throughout: the congruence *detects*
$\psi$.

**Corollary 4.3.**  Proving $\mathcal C(\eta)\equiv\psi(p)\eta\pmod p$ proves
$c'(pm)\equiv\psi(p)c'(m)\pmod p$ for all $m$ — the mod-$p$ half of the target.

**Where $\psi$ comes from.**  The $j=0$ slot is
$$\boxed{\;a_{p-1}\equiv\psi(p)\pmod p\;}$$
**[verified, every $p\le199$, all three rows]** (`07_lucas.log`, `11_cartier_wide.log`).  For $s_{18}$ this
reproduces $\chi_{-3}$ exactly, including $a_2\equiv0\pmod3$.  So on the $x$-side the
character is *not* imported from the modular side; it is the value at $x^{p-1}$ of one
explicit integral power series.

### 4.3 It is not a Dwork/Lucas factorisation **[verified]**

`07_lucas.log`.  The full Lucas congruence $a_{pj+r}\equiv a_ja_r\pmod p$ ($0\le r<p$)
**fails** massively: e.g. $s_7$, $p=11$: $566$ failures in $693$ tests; equivalently the
Dwork factorisation $\mathfrak A(x)\equiv\mathfrak A_{<p}(x)\mathfrak A(x^p)\pmod p$ fails
($\mathfrak A=\sum a_jx^j$).  Only the single slot $r=p-1$ survives, and it survives against
the *constant* $\psi(p)$ rather than against $a_{p-1}a_j$-type products.  Tested factor by
factor, none of $l(x)/x$, $1/\sqrt{P}$, $1/F$ satisfies the Cartier congruence on its own
(all three give $v_p=0$ at most primes); the property belongs to the product $\eta$.

So V8 is genuinely a statement that **$\eta$ spans a Cartier-stable line mod $p$ with
eigenvalue $\psi(p)$** — a unit-root/ordinary phenomenon, not a combinatorial Lucas law.

### 4.4 The exact ($\bmod p^2$) form: the Frobenius lift **[proved (equivalence); verified]**

Let $X_\sigma(x):=x\bigl(q(x)^p\bigr)$, so $q(X_\sigma)=q(x)^p$ exactly and
$X_\sigma\equiv x^p$ mod $p$.  **[verified, `09_frob.log`]** $X_\sigma\in\mathbf Z_p[[x]]$ for
$p\le7$ and all three rows (no coefficient has $p$ in its denominator, $\le120$ terms); e.g.
$s_7$, $p=3$: $X_\sigma=x^3+27x^4+639x^5+14778x^6+342324x^7+7994601x^8+\cdots$.  Let
$\sigma_*$ be the associated trace on $1$-forms,
$\operatorname{Res}_x(\omega\cdot Q\circ X_\sigma)=\operatorname{Res}_y(\sigma_*\omega\cdot Q)$;
then $\sigma_*\equiv\mathcal C\pmod p$ and, because
$c'(pm)=\operatorname{Res}_x\bigl(\eta\,q(x)^{-pm}\bigr)=\operatorname{Res}_x\bigl(\eta\,q(X_\sigma)^{-m}\bigr)$,

$$\boxed{\;(S)\ \Longleftrightarrow\ \sigma_*\eta\equiv\psi(p)\,\eta \pmod{p^2}. \;}$$

Equivalently, writing $q(x)^p=q(x^p)\bigl(1+p\,h(x)\bigr)$ with $h\in\mathbf Z_p[[x]]$ (this
is Frobenius on $\mathbf Z[[x]]$) and expanding $(1+ph)^{-m}$,
$$c'(pm)\equiv\operatorname{Res}_y\bigl(\mathcal C(\eta)q^{-m}\bigr)
-p\,m\operatorname{Res}_y\bigl(\mathcal C(\eta h)q^{-m}\bigr)\pmod{p^2},$$
so with $\mathcal C(\eta)=\psi(p)\eta+p\,\eta_1$ the missing second-order input is exactly
$$\eta_1\;\equiv\;D\bigl(\mathcal C(\eta h)\bigr)\pmod p$$
(the factor $m$ being $D=q\,d/dq$ under the residue pairing).  This is the precise shape of
the remaining gap; see §8.

**Remark (route B is not independent).**  Because $\Xi|U_{p^n}\equiv\psi(p)^n\Xi\pmod{p^2}$
(Prop. 2.1(4)), Hida's ordinary projector satisfies $e\Xi\equiv\Xi\pmod {p^2}$ whenever
$\operatorname{ord}\psi(p)\mid n!$; i.e. $\Xi$ is *already* its own ordinary projection mod
$p^2$ and computing $e\Xi$ returns no new information.  Also, the $p$-adic-modular-form
route can at best cover the primes for which the polar CM point is supersingular (inert or
ramified), whereas $(S)$ is verified at split primes too ($p=7,13,19,31,37,43$ for $s_7$).

---

## 5. The numerical map **[verified]**

`05_verify.log`, $M=1500$ coefficients, exact rational arithmetic.

* $m\mid c(m)$: holds, $m\le1500$, all rows.
* $n^2\mid\beta(n)$: holds, $n\le1500$, all rows.
* $(S)$: $e:=v_p\bigl(\gcd_m(c'(pm)-\psi(p)c'(m))\bigr)$ for every prime $p\le199$
  ($750$ tests at $p=2$ down to $7$ at $p=199$):

| row | $e$ |
|---|---|
| $s_7$ | $3$ at $p=2,5$; $\infty$ at $p=7$; $\;2$ at every other $p\le199$ |
| $s_{10}$ | $\infty$ at $p=5$; $\;2$ at every other $p\le199$ |
| $s_{18}$ | $2$ at every $p\le199$ |

* eq:magnetic $\Phi|U_{p^n}\equiv(\psi(p)p)^n\Phi\pmod{p^{e}}$, $p\le31$, $n\le4$: $e=n+2$
  exactly, except $e=n+3$ for $(s_7,p=2)$ and $(s_7,p=5)$, $e=\infty$ for $(s_7,7)$,
  $(s_{10},5)$, and $e=3n$ for $(s_{18},3)$.
* Refinements, both consequences of $n^2\mid\beta(n)$, both verified for $p\le19$,
  $n\le400/p$ (`10_exp.log`):
  $v_p\bigl(c'(pn)-\psi(p)c'(n)\bigr)\ \ge\ 2+2v_p(n)$ — no failures.
* $\beta$ is **not** multiplicative ($1238/1610$, $1072/1610$, $1610/1610$ coprime pairs fail).

Data files: `beta_s7.txt`, `beta_s10.txt`, `beta_s18.txt` ($\beta(n)$, $n\le600$),
`gamma_s*.txt` ($\gamma(n)=\beta(n)/n^2$, $n\le800$).

---

## 6. What is proved: the Atkin–Lehner cells

### 6.1 The Fricke sign **[proved]**

**Lemma 6.1 [proved].**  $u|W_N=1/(Cu)$ implies
$$F\big|_2W_N=-F,\qquad \rho(u)\big|W_N=-\rho(u),\qquad \Phi\big|_4W_N=-\Phi .$$

*Proof.*  For weight $0$, $D(h\circ\gamma)=(Dh)|_2\gamma$ with
$h|_k\gamma=\det(\gamma)^{k/2}(c\tau+d)^{-k}h(\gamma\tau)$, so
$F|_2W_N=D\log(u|W_N)=D\log\bigl(1/(Cu)\bigr)=-F$.  Next, $g(1/(Cu))=g(u)/(Cu^2)$ and
$1-C\cdot(1/(Cu))^2=-(1-Cu^2)/(Cu^2)$, whence
$$\rho\Bigl(\frac1{Cu}\Bigr)=\frac{1}{Cu}\cdot\frac{-(1-Cu^2)}{Cu^2}\cdot\frac{C^2u^4}{g(u)^2}
=-\frac{u(1-Cu^2)}{g(u)^2}=-\rho(u).$$
Finally $\Phi=F^2\rho(u)$ and $(F^2)|_4W_N=(F|_2W_N)^2=F^2$.  $\square$

(In passing this re-proves $x|W_N=x$, since $x=u/g(u)$ is fixed by $u\mapsto1/(Cu)$.)

**Lemma 6.2 [proved for $s_{10}$; num, 60 digits in general].**  For $s_{10}$,
$u|W_2=u$ and hence $u|W_5=1/(Cu)$, $\Phi|_4W_5=-\Phi$, $\Phi|_4W_2=+\Phi$.

*Proof.*  $W_2$ normalises $\Gamma_0(10)$ and permutes its cusps by $1\leftrightarrow2$,
$5\leftrightarrow10$; by the Ligozat divisor of $u$ (`cooper_sources` §2.2,
$[1,1,-1],[2,1,-1],[5,1,1],[10,1,1]$) this leaves $\operatorname{div}(u)$ invariant, so
$u|W_2=\epsilon u$ with $\epsilon^2=1$.  The fixed-point count for Atkin–Lehner involutions
gives $\nu(W_2)=h(-4)\prod_{p\mid5}\bigl(1+(\tfrac{-4}{p})\bigr)+h(-8)\prod_{p\mid5}\bigl(1+(\tfrac{-8}{p})\bigr)=1\cdot2+1\cdot0=2$:
$W_2$ has exactly two fixed points on $X_0(10)$, of discriminant $-4$ — i.e. precisely the two
poles of $\Phi$, where $u=u_\pm=(-3\pm4i)/25\ne0$.  If $\epsilon=-1$ then $u$ would vanish
there.  Hence $\epsilon=+1$, and $u|W_5=u|W_2|W_{10}=(1/(Cu))|W_2=1/(Cu)$.  $\square$

**[num, 60 digits]** `08_al.log` confirms all signs by direct evaluation of the eta quotients
at $\tau=0.13+0.37i$:

| row | $u|W_2$ | $u|W_5$ | $u|W_7$ | $u|W_9$ | $u|W_{18}$ | $\varepsilon_p$ on $\Phi$ |
|---|---|---|---|---|---|---|
| $s_7$ | – | – | $1/(Cu)$ | – | – | $\varepsilon_7=-1$ |
| $s_{10}$ | $u$ | $1/(Cu)$ | – | – | – | $\varepsilon_2=+1,\ \varepsilon_5=-1$ |
| $s_{18}$ | $1/(Cu)$ | – | – | $u$ | $1/(Cu)$ | $\varepsilon_2=-1,\ \varepsilon_9=+1$ |

(This confirms independently the correction of `cooper_sources` §6.2: $F|_2W_9=+F$ for $s_{18}$.)

### 6.2 The trace criterion **[proved]**

For $p\,\|\,N$ and $f$ of weight $k$ on $\Gamma_0(N)$, with
$(f|_kW_p)(\tau)$ the Atkin–Lehner involution at $p$,
$$\operatorname{Tr}^{N}_{N/p}(f)\;=\;f+p^{1-k/2}\,f\big|W_p\big|U_p . \tag{6.1}$$
*(Standard; verified here on $f=E_4$: at $(N,p)=(7,7)$ it gives $E_4+7^{-1}\cdot49E_4=8E_4$,
the index; at $(N,p)=(10,5)$ it gives $E_4+5^{-1}\cdot25E_4=6E_4$, again the index.)*

With $\Phi|W_p=\varepsilon_p\Phi$, $k=4$:
$$\boxed{\;\Phi\big|U_p=\varepsilon_p\,p\,\bigl(\operatorname{Tr}^N_{N/p}\Phi-\Phi\bigr),\qquad
\text{so}\qquad \Phi|U_p\equiv\psi(p)p\,\Phi\ (p^3)\iff
\operatorname{Tr}^N_{N/p}\Phi\equiv(1+\varepsilon_p\psi(p))\Phi\ (p^2).}$$

The right-hand side has **new content only when $1+\varepsilon_p\psi(p)=0$**, i.e.
$\varepsilon_p=-\psi(p)$: then it says the trace *vanishes mod $p^2$*, and rigidity upgrades
this to vanishing outright.  When $\varepsilon_p\psi(p)=+1$ the criterion is a tautological
rewriting of the congruence (verified as such: `04_fricke.log` shows
$\operatorname{Tr}^{10}_5\Phi\equiv2\Phi\pmod4$ for $s_{10}$ at $p=2$ and
$\operatorname{Tr}^{18}_9\Phi\equiv2\Phi\pmod4$ for $s_{18}$ at $p=2$, all $m\le150$;
in both cells $1+\varepsilon_2\psi(2)=2$).  Note the check also fixes the signs: with the
opposite sign of $\varepsilon_2$ the congruence fails at $m=1,5,9,\dots$

The cells with $\varepsilon_p\psi(p)=-1$ are exactly $(s_7,p=7)$ and $(s_{10},p=5)$.

### 6.3 THEOREM (the two exact cells) **[proved]**

**Theorem 6.3.**  $\Phi_{s_7}\big|U_7=7\,\Phi_{s_7}$ and $\Phi_{s_{10}}\big|U_5=5\,\Phi_{s_{10}}$,
identically.  Consequently $c(pm)=p\,c(m)$, $c'(pm)=c'(m)$, $\beta(pn)=0$ for all $n$, and
eq:magnetic holds with equality at $p=7$ for $s_7$ and at $p=5$ for $s_{10}$, for every
$n\ge1$ and every $m$.

*Proof.*  Both $\Phi$'s vanish at every cusp: $u$ has divisor supported on cusps with values
$0$ or $\infty$ there, hence $x=u/g(u)\to0$ at every cusp, hence $Dx$ vanishes at every cusp,
and $F$ is holomorphic there.  Both have, on $\mathbf H$, only double poles (in $\tau$) at
the CM points where $g(u)=0$.

*Case $s_7$, $p=N=7$.*  By Lemma 6.1, $\varepsilon_7=-1$, so (6.1) gives
$T:=\operatorname{Tr}^7_1\Phi=\Phi-\tfrac17\Phi|U_7$.  $T$ is a meromorphic weight-$4$ form on
$\mathrm{SL}_2(\mathbf Z)$: it is a finite sum of $\Phi|\gamma$, each of which is holomorphic
on $\mathbf H$ away from the $\mathrm{SL}_2(\mathbf Z)$-orbit of the two elliptic points of
order $3$ of $X_0(7)$ — all $\mathrm{SL}_2(\mathbf Z)$-equivalent to $\rho$ — where the pole
order in $\tau$ is at most $2$; and $T$ vanishes at $\infty$.  The valence formula for weight
$4$ on $\mathrm{SL}_2(\mathbf Z)$,
$$\operatorname{ord}_\infty T+\tfrac12\operatorname{ord}_iT+\tfrac13\operatorname{ord}_\rho T
+\sum{}'\operatorname{ord}_PT=\tfrac13,$$
with $\operatorname{ord}_\rho T\ge-2$, $\operatorname{ord}_iT\ge0$ and $\sum'\ge0$, forces
$\operatorname{ord}_\infty T\le\tfrac13+\tfrac23=1$.  But
$[q^1]T=c(1)-c(7)/7=1-7/7=0$ and $[q^2]T=c(2)-c(14)/7=-14-(-98)/7=0$ **[exact]**, so
$\operatorname{ord}_\infty T\ge3$ and therefore $T=0$, i.e. $\Phi|U_7=7\Phi$.

*Case $s_{10}$, $p=5$.*  By Lemma 6.2, $\varepsilon_5=-1$, so
$T:=\operatorname{Tr}^{10}_2\Phi=\Phi-\tfrac15\Phi|U_5$ is a meromorphic weight-$4$ form on
$\Gamma_0(2)$, holomorphic on $\mathbf H$ except for a pole of order $\le2$ in $\tau$ at the
image of the two disc-$(-4)$ points of $X_0(10)$.  On $X_0(2)$ there is exactly one Heegner
point of discriminant $-4$ (the number is $h(-4)\cdot\#\{\beta\bmod4:\beta^2\equiv-4\ (8)\}=1\cdot1$),
namely the elliptic point of order $2$; and $T$ vanishes at both cusps.  The valence formula
for $\Gamma_0(2)$ ($\mu=3$, $\nu_2=1$, $\nu_3=0$, two cusps, $k\mu/12=1$) gives
$$\operatorname{ord}_\infty T+\operatorname{ord}_0T+\tfrac12\operatorname{ord}_{e}T+\sum{}'=1,$$
so $\operatorname{ord}_\infty T\le1+1=2$.  **[exact]** $[q^m]T=c(m)-c(5m)/5=0$ for all
$m\le60$ (`04_fricke.log`), in particular for $m\le3$; hence $T=0$ and $\Phi|U_5=5\Phi$.
$\square$

*(The two finite verifications are exact rational computations; the proof needs only
$m\le2$, resp. $m\le3$.)*

### 6.4 Why exactly these cells **[proved $\Leftarrow$; verified $\Rightarrow$]**

The criterion of §6.2 is exact iff $\varepsilon_p=-\psi(p)$, which requires $p\,\|\,N$.  The
conceptual reason $U_p$ can preserve the polar divisor is that the pole is a Heegner point
of discriminant $D$ on $X_0(N)$ and $p\mid N$ splits in the order of discriminant $D$ with
class number one, so the Atkin–Lehner/$U_p$ correspondence permutes the disc-$D$ Heegner
points:

| row | $D$ | $p\mid N$ | $p$ in $\mathcal O_D$ | $h(D)$ | exact? |
|---|---|---|---|---|---|
| $s_7$ | $-3$ | $7$ | split ($-3\equiv2^2$) | $1$ | **yes** |
| $s_{10}$ | $-4$ | $2$ | ramified | $1$ | no ($\varepsilon_2=+1$) |
| $s_{10}$ | $-4$ | $5$ | split ($5=(2+i)(2-i)$) | $1$ | **yes** |
| $s_{18}$ | $-36$ | $2$ | ramified | $1$ | no ($\varepsilon_2=-1$, $\psi(2)=-1$) |
| $s_{18}$ | $-36$ | $3$ | divides the conductor | – | no (see §6.5) |

### 6.5 $s_{18}$ at $p=3$: $\psi(3)=0$ and a deeper degeneration **[verified]**

$\chi_{-3}(3)=0$, so $(S)$ reads $9\mid c'(3m)$.  In fact **[verified, $3^nm\le994$,
`09_frob.log`]**
$$v_3\bigl(c(3^nm)\bigr)\ \ge\ 3n,\qquad\text{i.e.}\qquad v_3\bigl(c'(3^nm)\bigr)\ \ge\ 2n,$$
sharp for $n\le5$ (minima exactly $3n$ over $331,110,36,12,4$ tests).  So the "$\psi(p)=0$"
cell is not merely $\Phi|U_3\equiv0\pmod{27}$ but $\Phi|U_{3^n}\equiv0\pmod{3^{3n}}$, i.e.
$\beta(3^je)\equiv0\pmod{3^{2j}}$ — which is precisely what $n^2\mid\beta(n)$ predicts.
$\Phi|U_3/27$ is not proportional to $\Phi$ ($2,-116,4266,\dots$ against $1,-10,54,\dots$).
$3^2\mid18$, so the trace identity (6.1) does not apply and this cell is not covered by §6.3.

### 6.6 The two anomalous cells $(s_7,p=2)$ and $(s_7,p=5)$ **[proved, given the values]**

`cooper_sources` §4.1 records $e=3$ (not $2$) there.  The reason is visible in $\gamma$:
$$\gamma_{s_7}(2)=-2,\qquad \gamma_{s_7}(5)=-5,\qquad\text{i.e.}\qquad
\beta_{s_7}(2)=-2^3,\quad \beta_{s_7}(5)=-5^3 .$$
Since $\gamma(p)\equiv0\pmod p$ in exactly these two cells (and in no other cell with
$p\le43$), $\beta(p)$ acquires one extra power of $p$ and the minimum of $v_p(\gamma(pe))$
over $e$ is $1$ rather than $0$ (`02_gamma.log`).  There is no additional structure: $\gamma$
is not a character and $\gamma(p)\ne-p$ at $p=3,11,13,\dots$ ($\gamma_{s_7}(3)=1$,
$\gamma_{s_7}(11)=-57$).  So this is an accident of small numbers on top of $n^2\mid\beta(n)$,
not an extra theorem.

---

## 7. The $q$-product, and Cartier's theorem as a sanity check

**Lemma 7.1 [proved].**  With $e_n:=\beta(n)/n$,
$$\sum_{m\ge1}\frac{c'(m)}{m}q^m=\sum_{n\ge1}e_n\sum_{k\ge1}\frac{\psi(k)}{k}q^{nk}.$$
For $\psi=\mathbf1$ this says $\Xi=D\log g$ with $g=\prod_{n\ge1}(1-q^n)^{-e_n}$; for
$\psi=\chi_{-3}$, $\Xi=\tfrac1{\sqrt{-3}}D\log G$ with
$G=\prod_{n\ge1}\bigl((1-\omega^2q^n)/(1-\omega q^n)\bigr)^{e_n}$, $\omega=e^{2\pi i/3}$.

**[verified, $n\le400$]** $n\mid\beta(n)$, so all $e_n\in\mathbf Z$, in all three rows; and
directly, $g=\exp\bigl(\sum c'(m)q^m/m\bigr)\in\mathbf Z[[q]]$ for $s_7$ and $s_{10}$
($n\le80$), while for $s_{18}$ it is *not* integral (first failure $n=3$) — exactly as the
character predicts.  Examples:
$$e^{s_7}=1,-4,3,8,-25,36,0,-144,351,-300,\dots;\quad
e^{s_{10}}=1,-2,-3,8,0,-18,21,32,-99,0,\dots;\quad
e^{s_{18}}=1,-2,6,-16,35,-84,217,-560,\dots$$
$$g_{s_7}=1+q-3q^2+0q^3+14q^4-23q^5-17q^6+125q^7-\cdots,\qquad
g_{s_{10}}=1+q-q^2-4q^3+5q^4+11q^5-\cdots$$

**Why this matters.**  Over $\mathbf F_p$, Cartier's theorem says $\mathcal C\omega=\omega$
iff $\omega$ is logarithmic, $\omega=dg/g$.  So for the $\psi=\mathbf1$ rows the mod-$p$
mechanism V8 is *equivalent* to the assertion that $\Xi\,dq/q$ is a logarithmic differential
mod $p$, and by Dieudonné–Dwork this is equivalent to the $p$-integrality of the $q$-product
$g$.  This is one more independent avatar of the same congruence, and it identifies the
statement to prove as an **integrality of a $q$-product** — the same flavour as the
integrality theorems for mirror maps.  ($n\mid\beta(n)$ and $\operatorname{rad}(n)^2\mid\beta(n)$
are independent; both follow from $n^2\mid\beta(n)$.)

---

## 8. The gap, stated precisely

**Proved unconditionally:** eq:magnetic at $(s_7,p=7)$ and $(s_{10},p=5)$ (Theorem 6.3);
the reduction of the general case to a single divisibility (Props. 2.1, 2.2); the exact
identification of the mechanism as a Cartier/Frobenius eigenvector statement (§4).

**Not proved, in decreasing order of strength.**

1. **(Master)** $n^2\mid\beta(n)$, i.e. $n^3\mid\sum_{d\mid n}\mu(d)\psi(d)\,d\,c(n/d)$.
   [verified $n\le1500$]  In Shimura–Borcherds language: $n^3\mid a(|D|n^2)$, three powers
   better than Paşol–Zudilin's $n\mid a(n^2)$, and *false* for their level-one examples (§3).
2. **(Target)** $(S)$, equivalently $\operatorname{rad}(n)^2\mid\beta(n)$, equivalently
   $$\sigma_*\eta\equiv\psi(p)\,\eta\pmod{p^2}$$
   for the Frobenius lift $\sigma:x\mapsto x(q(x)^p)$ of §4.4.  [verified $p\le199$]
3. **(Mod $p$)** $\mathcal C(\eta)\equiv\psi(p)\eta\pmod p$, i.e.
   $a_{pj+p-1}\equiv\psi(p)a_j\pmod p$ for $\eta=l(x)\,dx/(x\sqrt{P}F)$.  [verified $p\le101$]
   *This is the smallest missing brick*: it is a congruence for one explicit integral power
   series built from the Apéry-like recurrence, it implies $(S)\bmod p$, and $(2)\Rightarrow(3)$
   is a one-line expansion.  Its $j=0$ case, $a_{p-1}\equiv\psi(p)\pmod p$, is already the
   statement that determines the character.
4. Given (3), the lift to $\bmod\ p^2$ needs exactly one further identity mod $p$:
   $\eta_1\equiv D\bigl(\mathcal C(\eta h)\bigr)$, where
   $\mathcal C(\eta)=\psi(p)\eta+p\eta_1$ and $q(x)^p=q(x^p)(1+ph)$ (§4.4).

**Routes that are now excluded or downgraded.**

* Route (D) (Shimura–Borcherds / traces of singular moduli) cannot work on the polar data
  alone: §3 shows the congruence is invisible to "weight $4$, double pole at a CM point".
  It could still work if one *identifies the specific* weight-$5/2$ input $f$ at level $4N$
  and proves $n^3\mid a(|D|n^2)$ — but that divisibility is exactly the master congruence and
  is stronger than the Hecke-operator argument of Paşol–Zudilin delivers.
* Route (B) ($p$-adic modular forms, Hida): §4.4 Remark — the ordinary projection of $\Xi$ is
  $\Xi$ itself mod $p^2$, so it carries no new information; and the method only reaches
  inert/ramified $p$ whereas the congruence holds at split $p$ too.
* Route (C) (level-$Np$ identity, modular equation): usable only at $p\mid N$, where it gives
  Theorem 6.3; for $p\nmid N$ the trace $\operatorname{Tr}^{Np}_N\Phi=(p+1)\Phi$ is an
  identity with no content, and the mod-$p^3$ structure of the modular polynomial is exactly
  as hard as the target.
* Route (A) is the live one, in the sharp form (3) above.

---

## 9. Consequences for the repository

1. `paper/companions/main.tex`, Theorem `cooper`: hypothesis eq:magnetic can be replaced by
   the *single* congruence $(S)$ (Prop. 2.1) — the "one congruence per prime power" caveat in
   the current proof of (i) is unnecessary once one works with $\Xi$ rather than $\Phi$.  The
   verification block can be strengthened to $p\le199$, $m\le1500$, and can quote the sharp
   form $n^2\mid\beta(n)$.
2. The two exact cells are now **theorems**, not verifications: eq:magnetic holds
   unconditionally at $p=7$ for $s_7$ and $p=5$ for $s_{10}$ (Theorem 6.3), with the
   structural reason $\varepsilon_p=-\psi(p)$.  `cooper_sources/REPORT.md` §1 and §7's third
   bullet should be updated accordingly.
3. `cooper_sources/REPORT.md` §7 asks "why $p^3$ and not $p^4$, and why $p^4$ exactly at
   $(s_7,2),(s_7,5)$": answered in §6.6 — the sharp statement is $n^2\mid\beta(n)$ and the two
   cells are the accident $\beta_{s_7}(p)=-p^3$ at $p=2,5$.
4. `cooper_sources/REPORT.md` §7's first bullet proposes the Paşol–Zudilin route as "the
   natural attack".  §3 shows this cannot succeed at the level of the polar data alone.
5. The $W_9$-sign correction of `cooper_sources` §6.2 is confirmed independently (§6.1 table).

---

## 10. Files

| file | contents |
|---|---|
| `lib.gp` | rows, eta quotients, $F$, $x$, $\Phi$, $c$, $c'$, $\beta$, the recurrence |
| `01_beta.gp/.log` | $\beta=c'\star(\mu\psi)$; $\operatorname{rad}(n)^2\mid\beta(n)$; multiplicativity of $\beta$ |
| `02_gamma.gp/.log` | $n^2\mid\beta(n)$, sharpness, $\gamma=\beta/n^2$ |
| `03_pz.gp/.log` | the negative test on $\Delta/E_4^2$, $E_4\Delta/E_6^2$, $E_6\Delta/E_4^3$ |
| `04_fricke.gp/.log` | traces $\operatorname{Tr}^N_{N/p}\Phi$; the two exact cells; the mod-$4$ cells |
| `05_verify.gp/.log` | $M=1500$: magnetism, $n^2\mid\beta$, $(S)$ for $p\le199$, eq:magnetic $n\le4$, Lagrange–Bürmann |
| `06_cartier.gp/.log` | $\eta$, and $\mathcal C(\eta)\equiv\psi(p)\eta\ (p)$ with character controls |
| `07_lucas.gp/.log` | $a_{p-1}\equiv\psi(p)$; failure of the full Lucas/Dwork factorisation |
| `08_al.gp/.log` | Atkin–Lehner action on $u$, 60 digits |
| `09_frob.gp/.log` | $p$-integrality of the Frobenius lift $X_\sigma$; $s_{18}$ at $p=3$ |
| `10_exp.gp/.log` | $n\mid\beta(n)$, the $q$-product $g$, the refined $(S{+}{+})$ |
| `11_cartier_wide.gp/.log` | the Cartier congruence for all $p\le101$ ($\eta$ to $x^{1200}$); $a_{p-1}\equiv\psi(p)$ for $p\le199$ |
| `beta_s*.txt`, `gamma_s*.txt` | $\beta(n)$ ($n\le600$), $\gamma(n)=\beta(n)/n^2$ ($n\le800$) |

Run with `gp -q <file>`; total runtime under one minute.
