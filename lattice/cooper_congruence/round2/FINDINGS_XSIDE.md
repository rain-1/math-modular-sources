# The $x$-side attack on Cooper's congruence: routes (A)–(D)

*Round 2, 2026-09-02.  All computations in PARI/GP 2.15, exact over $\mathbf Q$ or over
$\mathbf F_p$/$\mathbf Z/p^k$; scripts `40_`–`47_` with matching `.log` files in this
directory.  Nothing outside `lattice/cooper_congruence/round2/` was touched.  Background:
`lattice/cooper_congruence/REPORT.md` (round 1) §§4, 7, 8.*

Notation as in round 1: rows $s_7,s_{10},s_{18}$ with $(B,C)=(13,49),(6,25),(14,1)$ and
$\psi=\mathbf1,\mathbf1,\chi_{-3}$; $A_n=[x^n]F$; $L_n=A_{n-1}/n\in\mathbf Z$;
$l(x)=\sum_{n\ge1}L_nx^n$; $P=1-2Bx+(B^2-4C)x^2=(1-\lambda_1x)(1-\lambda_2x)$;
$y:=\sqrt P\in\mathbf Z[[x]]$; $\Xi=l(x(q))$; $\eta=\Xi\,\frac{dq}q=a(x)\,dx$ with
$a=l/(xyF)$; $\mathcal C$ the Cartier operator.  **Target (3):**
$\mathcal C(\eta)\equiv\psi(p)\eta \pmod p$, i.e. $a_{pj+p-1}\equiv\psi(p)a_j$.
Throughout $p\nmid N$ unless said otherwise.  Write
$$\operatorname{comp}_0(S):=\sum_{n\ge0}[x^{pn}]S\cdot x^n .$$

---

## 0. Verdict table

| # | statement | status |
|---|---|---|
| **W1** | Route (A) **as posed in the task is empty**: for weight $0$ the geometric Hecke correspondence is $T_p=p\,U_p+V_p$ (not $U_p+V_p$), so the Kronecker congruence $\Phi_p(X,Y)\equiv(X^p-Y)(X-Y^p)$ says exactly $p\,U_pf\equiv0$ — a tautology.  The modular polynomial **mod $p$** carries no information about $U_p(x^b)\bmod p$; one needs it mod $p^2$ | **[proved]** §1 |
| **W2** | The consequence I first drew from the wrong normalisation, "$U_p(x^b)\equiv0$", is **false**, and so is the reduction it would give ($L_{pa}\equiv\psi(p)L_a$) | **[refuted, $p\le97$]** §1 |
| **W3** | **The correct $x$-line formula.** With $H:=(yF)^{p-1}\in\mathbf Z[[x]]$ and $f=\sum f_nx^n\in\mathbf Z_p[[x]]$, $$f\big|U_p\ \equiv\ f_0+\sum_{n\ge1}\bigl[x^{pn}\bigr]\bigl(f\,H\bigr)\,x^n \pmod p .$$ In particular $U_p(x^b)=\sum_{n\ge1}[x^{pn-b}]H\,x^n$ for $b\ge1$, and $U_p(1)=1$ $\iff$ $[x^{pn}]H=0$ $(n\ge1)$ | **[proved]** + **[verified]** §2 |
| **W4** | Hence **the target is** $\ \operatorname{comp}_0(l\,H)=\psi(p)\,l$, i.e. $[x^{pn}](l\,H)\equiv\psi(p)L_n$ $(n\ge1)$: $U_p$ in the $x$-coordinate is $f\mapsto f_0+\operatorname{comp}_0(fH)-\,(fH)_0$, it fixes $1$, and the claim is that $l$ is an eigenvector with eigenvalue $\psi(p)$ | **[proved]** (equivalence) + **[verified, $p\le29$]** §2 |
| **W5** | **$H$ is a rational function of $x$ mod $p$**, of degree $<p$ in numerator and denominator: the Lucas property $A_{pj+r}\equiv A_jA_r$ gives $F\equiv F_{<p}(x)F(x^p)$, hence $F^{p-1}\equiv1/F_{<p}$ and $$\boxed{H\ \equiv\ \frac{P(x)^{(p-1)/2}}{F_{<p}(x)}\pmod p},\qquad F_{<p}:=\sum_{n<p}A_nx^n .$$ | **[proved from Lucas]** + **[verified, $p\le23$]** §3 |
| **W6** | **The smallest missing brick, in fully explicit finite form:** $$\boxed{\ \bigl[x^{pn}\bigr]\Bigl(l(x)\,\frac{P(x)^{(p-1)/2}}{F_{<p}(x)}\Bigr)\ \equiv\ \psi(p)\,L_n\pmod p,\qquad n\ge1.\ }$$ Only $l$ is transcendental; the kernel is a ratio of two polynomials of degree $<p$ | **[proved equivalent to (3)]** + **[verified, all $p\le53$, $x^{1200}$, all rows]** §3 |
| **W7** | **The denominator is the supersingular polynomial.**  $\deg_xF_{<p}\bmod p=\bigl\lfloor \mu(N)(p-1)/(6\deg x)\bigr\rfloor$ exactly ($\mu=8,18,36$; $\deg x=2,4,8$) — i.e. the zeros of $F_{<p}$ are the supersingular points of $X_0(N)/\mathbf F_p$ seen on the $x$-line, doubled.  $F_{<p}$ is the Hasse invariant of the family (Igusa).  So (3) is an **ordinarity/unit-root statement at the supersingular divisor** | **[verified, $20$ primes $\times\,3$ rows, exact]** §3 |
| **W8** | **Intermediate ($\sqrt P$-free) form, not using Lucas.**  With $b:=l/(xF)\in\mathbf Z[[x]]$, $b_0=1$: $$[x^{pj+p-1}]\bigl(b\cdot P^{(p-1)/2}\bigr)\equiv\psi(p)\,b_j\pmod p .$$ The transcendental $1/\sqrt P$ is replaced by the **polynomial** $P^{(p-1)/2}$ of degree $p-1$ | **[proved equivalent]** + **[verified, all odd $p\le89$, $x^{1200}$]** §2 |
| **W9** | The same congruence **mod $p^2$ is false** | **[refuted, $p\le31$]** §2 |
| **W10** | The "conic" part alone is not the character: $[x^{p-1}]P^{(p-1)/2}=\bigl(\tfrac{B^2-4C}{p}\bigr)$, which is $\chi_{-3},\chi_{-4},\bigl(\tfrac3\cdot\bigr)$ for the three rows — **not** $\psi$.  The character comes from $b$ | **[proved]** + **[verified]** §2 |
| **W11** | **New supercongruence, and the source of $\psi$ on the combinatorial side:** $$\boxed{\ A_{p-1}\ \equiv\ \kappa\,\psi(p)\,p \pmod{p^2},\qquad \kappa=-2,\,-3,\,-3\ \ (s_7,s_{10},s_{18}),\ }$$ equivalently $\psi(p)\equiv\kappa^{-1}A_{p-1}/p$.  For $s_{18}$ this **manufactures $\chi_{-3}$ out of the Apéry-like recurrence alone** | **[verified, all $5\le p\le300$]**; **[proved for $s_{10}$, $p\ge5$]** §4 |
| **W12** | Refinement: $A_{p-1}\equiv-2p-6p^2$ $(s_7)$, $\equiv-3p-12p^2$ $(s_{10})$ mod $p^3$ | **[verified, $p\le59$]** §4 |
| **W13** | The bridge to the $j=0$ slot: $\kappa\,a_{p-1}\equiv A_{p-1}/p\equiv\kappa\psi(p)\pmod p$ | **[verified, $5\le p\le300$]** §4 |
| **W14** | No Lucas lift of W11: $A_{pa-1}\equiv\psi(p)pA_{a-1}\ (p^2)$ **fails** for $a\ge2$; equivalently $L_{pa}\not\equiv\psi(p)L_a$.  The congruence (3) is **not** reducible to a statement about $A$ alone | **[refuted, $p\le97$]** §1, §4 |
| **W15** | Lucas census: $A_{pj+r}\equiv A_jA_r$ **holds** ($0$ failures, $n\le600$, $p\le23$, all rows); it **fails** for $L_n$, for $a_j$ and for $b_j$ | **[verified]**; Lucas for $A$ also **[proved, Malik–Straub]** §4 |
| **W16** | Beukers/Dwork $A_{mp^s}\equiv A_{mp^{s-1}}\ (p^{3s})$ holds for $s_7,s_{10}$ but **fails for $s_{18}$**, which only reaches $p^{2s}$ | **[verified]** §4 |
| **W17** | No binomial-sum representation for $s_{18}$ was found (exhaustive search over products of $8$ standard blocks, exponents $\le5$, $20$ prefactor variants); $s_7,s_{10}$ formulas confirmed for $n\le40$ | **[verified / search failed]** §4 |
| **W18** | **(D)** The round-1 identity $\eta_1\equiv D\bigl(\mathcal C(\eta h)\bigr)\pmod p$ holds | **[verified, $p\le7$, all rows]** §5 |
| **W19** | **$\eta_1$ is exact mod $p$**: $\mathcal C(\eta_1)=0$.  Equivalently $a_{p^2n-1}\equiv\psi(p)a_{pn-1}\pmod{p^2}$.  This is forced by W18 (a $D$-image is exact), and is the $v_p=1$ case of W20 | **[proved from W18]** + **[verified, $p\le31$]** §5 |
| **W20** | The **$\psi$-twisted Dieudonné–Dwork tower** $$a_{pn-1}\equiv\psi(p)\,a_{n-1}\pmod{p^{\,1+v_p(n)}}\qquad(n\ge1)$$ holds for **all three rows**, extending round 1's V19 (which was $\psi=\mathbf1$ only) to the character row $s_{18}$ | **[verified, $p\le31$, $n\le1500/p$]** §5 |
| **W21** | $\eta_1$ has **no** independent description in the obvious sense: $\eta_1/\eta$ is not a rational function of $x$ mod $p$ (degree $\le12$ excluded, $p\le19$), and $\eta_1$ is not in the $\mathbf F_p$-span of $\{\eta,\ x\eta,\ F'/F,\ P'/P,\ xF/l,\ 1\}$ | **[refuted, stated ranges]** §5 |

**Headline.**  Route (A) as literally posed is empty (W1) — but the *right* version of it
is the live route, and it converts the smallest missing brick into a statement about one
**explicit rational function of degree $<p$**, namely $H=P^{(p-1)/2}/F_{<p}$, whose
denominator is the **supersingular polynomial** of $X_0(N)$ written in the coordinate $x$
(W5–W7).  Route (B) produced a new and clean supercongruence, $A_{p-1}\equiv\kappa\psi(p)p
\pmod{p^2}$, which shows that the character $\psi$ is already manufactured by the
Apéry-like recurrence itself (W11), proved for $s_{10}$.  Route (C) is closed as a
*reduction* (no Lucas law lifts the target, W14) but is the source of W5 and W11.
Route (D) turned the missing second-order input into the single clean assertion
"$\eta_1$ is exact", i.e. a $\psi$-twisted Dwork tower (W19, W20), verified for all three
rows including the character row.

---

## 1. Route (A) as posed: why it has no content, and what is true instead

### 1.1 The normalisation error **[proved]**

For a modular function $f$ on $\Gamma_0(N)$, $p\nmid N$, the geometric (degree $p+1$)
Hecke correspondence acts by
$$T_pf(\tau)=\sum_{i=0}^{p-1}f\Bigl(\frac{\tau+i}p\Bigr)+f(p\tau),\qquad
\sum_{i=0}^{p-1}f\Bigl(\frac{\tau+i}p\Bigr)=p\,(U_pf)(\tau),$$
because $\sum_{i}\zeta_p^{ni}=p$ or $0$.  Hence
$$T_pf=p\,U_pf+V_pf .$$
The Kronecker/Igusa congruence $\Phi_p(X,Y)\equiv(X^p-Y)(X-Y^p)\pmod p$ — equivalently
Eichler–Shimura $T_p\equiv \mathrm F+\mathrm F^{t}$ on $X_0(N)_{\mathbf F_p}$, the second
component acting by $0$ on functions because $\mathrm F$ is purely inseparable — says
exactly $T_pf\equiv V_pf\equiv f^p\pmod p$.  Substituting: $p\,U_pf\equiv0\pmod p$.
**Tautology.**

So the task's step "$U_p(x^b)=T_p(x^b)-V_p(x^b)=T_p(x^b)-x^{bp}\bmod p$" is wrong by a
factor $p$: the correct statement is $U_p(x^b)=\bigl(T_p(x^b)-x^{bp}\bigr)/p$, and its
reduction mod $p$ is the **first-order** term of the Kronecker congruence, i.e. it needs
the modular polynomial modulo $p^2$.  That is the weight-$0$ analogue of round 1's verdict
on route (C) at level $Np$, and it is the same obstruction: *the mod-$p^3$ (here mod-$p^2$)
structure of the modular polynomial is as hard as the target.*

For weight $k\ge2$ the same computation gives $T_p=U_p+p^{k-1}V_p$ and therefore
$U_p\equiv T_p\pmod{p^{k-1}}$, which is why the weight-$4$ form $\Phi$ satisfies
$\Phi|U_p\equiv\Phi|T_p\pmod{p^3}$ — but $\Phi$ is meromorphic and $T_p$ does not preserve
its polar divisor, so nothing is gained there either.

### 1.2 The false reduction, recorded **[refuted]** (`42_upx.log`)

Had "$U_p(x^b)\equiv0$" been true, then from $x^{pa+b}\equiv V_p(x^a)x^b$ and
$U_p(V_p(g)f)=g\,U_pf$ one would get $\Xi|U_p\equiv\sum_aL_{pa}x^a$ and hence
$$(S)\bmod p\iff L_{pa}\equiv\psi(p)L_a\iff A_{pa-1}\equiv\psi(p)\,p\,A_{a-1}\pmod{p^2}.$$
Both sides are **false**: `42_upx.log` records nonzero $U_p(x^b)\bmod p$ for essentially
every $b$ once $p\ge11$, and $L_{pa}\equiv\psi(p)L_a$ fails in $36/109$ tests already for
$s_7$, $p=11$ (and in $\ge12/13$ tests for $p\ge59$).  Recorded here because the same
mis-normalisation is easy to repeat.  *(The $a=1$ case of the false statement is
$A_{p-1}\equiv\psi(p)p$; what is actually true is W11, with $\kappa\ne1$.)*

### 1.3 What is true: $U_p(x^b)\bmod p$ is a rational function of $x$ **[verified]** (`42_upx3.log`)

Fitting $U_p(x^b)\bmod p$ against $\mathbf F_p(x)$ by exact linear algebra on $q$-expansions:

| row | $p$ | degree of the fit |
|---|---|---|
| $s_7$ | $3,5,11,13$ | $0,\ 1,\ 3,\ 3$ |
| $s_{10}$ | $3,7,11,13$ | $1,\ 2,\ 3,\ 4$ |

So route (A) *does* have content on the $x$-line — but the input is not the modular
polynomial mod $p$; it is the object identified in §2–§3.

---

## 2. The correct $x$-line formula for $U_p$ **[proved]**

**Proposition 2.1.**  Let $G:=1/(yF)\in\mathbf Z_p[[x]]$ and $H:=(yF)^{p-1}\in\mathbf Z[[x]]$.
For every $f=\sum_{n\ge0}f_nx^n\in\mathbf Z_p[[x]]$, regarded as the $q$-series $f(x(q))$,
$$f\big|U_p\ \equiv\ f_0\ +\ \sum_{n\ge1}\bigl[x^{pn}\bigr]\bigl(f\,H\bigr)\,x^n\pmod p .$$

*Proof.*  $Dx=xyF$ gives $\frac{dq}q=W\,dx$ with $W=\frac1{xyF}=x^{-1}G$.  Round 1
Prop. 4.2 gives $\mathcal C\bigl(f\frac{dq}q\bigr)=(f|U_p)\frac{dq}q$, i.e.
$f|U_p=\frac1W\sum_{j\ge-1}[x^{pj+p-1}](fW)\,x^j$.  Since $G\in\mathbf Z_p[[x]]$,
$G^p\equiv G(x^p)$, so $G=G(x^p)G^{-(p-1)}=G(x^p)H$ and $W\equiv x^{-1}H(x)G(x^p)$.
Then $[x^{pj+p-1}](fW)=[x^{p(j+1)}]\bigl(fH\,G(x^p)\bigr)=\sum_iG_i[x^{p(j+1-i)}](fH)$,
so $\sum_{j\ge-1}[x^{pj+p-1}](fW)x^j=G(x)\bigl(f_0x^{-1}+\sum_{m\ge0}[x^{p(m+1)}](fH)x^m\bigr)$.
Multiplying by $1/W=x/G$ gives the formula. $\square$

**Corollaries.**  $U_p(x^b)=\sum_{n\ge1}[x^{pn-b}]H\,x^n$ ($b\ge1$);
$U_p(1)=1+\sum_{n\ge1}[x^{pn}]H\,x^n$, and since $U_p(1)=1$,
$$\boxed{\ [x^{pn}]\bigl((yF)^{p-1}\bigr)\equiv0\pmod p\quad (n\ge1).\ }$$
**[verified]** (`42_upx4.log`, $p\le29$, $x^{420}$, all rows, $0$ nonzero coefficients);
and the formula for $U_p(x^b)$ was cross-checked against $q$-expansions
(`42_upx4.log`, $s_7$, $p=3,5,11,13$: *FORMULA CONFIRMED*).

**Corollary 2.2 (the target).**  Taking $f=l$ ($l_0=0$),
$$(S)\bmod p\iff \operatorname{comp}_0(l\,H)=\psi(p)\,l\iff [x^{pn}](l\,H)\equiv\psi(p)L_n\ (n\ge1).$$
**[verified, $p\le29$, $x^{420}$, all rows, $0$ failures]** (`42_upx4.log`).

So: *the operator $T(f):=f_0+\sum_{n\ge1}[x^{pn}](fH)x^n$ is $U_p$ in the $x$-coordinate;
it fixes $1$; the missing brick is that $l$ is an eigenvector of $T$ with eigenvalue
$\psi(p)$.*

### 2.1 The $\sqrt P$-free (hyperelliptic) form **[proved]** (`41_hyper.log`)

Independently of §3, one can strip the *transcendental* factor $1/\sqrt P$ out of $\eta$
without any input from the recurrence.  Write $\eta=b\,\frac{dx}{y}$ with
$b:=l/(xF)\in\mathbf Z[[x]]$, $b_0=1$.  Since $y^{p-1}=P^{(p-1)/2}$ and
$\mathcal C(h^p\omega)=h\,\mathcal C(\omega)$,
$$\mathcal C\Bigl(b\,\frac{dx}y\Bigr)=\frac1y\,\mathcal C_x\bigl(b\,P^{(p-1)/2}dx\bigr),$$
hence
$$\boxed{\ (3)\iff [x^{pj+p-1}]\bigl(b\cdot P^{(p-1)/2}\bigr)\equiv\psi(p)\,b_j\pmod p\quad(j\ge0).\ }$$
$P^{(p-1)/2}$ is a **polynomial of degree $p-1$**.  **[verified, every odd $p\le89$, all
three rows, $x^{1200}$: $0$ failures in $400$ tests at $p=3$ down to $13$ at $p=89$]**
(`41_hyper.log`; independently reconfirmed to $p\le31$ in `44_lucas.log`).
This extends round 1's verification of (3) from $p\le101$ (there done with the infinite
series $\sqrt P$) to a *finite-data* form.

**Mod $p^2$ it fails** for essentially every $p\ge5$ (`44_lucas.log`): the finite form is a
genuinely mod-$p$ statement.  **[refuted, $p\le31$]**

**Where $\psi$ is not.**  The $j=0$ slot reads
$\sum_{i=0}^{p-1}b_i\,[x^{p-1-i}]P^{(p-1)/2}\equiv\psi(p)$, and the "pure conic" term is
$$[x^{p-1}]P^{(p-1)/2}=(\lambda_1\lambda_2)^{(p-1)/2}=(B^2-4C)^{(p-1)/2}=\Bigl(\frac{B^2-4C}p\Bigr),$$
which equals $\chi_{-3}(p),\ \chi_{-4}(p),\ \bigl(\tfrac3p\bigr)$ for $s_7,s_{10},s_{18}$ —
**none** of which is $\psi$ (table in `41_hyper.log`).  The character is produced by the
interaction of $b$ with $P^{(p-1)/2}$, not by the branch points alone.

---

## 3. $H$ made explicit: the Lucas property and the supersingular polynomial

### 3.1 $H=P^{(p-1)/2}/F_{<p}$ **[proved from Lucas; verified]** (`46_hexplicit.log`)

The Lucas property $A_{pj+r}\equiv A_jA_r\pmod p$ ($0\le r<p$) holds for all three rows
(**[verified: $0$ failures in $601$ tests, every $p\le23$, $n\le600$]**, `44_lucas.log`; it
is a theorem of Malik–Straub type for Apéry-like sequences, already quoted in
`cooper_sources/REPORT.md` §5.4).  Equivalently
$$F(x)\equiv F_{<p}(x)\,F(x^p)\pmod p,\qquad F_{<p}:=\sum_{n=0}^{p-1}A_nx^n,$$
whence $F^{p-1}=F^p/F\equiv F(x^p)/F(x)\equiv 1/F_{<p}(x)$ and, with $y^{p-1}=P^{(p-1)/2}$,
$$\boxed{\ H=(yF)^{p-1}\ \equiv\ \frac{P(x)^{(p-1)/2}}{F_{<p}(x)}\pmod p\ }$$
**[verified: $F^{p-1}F_{<p}\equiv1$ to $x^{400}$, all rows, $p\le23$]**.

Substituting into Corollary 2.2 gives **W6**, the fully explicit form of the missing brick:
$$[x^{pn}]\Bigl(l(x)\,\frac{P(x)^{(p-1)/2}}{F_{<p}(x)}\Bigr)\equiv\psi(p)\,L_n\pmod p,\qquad n\ge1,$$
**[verified: all $p\le53$ with $p\nmid N$, all three rows, $x^{1200}$; $396$ tests at $p=3$
down to $20$ at $p=53$; $0$ failures]**.  Likewise the proved identity
$[x^{pn}]H\equiv0$ becomes the polynomial identity
$\operatorname{comp}_0\bigl(P^{(p-1)/2}/F_{<p}\bigr)=1$.

This is the sharpest available statement of the gap: *everything except $l$ is now a ratio
of two explicit polynomials of degree $<p$.*

### 3.2 The denominator is the supersingular divisor **[verified, exact]** (`47_bridge.log`)

$F_{<p}\bmod p$ drops degree (its top coefficient $A_{p-1}\equiv0$, W11), and the drop is
exactly what the geometry predicts.  $F^{p-1}=\bigl(F^{(p-1)/2}\bigr)^2$ and
$F^{(p-1)/2}/E_{p-1}$ is a modular *function* mod $p$ (weight $p-1$ over the Hasse
invariant, $E_{p-1}\equiv1$), so the polar divisor of $F^{p-1}$ is twice the supersingular
divisor of $X_0(N)_{\mathbf F_p}$; pushing to the $x$-line (of degree $\deg x=2,4,8$) and
using $\deg\operatorname{div}(E_{p-1})=\mu(N)(p-1)/12$ predicts
$$\deg_xF_{<p}\bmod p\ =\ \Bigl\lfloor \frac{\mu(N)\,(p-1)}{6\deg x}\Bigr\rfloor,\qquad
\mu=8,18,36,\quad \deg x=2,4,8 .$$
**[verified: exact equality for all $p\le79$ with $p\nmid N$, all three rows — $20$ primes
per row, no discrepancy]**.  E.g. $s_7$: $\deg F_{<p}=\lfloor2(p-1)/3\rfloor=6,8,10,12,14,\dots$
at $p=11,13,17,19,23$.

Independently, the modular-function fit $F^{(p-1)/2}=N(u)/D(u)$ over $\mathbf F_p$ succeeds
for $s_7$ at $p=5,7,11,13,17$ with $\deg D=2,0,6,8,10$ (`42_upx2.log`), consistent.

**Consequence.**  $F_{<p}$ is the Hasse invariant of the family (Igusa's theorem: the
truncated period is the Hasse polynomial), and the poles of the kernel $H$ are precisely the
supersingular points.  So (3) is an *ordinarity / unit-root* congruence attached to the
supersingular divisor — which fits round 1's V10 ("not a Dwork factorisation, a genuine
Cartier-eigenvector statement") and explains why no combinatorial Lucas law can deliver it.

---

## 4. Route (C): binomial sums, Lucas, and the new supercongruence

### 4.1 Closed forms **[verified $n\le40$ / search failed]** (`44_binom.log`, `44_s18search.log`)

* $s_{10}$: $A_n=\sum_k\binom nk^4$ — **MATCH**.
* $s_7$: $A_n=\sum_k\binom nk^2\binom{n+k}k\binom{2k}n$ — **MATCH**.
* $s_{18}$: **no formula found.**  Exhaustive search over all products of
  $\binom nk,\binom{2k}k,\binom{n+k}k,\binom{2k}n,\binom{2n-2k}{n-k},\binom{3k}k,\binom n{3k},\binom{n-k}k$
  with exponents $\le5$ times $20$ prefactor variants ($c^k,c^{n-k},c^{n-3k}$, $c\in\{2,3,4\}$,
  with and without $(-1)^k$) produced **no** match to $1,6,54,564,6390,76356,948276$.  The
  best near miss, $\sum_k\binom nk\binom{n+k}k\binom{2n-2k}{n-k}2^k$, agrees for $n\le5$ and
  fails at $n=6$ by $840$.  (Search machinery validated: it rediscovers $\sum\binom nk^4$
  instantly.)

### 4.2 The new supercongruence **[verified $5\le p\le300$; proved for $s_{10}$]** (`43_jslot.log`, `47_bridge.log`)

$$\boxed{\ A_{p-1}\ \equiv\ \kappa\,\psi(p)\,p\pmod{p^2},\qquad
\kappa_{s_7}=-2,\ \kappa_{s_{10}}=-3,\ \kappa_{s_{18}}=-3 .}$$
Verified for **every** prime $5\le p\le300$ in all three rows, with no exception
($p=3$ is exceptional in all rows).  Equivalently
$$\psi(p)\ \equiv\ \kappa^{-1}\,\frac{A_{p-1}}p\pmod p ,$$
so for $s_{18}$ the character $\chi_{-3}$ is manufactured by the Apéry-like recurrence
itself, with no modular or $x$-side input — a purely combinatorial companion to round 1's
V9 ($a_{p-1}\equiv\psi(p)$).  The two are linked by
$$\kappa\,a_{p-1}\ \equiv\ A_{p-1}/p\pmod p\qquad\textbf{[verified, }5\le p\le300\textbf{]}.$$

**Proof for $s_{10}$ ($p\ge5$).**  $A_{p-1}=\sum_{k=0}^{p-1}\binom{p-1}k^4$ and
$\binom{p-1}k\equiv(-1)^k(1-pH_k)\pmod{p^2}$ with $H_k=\sum_{i\le k}1/i$, so
$A_{p-1}\equiv\sum_k(1-4pH_k)=p-4p\sum_{k=0}^{p-1}H_k$.  Now
$\sum_{k=0}^{p-1}H_k=\sum_{i=1}^{p-1}\frac{p-i}i=pH_{p-1}-(p-1)\equiv1\pmod p$
(Wolstenholme: $p\mid H_{p-1}$).  Hence $A_{p-1}\equiv p-4p=-3p\pmod{p^2}$. $\square$
*(The same argument for $s_7$, $s_{18}$ is not available: $s_7$'s sum has the awkward
factor $\binom{2k}n$, and $s_{18}$ has no known binomial form.)*

**Refinement [verified, $p\le59$]:** $A_{p-1}\equiv-2p-6p^2\pmod{p^3}$ for $s_7$ and
$A_{p-1}\equiv-3p-12p^2\pmod{p^3}$ for $s_{10}$.

### 4.3 Lucas census **[verified]** (`44_lucas.log`)

| congruence | $s_7$ | $s_{10}$ | $s_{18}$ |
|---|---|---|---|
| $A_{pj+r}\equiv A_jA_r\ (p)$ | holds ($0/601$) | holds | holds |
| $L_{pj+r}\equiv L_jL_r\ (p)$ | fails | fails | fails |
| $a_{pj+r}\equiv a_ja_r\ (p)$ | fails (except $p\mid N$) | fails | fails (except $p=3$) |
| $b_{pj+r}\equiv b_jb_r\ (p)$ | fails (except $p=2$) | fails | fails (except $p=3$) |
| $A_{mp^s}\equiv A_{mp^{s-1}}\ (p^{3s})$ | holds ($p\ge3$) | holds ($p\ge5$) | **fails for all $p\ge5$** |
| best $e$: $A_{mp}\equiv A_m\ (p^e)$, $m\le40$ | $3$ | $3$ | **$2$** |
| $v_p(A_{p-1})$, $p\le43$ | $1$ (except $p=2$: $2$) | $1$ (except $p=3$: $2$) | $1$ (except $p=3$: $3$) |

Two things follow.  (i) The Lucas property of $A$ is the *only* combinatorial input that
survives, and it is exactly what makes $H$ rational (§3.1) — that is route (C)'s real
contribution.  (ii) $s_{18}$ is genuinely weaker at the Beukers/Dwork level ($p^{2s}$, not
$p^{3s}$), which is consistent with $\psi\ne\mathbf1$ and with round 1's V17.

---

## 5. Route (D): the lift to mod $p^2$

`45_eta1b.log`.  $\mathcal C(\eta)=\psi(p)\eta+p\,\eta_1$ with $\eta_1$ computed exactly
from the integers $a_j$; $h$ from $q(x)^p=q(x^p)(1+p\,h)$ via `serreverse`; for
$\omega=g\,dx$, $D(\omega):=d\bigl(g\cdot xyF\bigr)$ (the coefficient of $\frac{dq}q$).

1. **The round-1 identity holds.**  $\eta_1\equiv D\bigl(\mathcal C(\eta h)\bigr)\pmod p$:
   **[verified]** for $(s_7,p=2,3,5)$, $(s_{10},p=3,7)$, $(s_{18},p=5,7)$, on $193$, $128$,
   $76$, $128$, $53$, $76$, $53$ coefficients respectively.  This closes the round-1 §4.4
   loop numerically.

2. **$\eta_1$ is exact.**  $\mathcal C(\eta_1)=0$ **[verified, $p\le19$, all rows]** —
   note this is *implied* by 1, since a $D$-image is an exact differential and
   $\mathcal C(df)=0$.  Explicitly it says
   $$a_{p^2n-1}\equiv\psi(p)\,a_{pn-1}\pmod{p^2}\qquad\textbf{[verified, }p\le31\textbf{]} .$$

3. **The $\psi$-twisted Dwork tower.**  More generally
   $$a_{pn-1}\equiv\psi(p)\,a_{n-1}\pmod{p^{\,1+v_p(n)}}\qquad(n\ge1)$$
   **[verified, every $p\le31$, all $n\le1500/p$, all three rows, $0$ failures]**.  For
   $\psi=\mathbf1$ this is exactly Dieudonné–Dwork for
   $\mathcal G=\exp\int_0^xa$ (round 1 V19); the verification here shows the *same* tower
   holds verbatim for $s_{18}$ with the character inserted — i.e. the correct
   modularity-free statement covering all three rows is the twisted tower, not the
   integrality of $\mathcal G$ (which indeed fails for $s_{18}$).

4. **No independent description of $\eta_1$.**  $\eta_1/\eta\bmod p$ admits no rational
   representation in $x$ of degree $\le12$ for any $p\le19$ in any row (the sole exceptions
   are the degenerate cells $(s_7,p=3)$, where $\eta_1\equiv0$).  Nor is $\eta_1$ in the
   $\mathbf F_p$-span of $\{\eta,\ x\eta,\ F'/F,\ P'/P,\ xF/l,\ 1\}$: the kernel of the
   $6\times200$ system is $0$ for every row and every $p\ge5$ tested (the one relation found,
   $s_7$ at $p=5$, is $4\,x\eta+F'/F=0$, a degeneracy among the basis vectors that does not
   involve $\eta_1$).  **This line is a dead end as posed.**

---

## 6. Where the gap now stands

**Proved here.**  W1 (route (A) as posed is empty), W3, W4 (the $x$-line formula for $U_p$
and the reformulation of the target), the equivalence in W8, W5 modulo the Lucas property,
W19 given W18, and W11 for $s_{10}$.

**Reduced to.**  The smallest missing brick is now, equivalently:

* **(3-x)** $\ [x^{pn}]\bigl(l\cdot P^{(p-1)/2}/F_{<p}\bigr)\equiv\psi(p)L_n\pmod p$, $n\ge1$
  — a congruence whose kernel is a ratio of two polynomials of degree $<p$, the denominator
  being the supersingular polynomial of $X_0(N)$ on the $x$-line;
* **(3-b)** $\ [x^{pj+p-1}]\bigl(b\cdot P^{(p-1)/2}\bigr)\equiv\psi(p)b_j\pmod p$,
  $b=l/(xF)$ — the same without the Lucas input;
* **(3-T)** $\ l$ is an eigenvector, with eigenvalue $\psi(p)$, of the $\mathbf F_p$-linear
  operator $T(f)=f_0+\operatorname{comp}_0(fH)-(fH)_0$ which is $U_p$ in the $x$-coordinate
  and which fixes $1$.

**Still open, in order of size.**

1. The eigenvector statement (3-T).  A partial-fraction attack looks natural: for
   $\alpha\in\overline{\mathbf F}_p$, $\operatorname{comp}_0\bigl(1/(x-\alpha)\bigr)
   =\alpha^{p-1}/(x-\alpha^p)$, so on the rational part of $H$ the operator permutes the
   poles by Frobenius with multiplier $\alpha^{p-1}$ — and the poles are the supersingular
   points, which live in $\mathbf F_{p^2}$.  That is a plausible source for a $\pm1$
   eigenvalue.  What blocks it is that $l$ is *not* rational, so the pairing of $l$ against
   the partial fractions is not finite.  This is the concrete next step.
2. $A_{p-1}\equiv\kappa\psi(p)p\pmod{p^2}$ for $s_7$ and $s_{18}$ (W11).  For $s_7$ this is
   an elementary-looking binomial identity; for $s_{18}$ one first needs a closed form
   (W17 says none of the standard shapes works).
3. The mod-$p^2$ lift: given (3), the remaining input is $\eta_1\equiv D(\mathcal C(\eta h))$
   (W18, verified) — for which route (D) offers no structural handle (W21).  The cleanest
   equivalent target is the twisted Dwork tower W20, which packages the whole $p$-adic
   family into one statement.

**Routes now closed.**  Route (A) *as posed* (W1).  Any attempt to derive (3) from a Lucas
law for $A$, $L$, $a$ or $b$ alone (W14, W15).  Any attempt to read $\psi$ off the branch
points of $P$ (W10).  A description of $\eta_1$ inside the obvious $6$-dimensional space
(W21).

---

## 7. Files

| file | contents |
|---|---|
| `40_core.gp` | shared library: `AVEC/FSER/LSER/PSER/ETAS/AJ` on top of `lib.gp` |
| `40_check.gp/.log` | sanity: $A_n$, $a_j$ reproduce round 1 §4.1; $A_n=[x^n]F$ vs $q$-expansions |
| `41_hyper.gp/.log` | the $\sqrt P$-free form $[x^{pj+p-1}](bP^{(p-1)/2})\equiv\psi(p)b_j$, $p\le89$; the conic term |
| `42_upx.gp/.log` | the refuted claims $U_p(x^b)\equiv0$, $L_{pa}\equiv\psi(p)L_a$, $A_{pa-1}\equiv\psi(p)pA_{a-1}$; the discovery of $A_{p-1}\bmod p^2$ |
| `42_upx2.gp/.log` | $F^{(p-1)/2}=N(u)/D(u)$ over $\mathbf F_p$ (Hasse-invariant structure).  *(The second block of this script aborted on a series-precision error; it is superseded by `42_upx3.gp`.)* |
| `42_upx3.gp/.log` | $U_p(x^b)\bmod p$ is a rational function of $x$; degrees |
| `42_upx4.gp/.log` | the proved formula $U_p(x^b)=x\sum_jH_{pj+p-b}x^j$, cross-checked; $[x^{pn}]H=0$; $\operatorname{comp}_0(lH)=\psi(p)l$; $H$ rational |
| `43_jslot.gp/.log` | route (B): $a_{p-1}\bmod p,p^2$; the defect $d(p)$; $j=1,2$ slots; $A_{p-1}\bmod p^2,p^3$ |
| `44_binom.gp/.log`, `44_s18search.gp/.log`, `44_s18near.gp/.log`, `44_lucas.gp/.log` | route (C): closed forms, the failed $s_{18}$ search, the Lucas census |
| `45_eta1.gp` | first attempt at route (D); aborted on series-precision errors, **superseded by `45_eta1b.gp`** (kept only because it was swept into another agent's commit) |
| `45_eta1b.gp/.log` | route (D): $\eta_1=D(\mathcal C(\eta h))$; $\mathcal C(\eta_1)=0$; the negative lindep tests |
| `46_hexplicit.gp/.log` | $F^{p-1}=1/F_{<p}$; the explicit target (3-x); the twisted Dwork tower; $\mathcal C(\eta_1)=0$ |
| `47_bridge.gp/.log` | $\psi(p)=\kappa^{-1}A_{p-1}/p$; $\kappa a_{p-1}=A_{p-1}/p$; $\deg F_{<p}$ vs the supersingular count |

Run with `gp -q <file>`; total runtime a few minutes.
