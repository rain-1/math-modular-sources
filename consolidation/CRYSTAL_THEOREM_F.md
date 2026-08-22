# Theorem F and the tower trichotomy in Dwork-crystal language

*Fable (Opus 5), 2026-08-22.  Recasts `paper/sections/04_padic_euler_draft.tex`
(Theorem F) and `consolidation/GOOD_PRIME_TOWERS.md` (tower trichotomy) so that
the statements no longer mention modular forms.  New computations:
`lattice/crystal_thmF/` (`sep.gp`, `dwork.gp`, `towers2.gp`, `slopes.gp`; logs
alongside).  Tags: [proved] [verified] [lit] [conj] [open].*

---

## 0. Verdict, first

| claim | verdict |
|---|---|
| The pure half of the picture ($a_n$ integral, Dwork congruences, tower eigenvalue $1$ on the $A$-row) is a **theorem for every row with a Laurent-polynomial model**, at every $p$ including $p=2$ | **[proved, lit]** Mellit–Vlasenko / Samol–van Straten; Beukers–Vlasenko I Thm 5.3, Cor 5.9 |
| The Frobenius matrix at a MUM point is triangular with diagonal $(1,p,\dots,p^{r-1})$, and its off-diagonal entries are $p$-adic zeta values | **[proved, lit]** Beukers–Vlasenko, *Frobenius structure and $p$-adic zeta values* (arXiv:2302.09603), Prop. 1.3, Thms 1.4–1.5 |
| The existence criterion of Theorem F is **not** the Euler factor but a *Frobenius-slope separation* $\sigma_p>0$, and $\sigma_p=v_p(c)$ for every (R2)/(R3) row | **[proved]** (exact Casoratian) + **[verified]** 84 cells, `sep.gp` |
| The trichotomy is a statement about the *extension* crystal, with eigenvalue $\varepsilon(p)p^{-(w+1)}$ on the extension direction, $w+1$ = Hodge depth (denominator exponent), **not** the modular weight | **[verified]** 144 old cells + **20 new non-modular cells**, `towers2.gp` |
| **Zudilin's Catalan row has a Dwork crystal.** $R_m:=16^mQ_m\in\mathbf Z$ and satisfies the full Dwork tower $R_{mp^s}\equiv R_{mp^{s-1}}\pmod{p^s}$ — 2288 tests, $p=2,3,5,7$, 0 failures | **[verified]**, `dwork.gp` — *new*, and it contradicts the reading of `SOURCES_S18_ZUDILIN.md` §5 that the row has no geometric model |
| Zudilin's good-prime tower eigenvalue is $\chi_{-4}(p)\,p^{-2}$ **exactly** ($p=3,5,7,11$) — the character of Catalan's constant $G=L(2,\chi_{-4})$, read off a non-modular row | **[verified]**, `towers2.gp` — *new* |
| Brown–Zudilin's $\zeta(5)$ row also has the Dwork property ($Q_n$: 1716 tests, 0 failures) and has tower eigenvalue $p^{-5}$ **exactly** at $p=2,3,5,7$ (and $p^{-3}$ on the $\zeta(3)$ direction), trivial character | **[verified]**, `towers2.gp` — *new*; refines the "no tower limit" of `ZETA5_TWO_ROW.md` §4 |
| The crystal argument reproduces the **exact** Zudilin tail $v_2(\zeta_2(2)-P_m/Q_m)=8m-1-4s_2(m)$ | **[proved]** from the Casoratian + the denominator law, §5.3; re-verified to $m=200$ |
| The crystal argument produces the **value** $\xi_2^{\rm Zud}=\zeta_2(2)$ | **[open]** — the value is the syntomic regulator; crystal theory gives existence and rate, not the $L$-value.  Beukers' Padé proof (`ZUDILIN_2ADIC.md`) remains the only proof |
| The cuspidal row is **not** a Dwork row: $a_n$ fails $a_{mp^s}\equiv a_{mp^{s-1}}$ (531 failures, first at $p=2,s=2,m=1$) yet still has $\xi_2=0$ | **[verified]** — the two hypotheses are independent |

The one-line summary: **Theorem F is not about Eisenstein series; it is about the
$p$-adic separation of the Frobenius slopes of an extension of a MUM crystal by
$\mathbf Q_p(0)$.  The Eisenstein/Euler-factor formulation is the modular
*computation* of that separation, and it computes the *value*, which crystal
theory does not.**

---

## 1. The crystal set-up (what Beukers–Vlasenko actually give)

Throughout $R$ is a $p$-adically complete characteristic-zero ring with a
Frobenius lift $\sigma$, $\Lambda\in R[x_1^{\pm},\dots,x_n^{\pm}]$ a Laurent
polynomial with Newton polytope $\Delta$, and $A_m=\mathrm{CT}(\Lambda^m)$.

**(1.1) The module and the Cartier operator.** [lit: BV-I §2–3]
$\Omega_\Lambda(\mu)$ is spanned by $\omega_u=(u_0-1)!\,x^u/\Lambda^{u_0}$ over
the integral points of the cone $C(\Delta)$ supported in an *open* $\mu\subseteq\Delta$
(complement of a union of faces; $\mu=\Delta^\circ$ is the smallest).  The
Cartier operator is $C_p(\sum a_kx^k)=\sum a_{pk}x^k$; it satisfies
$C_p\theta_i=p\,\theta_iC_p$ and $C_p(g(x^p)h)=g(x)C_p(h)$
[BV-I Lemma 3.2].  **Caution:** $C_p(\Omega_\Lambda)\subseteq\widehat\Omega_{\Lambda^\sigma}$
is proved only for $p>2$ [BV-I Prop. 3.3].

**(1.2) Hasse–Witt and the unit-root crystal.** [lit: BV-I Def. 3.5, Thm 4.3, Rmk 4.4, Def. 4.9]
$(\beta_m)_{u,v}=$ coefficient of $x^{mv-u}$ in $\Lambda^{m-1}$; $\beta_p(\mu)$ is
the Hasse–Witt matrix.  *If $\beta_p(\mu)$ is invertible* then
$Q_\Lambda(\mu)=\widehat\Omega_\Lambda(\mu)/U_\Lambda(\mu)$ is free of rank
$h=\#\mu_{\mathbf Z}$, with $U_\Lambda(\mu)=\{\omega:C_p^s\omega\equiv0\ (p^s)\ \forall s\}$
the largest subcrystal on which $C_p$ is divisible by $p$; $C_p$ is invertible on
$Q_\Lambda(\mu)$, which is therefore the **unit-root quotient**.  Ordinarity is
*exactly* invertibility of $\beta_p(\mu)$ mod $p$ — there is no other hypothesis.

**(1.3) Dwork's congruences.** [lit: BV-I Thm 5.3, Thm 5.7; BV-II Thm 3.2, Cor. 3.1;
Mellit–Vlasenko IJNT 12 (2016) Thm 1; Samol–van Straten]
With $\Lambda_\sigma$ the matrix of $C_p$ on $Q_\Lambda(\mu)$,
$$\beta_{mp^s}(\mu)\equiv\Lambda_\sigma\,\beta^\sigma_{mp^{s-1}}(\mu)\pmod{p^s},
\qquad \Lambda_\sigma=\lim_s\beta_{p^s}\,\sigma(\beta_{p^{s-1}})^{-1},$$
and on expansion coefficients $a_{p^sk}\equiv\Lambda_\sigma a^\sigma_{p^{s-1}k}\pmod{p^s}$.
In the rank-one case ($0$ the unique interior lattice point of a reflexive
$\Delta$, $\mu=\Delta^\circ$, $h=1$) this is the scalar statement
$$f(t)/f(t^p)\equiv f_{p^s}(t)/f_{p^{s-1}}(t^p)\pmod{p^s},\qquad
f=\sum_mA_mt^m,$$
[BV-II Cor. 3.1 + Thm 3.2], and its coefficient shadow is the **Gauss/Dwork
tower** $A_{mp^s}\equiv A_{mp^{s-1}}\pmod{p^s}$ [BV-I Cor. 5.9].
*Mellit–Vlasenko's version has no hypothesis on $p$: $p=2$ is allowed.*

**(1.4) The MUM Frobenius matrix.** [lit: BV, arXiv:2302.09603 Def. 1.1, Prop. 1.3]
Let $L$ be MUM at $t=0$ of rank $r$, with the standard basis
$y_i=\sum_{j\le i}F_j(t)\log^{i-j}t/(i-j)!$, $F_0=f$.  A $p$-adic Frobenius
structure is $A=\sum_jA_j\theta^j$ with $A_0(0)=1$ such that $A(y_i(t^p))$ again
solves $L$.  **Proposition 1.3.**
$$A\bigl(y_i(t^p)\bigr)=p^i\sum_{j\le i}\alpha_j\,y_{i-j}(t),\qquad
\alpha_j=A_j(0)\in\mathbf Z_p,\ \alpha_0=1 .$$
So the Frobenius on the *pure* crystal at the MUM point is
$$\Phi_{\rm pure}=\begin{pmatrix}1&\ast&\ast&\cdots\\0&p&\ast&\\0&0&p^2&\\ &&&\ddots\end{pmatrix},$$
**upper triangular, unipotent-times-diagonal $(1,p,\dots,p^{r-1})$**, and the
off-diagonal entries $\alpha_j$ are, for the simplicial and hyperoctahedral
families, explicit generating series in $\Gamma_p$: $\alpha_j=[x^j]\exp(-\sum_{m\ge2}\zeta_p(m)x^m/m)$
[Thm 1.5], so $\alpha_1=\alpha_2=0$, $\alpha_3=-\zeta_p(3)/3$, $\alpha_5=-\zeta_p(5)/5$,
… .  Only **odd** $\zeta_p$ occur, because $\zeta_p$ vanishes at positive even
integers.  (Same structure conjectured by Candelas–de la Ossa–van Straten,
arXiv:2104.07816 §4.4, with $\gamma=(\chi(X)/y)\zeta_p(3)$.)

**(1.5) The extension, and where the row's $b_n$ lives.** [proved, elementary]
For (R2)/(R3) the two solution sequences of the recurrence give $A(t)=\sum a_nt^n$
with $LA=0$ and $B(t)=\sum b_nt^n$ with
$$L\,B=\rho\qquad(\rho\ \text{a nonzero constant, from the $n=0$ defect } b_1\ne \tfrac{b}{1}b_0).$$
So **$B$ is not a second solution of $L$: it is the period of an extension**
$$0\longrightarrow \mathcal M\longrightarrow\widetilde{\mathcal E}\longrightarrow\mathbf Q_p(0)\longrightarrow0$$
of the MUM crystal $\mathcal M$ by the trivial crystal, the class being $\rho$.
$B$ is pinned by $B(0)=0$ up to adding multiples of $f$ — and *that* ambiguity is
the whole subject: $\xi_p$ is the canonical choice.

Because $B$ carries no logarithm, Frobenius on $\widetilde{\mathcal E}$ in the basis
$(y_0,\dots,y_{r-1},B)$ has the shape
$$\widetilde\Phi=\begin{pmatrix}\Phi_{\rm pure}&\ \widetilde\xi\,e_1\\0&\varepsilon\end{pmatrix},
\qquad \Phi(B)=\varepsilon B+\widetilde\xi\,y_0 ,$$
the $y_1,\dots,y_{r-1}$ components of $\Phi(B)$ vanishing because $\Phi(B)$ is
log-free.  **Splitting criterion.** $\widetilde\xi$ can be removed by
$B\mapsto B-cy_0$ iff $\varepsilon\ne1$ (solve $c(\varepsilon-1)=\widetilde\xi$).
Hence:

> **the extension splits $p$-adically iff the Frobenius eigenvalue on the trivial
> quotient differs from the unit root $1$ of $\mathcal M$; when they coincide the
> Frobenius is unipotent and $\widetilde\xi$ is a genuine invariant — the
> $p$-adic period of the extension.**

This is the exact analogue, one step up, of the structure BV compute for the pure
crystal in Dwork crystals III §7 (Prop. 7.7, Lemma 7.13): there the off-diagonal
$\lambda_1=\frac{F(t)F(t^\sigma)}{W(t^\sigma)}\log\bigl(\gamma^{p-1}q^p/q^\sigma\bigr)$
vanishes exactly for the *excellent* Frobenius lift $q^\sigma=\gamma^{p-1}q^p$
(Def. 7.14), $q=t\exp(G/F)$ the canonical coordinate.

---

## 2. The measured eigenvalues, and what they are

The project measures $U_p$ on **coefficient sequences**, $(U_pX)(n)=X(pn)$, which
is dual to Frobenius on solutions; eigenvalues are inverted.  Write
$\rho^A_s=A(ap^{s+1})/A(ap^s)$, $\rho^B_s=B(ap^{s+1})/B(ap^s)$.

**(2.1) $\rho^A_s\to1$ is Dwork's congruence.** [proved for Laurent rows]
$A_{mp^{s+1}}\equiv A_{mp^s}\pmod{p^{s+1}}$ is exactly §1.3, hence holds for
*every* row with a Laurent-polynomial model **including $p=2$** (Mellit–Vlasenko).
This upgrades the "verified, $p=2,3$ outside every quoted theorem" caveat of
`GOOD_PRIME_TOWERS.md` §9 to a theorem for all rows in the census (all of which
are constant-term rows: `NOME_INTEGRALITY.md` §8, Straub).  Verified again here:
`dwork.gp`, 5737 tests per row, 0 failures, rows C, E, F, $\gamma$, $\alpha$, $s_{18}$.

**(2.2) The $B$-row eigenvalue is $\varepsilon(p)\,p^{-(w+1)}$, $w+1$ = Hodge depth.**
[verified: 144 old + 20 new cells]  Here $w+1$ is the *denominator exponent*
($d_n^{w+1}b_n\in\mathbf Z$), i.e. the number of integrations $D^{-(w+1)}$, i.e.
the Tate twist in $0\to\mathcal V(w+1)\to\mathcal E\to\mathbf Q(0)\to0$ —
**not** the rank of $L$.  New evidence, both rows non-modular:

| row | $w+1$ | $p$ | measured $\lim\rho^B_s\cdot p^{w+1}$ | $\varepsilon(p)$ |
|---|---|---|---|---|
| Zudilin Catalan | $2$ | $3,7,11$ | $-1$ (digitwise) | $\chi_{-4}(p)=-1$ |
| Zudilin Catalan | $2$ | $5$ | $+1$ | $\chi_{-4}(5)=+1$ |
| Brown–Zudilin $\zeta(5)$ | $5$ | $2,3,5,7$ | $+1$ | $\mathbf 1$ |
| Brown–Zudilin, $\widehat P$ ($\zeta(3)$ direction) | $3$ | $2,3$ | $+1$ | $\mathbf 1$ |

Convergence is linear in $s$ (e.g. BZ at $p=2$, $a=1$: $v_2(2^5\rho^B_s-1)=2,5,8,\dots,26$).
So **the character in the tower eigenvalue is the character of the $L$-value the
row computes**: $\chi_{-4}$ for Catalan's $G=L(2,\chi_{-4})$, trivial for
$\zeta(5)$ and $\zeta(3)$.  This is the first time the trichotomy has been seen
outside the modular census, and it is seen on a hypergeometric row and on a
cellular ($\mathcal M_{0,n}$) row.

*(Correction to `ZETA5_TWO_ROW.md` §4: "There is no tower limit" for Brown–Zudilin
is true of $P_n/Q_n$ but false of the correctly normalised tower — the rescaled
$\Lambda_a=\lim_sp^{5s}P(ap^s)/Q(ap^s)$ exists at $2,3,5,7$.  BZ is a **good
prime row at every prime tested**, which is precisely why it has no $p$-adic
Apéry limit anywhere.)*

**(2.3) Conjecture: $\varepsilon(p)=\gamma^{p-1}$.** [conj]
If $\gamma$ generates a quadratic field with $\gamma^2\in\mathbf Q$ then
$\gamma^{p-1}=(\gamma^2\!\mid\!p)$ is a quadratic character value.  BV-III
Def. 7.14 normalises the excellent Frobenius lift by $q^\sigma=\gamma^{p-1}q^p$
in the canonical coordinate.  This predicts $\varepsilon=\gamma^{p-1}$ with
$\gamma$ the normalisation constant relating the row's coordinate $t$ to the
canonical coordinate $q$.  Consistent with everything measured: $\varepsilon=\mathbf 1$
exactly for the rows with $\gamma\in\mathbf Q$ (Apéry $\zeta(3)$, whose mirror map
is the eta quotient $t(q)=(\eta_1\eta_6/\eta_2\eta_3)^{12}$ with rational
coefficients — `NOME_INTEGRALITY.md` §9 — and Brown–Zudilin), and $\varepsilon=\chi_{-4}$
for the Catalan row.  Not tested directly; the test is to compute $q^\sigma/q^p$
for one row at one prime.

---

## 3. The real existence criterion: Frobenius-slope separation

**Definition.** For a row put
$$\mu_1:=\lim_n\frac{v_p(a_n)}n,\qquad
\mu_{\rm det}:=\lim_n\frac{v_p\bigl(a_{n+1}b_n-a_nb_{n+1}\bigr)}n,\qquad
\boxed{\ \sigma_p:=\mu_{\rm det}-2\mu_1\ }$$
$\mu_1$ is the $p$-adic growth exponent of the holomorphic period (the unit-root
slope, $=0$ iff the row is $p$-integral with infinitely many $p$-unit terms);
$\mu_{\rm det}$ is the slope of the determinant of the crystal (the Wronskian);
$\sigma_p$ is the **gap between the two Frobenius slopes**.

**Proposition 3.1** [proved].  $v_p(b_{n+1}/a_{n+1}-b_n/a_n)=n\sigma_p+O(\log n)$.
Hence $b_n/a_n$ converges in $\mathbf Q_p$ iff $\sigma_p>0$, and then
$v_p(\xi_p-b_n/a_n)=n\sigma_p+O(\log n)$.
*Proof.* $b_{n+1}/a_{n+1}-b_n/a_n=-\mathrm{Cas}_n/(a_na_{n+1})$ and take valuations.
(For (R2)/(R3) $\mathrm{Cas}_n=-c^n/(n+1)^{2\ \rm or\ 3}$ exactly, `04_padic_euler_draft.tex`
Prop. `prop:dsplit`.)  $\square$

**Corollary 3.2** [proved for the ten (R2)/(R3) rows].  $\sigma_p=v_p(c)$.
Hence *the $p$-adic Apéry limit exists iff $p\mid c$*, with no reference to any
modular source.  Hypothesis (c$'$) of Theorem F is exactly $\sigma_p>0$.

**Verification** (`sep.gp`, $N=1500$, 84 cells = 14 rows $\times$ 6 primes):
$\sigma_p$ measured to $O(\log N/N)\approx0.01$ agrees with $v_p(c)$ in all
72 cells of the twelve normalised rows; the two rows whose recurrence is not in
(R2)/(R3) normalisation ($s_{18}$, cusp row) have a different Casoratian and
measured $\sigma_2=0$, $\sigma_2\approx2$ respectively — both matching the known
limit/no-limit verdicts of `EULER_CRITERION.md` §4.  Selected rows:

| row | $p$ | $\sigma_p$ (measured) | $v_p(c)$ |
|---|---|---|---|
| A $(7,2,-8)$ | 2 | $2.976$ | 3 |
| C $(10,3,9)$ | 3 | $1.993$ | 2 |
| E $(12,4,32)$ | 2 | $4.979$ | 5 |
| $\alpha$ Domb | 2 | $5.965$ | 6 |
| $\delta$ $(7,3,81)$ | 3 | $3.987$ | 4 |
| $\eta$ $(11,5,125)$ | 5 | $2.991$ | 3 |
| $\gamma$ (Apéry $\zeta(3)$) | all | $\le0.006$ | 0 |
| Zudilin Catalan | 2 | $7.96$ | — ($\mu_1=-3.98$) |

**Geometric reading.** $\sigma_p=-\max\{v_p(s):s\ \text{singular},\ s\ne0,\ |s|_p\ne1\}$:
the singular points of $L$ other than the $p$-adically dominant one lie *outside*
the closed unit disc.  Three regimes, which is the trichotomy:

| regime | singular locus vs. the MUM residue disc | consequence |
|---|---|---|
| $\sigma_p>0$ (**slope prime**) | all other singular fibres have $|t|_p>1$; the crystal is Frobenius-trivialisable on a strict neighbourhood of $\{|t|_p\le1\}$ | $\widetilde\Phi$ unipotent, $\xi_p=\widetilde\xi$ exists |
| $\sigma_p=0$ (**good prime**) | a singular fibre sits *on* $|t|_p=1$ (in the modular case: the supersingular discs, `THEORY_NOTES_05` §4) | $\widetilde\Phi$ has distinct eigenvalues $1$ and $\varepsilon p^{-(w+1)}$; extension splits; only rescaled tower limits $\Lambda_a$ survive |
| $\sigma_p<0$ or the slopes coincide irrationally | singular fibres *inside* $|t|_p<1$ | nothing: no limit, and no tower limit either |

Brown–Zudilin at $p=2$ was thought to be the third case (`ZETA5_TWO_ROW.md`: the
three roots of $4\lambda^3-2368\lambda^2-188\lambda+1$ all have $v_2=-\tfrac23$
and the cubic is irreducible over $\mathbf Q_2$).  The measurement of §2.2 shows
it is really the **second**: $Q_n$ is $2$-adically a unit, $\mu_1=0$, and the
tower is perfectly ordinary with eigenvalues $(1,2^{-5})$.  The irreducible
Newton polygon governs the *other* solutions, not the row.

---

## 4. Theorem F$'$

> **Theorem F$'$ (crystal form).** [status: (i) proved, (ii) proved, (iii) open — see §6]
> Let $L$ be a differential operator over $\mathbf Q$, MUM at $t=0$, with
> holomorphic solution $f=\sum_na_nt^n$, and let $B=\sum_nb_nt^n$ satisfy
> $LB=\rho\ne0$, $B(0)=0$; let $\widetilde{\mathcal E}$ be the corresponding
> extension of the Picard–Fuchs crystal $\mathcal M$ by $\mathbf Q_p(0)$.
> Fix a prime $p$ and assume
>
> **(F1) (crystal)** $L$ carries a $p$-adic Frobenius structure on a disc
> $|t|_p<\rho_0$ with $\rho_0>1$ — e.g. $a_n=\mathrm{CT}(\Lambda^n)$ for a Laurent
> polynomial with $0$ the unique interior lattice point of $\Delta$ and invertible
> Hasse–Witt matrix (BV-I Thm 4.3, BV-II Cor. 3.1);
>
> **(F2) (separation)** $\sigma_p>0$, i.e. the two Frobenius slopes of
> $\widetilde{\mathcal E}$ at the MUM point are distinct and the unit-root slope
> is attained by $a_n$ (equivalently: every singular fibre other than the
> $p$-adically dominant one lies at $|t|_p>1$).
>
> Then:
> **(i)** $\rho^A_s\to1$ and $\rho^B_s/\rho^A_s\to1$: the Frobenius on
> $\widetilde{\mathcal E}$ is **unipotent**;
> **(ii)** $\xi_p=\lim_nb_n/a_n$ exists in $\mathbf Q_p$, with
> $v_p(\xi_p-b_n/a_n)=n\sigma_p+O(\log n)$, and $\xi_p$ is the unique constant for
> which $B-\xi_pf$ converges on $|t|_p\le p^{\sigma_p}$; equivalently $\xi_p$ is
> the off-diagonal Frobenius entry $\widetilde\xi$ of $\widetilde{\mathcal E}$,
> the obstruction to splitting;
> **(iii)** *(not supplied by crystal theory)* the identification of $\widetilde\xi$
> with an $L$-value.
>
> If $\sigma_p=0$ the extension splits, $\xi_p$ does not exist, and the invariants
> are the rescaled tower limits $\Lambda_a=\lim_s\varepsilon(p)^sp^{(w+1)s}b_{ap^s}/a_{ap^s}$,
> the coordinates of the canonical Frobenius splitting.

**What replaces the modular hypotheses.**

| Theorem F (modular) | Theorem F$'$ (crystal) | status of the translation |
|---|---|---|
| (a) $\mathcal E_p(s)\mid P(s)$ | (F2) $\sigma_p>0$ | equivalent on the census: $\mathcal E_p\mid P\iff p\mid c\iff\sigma_p>0$, 36+72 tests, no exception — **the one unproved link** (§5.1) |
| (b) descent to $\mathbb P^1_t$, $\Gamma_t$-invariance | *absorbed*: the crystal lives on the $t$-line already | Theorem `thm:descentauto` becomes unnecessary |
| (c)/(c$'$) integral hauptmodul; $\rho=p^{\sigma_p}>1$ | (F2) again | identical content |
| (d) $v_p(a_n)=O(\log n)$ | (F1) via Dwork: $a_n\in\mathbf Z_p$ and $\mu_1=0$ | **upgraded from open to proved** for Laurent rows (Mellit–Vlasenko), the item `rem:remains`(1) |
| $\xi_p=-Q(w+1)\kappa_p$ | $\xi_p=\widetilde\xi$ | crystal gives existence and rate; the *value* stays modular/syntomic |

So Theorem F$'$ **proves more rows and fewer things**: it applies to every row
with a Laurent-polynomial model (hence to Zudilin, Cooper's $s_{18}$, Brown–Zudilin,
all the sporadic and CY rows), and it removes hypotheses (b), (c), (d); it does not
produce the Kubota–Leopoldt value.

---

## 5. Specialisations

### 5.1 The Euler-factor criterion in crystal terms

**Proposition 5.1** [proved except for one link].  For a modular row with source
$\Phi=P(V)E^{\psi,\varphi}_{w+2}$ the following are equivalent:
1. $1-\psi(p)p^{-s}$ divides $P(s)$ in $\mathbf Q[p^{-s}]$;  *(Theorem F (a))*
2. $p\mid c$;
3. $\sigma_p>0$;
4. $\lim_nb_n/a_n$ exists in $\mathbf Q_p$;
5. the Frobenius of $\widetilde{\mathcal E}$ at the MUM point is unipotent, i.e.
   $\rho^B_s/\rho^A_s\to1$ (equivalently $\lim\rho^B_s\in\{1,0\}$ rather than
   $\varepsilon(p)p^{-(w+1)}$).

$(2)\Leftrightarrow(3)$ is Cor. 3.2 [proved]; $(3)\Leftrightarrow(4)$ is
Prop. 3.1 [proved]; $(4)\Rightarrow(5)$ is immediate ($b_n=\xi_pa_n+O(\rho^{-n})$
with $\xi_p\ne0$ gives $\rho^B_s/\rho^A_s\to1$; $\xi_p=0$ gives $\lim\rho^B_s=0$,
the third branch), and $(5)\Rightarrow(4)$ under (F1).  **$(1)\Leftrightarrow(2)$
is verified (36 tests, `04_padic_euler_draft.tex` Remark, plus the 72 cells of
`sep.gp`) and not proved** — it is the bridge between the *source* (a statement
about the oldform polynomial) and the *host geometry* (a statement about where the
conifold points sit $p$-adically).  This is the sharpest form in which the Euler
criterion is still open, and it is now a clean statement about one modular curve
and its hauptmodul, with no $p$-adic analysis in it.

The Tate-curve reading of the eigenvalues is now literal: at $t=0$ the fibre is
nodal, $\Phi_{\rm pure}$ has diagonal $(1,p,\dots,p^{r-1})$ [BV 2302.09603
Prop. 1.3], the extension direction carries the Tate twist $(w+1)$, and the
quadratic twist $\varepsilon(p)=\psi(p)$ is the $\gamma^{p-1}$ of the excellent
Frobenius lift (§2.3).

### 5.2 Zudilin's Catalan row **has** a Dwork crystal

`SOURCES_S18_ZUDILIN.md` §5 concluded, correctly, that the row has no *modular*
parametrisation (order-4 operator, non-integral mirror map in the given
coordinate).  That conclusion does **not** survive rescaling.  With
$$R_m:=16^m\,Q_m:\qquad R_0,\dots,R_6=1,\ 28,\ 2596,\ 311536,\ 41759524,\ 5953698928,\ 882922368784,$$
`dwork.gp` finds
* $R_m\in\mathbf Z$ for all $m\le1200$ (no exception), and $v_2(R_m)=2s_2(m)$ **exactly** for $m\le12$ — the proved denominator law $v_2(Q_m)=-4m+2s_2(m)$ (`ZUDILIN_2ADIC.md` (5.4)) in normalised form;
* the full Dwork tower $R_{mp^s}\equiv R_{mp^{s-1}}\pmod{p^s}$ at $p=2,3,5,7$: **2288 instances, 0 failures.**

By BV-I Cor. 5.9 / Mellit–Vlasenko these congruences are exactly what a
constant-term sequence of a Laurent polynomial with $\mathbf 0$ the unique
interior lattice point produces.  So the row satisfies the *observable* content
of (F1); what is missing is the polynomial $\Lambda$ itself and a proof of the
Frobenius structure.  **This is a concrete, well-posed open problem: find
$\Lambda$ with $\mathrm{CT}(\Lambda^m)=16^mQ_m$.**  (The row is the moving-point
Padé specialisation $Q_m=q_m(\tfrac12-m)$, so a natural guess is a
one-variable-larger diagonal.)

**Where the $2$-adic limit comes from, in crystal terms.**  $\mu_1=-4$ exactly
(so after rescaling $\mu_1=0$: the row is *ordinary at $2$ after normalisation*),
$\mu_{\rm det}=0$, hence
$$\sigma_2=0-2(-4)=8>0,$$
measured $7.96$ at $N=1500$.  So (F2) holds at $p=2$ and fails at $3,5,7,11,13$
($\sigma_p=0$) — matching exactly the good/slope split of `EULER_CRITERION.md`.
Theorem F$'$(ii) therefore gives the existence of $\xi_2^{\rm Zud}$ **without any
modular input and without Beukers' Padé theory**.

### 5.3 The exact Zudilin tail, from the crystal side

**Proposition 5.2** [proved].  With $\mathrm{Cas}_m=(-1)^{m-1}(20m^2-8m+1)/(8m^2(2m-1)^2)$
(`ZUDILIN_2ADIC.md` (4.1)) and $v_2(Q_m)=-4m+2s_2(m)$ (ibid. (5.4)),
$$v_2\Bigl(\frac{P_{m+1}}{Q_{m+1}}-\frac{P_m}{Q_m}\Bigr)
= v_2(\mathrm{Cas}_{m+1})-v_2(Q_{m+1})-v_2(Q_m)
= 8m+3-4s_2(m{+}1)-4v_2(m{+}1),$$
and since $s_2(m+1)=s_2(m)+1-v_2(m+1)$ this equals $8m-1-4s_2(m)$; the increments
are strictly increasing, so
$$v_2\Bigl(\xi_2^{\rm Zud}-\frac{P_m}{Q_m}\Bigr)=8m-1-4s_2(m)\quad\text{exactly}.$$
Re-verified numerically at $m=5,10,17,32,33,64,100,127,128,200$ (`sep.gp`): exact
agreement in all ten.  **So the crystal route reproduces the exact tail of
`ZUDILIN_2ADIC.md` §5(iv) using only the Casoratian and the denominator law — the
Padé remainder estimate (5.3) is not needed for the rate.**  It is needed for the
*value*: nothing above identifies $\xi_2^{\rm Zud}$ with $\zeta_2(2)$.  That step
is Beukers' $\Theta_2(1/2)=8\zeta_2(2)$ and stays external.

### 5.4 Brown–Zudilin

$Q_n\in\mathbf Z$ and satisfies the Dwork tower ($p=2,3,5,7$; 1716 tests, 0
failures), so BZ is a Dwork row.  But $\mu_1=0$ and $\mu_{\rm det}=0$ at every
tested prime, so $\sigma_p=0$: **every prime is a good prime for Brown–Zudilin**,
the extension splits everywhere, and there is no $p$-adic Apéry limit at any $p$
— not because the row is bad, but because it is Apéry-perfect in the strongest
sense ($\prod\lambda_i$ a unit at every $p$).  The tower Frobenius is
$(1,p^{-5})$ on $(Q,P)$ and $(1,p^{-3})$ on $(Q,\widehat P)$, both with trivial
character, verified at $p=2,3,5,7$.  This is a **structural obstruction to the
two-row method for $\zeta(5)$ through BZ**, of exactly the kind
`THEORY_NOTES_06` §4(c) predicts: alignment primes are motivic invariants, and BZ
has none.  Changing the cellular integral will not create one unless it changes
the singular locus $p$-adically.

### 5.5 The tower trichotomy as a theorem for Laurent rows

> **Theorem (trichotomy, crystal form).** [(i) proved; (ii),(iii) verified]
> Let a row satisfy (F1).  Then $\lim_s\rho^A_s=1$ **(i, proved: Dwork's
> congruence, all $p$ including $2$)**, and
> $$\lim_s\rho^B_s=\begin{cases}
> \varepsilon(p)\,p^{-(w+1)}, & \sigma_p=0\ \text{(good)},\\
> 1, & \sigma_p>0,\ \xi_p\ne0,\\
> 0, & \sigma_p>0,\ \xi_p=0,\end{cases}$$
> with $w+1$ the Hodge depth and $\varepsilon$ the character of the row's $L$-value.
> The second and third lines are **proved** given the existence of $\xi_p$
> (Theorem F$'$(ii)); the first is verified on $144+20=164$ cells.

Census after this note: 12 modular rows $\times$ 6 primes $\times$ 2 bases (144,
`GOOD_PRIME_TOWERS.md`) $+$ Zudilin $\times$ 4 primes $\times$ 2 bases (8) $+$
Brown–Zudilin $\times$ 4 primes $\times$ 2 bases $\times$ 2 directions (12 usable)
$=$ **164 cells, no exception**.

---

## 6. Status ladder

**Proved, and now prime-uniform.**
1. $a_n\in\mathbf Z$, $\mu_1=0$, and $A_{mp^s}\equiv A_{mp^{s-1}}\pmod{p^s}$ for
   every constant-term row and every $p$ (Mellit–Vlasenko; Samol–van Straten).
   This closes `rem:remains`(1) of the draft *for Laurent rows* — the hypothesis
   $v_p(a_n)=O(\log n)$ is a theorem, not a measurement, and the
   Beukers–Coster-type input is no longer needed.
2. $\sigma_p=\mu_{\rm det}-2\mu_1$, and $b_n/a_n$ converges iff $\sigma_p>0$, with
   the exact rate (Prop. 3.1).  For (R2)/(R3), $\sigma_p=v_p(c)$.
3. The splitting criterion (§1.5): unipotent iff the two eigenvalues coincide;
   the obstruction $\widetilde\xi$ is then well defined.
4. The exact Zudilin tail (Prop. 5.2), from the Casoratian alone.

**Proved in the literature, imported.**
5. The MUM Frobenius matrix is triangular with diagonal $(1,p,\dots,p^{r-1})$ and
   $\Gamma_p$/odd-$\zeta_p$ off-diagonal entries (BV arXiv:2302.09603 Prop. 1.3,
   Thms 1.4–1.5) — for simplicial and hyperoctahedral families.
6. $p$-integrality of the canonical coordinate $q=t\exp(G/F)$ and the excellent
   lift (BV-III Cor. 7.11, Def. 7.14, Thm 7.15), for $p$ odd and $p\nmid\#\mathcal G$.

**Verified, not proved.**
7. $\varepsilon(p)p^{-(w+1)}$ at good primes: 164 cells.
8. Zudilin's row is a Dwork row (2288 congruence instances).
9. Brown–Zudilin is a Dwork row (1716 instances) and good at every prime.
10. $\mathcal E_p\mid P\iff p\mid c$ (36 + 72 tests).

**Open.**
11. **(the value)** Identify $\widetilde\xi$ with a $p$-adic $L$-value in general.
    This is exactly the syntomic-regulator item of `THEORY_NOTES_05` §5:
    Bannai–Kings for the modular case, and for a general Laurent row it would be a
    "$p$-adic Apéry limit = $p$-adic period of a mixed Tate/mixed elliptic motive"
    statement.  **No paper in the Beukers–Vlasenko series computes an extension
    period; they compute the pure-crystal $\alpha_j$.**  Note the shape mismatch
    worth understanding: their $\alpha_j$ involve only *odd* $\zeta_p$ (since
    $\zeta_p$ vanishes at even positive integers), while our $\xi_p$ are
    $\tfrac12L_p(w+1,\psi\omega^{-w})$ with $w+1$ even in half the census — so
    the extension periods are genuinely a different family of numbers from the
    pure-crystal ones, and even-weight $L_p$ values are exactly what the pure
    crystal cannot see.
12. **(the source–geometry bridge)** $\mathcal E_p\mid P\iff p\mid c$ (item 10).
13. **(the Laurent model for Zudilin)** find $\Lambda$ with $\mathrm{CT}(\Lambda^m)=16^mQ_m$;
    this would make Theorem F$'$ literally apply and make $\xi_2^{\rm Zud}$'s
    existence a corollary of published theory.
14. **($p=2$)** BV-I Prop. 3.3 and most of BV-III assume $p>2$; Mellit–Vlasenko
    does not.  Everything we need at $p=2$ (rows E, F, $\alpha$, $\varepsilon$,
    Zudilin) rests on the Mellit–Vlasenko half plus measurement.  A $p=2$ Cartier
    theory would remove the largest remaining caveat in the census.
15. **($\varepsilon=\gamma^{p-1}$)** §2.3.
16. **(higher rank)** for $r\ge3$ the extension can also be non-split against
    $y_1,\dots$; the $\log$-freeness argument of §1.5 rules this out for our
    rows, but not for a general extension of a MUM crystal.

---

## 7. Scripts and reproduction

```
lattice/crystal_thmF/sep.gp      # slope separation census, 84 cells; Zudilin tail  -> sep.log
lattice/crystal_thmF/dwork.gp    # Dwork tower congruences, all rows + Zudilin + BZ -> dwork.log
lattice/crystal_thmF/towers2.gp  # tower Frobenius eigenvalues, Zudilin + BZ        -> towers2.log
lattice/crystal_thmF/slopes.gp   # first (noisier) two-point rate estimate          -> slopes.log
```
All run with plain `gp -q` in a few minutes each; `lattice/euler_criterion/rows.gp`
supplies the row generators.  Nothing was run on **snake** for this note.

## 8. Literature, with the numbers used

* F. Beukers, M. Vlasenko, **Dwork crystals I**, arXiv:1903.11155, IMRN 2021 no. 12, 8807–8844.  Def. 3.5 (Hasse–Witt), Prop. 3.3 ($p>2$), Thm 4.3, Rmk 4.4, Def. 4.9, Thm 5.3, Thm 5.7, Cor. 5.9, Thm 6.1 (weight filtration), App. A (zeta functions).
* — **Dwork crystals II**, arXiv:1907.10390, IMRN 2021, 4427–4444.  Thm 2.3 ($p_v=\Lambda_\sigma\sigma(p_v)$), Cor. 3.1 ($\Lambda=q(t)/q(t^p)$), Thm 3.2 (= Mellit–Vlasenko), Rmk 4.5 (unit root as a limit), Thm 5.3 / Cor. 5.4 (A-hypergeometric).
* — **Dwork crystals III: from excellent Frobenius lifts towards supercongruences**, arXiv:2105.14841, IMRN 2023, 20433–20483.  Prop. 5.11–5.12, Thm 7.3, Prop. 7.7, Cor. 7.11, Lemma 7.13, Def. 7.14, Thm 7.15, Cor. 8.2.
* — **Frobenius structure and $p$-adic zeta values**, arXiv:2302.09603, Adv. Math. 480 (2025).  Def. 1.1, **Prop. 1.3** (the MUM Frobenius matrix), Thms 1.4–1.5.
* A. Mellit, M. Vlasenko, *Dwork's congruences for the constant terms of powers of a Laurent polynomial*, IJNT 12 (2016) 313–321, Thm 1.  K. Samol, D. van Straten, Ann. Math. Québec 39 (2015) 185–203.
* A. Straub, ANT 8 (2014) 1985–2008, Remark 1.4 (the Apéry $\Lambda$).
* F. Beukers, *Some congruences for the Apéry numbers*, JNT 21 (1985) 141–155 — proves $a_{mp^r-1}\equiv a_{mp^{r-1}-1}\bmod p^{3r}$, **not** the Dwork tower (see `NOME_INTEGRALITY.md` §8 for the citation correction).
* F. Beukers, *Irrationality of some $p$-adic $L$-values*, Acta Math. Sinica 24 (2008) 663–686 — Props 5.1, 6.1; the only route to the *value* $\zeta_2(2)$.
* P. Candelas, X. de la Ossa, D. van Straten, arXiv:2104.07816 §4.4 (the $\zeta_p(3)$ entry, conjectural).

**Correction of a premise in the brief.** arXiv:2004.04027 is not Dwork crystals II
(it is a paper of Chaika–Smillie–Weiss); the correct ids are as above.  There is
no Beukers–Vlasenko paper on "$p$-adic Apéry limits"; the closest is
arXiv:2302.09603, and it treats the pure crystal only.
