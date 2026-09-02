# The holonomy bound over a number field: the denominator rate is charged once per place

*Fable, 2026-09-02. Settles the "single load-bearing unproved input" of `CDT_FINDER.md` §7 (C2): whether the Calegari–Dimitrov–Tang bound, applied to functions defined over $K=\mathbf Q(\sqrt5)$, gives the sum of the two real places' budgets against one denominator rate $\tau$ (entry $+4.18$ for the $X_1(5)$ symmetric-square host) or the average (entry $-0.44$ at the hard ceiling, $-0.90$ with CDT's contour). It is the average. The argument is the one-variable Siegel–Cauchy–Liouville count of `CDT_UNPACKED.md` §1 and `TWO_VARIABLE_HOLONOMY.md` §2, run over $K$; the refined (Bost–Charles) form has the same shape because that method is adelic. A by-product is a Schwarz-lemma statement (§4) that explains where the inventory of a holonomy count comes from.*

## 1. Setting

$K$ a number field, $[K:Q]=n_K$, archimedean places $v$ with local degrees $d_v=[K_v:\mathbf R]\in\{1,2\}$, $\sum_vd_v=n_K$, embeddings $\sigma_v$. Power series $f_1,\dots,f_m\in K[[x]]$, linearly independent over $K(x)$, with
$$L(n)\,a_{i,n}\in\mathcal O_K\quad\text{for all }i,n,\qquad L(n)=\prod_j[1,\dots,b_jn]\in\mathbf Z,\qquad \log L(n)\sim\sigma n,$$
(the denominators are *rational* integers — least common multiples — as in every modular Apéry system, where $d_n^{\,r}b_n$ is an algebraic integer by Beukers' argument) and with $|\sigma_v(a_{i,n})|\le C^n$. For each $v$ a holomorphic $\varphi_v\colon(\mathbf D,0)\to(\mathbf C,0)$ such that every $f_i^{(v)}\circ\varphi_v$ ($f^{(v)}:=\sigma_v$ applied to the coefficients) is holomorphic on $\mathbf D$; put $M_v:=\max_{|z|=1}\log^+|\varphi_v(z)|$. Write $\overline{X}:=\frac1{n_K}\sum_vd_vX_v$ for the average over places.

## 2. Theorem

**Theorem 2.1 (naive level; proved at the same level of idealisation as the $\mathbf Q$-case in `CDT_UNPACKED.md` §1).** If $\sum_vd_v\log|\varphi_v'(0)|>n_K\sigma$ then
$$m\ \le\ \frac{\sum_vd_vM_v}{\sum_vd_v\log|\varphi_v'(0)|-n_K\sigma}\ =\ \frac{\overline{M}}{\overline{\log|\varphi'(0)|}-\sigma}.$$

**Theorem 2.2 (refined form, expected).** Replacing $M_v$ by the Bost–Charles integral $\mathrm{BC}_v$ of $\varphi_v$ and $\sigma$ by CDT's rearranged rate $\tau$ (their (6.0.4–6)):
$$m\ \le\ \frac{\overline{\mathrm{BC}}}{\overline{\log|\varphi'(0)|}-\tau}.$$
This is the shape produced by the slopes method over $\mathrm{Spec}\,\mathcal O_K$: the arithmetic degree of the evaluation map is a sum over *all* places, the archimedean ones contributing $d_v\mathrm{BC}_v$ and $d_v\log|\varphi_v'(0)|$, and the finite places above a rational prime $p$ contributing $\sum_{v\mid p}d_v\log|L(n)|_v=n_K\log|L(n)|_p$ — the denominators are charged $n_K$ times. (CDT's adelic Theorem 2.5 of arXiv:2510.04156 has exactly this bookkeeping.)

*Proof of 2.1.* Auxiliary function $F=\sum_iQ_i(x)f_i(x)$, $Q_i\in\mathcal O_K[x]$ of degree $<D$: $mD$ unknowns in $\mathcal O_K$, i.e. $mDn_K$ rational integers. Impose $[x^n]F=0$ for $n<T$: $T$ equations over $K$, i.e. $Tn_K$ over $\mathbf Z$ after choosing an integral basis; clear denominators by $L(n)$. Take $T=(1-\varepsilon)mD$; Siegel's lemma gives $Q_i\ne0$ (heights treated as in the one-variable heuristic). Since $F\ne0$ (linear independence over $K(x)$), let $n_0\ge T$ be the order of $F$ at $0$ and $\beta=[x^{n_0}]F\in K^\times$; $L(n_0)\beta\cdot(\text{height factor})$ is a nonzero algebraic integer, so
$$\sum_vd_v\log|\sigma_v(\beta)|\ \ge\ -n_K\log L(n_0)-o(n_0)\ \ge\ -n_K\sigma n_0-o(n_0)\qquad\text{(Liouville over }K\text{: the norm of a rational integer is its }n_K\text{-th power).}$$
At each place, $F^{(v)}\circ\varphi_v$ is holomorphic on $\mathbf D$ with a zero of order $n_0$ at $0$ and leading coefficient $\sigma_v(\beta)\varphi_v'(0)^{n_0}$, so Cauchy on the unit circle gives $\log|\sigma_v(\beta)|+n_0\log|\varphi_v'(0)|\le DM_v+o(D)$. Multiply by $d_v$, sum over $v$, and insert Liouville:
$$n_0\Bigl(\sum_vd_v\log|\varphi_v'(0)|-n_K\sigma\Bigr)\ \le\ D\sum_vd_vM_v+o(D).$$
With $n_0\ge(1-\varepsilon)mD$ and $\varepsilon\to0$ the claim follows. $\square$

## 3. Consequences

**3.1 (Mediant).** For a system defined over $\mathbf Q$ one may still choose different maps at the places of any $K\supset\mathbf Q$; but $\frac{\sum_vd_vM_v}{\sum_vd_v(\log|\varphi_v'(0)|-\sigma)}\ge\min_v\frac{M_v}{\log|\varphi_v'(0)|-\sigma}$, so nothing is gained: the number-field bound is never better than the best single-place bound. For a system genuinely defined over $K$ the conjugate places are not optional — every $f_i^{(v)}\circ\varphi_v$ must be holomorphic — and a place with poor geometry is a tax on the average.

**3.2 (The $X_1(5)$ symmetric-square host, `CDT_FINDER.md` §6).** The two real places have ceilings $\log 2839=7.951$ and $\log 23.08=3.139$, average $5.545=\log256$ (the fold product is a unit, $|N(t_2)|=1$). Against CDT's proportional inventory ($\tau=5.980$, $k=3$) the entry is $-0.435$ at the ceiling and $-0.898$ with CDT's contour loss; only the "best inventory" ($u_j=m/2$, $\tau=4.837$) gives $+0.244$, and a positive margin then needs $m\gtrsim49$. The row "unnormalised (sum over places): $+4.183$" of that table is not a bound that any version of the theorem provides. So the load-bearing input is settled, and it lands on the side that was already known to fail with CDT's architecture: the target $1,\zeta(3),L(3,\chi_5)$ needs an inventory with about half of its functions free of all three LCM layers.

**3.3 (Where the extra Eichler integration is paid).** The averaged form makes the cost of weight transparent: $\tau$ grows by about $1.7$ per extra integration while the archimedean average over places of a $K$-rational host is at best the $\mathbf Q$-value of the best place; there is no mechanism by which the second place *adds* budget.

## 4. A Schwarz-lemma lemma: why the inventory is finite and what it must consist of

Let $x$ be the Hauptmodul-type coordinate ($x=q+O(q^2)$) of a modular Apéry system. Every modular form $f$ on the group, and every Eichler integral, is holomorphic on $\mathbf H$ and at the cusp $i\infty$ in the variable $q$.

**Lemma 4.1 [proved].** If $\varphi\colon(\mathbf D,0)\to(\mathbf C,0)$ factors as $\varphi=x\circ\psi$ with $\psi\colon\mathbf D\setminus\{0\}\to\mathbf H$ holomorphic and $q\circ\psi$ holomorphic at $0$ (i.e. $\varphi$ lifts to the $q$-disc), then $|\varphi'(0)|\le1$, with equality only for $\varphi(z)=x(e^{i\theta}z)$.
*Proof.* $q\circ\psi\colon\mathbf D\to\mathbf D$ fixes $0$, so $|(q\circ\psi)'(0)|\le1$ by Schwarz, and $\varphi'(0)=x'(0)(q\circ\psi)'(0)=(q\circ\psi)'(0)$. $\square$

**Consequences.** (i) For any such lift *every* modular form and Eichler integral, and all products and powers of them, are holomorphic on $\mathbf D$ — infinitely many $\mathbf Q(x)$-independent integral functions — consistent with the bound only because $\log|\varphi'(0)|\le0<\tau$. (ii) A disc with $|\varphi'(0)|>1$ (CDT's has $161$) therefore does **not** lift: its image contains, in its interior, a point over which $x$ is not a covering — a cusp image (on the first or a further sheet); folds do not help, since an orbifold point of order two lies in $\mathbf H$ and passing through it is covered by the lemma. (iii) At a cusp image every modular form of weight $k\ge1$ has a $\log^k$ singularity in $x$, and the local system of a weight-$k$ row has exactly one holomorphic solution there. Hence **the inventory of a holonomy count is confined to the combinations that are regular at the interior cusp(s)** — the period-purified functions — and modular forms and their powers are excluded, which is why $m$ is small and why "denominator-free" functions are scarce: a denominator-free function regular at an interior cusp must be a non-modular combination, and the only known mechanism producing such objects, the free integration of a magnetic source, removes one LCM layer, not all of them.

## 5. Verdict

The number-field holonomy bound is the average over places. The $+4.18$ scenario for $1,\zeta(3),L(3,\chi_5)$ was an accounting artefact; the honest entries are $-0.44$ (ceiling) to $-0.90$ (CDT's contour), positive only with an inventory about half of whose functions are free of every LCM layer and $m\approx50$. By §4 such functions would have to be non-modular combinations regular at the interior cusp at both real places, and no source of them beyond the one-layer free integration is known. The mixed weight-4 targets remain the right ambition; the obstruction is the inventory, exactly as for Catalan, and not the normalisation.
