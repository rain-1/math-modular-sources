# Two-variable arithmetic holonomy bounds are multiplicative

*Fable, 2026-09-02. A closure note: what a holonomy bound for power series in two variables can and cannot do, at the level of the dimension count. The rigorous one-variable input is CDT's Theorem `basic main` (arXiv:2408.15403) with its Bost–Charles/slopes proof; the two-variable statements below are proved at the naive Cauchy–Liouville level and stated as expectations at the Bost–Charles level.*

## 1. Setting

$f_1,\dots,f_m\in\mathbf Q[[x,y]]$, linearly independent over $\mathbf Q(x,y)$, with coefficients $[x^ay^b]f_i\in\mathbf Z\cdot L_1(a)^{-1}L_2(b)^{-1}$ for denominator functions $L_1(a)=\prod_j[1..b_{1,j}a]$, $L_2(b)=\prod_j[1..b_{2,j}b]$ (rates $\sigma_1,\sigma_2$), and numerators of at most exponential growth. Maps $\varphi_1,\varphi_2\colon(\mathbf D,0)\to(\mathbf C,0)$ such that every $f_i(\varphi_1(z_1),\varphi_2(z_2))$ is holomorphic (or meromorphic) on the bidisc. Write $e_k:=\log|\varphi_k'(0)|-\sigma_k$ (naive entry; the rearranged $\tau$ replaces $\sigma$ in the refined version) and $M_k:=\max_{\mathbf T}\log^+|\varphi_k|$ (naive; the Bost–Charles integral $\mathrm{BC}_k$ in the refined version). Put $B_k:=M_k/e_k$, the naive one-variable bounds.

## 2. The naive two-variable bound

**Proposition 2.1 (naive form; proof at the Siegel–Liouville–Cauchy level, as in `CDT_UNPACKED.md` §1).** If $e_1,e_2>0$ then
$$m\ \le\ \frac{2\,M_1M_2}{e_1e_2}=2B_1B_2 .$$

*Argument.* Auxiliary function $F=\sum_iQ_i(x,y)f_i$, $Q_i\in\mathbf Z[x,y]$ of bidegree $<(D_1,D_2)$: $mD_1D_2$ unknowns. Impose $[x^ay^b]F=0$ on the weighted triangle $R_T=\{(a,b)\in\mathbf N^2:\ e_1a+e_2b<T\}$, which has $|R_T|\approx T^2/(2e_1e_2)$ points; take $|R_T|=(1-\varepsilon)mD_1D_2$. Since $F\ne0$, the set of $(a,b)$ with $\beta_{ab}:=[x^ay^b]F\ne0$ is nonempty and disjoint from $R_T$; let $(a_0,b_0)$ be an element of it that is *minimal for the product order* (such elements exist and form a finite antichain). Then all $\beta_{a'b'}$ with $a'\le a_0$, $b'\le b_0$, $(a',b')\ne(a_0,b_0)$ vanish, so the coefficient of $z_1^{a_0}z_2^{b_0}$ in $F(\varphi_1(z_1),\varphi_2(z_2))$ equals $\beta_{a_0b_0}\varphi_1'(0)^{a_0}\varphi_2'(0)^{b_0}$ exactly, and Cauchy's estimate on the torus gives
$$|\beta_{a_0b_0}|\,|\varphi_1'(0)|^{a_0}|\varphi_2'(0)|^{b_0}\le\max_{\mathbf T^2}|F\circ\varphi|\le e^{D_1M_1+D_2M_2+o(D)} .$$
Liouville: $\beta_{a_0b_0}$ is a nonzero rational with denominator dividing $L_1(a_0)L_2(b_0)\cdot(\text{height of }Q)$, so $|\beta_{a_0b_0}|\ge e^{-\sigma_1a_0-\sigma_2b_0-o(D)}$ (heights of the $Q_i$ are absorbed as in the one-variable heuristic). Hence $e_1a_0+e_2b_0\le D_1M_1+D_2M_2+o(D)$; but $(a_0,b_0)\notin R_T$ means $e_1a_0+e_2b_0\ge T$. So $T\le D_1M_1+D_2M_2$, i.e. $2(1-\varepsilon)mD_1D_2e_1e_2=T^2\le(D_1M_1+D_2M_2)^2$, and minimising over $D_2/D_1=M_1/M_2$ gives the claim. $\square$

The same argument in $d$ variables, with the simplex $\{\sum e_ka_k<T\}$ of volume $T^d/(d!\prod e_k)$, gives $m\le\frac{d^d}{d!}\prod_kB_k$.

**Remark 2.2 (Bost–Charles level).** In one variable the slopes method turns $M$ into the Bost–Charles integral and $\sigma$ into the rearranged $\tau$, with no loss of constants, because CDT's $d\to\infty$ "horizontal" replication uses uniformly distributed exponent vectors rather than a simplex. The two-variable version of that argument (a $d_1\times d_2$ array of replicated variables) is expected to give
$$m\ \le\ C\cdot\frac{\mathrm{BC}_1}{\log|\varphi_1'(0)|-\tau_1}\cdot\frac{\mathrm{BC}_2}{\log|\varphi_2'(0)|-\tau_2},\qquad C\in[1,2],$$
and the natural conjecture is $C=1$. Nothing below depends on the value of $C$ beyond $C\ge1$.

## 3. Inventories on product hosts are multiplicative

**Lemma 3.1.** Let $V_Z\subset\mathbf Q[[x]]$ and $V_N\subset\mathbf Q[[y]]$ be the admissible inventories of two one-variable hosts (pure functions plus the conditional functions produced by a hypothesis, e.g. $H_Z=bB_Z-aA_Z$ and $H_N=bB_N-aA_N$ under $G=a/b$). Then the admissible inventory on the product host with the product map $(\varphi_1,\varphi_2)$ is $V_Z\otimes V_N$, of dimension $m_Zm_N$: pure functions are tensor products of pure functions, and the conditional part is $H_Z\otimes V_N+V_Z\otimes H_N$ (both factors of a product must be holomorphic on the respective pulled-back discs).

**Corollary 3.2 (no gain).** Two variables give a contradiction iff $\frac{m_Z}{B_Z}\cdot\frac{m_N}{B_N}>C\ge1$. If both one-variable counts fail ($m_Z<B_Z$, $m_N<B_N$) the two-variable count fails; if one succeeds the one-variable proof already exists. In particular, for every Catalan host in the census (all with $m/B<1$) no product host can help, and this holds for the Hadamard-type functions $A_Z\odot B_N-A_N\odot B_Z$ as well: they are diagonal restrictions of elements of $V_Z\otimes V_N$, and `HADAMARD_HOST.md`'s $-11$ is the one-variable shadow of the same multiplicativity.

**Remark 3.3 (the $2$-adic bridge is one-dimensional).** The lattice construction's input $v_2(a^Z_mb^N_n-a^N_nb^Z_m)\ge24n$ holds along the ray $m=3n$ of the index lattice, not on a two-dimensional region; the $2$-adic radius of convergence of $W(x,y)=A_Z(x)B_N(y)-A_N(y)B_Z(x)$ in both variables is governed by the generic coefficients and equals that of the individual rows. So no two-variable $p$-adic overconvergence (and no multi-place gain of the `ADELIC_HOLONOMY.md` type) is available either.

## 4. What a gainful two-variable system would need

By Proposition 2.1/Remark 2.2 the only way to profit from two variables is a *super-multiplicative* inventory: a two-variable local system on $\mathbf P^1\times\mathbf P^1$ (or another surface) whose pure functions and conditional orbit together exceed the tensor product of what its restrictions to the two axes provide. For a tensor-product system (ranks $r_1r_2$) the inventory is exactly multiplicative. For the natural non-tensor systems the holonomic rank is *sub*-multiplicative (Appell $F_1$: rank $3<4$; $F_4$ is a twisted product by Bailey's formula), so their conditional orbits are smaller, not larger. A gainful system would need a hypothesis producing several independent conditional functions at once — e.g. a period that is constant along a whole fold *curve* while the pure inventory stays small — and no such system is known to us. This is the precise form of the obstruction; it is not a matter of constants.

## 5. Verdict

A two-variable arithmetic holonomy bound is provable (naively above; rigorously by transplanting the slopes method) and is multiplicative in the two one-variable bounds. It cannot help any problem whose one-variable counts fail on both factors, which includes every Catalan host in this repository, and it cannot see the $2$-adic alignment behind the two-row lattice construction. The genuinely open question it leaves is structural: find a non-tensor two-variable local system with a super-multiplicative conditional inventory.
