# Sol's directional proposal, evaluated quantitatively

**Verdict: the proposal fails, twice over, and the second failure is structural, not
quantitative.** (i) *Quantitatively*, the directional prime mass it can possibly supply is
$n^{-1}\log\Phi_n \approx 2.4$, whereas the coefficient box at the critical exponent has
$\beta_Z=E_2=14.3931\ldots$ — a shortfall of a factor $\approx 6$, and the shortfall is not
closable because $2.4$ is already the *absolute ceiling*
$n^{-1}\log\bigl(\gcd(V_n,U_n)_{\rm odd}/\gcd(X_n,Y_n)_{\rm odd}\bigr)$, measured with all
primes (including $\ell>6n$) for $n\le 60$. (ii) *Structurally*, the directional mass is
**entirely absorbed**: for every odd prime $\ell$ where the Nesterenko row has excess
valuation over the Zudilin row, the congruence $c_ZY_n+c_NU_n\equiv0\ (\mathrm{mod}\ S_n)$
**already forces** $\ell^{e_\ell}\mid c_Z$ for *every* vector of the lattice, kernel or not.
Net directional separation measured over $n=5,\dots,60$: $\mathbf{0}$ (exactly zero at every
$n$ tested except $n=10$, where it is $0.32/n$). This is Sol's own earlier absorption theorem
reappearing in the guise he thought had escaped it. The route gives no non-vanishing
information, and nothing here bears on the irrationality of $G$.

Scripts: `lattice/catalan_directional/{rows_common.gp, directional.gp, absorb.gp, detail.gp,
ceiling.gp}` (each run as `cat rows_common.gp X.gp > run.gp; gp -q run.gp`; total runtime
$<35$ s). Rows are rebuilt exactly, reusing `lattice/catalan_audit/{rows.gp,nest2.gp}`.

---

## 1. Normalisation, lattice, and the selection box

Fixed throughout, matching `NearCriticalAssembly.lean` / `NesterenkoCross.lean` and the
literal-window paper §§2–5:

$$X_n=D_{6n}^2\,2^{e_{3n}}Q_{3n},\quad Y_n=2^{e_{3n}}D_{6n}^2P_{3n},\qquad
V_n=4^{7n+1}D_{6n}^2B_n,\quad U_n=4^{7n}D_{6n}^2C_n,$$

with $e_m=\min(6m,\,4m+3+\lfloor\log_2(2m-1)\rfloor)$, $S_n=D_{6n}^2$ (Lean: `Sfac`),
$T_n=2^{\lfloor kn\rfloor}$ (Lean: `TexpNC`), $M_n=S_nT_n$, and
$\sigma(k)=n^{-1}\log M_n=12+k\log2$ (Lean: `sigmaNC`). The coefficient lattice is
$$\mathcal K_n=\{(c_Z,c_N)\in\mathbb Z^2:\ c_ZY_n+c_NU_n\equiv0\ (S_n),\quad
c_Za^Z_n+c_Na^N_n\equiv0\ (T_n)\},\qquad a^\bullet_n=X_n/S_n\ \text{resp.}\ V_n/S_n,$$
$q_n=(c_ZX_n+c_NV_n)/M_n$, $p_n=(c_ZY_n+c_NU_n)/M_n$. The second congruence is admissible
because $T_n\mid a^Z_nU_n-a^N_nY_n$; the Lean file proves $2^{24n}\mid a^Z_nY^N_n-a^N_ny^Z_n$
(`NesterenkoCross.lean:249`, `dvd_reduced_cross_N`), so $k<24$ is legitimate.

**Rates** (`NearCritical.lean`, paper (3.2)/(3.5)):
$A_1=27.5359435426\ldots$, $E_1=13.0995887908\ldots$, $A_2=29.3547738621\ldots$,
$E_2=14.3931452672\ldots$; $k_*=(E_1+E_2-12)/\log2=22.3512905953\ldots$;
$$F(k)=\tfrac12\bigl(E_1+E_2-12-k\log2\bigr)=\tfrac{\log2}{2}(k_*-k),\qquad
H(k)=F(k)+(A_2-E_2)=F(k)+14.9616285949\ldots$$

**The box.** Theorem 5.1's proof uses the anisotropic rectangle (paper (5.6))
$\mathcal C_n=\{|c_Z|\le e^{x_nn},\ |c_N|\le e^{(\kappa_n-x_n)n}\}$ with
$x_n=(\kappa_n+E_2-E_1)/2$ and $\kappa_n=n^{-1}\log[\mathbb Z^2:\mathcal K_n]$. In the
supercritical regime $\kappa_n\to\sigma(k)$ (the audit verified $[\mathbb Z^2:\mathcal K_n]$
really is $\asymp M_n$), so the **box exponents** are

$$\boxed{\ \beta_Z(k)=x=\tfrac12\bigl(12+k\log2+E_2-E_1\bigr),\qquad
\beta_N(k)=\sigma-x=\tfrac12\bigl(12+k\log2-E_2+E_1\bigr).\ }$$

| $k$ | $\sigma$ | $F(k)$ | $H(k)$ | $\beta_Z$ | $\beta_N$ |
|---|---|---|---|---|---|
| $k_*=22.35129$ | 27.49273 | $0$ | 14.96163 | **14.39315** ($=E_2$) | **13.09959** ($=E_1$) |
| 22.4 | 27.52650 | $-0.01688$ | 14.94475 | 14.41003 | 13.11647 |
| 23.0 | 27.94239 | $-0.22483$ | 14.73680 | 14.61797 | 13.32441 |
| 23.9 | 28.56622 | $-0.53674$ | 14.42489 | 14.92989 | 13.63633 |

Note $\beta_Z\to E_2$, $\beta_N\to E_1$ as $k\downarrow k_*$: at criticality the box is exactly
the "both source forms balanced" rectangle. **Taking $k$ only infinitesimally above $k_*$ does
not shrink the box** — Sol's phrase "any positive exponential directional mass is enough"
is simply wrong. $F$ shrinks to $0$, but $\beta_Z$ stays pinned at $E_2\approx14.39$. This is
the first and already fatal error in the proposal.

## 2. Valuation bookkeeping: when does $\ell\mid c_Z$ actually follow?

Under $G=a/b$, $(a,b)=1$, a kernel vector satisfies $c_Zr^Z_n+c_Nr^N_n=0$ with
$r^Z=aX_n-bY_n$, $r^N=aV_n-bU_n$. Hence for each prime $\ell$,
$$v_\ell(c_Z)\ \ge\ v_\ell(r^N)-v_\ell(r^Z)\ \ge\ \min(v_\ell V,v_\ell U)-v_\ell(r^Z).$$
The lower bound is *useless unless $v_\ell(r^Z)$ can be bounded above*, and $a,b$ are unknown.
The only unconditional upper bound is the ultrametric one: if $v_\ell(X_n)\ne v_\ell(Y_n)$ then
$v_\ell(r^Z)=\min(v_\ell X,v_\ell Y)$ exactly. If $v_\ell(X_n)=v_\ell(Y_n)=v$ then
$r^Z=\ell^v(aX'-bY')$ and $v_\ell(aX'-bY')$ can be arbitrarily large for suitable $a,b$ — no
bound. So the **rigorous** per-prime directional exponent is

$$e_\ell=\begin{cases}\max\bigl(0,\ \min(v_\ell V,v_\ell U)-\min(v_\ell X,v_\ell Y)\bigr),
& v_\ell(X_n)\ne v_\ell(Y_n),\\[2pt] 0,&\text{otherwise,}\end{cases}
\qquad \Phi_n=\prod_\ell \ell^{e_\ell}.$$

(Sol's "$\ell\mid V,U$ but $\ell\nmid r^Z$" never occurs literally, because $D_{6n}^2$ is
common to all four coordinates, so $\ell\mid r^Z$ for every window prime. The valuation form
above is the correct repair, and it is what his own hint "do not divide the common factors
out" amounts to.)

**Measured** for $n=5,\dots,60$ over all odd $\ell\le6n$ (`directional.gp`): the two columns
coincide — whenever there is excess, $v_\ell(X_n)\ne v_\ell(Y_n)$ holds — so the rigorous and
optimistic masses are equal, and
$$n^{-1}\log\Phi_n \in [1.62,\ 3.66],\quad \text{mean}\approx2.4,\ \text{no visible drift}.$$

**Absolute ceiling** (`ceiling.gp`, includes primes $>6n$):
$n^{-1}\log\bigl(\gcd(V_n,U_n)_{\rm odd}\bigr)-n^{-1}\log\bigl(\gcd(X_n,Y_n)_{\rm odd}\bigr)
= 2.16\ \text{–}\ 3.65$ for $n\le40$; primes $\ell>6n$ contribute $\le0.73/n$ sporadically and
$0$ typically. So $2.4$ is not an artefact of the scan range: **no** refinement of the
prime-window theory can raise $\Phi_n$ above $\approx e^{3n}$.

**The window $4n<\ell<6n$ contributes essentially nothing.** At $n=20$ the excess-valuation
primes are $\ell=3,5,7,11,29,41,43,47,53,61,67,71,73,79$ — scattered, with $e_\ell=1$, and the
window $(4n,6n)=(80,120)$ contributes $e_\ell=0$ at every prime. Sol's identification of
$(4n,6n)$ as the source of directional mass is not borne out.

## 3. The absorption check — this is what kills it

For each odd $\ell$, write $s=v_\ell(S_n)=2v_\ell(D_{6n})$, $y=v_\ell(Y_n)$, $u=v_\ell(U_n)$.
The congruence $c_ZY_n+c_NU_n\equiv0\ (\ell^{s})$ is solvable for $c_N$ iff
$v_\ell(c_ZY_n)\ge\min(u,s)$, i.e. **every** $c\in\mathcal K_n$ obeys
$$v_\ell(c_Z)\ \ge\ f_\ell:=\max\bigl(0,\ \min(u,s)-y\bigr).$$
Only $e_\ell-f_\ell$ is information about the *kernel* specifically. Measured
(`absorb.gp`, $n\le60$):

| $n$ | $n^{-1}\log\Phi_n$ | $f$-mass $/n$ | **net** $(e-f)$-mass$/n$ |
|---|---|---|---|
| 8 | 3.65211 | 3.65211 | **0** |
| 12 | 2.63299 | 2.63299 | **0** |
| 20 | 2.40256 | 2.40256 | **0** |
| 32 | 2.88001 | 2.88001 | **0** |
| 44 | 2.40908 | 2.40908 | **0** |
| 60 | 2.46653 | 2.46653 | **0** |

(only exception in the whole range: $n=10$, net $=0.32/n$, from a single prime with $u>s$.)

**Why, structurally.** The $n=20$ table (`detail.gp`) shows the mechanism: for every odd
$\ell\le6n$ one has $v_\ell(Y_n)=0$ (the factor $D_{6n}^2$ exactly clears the odd denominator
of $P_{3n}$; verified with at most two small-prime exceptions per $n$), while
$v_\ell(U_n)\le v_\ell(S_n)$ (verified for every $\ell$ and every $n\le60$ except one prime at
$n=10$). Then
$$e_\ell=\min(v_\ell V,v_\ell U)-0\le u=\min(u,s)-y=f_\ell .$$
So $e_\ell\le f_\ell$ identically: **the excess valuation of the Nesterenko row is exactly the
divisibility the modulus $S_n$ already imposes on $c_Z$ for all lattice vectors.** An
$\ell$ that divides $c_Z$ for every vector of $\mathcal K_n$ separates nothing.

The symmetric statement (Zudilin excess forcing $\ell\mid c_N$) is worse: the rigorous mass is
$0$ at almost every $n$ (max $0.245/n$ at $n=12$), and its net after subtracting
$f^N_\ell=\max(0,\min(y,s)-u)$ is $0$ at every $n$ tested.

**Variant.** With the *literal-window* modulus $S_n^{\rm w}$ of paper (3.6) instead of
$D_{6n}^2$, absorption is incomplete: net mass $\approx1.6$–$2.2$ per $n$. But
$n^{-1}\log S^{\rm w}_n\approx5.6$, giving $\sigma\approx5.6+k\log2$ — deeply *subcritical*
($F>0$), the regime where the papers already conclude $\delta<1$ and where non-vanishing is
not the obstruction. To reach $\sigma\approx27.5$ one must reinstate $D_{6n}^2$, and the
absorption returns with it. There is no normalisation in which one has both the supercritical
$\sigma$ and un-absorbed directional mass.

## 4. Comparison, and what would remain even if the numbers had worked

$$n^{-1}\log\Phi_n^{\rm net}=0\quad\ll\quad n^{-1}\log\Phi_n^{\rm gross}\approx2.4
\quad\ll\quad \beta_Z(k_*)=14.3931.$$
Shortfall on the gross figure alone: $\approx11.99$ per $n$, a factor $e^{12n}$; on the net
figure the comparison is vacuous. Nothing marginal is happening.

Had the numbers come out favourable, the following would still all be open. They are listed
because the first two are load-bearing for *any* supercritical attempt, not just this one.

* **(a) Exact-valuation window laws for all $n$.** Everything in §§2–3 is a finite
  computation, $n\le60$. The statements "$v_\ell(Y_n)=0$", "$v_\ell(U_n)\le v_\ell(S_n)$",
  "$v_\ell(X_n)\ne v_\ell(Y_n)$ on the excess set" have no proof for general $n$; the paper's
  (3.7) only gives $S_n\mid X_n,V_n$.
* **(b) Supercritical 2-adic divisibility.** `dvd_reduced_cross_N`
  (`NesterenkoCross.lean:249`) *is* proved for all $n$, `sorry`-free, giving $2^{24n}\mid
  a^Z_nY^N_n-a^N_ny^Z_n$ — but **relative to the hypothesis bundle
  `NesterenkoPadicInputs`**, whose fields `val_four_BN`, `val_tail`/`val_Jform` (the exact
  2-adic slopes $-14n+3+s_2(n)+s_2(3n)$ and $28n-1-2s_2(n)-2s_2(3n)$) are *assumed*, not
  proved. So $k<24$ rests on two unproved exact-valuation laws. Note also that the audit's
  $T_n=2^{v_2(h_n)}$ with measured $v_2(h_n)/n\approx24.06\downarrow$ is a *different*,
  purely empirical modulus; only $\lfloor kn\rfloor$ with $k<24$ is covered.
* **(c) Minkowski selection at $F<0$.** `lattice_selection`
  (`GeometrySelection.lean:360`) needs no sign condition, but the conversion
  `worthiness_ge_of_balancedMinimaData` (`Final.lean:120`) and everything below it
  (`TwoRow.lean:380,498`, `Worthiness.lean:89`) require $0<F$. **Nothing in the repository
  covers $F\le0$**, and `worthiness_nearcritical` is explicitly gated on $k<k_*$. A
  supercritical selection theorem does not exist even in outline: at $F<0$ Minkowski's second
  theorem gives a short vector but no guarantee it is not the kernel vector — exactly the
  circularity the audit identified (`CATALAN_AUDIT.md` §4(c)).
* **(d) PNT inputs.** $\log D_m=m+o(m)$ and $\vartheta(6n)-\vartheta(4n)=2n+o(n)$ are used
  throughout; standard, but the $o(n)$'s must be uniform enough to survive the $F\to0^-$
  limit, where the whole margin *is* the $o(n)$.
* **(e)** Sol's own §"no amount of prime-factor information about $h_n$ alone can eliminate
  the supercritical zero branch" (his identity $\gamma_h=k_*\log2+H_*$) is correct and, note,
  applies verbatim here: $\Phi_n$ is a divisor of $\gcd(V_n,U_n)$, hence of the same
  arithmetic pool, and the required $e^{\beta_Zn}=e^{E_2n}$ exceeds anything that pool holds.

## 5. Verdict

The directional proposal is refuted on its own terms. Its central hope — that a prime dividing
both Nesterenko coordinates but not the Zudilin residual forces $\ell\mid c_Z$ — is correct as
a valuation statement but **carries zero information**, because the same primes divide $c_Z$
for every vector of the congruence lattice, kernel or not. Even ignoring that, the available
mass ($\approx2.4$, ceiling $\approx3.7$) is smaller than the box exponent
($\beta_Z=E_2\approx14.39$) by a factor of six, and $\beta_Z$ does **not** shrink as
$k\downarrow k_*$ — only $F$ does. The proposal's decisive premise, "we can take $k$ only
infinitesimally above $k_*$, so any positive directional mass suffices", is false.

No claim about the irrationality of Catalan's constant is made or supported here. The
conclusions of `CATALAN_AUDIT.md` stand unchanged: the proved results are $\delta=0.6589$
(LCM) and $\delta=0.9025$ (5:8), both $<1$, and the supercritical regime remains blocked by
non-vanishing, now with one more escape route closed.
