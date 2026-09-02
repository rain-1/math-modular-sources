# Arithmetic of the companion sequences from the companion formula

*Fable, 2026-09-02. Theory in §§1–3 and §5, census in §4. The census was run by an agent in exact rational arithmetic (report `REPORT.md` and 28 GP scripts in `lattice/companion_arithmetic/`; all fifteen sporadic rows regenerated from their recurrences, the repository used only as a cross-check). Companion to `ASYMPTOTIC_CONSTANTS.md` (the analytic side of the same fold lemma).*

**Summary.** The companion formula (1.1) makes the companion $b_n$ of an Apéry-like row an explicit finite object. Over the fifteen sporadic rows (Zagier's six, the Almkvist–Zudilin six, Cooper's three): (i) the formula holds exactly; (ii) the refined denominator $R_n$ of Theorem 1.1 is exact at every bad prime but only an upper bound at good primes, where the exact mechanism is $v_p(\mathrm{den}\,b_n)=\max(0,r-v_p(c(p)e_{n,p}))$ for $n/2<p\le n$; (iii) the valuation transfer $v_p(b_n)-v_p(a_n)=v_p(\xi_p)$ holds without exception wherever $\xi_p\ne0$, and for the two rows with $\xi_p=0$ the difference grows like $\sigma_pn$; (iv) two new exact digit laws, $v_2(a_n)=s_2(n)$ for Cooper's $s_{10}$ and $s_{18}$, with the companion law $v_2(b_n)=s_2(n)-2\lfloor\log_2n\rfloor-1$ — Cooper's free integration made exact; (v) a **Lucas congruence for the companions** at every good prime,
$$\beta_{np+m}\equiv\psi(p)\,\beta_n\,a_m\pmod p,\qquad\beta_n:=p^{\,r\lfloor\log_pn\rfloor}b_n,$$
observed in all $199$ (row, prime) cells and **proved in §5** for the twelve Eisenstein rows from (1.1) and the Malik–Straub Lucas congruences. It determines the character $\psi$ of Cooper's three meromorphic sources ($\mathbf 1,\mathbf 1,\chi_{-3}$; the last previously conjectural).

## 1. The formula as an arithmetic statement

For a modular Apéry system $(\Gamma,t,F)$ with exactly rectified operator (Theorem inversion of `paper/sections/02_sources.tex`), the normalised companion is $B=F\cdot D^{-r}\Phi$ with source $\Phi=F\,\theta_qt=\sum_{m\ge1}c(m)q^m$, $c(m)\in\mathbf Z$, $D=q\,d/dq$, and $r=w+1$. Writing $e_{n,m}:=[t^n]\bigl(F(q)\,q^m\bigr)\in\mathbf Z$ (integrality of $F$ and of $q(t)\in t+t^2\mathbf Z[[t]]$, the latter automatic from $t(q)\in q+q^2\mathbf Z[[q]]$ by Lagrange inversion),
$$\boxed{\ b_n=\sum_{m=1}^{n}\frac{c(m)}{m^{\,r}}\,e_{n,m}\ }\qquad(e_{n,m}=0\text{ for }m>n,\ e_{n,n}=1).\tag{1.1}$$
This is a finite formula for the numerators of an Apéry-like system in terms of the divisor sums $c(m)$ of Table tab:sources and the integer matrix $e_{n,m}$, the change of basis between $\{Fq^m\}$ and $\{t^n\}$. Normalisation throughout: $b_0=0$, $b_1=1$ (the companion solution $B_0=0,B_1=1$ of the recurrence).

**Theorem 1.1 (uniform integrality) [proved].** For every row with an integral source, $d_n^{\,r}b_n\in\mathbf Z$, $d_n=\mathrm{lcm}(1,\dots,n)$; more precisely
$$\mathrm{den}(b_n)\ \Bigm|\ R_n:=\operatorname{lcm}_{m\le n}\frac{m^r}{\gcd(m^r,c(m))}.$$
*Proof.* Each term of (1.1) has denominator dividing $m^r/\gcd(m^r,c(m))$. $\square$

This recovers in one line the case-by-case integrality statements for Zagier's six rows ($r=2$) and the six Almkvist–Zudilin rows ($r=3$), and it locates Cooper's free integration: $d_n^{\,r-1}b_n\in\mathbf Z$ as soon as $m\mid c(m)$ for all $m$ (a *magnetic* source in the sense of Paşol–Zudilin), which never happens for a holomorphic Eisenstein source ($c(p)=\psi(p)+\chi(p)p^{r}\not\equiv0\bmod p$ when $\psi(p)\ne0$) and does happen for the meromorphic sources of $s_7,s_{10},s_{18}$. The bound $R_n$ is sharp at the bad primes and not at the good ones (§4.1).

**Remark 1.2 (Lambert form).** For an Eisenstein source with $c(m)=\sum_{d\mid m}\chi(d)\psi(m/d)d^{r}$ one has $c(m)/m^r=\sum_{e\mid m}\psi(e)\chi(m/e)e^{-r}$, hence
$$D^{-r}\Phi=\sum_{e\ge1}\frac{\psi(e)}{e^{r}}\,g_\chi(q^e),\qquad g_\chi(q)=\sum_{k\ge1}\chi(k)q^k\in\mathbf Q(q),$$
a rational function of $q$ twisted by $\psi(e)e^{-r}$: the companion is $F$ times a twisted Lambert series, and every denominator of $b_n$ is an $e^{-r}$ with $e\le n$. Note the asymmetry of the two characters: $\psi$ (attached to $e^{-r}$) is the one that survives modulo $p$ in §5; $\chi$ is invisible there.

## 2. Valuations at the bad primes

Let $p\mid c$ (the constant term of the recurrence, $c=\lambda_1\lambda_2$), with $p$-adic slope $\sigma_p=v_p(c)+2\kappa_p>0$ and $p$-adic Apéry limit $\xi_p$ (Proposition C and Theorem F of the paper: $b_n/a_n\to\xi_p$ with $v_p(b_n-\xi_pa_n)\ge\sigma_pn-O(\log n)$).

**Proposition 2.1 [proved].** If $\xi_p\ne0$ then for all large $n$ (as soon as $v_p(a_n)+v_p(\xi_p)<\sigma_pn-O(\log n)$, which holds since $v_p(a_n)=O(\log n)$)
$$v_p(b_n)=v_p(a_n)+v_p(\xi_p).$$
If $\xi_p=0$ then $v_p(b_n)-v_p(a_n)\ge\sigma_pn-O(\log n)$.
*Proof.* $b_n=\xi_pa_n+(b_n-\xi_pa_n)$ and the second term has strictly larger valuation. $\square$

So at a bad prime the companion's valuation law is the row's law shifted by the constant $v_p(\xi_p)$, a Kubota–Leopoldt valuation. Example: Zagier's Catalan row $\mathbf E$ at $p=2$: $v_2(a_n)=2s_2(n)$ (exact digit law) and $\xi_2=\tfrac12\zeta_2(2)$ with $v_2(\xi_2)=-2$ (measured exactly in §4.2, so $v_2(\zeta_2(2))=-1$) give $v_2(b_n)=2s_2(n)-2$, exact for $1\le n\le6000$; in the normalisation of `catalan/PTPAH_K_final_resolution_2026-08-27.md`, whose companion has $v_2(b_1)=1$, the same law reads $2s_2(n)-1$. The census (§4.2) confirms the transfer on every row and every bad prime, and records which rows admit an exact digit law for $a_n$ (the transferable input) as opposed to a bound.

## 3. What is new here and what is not

Lucas congruences and $p$-adic valuations of the *rows* $a_n$ are studied in the literature (Malik–Straub, Delaygue, Straub); the *companions* $b_n$ are usually treated only through $d_n^rb_n\in\mathbf Z$. Formula (1.1) makes the companion an explicit finite object: the denominator bound $R_n$ (Theorem 1.1), the valuation transfer (Prop. 2.1), the good-prime tower structure of `GOOD_PRIME_TOWERS.md` ($b_n/a_n$ converges $p$-adically along $n\mapsto np^k$ to $\Gamma_p$-assemblies), and now a Lucas law for $b_n$ itself (Theorem 5.1) are all consequences of the same three-line identity.

## 4. Census

Fifteen rows: $\mathbf A$–$\mathbf F$ ($r=2$), $\alpha,\gamma,\delta,\varepsilon(=\mathbf T),\zeta,\eta$ ($r=3$; note $\zeta=(9,3,-27)$, with $c=+27$ the sequence is not even integral), Cooper $s_7=(13,4,-27,3)$, $s_{10}=(6,2,-64,4)$, $s_{18}=(14,6,192,-12)$ in the $(a,b,c,d)$ normalisation of (R3). All entries are exact computations over the stated range.

### 4.1 Denominators

* Least $k$ with $d_n^{\,k}b_n\in\mathbf Z$ for $n\le3000$: $k=2$ for $\mathbf A$–$\mathbf F$, $k=3$ for the six third-order Eisenstein rows, and $k=2$ for all three Cooper rows (the free integration; sharp, since $d_nb_n\notin\mathbf Z$ at $n=2$).
* The companion formula (1.1) holds **exactly as rationals** for all twelve Eisenstein rows, $n\le60$, with $q(t)\in\mathbf Z[[t]]$ and $e_{n,m}\in\mathbf Z$ — an independent re-derivation of the twelve-family source table from the recurrences and the $c(m)$ alone.
* $\mathrm{den}(b_n)\mid R_n\mid d_n^{\,r}$ always, but $\mathrm{den}(b_n)=R_n$ for at most $306$ of the $n\le3000$ (as few as $3$ for $\varepsilon$); $R_n=d_n^{\,r}$ identically exactly for the Apéry-perfect rows $\mathbf D,\gamma$ ($c=\pm1$). The failure is sharply localised:
  * **bad primes are exact**: $v_p(\mathrm{den}\,b_n)=v_p(R_n)$ for every bad $p\mid c$, every row, all $n\le3000$, with the single exception of $\varepsilon$ at $p=2$, where $v_2(b_n)=-1$ exactly when $n=2^k$ and $R_n$ over-predicts one factor of $2$ otherwise (explained by the digit law of §4.3);
  * **good primes, exact mechanism**: for $n/2<p\le n$ only $m=p$ contributes a $p$ to (1.1), and $v_p(\mathrm{den}\,b_n)=\max\bigl(0,\ r-v_p(c(p)\,e_{n,p})\bigr)$ [all twelve rows, all $3\le n\le120$, all such $p$: $859$ cases per row, no exception]. So $R_n$ fails exactly when $p\mid e_{n,p}$ — a property of the nome $q(t)$, not of the source.

### 4.2 Bad primes: valuations, transfer, and $\xi_p$

* $v_p(a_n)\le\lambda_ps_p(n)$ with the paper's integer $\lambda_p$ reproduced exactly ($4$ for $\mathbf A,\mathbf F$ at $2$; $5$ for $\alpha$; $3$ for $\varepsilon,\delta,s_{18}@3$; $2$ for $\mathbf C,\mathbf F@3,\zeta,\eta$; $1$ for $\mathbf B,s_7,s_{10},s_{18}@2$), and the sharp rational $\lambda_p^*=\max_nv_p(a_n)/s_p(n)$ recorded ($\mathbf A$: $24/7$; $\mathbf C$: $13/7$; $\mathbf F$: $10/3$, $13/7$; $\alpha$: $9/2$; $\delta$: $12/5$; $\zeta$: $27/14$; $\eta$: $11/10$; $s_7$: $1/2$; $s_{18}@3$: $13/6$). No affine-in-$(n,s_p(n),\lfloor\log_pn\rfloor)$ law closes any of the non-exact cells (residual widths $6$–$16$).
* **Transfer (Prop. 2.1) confirmed**: $v_p(b_n)-v_p(a_n)$ is constant on $1\le n\le3000$ (at worst one exceptional small $n$) and equals $v_p(\xi_p)$ wherever $\xi_p\ne0$: $-1$ for $\mathbf B,\mathbf C,\mathbf F@3,\delta,\eta,s_{18}@3$; $-2$ for $\mathbf E,\alpha$; $-4$ for $\varepsilon$; $0$ for $\mathbf F@2$. Consistency across rows: $v_2(\zeta_2(3))=-2$ read off both $\alpha$ ($\xi_2=\frac13\zeta_2(3)$) and $\varepsilon$ ($\xi_2=\frac14\zeta_2(3)$); $v_3(\zeta_3(2))=-1$ read identically off $\mathbf B,\mathbf C,\mathbf F,s_{18}$, which even share the same eight leading $3$-adic digits.
* The two rows with $\xi_p=0$ ($\mathbf A$ at $2$, $\zeta$ at $3$ — exactly the rows with non-trivial co-divisor character) have $v_p(b_n/a_n)=3n+O(1)$, i.e. $\sigma_p=v_p(c)=3$, with $\xi_p=0$ to $p^{8954}$.
* Cooper's rows at their nominal bad primes have no limit: $v_p(b_n/a_n)=-2\lfloor\log_pn\rfloor-1$ exactly for $s_{10},s_{18}$ at $p=2$, and $-2\lfloor\log_3n\rfloor+[0,1]$ for $s_7$ at $p=3$ — the "slope $0$" of `ROW_LEDGER.md` made exact.
* All four Conjecture-D ratios ($\xi^{\mathbf B}=\xi^{\mathbf C}=\xi^{s_{18}}$, $4\xi^{\mathbf F}=5\xi^{\mathbf B}$, $3\xi^\alpha=4\xi^\varepsilon$) hold to the full Cauchy precision available at $N=3000$ ($3^{5991}$, $3^{2994}$, $3^{5991}$, $2^{11956}$).

### 4.3 Exact digit laws

A scan of all fifteen rows and all $p\le43$ for $v_p(a_n)=\lambda s_p(n)+\mu$ finds exactly four laws with $\lambda>0$ plus one near-law, all verified to $n\le10^4$:
$$v_3(a^{\mathbf B}_n)=s_3(n),\qquad v_2(a^{\mathbf E}_n)=2s_2(n),\qquad v_2(a^{s_{10}}_n)=v_2(a^{s_{18}}_n)=s_2(n),\qquad v_2(a^{\varepsilon}_n)=3s_2(n)-[n\text{ odd}].$$
The two Cooper laws are new, and the two Cooper valuation vectors are *identical*, for $a_n$ and for $b_n$ (matching the identical vectors noticed for $v_2(\mathrm{den}[q^m]\sqrt F)$ in `WEIGHT_DROP.md`). The implied companion laws, verified directly to $n\le6000$:
$$v_3(b^{\mathbf B}_n)=s_3(n)-1,\quad v_2(b^{\mathbf E}_n)=2s_2(n)-2,\quad v_2(b^{\varepsilon}_n)=3s_2(n)-[n\text{ odd}]-4\ (n\ge2),\quad v_2(b^{s_{10}}_n)=v_2(b^{s_{18}}_n)=s_2(n)-2\lfloor\log_2n\rfloor-1 .$$
Corollaries: $b^{\mathbf B}_n,b^{\mathbf E}_n$ are $p$-integral with $v_p=0$ exactly at the $p$-power indices; $b^{\varepsilon}_n$ has $v_2=-1$ iff $n=2^k$ (the one exception of §4.1); for $s_{10},s_{18}$, $v_2(d_n^2b_n)=s_2(n)-1\ge0$ with equality at $n=2^k$, so the free integration is exactly the coefficient $2$ (not $3$) of $\lfloor\log_2n\rfloor$.

### 4.4 Good primes: towers and the Lucas law

* **Towers.** For $\gamma$ and $\mathbf D$ (no bad primes) the literal congruence $b_{np^{k+1}}/a_{np^{k+1}}\equiv b_{np^k}/a_{np^k}$ fails by exactly $p^{r(k+1)}$; the normalised $\Lambda_k=p^{rk}b_{np^k}/a_{np^k}$ satisfies $\Lambda_{k+1}\equiv\Lambda_k\pmod{p^{e(k)}\Lambda_k}$ with $e(k)=3(k+1)$ in $18$ of $24$ cells (and $3(k+1)-|v_p(\Lambda_0)|$ in the six non-unit cells): the gain is the Beukers–Coster $3$ per level in both rows, although $r=3$ for $\gamma$ and $r=2$ for $\mathbf D$.
* **Lucas for $d_n^rb_n$: no.** $N_n:=d_n^rb_n$ satisfies neither $N_{np+m}\equiv N_nN_m$ (hit rates at chance) nor an Atkin–Swinnerton-Dyer three-term law. (Control: $a_{np+m}\equiv a_na_m\pmod p$ holds for all fifteen rows, $p\le13$, $np+m\le900$, no exception.)
* **Lucas for the companion: yes.** With $\beta_n:=p^{k_p\lfloor\log_pn\rfloor}b_n\in\mathbf Z_p$, $k_p=r$ for the Eisenstein rows and $k_p=2$ for Cooper's,
$$\beta_{np+m}\equiv\psi(p)\,\beta_n\,a_m\pmod p\qquad(n\ge1,\ 0\le m<p)$$
at every good prime, with $\psi$ the character of the first $L$-factor of $L(\Phi,s)$: verified as a search for $u\in\mathbf F_p$ over all $(n,m)$ with $np+m\le2000$ — a unique $u$ in every one of the $199$ cells ($15$ rows, $p\le43$), and $u=\psi(p)$ in every cell: $u\equiv1$ for $\mathbf A,\mathbf D,\alpha,\gamma,\delta,\varepsilon$; $u\equiv\chi_{-3}(p)$ for $\mathbf B,\mathbf C,\mathbf F,\zeta$; $u\equiv\chi_{-4}(p)$ for $\mathbf E$; $u\equiv\chi_5(p)$ for $\eta$. The iterated form $\beta_n\equiv\psi(p)^j\,b_{n_j}\prod_{i<j}a_{n_i}$ ($n=(n_j\cdots n_0)_p$) was checked directly for all $n\le2000$ in all cells. The congruence is modulo $p$ and no better (minimum defect valuation $1$ in $187$ cells, $2$ in twelve).
* **Cooper's characters.** The same search on Cooper's rows returns $\psi_{s_7}=\mathbf 1$, $\psi_{s_{10}}=\mathbf 1$, $\psi_{s_{18}}=\chi_{-3}$ (for $s_{18}$: twelve unramified determinations $5\le p\le43$, every one equal to $\chi_{-3}(p)$, plus $u=0$ at $p=3$), confirming by exact computation the $(\psi,\varphi)$ entries of `WEIGHT_DROP.md` §V2 and the conjectural $s_{18}$ line of `EULER_CRITERION.md` §4.1, without ever touching the meromorphic source. For $s_7,s_{10}$ the normalisation $k_p=2$ works at *every* prime including the nominally bad $3$, resp. $2$, with $u=1$: the exact form of "no $p$-adic slope anywhere" for these two rows. At the bad primes of the Eisenstein rows the unnormalised search gives $u=1$ in the unipotent cells and $u=0$ exactly for the two $\xi_p=0$ rows, reproducing the trichotomy of `obs:tower` modulo $p$.


### 4.5 Cooper's three sources: meromorphic, magnetic, and one congruence for everything

Second agent census (`lattice/cooper_sources/`, exact over $\mathbf Q$ to $m\le400$, Ligozat divisors, 60-digit numerics for the CM data). For the weight-two rows the canonical source is $\Phi=F\,Dx=x\sqrt{P(x)}\,F^2=u(1-Cu^2)F^2/(1+Bu+Cu^2)^2$ (Theorem 3.4 of `ASYMPTOTIC_CONSTANTS.md`), and for Cooper's rows:

* **Meromorphy.** $\Phi$ has a double pole (non-zero residue) at the points where $1+Bu+Cu^2=0$, which are CM points and never cusps: the two order-3 elliptic points $(5\pm\sqrt{-3})/14$ of $X_0(7)$ (disc $-3$) for $s_7$; the two order-2 elliptic points $(3+i)/10,(7+i)/10$ of $X_0(10)$ (disc $-4$) for $s_{10}$; the fixed point $(3+i)/6$ of the Atkin–Lehner involution $W_9$ (disc $-36$; $X_0(18)$ has no elliptic points) for $s_{18}$. Exact principal part $A_2=\nu^2/(4\pi^2g'(u_0))$ with $g'(u_0)^2=B^2-4C=c$, and $c(m)=-\nu^2\bigl(q_0^{-m}/g'(u_0)+\overline{\cdots}\bigr)\bigl(m-\tfrac1{2\pi\,\mathrm{Im}\,\tau_0}\bigr)+O(R_2^m)$, matching the exact coefficients to 60 digits at $m\approx400$; growth rates $e^{2\pi\,\mathrm{Im}\tau_0}=2.175682\ldots,\ 1.874456\ldots,\ 2.849653\ldots$. So $\Phi\notin M_4(\Gamma_0(N))$ for all three (also checked by `mftobasis`), and no Eisenstein-plus-cusp identification exists. What is exact is the factorisation above, with $F=\frac16(7E_2(7\tau)-E_2(\tau))$, $\frac1{12}(10E_2(10\tau)+5E_2(5\tau)-2E_2(2\tau)-E_2(\tau))$, $\frac14(18E_2(18\tau)-9E_2(9\tau)-12E_2(6\tau)+6E_2(3\tau)+2E_2(2\tau)-E_2(\tau))$ and $F^2$ = Eisenstein part $\bigl(\sum_{d\in S}d^2E_4(d\tau)\bigr)/\sum_{d\in S}d^2$ ($S=\{1,7\},\{1,2,5,10\},\{3,6\}$) plus $\frac{16}5f_{7.4.a.a}$, $\frac{28}{13}(f_{5.4.a.a}(\tau)+4f_{5.4.a.a}(2\tau))$, $12(\eta(3\tau)^8+4\eta(6\tau)^8)$ respectively — for $s_{18}$ the CM newform of level 9 for $\mathbf Q(\sqrt{-3})$, the visible origin of $\psi=\chi_{-3}$. The critical slot is still Eisenstein: $L(D^{-1}\Phi,2)=\zeta(2)/7,\ \zeta(2)/5,\ \tfrac12L(2,\chi_{-3})$ (60 digits).
* **The congruence.** For every prime $p\le43$, good and nominally bad alike, with $c'(m)=c(m)/m$,
$$c'(p^im)\equiv\psi(p)^i\,c'(m)\pmod{p^2}\quad(p\nmid m,\ i=1,2,3),\qquad\text{equivalently}\qquad \Phi\,|\,U_p\equiv\psi(p)\,p\,\Phi\pmod{p^3},$$
sharp at $p^2$ (resp. $p^3$) in 38 of the 42 cells ($p^3$ for $s_7$ at $p=2,5$; identically zero at $p=7$ for $s_7$ and $p=5$ for $s_{10}$, where $c(pm)=p\,c(m)$ exactly). Exact multiplicativity of $c'$ fails for almost all coprime pairs, and the next digit $(c'(p)-\psi(p))/p^2\bmod p$ is not a character: the mod-$p^2$ shape is exactly as much Eisenstein structure as the source carries.
* **One congruence, two consequences.** Dividing by $m$, the congruence is hypothesis (a) of Theorem 5.1 with $r-1=2$, hence the companion Lucas law $\beta_{np+m}\equiv\psi(p)\beta_na_m$ with $k_p=2$ — at bad primes too, which is what the first census observed. Iterating, $\Phi|U_{p^n}\equiv(\psi(p)p)^n\Phi\pmod{p^{n+2}}$ gives $p^n\mid m\Rightarrow p^n\mid c(m)$, i.e. Paşol–Zudilin's *strong $p$-magnetic property* (with two extra powers of $p$ and a character), hence $m\mid c(m)$, i.e. $D^{-1}\Phi\in\mathbf Z[[q]]$, i.e. Cooper's free integration $d_n^2b_n\in\mathbf Z$ (equivalently Bogner's $(n+1)\mid a_n$, re-verified to $n\le3000$). The free integration and the companion Lucas law are the same congruence. Cooper's sources are the level-$N$ members of the family of *magnetic modular forms* of Paşol–Zudilin (Nagoya Math. J. 248, 2022), whose level-one examples $\Delta/E_4^2$ and $E_4\Delta/E_6^2$ have double poles at the CM points of discriminant $-3$ and $-4$; their folklore conjecture that no holomorphic form is magnetic is the conceptual reason these sources have to be meromorphic. Proving the congruence (their route is a Shimura–Borcherds lift) would make the Cooper Lucas law unconditional.


## 5. The companion Lucas law is a consequence of the companion formula

**Theorem 5.1 [proved, given (a) and (b)].** Let the row have an integral modular parametrisation $t(q)\in q+q^2\mathbf Z[[q]]$, $F\in1+q\mathbf Z[[q]]$, with source coefficients $c(m)$ such that, at the good prime $p$,
$$c(p^im)\equiv\psi(p)^i\,c(m)\pmod{p^r}\qquad(p\nmid m,\ i\ge1)\tag{a}$$
(this congruence is all the proof uses; it holds with equality of multiplicativity for $\Phi=\sum_{d\mid N}\lambda_dE_{\psi,\chi,r}(d\tau)$ with $c_E(m)=\sum_{e\mid m}\chi(e)\psi(m/e)e^r$ and $p\nmid N$: multiplicativity, and $c_E(p^i)=\sum_{a\le i}\chi(p)^a\psi(p)^{i-a}p^{ra}$), and suppose the row satisfies the Dwork–Lucas congruence
$$F(t)\equiv F_{<p}(t)\,F(t^p)\pmod p,\qquad F_{<p}:=\sum_{m<p}a_mt^m,\tag{b}$$
which is equivalent to $a_{np+m}\equiv a_na_m\pmod p$ and is the Malik–Straub theorem for the sporadic rows. Then with $\beta_n=p^{r\lfloor\log_pn\rfloor}b_n$ and $n=(n_j\cdots n_0)_p$,
$$\beta_n\equiv\psi(p)^{\,j}\,b_{n_j}\prod_{i<j}a_{n_i}\pmod p,\qquad\text{hence}\qquad\beta_{np+m}\equiv\psi(p)\,\beta_na_m\pmod p .$$

*Proof.* Group the terms of $\Theta=D^{-r}\Phi=\sum_{m\ge1}c(m)m^{-r}q^m$ by the exact power of $p$ dividing the index: $\Theta=\sum_{i\ge0}p^{-ri}\Theta_i(q^{p^i})$ with $\Theta_i(q):=\sum_{p\nmid m}c(p^im)m^{-r}q^m\in\mathbf Z_p[[q]]$. By (a), $\Theta_i\equiv\psi(p)^i\Theta_0\pmod{p^r}$. In $b_n=[t^n](F\Theta)$ only $i\le j=\lfloor\log_pn\rfloor$ contribute (since $p^i\le p^im\le n$), so
$$\beta_n=p^{rj}b_n=\sum_{i\le j}p^{r(j-i)}[t^n]\bigl(F\,\Theta_i(q^{p^i})\bigr)\equiv[t^n]\bigl(F\,\Theta_j(q^{p^j})\bigr)\equiv\psi(p)^j\,[t^n]\bigl(F(t)\,\Theta_0(q(t)^{p^j})\bigr)\pmod{p^r},$$
because the $i<j$ terms carry the factor $p^{r(j-i)}$ in front of a $p$-integral series. Now reduce modulo $p$: $q(t)^{p^j}\equiv q(t^{p^j})$ (Frobenius on $\mathbf Z[[t]]$), and iterating (b), $F(t)\equiv\prod_{i<j}F_{<p}(t^{p^i})\cdot F(t^{p^j})$. Hence
$$\beta_n\equiv\psi(p)^j\,[t^n]\Bigl(\prod_{i<j}F_{<p}(t^{p^i})\cdot\bigl(F\Theta_0\bigr)(t^{p^j})\Bigr)\pmod p .$$
The first product has degree $<p^j$ in $t$ and the second factor is a series in $t^{p^j}$, so the coefficient of $t^n$ factors as $[t^{n\bmod p^j}]\prod_{i<j}F_{<p}(t^{p^i})\cdot[s^{n_j}](F\Theta_0)(s)=\prod_{i<j}a_{n_i}\cdot[s^{n_j}](F\Theta_0)$, and $[s^{n_j}](F\Theta_0)=b_{n_j}$ because every index $m\le n_j<p$ is prime to $p$. The two-step form follows by comparing the digit expansions of $n$ and $np+m$. $\square$

**Remarks.** (1) The intermediate statement $\beta_n\equiv\psi(p)^j[t^n](F\Theta_0(q^{p^j}))$ holds modulo $p^r$; only the two Frobenius reductions cost precision, which is why the law is modulo $p$ and no better in the census. A supercongruence would need Dwork's $F(t)/F(t^p)\equiv F_{<p^{s}}(t)/F_{<p^{s-1}}(t^p)\pmod{p^s}$ together with a lift of $q(t)^{p}\equiv q(t^p)$; the right modulo-$p^2$ correction is not identified. (2) The character that survives is $\psi$, the one attached to $e^{-r}$ in Remark 1.2 — the same $\psi$ whose value $\psi(p)p^{-w}$ governs the good-prime towers of `GOOD_PRIME_TOWERS.md`; Theorem 5.1 extends that statement from $n=ap^s$ to all $n$. (3) For Cooper's rows the census law has $k_p=2=r-1$ and $\psi=\mathbf 1,\mathbf 1,\chi_{-3}$; the proof goes through verbatim if $c(m)=m\,c'(m)$ with $c'$ satisfying (a) with $r-1=2$ in place of $r$. This is now **verified exactly** for every prime $p\le43$ (§4.5): $c'(p^im)\equiv\psi(p)^ic'(m)\pmod{p^2}$, sharp at $p^2$, at good and bad primes alike — so the companion Lucas law for Cooper's rows follows from the proof above. The stronger guess of the first version, that $c'$ is literally the coefficient system of a weight-$(r-1)$ Eisenstein series, is false: the sources are meromorphic (§4.5), $c'$ is not multiplicative and grows exponentially; it is Eisenstein-like modulo $p^2$ only.

## 6. Open

* No exact digit law for $v_p(a_n)$ in the ten remaining (row, bad prime) cells; whether a carry-counting law exists is open.
* The modulo-$p^2$ correction to Theorem 5.1, and the $e(k)=3(k+1)$ tower gain (fitted on $24$ cells, not derived).
* The good-prime denominator mechanism of §4.1 was checked to $n\le120$ (it needs the matrix $e_{n,m}$); extending it is computation, not ideas.

## 7. Reproduction

`lattice/companion_arithmetic/`: `lib.gp` (the fifteen rows, exact $(a_n,b_n)$), `src.gp` (the twelve $c(m)$), `t1*` (denominators, the companion formula, bad-prime exactness, the large-prime mechanism), `t2*` (valuations, $\xi_p$, Conjecture-D ratios), `t3*` (digit-law scan and verification to $n\le10^4$), `t4*` (towers, Lucas/ASD searches, the $\psi(p)$ law and its iterated digit form); full report in `REPORT.md`. `gp -q <file>`.
