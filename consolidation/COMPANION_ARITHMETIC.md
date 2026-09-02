# Arithmetic of the companion sequences from the companion formula

*Fable, 2026-09-02. Theory here; census in `scratchpad/companions/REPORT.md` (agent, pending), to be merged into §4. Companion to `ASYMPTOTIC_CONSTANTS.md` (the analytic side of the same fold lemma).*

## 1. The formula as an arithmetic statement

For a modular Apéry system $(\Gamma,t,F)$ with exactly rectified operator (Theorem inversion of `paper/sections/02_sources.tex`), the normalised companion is $B=F\cdot D^{-r}\Phi$ with source $\Phi=F\,\theta_qt=\sum_{m\ge1}c(m)q^m$, $c(m)\in\mathbf Z$, $D=q\,d/dq$, and $r=w+1$. Writing $e_{n,m}:=[t^n]\bigl(F(q)\,q^m\bigr)\in\mathbf Z$ (integrality of $F$ and of $q(t)\in t+t^2\mathbf Z[[t]]$),
$$\boxed{\ b_n=\sum_{m=1}^{n}\frac{c(m)}{m^{\,r}}\,e_{n,m}\ }\qquad(e_{n,m}=0\text{ for }m>n,\ e_{n,n}=1).\tag{1.1}$$
This is a finite formula for the numerators of an Apéry-like system in terms of the divisor sums $c(m)$ of Table tab:sources and the integer matrix $e_{n,m}$, the change of basis between $\{Fq^m\}$ and $\{t^n\}$.

**Theorem 1.1 (uniform integrality) [proved].** For every row with an integral source, $d_n^{\,r}b_n\in\mathbf Z$, $d_n=\mathrm{lcm}(1,\dots,n)$; more precisely
$$\mathrm{den}(b_n)\ \Bigm|\ R_n:=\operatorname{lcm}_{m\le n}\frac{m^r}{\gcd(m^r,c(m))}.$$
*Proof.* Each term of (1.1) has denominator dividing $m^r/\gcd(m^r,c(m))$. $\square$

This recovers in one line the case-by-case integrality statements for Zagier's six rows ($r=2$) and the six Almkvist–Zudilin rows ($r=3$), and it locates Cooper's free integration: $d_n^{\,r-1}b_n\in\mathbf Z$ as soon as $m\mid c(m)$ for all $m$ (a *magnetic* source in the sense of Paşol–Zudilin), which never happens for a holomorphic Eisenstein source ($c(p)=\psi(p)+\chi(p)p^{r}\not\equiv0\bmod p$ when $\psi(p)\ne0$) and does happen for the meromorphic sources of $s_7,s_{10},s_{18}$.

**Remark 1.2 (Lambert form).** For an Eisenstein source with $c(m)=\sum_{d\mid m}\chi(d)\psi(m/d)d^{r}$ one has $c(m)/m^r=\sum_{e\mid m}\psi(e)\chi(m/e)e^{-r}$, hence
$$D^{-r}\Phi=\sum_{e\ge1}\frac{\psi(e)}{e^{r}}\,g_\chi(q^e),\qquad g_\chi(q)=\sum_{k\ge1}\chi(k)q^k\in\mathbf Q(q),$$
a rational function of $q$ twisted by $\psi(e)e^{-r}$: the companion is $F$ times a twisted Lambert series, and every denominator of $b_n$ is an $e^{-r}$ with $e\le n$.

## 2. Valuations at the bad primes

Let $p\mid c$ (the constant term of the recurrence, $c=\lambda_1\lambda_2$), with $p$-adic slope $\sigma_p=v_p(c)+2\kappa_p>0$ and $p$-adic Apéry limit $\xi_p$ (Proposition C and Theorem F of the paper: $b_n/a_n\to\xi_p$ with $v_p(b_n-\xi_pa_n)\ge\sigma_pn-O(\log n)$).

**Proposition 2.1 [proved].** If $v_p(a_n)+v_p(\xi_p)<\sigma_pn-O(\log n)$ (which holds for all large $n$ since $v_p(a_n)=O(\log n)$), then
$$v_p(b_n)=v_p(a_n)+v_p(\xi_p).$$
*Proof.* $b_n=\xi_pa_n+(b_n-\xi_pa_n)$ and the second term has strictly larger valuation. $\square$

So at a bad prime the companion's valuation law is the row's law shifted by the constant $v_p(\xi_p)$, a Kubota–Leopoldt valuation. Example: Zagier's Catalan row $\mathbf E$ at $p=2$: $v_2(a_n)=2s_2(n)$ (exact digit law, `THEOREM_F_HYPOTHESES.md`) and $\xi_2=\tfrac12\zeta_2(2)$ with $v_2=-1$ give $v_2(b_n)=2s_2(n)-1$, the exact law found in `catalan/PTPAH_K_final_resolution_2026-08-27.md`. The census below tests this on every row and every bad prime, and records which rows admit an exact digit law for $a_n$ (the transferable input) as opposed to a bound.

## 3. What is new here and what is not

Lucas congruences and $p$-adic valuations of the *rows* $a_n$ are studied in the literature (Malik–Straub, Delaygue, Straub); the *companions* $b_n$ are usually treated only through $d_n^rb_n\in\mathbf Z$. Formula (1.1) makes the companion an explicit finite object: the exact denominator $R_n$ (Theorem 1.1), the valuation transfer (Prop. 2.1), and, at good primes, the tower structure of `GOOD_PRIME_TOWERS.md` ($b_n/a_n$ converges $p$-adically along $n\mapsto np^k$ to $\Gamma_p$-assemblies) are all consequences of the same three-line identity.

## 4. Census

*(to be filled from the agent report: exact denominators vs $R_n$, digit laws, $v_p(b_n)-v_p(a_n)$ vs $v_p(\xi_p)$, good-prime tower congruences for $\gamma$ and $\mathbf D$.)*
