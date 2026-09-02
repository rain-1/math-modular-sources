# Cooper's magnetic congruence: reduction, two proved cells, and the master conjecture

*Fable, 2026-09-02. Full report with 13 scripts: `lattice/cooper_congruence/REPORT.md`. Status: the congruence of `paper/companions` Theorem "cooper" (eq:magnetic) is proved at two of its cells and reduced everywhere else to one divisibility; the numerics show that divisibility is the shadow of a structural statement about the Eichler integral of the source, stated as Conjecture 4.1 below.*

## 1. The objects

Cooper's rows $s_7,s_{10},s_{18}$: $u$ an eta quotient with $u|W_N=1/(Cu)$, $x=u/(1+Bu+Cu^2)$, $F=D\log u\in M_2(\Gamma_0(N))$, $A_n=[x^n]F$, source $\Phi=F\cdot Dx=x\sqrt{P(x)}F^2=\sum_{m\ge1}c(m)q^m$, a meromorphic weight-4 form with a double pole at a CM point. Magnetism ($m\mid c(m)$) is proved (Bogner: $(n+1)\mid A_n$), so
$$\Xi:=D^{-1}\Phi=\sum_{m\ge1}c'(m)q^m\in\mathbf Z[[q]],\qquad c'(m)=c(m)/m,\qquad \Xi=l(x(q)),\ \ l(x)=\sum_n\frac{A_n}{n+1}x^{n+1}\in\mathbf Z[[x]].$$
Characters $\psi=\mathbf 1,\mathbf 1,\chi_{-3}$.

## 2. Reduction [proved]

**Theorem 2.1.** The following are equivalent:
(i) eq:magnetic: $\Phi|U_{p^n}\equiv(\psi(p)p)^n\Phi\pmod{p^{n+2}}$ for all $n\ge1$;
(ii) $(S)$: $\Xi|U_p\equiv\psi(p)\,\Xi\pmod{p^2}$, i.e. $c'(pm)\equiv\psi(p)c'(m)\pmod{p^2}$ for all $m$;
(iii) $p^2\mid\beta(n)$ for every $n$ divisible by $p$, where $\beta:=c'\star(\mu\psi)$ is the Möbius inverse, $\beta(n)=\sum_{d\mid n}\mu(d)\psi(d)c'(n/d)$, so that $\Xi=\sum_{e\ge1}\beta(e)\Lambda_\psi(q^e)$, $\Lambda_\psi(t)=\sum_d\psi(d)t^d$.
Over all $p$: (ii) $\iff\operatorname{rad}(n)^2\mid\beta(n)$ for all $n$.

*Proof.* (i)$\Rightarrow$(ii) is the case $n=1$ divided by $pm$; (ii)$\Rightarrow$(i) by induction, since $U_p$ preserves $p^2\mathbf Z[[q]]$ and $c(p^nm)=p^nm\,c'(p^nm)$. (ii)$\iff$(iii): $c'(pn)-\psi(p)c'(n)=\sum_{de=m}\psi(d)\beta(p^{v+1}e)$ for $n=p^vm$, $p\nmid m$, and the $\psi$-convolution is unitriangular. $\square$

**Theorem 2.2 (two cells).** $\Phi_{s_7}|U_7=7\Phi_{s_7}$ and $\Phi_{s_{10}}|U_5=5\Phi_{s_{10}}$ identically. *Proof sketch.* $\Phi|_4W_N=-\Phi$ for all three rows (from $u|W_N=1/(Cu)$, $F|W_N=-F$, $Dx|W_N=+Dx$); for $p\,\|\,N$ with Atkin–Lehner sign $\varepsilon_p=-1$ the trace to level $N/p$ is $\mathrm{Tr}\,\Phi=\Phi-\tfrac1p\Phi|U_p$. For $s_7$ this trace is a meromorphic weight-4 form on $\mathrm{SL}_2(\mathbf Z)$ with at most a double pole at $\rho$, vanishing at $\infty$; the valence formula bounds its order at $\infty$ by $1$, and the exact values $c(7)=7$, $c(14)=-98$ make it vanish to order $\ge3$, so it is zero. For $s_{10}$ the same on $\Gamma_0(2)$ with the single disc-$(-4)$ point and $c(5m)=5c(m)$ for $m\le3$. $\square$ These are exactly the cells where eq:magnetic was observed to be an identity; they are the cells with $\varepsilon_p=-\psi(p)$.

## 3. Where the character lives [proved equivalence, verified to $p\le101$]

By Lagrange inversion, $c(m)=[x^{m-1}]F(x)G(x)^m$ with $G=x/q(x)$, hence $c'(m)=\mathrm{Res}_x\bigl(\eta\,q(x)^{-m}\bigr)$ with
$$\eta=\frac{l(x)\,dx}{x\sqrt{P(x)}\,F(x)}=\Xi\,\frac{dq}{q}=\sum_ja_jx^jdx\in\mathbf Z[[x]]dx .$$
The Cartier operator is intrinsic, so $(S)\bmod p\iff\mathcal C(\eta)\equiv\psi(p)\eta\pmod p\iff a_{pj+p-1}\equiv\psi(p)a_j\pmod p$. The $j=0$ slot, $a_{p-1}\equiv\psi(p)\pmod p$, holds for all $p\le199$ and *manufactures* $\chi_{-3}$ for $s_{18}$ on the $x$-side, with no modular input. For $\psi=\mathbf 1$, Cartier's theorem turns this into: $\eta$ is logarithmic mod $p$, i.e. $\exp\int_0^x\eta\in\mathbf Z_p[[x]]$, a mirror-map-type integrality (verified to $x^{200}$). The full Lucas factorisation $a_{pj+r}\equiv a_ja_r$ fails; only the $r=p-1$ slot holds, so this is a unit-root statement, not a Dwork product.

**Negative control.** $(S)$ fails for Paşol–Zudilin's level-one magnetic forms $\Delta/E_4^2$, $E_4\Delta/E_6^2$, $E_6\Delta/E_4^3$ already at $n=2,3$. The congruence is therefore not a property of "weight 4 with a double pole at a CM point"; it is a property of Cooper's sources specifically.

## 4. The master conjecture

Numerically (all three rows, $n\le1500$, sharp) $n^2\mid\beta(n)$. Write $\beta(n)=n^2\gamma(n)$, $\gamma\in\mathbf Z$, $\gamma(1)=1$; e.g. $\gamma_{s_7}=1,-2,1,2,-5,6,0,-18,39,-30,-57,234,\ldots$, non-multiplicative, exponentially growing. Unwinding the Möbius inversion, this is:

**Conjecture 4.1 (Eichler-integral form).** With $K:=\sum_d\gamma(d)q^d\in q\mathbf Z[[q]]$,
$$\Xi=\sum_{e\ge1}\psi(e)\,(D^2K)(q^e)=\sum_{e,d}\psi(e)\,\gamma(d)\,d^2\,q^{ed},\qquad
\Phi=D^3\Bigl[\sum_{e\ge1}\psi(e)\,e^{-2}\,K(q^e)\Bigr],\qquad
L(\Phi,s)=L(s-1,\psi)\cdot\sum_d\gamma(d)d^{3-s}.$$
Equivalently: the Eichler integral $D^{-3}\Phi=\sum c(n)n^{-3}q^n$ (Bol, weight $-2$) equals $\sum_e\psi(e)e^{-2}K(q^e)$ with $K$ integral; or $c(n)=n\sum_{d\mid n}\psi(n/d)\,d^2\gamma(d)$.

This says that $\Xi$ has exactly the shape of an inner-orientation weight-three Eisenstein series, $E_3^{\psi,\mathbf 1}=\sum_e\psi(e)(D^2\Lambda_{\mathbf 1})(q^e)$, with the Lambert kernel $\Lambda_{\mathbf 1}=q/(1-q)$ (i.e. $\gamma\equiv1$) replaced by an integral series $K$; and that $\Phi=D$ of it. That is the corrected form of `COMPANION_ARITHMETIC.md` Remark 5.1(3): the source *is* $D$ of a weight-three Eisenstein-type series, generalised to a non-multiplicative kernel. Consequences of 4.1 (all verified in the stated ranges):
* $(S)$ and eq:magnetic at every prime (the $p\mid d$ terms contribute $0\bmod p^2$ and the others carry $\psi(p)$);
* the supercongruences $c'(p^2)\equiv\psi(p)c'(p)\pmod{p^4}$, $c'(p^k)\equiv\psi(p)c'(p^{k-1})\pmod{p^{2k}}$, $c'(pq)-\psi(p)c'(q)-\psi(q)c'(p)+\psi(pq)\equiv0\pmod{p^2q^2}$;
* the $q$-product $\Xi=D\log\prod_n(1-q^n)^{-n\gamma(n)}$ for $\psi=\mathbf 1$ (only $n\mid\beta(n)$ is needed for this), and its $\mathbf Z[\omega]$-analogue for $s_{18}$.

**Shimura–Borcherds reading.** If $\Phi$ is the lift of a weight-$5/2$ form $f=\sum a(n)q^n$ at level $4N$ with discriminant $D$, then $c(n)=\sum_{d\mid n}\psi(d)d\,a(|D|n^2/d^2)$ and Conjecture 4.1 reads $a(|D|m^2)=m^3\gamma(m)$: the square-index coefficients of the half-integral-weight input are divisible by $m^3$, three powers more than Paşol–Zudilin's level-one forms have, which is why those fail $(S)$. Divisibility of coefficients at $p^{2k}$-multiples of a discriminant by powers of $p$ is the Ahlgren–Ono/Edixhoven phenomenon for traces of singular moduli, proved there by the $p$-adic geometry of CM points; the weight-$5/2$ (rather than $3/2$) input would account for the third power. This is the route now being tried.

## 5. Status

| statement | status |
|---|---|
| eq:magnetic $\iff(S)\iff\operatorname{rad}(n)^2\mid\beta(n)$ | proved |
| eq:magnetic at $(s_7,7)$, $(s_{10},5)$, all $n$, all $m$, with equality | **proved** (Thm 2.2) |
| $(S)$ at all $p\le199$; $n^2\mid\beta(n)$ for $n\le1500$ | verified |
| $\mathcal C(\eta)\equiv\psi(p)\eta\pmod p$, i.e. $a_{pj+p-1}\equiv\psi(p)a_j$ | verified $p\le101$, **the smallest missing brick** |
| lift from mod $p$ to mod $p^2$ | needs one identity mod $p$ for the Frobenius lift $q(x)^p=q(x^p)(1+ph)$ |
| Conjecture 4.1 | open; identification of $K$ (traces? Shimura preimage?) is the lever |

Routes closed by the report: Hida/ordinary projection carries no information (the projection of $\Xi$ is $\Xi$ mod $p^2$) and only reaches inert primes; the level-$Np$ modular-polynomial route works only at $p\mid N$ (where it gives Theorem 2.2); "generic magnetism" is refuted by the negative control.
