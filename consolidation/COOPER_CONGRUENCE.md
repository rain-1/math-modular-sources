# Cooper's magnetic congruence: reduction, proved cells, the level-one theorem, and the trace identification

*Fable, 2026-09-02, version 2 (round 2 superseded several statements of version 1; the corrections are marked). Reports: `lattice/cooper_congruence/REPORT.md` (round 1) and `lattice/cooper_congruence/round2/REPORT.md` with `FINDINGS_PZ.md`, `FINDINGS_ROWS.md`, `FINDINGS_MODP.md`, `FINDINGS_XSIDE.md` (round 2). Status: the congruence of `paper/companions` Theorem "cooper" (eq:magnetic) is proved at two cells, is proved in full for every level-one magnetic form of Paşol–Zudilin, and for Cooper's rows is reduced to one finite-dimensional statement about vector-valued weight-5/2 forms.*

## 1. The objects

Cooper's rows $s_7,s_{10},s_{18}$: $u$ an eta quotient with $u|W_N=1/(Cu)$, $x=u/g(u)$, $g(u)=1+Bu+Cu^2$, $F=D\log u\in M_2(\Gamma_0(N))$, $A_n=[x^n]F$, source $\Phi=F\cdot Dx=x\sqrt{P(x)}F^2=\sum c(m)q^m$, meromorphic of weight 4 with a double pole at a CM point $\tau_0$ of discriminant $D_0=-3,-4,-36$. Magnetism $m\mid c(m)$ is proved (Bogner), so $\Xi:=D^{-1}\Phi=\sum c'(m)q^m\in\mathbf Z[[q]]$, $c'=c(m)/m$. Characters $\psi=\mathbf 1,\mathbf 1,\chi_{-3}$. Möbius inverse $\beta:=c'\star(\mu\psi)$, so $\Xi=\sum_e\beta(e)\Lambda_\psi(q^e)$.

## 2. Reduction and the two proved cells [proved]

**Theorem 2.1.** eq:magnetic ($\Phi|U_{p^n}\equiv(\psi(p)p)^n\Phi\bmod p^{n+2}$ for all $n$) $\iff$ $\Xi|U_p\equiv\psi(p)\Xi\pmod{p^2}$ $\iff$ $p^2\mid\beta(n)$ whenever $p\mid n$; over all $p$, $\iff\operatorname{rad}(n)^2\mid\beta(n)$.

**Theorem 2.2.** $\Phi_{s_7}|U_7=7\Phi_{s_7}$ and $\Phi_{s_{10}}|U_5=5\Phi_{s_{10}}$ identically ($\Phi|_4W_N=-\Phi$, Atkin–Lehner trace to level $N/p$, valence formula, two exact coefficients). These are the cells with $\varepsilon_p=-\psi(p)$, and they are the only ones where eq:magnetic is an identity.

## 3. The master conjecture and its Eichler-integral form

Verified for $n\le12000$ in all three rows, sharp except at $(s_7,3)$ and $(s_{10},2)$ where one extra power appears:
$$n^2\mid\beta(n).$$
With $\gamma(n)=\beta(n)/n^2$ and $K=\sum\gamma(d)q^d$ this is
$$\Xi=\sum_{e\ge1}\psi(e)(D^2K)(q^e),\qquad D^{-3}\Phi=\sum_e\psi(e)e^{-2}K(q^e)=\sum_n\gamma(n)\,\mathrm{Li}_2^{\psi}(q^n),\qquad L(\Phi,s)=L(s-1,\psi)\sum_d\gamma(d)d^{3-s}:$$
$\Xi$ has the shape of an inner weight-three Eisenstein series whose Lambert kernel $q/(1-q)$ is replaced by the integral kernel $K$, and $\Phi=D$ of it. Consequences: eq:magnetic at every prime; $c'(p^k)\equiv\psi(p)c'(p^{k-1})\pmod{p^{2k}}$; the $q$-product $\Xi=D\log\prod_n\bigl((1-q^n)^n\bigr)^{-\gamma(n)}$ for $\psi=\mathbf 1$.

**Correction to version 1.** Version 1 claimed the congruence fails for Paşol–Zudilin's level-one magnetic forms and concluded that it is special to Cooper's rows. That test used the wrong character. With $\psi=\chi_{D_0}$, $D_0$ the discriminant of the pole, and with their integralising constants ($64\Delta/E_4^2$, $108E_4\Delta/E_6^2$), $n^2\mid\beta_\psi(n)$ holds at every prime, and likewise for every entry of their Table 1 (verified $n\le500$, resp. $150$, sharp). The conjecture is a general statement about Shimura–Borcherds lifts; Cooper's rows are the level-$4N$ genus-character-twisted case, and the twist is exactly what makes $\psi$ trivial (resp. $\chi_{-3}$).

## 4. The level-one theorem [proved]

Let $f=c\,q^{-m_0}+O(q)\in S^{!,+}_{5/2}(\Gamma_0(4))$ with $-m_0$ a discriminant, $\Phi=\Psi(f)$ its Shimura–Borcherds lift, $\psi=\chi_{-m_0}$. Because $\dim S_4(\mathrm{SL}_2(\mathbf Z))=0$ a weakly holomorphic plus form is determined by its principal part, which forces the exact Hecke identities $f|T_{p^2}=\chi_{-m_0}(p)\,p\,f+p^3c\,g_{m_0p^2}$ and $g_{m_0p^{2r}}|T_{p^2}=g_{m_0p^{2r-2}}+p^3g_{m_0p^{2r+2}}$ on Borcherds–Zagier's integral basis $\{g_M\}$. A two-step induction gives, for $p\nmid m$,
$$a(p^{2a}m^2)+\bigl(1-\psi(p)\bigr)\sum_{j<a}p^{a-j}a(p^{2j}m^2)=p^{3a}\,c\,g_{m_0p^{2a}}(m^2),$$
uniformly for $\psi(p)=+1,-1,0$, and the left side is exactly what the $\psi$-twisted Möbius inversion assembles; hence $v_p(\beta_\psi(n))\ge2v_p(n)$ at every prime with $c\in\mathbf Z_p$, i.e. $n^2\mid\beta_\psi(n)$. This covers all of Paşol–Zudilin's Table 1 and is sharp. It shows *why* $\beta_\psi$ rather than $a(m^2)/m$ is the right object: at inert and ramified primes the correction terms are precisely the twisted Möbius terms.

## 5. The identification of the kernel [verified to 70–135 digits]

Let $f:=1/(xF)=g(u)/(uF)$. It is a weakly holomorphic form of weight $-2$ on $\Gamma_0(N)$, holomorphic on $\mathbf H$, vanishing simply at the polar CM points, $f|_{-2}W_N=-f$, with integral expansion ($f_{s_7}=q^{-1}+5+19q+52q^2+\cdots$), canonical (unique with pole order $\le1$ at both cusps) [proved]. Let $\widehat f:=Df+f/(2\pi y)=-R_{-2}f/4\pi$ be its Shimura–Maass derivative, a real-analytic $\Gamma_0(N)$-invariant function with $\widehat f|W_N=-\widehat f$ and algebraic CM values. Then, with the sum over $\Gamma_0(N)$-classes of Heegner forms $[A,B,C]$ of discriminant $D_0m^2$, $N\mid A$, $B\equiv\delta_0m\pmod{2N}$, twisted by the genus character:
$$\beta_{s_7}(m)=\sqrt{-3}\sum_Q\chi_{-3}(Q)\,\frac{\widehat f(\alpha_Q)}{\omega_Q}\qquad(7\nmid m),$$
and analogous laws for $s_{10}$ (twist $\chi_{-4}$, an Atkin–Lehner factor $2$ at even $m$) and $s_{18}$ (content-corrected character, odd $m$). The constant is the polar coefficient, $\lambda=4\pi^2A_2$. At level one the same recipe with $f=E_4E_6/\Delta$ reproduces Paşol–Zudilin's $f_{4a}$, $f_{4b}$ coefficient by coefficient at *all* admissible $d$, i.e. the weight-$5/2$ input is the generating function of the twisted CM traces of $\widehat f$ (the theta lift, checked). So

> **the master conjecture is an Ahlgren–Ono/Edixhoven statement: $n^2$ divides the genus-character-twisted Heegner trace at discriminant $D_0n^2$ of the first Shimura–Maass derivative of the weight-$(-2)$ form $1/(xF)$.**

Three routes are closed by this: $\beta(m)$ is not a trace of a modular *function* (it carries a $1/m$ tail with an irrational coefficient $1/(2\pi\operatorname{Im}\tau_0)$); $K$ is not modular and not even holonomic (its radius of convergence is transcendental by Gelfond–Schneider); and $\Phi^\flat=\sum_e\mu(e)\psi(e)e\,\Phi(e\tau)=D^3K$ is not modular.

## 6. The clean form, and the gap

For $s_7$ the weight-$5/2$ input is computed constructively: $c(d)=\sqrt{-3}\operatorname{Tr}_{-3d}(\widehat f)\in\mathbf Z$ for every admissible $d$ ($-3d$ a square mod $28$), with $c(m^2)=\beta_{s_7}(m)$, and
$$d\mid c(d)\quad\text{for every admissible }d\qquad\text{[verified: all 690 admissible }d\le2500\text{, no failure]}.$$
The master conjecture is the case $d=m^2$. The quotient $G(d)=c(d)/d$ extends $\gamma$ from the squares to all discriminants, is sharp (minimal valuation $0$ at $p=3,11,\ldots,29$, $1$ at $p=2,5$, exactly round 1's anomalous cells), and vanishes exactly when $7\mid d$ (the Atkin–Lehner cell off the squares). For $s_{10}$ it holds as well, with a single uniform constant $\lambda=2i=-\nu^2/g'(u_0)$ (420 admissible $d\le1500$, no failure, recomputed with minimal-$A$ Heegner representatives; the degenerate $d$ are exactly $25\mid d$, and $c(d)=0$ exactly when $5\mid d$ or $d\equiv5\pmod 8$; sharp at every $p\le59$), which dissolves the parity split of the earlier calibration. For $s_{18}$ the statement is vacuous: off the squares the twisted trace is $0$ for $d\equiv2\pmod3$ and irrational otherwise, so no integer $c(d)$ exists there; the clean form holds exactly for the rows whose CM discriminant is fundamental ($-3,-4$) and must be stated on the squares only for the conductor-$3$ row ($-36$).

Of the three ingredients of the level-one proof, two are in place at level $28$: uniqueness, because $\dim S_4(\Gamma_0(7))^{W_7=-1}=0$ (the eigenform has $a_7=-7$, Fricke sign $+1$) and $\Phi$ lives in the $W_7=-1$ part and vanishes at both cusps; and the eigenvalue, $c(p^2d)+p\,c(d)+p^3c(d/p^2)=\lambda c(d)$ with $\lambda\equiv p\pmod{p^2}$ constant in $d$ and a trivial middle character (verified $p\le19$), which in this normalisation is exactly what the target exponent $2a$ requires. **The single missing ingredient is an integral basis of the vector-valued space $M^!_{5/2}(\rho_L)$ for the Weil representation of the discriminant form $\mathbf Z/2N$, indexed by $(\beta,d)$, together with the tower relation written for $\rho_L$.** The support law of $c(d)$ is the index condition of that vector-valued form, not a scalar Kohnen condition (a scalar search in $M_{29/2}(\Gamma_0(28))/\Delta(4\tau)$ finds nothing, while the same pipeline recovers $f_{4a}$ at level 4).

## 7. Other results of round 2

* On the $x$-line the round-1 brick becomes $[x^{pn}]\bigl(l\,P^{(p-1)/2}/F_{<p}\bigr)\equiv\psi(p)L_n\pmod p$ with $F_{<p}$ the supersingular polynomial of $X_0(N)/\mathbf F_p$ in the coordinate $x$; the naive Kronecker route in weight $0$ is empty ($T_p=pU_p+V_p$).
* A new supercongruence $A_{p-1}\equiv\kappa\psi(p)p\pmod{p^2}$, $\kappa=-2,-3,-3$ (verified $p\le300$, proved for $s_{10}$): the recurrence alone manufactures $\chi_{-3}$.
* $\gamma\bmod2$ is multiplicative and $K\bmod2$ is a rational function of $q$; $\Xi\bmod2$ is Eisenstein ($c'(m)$ odd iff the $2P$-free part of $m$ is a square); outside four cells $K\bmod p$ has no structure and $\gamma(p)\bmod p$ is equidistributed.
* Magnetism is a genuine constraint: $E_4/H_D(j)^2$ is not magnetic for any tested $D$.
* `pz.txt` mangles Paşol–Zudilin's Lemma 1: $f_{4a}=\tfrac1{64}q^{-3}+q-506q^4+\cdots$, $f_{4b}=-\tfrac1{108}q^{-4}+q+\cdots$.

## 8. Status

| statement | status |
|---|---|
| eq:magnetic $\iff(S)\iff\operatorname{rad}(n)^2\mid\beta(n)$ | proved |
| eq:magnetic at $(s_7,7)$, $(s_{10},5)$ | **proved** |
| $n^2\mid\beta_\psi(n)$ for every level-one Shimura–Borcherds lift with single-index principal part (all of Paşol–Zudilin's Table 1) | **proved** |
| $\gamma$ = twisted CM traces of $\widehat f$, $f=1/(xF)$ | verified 70–135 digits, all three rows |
| $n^2\mid\beta(n)$ for Cooper's rows; $d\mid c(d)$ for $s_7$ | verified $n\le12000$, $d\le2500$ |
| level-$28$ proof: uniqueness and eigenvalue | proved / verified |
| level-$28$ proof: integral basis of $M^!_{5/2}(\rho_L)$ and the tower | **open, the whole gap** |
