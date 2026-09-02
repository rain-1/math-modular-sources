# Cooper's magnetic congruence, round 2: $\gamma$ identified as a twisted CM trace, the master conjecture shown to be a general property of Shimura–Borcherds lifts (and proved at level one), and its clean form $d\mid c(d)$

*Working note, 2026-09-02.  PARI/GP 2.15.4; all divisibility claims are exact integer
arithmetic, all CM evaluations are floating point with the precision quoted.  Scripts
`01_*`–`70_*` with matching `.log` files in this directory; four companion verdict files
`FINDINGS_PZ.md`, `FINDINGS_MODP.md`, `FINDINGS_XSIDE.md`, `FINDINGS_ROWS.md`; two working notes
`NOTES_IDENT.md`, `NOTES_TRACE.md`.  Background: `consolidation/COOPER_CONGRUENCE.md`,
`lattice/cooper_congruence/REPORT.md` (round 1), `lattice/cooper_sources/REPORT.md` §2,
`lattice/cooper_sources/pz.txt` (Paşol–Zudilin).  Nothing outside
`lattice/cooper_congruence/round2/` was created or modified.*

Claims are tagged **[proved]**, **[verified, range]** (exact computation over the stated
range), **[num, digits]**, **[refuted, range]**, **[conjectural]**.

---

## 0. Verdict table

| # | statement | status |
|---|---|---|
| **A0** | **THEOREM (level one).**  Let $f=c\,q^{-m_0}+O(q)\in S^{!,+}_{5/2}(\Gamma_0(4))$ with $-m_0$ a discriminant, $\Phi=\Psi(f)$, $\psi=\chi_{-m_0}$, $\beta_\psi=(A/n)\star(\mu\psi)$.  Then $n^2\mid\beta_\psi(n)$ for every $n$, at every prime $p$ with $c\in\mathbf Z_p$ — split, inert **and** ramified.  Proof: an exact $T_{p^2}$ eigen-identity plus a tower, giving $$a(p^{2a}m^2)+(1-\psi(p))\sum_{j<a}p^{a-j}a(p^{2j}m^2)=p^{3a}c\,g_{m_0p^{2a}}(m^2)$$ | **[proved]** `FINDINGS_PZ.md` §2.7, §5.3 |
| **A1** | **THE MASTER CONJECTURE IS NOT SPECIAL TO COOPER'S ROWS.**  Round 1's negative control (V6, `COOPER_CONGRUENCE.md` §3) used the *wrong character*.  With $\psi=\chi_{D_0}$, $D_0$ the discriminant of the pole, and with Paşol–Zudilin's own integralising constants, $n^2\mid\beta_\psi(n)$ holds for $64\,\Delta/E_4^2$ ($\psi=\chi_{-3}$) and $108\,E_4\Delta/E_6^2$ ($\psi=\chi_{-4}$) at **every** prime | **[verified, $n\le500$]** §5.1 |
| **A2** | The same holds across **all** of Paşol–Zudilin's Table 1 of strong magnetic weight-$4$ forms $\Psi(f_m)$: $m=7,8,11,19,43,67,163,15$ and the two $T_4$-twisted entries — magnetism $m\mid c(m)$ and $n^2\mid\beta_{\chi_{-m}}(n)$ in every case, and **sharp** ($\min_n\bigl(v_p(\beta(n))-2v_p(n)\bigr)=0$ for $p\ge7$) | **[verified, $n\le150$]** §5.2 |
| **A3** | So Conjecture 4.1 should be restated as a statement about **Shimura–Borcherds lifts**: if $\Phi=\Psi(f)$ with $f\in S^{!,+}_{5/2}$ having a single-index principal part $c\,q^{-m_0}$, then with $\psi=\chi_{-m_0}$ one has $n^2\mid\beta_\psi(n)$.  Route (D), closed by round 1 on the strength of the negative control, is **reopened** and is now the main route | **[conjectural]**, mechanism in `FINDINGS_PZ.md` §5.3 |
| **A4** | Magnetism is a genuine constraint: the naive family $E_4/H_D(j)^2$ (weight $4$, double poles exactly on the disc-$D$ Heegner divisor) is **not** magnetic for $D=-7,-8,-11,-19,-43,-67,-163,-15$ | **[refuted, $m\le260$]** §5.2 |
| **R1** | **THE IDENTIFICATION (row $s_7$).**  Let $f=1/(xF)=q^{-1}+5+19q+52q^2+\cdots$, the weakly holomorphic modular form of **weight $-2$** on $\Gamma_0(7)$ with a simple pole at the cusps, $f|_{-2}W_7=-f$; let $\widehat f:=Df+\dfrac{f}{2\pi y}=-\dfrac{R_{-2}f}{4\pi}$ ($R_{-2}$ = Maass raising), a real-analytic $\Gamma_0(7)$-invariant weight-$0$ function.  Then $$\beta_{s_7}(m)=\sqrt{-3}\sum_{Q}\chi_{-3}(Q)\,\frac{\widehat f(\alpha_Q)}{\omega_Q},$$ the sum over the $\Gamma_0(7)$-classes of Heegner forms $[A,B,C]$ of discriminant $-3m^2$ with $7\mid A$, $B\equiv5m\ (14)$ | **[verified, every $m\le45$ with $7\nmid m$, $\ge70$ significant digits]** §2 |
| **R2** | $f$ is weight $-2$, weakly holomorphic, holomorphic on $\mathbf H$ (vanishing to order $1$ at the polar CM points), $f|W_N=-f$; unique with pole order $\le1$ at both cusps.  The construction $f=1/(xF)$ is uniform across the three rows | **[proved]** §2.1 |
| **R3** | $\widehat f$ is $\Gamma_0(N)$-invariant of weight $0$; computed two independent ways (closed form via $E_2^*,E_4$ and $\mathrm{SL}_2(\mathbf Z)$-reduction; $q$-series with $\Gamma_0(N)+N$-reduction) agreeing to $50$ digits | **[proved]** (invariance) + **[num, 50]** §2.2 |
| **R4** | **The trace formula is general.**  With $f=E_4E_6/\Delta=q^{-1}-240-\cdots$ on $\mathrm{SL}_2(\mathbf Z)$ and the same $\widehat f$: $\beta_{F_{4a}}(m)=\mathrm{Tr}_{-3m^2}(\widehat f)/192$ and $\beta_{F_{4b}}(m)=-\mathrm{Tr}_{-4m^2}(\widehat f)/432$ (both with $\psi=\mathbf1$; the untwisted trace gives $\beta_{\chi_{D_0}}$ instead, so the genus character in the trace and the character $\psi$ are dual) | **[verified, $m\le12$, 50 digits]** §2.5 |
| **R1b** | The same law for $s_{10}$: $\beta(m)=i\,T(m)$ ($m$ odd), $2i\,T(m)$ ($m$ even), $\beta$-class $6m\ (20)$, twist $\chi_{-4}$; the factor $2$ is the Atkin–Lehner $W_2$, under which $\widehat f$ is invariant and which double-counts the classes at odd $m$ | **[verified, $m\le30$, 130 digits]** §2.4 |
| **R1c** | And for $s_{18}$: $\beta(m)=-\tfrac1{4\sqrt3}T^*(m)$ ($3\nmid m$), $-\tfrac1{2\sqrt3}T^*(m)$ ($3\mid m$), $m$ odd, with $\beta$-class $18m\ (36)$, the content-corrected character $\chi^*=\chi_{-3}\cdot(\tfrac{-4}{\operatorname{cont}})$, and the classes with $3\mid\operatorname{cont}$ excluded (**[proved]**: they have no Heegner representative).  Even $m$ degenerate | **[verified, odd $m\le29$, 135 digits]** §2.4 |
| **R4b** | **The theta lift in full.**  $a(d)=\dfrac{\sqrt d}{192}\operatorname{Tr}_{-3d}(\widehat f)$ for **every** $d\equiv0,1\ (4)$, checked against Paşol–Zudilin's explicitly constructed $f_{4a}=\tfrac1{64}q^{-3}+q-506q^4+\tfrac{131565}{64}q^5-\cdots$, including the non-square $d$.  So the weight-$5/2$ Shimura–Borcherds input **is** the generating function of the twisted CM traces of $\widehat f$ | **[verified, $d\le40$, 45 digits]** §2.7 |
| **R4c** | **The $s_7$ weight-$5/2$ input, computed.**  $c(d)=\sqrt{-3}\operatorname{Tr}_{-3d}(\widehat f)$ is an **integer** for every admissible $d$ (those with $-3d$ a square mod $28$: $d\equiv0,1\ (4)$ and $d\equiv0,1,2,4\ (7)$), with $c(m^2)=\beta_{s_7}(m)$; table in §2.8.  This gives the level-$28$ input constructively | **[verified, $d\le88$, 50 digits]** §2.8 |
| **R5** | **The constant is the polar coefficient.**  If $\Phi=A_2(\tau-\tau_0)^{-2}+\cdots$ then $\lambda=4\pi^2A_2$ up to sign: $4\pi^2A_2=1/192,\ -1/432$ for $F_{4a},F_{4b}$ (exactly the observed constants), and $|4\pi^2A_2|=\nu^2/|g'(u_0)|=\sqrt3,\,2,\,2/\sqrt3$ for $s_7,s_{10},s_{18}$, matching $|\sqrt{-3}|=\sqrt3$ | **[proved]** (level one) + **[verified]** §2.5 |
| **R6** | **The master conjecture is a divisibility of traces of singular moduli:** $n^2\mid\beta(n)\iff n^2\mid\lambda\,\mathrm{Tr}_{D_0n^2}(\widehat f)$ — the Ahlgren–Ono/Edixhoven setting, for the first Shimura–Maass derivative of a weight $-2$ form rather than for a modular function.  Route (B) is now precisely stated | **[proved, given R1]** §2.6 |
| **R7** | **Why the naive route (A)(1) fails.**  $\beta(m)$ is **not** a $\mathbf Q$-linear combination of twisted CM traces of modular *functions* at the discriminants $D_0m^2$: those are exponential polynomials with constant coefficients, whereas $\beta(m)$ carries a $1/m$ tail, $\beta(m)=\lambda_0(1-\kappa/m)R^m\sin(\cdot)+O(R^{\theta m})$, $\kappa=1/(2\pi\operatorname{Im}\tau_0)\notin\mathbf Q$.  Numerically $m\beta(m)/((m-\kappa)\tau_H(m))\to3$ to $11$ digits at $m=33$ | **[proved]** + **[verified, $m\le40$, all rows]** §3.1 |
| **R8** | Route (A)(3) is closed: $K=\sum\gamma(n)q^n$ is not a weakly holomorphic, holomorphic, quasi- or meromorphic modular form of any weight and level, and is **not holonomic** over $\mathbf Q(q)$ — its radius of convergence $e^{-\pi\sqrt{|D_0|}/(2N)}$ is transcendental (Gelfond–Schneider) while a holonomic series has algebraic singularities.  Same for $DK$, $D^2K$ | **[proved]** §3.2 |
| **R9** | $\Phi^\flat:=\sum_e\mu(e)\psi(e)e\,\Phi(e\tau)=D^3K$, whose coefficients are $n\beta(n)=a(n^2)$, is **not** a modular form (poles at $\tau_0/e$ for every squarefree $e$) | **[proved]**; $\Phi^\flat/\Phi$ not rational in $u$ of degree $\le6$ **[refuted]** §3.3 |
| **R10** | **Dictionary, calibrated.**  $\beta(m)=a(m^2)/m$ for the weight-$5/2$ Shimura–Borcherds input, so the master conjecture is $m^3\mid a(m^2)$.  Checked against the literature: for $F_{4a}$ this gives $a(1)=1$, $a(4)=-506$, matching $64f_{4a}=q^{-3}+q-506q^4+\cdots$ | **[proved]** + **[verified]** §4.1 |
| **R11** | Equivalent shapes: $\Xi=\sum_e\psi(e)(D^2K)(q^e)$; the **dilogarithm ladder** $D^{-3}\Phi=\sum_n\gamma(n)\mathrm{Li}_2^\psi(q^n)$; $\mathfrak g=\prod_n((1-q^n)^n)^{-\gamma(n)}$; $\Phi^\flat$ **triply magnetic**; $c(pn)\equiv\psi(p)p\,c(n)\pmod{p^{3(1+v_p(n))}}$; $m^3\mid a(m^2)$; $n^2\mid\mathrm{Tr}$ | **[proved]** (equivalences) §4.2 |
| **R12** | $n^2\mid\beta(n)$, all three rows | **[verified, $n\le12000$]** (round 1 had $1500$) |
| **M1–M3** | **A genuine mod-$p$ structure, only at $p=2$ (all rows) and $p=3$ ($s_{18}$).**  $\gamma(n)$ is odd $\iff$ $2^{1+v_2(N)}\nmid n$ and $P\nmid n$ ($P=7,5,3$); so $\gamma\bmod2$ is **multiplicative** and $K\bmod2$ is a **rational** function of $q$ (periods $14,20,12$); $K_{s_{18}}\equiv(q+2q^2+2q^3)/(1+q)^3\pmod3$ | **[verified, $n\le12000$]** `FINDINGS_MODP.md` |
| **M10** | Mod $2$ the whole conjecture degenerates to a weight-$3$ Eisenstein/theta identity: $c'(m)$ is odd $\iff$ the $2P$-free part of $m$ is a perfect square.  A concrete finite-looking target whose proof gives eq:magnetic at $p=2$ unconditionally | **[verified, $m\le6000$]** `FINDINGS_MODP.md` |
| **M4, L1–L6, N1** | Outside those four cells $K\bmod p$ is neither rational nor algebraic; all Lucas/Dwork laws for $\gamma$ fail massively; $\gamma(p)\bmod p$ matches nothing in a sweep of several hundred arithmetic candidates over $42$ primes (best $5/42$ = noise) and is equidistributed | **[refuted, stated ranges]** `FINDINGS_MODP.md` |
| **M9** | New unbounded divisibilities: $v_3(\gamma_{s_7}(n))\ge v_3(n)-1$ and $v_2(\gamma_{s_{10}}(n))\ge v_2(n)-1$, both sharp — $n^2\mid\beta(n)$ is **not** sharp at $(s_7,3)$, $(s_{10},2)$ | **[verified, $n\le12000$]** `FINDINGS_MODP.md` |
| **W1** | The task's $x$-side route (A) is **empty as posed**: in weight $0$, $T_p=p\,U_p+V_p$, so the Kronecker congruence mod $p$ reads $p\,U_pf\equiv0$ | **[proved]** `FINDINGS_XSIDE.md` |
| **W3–W7** | Its correct form is a real reduction: $f|U_p\equiv f_0+\sum_{n\ge1}[x^{pn}](fH)x^n\pmod p$ with $H=(\sqrt PF)^{p-1}\equiv P^{(p-1)/2}/F_{<p}$, where $F_{<p}=\sum_{n<p}A_nx^n$ is the **supersingular polynomial** of $X_0(N)/\mathbf F_p$ in the coordinate $x$ ($\deg F_{<p}=\lfloor\mu(N)(p-1)/(6\deg x)\rfloor$ exactly).  Round 1's smallest missing brick becomes $[x^{pn}]\bigl(l\,P^{(p-1)/2}/F_{<p}\bigr)\equiv\psi(p)L_n\pmod p$ | **[proved]** (equivalence) + **[verified, $p\le53$]** `FINDINGS_XSIDE.md` |
| **W11** | **New supercongruence** $A_{p-1}\equiv\kappa\,\psi(p)\,p\pmod{p^2}$, $\kappa=-2,-3,-3$; for $s_{18}$ it manufactures $\chi_{-3}$ from the Apéry-like recurrence alone | **[verified, $5\le p\le300$]**, **[proved for $s_{10}$]** `FINDINGS_XSIDE.md` |
| **W20** | The **$\psi$-twisted Dieudonné–Dwork tower** $a_{pn-1}\equiv\psi(p)a_{n-1}\pmod{p^{1+v_p(n)}}$ holds for all three rows, extending round 1's V19 to the character row | **[verified, $p\le31$]** `FINDINGS_XSIDE.md` |
| **P6** | **THE CLEAN FORM.**  With $c(d)$ the twisted CM trace coefficients of §2.8, $$\boxed{\ d\ \bigm|\ c(d)\ }\qquad\text{for every admissible }d,$$ verified for **all $690$ admissible $d\le2500$, zero failures**.  Conjecture 4.1 is the case $d=m^2$, since $c(m^2)=\beta_{s_7}(m)$ and $m^2\mid\beta(m)$ is $m^2\mid c(m^2)$; and $G(d):=c(d)/d$ extends $\gamma$ from the squares to all discriminants, $G(m^2)=\gamma_{s_7}(m)$ **[verified]**.  $G=1,-2,2,1,2,0,-5,0,5,6,6,4,4,15,0,14,20,-18,10,14,\dots$ at $d=1,4,8,9,16,21,25,28,29,32,36,37,44,53,\dots$  **Sharp**: $\min_{p\mid d}v_p(G(d))=0$ at $p=3,11,13,17,19,23,29$ and $=1$ at $p=2,5$ — exactly round 1's two anomalous cells $(s_7,2)$, $(s_7,5)$ | **[verified, $d\le2500$]** §5.5 |
| **P3** | Equivalently, on the $p$-parts: $p^{2a}\mid c(p^{2a}d)$ for every prime $p$, $a\ge1$, admissible $d$ — **295 tests, $0$ failures**, $p\le29$, $a\le4$, $p^{2a}d\le2500$, including the $118$ cells with $p\mid d$ | **[verified]** §5.5 |
| **P4** | **Uniqueness survives at level $28$.**  $\dim S_4(\Gamma_0(7))=1$, its eigenform has $a_7=-7$ hence Fricke eigenvalue $+1$, so $\dim S_4(\Gamma_0(7))^{W_7=-1}=0$ — and $\Phi|_4W_7=-\Phi$ lives exactly there.  This replaces $\dim S_4(\mathrm{SL}_2(\mathbf Z))=0$ in the level-one proof.  ($\dim M_4(\Gamma_0(7))=3$, so one also needs that $\Phi$ vanishes at both cusps, which it does) | **[proved]** §5.5 |
| **P5** | **The Hecke eigenvalue is right on the nose.**  Fitting $c(p^2d)+p\,c(d)+p^3c(d/p^2)=\lambda\,c(d)$ over all admissible $d$ gives $\lambda\equiv p\pmod{p^2}$, **constant in $d$**, with a **trivial** middle character ($\chi_{-3}(p)^2=1$), for $p=2,3,5,11,13,17,19$.  In the $c$-normalisation the target exponent is $2a$ (not $3a$), for which $\lambda\equiv p\ (p^2)$ is exactly what is needed | **[verified, $p\le19$]** `FINDINGS_PZ.md` §4.5 |
| **P2** | An **exact** $T_{p^2}$ eigen-identity behind Paşol–Zudilin's Theorem 1, with a tower: $f|T_{p^2}=\chi_{-m_0}(p)pf+p^3c\,g_{m_0p^2}$ and $g_{m_0p^{2r}}|T_{p^2}=g_{m_0p^{2r-2}}+p^3g_{m_0p^{2r+2}}$; and $v_p(a(p^{2j}m^2))\ge3j+v_p(a(m^2))$ $\iff$ $\lambda\equiv p\pmod{p^3}$ | **[proved / verified $p\le17$]** `FINDINGS_PZ.md` |
| **GAP** | $n^2\mid\beta(n)$: still open, but now in three interlocking precise forms — (i) $n^2$ divides the twisted CM trace $\mathrm{Tr}_{D_0n^2}(\widehat f)$; (ii) the $T_{p^2}$ eigenvalue of the weight-$5/2$ input is $\equiv p\pmod{p^3}$ together with the tower relation; (iii) the $x$-line congruence W6 | **open** §8 |

**Headline.**  Three things changed.  *First*, the kernel $K$ of Conjecture 4.1 is identified:
$\gamma(m)=\lambda m^{-2}\mathrm{Tr}_{D_0m^2}(\widehat f)$, where $\widehat f=-R_{-2}f/4\pi$
is the first Shimura–Maass derivative of the weight-$(-2)$ weakly holomorphic form
$f=1/(xF)$ and $\mathrm{Tr}$ is the genus-character-twisted Heegner trace; the master
conjecture becomes an Ahlgren–Ono/Edixhoven divisibility of traces of singular moduli.
*Second, and more consequential*: the congruence is **not** special to Cooper's rows, **and
at level one it is now a theorem**.
Round 1's negative control tested the wrong character; with $\psi=\chi_{D_0}$ the master
congruence holds for Paşol–Zudilin's level-one magnetic forms and for every entry of their
Table 1 that we could transcribe.  Conjecture 4.1 is therefore a general statement about
Shimura–Borcherds lifts of weight-$5/2$ forms; route (D) is reopened; and for
$f\in S^{!,+}_{5/2}(\Gamma_0(4))$ with a single-index principal part the statement is
**proved** (A0), covering every entry of Paşol–Zudilin's Table 1.  Cooper's rows are the
level-$4N$, genus-character-twisted analogue: the $\chi_{D_0}(Q)$ twist inside the trace is
exactly what turns the Hecke eigenvalue $\chi_{D_0}(p)p$ into $p$, hence $\psi=\mathbf1$
and divisibility at *every* prime.  Two of the three inputs that the level-one proof needs
are now in place at level $28$ as well — uniqueness, from
$\dim S_4(\Gamma_0(7))^{W_7=-1}=0$ (P4), and the eigenvalue $\lambda\equiv p\pmod{p^2}$
with a trivial middle character (P5) — so the single genuinely missing ingredient is an
integral basis of $M^!_{5/2}(\rho_L)$.  *Third*, the statement to aim at is cleaner than
Conjecture 4.1: with $c(d)$ the twisted CM trace coefficients,
$$d\mid c(d)\qquad\text{for every admissible }d$$
(all $690$ admissible $d\le2500$, zero failures), of which Conjecture 4.1 is the case
$d=m^2$; the quotient $G(d)=c(d)/d$ extends Cooper's kernel $\gamma$ from the squares to all
discriminants, and its sharpness fingerprint picks out exactly round 1's two anomalous cells
$(s_7,2)$ and $(s_7,5)$ (P6).

---

## 1. Setting, and what round 1 left

Cooper's rows $s_7,s_{10},s_{18}$: $N=7,10,18$; $u$ the eta quotient with $u|W_N=1/(Cu)$;
$g(u)=1+Bu+Cu^2$ with $(B,C)=(13,49),(6,25),(14,1)$; $x=u/g(u)$; $F=D\log u\in
M_2(\Gamma_0(N))$; $\Phi=F\,Dx=F^2u(1-Cu^2)/g(u)^2=\sum_{m\ge1}c(m)q^m$, a meromorphic
weight-$4$ form, holomorphic at all cusps, with a **double** pole at the CM points
$\tau_0$ where $g(u)=0$, of discriminants $D_0=-3,-4,-36$; $\psi=\mathbf1,\mathbf1,\chi_{-3}$;
$c'(m)=c(m)/m\in\mathbf Z$ (magnetism); $\beta=c'\star(\mu\psi)$; $\gamma(n)=\beta(n)/n^2$.

Round 1 proved eq:magnetic $\iff$ $(S)$ $\iff$ $\operatorname{rad}(n)^2\mid\beta(n)$, proved
two Atkin–Lehner cells outright, identified the mechanism as a Cartier-eigenvector statement
on the $x$-line, and left two things open: **(i)** what $\gamma$ *is*, and **(ii)** the
mod-$p$ brick $\mathcal C(\eta)\equiv\psi(p)\eta$.  This note answers (i) and reduces (ii).

Numerical baseline recomputed from scratch: $c',\beta,\gamma$ for $n\le12000$, all rows,
exact; $n^2\mid\beta(n)$ throughout (`01_data.gp`, `20_extend12000.gp`).

---

## 2. The identification of $\gamma$

### 2.1 The weight $-2$ form **[proved]**

**Lemma 2.1.**  $f:=\dfrac1{x\,F}=\dfrac{g(u)}{u\,F}$ is a weakly holomorphic modular form
of weight $-2$ on $\Gamma_0(N)$, holomorphic on $\mathbf H$, vanishing to order $1$ (in
$\tau$) at each polar CM point, with $f|_{-2}W_N=-f$ and integral $q$-expansion
$$f_{s_7}=q^{-1}+5+19q+52q^2+137q^3+316q^4+695q^5+1440q^6+2887q^7+5564q^8+\cdots$$

*Proof.*  $x|W_N=x$ and $F|_2W_N=-F$ (round 1, Lemma 6.1) give the weight and the sign.
$x$ has a pole of order $\nu=\operatorname{ord}_\tau(u-u_0)$ at $\tau_0$ ($\nu=3,4,4$), because
$g(u_0)=0$ simply.  $F=Du/u$ vanishes there to order $\nu-1$.  Hence $f=(1/x)(1/F)$ has
order $\nu-(\nu-1)=1$.  Away from $\tau_0$, $F$ does not vanish: the divisor of $F$ on
$X_0(N)$ has degree $2\mu(N)/12=4/3,\,3,\,6$, is supported off the cusps (constant term $1$
at $\infty$, and $F|W_N=-F$), and the polar CM points already account for all of it
($2\cdot\frac13\cdot 2$ for $s_7$; $3\cdot\frac12\cdot2$ for $s_{10}$; $3+3$ for $s_{18}$).
Zeros of $x$ occur only at the zeros of $u$, which are cusps. $\square$

For weight $-2$ with pole order $\le1$ at both cusps the space is one-dimensional
(any such form is $G/F$ with $G$ a modular function vanishing at $u=u_\pm$, i.e. divisible
by $g(u)$, and the pole condition forces $G=cg(u)/u$), so $f$ is canonical.

### 2.2 The Shimura–Maass derivative **[proved + num, 50 digits]**

Let $R_k=2i\,\partial_\tau+k\,y^{-1}$ be the Maass raising operator ($y=\operatorname{Im}\tau$),
so that $R_{-2}$ carries weight $-2$ to weight $0$, and set
$$\widehat f:=-\frac{R_{-2}f}{4\pi}=Df+\frac{f}{2\pi y},\qquad D=q\frac d{dq}.$$
$\widehat f$ is a **real-analytic $\Gamma_0(N)$-invariant function**, with
$\widehat f|W_N=-\widehat f$.  Its CM values are algebraic (Shimura–Maass); e.g.
$\widehat f\bigl(\frac{-5+\sqrt{-3}}{14}\bigr)=-\sqrt{-3}$.

Two independent evaluations agree to $50$ digits (`maass.gp`, `maass2.gp`, `t3.gp`, `t5.gp`):

* closed form $\widehat f=-f^2\bigl(F\,Dx+x\,DF\bigr)+\dfrac f{2\pi y}$ with
  $Dx=Fu(1-Cu^2)/g(u)^2$ and $DF=\sum_d\lambda_d\,d\,(E_2(d\tau)^2-E_4(d\tau))/12$, the
  $E_2^*,E_4$ evaluated by reduction to the $\mathrm{SL}_2(\mathbf Z)$ fundamental domain;
* the $q$-expansions of $f$ and $Df$ after reduction in $\Gamma_0(N)+N$ (used at the zeros
  of $F$, i.e. exactly at the polar CM points, where the closed form is $0\cdot\infty$).

### 2.3 The trace identity **[verified, $m\le45$, $\ge70$ significant digits]**

For $d=D_0m^2$ let $\mathcal Q_m$ be the set of $\Gamma_0(N)$-classes of Heegner forms
$[A,B,C]$ with $N\mid A$, $B\equiv\delta_0m\pmod{2N}$ ($\delta_0^2\equiv D_0\bmod 4N$;
$\delta_0=5$ for $s_7$), $B^2-4AC=D_0m^2$; let $\chi_{D_0}$ be the genus character and
$\omega_Q=|\overline{\Gamma_0(N)}_{\alpha_Q}|$.  Put
$$\operatorname{Tr}_m(\widehat f)=\sum_{Q\in\mathcal Q_m}\chi_{D_0}(Q)\,
\frac{\widehat f(\alpha_Q)}{\omega_Q},\qquad \alpha_Q=\frac{-B+\sqrt{D_0}\,m}{2A}.$$

$$\boxed{\ \beta_{s_7}(m)=\sqrt{-3}\;\operatorname{Tr}_m(\widehat f)\ }\qquad(7\nmid m)$$

**[verified]** for every $m\le45$ with $7\nmid m$: the worst absolute error is
$4\cdot10^{-61}$ against values up to $2.3\cdot10^{15}$, and the number of $\Gamma_0(7)$
Heegner classes found equals the number of $\mathrm{SL}_2(\mathbf Z)$-classes of
discriminant $-3m^2$ in every single case (`13_ident7.gp/.log`, `17_ident7wide.gp/.log`).
Equivalently
$$\gamma_{s_7}(m)=\frac{\sqrt{-3}}{m^2}\operatorname{Tr}_m(\widehat f),\qquad
a(m^2)=m\,\sqrt{-3}\operatorname{Tr}_m(\widehat f).$$
For $7\mid m$ the two $\beta$-classes $\pm5m$ merge, the trace becomes real, and
$\beta(m)=0$ by round 1's Theorem 6.3.

### 2.4 The other two rows **[verified, 130–135 digits]**

`FINDINGS_ROWS.md` (ADDENDUM), scripts `58_`, `59_`.  The same law holds, with the
Atkin–Lehner group of the level intervening.

**$s_{10}$** ($N=10$, $D_0=-4$, $\beta$-class $6m\bmod20$ — *not* $3m$, since
$3^2\not\equiv-4\ (40)$ — genus character $\chi_{-4}$, $f=1/(xF)$).  With
$T(m)=\sum_Q\chi_{-4}(Q)\widehat f(\alpha_Q)/\omega_Q$,
$$\beta_{s_{10}}(m)=i\,T(m)\ (m\text{ odd}),\qquad \beta_{s_{10}}(m)=2i\,T(m)\ (m\text{ even}),$$
**[verified to 130 digits, every $m\le30$ with $5\nmid m$]**.  The factor $2$ is the extra
Atkin–Lehner involution: $\Phi|W_2=+\Phi$ gives $f|W_2=+f$, so $\widehat f$ is
$W_2$-invariant; for **odd** $m$, $W_2$ preserves the $\beta$-class and both $\chi_{-4}$ and
$\widehat f$ termwise, acts as an involution on the classes with exactly one fixed class, so
$T$ double-counts every free orbit (the trace over a $W_2$-transversal gives the ratio
$2.0137,2.00225,2.00061,2.00028,1.999978,1.9999962$ at $m=7,9,11,13,17,19$).  For **even**
$m$, $\chi_{-4}$ vanishes exactly on the even-content forms, which are exactly the classes
whose $\Gamma_0(10)$-Heegner fibre has size $3$ rather than $1$ — so the GKZ bijection
genuinely fails there, the $W_2$-doubling does not happen, and $|\lambda|=2$ as §2.5
predicts.

**$s_{18}$** ($N=18$, $D_0=-36$, $\beta$-class $18m\bmod36$).  Three corrections were needed:
1. **[proved]** forms of discriminant $-36m^2$ with $3\mid\operatorname{cont}(Q)$ have **no**
   $\Gamma_0(18)$-Heegner representative at all ($Q=3Q'$ needs $6\mid Q'(p,r)$, and
   $3\mid p^2+r^2$ forces $3\mid p,3\mid r$ against $\gcd(p,r)=1$); restricting to
   $3\nmid\operatorname{cont}$ makes the class counts match for every $m\le30$;
2. a **bug in `heeg.gp`**: `heegrep` tries only one completion $(q,s)$ of a first column
   $(p,r)$, whereas the general matrix is $\binom{p\ \ q+tp}{r\ \ s+tr}$, shifting $B$ by
   $2At$ — so valid, in particular *smallest-$A$*, representatives were silently discarded
   (it returned $A=5634$ instead of $A=18$ for $[1,0,9]$, making $\operatorname{Im}\alpha$
   tiny and both evaluations of $\widehat f$ fail).  Fixed in `59_final3.gp` (`heegmin2`),
   which solves $B_0+2At\equiv\beta$ and always takes the minimal-$A$ representative.  *This
   affected none of the $s_7$ results*, where `#found = #classes` throughout and the
   closed-form evaluation of $\widehat f$ is insensitive to the height;
3. the genus character needs a **content factor**: with $\chi_{-3}$ alone the identity holds
   at $m=1,3,5,9,13,15,17,25$ and fails at $m=7,11,19,21,23$ — exactly the $m$ with a prime
   factor $\equiv3\ (4)$ other than $3$, i.e. inert in $\mathbf Q(i)$ — by exactly
   $\pm8\sqrt3$.  Since $-36$ is not fundamental, the character that works is
   $$\chi^*(Q)=\chi_{-3}(Q)\cdot\Bigl(\tfrac{-4}{\operatorname{cont}Q}\Bigr),$$
   trivial on primitive forms, which is why it is invisible for small $m$.

With those, and $T^*=\sum_Q\chi^*(Q)\widehat f(\alpha_Q)/\omega_Q$,
$$\beta_{s_{18}}(m)=-\tfrac1{4\sqrt3}T^*(m)\ \ (m\text{ odd},\ 3\nmid m),\qquad
\beta_{s_{18}}(m)=-\tfrac1{2\sqrt3}T^*(m)\ \ (m\text{ odd},\ 3\mid m),$$
**[verified to 135 digits, $m=1,3,5,\dots,29$]**.  For **even** $m$, $18m\equiv0\ (36)$, the
$\pm\beta$-classes merge and the family degenerates — the $s_{18}$ analogue of $7\mid m$ for
$s_7$ and $5\mid m$ for $s_{10}$; that half remains open.

### 2.5 The formula is general; the constant is the polar coefficient **[verified]**

The same construction at level one (`19_pzcontrol2.gp/.log`): $f=E_4E_6/\Delta=
q^{-1}-240-141444q-\cdots$ is *the* weight-$(-2)$ weakly holomorphic form on
$\mathrm{SL}_2(\mathbf Z)$ with pole order $1$, and with the same $\widehat f$ and the plain
twisted traces,
$$\beta_{F_{4a}}(m)=\frac{\operatorname{Tr}_{-3m^2}(\widehat f)}{192},\qquad
\beta_{F_{4b}}(m)=-\frac{\operatorname{Tr}_{-4m^2}(\widehat f)}{432}
\qquad\textbf{[verified, }m\le12\textbf{, 50 digits]}$$
for Paşol–Zudilin's $F_{4a}=\Delta/E_4^2$ (double pole at $\rho$) and $F_{4b}=E_4\Delta/E_6^2$
(at $i$).  A third discriminant confirms it: for $\Psi(f_7)$ of Table 1 (pole discriminant
$-7$) one gets $\beta_{\mathbf1}(m)=\tfrac17\operatorname{Tr}^{\chi_{-7}}_{-7m^2}(\widehat f)$
**and** $\beta_{\chi_{-7}}(m)=\tfrac17\operatorname{Tr}^{\mathbf1}_{-7m^2}(\widehat f)$,
both **[verified, $m\le12$, 45 digits]** (`65_trace7disc.gp`) — so the genus character in
the trace and the character $\psi$ in the Möbius inversion are exactly dual, and the
*master* divisibility (which needs $\psi=\chi_{D_0}$) is the divisibility of the
**untwisted** CM trace by $n^2$.  Writing $\Phi=A_2(\tau-\tau_0)^{-2}+\cdots$, one computes
$A_2(F_{4a})=\Delta(\rho)/E_4'(\rho)^2=1/(768\pi^2)$ and
$A_2(F_{4b})=E_4(i)\Delta(i)/E_6'(i)^2=-1/(1728\pi^2)$, so in both cases
$$\lambda=4\pi^2A_2 .$$
For Cooper's rows `cooper_sources` §2.3 gives $A_2=\nu^2/(4\pi^2g'(u_0))$, hence
$$|\lambda|=\frac{\nu^2}{|g'(u_0)|}=\sqrt3\ (s_7),\qquad 2\ (s_{10}),\qquad
\tfrac2{\sqrt3}\ (s_{18}),$$
**[verified]** the rule is exact for $s_7$ ($|\lambda|=\sqrt3$) and for $s_{10}$ at even $m$
($|\lambda|=2$); at odd $m$ for $s_{10}$ it is off by the factor $2$ explained in §2.4 by
$W_2$.  For $s_{18}$ it is off by $8$ (resp. $4$) — the observed constants are
$1/(4\sqrt3)$ and $1/(2\sqrt3)$ against the predicted $2/\sqrt3$.  All the discrepancies are
powers of $2$ and follow the same $2\mid m$ / $3\mid m$ splits, so they are presumably the
same Atkin–Lehner phenomenon compounded by the non-fundamental conductor of $-36$, but this
has **not** been derived; it is recorded as an open point.  (The sign depends on which of
the two Fricke-conjugate CM points the $\beta$-class selects.)

**This is the decisive structural point.**  The trace identity is a *theta lift* — it is the
Bruinier–Funke/Katok–Sarnak correspondence "integral weight $2-2k$ $\to$ half-integral weight
$k+1/2$" at $k=2$, whose input weight is $-2$ and output weight $5/2$, exactly the
Shimura–Borcherds input of Paşol–Zudilin.  It holds for Cooper's rows and for
Paşol–Zudilin's level-one forms alike.  What is special about Cooper's rows is therefore
**not** the shape of $\gamma$ but the *arithmetic of the traces*: $m^2\mid\beta(m)$ holds
for them and fails at once for the level-one forms ($\beta_{F_{4a}}(2)=-253$).

### 2.6 The master conjecture as a divisibility of traces of singular moduli

$$\boxed{\ n^2\mid\beta(n)\ \iff\ n^2\ \Bigm|\ \lambda\operatorname{Tr}_{D_0n^2}(\widehat f)\ }$$
i.e. exactly the Ahlgren–Ono/Edixhoven phenomenon — divisibility of traces of singular
moduli at discriminants $D_0p^{2k}$ by powers of $p$ — but for the **first Shimura–Maass
derivative of a weight $-2$ form** in place of a modular function.  Route (B) of the task is
now a precisely stated problem, and its $p$-adic ingredients are available: the disc-$D_0n^2$
CM points lie in the $p$-adic residue discs of the disc-$D_0$ points, and on the ordinary
locus $R_{-2}$ agrees at CM points with the Katz $\theta$-operator under the Katz $p$-adic
period comparison, so the "one derivative" of $R_{-2}$ is expected to contribute the extra
factor $p$ on top of the factor $p$ that Edixhoven's argument gives for a weight-$0$ input —
which is precisely the observed $p^2$.

### 2.7 The theta lift, verified in full **[verified, $d\le40$, 45 digits]**

`66_thetalift.gp`.  Building Paşol–Zudilin's weight-$5/2$ form
$f_{4a}=\tfrac78g_0+\tfrac1{768}g_1-\tfrac1{768}g_2=\tfrac1{64}q^{-3}+q-506q^4+
\tfrac{131565}{64}q^5-66516q^8+180249q^9-\cdots$ from their §4 and comparing with the
twisted traces of $\widehat f$ at **all** admissible discriminants:
$$\boxed{\ a(d)=\frac{\sqrt d}{192}\,\operatorname{Tr}_{-3d}(\widehat f)\qquad
\text{for every }d\equiv0,1\pmod4,\ d>0.\ }$$
Verified for $d=1,4,5,8,9,12,13,16,17,20,21,24,25,28,29,32,33,36,37,40$ to $45$ decimals,
including the *non-square* $d$ (where both sides are irrational multiples of $\sqrt d$:
e.g. $a(5)=131565/64$, $a(13)=189072723/32$).

This is the Bruinier–Funke/Katok–Sarnak theta lift "integral weight $2-2k\to$
half-integral weight $k+1/2$" at $k=2$, written out explicitly and checked against the
literature: **the weight-$5/2$ Shimura–Borcherds input of a magnetic weight-$4$ form is the
generating function of the twisted CM traces of the Maass-raised weight-$(-2)$ form.**  Two
consequences:

* the square-index case $d=m^2$ is exactly $a(m^2)=m\beta(m)$ of §4.1, so the master
  conjecture $m^3\mid a(m^2)$ **is** $m^2\mid\operatorname{Tr}_{D_0m^2}(\widehat f)$;
* it gives a *constructive* handle on the level-$4N$ weight-$5/2$ input for Cooper's rows,
  which `FINDINGS_PZ.md` found hard to reach with PARI's `mf` package directly: its
  coefficients are $\sqrt d\,\lambda\operatorname{Tr}_{D_0d}(\widehat f)$, computable from
  $\widehat f=Df+f/(2\pi y)$, $f=1/(xF)$.

### 2.8 The weight-$5/2$ input of $s_7$, computed **[verified, $d\le88$]**

`68_f52_s7.gp`, `69_f52_s7b.gp`.  Running the trace at *all* admissible discriminants rather
than only at $D_0m^2$ produces the weight-$5/2$ input of Cooper's row $s_7$ directly.  Put
$$c(d):=\sqrt{-3}\sum_Q\chi_{-3}(Q)\,\frac{\widehat f(\alpha_Q)}{\omega_Q},$$
the sum over $\Gamma_0(7)$-classes of Heegner forms of discriminant $-3d$ with $7\mid A$ and
$B\equiv\beta\ (14)$, $\beta$ the class $\equiv5\sqrt d$ when $d$ is a square.  Then
**$c(d)\in\mathbf Z$ for every admissible $d$** and $c(m^2)=\beta_{s_7}(m)$:

| $d$ | 1 | 4 | 8 | 9 | 16 | 21 | 25 | 28 | 29 | 32 | 36 | 37 | 44 | 49 | 53 | 56 | 57 | 60 | 64 | 65 | 72 | 77 | 81 | 84 | 85 | 88 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| $\pm c(d)$ | 1 | 8 | 16 | 9 | 32 | 0 | 125 | 0 | 145 | 192 | 216 | 148 | 176 | $\ast$ | 795 | 0 | 798 | 1200 | 1152 | 650 | 1008 | 0 | 3159 | 0 | 2975 | 4224 |

($\ast$: at $d=49$ the two $\beta$-classes merge and the trace degenerates, as at $7\mid m$
in the square case, where $\beta=0$ by the Atkin–Lehner cell; the other multiples of $7$ give
$0$.)  The admissible $d$ are exactly $d\equiv0,1\pmod4$ **and** $d\equiv0,1,2,4\pmod7$
(i.e. $-3d$ a square mod $28$) — the extra condition mod $7$ being the level-$28$ analogue of
the Kohnen plus condition.  This gives a *constructive* route to the level-$4N$ input, which
`FINDINGS_PZ.md` found inaccessible to a blind search in $\{g/\Delta(4\tau)^r\}$.

---

## 3. What $\gamma$ is **not**: three routes closed

### 3.1 Not a CM trace of a modular function **[proved]**

`NOTES_TRACE.md`, `02_trace7.gp`–`07_obstr7.gp`.  For a modular function $G$ on $X_0(N)$
with poles only at the cusps, the Heegner point of largest height at discriminant $D_0m^2$
is $\alpha=-m\bar\tau_0$ (height $m\operatorname{Im}\tau_0$), together with its Fricke
partner, all other classes having height $\le\frac12m\operatorname{Im}\tau_0$; hence
$$\operatorname{Tr}^{\mathrm{tw}}_m(G)=\sum_{j\ge1}2\operatorname{Re}\bigl(c_jq_0^{-jm}\bigr)
+O(R^{\theta m}),\qquad \theta<1,$$
with **constant** $c_j$ (the principal-part coefficients of $G$).  By `cooper_sources` §2.4
the polar data of $\Phi$ give
$\beta(m)=-\nu^2\,2\operatorname{Re}\bigl(q_0^{-m}/g'(u_0)\bigr)\bigl(1-\kappa/m\bigr)
+O(R^{\theta m})$ with $\kappa=1/(2\pi\operatorname{Im}\tau_0)$.  The $1/m$ tail has size
$R^m/m\gg R^{\theta m}$ and is not an exponential polynomial with constant coefficients, so
$\beta$ is not a finite $\mathbf Q$-linear combination of such traces.  $\square$

Numerically, with $H=1/u-49u$ the unique Fricke-anti-invariant Hauptmodul combination with
poles of order $\le1$, $\operatorname{Tr}^{\mathrm{tw}}_m(H)=i\sqrt3\,\tau_H(m)$ with
$\tau_H(m)\in\mathbf Z$ ($1,-6,5,18,-57,90,\ast,-458,1227,-1138,\dots$) and
$$\frac{m\,\beta(m)}{(m-\kappa)\,\tau_H(m)}\longrightarrow3
\quad(3.000000000012\ \text{at }m=33),$$
while an exact three-parameter fit against the full basis $\{1,1/u-49u,1/u+49u\}$ solved at
$m=22,23,24$ misses $m=25,\dots,34$ by $\approx1\%$ — a relative error of size $\kappa/m$, not
$R^{-m/2}$.  The $E_2^*$-completed control also fails: numerically
$\sum_Q\chi\omega^{-1}(E_2^*H/F)(\alpha_Q)=-3\operatorname{Tr}^{\mathrm{tw}}_m(H)+O(1)$,
producing no new sequence (`06_e2trace7.log`).

The moral is exactly §2: the missing $1/m$ is supplied by the $f/(2\pi y)$ term of the Maass
raising operator, i.e. by the fact that $\Phi$ has a **double** pole and $\Xi=D^{-1}\Phi$ is
an Eichler integral.

### 3.2 $K$ is not modular and not holonomic **[proved]**

$|\gamma(n)|\sim c\,R^n/n^2$ with $R=e^{2\pi\operatorname{Im}\tau_0}=e^{\pi\sqrt{|D_0|}/(2N)}
>1$ (`09_gamma_probe.log`, and $|\gamma(n)|^{1/n}\to R$ from below).  Hence:

* $K$ is not a weakly holomorphic modular form of any weight and level: those have
  coefficients $O(e^{c\sqrt n})$.  A fortiori not holomorphic, not quasi-modular.
* $K$ is not a *meromorphic* modular form either: those have coefficients
  $\mathrm{poly}(n)\cdot R^n$, with no $1/n^2$ factor.  Same for $DK$ ($\sim R^n/n$) and
  $D^2K=B=\sum\beta(n)q^n$ ($\sim R^n(1-\kappa/n)$).
* $K$ is not holonomic over $\mathbf Q(q)$: a $D$-finite series has its singularities among
  the roots of the leading polynomial, hence an algebraic radius of convergence, but
  $R^{-1}=e^{-\pi\sqrt{|D_0|}/(2N)}=e^{-\pi\sqrt d}$ with $d=3/49,\,1/25,\,1/9$ rational
  positive, which is transcendental by Gelfond–Schneider ($e^{\pi\sqrt d}=(-1)^{-i\sqrt d}$).

### 3.3 $\Phi^\flat$ is not modular **[proved]**

$\Phi^\flat:=\sum_{e\ge1}\mu(e)\psi(e)\,e\,\Phi(e\tau)=D^3K=DB$ has $[q^n]\Phi^\flat=n\beta(n)
=a(n^2)$ **[verified, $n\le400$, `08_struct.gp`]**, and $\Phi^\flat\approx\Phi$
asymptotically.  But $\Phi(e\tau)$ has poles at $\tau_0/e$, so $\Phi^\flat$ is meromorphic on
$\mathbf H$ with poles at infinitely many $\Gamma_0(N)$-inequivalent points (accumulating at
the cusp $0$), hence is not a modular form.  Consistently, $\Phi^\flat/\Phi$ is not a
rational function of $u$ of degree $\le6$ (`10b_flat.log`, $60$ exact equations, kernel $0$).

---

## 4. The dictionary and the equivalent shapes of the conjecture

### 4.1 $\beta(m)=a(m^2)/m$ **[proved; calibrated]**

If $\Phi=\Psi(f_{5/2})$ is the Shimura–Borcherds lift (Paşol–Zudilin (5), $k=2$, $D=1$) of
$f_{5/2}=\sum a(n)q^n$, i.e. $c(n)=\sum_{d\mid n}\psi(d)\,d\,a(n^2/d^2)$, then
$$c'(n)=\sum_{m\mid n}\frac{a(m^2)}m\quad\Longrightarrow\quad \beta(m)=\frac{a(m^2)}m,$$
so **the master conjecture is $m^3\mid a(m^2)$**, three powers more than Paşol–Zudilin's
theorem $m\mid a(m^2)$ delivers.  Calibration against the literature (`18_pzcontrol.gp`):
for $F_{4a}$ one gets $\beta=1,-253,60083,-14090996,\dots$, hence
$a(m^2)=m\beta(m)=1,-506,180249,\dots$, and $a(1)=1$, $a(4)=-506$ agree with Paşol–Zudilin's
$64f_{4a}=q^{-3}+q-506q^4+\cdots$ (Lemma 1(a)).  See also `FINDINGS_PZ.md`.

*Caveat on the character.*  The identity $\beta(m)=a(m^2)/m$ is for the untwisted lift
($\psi=\mathbf1$, $D=1$).  For the $\psi$-twisted lift
$c(n)=\sum_{d\mid n}\psi(d)\,d\,a(|D|n^2/d^2)$ with $\chi_D=\psi$ the same computation
gives $\beta_\psi(m)=a(|D|m^2)/m$, so with $\psi=\chi_{D_0}$ (the character for which the
master congruence holds at level one, §5) the conjecture reads $m^3\mid a(|D_0|m^2)$ — the
form stated in `COOPER_CONGRUENCE.md` §4.

### 4.2 Equivalent shapes **[proved]**

With $K=\sum\gamma(d)q^d$ and $\Lambda_\psi$ the $\psi$-Lambert kernel, the following are
equivalent to $\gamma\in\mathbf Z$ (i.e. to $n^2\mid\beta(n)$):

1. $\Xi=\sum_{e\ge1}\psi(e)\,(D^2K)(q^e)$, i.e. $c(n)=n\sum_{d\mid n}\psi(n/d)d^2\gamma(d)$
   — $\Xi$ is a **$\gamma$-twisted weight-3 Eisenstein series** $E_3^{\psi,\mathbf1}$ with the
   Lambert kernel replaced by $K$;
2. the **dilogarithm ladder** $\;D^{-3}\Phi=\sum_{n\ge1}\gamma(n)\,\mathrm{Li}_2^\psi(q^n)$,
   $\mathrm{Li}_2^\psi(t)=\sum_k\psi(k)t^k/k^2$; one derivative down,
   $D^{-2}\Phi=\log\mathfrak g$ with $\mathfrak g=\prod_n\bigl((1-q^n)^{n}\bigr)^{-\gamma(n)}$;
3. $\Phi^\flat=\sum_e\mu(e)\psi(e)e\Phi(e\tau)$ is **triply magnetic**
   ($D^{-3}\Phi^\flat\in\mathbf Z[[q]]$);
4. $c(pn)\equiv\psi(p)\,p\,c(n)\pmod{p^{3(1+v_p(n))}}$ for every prime $p$ and every $n$
   — the "eq:magnetic with multiplicity" form;
5. $m^3\mid a(m^2)$ for the weight-$5/2$ Shimura–Borcherds input;
6. $n^2\mid\lambda\operatorname{Tr}_{D_0n^2}(\widehat f)$ (§2.6).

(All six are formal rewritings; (1)–(3) verified to $n\le400$ in `08_struct.gp`, (4)
follows from (1) by the computation of round 1 Prop. 2.2 applied at each prime power.)

---

## 5. The master conjecture is not special to Cooper's rows

### 5.1 Round 1's negative control was the wrong character **[verified, $n\le500$]**

`60_pzcheck.gp`, `61_pznorm.gp`.  Round 1 §3 (= `COOPER_CONGRUENCE.md` §3, V6) tested
$\operatorname{rad}(n)^2\mid\beta(n)$ for Paşol–Zudilin's $F_{4a}=\Delta/E_4^2$ and
$F_{4b}=E_4\Delta/E_6^2$ **with $\psi=\mathbf1$** and found immediate failure
($\beta_{F_{4a}}(2)=-253$).  But the character attached to these forms is the genus
character of the **pole discriminant**: $\chi_{-3}$ for $F_{4a}$ (double pole at $\rho$)
and $\chi_{-4}$ for $F_{4b}$ (double pole at $i$).  With those characters:

| form | $\psi$ | scale | $n^2\mid\beta(n)$, $n\le500$ | $\min_n\bigl(v_p(\beta(n))-2v_p(n)\bigr)$, $p=2,3,5,7,\dots$ |
|---|---|---|---|---|
| $\Delta/E_4^2$ | $\chi_{-3}$ | $1$ | **no**, only at even $n$ | $-2,\,0,\,1,\,0,\,0,\dots$ |
| $\Delta/E_4^2$ | $\chi_{-3}$ | $64$ | **YES** | $4,\,0,\,1,\,0,\,0,\dots$ |
| $E_4\Delta/E_6^2$ | $\chi_{-4}$ | $1$ | **no**, only at multiples of $3$ | $0,\,-2,\,1,\,0,\dots$ |
| $E_4\Delta/E_6^2$ | $\chi_{-4}$ | $108$ | **YES** | $2,\,1,\,1,\,0,\dots$ |

The scales $64$ and $108$ are exactly Paşol–Zudilin's own constants — the ones making the
weight-$5/2$ input integral ($f_{4a}\in\frac1{64}q^{-3}\mathbf Z[[q]]$,
$f_{4b}\in\frac1{108}q^{-4}\mathbf Z[[q]]$, their Lemma 1) — and they cure exactly the bad
prime.  The divisibility is **sharp**: the minimum of $v_p(\beta(n))-2v_p(n)$ is $0$ for
every $p\ge7$.

**Therefore round 1's V6 and the sentence "no argument that only uses meromorphic weight $4$
with a double pole at a CM point can prove it" are wrong, and route (D) is reopened.**

### 5.2 The whole of Paşol–Zudilin's Table 1 **[verified, $n\le150$]**

`63_table1.gp`, `64_table1b.gp`.  Table 1 of `pz.txt` lists strong magnetic weight-$4$ forms
as $\Psi(c^{-1}f_m)$ with explicit expressions $E_4\cdot(\text{polynomial in }j)/H_{-m}(j)^2$.
Rescaling to $\Psi(f_m)$ and taking $\psi=\chi_{-m}$:

| form | $m$ | $h(-m)$ | $m\mid c(m)$ | $n^2\mid\beta(n)$ | sharp at $p\ge7$ |
|---|---|---|---|---|---|
| $\Psi(f_7)$ | $7$ | 1 | yes | **YES** | yes |
| $\Psi(f_8)$ | $8$ | 1 | yes | **YES** | yes |
| $\Psi(f_{11})$ | $11$ | 1 | yes | **YES** | yes |
| $\Psi(f_{19})$ | $19$ | 1 | yes | **YES** | yes |
| $\Psi(f_{43})$ | $43$ | 1 | yes | **YES** | yes |
| $\Psi(f_{67})$ | $67$ | 1 | yes | **YES** | yes |
| $\Psi(f_{163})$ | $163$ | 1 | yes | **YES** | yes |
| $\Psi(f_{15})$ | $15$ | 2 | yes | **YES** | yes |
| $\Psi(f_3|T_4)$ | $3$ | 1 | yes | **YES** | yes |
| $\Psi(f_4|(1-\tfrac12T_4))$ | $4$ | 1 | yes | **YES** | yes |

($\Psi(f_{20})$ and $\Psi(f_{23})$ already fail magnetism at $m=7$ in our transcription; the
numerator coefficients of those two entries are ambiguous in the plain-text dump of the
paper, so we treat them as untested, not as counterexamples.)

**A control.**  Magnetism is a real constraint, not automatic for forms with the right
polar divisor: the naive family $\Phi_D=E_4/H_D(j)^2$ — weight $4$, holomorphic at the
cusp, double poles exactly on the disc-$D$ Heegner divisor — is **not** magnetic for
$D=-7,-8,-11,-19,-43,-67,-163,-15$ (already $2\nmid c(2)=1$), `62_family.gp`.  So the
statement really is about Shimura–Borcherds lifts, i.e. about the polar divisor *together
with* the correct linear combination.

### 5.3 The mechanism, and why Cooper's rows have $\psi=\mathbf1$

`FINDINGS_PZ.md` establishes, for $f\in S^{!,+}_{5/2}$ at level $4$ with single-index
principal part $c\,q^{-m_0}$ (so that $\dim S_{5/2}=0$ forces uniqueness):
$$f|T_{p^2}=\chi_{-m_0}(p)\,p\,f+p^3c\,g_{m_0p^2},\qquad
g_{m_0p^{2r}}|T_{p^2}=g_{m_0p^{2r-2}}+p^3g_{m_0p^{2r+2}}\ (r\ge1),$$
**exactly** (not merely as a congruence), and deduces
$$v_p\bigl(a(p^{2j}m^2)\bigr)\ \ge\ 3j+v_p\bigl(a(m^2)\bigr)\ \text{ for all }j
\quad\Longleftrightarrow\quad \lambda\equiv p\pmod{p^3},$$
with the exact formula $a(p^{2j}m^2)=p^{3j}h_j(m^2)$ under that hypothesis.  Since
$\lambda=\chi_{-m_0}(p)p$, a *single-index untwisted* input has $n^3\mid a(n^2)$ exactly at
the split primes — which is precisely the observation that $\psi=\chi_{D_0}$ is the right
character for the level-one forms, since the $\psi$-twisted Möbius inversion contributes a
second factor $\chi_{D_0}(p)$ and turns the eigenvalue into $\chi_{D_0}(p)^2p=p$.

For Cooper's rows $\psi=\mathbf1,\mathbf1,\chi_{-3}$ while $D_0=-3,-4,-36$, so the character
is *not* $\chi_{D_0}$; consistently, the genus character $\chi_{D_0}$ appears instead **inside
the trace** (R1) — and the two are dual (R4).  In other words, the $\chi_{-3}(Q)$ twist in
$$\beta_{s_7}(m)=\sqrt{-3}\sum_Q\chi_{-3}(Q)\,\omega_Q^{-1}\widehat f(\alpha_Q)$$
is not cosmetic: it is exactly what converts the eigenvalue $\chi_{-3}(p)p$ into $p$ and so
delivers the divisibility at *every* prime rather than only at the split ones.

### 5.4 The level-one case is a theorem **[proved]**

`FINDINGS_PZ.md` §2.7, numerics `37_thm.gp/.log`.  Let
$f=c\,q^{-m_0}+O(q)\in S^{!,+}_{5/2}(\Gamma_0(4))$ with $-m_0$ a discriminant, let
$\{g_M=q^{-M}+O(q)\}$ be Borcherds–Zagier's integral basis of the Kohnen plus space, put
$\Phi=\Psi(f)=\sum A(n)q^n$, $\psi=\chi_{-m_0}$, $\beta_\psi=(A/n)\star(\mu\psi)$.  Fix a
prime $p$ with $c\in\mathbf Z_p$ and $p^2>m_0$, and $m$ prime to $p$; write
$x_{r,j}=c\,g_{m_0p^{2r}}(p^{2j}m^2)$.  Because $\dim S_4(\mathrm{SL}_2(\mathbf Z))=0$, a
weakly holomorphic plus form is determined by its principal part, which forces the two
**exact** identities
$$f|T_{p^2}=\chi_{-m_0}(p)\,p\,f+p^3c\,g_{m_0p^2},\qquad
g_{m_0p^{2r}}|T_{p^2}=g_{m_0p^{2r-2}}+p^3g_{m_0p^{2r+2}}\ (r\ge1).$$
A two-step induction (the key lemma being
$x_{r,k}-x_{r-1,k-1}=p^{3k}x_{r+k,0}-p^{3k-2}x_{r+k-1,0}$) yields the closed formula
$$\boxed{\ a(p^{2a}m^2)+\bigl(1-\psi(p)\bigr)\sum_{j=0}^{a-1}p^{a-j}a(p^{2j}m^2)
= p^{3a}\,c\,g_{m_0p^{2a}}(m^2)\ }$$
uniformly in $\psi(p)=+1,-1,0$.  Since $A/n=\mathbf1\star\alpha$ with $\alpha(e)=a(e^2)/e$
and $\beta_\psi=\alpha\star\nu$, $\nu=\mathbf1\star\mu\psi$, $\nu(e)=\prod_{p\mid e}(1-\psi(p))\in\mathbf Z$,
the left-hand side is exactly what $\beta_\psi$ assembles, and $v_p\ge3a$ gives
$v_p(\beta_\psi(n))\ge2v_p(n)$; over all $p$, $n^2\mid\beta_\psi(n)$.  $\square$

At split $p$ the correction terms vanish and one gets $p^{3a}\mid a(p^{2a}m^2)$ outright; at
inert $p$ ($1-\psi=2$) and ramified $p$ ($1-\psi=1$) the correction terms are precisely what
the $\psi$-twisted Möbius inversion supplies — which is *why* $\beta_\psi$ and not
$a(m^2)/m$ is the right object.  **[verified]** $v_p(S_a)\ge3a$ for all $a\ge1$ and all
$p^{2a}m^2\le3000$, every $2\le p\le53$, except the single prime dividing the denominator of
$c$ ($p=2$ for $f_{4a}$, $c=1/64$; $p=3$ for $f_{4b}$, $c=-1/108$) — exactly the constants
$64$ and $108$ of §5.1; and the excess is $0$, so the result is sharp.

The theorem covers every entry of Table 1 (all are $\Psi(\kappa f_{m_0})$), so §5.2 is
explained.  It does **not** cover Cooper's rows: at level $4N$ (i) $\dim S^+_{5/2}(\Gamma_0(28))=1$,
so the principal part no longer determines the form; (ii) there are six cusps, so the
principal-part datum is a vector, not a single index, and $T_{p^2}$ mixes them; (iii)
Borcherds–Zagier's integral basis $\{g_M\}$ is a level-$4$ statement.  See §8.

### 5.5 The level-$28$ inputs: two of the three obstacles removed

`72_check.gp/.log`, `FINDINGS_PZ.md` §§4.4, 4.5, 8.

**Uniqueness [proved].**  $\dim S_4(\Gamma_0(7))=1$; PARI gives the eigenform
$q-q^2-2q^3-7q^4+16q^5+2q^6-7q^7+\cdots$ with $a_7=-7$, hence Fricke eigenvalue
$-a_7/7^{k/2-1}=+1$ (confirmed by `mfatkineigenvalues`).  Since $\Phi|_4W_7=-\Phi$ — and
$\widehat f|W_7=-\widehat f$ is precisely what forces the trace onto one $\beta$-class — the
object lives in the $W_7=-1$ eigenspace, where
$$\dim S_4(\Gamma_0(7))^{W_7=-1}=0 .$$
So a weakly holomorphic form is determined by its principal part there, exactly as
$\dim S_4(\mathrm{SL}_2(\mathbf Z))=0$ does at level one, and Paşol–Zudilin's uniqueness step
survives at level $28$ with no vector-valued dimension formula.  (Note
$\dim M_4(\Gamma_0(7))=3$, so the anti-invariant part of $M_4$ is *not* zero; the argument
also uses that $\Phi$ vanishes at both cusps, which it does.)

**The eigenvalue [verified, $p\le19$].**  Fitting
$c(p^2d)+\tau\,p\,c(d)+p^3c(d/p^2)=\lambda\,c(d)$ against the table of §2.8 with $\tau=+1$
(no Kronecker symbol) returns $\lambda\equiv p\pmod{p^2}$, **constant in $d$**, at
$p=2,3,5,11,13,17,19$ — the genus character squaring away, $\chi_{-3}(p)^2=1$, exactly the
predicted mechanism.  In the $c$-normalisation $c(m^2)=\beta(m)$ (no extra factor $m$, unlike
the level-one $a(m^2)=m\beta(m)$), so the target is $v_p(c(p^{2a}m^2))\ge2a$, for which
$\lambda\equiv p\ (p^2)$ is the exact requirement.

**The sharpest statement.**  Independently of any Hecke normalisation, the table gives
$$p^{2a}\ \bigm|\ c\bigl(p^{2a}d\bigr)\qquad p\ \text{prime},\ a\ge1,\ d\ \text{admissible},$$
**[verified: $295$ tests, $0$ failures, $p\le29$, $a\le4$, $p^{2a}d\le2500$]**, sharp
(minimal excess $0$ at $p=3,11,13,17,19,23,29$), and including the $118$ cells with
$p\mid d$ — which `FINDINGS_PZ.md` could not reach — which pass with minimal excess $3$ at
$p=2$ and $2$ at $p=3,5$.  So the divisibility is not an artefact of restricting to $d$
prime to $p$.

Collecting the primes, this is the single clean assertion
$$\boxed{\ d\ \bigm|\ c(d)\ \text{ for every admissible }d\ }
\qquad\textbf{[verified: all 690 admissible }d\le2500\textbf{, 0 failures]}$$
(`75_dcd.gp/.log`), of which **Conjecture 4.1 is the case $d=m^2$**: $c(m^2)=\beta_{s_7}(m)$,
and $m^2\mid\beta(m)$ *is* $m^2\mid c(m^2)$.  The quotient
$$G(d):=c(d)/d\in\mathbf Z$$
extends Cooper's kernel from the squares to all discriminants, $G(m^2)=\gamma_{s_7}(m)$
**[verified, $m\le20$]**:
$$G=1,\,-2,\,2,\,1,\,2,\,0,\,-5,\,0,\,5,\,6,\,6,\,4,\,4,\,15,\,0,\,14,\,20,\,-18,\,10,\,14,\,0,\,39,\dots$$
at $d=1,4,8,9,16,21,25,28,29,32,36,37,44,53,56,57,60,64,65,72,77,81,\dots$  It is **sharp**:
$\min_{p\mid d}v_p(G(d))=0$ for $p=3,11,13,17,19,23,29$ and $=1$ for $p=2,5$ — precisely
round 1's two anomalous cells $(s_7,p=2)$ and $(s_7,p=5)$, now explained as a property of
the whole trace sequence rather than an accident at $\beta(2),\beta(5)$.

*Caveat on generality.*  "$d\mid c(d)$" is a statement about the **one-$\beta$-class**
(vector-valued component) trace.  At level one the full trace at non-square $d$ is an
irrational multiple of $\sqrt d$, so the statement does not even parse there; it is the
$\Gamma_0(N)$, $N>1$ phenomenon.  Whether it holds for $s_{10}$ and $s_{18}$ has not been
tested and is the natural next check.

---

## 6. The mod-$p$ structure of $\gamma$

Full details in `FINDINGS_MODP.md` (scripts `20_`–`29_`).  Summary of what is new:

* **Mod-$2$ law.**  $\gamma(n)$ is odd $\iff$ $2^{1+v_2(N)}\nmid n$ and $P\nmid n$, with
  $P=7,5,3$ the distinguished prime of the row.  Hence $\gamma\bmod2$ is *multiplicative*
  (although $\gamma$ is not) and $K\bmod2$ is a *rational* function of $q$ with denominator
  $1-q^{14},1-q^{20},1-q^{12}$. **[verified, $n\le12000$]**
* **Mod-$3$ law for $s_{18}$** ($\psi(3)=0$): $K_{s_{18}}\equiv(q+2q^2+2q^3)/(1+q)^3\pmod3$,
  period $6$; in particular $3\nmid\gamma_{s_{18}}(n)$ ever, and
  $\gamma(n)\gamma(n+1)\equiv-1\pmod9$ for $n\equiv1\ (3)$, sharp.
* Those four cells are the **only** ones: elsewhere $K\bmod p$ is neither rational
  (degree $\le250$) nor algebraic (boxes $(40,6),(20,20),(60,4)$), and the $p$-kernel is
  infinite, so Christol fails. **[refuted, $p\le43$]**
* Every Lucas/Dwork law for $\gamma$ fails massively; round 1's $(S)$ does **not** descend
  to $\gamma$.  $\gamma(p)\bmod p$ (= round 1's $\lambda(p)$) matched nothing in a sweep of
  several hundred arithmetic candidates over $42$ primes, and is equidistributed.
* **Mod $2$ the whole conjecture degenerates**: $c'(m)$ is odd $\iff$ the $2P$-free part of
  $m$ is a perfect square — a weight-$3$ Eisenstein/theta identity.  Proving it gives
  eq:magnetic at $p=2$ for all three rows unconditionally.  This is the most concrete
  finite-looking target produced by this round.
* New unbounded divisibilities $v_3(\gamma_{s_7}(n))\ge v_3(n)-1$,
  $v_2(\gamma_{s_{10}}(n))\ge v_2(n)-1$ (both sharp): $n^2\mid\beta(n)$ is not sharp at
  $(s_7,3)$ and $(s_{10},2)$, neither of which is an Atkin–Lehner cell.
* $\gamma_{s_7}(n)=0\iff7\mid n$ and $\gamma_{s_{10}}(n)=0\iff5\mid n$ — the "$\Leftarrow$"
  half **proved** from round 1's exact cells; $\gamma_{s_{18}}$ never vanishes and has sign
  $(-1)^{n-1}$.

### 6.1 Further verifications (all three rows, `FINDINGS_ROWS.md`)

* $n^2\mid\beta(n)$ exactly to $n\le12000$; $\gamma(n)=0$ exactly for $7\mid n$ ($s_7$),
  $5\mid n$ ($s_{10}$), never ($s_{18}$).
* Sharpness: $\min_n v_p(\gamma(n))=0$ for every $p\le60$ in every row; the minimum over
  $p\mid n$ is $1$ only at $(s_7,p=2)$ and $(s_7,p=5)$ — round 1's two anomalous cells.
* The refined supercongruences of `COOPER_CONGRUENCE.md` §4:
  $c'(p^k)\equiv\psi(p)c'(p^{k-1})\pmod{p^{2k}}$, **no failures** over $1493$ pairs $(p,k)$
  with $p^k\le12000$ per row; and
  $c'(pq)-\psi(p)c'(q)-\psi(q)c'(p)+\psi(pq)\equiv0\pmod{p^2q^2}$, **no failures** over
  $3101$ pairs $pq\le12000$ per row.  The modulus $p^{2k}$ is sharp for every $p$ **except**
  $p=2,3,5$, where $v_p$ of the difference grows like $3k-1$ or $2k+1$ uniformly up to
  $k=13$ — a stronger congruence at the small primes, not previously recorded.
* **Borcherds-product test refuted**: $\prod_n(1-q^n)^{\tau_\chi(n)}$, with $\tau_\chi$ the
  weight-$0$ twisted trace of §3.1, is **not** a meromorphic modular form on $\Gamma_0(N)$
  with divisor supported on the cusps and the CM points, in any of the three rows (exact
  linear algebra on its logarithmic derivative against
  $\{1\}\cup\{E_2(d\tau)\}\cup\{Fu^i/g^j\}$, with a positive control that recovers $E_2-1$).
* The weight-$0$ obstruction of §3.1 is measured independently in **all three** rows: the
  residual $\sup_m|(r-L)R^{m/2}|$ is $5.4\cdot10^5,\,2.6\cdot10^4,\,3.2\cdot10^4$ **without**
  the $\kappa$-correction against $0.29,\,0.39,\,0.0028$ **with** it, and the $1/m$
  coefficient equals $-L\kappa$ to $11$–$13$ significant digits.

---

## 7. The $x$-side: what became of round 1's "smallest missing brick"

Full details in `FINDINGS_XSIDE.md` (scripts `40_`–`47_`).

* Route (A) **as posed in the task is empty**: in weight $0$ the Hecke correspondence is
  $T_p=p\,U_p+V_p$, so the Kronecker congruence $\Phi_p(X,Y)\equiv(X^p-Y)(X-Y^p)\bmod p$
  says $p\,U_pf\equiv0$, a tautology.  One needs the modular polynomial mod $p^2$. **[proved]**
* Its **correct** form is a genuine reduction.  With $H=(\sqrt PF)^{p-1}$,
  $f|U_p\equiv f_0+\sum_{n\ge1}[x^{pn}](fH)x^n\pmod p$, and the Lucas law
  $A_{pj+r}\equiv A_jA_r$ (proved, Malik–Straub) gives $H\equiv P^{(p-1)/2}/F_{<p}$.  So the
  brick is the **finite explicit** congruence
  $$\bigl[x^{pn}\bigr]\Bigl(l(x)\,P(x)^{(p-1)/2}\big/F_{<p}(x)\Bigr)\equiv\psi(p)L_n\pmod p .$$
  **[verified, all $p\le53$, $x^{1200}$, all rows]**
* $\deg F_{<p}=\lfloor\mu(N)(p-1)/(6\deg x)\rfloor$ **exactly** for all $20$ primes
  $p\le79$ and all three rows: $F_{<p}$ is the **supersingular polynomial** of
  $X_0(N)/\mathbf F_p$ in the coordinate $x$ (the Hasse invariant of the family).  The brick
  is an ordinarity/unit-root statement at the supersingular divisor.
* New supercongruence $A_{p-1}\equiv\kappa\psi(p)p\pmod{p^2}$, $\kappa=-2,-3,-3$,
  **[verified $5\le p\le300$]**, **[proved for $s_{10}$]** — for $s_{18}$ it manufactures
  $\chi_{-3}$ from the Apéry-like recurrence alone, with no modular input.
* The $\psi$-twisted Dieudonné–Dwork tower $a_{pn-1}\equiv\psi(p)a_{n-1}\pmod{p^{1+v_p(n)}}$
  holds for all three rows (round 1 had only the $\psi=\mathbf1$ case).
* Negatives: no binomial-sum representation of $s_{18}$ was found; Beukers' congruence
  $A_{mp^s}\equiv A_{mp^{s-1}}\ (p^{3s})$ holds for $s_7,s_{10}$ but only to $p^{2s}$ for
  $s_{18}$; Lucas fails for $L_n$, $a_j$, $b_j$; $\eta_1$ has no description in the span of
  the obvious differentials.

---

## 8. The gap, stated three ways

$n^2\mid\beta(n)$ is still open.  It is now equivalent to each of:

1. **(Traces)** $n^2$ divides $\lambda\operatorname{Tr}_{D_0n^2}(\widehat f)$, the
   genus-character-twisted Heegner trace of the first Shimura–Maass derivative of the
   weight-$(-2)$ form $f=1/(xF)$.  This is an Ahlgren–Ono/Edixhoven statement; what is
   needed beyond the literature is the **weight-$5/2$, twisted, exact-eigenvalue** version
   (the published results are weight $3/2$, untwisted, and give a congruence mod $p$
   conditional on splitting).
2. **(Hecke)** the weight-$5/2$ input $f_{5/2}$ satisfies $f|T_{p^2}=\lambda f+p^3h_1$ with
   $\lambda\equiv p\pmod{p^3}$, together with the tower $h_r|T_{p^2}=h_{r-1}+p^3h_{r+1}$
   and $p$-integrality of the tower.  At level $4$ with a single-index principal part this
   is now a **theorem** (§5.4), so the *entire* remaining gap for Cooper's rows is the
   twisted analogue of those two identities.  §2.8 pins down the setting: the support law
   for $c(d)$ ($-3d$ a square mod $28$, i.e. $d\equiv0,1\ (4)$ **and** $d\equiv0,1,2,4\ (7)$)
   is *not* a Kohnen plus condition inside a scalar space but the index condition of a
   **vector-valued** form for the Weil representation of the discriminant form $\mathbf Z/2N$
   — the $\beta$-class is the component index, and $\widehat f|W_N=-\widehat f$ is why the
   sum over both classes vanishes.  Consistently, a direct search **[refuted]** finds no
   scalar element of $\{g/\Delta(4\tau)\in M_{29/2}(\Gamma_0(28))\}$ (dimension $17$) with
   $a(m^2)=m\beta_{s_7}(m)$ for $m\le40$, while the same pipeline recovers $f_{4a}$ exactly
   at level $4$.  So what is needed is: (i) vanishing of the space of holomorphic
   vector-valued weight-$5/2$ cusp forms for that Weil representation (the analogue of
   $\dim S_4=0$) — **done**, see P4: the Atkin–Lehner sign gives it, no dimension formula
   needed; (ii) the principal-part *vector*, which §2.8 makes explicit; and (iii) an
   **integral basis of $M^!_{5/2}(\rho_L)$ indexed by $(\beta,d)$** together with the tower
   relation (T2) written for $\rho_L$.  Only (iii) is genuinely missing.  Everything after
   that is the induction of §5.4 with target exponent $2a$ in place of $3a$.
3. **($x$-line, mod $p$ only)**
   $[x^{pn}]\bigl(l\,P^{(p-1)/2}/F_{<p}\bigr)\equiv\psi(p)L_n\pmod p$, with $F_{<p}$ the
   supersingular polynomial; plus the round-1 lift to mod $p^2$, which is now the single
   assertion "$\eta_1$ is exact".

Of these, (2) is the most likely to close: it is a finite-dimensional statement about a
space of half-integral weight forms, and the level-one case is already a theorem.

---

## 8bis. New phenomena recorded this round

1. **The trace identification itself** (R1, R1b, R1c) and its level-one generality (R4).
2. **The theta lift written out**: $a(d)=\sqrt d\,\operatorname{Tr}_{-3d}(\widehat f)/192$ at
   level one for *every* admissible $d$ (R4b) — the first explicit check we know of that
   Paşol–Zudilin's $f_{4a}$ is the trace generating function.
3. **The correction of the negative control** (A1, A2) and the resulting **theorem** at level
   one (A0).
4. **The mod-$2$ law** $\gamma(n)$ odd $\iff2^{1+v_2(N)}\nmid n$ and $P\nmid n$, and the
   mod-$3$ law for $s_{18}$: $\gamma\bmod2$ is multiplicative and $K\bmod p$ rational in
   exactly four cells and nowhere else (M1–M5).
5. **$\Xi\bmod2$ is Eisenstein**: $c'(m)$ odd $\iff$ the $2P$-free part of $m$ is a square
   (M10).
6. **Two extra unbounded divisibilities**: $v_3(\gamma_{s_7}(n))\ge v_3(n)-1$ and
   $v_2(\gamma_{s_{10}}(n))\ge v_2(n)-1$, sharp (M9).
7. **$\gamma(p)\bmod p$ is structureless**: equidistributed, matching none of several hundred
   arithmetic candidates (N1, N3) — round 1's $\lambda(p)$ question is answered negatively.
8. **The supersingular polynomial appears** as the denominator of the $x$-line kernel (W7).
9. **A new supercongruence** $A_{p-1}\equiv\kappa\psi(p)p\pmod{p^2}$ manufacturing $\chi_{-3}$
   from the recurrence alone (W11), proved for $s_{10}$.
10. **A stronger congruence at small primes**: the modulus $p^{2k}$ in
    $c'(p^k)\equiv\psi(p)c'(p^{k-1})$ is sharp for every $p$ except $p=2,3,5$, where the
    valuation grows like $3k-1$ or $2k+1$ up to $k=13$.
11. **Two negatives worth recording**: the Borcherds-product test on the weight-$0$ traces is
    refuted in all three rows; and the naive family $E_4/H_D(j)^2$ is not magnetic, so
    magnetism is a genuine constraint (A4).

---

## 9. Consequences for the repository

1. **`consolidation/COOPER_CONGRUENCE.md` §3 (the negative control) must be corrected.**
   $(S)$ does **not** fail for Paşol–Zudilin's level-one magnetic forms; it fails only for
   $\psi=\mathbf1$, which is the wrong character for them.  With $\psi=\chi_{D_0}$ and their
   own integralising constants it holds at every prime, and the same is true for all of
   their Table 1.  The sentence "no argument that only uses meromorphic weight $4$ with a
   double pole at a CM point can prove it" is false, and the "Routes closed" line about
   route (D) must be withdrawn.
2. Conjecture 4.1 should be **restated in the generality it actually has**: for
   $\Phi=\Psi(f)$ a Shimura–Borcherds lift of $f\in S^{!,+}_{5/2}$ with principal part
   $c\,q^{-m_0}$, one has $n^2\mid\beta_{\chi_{-m_0}}(n)$; Cooper's rows are the level-$N$
   twisted case, where the genus character makes $\psi=\mathbf1$ (resp. $\chi_{-3}$).
3. Round 1's V4 verification extends from $n\le1500$ to $n\le12000$.
4. Round 1's V19 (the Dieudonné–Dwork form of the brick) extends to the character row
   $s_{18}$ in the $\psi$-twisted form W20.
5. `cooper_sources/REPORT.md` §7's proposal of the Paşol–Zudilin route as "the natural
   attack" is **reinstated**: §3's objection was based on the wrong character.
6. The identification R1 supplies what `COOPER_CONGRUENCE.md` §5 calls "the lever":
   $K$ is the generating function of twisted CM traces of $\widehat f$.
7. `COOPER_CONGRUENCE.md` §4's Shimura–Borcherds reading "$a(|D|m^2)=m^3\gamma(m)$" is
   correct in the $|D|=D_k$ (i.e. $D=1$ for $k$ even) reading, **not** with $|D|$ the pole
   discriminant.
8. The plain-text dump `pz.txt` mangles Paşol–Zudilin's Lemma 1: the correct normalisations
   are $f_{4a}=\tfrac1{64}q^{-3}+q-506q^4+\cdots$ and $f_{4b}=-\tfrac1{108}q^{-4}+q+\cdots$
   (the $1/64$ and $1/108$ sit on the denominator line of the two-line display and read as
   part of the linear combination).
9. Conjecture 4.1 should be superseded by P6: **$d\mid c(d)$** for the twisted CM trace
   coefficients $c(d)$ of §2.8, of which it is the case $d=m^2$; the quotient $G(d)=c(d)/d$
   extends $\gamma$ to all discriminants.

---

## 10. Files

| file | contents |
|---|---|
| `lib.gp` | rows, eta quotients, $F$, $x$, $\Phi$, $c$, $c'$, $\beta$ (copied from round 1) |
| `heeg.gp` | reduced forms of a discriminant, $\Gamma_0(N)$-Heegner representatives with a fixed $\beta$-class, stabiliser orders, genus characters.  **Caveat:** `heegrep` tries only one completion $(q,s)$ per first column $(p,r)$ and so can return a non-minimal $A$; harmless for $s_7$ and $s_{10}$ (class counts match and the closed-form $\widehat f$ is height-insensitive) but fatal for $s_{18}$ — use `heegmin2` of `59_final3.gp` instead |
| `e2.gp` | $E_2^*$ by reduction to the $\mathrm{SL}_2(\mathbf Z)$ fundamental domain |
| `wt2.gp`, `wt2b.gp` | first versions of $f=1/(xF)$ and $\widehat f$ for $s_7$ (closed form; $q$-series) |
| `maass.gp`, `maass2.gp` | $\widehat f$ for all three rows, two independent evaluations, with the $q$-series fallback at the zeros of $F$ |
| `01_data.gp` | $c,c',\beta,\gamma$ for $n\le600$; $n^2\mid\beta(n)$ |
| `02`–`07` | twisted CM traces of modular *functions*; the double-pole obstruction (`NOTES_TRACE.md`) |
| `08_struct.gp` | the structural reformulations, exact, $n\le400$ |
| `09_gamma_probe.gp` | Dirichlet inverse and growth of $\gamma$ |
| `10`, `10b` | $\Phi^\flat$ is not modular |
| `11`–`13`, `17` | **the identification** for $s_7$ (`NOTES_IDENT.md`) |
| `14`, `15`, `16` | first attempts for $s_{10}$, $s_{18}$ (superseded by `FINDINGS_ROWS.md`) |
| `18`, `19` | the level-one control: the trace formula for $F_{4a}$, $F_{4b}$ |
| `20`–`29` | mod-$p$ structure of $\gamma$ (`FINDINGS_MODP.md`) |
| `30`–`39` | Paşol–Zudilin calibration, the $T_{p^2}$ eigen-identity and tower, level-$28$ feasibility (`FINDINGS_PZ.md`) |
| `40`–`47` | the $x$-side (`FINDINGS_XSIDE.md`) |
| `50`–`59` | the three rows: traces, obstruction, master verification to $n\le12000$ (`FINDINGS_ROWS.md`) |
| `60`–`64` | **the correction**: the master congruence for Paşol–Zudilin's forms and for their Table 1 |
| `beta_s*.txt`, `gamma_s*.txt`, `cp_s*.txt` and the `20*`, `22_`, `50_` variants | data |

Run with `gp -q <file>`; the CM evaluations (`11`–`19`, `50`–`59`) are the slow ones.
