> **Attribution (2026-08-22, found by River):** this theorem is **Beukers, "Irrationality proofs using modular forms", Astérisque 147–148 (1987), Theorem 3**: $\sum a_nq^n=(\eta_1^9\eta_6^9/\eta_2^3\eta_3^3)^{1/2}$ ($=4\Psi$), and $\sum a_n/n^2$ is irrational, proved via the $4^n$ denominators of $\sqrt E$ and $4e<(\sqrt2+1)^4$. Our contribution is the framing (square root of Apéry's recurrence), the explicit recurrence/Casoratian/source, the measure $\mu\le50.66$, and the formalization blueprint. PDF: papers/beukers1987.pdf.

# The square root of Apéry's recurrence: a formalization blueprint

*Standalone statement and proof, organised for Lean 4 / Mathlib autoformalization. Fable, 2026-08-22. Source of the mathematics: `consolidation/SQRT_APERY.md`, `paper/sections/02_sqrt_rows.tex`; every identity below was re-verified exactly (`lattice/sqrt_apery_*.gp`).*

## 0. Scope and design

**Main Theorem.** Define integers $a_n$ and rationals $b_n$ by
$$a_0=1,\ a_1=10,\qquad b_0=0,\ b_1=1,\qquad (n+1)^2u_{n+1}=(136n^2+68n+10)\,u_n-4(2n-1)^2\,u_{n-1}\quad(n\ge1).$$
Then $b_n/a_n$ converges to a real number $\xi=0.100187449229339406167758682133\ldots$, and $\xi$ is irrational.

The proof is elementary except for **one external input** (§4, integrality of the nome of Apéry's operator), which is a classical theorem with two independent published proofs. The modular-forms *interpretation* ($\xi=L(\Psi,2)$) is not needed and is not part of the blueprint.

Design choices for the formalization:
- Everything lives in $\mathbb Q[[u]]$ (formal power series, `PowerSeries ℚ`), with $\theta=u\,\frac{d}{du}$ (`PowerSeries.derivative` composed with multiplication by `X`). No analysis until §7, and there only real limits of rational sequences.
- Integrality statements are statements about $\mathbb Z[[u]]\subset\mathbb Q[[u]]$ and about $2$-adic valuations of rationals.
- The square root $f=\sqrt F$ is *defined* by coefficient recursion (§2), not by a binomial series; the binomial series is used only in the integrality lemma, as an identity of formal series.
- The irrationality criterion is the standard one (§7.4).

Notation: $\theta=u\frac d{du}$ on $\mathbb Q[[u]]$; $[u^n]g$ the $n$-th coefficient; $d_n=\operatorname{lcm}(1,\dots,n)$; $C_k=\frac1{k+1}\binom{2k}{k}$ the Catalan numbers; $v_2$ the $2$-adic valuation.

## 1. Apéry's numbers and operator (inputs from the literature, both formalizable)

**Definition 1.1.** $A_n=\sum_{k=0}^n\binom nk^2\binom{n+k}k^2$.

**Theorem 1.2 (Apéry's recurrence).** $(n+1)^3A_{n+1}=(2n+1)(17n^2+17n+5)A_n-n^3A_{n-1}$ for $n\ge1$, with $A_0=1,A_1=5$.
*Status:* classical; proved by creative telescoping (Zeilberger certificate; e.g. van der Poorten's account, or the WZ proof formalised in Coq by Chyla–Mahboubi–Sibut-Pinote–Tassi, "A formal proof of Apéry's theorem", 2014). A Lean proof can follow the Coq one: exhibit the rational certificate $R(n,k)$ and verify the telescoping identity.

**Corollary 1.3 (operator form).** Put $F=\sum_{n\ge0}A_n(4u)^n\in1+u\mathbb Z[[u]]$ and
$$L_{\rm Ap}:=\theta^3-4u\,(2\theta+1)(17\theta^2+17\theta+5)+16u^2(\theta+1)^3 .$$
Then $L_{\rm Ap}F=0$ in $\mathbb Q[[u]]$.
*Proof.* Coefficient extraction: $[u^{n+1}]$ of $L_{\rm Ap}F$ is $(n+1)^3A_{n+1}4^{n+1}-4\cdot4^n(2n+1)(17n^2+17n+5)A_n+16\cdot4^{n-1}n^3A_{n-1}=4^{n+1}\big((n+1)^3A_{n+1}-(2n+1)(17n^2+17n+5)A_n-\ldots\big)$; check the sign of the last term: $16u^2(\theta+1)^3$ applied to $A_{n-1}(4u)^{n-1}$ contributes $16\cdot4^{n-1}n^3A_{n-1}=4^{n+1}n^3A_{n-1}$ at $u^{n+1}$, and Theorem 1.2 reads $(n+1)^3A_{n+1}-(2n+1)(\ldots)A_n+n^3A_{n-1}=0$. $\square$

(The variable $u=t/4$ is used throughout; $t$ itself never appears in the proof.)

## 2. The square root and its operator

**Definition 2.1.** $f\in\mathbb Q[[u]]$ is the unique series with $f^2=F$ and $[u^0]f=1$: $f_0=1$, $f_n=\tfrac12\big([u^n]F-\sum_{k=1}^{n-1}f_kf_{n-k}\big)$.

**Definition 2.2.** $L_1:=P\theta^2+Q\theta+R$ with $P=1-136u+16u^2$, $Q=-68u+16u^2$, $R=-10u+4u^2$.

**Lemma 2.3 (symmetric-square identity).** For every $g\in\mathbb Q[[u]]$, with $w:=L_1g$,
$$L_{\rm Ap}(g^2)=2g\,\theta w+6\,(\theta g)\,w .$$
*Proof.* Polynomial identity in $g,\theta g,\theta^2g,\theta^3g$ after substituting $\theta^2g=(w-Q\theta g-Rg)/P$ and $\theta^3g=\theta$ of that; the denominators $P$ cancel identically. Verified by computer algebra (sympy; `lattice/sqrt_apery_ops.gp` for the numerical form). For Lean: expand both sides as polynomials in $g_0=g,g_1=\theta g,g_2=\theta^2g,g_3=\theta^3g$ with $w=Pg_2+Qg_1+Rg_0$, $\theta w=\theta P\cdot g_2+Pg_3+\theta Q\cdot g_1+Qg_2+\theta R\cdot g_0+Rg_1$ and compare; `ring` should close it. $\square$

**Proposition 2.4.** $L_1f=0$.
*Proof.* Let $w=L_1f$. By Lemma 2.3 and Corollary 1.3, $2f\theta w+6(\theta f)w=0$. Hence $\theta(wf^3)=f^2\big(f\theta w+3(\theta f)w\big)=\tfrac{f^2}2\big(2f\theta w+6(\theta f)w\big)=0$, so $wf^3$ is constant. Its constant term is $w_0f_0^3=w_0=[u^0](L_1f)=R(0)f_0=0$. Thus $wf^3=0$; $f$ is a unit ($f_0=1$), so $w=0$. $\square$

**Corollary 2.5 (the recurrence).** Define $a_n:=[u^n]f$ (equivalently $a_n=4^n[t^n]\sqrt{F}$ in the variable $t=4u$). Then $(n+1)^2a_{n+1}=(136n^2+68n+10)a_n-4(2n-1)^2a_{n-1}$, $a_0=1$, $a_1=10$.
*Proof.* $[u^{n+1}](L_1f)=0$: $P\theta^2$ gives $(n+1)^2a_{n+1}-136n^2a_n+16(n-1)^2a_{n-1}$; $Q\theta$ gives $-68na_n+16(n-1)a_{n-1}$; $R$ gives $-10a_n+4a_{n-1}$. Sum: $(n+1)^2a_{n+1}-(136n^2+68n+10)a_n+(16n^2-32n+16+16n-16+4)a_{n-1}$, and $16n^2-16n+4=4(2n-1)^2$. $a_1=[u^1]f=\tfrac12[u^1]F=\tfrac12\cdot20=10$. $\square$

So the sequence $a_n$ of the Main Theorem equals $[u^n]f$.

## 3. Integrality of $a_n$

**Lemma 3.1.** $4^k\binom{1/2}{k}=(-1)^{k-1}\,2\,C_{k-1}$ for $k\ge1$. *(Induction on $k$ from $\binom{1/2}{k+1}=\binom{1/2}{k}\frac{1/2-k}{k+1}$ and $C_k=C_{k-1}\frac{2(2k-1)}{k+1}$.)*

**Theorem 3.2.** For every $G\in1+x\mathbb Z[[x]]$, the series $\sqrt G$ (constant term $1$) satisfies $4^n[x^n]\sqrt G\in\mathbb Z$ for all $n$.
*Proof.* $\sqrt G=\sum_{k\ge0}\binom{1/2}{k}S^k$ with $S=G-1\in x\mathbb Z[[x]]$ (formal binomial series; this is the unique square root with constant term $1$, by uniqueness in Definition 2.1). $[x^n]S^k=0$ for $k>n$, and for $k\le n$, $4^n\binom{1/2}{k}[x^n]S^k=4^{n-k}\cdot4^k\binom{1/2}{k}\cdot[x^n]S^k\in\mathbb Z$ by Lemma 3.1. $\square$

**Corollary 3.3.** $a_n\in\mathbb Z$ for all $n$.
*Proof.* $f(u)=\sqrt{F}$ and $F=\sum A_n(4u)^n=G(4u)$ with $G(x)=\sum A_nx^n\in1+x\mathbb Z[[x]]$; so $[u^n]f=4^n[x^n]\sqrt G\in\mathbb Z$. $\square$

(Remark: the sharper statement $\lambda^n[x^n]\sqrt G\in\mathbb Z$ with $\lambda=\max(1,2^{2-e})$, $e=\min_{n\ge1}v_2([x^n]G)$, is not needed.)

## 4. The nome: the single external input

Work now in the variable $t=4u$ as well as $u$; $F\in1+t\mathbb Z[[t]]$ as a series in $t$ (Definition 1.1), and $f=\sqrt F$. Let $\sigma:=\sqrt{1-34t+t^2}=\sqrt{P}\in1+u\mathbb Q[[u]]$ (constant term $1$), so that $\sigma^2=P$ and, since $Q=\tfrac12\theta_uP$,
$$\theta_u\sigma=\frac{Q}{\sigma}\qquad\text{i.e.}\qquad \frac{\theta_u\sigma}{\sigma}=\frac QP. \tag{4.0}$$

**Definition 4.1 (nome).** $q\in t\,\mathbb Q[[t]]$ is the unique series $q=t+O(t^2)$ with
$$\theta_uq=\frac{q}{F\sigma}=\frac{q}{f^2\sigma}, \tag{4.1}$$
i.e. $\lambda:=\theta_u\log q=1/(f^2\sigma)$. (Existence and uniqueness: $\frac1{F\sigma}-1\in t\mathbb Q[[t]]$ determines $\log(q/t)$ coefficientwise; exponentiate. In Lean, define $q$ by the coefficient recursion that (4.1) gives. Note (4.1) is invariant under $t\leftrightarrow u=t/4$ since $\theta$ is; the normalisation $q=t+O(t^2)$ is the one under which Theorem 4.2 holds.) Let $u(q)$ (equivalently $t(q)=4u(q)$) denote the compositional inverse.

**Lemma 4.1a (elementary).** $\sigma=\sqrt{1-34t+t^2}\in1-17t+16t^2\mathbb Z[[t]]$; hence $F\sigma\in\mathbb Z[[t]]$ and $1/(F\sigma)\in\mathbb Z[[t]]$.
*Proof.* $1-34t+t^2=(1-17t)^2-288t^2$ and $288=2^5\cdot3^2$, so $\sigma=(1-17t)\sqrt{1+y}$ with $y=-288t^2/(1-17t)^2\in2^5t^2\mathbb Z[[t]]$, and $4^k\binom{1/2}k=\pm2C_{k-1}$ (Lemma 3.1) absorbs the $4^k$ in $\binom{1/2}ky^k$. $F\sigma$ is then a product of integral series with constant term $1$, hence a unit of $\mathbb Z[[t]]$. $\square$

**Theorem 4.2 (integrality of the nome).** $q(t)\in t\,\mathbb Z[[t]]$, its inverse $t(q)\in q\,\mathbb Z[[q]]$, and $F(t(q))\in1+q\,\mathbb Z[[q]]$.
*Status and proof structure* (`consolidation/NOME_INTEGRALITY.md`): write $1/(F\sigma)=\sum_{n\ge0}r_nt^n\in\mathbb Z[[t]]$ (Lemma 4.1a), so that $\log(q/t)=\sum_{n\ge1}r_nt^n/n$. By Dwork's lemma (proved in full there), $q/t\in1+t\mathbb Z_p[[t]]$ for all $p$ **if and only if** the Gauss congruences
$$r_{mp^a}\equiv r_{mp^{a-1}}\pmod{p^a}\qquad(\text{all primes }p,\ a\ge1,\ m\ge1)$$
hold. These are supplied by Krattenthaler--Rivoal, *S\'emin.\ Congr.* 27 (2011), arXiv:0804.3049, Theorem 2, specialised to Ap\'ery's numbers ($d=k=2$, $N^{(1)}=N^{(2)}=(2,1)$), whose proof uses only Dwork-type congruences for $A_n$ (Samol--van Straten / Mellit--Vlasenko for the Laurent polynomial $\Lambda=(x+y)(z+1)(x+y+z)(y+z+1)/xyz$ of which $A_n$ are the constant terms of powers); the bridge between their mirror map and Definition 4.1 is the identity $\theta(G/F)=1/(F\sigma)-1$ for the log-solution $G$ of Ap\'ery's operator, a WZ-certifiable recurrence of the same status as Theorem 1.2. Parts (ii),(iii) follow from (i) by Lagrange inversion and composition. Numerically verified: $q(t)$ integral to $t^{4000}$, Gauss congruences for $p\le200$, $n\le4000$ (10\,832 instances), $t(q),F(t(q))$ to $q^{1200}$.
*Formalization options:* (a) take the Gauss congruences for $(r_n)$ as the single axiom (each instance a finite integer statement); (b) formalize Krattenthaler--Rivoal's proof; (c) cite Beukers' eta-product identities.

**Consequences.** (4.3) $F\,t\in q\,\mathbb Z[[q]]$ and $F\in 1+q\mathbb Z[[q]]$ as series in $q$. (4.4) For every $m\ge1$, $q^m\in t^m\mathbb Z[[t]]=(4u)^m\mathbb Z[[4u]]$; in particular $[u^k]q^m\in4^k\mathbb Z$ and $=0$ for $k<m$.

## 5. The companion and its source

**Definition 5.1.** $B:=\sum_{n\ge1}b_nu^n\in\mathbb Q[[u]]$ with $b_n$ as in the Main Theorem.

**Lemma 5.2.** $L_1B=u$, and $B$ is the unique series with $L_1B=u$, $[u^0]B=0$, $[u^1]B=1$.
*Proof.* Coefficient extraction as in Corollary 2.5: the coefficient of $u^{n+1}$ ($n\ge1$) is the recurrence; the coefficient of $u^1$ is $P(0)b_1=1$; the coefficient of $u^0$ is $0$. Uniqueness: the recurrence determines $b_{n+1}$ from $b_n,b_{n-1}$. $\square$

**Lemma 5.3 (key identity).** For every $h\in\mathbb Q[[q]]$, regarded as a series in $u$ via $q=q(u)$,
$$L_1(f\,h)=\frac{\theta_q^2h}{f^3}.$$
*Proof.* Since $L_1f=0$: $L_1(fh)=Pf\,\theta_u^2h+(2P\theta_uf+Qf)\,\theta_uh$. Chain rule: $\theta_uh=\lambda\,\theta_qh$ with $\lambda=\theta_u\log q$, hence $\theta_u^2h=\lambda^2\theta_q^2h+(\theta_u\lambda)\theta_qh$. Therefore
$$L_1(fh)=Pf\lambda^2\,\theta_q^2h+\Big[Pf\,\theta_u\lambda+(2P\theta_uf+Qf)\lambda\Big]\theta_qh .$$
By (4.1), $\lambda=f^{-2}\sigma^{-1}$, so $\theta_u\lambda=-\lambda\big(2\theta_uf/f+\theta_u\sigma/\sigma\big)=-\lambda\big(2\theta_uf/f+Q/P\big)$ by (4.0); substituting, the bracket is $Pf\lambda\big(-2\theta_uf/f-Q/P\big)+2P\lambda\theta_uf+Q\lambda f=0$. Finally $Pf\lambda^2=\sigma^2f/(f^4\sigma^2)=f^{-3}$. $\square$

**Definition 5.4.** $\Psi:=f^3u=F\,f\,u\in q\,\mathbb Q[[q]]$, $\Psi=\sum_{m\ge1}\psi(m)q^m$, and $\Theta:=\theta_q^{-2}\Psi:=\sum_{m\ge1}\frac{\psi(m)}{m^2}q^m$ (so $\theta_q^2\Theta=\Psi$).

**Theorem 5.5.** $B=f\cdot\Theta$; i.e. $b_n=[u^n](f\,\Theta)$ for all $n$.
*Proof.* By Lemma 5.3, $L_1(f\Theta)=\Psi/f^3=u$. Also $f\Theta=\psi(1)\,q+O(q^2)=u+O(u^2)$ since $\psi(1)=[q^1](f^3u)=1$. By the uniqueness in Lemma 5.2, $f\Theta=B$. $\square$
(Exact check: $[u^n](f\Theta)=b_n$ for all $n\le256$, `lattice/sqrt_apery_verify.gp`.)

## 6. Integrality of $d_n^2b_n$

**Lemma 6.1.** $4^m\psi(m)\in\mathbb Z$ for all $m\ge1$.
*Proof.* $\Psi=f\cdot\tfrac14(Ft)$ as series in $q$. By (4.3), $Ft\in q\mathbb Z[[q]]$ and $F\in1+q\mathbb Z[[q]]$; by Theorem 3.2 applied to $G=F(t(q))$, $4^jf_j\in\mathbb Z$ for $f_j:=[q^j]f$. Hence
$4^m\psi(m)=\tfrac14\sum_{j=0}^{m-1}\big(4^jf_j\big)\big(4^{m-j}[q^{m-j}](Ft)\big)$, and each term is an integer divisible by $4$ because $m-j\ge1$. $\square$

**Lemma 6.2.** For $1\le m\le n$, $e_{n,m}:=[u^n]\big(f\,q^m\big)$ is an integer divisible by $4^m$; and $e_{n,m}=0$ for $m>n$.
*Proof.* $e_{n,m}=\sum_{i=0}^{n}a_i\,[u^{n-i}]q^m$ with $a_i=[u^i]f\in\mathbb Z$ (Corollary 3.3), and by (4.4) $[u^{n-i}]q^m\in4^{n-i}\mathbb Z$, vanishing unless $n-i\ge m$; so every nonzero term is divisible by $4^m$. $\square$

**Theorem 6.3.** $d_n^2\,b_n\in\mathbb Z$ for all $n\ge1$.
*Proof.* By Theorem 5.5, $b_n=[u^n](f\Theta)=\sum_{m=1}^n\frac{\psi(m)}{m^2}e_{n,m}=\sum_{m=1}^n\frac1{m^2}\big(4^m\psi(m)\big)\frac{e_{n,m}}{4^m}$, each summand in $\frac1{m^2}\mathbb Z$ by Lemmas 6.1–6.2, and $m\mid d_n$ for $m\le n$. $\square$
(Verified exactly for $n\le1500$; $d_nb_n\notin\mathbb Z$ already at $n=2$: $b_2=107/2$.)

## 7. The linear forms: positivity, decay, irrationality

**Lemma 7.1 (Casoratian).** $a_nb_{n-1}-a_{n-1}b_n=-C_{n-1}^2$ for $n\ge1$.
*Proof.* Induction: if $W_n:=a_nb_{n-1}-a_{n-1}b_n$, the recurrence gives $(n+1)^2W_{n+1}=4(2n-1)^2W_n$; $W_1=-1$; and $C_n/C_{n-1}=\frac{2(2n-1)}{n+1}$ gives $C_n^2/C_{n-1}^2=\frac{4(2n-1)^2}{(n+1)^2}$. $\square$

**Lemma 7.2 (growth).** $a_n>0$ and $a_{n+1}\ge\rho\,a_n$ for all $n\ge n_0$, for explicit $(\rho,n_0)$ with $\rho>16e^{2}=118.23$; e.g. $\rho=120$, $n_0=20$ (numerically $a_{n+1}/a_n$ increases from $53.4$ at $n=1$ through $118.2$ near $n=13$ towards $\lambda_1=135.88$).
*Proof sketch (Lean-friendly).* Let $r_n=a_{n+1}/a_n$. The recurrence gives $r_n=\frac{136n^2+68n+10}{(n+1)^2}-\frac{4(2n-1)^2}{(n+1)^2}\cdot\frac1{r_{n-1}}$. Show by induction that $r_n\in[\rho,136]$ for $n\ge n_0$: if $r_{n-1}\ge\rho$ then $r_n\ge\frac{136n^2+68n+10}{(n+1)^2}-\frac{4(2n-1)^2}{(n+1)^2\rho}$, which is $\ge\rho$ for $n\ge n_0$ (a rational inequality in $n$, checkable for the base case and monotone in $n$). The base case $r_{n_0}\ge\rho$ is a finite computation. $\square$

**Proposition 7.3.** The sequence $b_n/a_n$ is increasing and bounded; let $\xi=\lim b_n/a_n$. Then for all $n$,
$$a_n\xi-b_n=a_n\sum_{m\ge n}\frac{C_m^2}{a_ma_{m+1}}>0,\qquad\text{and}\qquad 0<a_n\xi-b_n\le\frac{16^n}{a_n\,\rho}\cdot\frac{1}{1-16/\rho^2}\quad(n\ge n_0).$$
*Proof.* $\frac{b_{m+1}}{a_{m+1}}-\frac{b_m}{a_m}=\frac{C_m^2}{a_ma_{m+1}}>0$ by Lemma 7.1; the sum converges by Lemma 7.2 ($a_m\ge a_{n_0}\rho^{m-n_0}$ and $C_m^2\le16^m$, $16/\rho^2<1$). Telescoping gives the first display; the bound follows from $C_m^2\le16^m$, $a_m\ge a_n\rho^{m-n}$, $a_{m+1}\ge\rho a_m$. $\square$

**Lemma 7.4 (lcm growth; EXTERNAL but available in Lean).** $\log d_n=n+o(n)$ (prime number theorem), or the weaker explicit $d_n\le e^{1.03883\,n}$ for all $n$ (Rosser–Schoenfeld). *Lean:* PNT+ provides `rate_Dlcm` (used in the Catalan worthiness project). Note that Mathlib's $d_n\le4^n$ is **not** sufficient here: one needs $\log d_n<\tfrac12\log(\rho/16)\cdot\ldots$, precisely $2\log d_n+n\log16-n\log\rho\to-\infty$, i.e. $\limsup\frac{\log d_n}n<\tfrac12\log\frac{\rho}{16}=1.0074$ for $\rho=120$ (or $1.0696$ using $\lambda_1$ itself); any $\varepsilon<0.007$ in $\log d_n\le(1+\varepsilon)n$ works, so PNT (or Rosser–Schoenfeld with $\rho\ge 16e^{2\cdot1.03883}=127.7$, which holds for $n_0$ a little larger) suffices.

**Theorem 7.5 (Main Theorem).** $\xi\notin\mathbb Q$. Moreover $\mu(\xi)\le1+\frac{\log\lambda_1+2}{\log(1/\lambda_2)-2}=50.66$, $\lambda_{1,2}=4(17\pm12\sqrt2)$.
*Proof.* Put $q_n=d_n^2a_n\in\mathbb Z$, $p_n=d_n^2b_n\in\mathbb Z$ (Theorem 6.3). Then $q_n\xi-p_n=d_n^2(a_n\xi-b_n)>0$ (Prop. 7.3) and $q_n\xi-p_n\le d_n^2\cdot\frac{16^n}{a_n\rho}\cdot c\le c'\,\frac{(d_n^2\,16^n)}{\rho^{n}}\to0$ by Lemma 7.4 and $16e^{2(1+\varepsilon)}<\rho$. A real number $\xi$ with integers $q_n,p_n$, $q_n\xi-p_n\ne0$, $q_n\xi-p_n\to0$ is irrational: if $\xi=A/B$ then $|q_n\xi-p_n|\ge1/B$. The measure follows from the standard lemma with $\log|q_n|\sim(\log\lambda_1+2)n$ and $\log|q_n\xi-p_n|\sim-(\log(1/\lambda_2)-2)n$ (Poincaré's theorem for the exact rates; not needed for irrationality). $\square$

## 8. Numerical certificates (for sanity checks in the formalization)

- $a_0,\dots,a_{10}$: $1,10,534,40900,3672550,360764460,37546578300,4068324634440,454111557020550,51858065465743900,6029992129304117524$.
- $b_1,\dots,b_8$: $1,\ \tfrac{107}2,\ \tfrac{12293}3,\ \tfrac{4415321}{12},\ \tfrac{1807203551}{50},\ \tfrac{22570175443}6,\ \tfrac{1497911874023249}{3675},\ \tfrac{17834541196873427}{392}$.
- $d_n^2b_n$ for $n=1..8$: $1,214,147516,52983852,130118655672,13542105265800,71899769953115952,32102174154372168600$.
- $a_nb_{n-1}-a_{n-1}b_n=-C_{n-1}^2$: $-1,-1,-4,-25,-196,-1764$ for $n=1..6$.
- $\xi=0.10018744922933940616775868213306112163\ldots$; $b_{60}/a_{60}$ agrees with $\xi$ to $>160$ digits.
- $F\sqrt{1-34t+t^2}=1-12t-156t^2-2964t^3-66300t^4-\cdots$ (integral; Remark 4.3).

## 9. What is *not* part of the blueprint
The identification $\xi=L(\Psi,2)$ with $\Psi=\tfrac14(\eta_1\eta_6)^{9/2}(\eta_2\eta_3)^{-3/2}$ (a non-congruence weight-three form), the symmetric-square/Beukers modular interpretation, and the exact asymptotics $|a_n\xi-b_n|^{1/n}\to\lambda_2$. The irrationality proof needs only §§1–7.

## 10. Dependency graph
```
Thm 1.2 (Apéry recurrence) ──▶ Cor 1.3 (L_Ap F = 0)
                                      │
Def 2.1 (f = √F) + Lemma 2.3 ────────▶ Prop 2.4 (L_1 f = 0) ──▶ Cor 2.5 (recurrence for a_n)
Lemma 3.1 ──▶ Thm 3.2 ──▶ Cor 3.3 (a_n ∈ ℤ)
Thm 4.2 (nome integral; EXTERNAL) ──▶ (4.3),(4.4)
Def 4.1 (nome) + (4.0) ──▶ Lemma 5.3 (L_1(fh) = θ_q²h/f³) ──▶ Thm 5.5 (B = f·Θ)
Thm 3.2 + (4.3) ──▶ Lemma 6.1 ;  Cor 3.3 + (4.4) ──▶ Lemma 6.2 ;  + Thm 5.5 ──▶ Thm 6.3 (d_n² b_n ∈ ℤ)
Lemma 7.1 (Casoratian) + Lemma 7.2 (growth) ──▶ Prop 7.3 (positivity, bound)
Prop 7.3 + Thm 6.3 + Lemma 7.4 (PNT+) ──▶ Thm 7.5 (ξ ∉ ℚ)
```
