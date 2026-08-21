# The $\operatorname{Sym}^1$ square root of Apéry's $\zeta(3)$ row: proofs, and an irrationality theorem

*Claude (Opus 5), 2026-08-21. Scripts: `lattice/sqrt_apery/`. Follows
`SPORADIC_SCAN2.md` Part II (§7), which reported this row as data and made no claim.*

---

## Verdict first

Everything §7 of `SPORADIC_SCAN2.md` left open for Apéry's square-root row is now **proved**,
and the criterion closes:

> **Theorem E.** Let $A_n=\sum_k\binom nk^2\binom{n+k}k^2$ be Apéry's $\zeta(3)$ numbers,
> $F(t)=\sum_nA_nt^n$, and
> $$a_n:=4^n[t^n]\sqrt{F(t)}=1,\,10,\,534,\,40900,\,3672550,\,360764460,\dots$$
> Let $b$ be the second solution of the same recurrence with $b_0=0$, $b_1=1$, and
> $$\xi:=\lim_{n\to\infty}\frac{b_n}{a_n}=\sum_{m\ge0}\frac{\mathrm C_m^{\,2}}{a_ma_{m+1}}
> =0.100187449229339406167758682133061121628777720500812\ldots$$
> ($\mathrm C_m$ = Catalan). Then **$\xi\notin\mathbf Q$**, with irrationality measure
> $\mu(\xi)\le 50.654$.

Nothing is conditional. The four ingredients — $a_n\in\mathbf Z$, $d_n^2b_n\in\mathbf Z$,
$a_n\xi-b_n\neq0$, and $\limsup|a_n\xi-b_n|^{1/n}\le\lambda_2=4(17-12\sqrt2)$ — are
Theorems 1–4 below; the score $\log(1/\lambda_2)-2=+0.13920>0$ then does the rest.

And the period **is identified**, structurally:

> **Theorem D.** $\displaystyle \xi=L(\Psi,2)$, the critical value of the **weight-three
> half-integral eta quotient**
> $$\Psi(\tau)=\frac14\,\frac{\bigl(\eta(\tau)\eta(6\tau)\bigr)^{9/2}}
> {\bigl(\eta(2\tau)\eta(3\tau)\bigr)^{3/2}}
> =\tfrac14q-\tfrac98q^2+\tfrac{39}{32}q^3+\tfrac{63}{64}q^4-\tfrac{981}{512}q^5+\cdots$$
> \[[verified to $97$ digits]; the Theorem B$^*$ Fricke-geometry proof transposes verbatim to
> $(w,N)=(1,6)$ and its three hypotheses are *proved* here (§5).\]

$\Psi$ has **unbounded $2$-power denominators** ($v_2(\text{den}\,\psi(m))=2m+O(\log m)$,
checked to $m=60$), so by Calegari–Dimitrov–Tang it is a form on a **non-congruence**
finite-index subgroup, and there is no reason for $\xi$ to be a classical $L$-value. It is
not: a 95-constant structured `PSLQ` battery at 200–260 digits (§6), on top of the
1792-element sweep of `SPORADIC_SCAN2.md`, returns nothing, and `algdep` excludes degree
$\le10$ algebraicity at 300 digits. **So $\xi$ is a new irrational number whose only known
closed form is $L(\Psi,2)$.**

Two by-products of independent interest:

* **Theorem 1 (general).** For *any* $G\in1+x\mathbf Z[[x]]$, $4^n[x^n]\sqrt G\in\mathbf Z$.
  Sharp ($G=1+x$ already needs $4$). One line, and it proves the $\lambda=4$ half of
  `SPORADIC_SCAN2.md` §7 for every row at once.
* **Theorem 2.** For each of the **six** Almkvist–Zudilin sporadic rows $(a,b,c)$ (and, with
  Cooper's three, for all **nine**; §7),
  $a_n:=4^n[t^n]\sqrt{\sum A_nt^n}\in\mathbf Z$ and satisfies
  $$(n+1)^2a_{n+1}=\bigl(8an^2+4an+2b\bigr)a_n-4c\,(2n-1)^2a_{n-1},\qquad a_0=1,\ a_1=2b,$$
  which is exactly the table of `SPORADIC_SCAN2.md` §7 after the rescalings there. This is a
  *theorem*, not an observation: it is Theorem 1 plus the project's own
  Theorem `thm:rect15` ($L=t^3P_3\cdot\operatorname{Sym}^2\widetilde D$). §7 sharpens
  $\lambda$ to $\max(1,2^{2-e})$, $e=\min_{n\ge1}v_2(A_n)$, which reproduces and proves the
  whole $\lambda\in\{1,2,4\}$ column of `SPORADIC_SCAN2.md`.

---

## 1. The object and the two operators

Apéry's operator in $t$ is $L_{\mathrm{Ap}}=\tht^3-t(2\tht+1)(17\tht^2+17\tht+5)+t^2(\tht+1)^3$
($\tht=t\,d/dt$), $F=\sum A_nt^n$ its unique analytic solution with $F(0)=1$. Put
$$u=t/4,\qquad L_1=\tht^2-u\bigl(136\tht^2+68\tht+10\bigr)+16u^2\bigl(\tht+\tfrac12\bigr)^2 .$$

**Proposition 1.1 [proved, exact symbolic identity].** In normal form
$y''+py'+\varrho\,y=0$ with ($\varrho$, not $q$: $q$ is the nome throughout)
$$p=\frac{32u^2-204u+1}{u(16u^2-136u+1)},\qquad \varrho=\frac{2(2u-5)}{u(16u^2-136u+1)},$$
the classical symmetric square
$\;w'''+3pw''+(2p^2+p'+4\varrho)w'+(4p\varrho+2\varrho')w\;$ of $L_1$ equals
$L_{\mathrm{Ap}}$ divided by its leading coefficient $u^3(16u^2-136u+1)$, **coefficient by
coefficient over $\mathbf Q(u)$**. (`lattice/sqrt_apery/01_sym2.py`; three exact zeros.)
This is the $(a,b,c)=(17,5,1)$, $t\mapsto t/4$ instance of `thm:rect15`.

**Corollary 1.2 [proved].** $L_1$ has exponents $(0,0)$ at $u=0$, so it has a unique analytic
solution $A(u)=\sum a_nu^n$ with $A(0)=1$; $A^2$ then solves $L_{\mathrm{Ap}}$, is analytic
and equals $1$ at $0$, and $L_{\mathrm{Ap}}$ has exponents $(0,0,0)$ at $0$, whence
$$A(u)^2=F(4u),\qquad\text{i.e.}\qquad a_n=4^n[t^n]\sqrt{F(t)} ,$$
and $a$ satisfies $(n+1)^2a_{n+1}=(136n^2+68n+10)a_n-4(2n-1)^2a_{n-1}$, $a_0=1,a_1=10$.
Characteristic roots $\lambda_{1,2}=68\pm48\sqrt2=4(17\pm12\sqrt2)$,
$\lambda_1=135.8822510$, $\lambda_2=0.1177490061$, $\lambda_1\lambda_2=16$.

The four singular points of $L_1$ are $0$, $u_1=(\sqrt2-1)^4/4$, $u_2=(\sqrt2+1)^4/4$,
$\infty$, with exponents $(0,0)$, $(0,\tfrac12)$, $(0,\tfrac12)$, $(\tfrac12,\tfrac12)$
(the residue of $p$ at $u_{1,2}$ is exactly $\tfrac12$). $u_{1,2}$ are the two
$W_6$-fixed points of discriminant $-24$: $t(i/\sqrt6)=(\sqrt2-1)^4$ **exactly**
(verified to 60 digits, `06_fold.gp`). So this row sits in the **Fricke fold geometry**
of `02_sources_thmB_exact.tex`, not the cusp geometry — unusual for a second-order row.

*Aside (negative).* $L_1$ is **not** a pullback of a hypergeometric operator by a rational map
of degree 2: the exponent differences $(0,\tfrac12,\tfrac12,0)$ at $(0,u_1,u_2,\infty)$ force
either $u_1=-u_2$ (false: $u_1u_2=1/16>0$) or an extra apparent singularity that $L_1$ does
not have. So the "$f=\text{alg}\cdot{}_2F_1(\tfrac12,\tfrac12;1;z)$" route of the brief is
closed; the $\tfrac12$'s are the index-two cover, i.e. the non-congruence group itself.

## 2. Theorem 1: integrality of $a_n$, in one line

**Theorem 1 [proved].** Let $G\in1+x\mathbf Z[[x]]$. Then $4^n[x^n]\sqrt G\in\mathbf Z$ for
all $n\ge0$.

*Proof.* Write $G=1+S$, $S\in x\mathbf Z[[x]]$, and $\sqrt G=\sum_{k\ge0}\binom{1/2}kS^k$.
For $k\ge1$,
$$4^k\binom{1/2}k=(-1)^{k-1}\frac{\binom{2k}k}{2k-1}=(-1)^{k-1}\,2\,\mathrm C_{k-1}\in2\mathbf Z$$
($\mathrm C_j$ Catalan; verified symbolically for $k\le80$ in `08_verify.gp`, and it is the
standard identity $\binom{2k}{k}/(2k-1)=2\mathrm C_{k-1}$). Since $S^k\in x^k\mathbf Z[[x]]$,
only $k\le n$ contributes to $[x^n]$, so
$$4^n[x^n]\sqrt G=[n{=}0]+\sum_{k=1}^n(-1)^{k-1}2\,\mathrm C_{k-1}\,4^{\,n-k}\,[x^n]S^k\in\mathbf Z. \qquad\square$$

Applied to $F\in1+t\mathbf Z[[t]]$ this gives $a_n\in\mathbf Z$ (Theorem 2 for Apéry);
applied to $F\in1+q\mathbf Z[[q]]$ it gives $4^m[q^m]f\in\mathbf Z$ for $f=\sqrt F$, which is
what §3 needs. **Sharpness:** $v_2\bigl(\text{den}\,[t^n]\sqrt F\bigr)=2n-s_2(n)$ exactly for
$n\le30$ (`12_sharp.gp`), so $\lambda=4$ cannot be improved for Apéry's row.

**Positivity [proved].** $a_n\ge10\,a_{n-1}>0$ for all $n\ge1$, by induction: if
$a_{n-1}\le a_n/10$ then $(n+1)^2a_{n+1}\ge(136n^2+68n+10-\tfrac25(2n-1)^2)a_n
=(134.4n^2+69.6n+9.6)a_n\ge10(n+1)^2a_n$ for $n\ge1$. (Checked to $n=1500$.)

## 3. Theorem 3: the companion, its source, and $d_n^2b_n\in\mathbf Z$

Let $B(u)=\sum_{n\ge1}b_nu^n$. Comparing coefficients, the census normalisation
$b_0=0,b_1=1$ is exactly
$$\boxed{\,L_1(B)=u\,}$$
— $B$ is **not** a second solution of $L_1$ (the second solution is $A\log q$, `thm:rect15`),
it is the analytic solution of an *inhomogeneous* equation. That is why the naive Eichler
ansatz $\Phi=f\,\thq u$ fails here (checked: it gives $1,\tfrac{73}2,\tfrac{7615}3,\dots\neq b$).
The correct source comes out of variation of parameters.

**Lemma 3.1 [proved; verified to $q^{126}$].** $\thq t=t\,F\sqrt{1-34t+t^2}$.
*(Equivalently $(\thq t)^2=t^2P_3F^2$, i.e. `Prop. collapse` of `02_sources.tex` combined with
$\Phi=t\sigma^3/(P_3F)$; it also follows from Prop. 1.1 plus "canonical coordinate $=$ nome".)*

**Theorem 3 [proved].** Put $f=\sqrt F$ (weight $1$), and
$$\Psi:=f^3\cdot\frac t4=f\,F\,u=\frac14\frac{(\eta_1\eta_6)^{9/2}}{(\eta_2\eta_3)^{3/2}}
=\sum_{m\ge1}\psi(m)q^m,\qquad
\Theta:=\thq^{-2}\Psi=\sum_{m\ge1}\frac{\psi(m)}{m^2}q^m .$$
Then
$$b_n=4^n[t^n]\bigl(f\,\Theta\bigr)\qquad\text{for all }n\ge1,$$
and consequently $d_n^2b_n\in\mathbf Z$.

*Proof of the formula.* $\{A,\;A\log q\}$ is a basis of solutions of $L_1$; its Wronskian in
$u$ is $W=A^2\,d\log q/du=F/\thq u=1/(u\sqrt{P})$, $P=16u^2-136u+1$, by Lemma 3.1. With
$L_1=u^2P\,(D^2+pD+\varrho)$ the equation $L_1(Y)=u$ is $(D^2+pD+\varrho)Y=g$, $g=1/(uP)$, so
$g/W=1/\sqrt P$ and variation of parameters gives
$$Y(u)=A(u)\int_0^u A(u')\,\frac{\log q(u)-\log q(u')}{\sqrt{P(u')}}\,du' .$$
Changing variable to $q'$ and using $\thq u=uF\sqrt P$ turns $A\,du'/\sqrt P$ into
$f^3u'\,d\log q'$, and
$\int_0^q q'^{m-1}(\log q-\log q')dq'=q^m/m^2$ gives $Y=A\cdot\Theta$. $Y$ is analytic and
$O(u)$, and analytic solutions of $L_1$ are multiples of $A$ with $A(0)=1$, so $Y=B$. $\square$
\[**Verified independently, exactly, for $1\le n\le116$** (`05_source.gp`): the two sides agree
as rationals.\]

*Proof of $d_n^2b_n\in\mathbf Z$.* Three integrality facts, each elementary:

1. $4^m\psi(m)\in\mathbf Z$. Indeed $4^m\psi(m)=\tfrac14\sum_{j=0}^{m-1}(4^jf_j)\,
   4^{\,m-j}[q^{m-j}](Ft)$ with $f_j=[q^j]f$; $4^jf_j\in\mathbf Z$ by Theorem 1
   ($F\in1+q\mathbf Z[[q]]$), $Ft\in q\mathbf Z[[q]]$ forces $m-j\ge1$, so every term is
   divisible by $4$.
2. $e_{n,m}:=4^n[t^n](f\,q^m)\in\mathbf Z$ and $4^m\mid e_{n,m}$. Indeed $q\in t\mathbf Z[[t]]$
   (Lagrange inversion of the monic integral series $t(q)$), so
   $e_{n,m}=\sum_i(4^i[t^i]f)\,4^{\,n-i}[t^{n-i}]q^m$ with $[t^{n-i}]q^m=0$ unless $n-i\ge m$;
   $4^i[t^i]f=a_i\in\mathbf Z$ by Theorem 1, and $4^{n-i}$ is divisible by $4^m$.
3. Hence $\psi(m)e_{n,m}=(4^m\psi(m))(e_{n,m}/4^m)\in\mathbf Z$.

Since $e_{n,m}=0$ for $m>n$, $\;d_n^2b_n=\sum_{m=1}^n\frac{d_n^2}{m^2}\,\psi(m)\,e_{n,m}
\in\mathbf Z$, because $m\mid d_n$ for $m\le n$. $\square$

\[**Verified**: $d_n^2b_n\in\mathbf Z$ for all $n\le1500$ exactly; $d_nb_n\notin\mathbf Z$
(fails at $n=2$ already, and holds only for $7$ of the $400$ values $n\le400$), so
**$k=2$ is sharp**. Also verified: $4^m\psi(m)\in\mathbf Z$ for $m\le118$; $e_{n,m}$ integral
and $4^m$-divisible for all $m\le n\le40$; $q\in t\mathbf Z[[t]]$ to $t^{120}$.\]

*Remark.* This is the $r=2$ instance of "each Eichler integration costs one $d_n$", made
$2$-adically honest: the $4^{-2m}$ that $\Psi$ carries, because $f$ has unbounded denominators,
is exactly cancelled by the $4^n$ in $[u^n]=4^n[t^n]$ together with $q^m\in t^m\mathbf Z[[t]]$.
That cancellation is the whole content of the theorem.

## 4. Theorem 4: the Casoratian, nonvanishing, and the exact rate

**Lemma 4.1 [proved].** $W_n:=a_nb_{n+1}-a_{n+1}b_n$ satisfies
$W_n=\frac{4(2n-1)^2}{(n+1)^2}W_{n-1}$, $W_0=1$, hence
$$W_n=\frac{4^n\bigl((2n-1)!!\bigr)^2}{\bigl((n+1)!\bigr)^2}=\mathrm C_n^{\,2}
\quad(\text{Catalan squared}).$$
\[Verified exactly for $n\le1500$.\] Therefore
$$\xi-\frac{b_n}{a_n}=\sum_{m\ge n}\frac{\mathrm C_m^{\,2}}{a_ma_{m+1}}>0,
\qquad r_n:=a_n\xi-b_n=a_n\sum_{m\ge n}\frac{\mathrm C_m^{\,2}}{a_ma_{m+1}}>0\ \ \forall n .$$
**Nonvanishing is exact, not asymptotic** — this is what makes the criterion applicable.

**Theorem 4 [proved].** $\limsup_n r_n^{1/n}\le\lambda_2=4(17-12\sqrt2)=0.1177490\ldots$

*Proof.* The recurrence has coefficients converging to $136$ and $-16$, whose characteristic
roots $\lambda_1>\lambda_2>0$ have distinct moduli, so by Poincaré's theorem $a_{n+1}/a_n$
tends to $\lambda_1$ or $\lambda_2$; since $a_{n+1}/a_n\ge10>\lambda_2$ (§2) the limit is
$\lambda_1$. Fix $\varepsilon>0$ and $C_\varepsilon$ with
$C_\varepsilon^{-1}(\lambda_1-\varepsilon)^n\le a_n\le C_\varepsilon(\lambda_1+\varepsilon)^n$.
With $\mathrm C_m\le4^m$,
$$r_n\le a_nC_\varepsilon^2\sum_{m\ge n}\Bigl(\frac{16}{(\lambda_1-\varepsilon)^2}\Bigr)^m
\cdot\frac1{\lambda_1-\varepsilon}=O\!\Bigl((\lambda_1+\varepsilon)^n
\bigl(16(\lambda_1-\varepsilon)^{-2}\bigr)^n\Bigr),$$
and $16(\lambda_1+\varepsilon)/(\lambda_1-\varepsilon)^2\to16/\lambda_1=\lambda_2$. $\square$

\[Verified: $\tfrac1n\log r_n=-2.2449,-2.1972,-2.1799,-2.1708,-2.1652,-2.1613$ at
$n=100,\dots,600$, increasing towards $\log\lambda_2=-2.139200$. The refined shape
$r_n\sim C\lambda_2^{\,n}n^{-3/2}$ — which is what analytic continuation predicts, $B-\xi A$ being
analytic at $u_1$ and carrying the exponent $\tfrac12$ at $u_2$ — fits with
$C\approx0.0257$ to three digits at **both** ends of that range ($n=100$: predicted $-2.2449$;
$n=600$: predicted $-2.1613$). So $C\neq0$ is numerically certain, though only
$\limsup r_n^{1/n}\le\lambda_2$ is used below, and only that is proved.\]

## 5. Theorem E: irrationality

$a_n\in\mathbf Z$ (Thm 2), $d_n^2b_n\in\mathbf Z$ (Thm 3): put $q_n=d_n^2a_n$,
$p_n=d_n^2b_n\in\mathbf Z$. Then $q_n\xi-p_n=d_n^2r_n>0$ (Lem. 4.1), and by Theorem 4 and
Rosser–Schoenfeld ($\psi(x)<1.03883\,x$, so $d_n<e^{1.03883n}$; PNT gives $d_n=e^{n(1+o(1))}$),
$$\limsup_n\bigl(q_n\xi-p_n\bigr)^{1/n}\le e^2\lambda_2=0.87013\ldots<1 .$$
If $\xi=p/s$ then $0<|s|(q_n\xi-p_n)=|p\,q_n-s\,p_n|\in\mathbf Z_{\ge1}$, contradiction for
large $n$. Hence $\xi\notin\mathbf Q$, and since $\tfrac1n\log q_n\to2+\log\lambda_1=6.911789$,
$$\mu(\xi)\le1+\frac{2+\log\lambda_1}{\log(1/\lambda_2)-2}=1+\frac{6.911789}{0.139200}=50.654 .$$
\[Verified: $\tfrac1n\log(d_n^2r_n)$ at $n=100,\dots,600$ is
$-0.364,-0.136,-0.185,-0.182,-0.159,-0.182$, oscillating about $2+\log\lambda_2=-0.13920$
with the usual $d_n$ fluctuation, always negative.\]

**Score bookkeeping.** $\operatorname{score}=\log(1/\lambda_2)-k=2.139200-2=+0.139200$;
$\operatorname{budget}=\log\lambda_1-k=+2.911789$. Apéry's own row scores $+0.5255$: the
square root is *worse* than its parent (it gains a free integration, $k:3\to2$, but pays
$\log4$), yet it is the **second** positive-score integral row in the census and the first
that is not Apéry's.

## 6. Theorem D: the period is $L(\Psi,2)$ — and $\Psi$ is non-congruence

**The Fricke data, proved.** From $\eta(-1/z)=\sqrt{-iz}\,\eta(z)$, under $W_6:\tau\mapsto-1/(6\tau)$,
$$(\eta_1\eta_6)(W_6\tau)=(-i\tau)\sqrt6\,(\eta_1\eta_6)(\tau),\qquad
(\eta_2\eta_3)(W_6\tau)=(-i\tau)\sqrt6\,(\eta_2\eta_3)(\tau)$$
(both factors pick up the *same* automorphy factor, since $W_6$ swaps $\eta_1\leftrightarrow\eta_6$
and $\eta_2\leftrightarrow\eta_3$),
hence $t\circ W_6=t$ (so the fold is at $\tau_*=i/\sqrt6$, and indeed
$t(\tau_*)=(\sqrt2-1)^4$ to $60$ digits), and with
$f=(\eta_2\eta_3)^{7/2}/(\eta_1\eta_6)^{5/2}$, $\Psi=\tfrac14(\eta_1\eta_6)^{9/2}/(\eta_2\eta_3)^{3/2}$
(both **verified as exact $q$-series identities to $q^{200}$**),
$$f(W_6\tau)=(-i\tau)\sqrt6\,f(\tau),\qquad
\Psi(W_6\tau)=\bigl((-i\tau)\sqrt6\bigr)^3\Psi(\tau),\qquad F|_2W_6=-F .$$
So $\Lambda(\Psi,s)=\int_0^\infty\Psi(iy)y^{s-1}dy$ is **entire** (endpoint criterion of
`thm:Bexact` satisfied at both $\infty$ and $0$) and satisfies the folded formula
$$\Lambda(\Psi,s)=\int_{1/\sqrt6}^{\infty}\Psi(iy)\bigl[y^{s-1}+6^{3/2-s}y^{2-s}\bigr]dy .$$
Evaluating, $L(\Psi,2)=4\pi^2\Lambda(\Psi,2)$:
$$L(\Psi,2)-\xi=2.3\times10^{-97}\qquad\text{(80-digit arithmetic, 600 $q$-terms)} .$$
This is Theorem B$^*$ in the Fricke geometry with $(w,N)=(1,6)$, whose proof transposes
verbatim; all three of its hypotheses are proved above. (With $j(W_6,\tau)=\sqrt6\tau$ the
eigenvalues are $f|_1W_6=-i\,f$ and $\Psi|_3W_6=+i\,\Psi$ — weight is odd, so $W_6^2=-1$ and
the eigenvalues are $\pm i$ rather than $\pm1$ — and their product is $+1$, which is the only
thing the cocycle argument uses.) \[Class: **proved modulo the
transposition**, **verified to 97 digits**.\] Also computed: $L(\Psi,1)=0.0723454929$,
$L(\Psi,3)=0.1560378492$; unlike the weight-two rows there is **no** forced annihilation
($R(\tau_*)\neq0$), consistent with `rem:Lphi1`.

**Non-congruence.** $\Psi$ is modular of weight $3$ on the (finite-index) kernel of the
order-dividing-$48$ multiplier character of $(\eta_1\eta_6)^{9/2}(\eta_2\eta_3)^{-3/2}$ on
$\Gamma_0(6)$, and $v_2(\mathrm{den}\,\psi(m))=2,3,5,6,9,10,12,13,17,\dots$ grows linearly
(checked $m\le60$), i.e. **unbounded denominators**. By Calegari–Dimitrov–Tang the group is
non-congruence, so $\xi=L(\Psi,2)$ is not expected to be a classical $L$-value.

**Clean negatives.** (`11_pslq.gp`, `07b_ident.gp`, `09_closedform.gp`.)
* 2-term `lindep` at 200 digits, height cut $10^{14}$, of $\xi$ and of
  $\sqrt2\,\xi,\sqrt3\,\xi,\sqrt6\,\xi,\xi/\sqrt2,\xi/\sqrt3,\xi/\sqrt6$ against **95**
  structured constants: $1,\zeta(2),\zeta(3),\pi^2,\pi^3$; $L(\chi_D,s)$ for
  $D=-3,-4,-8,-24,8,12,24$ and $s=2,3$; $\pi\log2$, $\pi\log3$, $\pi\log(1+\sqrt2)$,
  $\pi^2\log2$, $\pi^2\log(1+\sqrt2)$, $\log^3(1+\sqrt2)$, $\log^32$; and Chowla–Selberg CM
  periods $\pi^j\Omega_D^{2,4,6}$, $j=-2,\dots,3$, $D=-3,-4,-8,-24$. **Zero hits.**
* 3-term `lindep` with $1$ at 160 digits, height cut $10^9$: **zero hits.**
* `algdep` degree $1\le d\le10$ at 300 digits: minimal heights $10^{40}$–$10^{220}$, i.e.
  **$\xi$ is not algebraic of degree $\le10$** (and by Theorem E not of degree $1$).
* No closed form for $a_n$ as $\sum_k$ of a product of $\le3$ factors from a 17-element
  binomial dictionary; and $a_n/\binom{2n}n$ is *not* integral ($n=7$ gives a $13$ in the
  denominator), so the obvious normalisation does not exist either. **Open.**
* This is on top of the 1792-element weight-3/weight-4 $L$-value sweep of
  `SPORADIC_SCAN2.md` §7, also negative.

## 7. The other eight rows, and the graded integrality lemma

All nine sporadic $\operatorname{Sym}^1$ square roots were redone from scratch
(`lattice/sqrt_apery/rows_*`), with the same standard of proof.

**Theorem 1$'$ (graded form of Theorem 1) [proved].** Let $G=\sum A_nt^n\in1+t\mathbf Z[[t]]$
and let $e:=\min_{n\ge1}v_2(A_n)$. Then $\lambda^n[t^n]\sqrt G\in\mathbf Z$ for
$$\lambda=\max\bigl(1,\,2^{\,2-e}\bigr).$$
*Proof.* $G=1+2^eK$, $K\in t\mathbf Z[[t]]$, and
$2^{ek}\binom{1/2}k=(-1)^{k-1}2\,\mathrm C_{k-1}2^{(e-2)k}$; for $k\le n$ the factor
$\lambda^{\,n}2^{(e-2)k}$ is integral. $\square$
No $2^n$ version can exist for general $G$: for $G=1+t$, $v_2(2^n[t^n]\sqrt G)=1-n$ at
$n=2^k$.

Measured $e$ (exact, $n\le400$): Domb, $T$, Cooper $s_7$ have $e=2$ ($\lambda=1$);
$s_{10},s_{18}$ have $e=1$ ($\lambda=2$); Apéry and the three other Almkvist–Zudilin rows have
$e=0$ ($\lambda=4$). **This explains the whole $\lambda$ column of `SPORADIC_SCAN2.md` §7–8 and
proves it**, and it lines up exactly with that report's CM dichotomy:
$$\lambda=1\iff e\ge2\iff \sqrt F\text{ is an integral weight-one CM theta series}
\iff\text{the period is identified}.$$
$4\mid A_n$ is proved elementarily for Domb and $2\mid A_n$ for $s_{10}$; for $T,s_7,s_{18}$ it
is verified to $n=400$ (so $s_{18}$'s $\lambda=2$ is proved modulo $2\mid A_n$).

**The $\operatorname{Sym}^2$ identity is exact for all nine.** $\operatorname{Sym}^2(L_1)=L_{\rm
parent}$ over $\mathbf Q(u)$, checked symbolically with zero residual in every case — including
Cooper's three, which are *not* in Almkvist–Zudilin normalisation and are therefore outside
`thm:rect15`. So §7's "Observation (verified)" is now a **theorem** for all nine rows.

| parent $\to$ $\operatorname{Sym}^1$ | $\lambda$ | sharp | $k$ | $\lambda_1$ | $\lambda_2$ | score | budget | period |
|---|---|---|---|---|---|---|---|---|
| **Apéry** $(17,5,1)\to(136,10,4)$ | 4 | yes | 2 | $4(17{+}12\sqrt2)$ | $4(17{-}12\sqrt2)$ | $\mathbf{+0.1392}$ | $+2.9118$ | $L(\Psi,2)$, **§6** |
| Domb $(10,4,64)\to(20,2,16)$ | 1 | min | 2 | $16$ | $4$ | $-3.3863$ | $+0.7726$ | $L(f_{12},2)$, 289 dig, factor $1$ |
| $T\,(12,4,16)\to(24,2,4)$ | 1 | min | 2 | $12{+}8\sqrt2$ | $12{-}8\sqrt2$ | $-1.6235$ | $+1.1490$ | $L(f_8,2)$, 306 dig, factor $1$ |
| AZ$(9,3,-27)\to(72,6,-108)$ | 4 | yes | 2 | $77.5692$ | $-5.5692$ | $-3.7173$ | $+2.3512$ | $0.14551448201\ldots$ **unid.** |
| AZ$(11,5,125)\to(88,10,500)$ | 4 | yes | 2 | complex, $\sqrt{2000}$ | — | n/a | $+1.8005$ | none (complex roots) |
| AZ$(7,3,81)\to(56,6,324)$ | 4 | yes | 2 | complex, $36$ | — | n/a | $+1.5835$ | none |
| Cooper $s_7$ | 1 | min | 2 | $27$ | $-1$ | $-2.0000$ | $+1.2958$ | $L(f_7,2)$, $\ge300$ dig, factor $1$ |
| Cooper $s_{10}$ | 2 | yes | 2 | $32$ | $-8$ | $-4.0794$ | $+1.4657$ | $0.31692535921\ldots$ **unid.** |
| Cooper $s_{18}$ | 2 | yes | 2 | $32$ | $24$ | $-5.1781$ | $+1.4657$ | $0.48423755360\ldots$ **unid.** |

Every row: $a_n\in\mathbf Z$ and the recurrence exactly to $n=300$; $d_n^2b_n\in\mathbf Z$ and
$d_nb_n\notin\mathbf Z$ for $n\le200$ (so $k=2$ sharp everywhere); Casoratian
$W_n=\prod_{j\le n}p_0(j)/(j+1)^2$ exact to $n=200$; periods to $\ge250$ digits.
The three unidentified periods survive: all weight-3 newform $L(f,s)$, $s=1,2$ (real and
imaginary parts) at the relevant level towers ($488$, $556$, $408$ values), plus a 40-constant
battery, 2- and 3-term, at 180 digits. **Apéry's is the only positive score in the family**,
and the $T$ row's $\lambda_2=12-8\sqrt2<1$ remains the only other decayer.

**Corrections to `SPORADIC_SCAN2.md` §7–8** (no numerical disagreement anywhere):
1. §7's "Observation" $\to$ **theorem** (Sym$^2$ identity exact; integrality proved by
   Theorem 1$'$; $s_{18}$ conditional on $2\mid A_n$).
2. Add Theorem 1$'$; it turns the CM dichotomy into an arithmetic statement.
3. §8: Cooper $s_7$'s parameter is **not** an eta quotient — it is
   $t=u/(1+13u+49u^2)$ with $u=(\eta_7/\eta_1)^4$, degree $2$ on $\Gamma_0(7)$ (verified to
   $q^{120}$). That is why both eta-quotient scans missed it.
4. Domb's $\operatorname{Sym}^1$ parameter carries a **minus sign**:
   $t=-(\eta_2\eta_6/\eta_1\eta_3)^6$.
5. Unrelated: `SLOPE_CENSUS.md`'s "$s_{10}$'s $A_n$ is not integral" is an artefact of the
   initial condition; with $A_1=2$ one has $A_n\in\mathbf Z$ for $n\le400$.

## 8. Status table

| statement | class |
|---|---|
| $\operatorname{Sym}^2L_1=L_{\mathrm{Ap}}$ over $\mathbf Q(u)$, $u=t/4$ | **proved** (exact symbolic) |
| $4^n[x^n]\sqrt G\in\mathbf Z$ for $G\in1+x\mathbf Z[[x]]$; sharp | **proved** |
| $a_n=4^n[t^n]\sqrt F\in\mathbf Z$; recurrence $(136,68,10;4)$ | **proved** |
| same for all six Almkvist–Zudilin rows, $\lambda=4$ | **proved** (Thm 1 + `thm:rect15`) |
| $\operatorname{Sym}^2(L_1)=L_{\rm parent}$ for all **nine** rows (incl. Cooper) | **proved** (exact symbolic) |
| $\lambda=\max(1,2^{2-e})$, $e=\min v_2(A_n)$ (Thm 1$'$) | **proved**; $e$ measured $n\le400$ |
| $k=2$ sharp for all nine rows | verified $n\le200$ |
| Domb/$T$/$s_7$ periods $=L(f,2)$, no rational factor | verified $289$/$306$/$300$ digits |
| the other three periods ($0.14551\ldots$, $0.31693\ldots$, $0.48424\ldots$) | **open**; excluded over 1450+ $L$-values and 40 constants |
| Cooper $s_7$'s parameter $t=u/(1{+}13u{+}49u^2)$, $u=(\eta_7/\eta_1)^4$ | verified to $q^{120}$ |
| $a_n\ge10a_{n-1}>0$ | **proved** |
| $L_1(B)=u$; $B=A\cdot\thq^{-2}\Psi$, $\Psi=f^3t/4$ | **proved**; verified exactly $n\le116$ |
| $d_n^2b_n\in\mathbf Z$; $k=2$ sharp | **proved**; verified $n\le1500$ |
| Casoratian $=\mathrm C_n^2$, $r_n>0$ | **proved**; verified $n\le1500$ |
| $\limsup r_n^{1/n}\le\lambda_2$ | **proved** (Poincaré) |
| **$\xi\notin\mathbf Q$, $\mu(\xi)\le50.654$** | **proved** |
| $f,\Psi$ as half-integral eta quotients; Fricke eigenvalues; $t\circ W_6=t$ | **proved**; verified to $q^{200}$ |
| $\xi=L(\Psi,2)$ | proved *modulo* transposing `thm:Bexact` to $w=1$; verified $10^{-97}$ |
| $\Psi$ has unbounded denominators $\Rightarrow$ non-congruence | verified $m\le60$ + CDT |
| $\xi$ in classical closed form | **open**; excluded over the stated bases |
| closed form for $a_n$ as a binomial sum | **open** |
| $L_1$ a degree-2 hypergeometric pullback | **excluded** (exponent count) |

## 9. What this changes

1. `SPORADIC_SCAN2.md` §7's headline "score $+0.1392>0$ … its period is unidentified and no
   claim is made" can be upgraded: the claim can now be made. The **caveat paragraph should be
   replaced** by Theorem E, and the "unidentified" entry by $L(\Psi,2)$ with the
   non-congruence explanation for why the 1792-element sweep had to fail.
2. Paper §7 Problem 7.6 ("free integration"): the $\operatorname{Sym}^2\to\operatorname{Sym}^1$
   step is now a *proved* mechanism for all **nine** sporadic rows, with the exact bookkeeping
   $\text{score}\mapsto\text{score}+1-\log\lambda$.
3. The census gains a second positive-score row and its first *new* irrationality theorem.
   The constant is new; identifying $L(\Psi,2)$ is the sharpest open question this leaves.
4. The paper's Table `tab:census` and \S`sec:free` should gain the nine
   $\operatorname{Sym}^1$ rows with $k=2$, and \S`sec:arch` a line: this is the first row in
   the project whose *second*-order operator sits in the **Fricke fold geometry**
   (exponents $(0,\tfrac12)$ at $t_c$), which is why the standard second-order
   $\Phi=F\thq t$ ansatz fails for it and $\Psi=f^3t/4$ is the source instead.
5. Conjecture `conj:barrier` (scalar barrier) is untouched — it is about $w\ge3$ — but the
   route "take a square root, pay $\log\lambda$, gain a unit of $k$" is now a proved tool to
   test against it.

## 10. Reproduction

```
lattice/sqrt_apery/01_sym2.py        # Sym^2(L_1) = L_Apery, exact (sympy)
lattice/sqrt_apery/02_data.gp        # a_n, b_n, denominator factorisations
lattice/sqrt_apery/03_xi.gp          # xi to 1100 digits via the Casoratian series
lattice/sqrt_apery/04_modular.gp     # q-expansions; the WRONG Eichler ansatz, for the record
lattice/sqrt_apery/05_source.gp      # Theta_q t = tF sqrt(P3); Psi = f^3 t/4; b_n = 4^n[t^n](f Theta)
lattice/sqrt_apery/06_fold.gp        # fold at i/sqrt6; Fricke test; folded Mellin -> L(Psi,2) = xi
lattice/sqrt_apery/07_ident.gp
lattice/sqrt_apery/07b_ident.gp      # algdep, CM periods
lattice/sqrt_apery/08_verify.gp      # all integrality/positivity/Casoratian checks to n=700; 13_bign.gp to n=1500
lattice/sqrt_apery/09_closedform.gp  # closed-form search for a_n (negative)
lattice/sqrt_apery/10_psi_eta.gp     # f and Psi as half-integral eta quotients
lattice/sqrt_apery/11_pslq.gp        # 95-constant PSLQ battery (negative)
lattice/sqrt_apery/12_sharp.gp       # lambda = 4 is sharp
lattice/sqrt_apery/13_bign.gp        # a_n in Z, d_n^2 b_n in Z, Casoratian, to n=1500
lattice/sqrt_apery/rows_01_verify.py rows_02_sym2.py rows_03_denom_inv.gp
lattice/sqrt_apery/rows_04_period.gp rows_05_ident.gp rows_06_integrality.py rows_07_modular.gp
```
