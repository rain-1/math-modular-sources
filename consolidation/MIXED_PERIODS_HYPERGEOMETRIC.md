# Mixed dilogarithm periods on the hosts $1/\sqrt{1\mp 4mx}$: a family of linear independence theorems from the CDT holonomy bound

*Fable, 2026-09-01. Companion computations: session scratchpad `mixed/` (REPORT.md by an Opus agent, `dens*.gp`, `periods.py`). Literature: Calegari–Dimitrov–Tang, arXiv:2408.15403 v2 (`papers/cdt/cdt2/L2chi.tex`), Theorem "three elements" (= Theorem `noextrassmall`), Lemma "around 1/4", Theorem "mixed".*

*Tags: **[proved]** = complete argument in this note (modulo the cited CDT theorem); **[verified]** = exact computation over the stated range; **[numerical]** = high-precision numerics with PSLQ; **[open]**.*

---

## 0. Summary

CDT's Theorem `noextrassmall` says: if $f\in\mathbf Q[[x]]$ has denominator type $[1..n][1..n/2]$, is holomorphic on $\mathbf C\setminus[1,\infty)$, and continues meromorphically along all paths in $\mathbf P^1\setminus\{0,\delta,1,\infty\}$ for some $\delta\in(-\infty,1)$, then $f$ has **no monodromy at $\delta$ on any sheet** (part $(\ast)$; the full form says $f\in\mathbf Q[x,\tfrac1{1-x}]\langle 1,\log(1-x),\log^2(1-x)\rangle$). CDT use it once, with $\delta=\tfrac14$ on the host $H_A=1/\sqrt{1-4x}$, to get their "mixed" theorem for $L(2,\chi_{-3})-L(1,\chi_{-3})\log3$.

The point of this note is that **the fourth puncture $\delta$ is arbitrary**, so the same theorem applies to every host
$$H_A^{(m)}(x)=\frac1{\sqrt{1-4mx}}=\sum_{n\ge0}\binom{2n}{n}m^nx^n\in\mathbf Z[[x]],\qquad \delta=\frac1{4m},$$
and to $1/\sqrt{1+4mx}$ ($\delta=-\tfrac1{4m}$), always with the puncture at $1$. Two ingredients are needed and both are elementary:

* **Half-integration lemma [proved, §2].** If $h\in\mathbf Z[[t]]$ and $k=\sum k_jt^j$ has $[1..j]\,k_j\in\mathbf Z$, then $\int_0^x hk\,dt$ and $h\int_0^x hk\,dt$ have denominator type $[1..n][1..\lfloor n/2\rfloor]$ **exactly** (no constant). This is why $H_D=H_A\int H_A\log(1-t)/(1-t)$ has the $\tau=\tfrac32$ type on every host, and it also shows that CDT's $H_C=H_A\int H_A\log(1-t)/t$ does **not** (it has type $[1..n]^2$ with unbounded excess: its coefficient of $x^2$ is $-13/4$) — see §6, an erratum to CDT's Lemma "around 1/4".
* **Fold lemma [proved, §3].** With $u=\sqrt{1-4mx}$, every $H[k]:=H_A\int_0^xH_Ak\,dt$ decomposes as $c[k]\,H_A+E_k(x)$ with $E_k$ analytic on the disc $|x-\delta|<1-\delta$ (which contains $0$ and reaches the puncture $1$), and $c[k]=\int_0^\delta H_Ak\,dt$; so the fold monodromy is $T_\delta H[k]=H[k]-2c[k]H_A$, and a rational relation among $1$ and the $c[k]$ produces a fold-regular function of type $[1..n][1..n/2]$ on $\mathbf P^1\setminus\{0,\delta,1,\infty\}$, which then has to have no monodromy at $\delta$ on the sheets reached by looping around $1$; a two-sheet computation (§3.3) shows it does, unless the relation is trivial.

**Theorem A [proved modulo CDT $(\ast)$].** For every integer $m\ge1$ put $D'=4m-1$ and $\theta_m=2\arctan(1/\sqrt{D'})$, so that $\pi-\theta_m=2\arctan\sqrt{4m-1}$. Then
$$1,\qquad \frac{\theta_m}{\sqrt{D'}},\qquad \frac{\theta_m\log\frac{4m-1}{m}-2\,\mathrm{Cl}_2(\pi-\theta_m)}{\sqrt{D'}}$$
are linearly independent over $\mathbf Q$. Equivalently
$$2\,\mathrm{Cl}_2\!\bigl(2\arctan\sqrt{4m-1}\bigr)-2\arctan\!\tfrac1{\sqrt{4m-1}}\cdot\log\tfrac{4m-1}{m}\ \notin\ \mathbf Q\sqrt{4m-1}+\mathbf Q\,\arctan\tfrac1{\sqrt{4m-1}}.$$
$m=1$ is CDT's Theorem "mixed" (three-period version): $\theta_1=\pi/3$, $2\mathrm{Cl}_2(2\pi/3)=\sqrt3L(2,\chi_{-3})$. For $m\ge2$ the angle $\theta_m$ is not a rational multiple of $\pi$; $\mathrm{Cl}_2(\pi-\theta_m)=D(-e^{-i\theta_m})$ is the Bloch–Wigner dilogarithm (volume of an ideal hyperbolic tetrahedron) at the point $-\bar w_m$, $w_m=\frac{(4m-2)+2\sqrt{-(4m-1)}}{4m}\in\mathbf Q(\sqrt{-(4m-1)})$, $|w_m|=1$. The discriminants $-(4m-1)$ run through $-3,-7,-11,-15,-19,-23,\dots$; $m=2,3,5,11,17,41$ give the Heegner fields.

**Theorem B [proved modulo CDT $(\ast)$].** For every integer $m\ge1$ put $c=\sqrt{4m+1}$ and $\varepsilon=\frac{c+1}{c-1}$. Then
$$1,\qquad \frac{\log\varepsilon}{c},\qquad \frac{R_m}{c},\qquad R_m:=2\,\mathrm{Li}_2\!\Bigl(\frac{c-1}{2c}\Bigr)-\frac{\pi^2}6+\log(2c)\log\varepsilon+\log\frac{c+1}{2c}\log\frac{c-1}{2c}-\tfrac12\log^2(c+1)+\tfrac12\log^2(c-1)$$
are linearly independent over $\mathbf Q$. When $4m+1=c^2$ is a square ($c=3,5,7,\dots$, $m=2,6,12,\dots$) the dilogarithm sits at the rational point $\frac{c-1}{2c}=\frac13,\frac25,\frac37,\dots$; in particular ($c=3$)
$$2\,\mathrm{Li}_2\!\bigl(\tfrac13\bigr)-\frac{\pi^2}6+\log^23-\tfrac12\log^22\ \notin\ \mathbf Q+\mathbf Q\log2 .$$
(The irrationality of $\mathrm{Li}_2(1/3)$ itself is open; CDT list $n\in\{-4,-3,-2,2,3,4,5\}$ as the values for which $\mathrm{Li}_2(1/n)\notin\mathbf Q$ is unknown.)

**What this is and is not.** These are "mixed" statements of exactly CDT's showcase type, obtained at the cheapest level ($\tau=\tfrac32$, five functions, a single bivalent map) — no new holonomy bound is used. They are not pure $L$-value results, and §5 explains structurally why the $\log\frac{4m-1}m$ mixing cannot be removed inside the $\tau=\tfrac32$ layer: every period of that layer lies in $\mathbf Q+\mathbf Q\,\theta_m/\sqrt{D'}+\mathbf Q\sqrt{D'}\theta_m+\mathbf Q\,Q_m/\sqrt{D'}$ with one fixed combination $Q_m=2\mathrm{Cl}_2(\pi-\theta_m)-\theta_m\log\frac{D'}m$.

---

## 1. The CDT input

Throughout $[1..n]=\mathrm{lcm}(1,\dots,n)$ and $[1..n/2]:=[1..\lfloor n/2\rfloor]$. A series $\sum a_nx^n\in\mathbf Q[[x]]$ *has type* $T(n)$ if $T(n)a_n\in\mathbf Z$ for all $n\ge0$ ($T(0)=1$).

> **Theorem $(\ast)$** (CDT, arXiv:2408.15403, Theorem `noextrassmall`, part $(\ast)$; proved in their §"bivalent app" from Theorem `basic main` with the bivalent map $\varphi(z)=8(z+z^3)/(1+z)^4$, $\varphi'(0)=8$, Bost–Charles integral $\log8+4G/\pi$, and the five functions $1,\log(1-x),\log^2(1-x),f,f(x/(x-1))$ with $\tau=69/50$; the bound reads $5\le4.640$.)
> Let $f\in\mathbf Q[[x]]$ have type $[1..n][1..n/2]$, be holomorphic in $\mathbf C\setminus[1,\infty)$, and continue as a meromorphic function along all paths in $\mathbf P^1\setminus\{0,\delta,1,\infty\}$ for some $\delta\in(-\infty,1)$. Then $f$ continues meromorphically along all paths in $\mathbf P^1\setminus\{0,1,\infty\}$.

Nothing in the proof uses the value of $\delta$ (only $\delta\ne\delta/(\delta-1)$, i.e. $\delta\ne2$, which is automatic). Here is the proof, from the source (`L2chi.tex`, §"bivalent app"), reduced to its inputs, so that the only remaining citation is the holonomy bound itself.

**Input 1 — the basic holonomy bound** (CDT Theorem `basic main`). Let $\mathbf b=(b_{i,j})$ be an $m\times r$ array of nonnegative reals whose columns have the form $0=b_{1,j}=\dots=b_{u_j,j}<b_{u_j+1,j}=\dots=b_{m,j}=:b_j$, with row sums $\sigma_i$, and put $\tau(\mathbf b)=\frac1{m^2}\sum_i(2i-1)\sigma_i$. Let $\varphi\colon(\mathbf D,0)\to(\mathbf C,0)$ be holomorphic with $|\varphi'(0)|>e^{\sigma_m}$. If $f_1,\dots,f_m\in\mathbf Q[[x]]$ are $\mathbf Q(x)$-linearly independent, $f_i=\sum_na_{i,n}x^n/\prod_j[1..b_{i,j}n]$ with $a_{i,n}\in\mathbf Z$, and every $f_i\circ\varphi$ is meromorphic on $\mathbf D$, then
$$m\le\frac{\iint_{\mathbf T^2}\log|\varphi(z)-\varphi(w)|\,d\mu(z)d\mu(w)}{\log|\varphi'(0)|-\tau(\mathbf b)}.$$
(The repository's `CDT_UNPACKED.md` §1 and `ADELIC_HOLONOMY.md` §2 rederive this bound: Siegel's lemma with $mD$ unknowns and $(1-\varepsilon)mD$ vanishing conditions, Liouville with the rearranged denominator rate $\tau$, and the Bost–Charles capacity estimate on the circle.)

**Input 2 — the bivalent map** (CDT Basic Remark `bivalent clean`, Lemma `bivalent BC integral`). $\varphi(z)=8(z+z^3)/(1+z)^4=1-\bigl(\frac{1-z}{1+z}\bigr)^4$ maps $\mathbf D$ onto $\mathbf C\setminus\{1\}$, is a bijection $(-1,1)\to(-\infty,1)$, maps each open half-disc conformally onto $\mathbf C\setminus(-\infty,1]$, has $\varphi'(z)=8(1-z)^3/(1+z)^5\ne0$ on $\mathbf D$, $\varphi'(0)=8$, and its Bost–Charles integral is $\log8+4G/\pi=3.2452\ldots$ (Smyth's Mahler measure $m(1+x+y-xy)=2G/\pi$ appears as the overflow term).

**Input 3 — overconvergent pullback** (CDT Proposition `overconvergence`). If $g$ is a germ at $0$ continuing as a multivalued holomorphic function on $\mathbf P^1\setminus\Sigma$, $\Sigma=\Sigma^0\sqcup\Sigma^1$, $\varphi\colon\mathbf D\to\mathbf P^1\setminus\Sigma^1$ with $\varphi(0)=0$ restricts univalently to a contractible $\Omega\ni0$ with $\varphi^{-1}(\Sigma^0)\subset\Omega$, and $g$ is holomorphic on $\varphi(\Omega)$, then $g\circ\varphi$ converges on all of $\mathbf D$. For the bivalent map and $\Sigma^0\subset(-\infty,1)$ one takes $\Omega$ a thin neighbourhood of the diameter $(-1,1)$: real values $w<1$ have exactly one preimage in $\mathbf D$ (on the diameter), so $\varphi$ is injective near the diameter and $\varphi^{-1}(\varphi(\Omega))=\Omega$ for $\Omega$ thin enough.

**Proof of Theorem $(\ast)$.** Suppose $f$ as in the statement has nontrivial monodromy at $\delta$ on some sheet. The involution $\iota(x)=x/(x-1)$ fixes $0$, exchanges $1\leftrightarrow\infty$, maps $(-\infty,1)$ to itself and $\delta$ to $\delta'=\delta/(\delta-1)\ne\delta$, and preserves the type $[1..n][1..n/2]$ (it is $x\mapsto -x-x^2-x^3-\dots$, an integral substitution). So $f\circ\iota\in\mathbf Q[[x]]$ has the same properties with $\delta'$ in place of $\delta$. The five functions $1,\log(1-x),\log^2(1-x),f,f\circ\iota$ are $\mathbf C(x)$-linearly independent (the first three have no monodromy away from $\{0,1,\infty\}$, $f$ has monodromy at $\delta$ but not at $\delta'$, $f\circ\iota$ the reverse), of types $1,[1..n],[1..n][1..n/2],[1..n][1..n/2],[1..n][1..n/2]$; the array has two columns $(0,1,1,1,1)$, $(0,0,\tfrac12,\tfrac12,\tfrac12)$ and $\tau=\frac{1\cdot0+3\cdot1+(5+7+9)\cdot\frac32}{25}=\frac{69}{50}$. Apply Input 3 with $\Sigma^1=\{1,\infty\}$, $\Sigma^0=\{0,\delta,\delta'\}\subset(-\infty,1)$: after multiplying by a polynomial $Q$ clearing poles, all five pullbacks are holomorphic on $\mathbf D$ (the hypothesis "holomorphic on $\mathbf C\setminus[1,\infty)$" is exactly holomorphy on $\varphi(\Omega)$). Input 1 gives $5\le\frac{\log8+4G/\pi}{\log8-69/50}=4.6404$, a contradiction. $\square$

**Two free generalisations [proved by the same argument].** (a) *Several overconvergent punctures.* If $S\subset(-\infty,1)$ is finite and $f$ is holomorphic on $\mathbf C\setminus[1,\infty)$ with meromorphic continuation on $\mathbf P^1\setminus(\{0,1,\infty\}\cup S)$, then $f$ has no monodromy at any point of $S$ on any sheet: take $\Sigma^0=\{0\}\cup S\cup\iota(S)$. This is what makes the two-fold (elliptic) hosts of §8 admissible. (b) *Slightly larger denominators.* The contradiction $5\le3.2452/(\log8-\tau)$ persists for $\tau<1.4303$, i.e. for $f,f\circ\iota$ of any type of rate $\sigma\le1.579$ (for instance $[1..n][1..n/2][1..n/13]$), with the same three pure functions.

---

## 2. Denominators: the half-integration lemma

**Lemma 2.1 [proved].** Let $h=\sum_ih_it^i\in\mathbf Z[[t]]$ and $k=\sum_{j\ge0}k_jt^j\in\mathbf Q[[t]]$ with $[1..j]\,k_j\in\mathbf Z$ for all $j$. Then $J(x):=\int_0^xh(t)k(t)\,dt$ has type $[1..n][1..n/2]$, and so does $h\cdot J$.

*Proof.* $[x^n]J=\frac1n\sum_{j=0}^{n-1}h_{n-1-j}k_j$. Fix a prime $p$ and write $e=\lfloor\log_pn\rfloor$, $e'=\lfloor\log_p(n/2)\rfloor$, $f=v_p(n)$. Since $j\le n-1$, $v_p(k_j)\ge-\lfloor\log_p(n-1)\rfloor=:-g$, hence $v_p([x^n]J)\ge-(f+g)$, and we must show $f+g\le e+e'$. Note $e'\ge e-1$ always ($n/2\ge p^e/2\ge p^{e-1}$).
*Case $n=p^e$:* $f=e$, $g=e-1$, $e'=e-1$; equality.
*Case $n\ne p^e$:* then $n-1\ge p^e$, so $g=e$. If $f=e$ then $p^e\mid n$ and $n\ne p^e$ force $n\ge2p^e$, so $e'=e$ and $f+g=2e=e+e'$. If $f\le e-1$ then $f+g\le2e-1\le e+e'$.
Multiplying by $h\in\mathbf Z[[t]]$ does not change the type because $[1..n'][1..n'/2]$ divides $[1..n][1..n/2]$ for $n'\le n$. $\square$

**Corollary 2.2 [proved].** On every host $H_A=1/\sqrt{1\mp4mx}$:
* $H_A$ has type $1$;
* $H_B:=H_A\int_0^xH_A\,\frac{dt}{1-t}$ has type $[1..n]$ (the coefficient of $x^n$ in $\int_0^xH_A/(1-t)$ is $\frac1n\sum_{i<n}\binom{2i}im^i$);
* $H_D:=H_A\int_0^xH_A\,\frac{\log(1-t)}{1-t}\,dt$ has type $[1..n][1..n/2]$, because $\frac{\log(1-t)}{1-t}=-\sum_{j\ge1}\mathcal H_jt^j$ with harmonic numbers $\mathcal H_j$, and $[1..j]\mathcal H_j\in\mathbf Z$;
* more generally $H[k]:=H_A\int_0^xH_Ak\,dt$ has type $[1..n][1..n/2]$ for every $k=r(t)\log(1-t)$ with $r\in\mathbf Z[t,\frac1{1-t}]$, and for $k=r(t)\in\mathbf Z[t,\frac1{1-t}]$ it has type $[1..n]$.

**[verified exact, $n\le400$, $m\in\{1,2,3,4,5,6,7,11,17,41\}$]** $H_A$ integral, $H_B$ exactly type $[1..n]$, $H_D$ exactly type $[1..n][1..n/2]$ (excess $1$ throughout), $H_C$ exactly type $[1..n]^2$ with excess over $[1..n][1..n/2]$ equal to $\bigl(\prod_{2n/3<p\le n}p\bigr)\cdot(\text{small cofactor})=e^{n/3+o(n)}$ (agent report `mixed/REPORT.md` §1). The kernel scan (§2 there, 88 kernels $t^i(1-t)^{-j}\log^\ell(1-t)$, outer factors $H_A,H_A^3,H_A/(1-x)$, $m=1,2$, $n\le200$–$400$) shows an exact ladder: $H[t^i(1-t)^{-j}\log^\ell(1-t)]$ has exact type $[1..n][1..n/2]\cdots[1..n/(\ell+1)]$ for $j\ge1$, the outer factor never changes the type, and the $1/t$ kernels are the only ones that fall off the ladder (to $[1..n]^2$). Two facts beyond Lemma 2.1: $H[\log^2(1-t)]$ (no $(1-t)^{-1}$) is exactly $[1..n][1..n/2]$ although $\log^2(1-t)$ is not of type $[1..j]$, and its period is not new: $c[\log^2(1-t)]=(4-2D'c_B+D'c_D)/m$ (PSLQ, $10^{-150}$, six values of $m$).

**Why $H_C$ fails.** For $k=\log(1-t)/t=-\sum_{j\ge0}t^j/(j+1)$ the hypothesis $[1..j]k_j\in\mathbf Z$ fails at $j=p-1$: the $1/p$ of $k_{p-1}$ meets the $1/p$ of the outer $\frac1n$ at $n=p$. Indeed $H_C=-x-\tfrac{13}4x^2-\tfrac{197}{18}x^3-\cdots$ on the $m=1$ host, and its excess over $[1..n][1..n/2]$ is unbounded (it is $2\cdot7$ at $n=8$, $2\cdot11\cdot13$ at $n=16$, …; it contains the primes in $(n/2,n]$), while its excess over $[1..n]^2$ is $1$. So $H_C$ is a genuine $\tau=2$ function. **[verified, $n\le250$]**

---

## 3. The fold and the sheets

Fix $m\ge1$; the case $1/\sqrt{1+4mx}$ is identical with $m\to-m$ and $D'\to-(4m+1)$, and is written out in §4.2. Put $\delta=\tfrac1{4m}$, $D'=4m-1$, and for $x$ near $[0,\delta]$ let $u=\sqrt{1-4mx}$, so $x=\frac{1-u^2}{4m}$, $H_A=\frac1u$, $1-x=\frac{D'+u^2}{4m}$.

**Lemma 3.1 (fold decomposition) [proved].** Let $k\in\mathbf Q[[t]]$ extend holomorphically to $\mathbf C\setminus\{1\}$ (all kernels of Corollary 2.2 do) and set $\kappa(v):=k\bigl(\frac{1-v^2}{4m}\bigr)$, an even function of $v$, holomorphic on $|v|<\sqrt{D'}$. Then
$$H[k](x)=\frac1{2mu}\int_u^1\kappa(v)\,dv=c[k]\,H_A(x)-\frac1{2m}\cdot\frac1u\int_0^u\kappa(v)\,dv,\qquad c[k]:=\frac1{2m}\int_0^1\kappa(v)\,dv=\int_0^\delta H_A(t)k(t)\,dt .$$
The last term, $E_k(x)$, is an even holomorphic function of $u$ on $|u|<\sqrt{D'}$, hence a holomorphic function of $x$ on the disc $|x-\delta|<D'/(4m)=1-\delta$; this disc contains $0$ and has the puncture $1$ on its boundary.

*Proof.* Substitute $t=\frac{1-v^2}{4m}$ in $\int_0^xH_Ak\,dt$: $H_A(t)dt=\frac1v\cdot\bigl(-\frac{v\,dv}{2m}\bigr)$ and $t\colon0\to x$ is $v\colon1\to u$. Since $\kappa$ is even, $\int_0^u\kappa$ is odd in $u$ and $\frac1u\int_0^u\kappa$ is even. $\square$

Consequently $T_\delta H[k]=H[k]-2c[k]H_A$ for the local monodromy $T_\delta$ at the fold ($T_\delta H_A=-H_A$), and

**Corollary 3.2 [proved].** If $a_0+\sum_ia_ic[k_i]=0$ with $a_i\in\mathbf Q$, then $f:=a_0H_A+\sum_ia_iH[k_i]=\sum_ia_iE_{k_i}$ is holomorphic on the disc $|x-\delta|<1-\delta$ and hence, being single-valued on the simply connected slit plane and singular there only possibly at $1$, holomorphic on $\mathbf C\setminus[1,\infty)$; it continues meromorphically (indeed holomorphically) along all paths in $\mathbf P^1\setminus\{0,\delta,1,\infty\}$; and it has type $[1..n][1..n/2]$ by Corollary 2.2 (clear denominators of the $a_i$). So Theorem $(\ast)$ applies to $f$.

### 3.3 The sheets reached by looping around $1$

Let $T_1$ denote continuation along a loop based near $0$ that encircles $1$ once (and not $\delta$), and $\epsilon=\pm1$ its orientation. Write $H_A(1)$ for the value of the (unchanged, since the loop does not enclose $\delta$) branch of $H_A$ at $1$: $H_A(1)^2=1/(1-4m)=-1/D'$, so $H_A(1)=\pm i/\sqrt{D'}\ne0$.

**Lemma 3.3 [proved].** With $I_B=\int_0^xH_A/(1-t)$ and $I_D=\int_0^xH_A\log(1-t)/(1-t)$:
$$T_1H_A=H_A,\qquad T_1H_B=H_B+r_BH_A,\quad r_B=-2\pi i\epsilon H_A(1)\ne0,\qquad T_1H_D=H_D+2\pi i\epsilon\,H_B+c_1H_A$$
for some constant $c_1$.

*Proof.* Near $t=1$, $H_A(t)=H_A(1)+(t-1)\tilde H(t)$ with $\tilde H$ analytic. So $I_B'=H_A/(1-t)$ has a simple pole with residue $-H_A(1)$, $I_B=-H_A(1)\log(1-x)+(\text{analytic})$, and $\log(1-x)\mapsto\log(1-x)+2\pi i\epsilon$ gives the first formula. Likewise $I_D=-\tfrac12H_A(1)\log^2(1-x)+B(x)\log(1-x)+A(x)$ with $A,B$ analytic at $1$, so $T_1I_D=I_D-2\pi i\epsilon H_A(1)\log(1-x)+(\text{analytic})=I_D+2\pi i\epsilon I_B+(\text{analytic})$. But $\{1,I_B,I_D\}$ is a basis of solutions of the third-order operator annihilating $I_D$ (whose homogeneous part annihilates $H_A/(1-t)$ and $H_A\log(1-t)/(1-t)$), and monodromy preserves the solution space with constant coefficients; so the analytic remainder is a constant $c_1$. Multiply by $H_A$. $\square$

**Proposition 3.4 [proved].** Let $f=aH_A+bH_B+dH_D$ with $a,b,d\in\mathbf Q$ and $a+bc_B+dc_D=0$ (so $T_\delta f=f$). If $(a,b,d)\ne0$ then some sheet of $f$ has nontrivial monodromy at $\delta$; precisely, for $k\ge1$,
$$(T_\delta-1)\,T_1^kf=-2H_A\Bigl[k\bigl(b\,r_B+d\,c_1+2\pi i\epsilon\,d\,c_B\bigr)+\pi i\epsilon\,d\,r_B\,k(k-1)\Bigr],$$
and the vanishing of this for $k=1,2$ forces $d=0$, then $b=0$, then $a=0$.

*Proof.* Iterating Lemma 3.3, $T_1^kH_D=H_D+2\pi i\epsilon k\,H_B+\bigl(kc_1+2\pi i\epsilon r_B\tbinom k2\bigr)H_A$ and $T_1^kH_B=H_B+kr_BH_A$. Apply $T_\delta-1$, which kills $f$ by hypothesis and sends $H_A\mapsto-2H_A$, $H_B\mapsto-2c_BH_A$. The coefficient of $k(k-1)$ is $\pi i\epsilon\,d\,r_B$ with $r_B\ne0$. $\square$

Since $H_A\ne0$ has a genuine branch point at $\delta$, a nonzero multiple of $H_A$ as $\delta$-monodromy contradicts the meromorphic continuability of $f$ along all paths in $\mathbf P^1\setminus\{0,1,\infty\}$ asserted by Theorem $(\ast)$. Hence:

**Theorem 3.5 [proved modulo $(\ast)$].** On each host $1/\sqrt{1\mp4mx}$ the numbers $1$, $c_B=c[\tfrac1{1-t}]$, $c_D=c[\tfrac{\log(1-t)}{1-t}]$ are $\mathbf Q$-linearly independent.

The same argument works for any finite family of kernels $k_i$ from Corollary 2.2 once the $T_1$-structure is written down (it is always unipotent with the $\log^2$-coefficient proportional to $H_A(1)\cdot(\text{leading coefficient of }k_i\text{ at }1)$); but by §5 no further kernel produces a period outside the span already controlled.

---

## 4. The periods in closed form

### 4.1 The imaginary family $1/\sqrt{1-4mx}$ (Theorem A)

With $a=\sqrt{D'}$, $\varphi=\arctan(1/a)=\theta_m/2$:
$$c_B=\frac1{2m}\int_0^1\frac{4m}{D'+v^2}dv=\frac2a\arctan\frac1a=\frac{\theta_m}{a},$$
$$c_D=2\int_0^1\frac{\log\frac{D'+v^2}{4m}}{D'+v^2}dv=\frac2a\Bigl[\varphi\log\frac{D'}{4m}+\int_0^\varphi\log(1+\tan^2\psi)\,d\psi\Bigr]=\frac2a\Bigl[\varphi\log\frac{D'}{4m}+2\varphi\log2-\mathrm{Cl}_2(\pi-2\varphi)\Bigr]=\frac{\theta_m\log\frac{D'}m-2\,\mathrm{Cl}_2(\pi-\theta_m)}{\sqrt{D'}},$$
using $\int_0^\varphi\log(2\cos\psi)\,d\psi=\tfrac12\mathrm{Cl}_2(\pi-2\varphi)$. **[numerical, 50 digits, $m=1,2,3,5,11$: PSLQ returns exactly $c_D=\theta\log D'/a-2\mathrm{Cl}_2(\pi-\theta)/a$ at $m=1$ and the displayed formula in general; see REPORT.md.]**

Theorem A is Theorem 3.5 with these values. Note $\tan\frac{\pi-\theta_m}2=\sqrt{D'}$ and $e^{i\theta_m}=\frac{(\sqrt{D'}+i)^2}{D'+1}=w_m$.

**[numerical]** $\mathrm{Cl}_2(\theta_m)$ and $\mathrm{Cl}_2(\pi-\theta_m)$ are **not** rational multiples of $\sqrt{D'}L(2,\chi_{-D'})$ for $m=2,3,5,11$ (PSLQ to 25 digits finds only garbage relations), as expected: $w_m\wedge(1-w_m)\ne0$ in $\Lambda^2\mathbf Q(\sqrt{-D'})^\times\otimes\mathbf Q$ for $m\ge2$ (for $m=2$: $w_2=-\bar\pi/\pi$ and $1-w_2=1/\pi$ with $\pi=\frac{1+\sqrt{-7}}2$, so $w_2\wedge(1-w_2)=-\bar\pi\wedge\pi\ne0$), so $w_m$ is not in the Bloch group and its dilogarithm is not a $\zeta_K(2)$-multiple. Only $m=1$ (roots of unity) gives an $L$-value.

### 4.2 The real family $1/\sqrt{1+4mx}$ (Theorem B)

Here $\delta=-\tfrac1{4m}$, $u=\sqrt{1+4mx}$, $x=\frac{u^2-1}{4m}$, $1-x=\frac{c^2-u^2}{4m}$ with $c^2=4m+1$; Lemma 3.1 and Proposition 3.4 go through verbatim (now $H_A(1)=1/c$ is real, $r_B=-2\pi i\epsilon/c\ne0$). The periods:
$$c_B=2\int_0^1\frac{dv}{c^2-v^2}=\frac1c\log\frac{c+1}{c-1},\qquad c_D=2\int_0^1\frac{\log\frac{c^2-v^2}{4m}}{c^2-v^2}dv=\frac{R_m}c,$$
with $R_m$ as in §0; the evaluation uses $\frac1{c^2-v^2}=\frac1{2c}\bigl(\frac1{c-v}+\frac1{c+v}\bigr)$, $\int\frac{\log(2c-w)}w\,dw=\log(2c)\log w-\mathrm{Li}_2\bigl(\frac w{2c}\bigr)$ and the reflection $\mathrm{Li}_2(z)+\mathrm{Li}_2(1-z)=\frac{\pi^2}6-\log z\log(1-z)$ to merge $\mathrm{Li}_2(\frac{c+1}{2c})$ into $\mathrm{Li}_2(\frac{c-1}{2c})$. **[numerical, 50 digits: PSLQ confirms, e.g. $c=3$: $3c_D=2\mathrm{Li}_2(\tfrac13)-\tfrac{\pi^2}6-\tfrac12\log^22+\log^23$; $c=5$: $5c_D=2\mathrm{Li}_2(\tfrac25)-\tfrac{\pi^2}6+\ldots$; $c=7$: $2\mathrm{Li}_2(\tfrac37)$; $c=9$: $2\mathrm{Li}_2(\tfrac49)$.]**

For non-square $4m+1$ the point $\frac{c-1}{2c}=\frac12-\frac1{2\sqrt{4m+1}}$ is a real quadratic irrationality; for $m=1$ it is $\frac{5-\sqrt5}{10}$ and $\mathrm{Li}_2$ there is not in the span of $\pi^2$ and products of $\log\varphi_{\rm golden},\log2,\log5$ **[numerical]**.

---

## 5. The whole $\tau=\tfrac32$ layer has a three-dimensional period space

**[proved]** Let $k=r(t)\log(1-t)$ or $k=r(t)$ with $r\in\mathbf Q[t,\frac1{1-t}]$ (Corollary 2.2). Under $t=\frac{1-v^2}{4m}$, $1-t=\frac{D'+v^2}{4m}$, the period is
$$c[k]=\frac1{2m}\int_0^1\tilde r(v^2)\bigl[\log(D'+v^2)-\log(4m)\bigr]dv,\qquad\tilde r(v^2)=r\bigl(\tfrac{1-v^2}{4m}\bigr)\in\mathbf Q\bigl[v^2,\tfrac1{D'+v^2}\bigr].$$
Split $\tilde r$ into a polynomial part and pole parts $(D'+v^2)^{-j}$.
* Polynomial part: $\int_0^1v^{2i}\log(D'+v^2)dv=\frac{\log(D'+1)}{2i+1}+(\mathbf Q+\mathbf Q\,D'^{\,i}\sqrt{D'}\theta_m)$, and the $\log(D'+1)=\log(4m)$ cancels against $-\log(4m)\int_0^1v^{2i}dv$. Contribution: $\mathbf Q+\mathbf Q\sqrt{D'}\theta_m$.
* $j=1$: $c_D$-type, contribution $\mathbf Q\,Q_m/\sqrt{D'}+\mathbf Q\,\theta_m/\sqrt{D'}$ after removing $\log(4m)c_B$ — this is where the $\log\frac{D'}{m}$ is born ($\log D'$ from $\int\log(D'+v^2)/(D'+v^2)$, $-\log4m$ from the kernel, $+2\log2$ from $\int\log\cos$).
* $j\ge2$: differentiate $G(s):=\int_0^1\frac{\log(s+v^2)}{s+v^2}dv=\frac1{\sqrt s}\bigl[\varphi\log s+2\varphi\log2-\mathrm{Cl}_2(\pi-2\varphi)\bigr]$ ($\varphi=\arctan s^{-1/2}$) in $s$: since $\frac{d}{ds}\mathrm{Cl}_2(\pi-2\varphi)=2\varphi'\log(2\cos\varphi)$ and $\log s-2\log\cos\varphi=\log(s+1)$, one finds $\partial_sG=-\tfrac12s^{-3/2}[\cdots]+s^{-1/2}\bigl(\varphi'\log(s+1)+\varphi/s\bigr)$, and again the stand-alone $\log(s+1)=\log(4m)$ cancels against the $-\log(4m)\int(D'+v^2)^{-j}$ term. Contribution: the same bracket $\varphi\log\frac{D'}m-\mathrm{Cl}_2(\pi-\theta_m)$ times rationals, plus $\mathbf Q\theta_m/\sqrt{D'}+\mathbf Q$.

Hence **every period of the layer lies in $V_m:=\mathbf Q\oplus\mathbf Q\frac{\theta_m}{\sqrt{D'}}\oplus\mathbf Q\sqrt{D'}\theta_m\oplus\mathbf Q\frac{Q_m}{\sqrt{D'}}$**, and since $\mathbf Q\frac{\theta_m}{\sqrt{D'}}+\mathbf Q\sqrt{D'}\theta_m=\theta_m\cdot\mathbf Q(\sqrt{D'})$, the strongest statement available from the layer is Theorem A together with the trivial $\theta_m\notin\mathbf Q(\sqrt{D'})$. The kernel $\log(1-t)$ itself realises $\sqrt{D'}\theta_m$: $c[\log(1-t)]=\frac{\sqrt{D'}\theta_m-2}{2m}$. Kernels with a factor $H_A^{2s}$, $s\ge1$, give divergent fold integrals; their $H[k]$ has a pole at $\delta$ after subtracting the finite-part multiple of $H_A$, and $(1-4mx)^sH[k]$ is again admissible; the periods are Hadamard finite parts and stay in $V_m$ (agent report, Task 2–3).

The pure Clausen value therefore cannot be reached at $\tau=\tfrac32$; reaching it at $\tau=2$ is exactly CDT's main theorem for $m=1$ (host $X_0(6)$, fourteen functions, symmetrisation), and for $m\ge2$ nothing analogous is available since the hosts $1/\sqrt{1-4mx}$ carry no second Eisenstein-type class.

---

## 6. Remark on CDT's Lemma "around 1/4" and Theorem "mixed"

CDT (arXiv v2 numbering: Lemma 2.11.13, Theorem 2.11.17; the type condition is (2.5.3) in Theorem 2.8.4, with no $A^{n+1}$ slack) define $H_C=\frac1{\sqrt{1-4x}}\int_0^x\frac{\log(1-t)}{t\sqrt{1-4t}}dt$ and assert it has type $[1..n][1..n/2]$. It does not: $H_C=-x-\tfrac{13}4x^2-\tfrac{197}{18}x^3-\cdots$, the excess over $[1..n][1..n/2]$ is unbounded (§2), and the true type is $[1..n]^2$ **[verified, $n\le250$]**. Their Theorem "mixed" claims the $\mathbf Q$-linear independence of $1,\ \pi/\sqrt3,\ \pi^2,\ 3L(2,\chi_{-3})-\frac\pi{\sqrt3}\log3$; the argument via $(\ast)$ only supports the three-element statement without $\pi^2$ (which is Theorem A at $m=1$; note $c_B=\pi/(3\sqrt3)$, $c_D=\frac\pi{3\sqrt3}\log3-L(2,\chi_{-3})$). By §5 no $[1..n][1..n/2]$-function of the layer has fold period $\pi^2$ (the value $-\pi^2/18$ of $\int_0^{1/4}\frac{\log(1-t)}{t\sqrt{1-4t}}dt$ is correct **[numerical]**, but its function is a $\tau=2$ object), so the $\pi^2$ clause needs a different input. This should be checked against the published version and, if confirmed, communicated to the authors.

---

## 7. Remarks and outlook

1. **Mahler measures [proved identity, verified $10^{-121}$].** $\mathrm{Cl}_2(\pi-\theta_m)=D(e^{i(\pi-\theta_m)})$ is the volume of the ideal tetrahedron with shape $-\bar w_m$. The Cassaigne–Maillot formula $\pi\,m(a+bx+cy)=\alpha\log a+\beta\log b+\gamma\log c+D(\tfrac ba e^{i\gamma})$ (triangle with sides $a,b,c$, opposite angles $\alpha,\beta,\gamma$) applied to the isosceles triangle $(1,1,\sqrt{D'/m})$, whose apex angle has cosine $1-\frac{D'}{2m}=-\frac{2m-1}{2m}$, i.e. equals $\pi-\theta_m$ exactly, gives
$$Q_m=2\pi\,m\bigl(1+x+\sqrt{D'/m}\,y\bigr)-\pi\log\frac{D'}m=\pi\Bigl[m\bigl(m(1+x)^2-(4m-1)y^2\bigr)-\log(4m-1)\Bigr],$$
using $2m(1+x+ty)=m((1+x)^2-t^2y^2)$ and $m(\lambda P)=\log\lambda+m(P)$. Hence **Theorem A is a statement about the Mahler measure of the integer polynomial $P_m=m(1+x)^2-(4m-1)y^2$**: $\pi\bigl(m(P_m)-\log(4m-1)\bigr)\notin\mathbf Q\sqrt{4m-1}+\mathbf Q\,\theta_m$, i.e. $1$, $\theta_m/\sqrt{4m-1}$ and $\pi\bigl(m(P_m)-\log(4m-1)\bigr)/\sqrt{4m-1}$ are $\mathbf Q$-linearly independent. For $m=1$, $P_1=(1+x)^2-3y^2$ and the identity is CDT's $m((1+x+y)^4/3)=4m(1+x+y)-\log3$ (Smyth); for $m\ge2$ the apex angle is irrational and no integer-sided triangle exists (a scan of all $31\,079$ triples $a\le b\le c\le60$ found nothing), which is why the algebraic triple is the right one. Fold-regularity was also confirmed numerically: the coefficients of $H_D-c_DH_A$ satisfy $n|a_n|\asymp\log n$ (radius exactly $1$, a $\log^2$ singularity at $1$), those of $H_B-c_BH_A$ have $n|a_n|$ constant; the cancellation is $\approx430$ digits at $n=400$, $m=3$.
2. **Where this sits in the programme, and why modular hosts do not join.** The modular Apéry companions $B=F\,\theta_q^{-r}\Phi$ are $\tau=w+1\ge2$ objects; the whole difficulty of the Catalan and $L(2,\chi)$ hosts in this repository is the $\tau=2$ count. The functions here are the $\tau=\tfrac32$ shadow: a *single* integration of $H_A$ against a type-$[1..j]$ kernel. Lemma 2.1 holds verbatim with $h=F$ a modular form with integral $t$-expansion, so $F\int F\,r(t)\log(1-t)\,dt$ has type $[1..n][1..n/2]$ on any modular host; but the fold mechanism does **not** transfer. What makes §3 work is that $H_A$ has *scalar* monodromy $-1$ at the fold, so that $T_\delta H[k]-H[k]$ is a multiple of one fixed function. At a cusp fold (all second-order rows: A, C, E, F) the monodromy of $F$ is unipotent, $T F=F+2\pi i\beta y_1$, and $T(F\!\int\!Fk)-F\!\int\!Fk$ contains the two new functions $y_1\!\int_{\delta}^xFk$ and $F\!\int_\delta^xy_1k$; no rational combination of $F$ and finitely many $F\!\int\!Fk_i$ is fold-regular. At an order-2 elliptic fold (third-order rows, Fricke point) a Fricke-odd $F$ does have scalar monodromy $-1$, but then the second singular point is the Galois conjugate of the fold (Apéry's $\zeta(3)$ host: $(\sqrt2\pm1)^4$) or sits at $\tfrac14$ with the fold at $\tfrac1{16}$ (Domb), and moving it to $1$ costs $4^n$ in the denominators — the width law once more. So $1/\sqrt{1\mp4mx}$ (and the equivalent two-sign-point hosts $1/\sqrt{(1-x)(1-c^2x)}$, which reproduce the real family) are essentially the only hosts on which Theorem $(\ast)$ bites.
3. **Not attempted:** quantitative versions; the $\tau=2$ layer (would need CDT's fine theorems and a symmetrisation of the host, which $\mathbf P^1\setminus\{0,\delta,1,\infty\}$ has only for $\delta\in\{-1,\tfrac12,2\}$, i.e. $m\in\{-\tfrac14,\tfrac12\}$ — not integers).

---

## 8. Two-fold (elliptic) hosts — the extension beyond CDT's setting

*Theory [proved]; numerics pending (agent report `scratchpad/elliptic/REPORT.md`).*

Let $H=1/\sqrt{P(x)}$ with $P\in\mathbf Z[x]$, $P(0)=1$, $P$ squarefree with all roots real, $H\in\mathbf Z[[x]]$ (e.g. $P=(1-x)(1-9x)(1-4mx)$, a product of factors each of which has integral $P_i^{-1/2}$). Let $\delta_1<\dots<\delta_s$ be the roots of $P$ in $(0,1)$ and suppose $1$ is either a root of $P$ or not; the curve $y^2=P(x)$ has genus $\lfloor(\deg P-1)/2\rfloor$ — genus $1$ for $\deg P=3,4$. With kernels $k$ as in Corollary 2.2, $H[k]=H\int_0^xHk$ has type $[1..n][1..n/2]$ (Lemma 2.1, unchanged), $c_i[k]:=\int_0^{\delta_i}Hk\,dt$.

**Lemma 8.1 (local fold decomposition) [proved].** Near $\delta_i$ put $t=\delta_i+s^2$. Then $H$ is an odd function of $s$, $Hk\,dt=(\text{even in }s)\,ds$, $\int_{\delta_i}^xHk$ is odd in $s$, and $H[k]=c_i[k]\,H+(\text{even in }s)$; so $T_{\delta_i}H[k]=H[k]-2c_i[k]H$ at every fold, with the same scalar structure as in §3.

Hence $f=a_0H+\sum_ja_jH[k_j]$ is holomorphic on $\mathbf C\setminus[1,\infty)$ iff $a_0+\sum_ja_jc_i[k_j]=0$ for **every** $i=1,\dots,s$ ($s$ linear conditions), and then Theorem $(\ast)$ in its generalised form (a) applies with $S=\{\delta_1,\dots,\delta_s\}$, and the sheet argument of §3.3 (which only used the puncture $1$ and one fold) shows $f$ has monodromy at $\delta_1$ on the sheets reached by looping around $1$ unless $(a_j)=0$ — when $1$ is itself a root of $P$ the loop around $1$ acts by a sign on $H$ and the bookkeeping changes but the conclusion persists (to be written out with the numerics).

**The $P(1)=0$ bookkeeping [proved].** Near $t=1$ write $H=(1-t)^{-1/2}g(t)$, $g$ analytic and nonvanishing. Then $I_B=\int_0^xH/(1-t)=(1-x)^{-1/2}A_1(x)+A_2(x)$ and $I_D=\int_0^xH\log(1-t)/(1-t)=(1-x)^{-1/2}[A_3(x)\log(1-x)+A_4(x)]+A_6(x)$ with $A_i$ analytic at $1$ (all exponents produced by integrating $(1-t)^{-3/2+j}\log^{0,1}(1-t)$ are half-integral, so no term $\log(1-x)$ without $(1-x)^{-1/2}$ occurs), $A_1(1)=A_3(1)=2g(1)$. The loop $T_1$ acts by $H\mapsto-H$, $(1-x)^{-1/2}\mapsto-(1-x)^{-1/2}$, $\log\mapsto\log+2\pi i\epsilon$, so
$$T_1H_B=H_B-2A_2H,\qquad T_1H_D=H_D+2\pi i\epsilon\frac{A_3}{A_1}H_B-2A_6H+\dots$$
and since $T_1$ preserves the solution spaces $\{H_B+cH\}$ and $\{H_D+\alpha H_B+\beta H\}$ of the inhomogeneous first-order equations, the coefficients are constants: $T_1H_B=H_B+r_B'H$ with $r_B'=-2A_2(1)=-2\,\mathrm{FP}\!\int_0^1\frac{H\,dt}{1-t}$ (Hadamard finite part), and $T_1H_D=H_D+2\pi i\epsilon H_B+c_1H$. This is the structure of Lemma 3.3 with $r_B$ replaced by $r_B'$, and Proposition 3.4 goes through verbatim provided $r_B'\ne0$, i.e. **the finite part of the divergent third-kind integral $\int_0^1\frac{dt}{(1-t)y}$ is nonzero** — a number to be computed per host (it is $-\tfrac12r_B'$; for the genus-$0$ hosts of §4.2 with $P=(1-x)(1-c^2x)$ one has $I_B=\tfrac14(1-u)$, $u=\sqrt{(1-c^2x)/(1-x)}$, so $A_2\equiv\tfrac14\ne0$).

**Theorem 8.2 [proved modulo $(\ast)$ and $\mathrm{FP}\int_0^1\frac{dt}{(1-t)y}\ne0$].** For such a host and kernels $k_1,\dots,k_r$, the $s\times(r+1)$ matrix $M=\bigl(1\mid c_i[k_j]\bigr)_{i,j}$ has no nonzero rational vector in its kernel.

For $s=2$, $r=2$ ($k_1=\frac1{1-t}$, $k_2=\frac{\log(1-t)}{1-t}$) the kernel of $M$ is spanned by the vector of $2\times2$ minors
$$\bigl(c_1[k_1]c_2[k_2]-c_2[k_1]c_1[k_2],\ \ c_1[k_2]-c_2[k_2],\ \ c_2[k_1]-c_1[k_1]\bigr),$$
so the statement is: **these three numbers are not $\mathbf Q$-proportional to a rational vector**, i.e. they span a $\mathbf Q$-space of dimension $\ge2$. The last two are integrals *between two branch points* of the curve $y^2=P(x)$:
$$\Delta_B=\int_{\delta_2}^{\delta_1}\frac{dt}{(1-t)y},\qquad \Delta_D=\int_{\delta_2}^{\delta_1}\frac{\log(1-t)\,dt}{(1-t)y},$$
a period of a differential of the third kind on the elliptic curve, and a "log-twisted" period — the kind of integral that computes Mahler measures of polynomials defining elliptic curves (Boyd, Rodriguez-Villegas) and elliptic dilogarithms. So the two-fold hosts move the method from $\mathbf P^1\setminus\{0,1,\infty\}$-periods (polylogarithms at algebraic points) to genuinely elliptic mixed periods, at no extra cost in the holonomy count. What these numbers are for the concrete hosts $E_m\colon y^2=(1-x)(1-9x)(1-4mx)$ is the content of the elliptic report.

**Genus 0 is vacuous [proved + numerical].** For the two-fold host $H=1/\sqrt{(1-4m_1x)(1-4m_2x)}$ (branched at two finite points only, so $y^2=P$ is rational) the difference periods are elementary: with $\delta_1=\frac1{4m_2}<\delta_2=\frac1{4m_1}$, the substitution $t=\frac{\delta_1+\delta_2}2-\frac{\delta_2-\delta_1}2\cos\phi$ gives $dt/|y|=d\phi/(4\sqrt{m_1m_2})$ and $1-t=a+b\cos\phi$ with $a=1-\frac{m_1+m_2}{8m_1m_2}$, $b=\frac{\delta_2-\delta_1}2$, $s:=\sqrt{a^2-b^2}=\frac{\sqrt{D_1D_2}}{4\sqrt{m_1m_2}}$ ($D_i=4m_i-1$), and the classical integrals $\int_0^\pi\frac{d\phi}{a+b\cos\phi}=\frac\pi s$, $\int_0^\pi\frac{\log(a+b\cos\phi)}{a+b\cos\phi}d\phi=\frac\pi s\log\frac{2s^2}{a+s}$ yield
$$\Delta_B=\int_{\delta_1}^{\delta_2}\frac{dt}{(1-t)|y|}=\frac{\pi}{\sqrt{D_1D_2}},\qquad \Delta_D=\int_{\delta_1}^{\delta_2}\frac{\log(1-t)\,dt}{(1-t)|y|}=\frac{\pi}{\sqrt{D_1D_2}}\log\frac{2s^2}{a+s}$$
(PSLQ confirms the second for $(m_1,m_2)=(1,2),(1,3),(2,3),(1,5)$; the values are $\Delta_B=i\pi/\sqrt{21}$ etc. on the imaginary branch). The second-fold conditions of Theorem 8.2 then only involve $\Delta_D/\Delta_B=\log(\text{algebraic})\notin\mathbf Q$ (Lindemann), so the two-fold statement is empty in genus $0$: the dilogarithm lives entirely in the individual fold periods $c_i$, and the *difference* of two folds sees only the residue structure at the puncture. Genus $\ge1$ is where $\Delta_B,\Delta_D$ become genuine periods.

### 8.1 Genus one: what the two-fold theorem says [numerical identifications at $10^{-90}$, agent report `lattice/mixed_periods/agent_reports/REPORT_elliptic_hosts.md`]

Hosts $E_m\colon P=(1-x)(1-9x)(1-4mx)$ ($m=1,\dots,6,12$), $F_m\colon P=(1-x)(1-25x)(1-4mx)$ ($m=1,2,3$), and $G\colon P=(1-4x)(1-8x)(1-12x)$. All eleven have $H\in\mathbf Z[[x]]$, $H_B$ of exact type $[1..n]$, $H_D$ and $H_L=H[\log(1-t)]$ of exact type $[1..n][1..n/2]$ ($n\le300$). Conductors $30,42,198,210,6270,2070,10998$ ($E_m$), $630,3570,4290$ ($F_m$), $24$ ($G$, isogenous to Boyd's $k=2,8$ curves). On the real branch continued from the upper half-plane, $H(u+i0)=i^{\#\{\text{branch points}<u\}}|P(u)|^{-1/2}$, so the near-fold period is real, $c_{\rm near}[k]=A[k]$, and the far-fold period is $c_{\rm far}[k]=A[k]+iB[k]$ with $B[k]=\int_{r_1}^{r_2}k\,dt/|y|$.

**Identities for $E_m,F_m$** (write $\{a,b\}=\{9\text{ or }25,\,4m\}$, $\Delta=(a-1)(b-1)$, $f=ab/\Delta$, $c=1/f=(a-1)(b-1)/(ab)$, and $\Omega=\int_{r_1}^{r_2}dx/y$, $\eta=\int_{r_1}^{r_2}x\,dx/y$ the period and quasi-period over the fold cycle):
$$B_B=f(\Omega-\eta),\qquad B_L=\tfrac12\Omega\log c,\qquad B_D=\frac{(a+b)\Omega-2ab\,\eta}{\Delta}+\tfrac f2(\Omega-\eta)\log c,\qquad A_B=f(\Omega_{01}-\eta_{01})-\tfrac2\Delta,$$
all ten hosts, residuals $\lesssim10^{-90}$. Mechanism: $x=1$ is a Weierstrass point, translation by the $2$-torsion point $(1,0)$ is an involution $\sigma$ swapping the two folds with $(1-x)(1-\sigma x)=c$, so averaging $\log(1-x)$ over $\{x,\sigma x\}$ gives the second identity; the first is the reduction of $dx/((1-x)y)$ — a second-kind differential, double pole without residue at the Weierstrass point — to $\{dx/y,x\,dx/y\}$. The near-fold periods $A_D,A_L$ are **not** in $\mathbf Q\langle1,\Omega_{01},\eta_{01},\Omega,\eta,\Omega\log c,\eta\log c\rangle$ (PSLQ, $70$ digits): they are the genuinely new mixed elliptic transcendentals. For $G$ (puncture not a branch point, $dx/((1-x)y)$ of the third kind) nothing identifies linearly; the Riemann bilinear relation $B_B\Omega_{23}-C_B\Omega_{12}=-2\pi u/\sqrt{231}$ holds exactly, $u$ the elliptic logarithm of the point over $x=1$.

**Theorem C (two-fold dichotomy) [proved modulo $(\ast)$; $\mathrm{FP}\ne0$ verified on all hosts].** For $E_m,F_m$ a rational $(\alpha,\beta,\delta)\ne0$ with $\alpha+\beta c_i[B]+\delta c_i[D]=0$ at both folds would need $\beta B_B+\delta B_D=0$ and $\alpha+\beta A_B+\delta A_D=0$; by the identities the first reads $\log c=(p\Omega+4\eta)/(\Omega-\eta)$ with $p=-(2\beta/\delta+2(a+b)/(ab))\in\mathbf Q$ (and $\delta=0$ forces $\beta=\alpha=0$). Hence: **for each host, either $\log\frac{(a-1)(b-1)}{ab}\ne\frac{p\,\Omega+4\eta}{\Omega-\eta}$ for every rational $p$, or, for the unique $p_0=((\Omega-\eta)\log c-4\eta)/\Omega$ (if rational), $A_D-\bigl(\tfrac{p_0}2+\tfrac{a+b}{ab}\bigr)A_B\notin\mathbf Q$.** Numerically $p_0$ is not a rational of height $\le10^{12}$ on any host, so the first alternative is what is (presumably) true, and Theorem C then says nothing about $A_D$. One cannot force the far-fold condition rationally: among the second-kind kernels $(1-t)^{-j}$, $j=1,2,3$, the unique rational combination with vanishing far-fold difference is the exact differential $d(R\sqrt P)$, whose near-fold period is algebraic. So in genus one the method yields a dichotomy, not an unconditional statement about the new periods; the structural content is the identities above.

---

## 9. From square roots to $k$-th roots: Kummer hosts $(1-k^kmx)^{-1/k}$

*Theory [proved]; numerics: agent report `scratchpad/kummer/REPORT.md` (pending).*

The mechanism of §3 used only that $H_A$ has rank one with scalar fold monodromy. The same is true of
$$H(x)=(1-k^kmx)^{-1/k},\qquad k\ge2,\ m\in\mathbf Z\setminus\{0\},\qquad\delta=\frac1{k^km},\qquad D:=k^km-1,$$
with monodromy $\zeta_k^{-1}$ at $\delta$.

**Lemma 9.1 (integrality) [proved].** $(1-Nx)^{-j/k}\in\mathbf Z[[x]]$ for $N=k^km$, $m\in\mathbf Z$, and all $j\ge0$. *Proof.* $[x^n]=\frac{(j/k)_n}{n!}N^n=\frac{\prod_{i<n}(j+ki)}{n!}\,k^{(k-1)n}m^n$. For $p\nmid k$ the arithmetic progression $j,j+k,\dots,j+(n-1)k$ has difference prime to $p$, so among any $p^e$ consecutive terms exactly one is divisible by $p^e$, whence $v_p(\prod)\ge\sum_e\lfloor n/p^e\rfloor=v_p(n!)$. For $p\mid k$, $v_p(n!)<n/(p-1)\le n\le(k-1)n\,v_p(k)$. $\square$ **[verified exact, $n\le150$, $k\le6$, $j<k$]**

**Sectors.** For $a\in\{1,\dots,k-1\}$ and a kernel $\kappa$ put $H^{(a)}[\kappa]:=H^a\int_0^xH^{k-a}\kappa\,dt$ and $c^{(a)}[\kappa]:=\int_0^\delta H^{k-a}\kappa\,dt$. By Lemma 2.1 (with $h=H^{k-a}\in\mathbf Z[[t]]$ and outer factor $H^a\in\mathbf Z[[t]]$), $H^{(a)}[\tfrac1{1-t}]$ has type $[1..n]$ and $H^{(a)}[\tfrac{\log(1-t)}{1-t}]$, $H^{(a)}[\log(1-t)]$ have type $[1..n][1..n/2]$ **[verified exact, excess $1$, $n\le150$, $k\le6$, $a=1,2$]**.

**Lemma 9.2 (fold decomposition) [proved].** With $u=(1-k^kmx)^{1/k}$, $x=\frac{1-u^k}{k^km}$, $dt=-\frac{u^{k-1}du}{k^{k-1}m}$, $H^{k-a}=u^{-(k-a)}$, and $\tilde\kappa(v):=\kappa\bigl(\frac{1-v^k}{k^km}\bigr)$ (a function of $v^k$, holomorphic for $|v|<D^{1/k}$):
$$H^{(a)}[\kappa]=\frac{u^{-a}}{k^{k-1}m}\int_u^1v^{a-1}\tilde\kappa(v)\,dv=c^{(a)}[\kappa]\,H^a-\frac{u^{-a}}{k^{k-1}m}\int_0^uv^{a-1}\tilde\kappa(v)\,dv,\qquad c^{(a)}[\kappa]=\frac1{k^{k-1}m}\int_0^1v^{a-1}\tilde\kappa(v)\,dv,$$
and $u^{-a}\int_0^uv^{a-1}\tilde\kappa(v)dv$ is a holomorphic function of $u^k$, i.e. of $x$, on $|x-\delta|<1-\delta$. Hence $T_\delta H^{(a)}[\kappa]=H^{(a)}[\kappa]+(\zeta_k^{-a}-1)c^{(a)}[\kappa]H^a$, and a rational relation $a_0+\sum_ia_ic^{(a)}[\kappa_i]=0$ makes $f=a_0H^a+\sum_ia_iH^{(a)}[\kappa_i]$ holomorphic on $\mathbf C\setminus[1,\infty)$ with type $[1..n][1..n/2]$. Explicitly
$$c^{(a)}_B=k\int_0^1\frac{v^{a-1}dv}{D+v^k},\qquad c^{(a)}_D=k\int_0^1\frac{v^{a-1}\log\frac{D+v^k}{k^km}}{D+v^k}dv,\qquad c^{(a)}[\log(1-t)]=\frac1{k^{k-1}m}\int_0^1v^{a-1}\log\frac{D+v^k}{k^km}\,dv .$$

**Lemma 9.3 (sheets) [proved].** $H$ is analytic at $1$ with $H^{k-a}(1)=(-D)^{-(k-a)/k}\ne0$; Lemma 3.3 holds verbatim with $H_A$ replaced by $H^a$: $T_1H^{(a)}_B=H^{(a)}_B+rH^a$, $r=-2\pi i\epsilon H^{k-a}(1)\ne0$, and $T_1H^{(a)}_D=H^{(a)}_D+2\pi i\epsilon H^{(a)}_B+c_1H^a$. Since $T_\delta-1$ multiplies $H^a$ by $\zeta_k^{-a}-1\ne0$, the computation of Proposition 3.4 goes through unchanged.

**Theorem D [proved modulo $(\ast)$].** For every $k\ge2$, every $m\in\mathbf Z\setminus\{0\}$ and every sector $a\in\{1,\dots,k-1\}$, the three numbers
$$1,\qquad k\int_0^1\frac{v^{a-1}\,dv}{k^km-1+v^k},\qquad k\int_0^1\frac{v^{a-1}\log\frac{k^km-1+v^k}{k^km}}{k^km-1+v^k}\,dv$$
are linearly independent over $\mathbf Q$. ($k=2$, $a=1$, $m\gtrless0$ are Theorems A and B.) By partial fractions over the roots $\alpha_j=D^{1/k}e^{i\pi(2j+1)/k}$ of $v^k=-D$, the second number is a $\mathbf Q(\zeta_{2k},D^{1/k})$-combination of $\pi$ and logarithms of elements of that field, and the third a combination of dilogarithms at the points $\alpha_j/(\alpha_j-\alpha_l)$, $(\alpha_j-1)/(\alpha_j-\alpha_l)$ and products of such logarithms. For $k=3$, $m=1$ the field is $\mathbf Q(\omega,\sqrt[3]{26})$: the statements leave the quadratic world for pure cubic fields. The $k-1$ sectors of one host give $k-1$ separate statements (their fold monodromies are multiples of the independent functions $H^a$, so a relation mixing sectors is not fold-regular); whether periods of different sectors are $\mathbf Q$-related is being tested.

**Refinement (integrality threshold) [proved].** The proof of Lemma 9.1 only needs $v_p(N/k)\ge1$ for the primes $p\mid k$ (then $v_p((N/k)^n)\ge n\ge v_p(n!)$), so $(1-Nx)^{-j/k}\in\mathbf Z[[x]]$ for all $N\in k\,\mathrm{rad}(k)\,\mathbf Z$: hosts $(1\mp9mx)^{-1/3}$, $(1\mp8mx)^{-1/4}$, $(1\mp25mx)^{-1/5}$, $(1\mp36mx)^{-1/6}$, with $D=N\mp1$ in Theorem D (the hypothesis $k^km$ there can be replaced by $k\,\mathrm{rad}(k)\,m$).

**Cyclotomic collapse [numerics pending].** When $N\mp1=n^k$ is a perfect $k$-th power the roots of $D\pm v^k$ are $n\cdot(\text{roots of unity})$, the dilogarithm points $c_j/(c_j-c_l)=1/(1-\zeta)$ and $(c_j-1)/(c_j-c_l)$ lie in $\mathbf Q(\zeta_{2k})$, and the fold-type points $1/(1-\zeta)$ are Bloch-group elements whose $D$-values are $L$-values: $k=3$ ($N=9$: $D=8$; $N=-63$: $64$; $N=126$: $125$) gives $\sqrt3L(2,\chi_{-3})$, $k=4$ ($N=-80$: $81$; $N=-624$: $625$) gives **Catalan's constant** $G=D((1+i)/2)$, $k=5$ ($N=-1025$: $4^5$; $N=7775$: $6^5$) gives $L(2,\psi)$ for the quartic character mod $5$, $k=6$ ($N=-15624$) gives $\mathbf Q(\zeta_{12})$-values. The base-point ($t=0$, $u=1$) contributions $(c_j-1)/(c_j-c_l)$ are generically not Bloch elements, so the periods are mixtures of an $L$-value with generic cyclotomic dilogarithms and logarithms; e.g. for $(1+80x)^{-1/4}$, sector $1$: $27c_B=2\arctan\tfrac13+\log2$ **[numerical, 50 digits]**, and $c_D$ involves $G$, $\mathrm{Cl}_2(\pi-2\arctan\tfrac13)$ and $\mathrm{Li}_2$ at $(1+i)/3$-type points (exact decomposition being computed).

### 9.1 Results of the Kummer computations [agent report `lattice/mixed_periods/agent_reports/REPORT_kummer_hosts.md`]

Conventions: $H=(1-\sigma Nx)^{-1/k}$ with $\sigma=\pm1$ ($\sigma=+1$ the "imaginary" family, $\sigma=-1$ the real one), $D=N-\sigma$, $Q(u)=u^k+\sigma D$, roots $\alpha_j$ of $Q$; then $1-t=|Q(u)|/N$ and $c^{(a)}[\kappa]=k\int_0^1u^{a-1}\tilde\kappa\,du/Q(u)$ up to an overall sign in the real family (the periods of the real family are negative in the normalisation of §9).

* **Integrality is exactly $k\,\mathrm{rad}(k)\mid N$** [verified exact, $N\le400$, $k\le7$, both signs]; all types exact (excess $1$) for $k\le7$, $m\le3$, $n\le300$, all sectors.
* **Closed form** [proved by partial fractions; verified 111 digits]: since the power sums $\sum_j\alpha_j^i$ vanish for $1\le i\le k-1$,
$$c^{(a)}_B=-\frac\sigma D\sum_{j=0}^{k-1}\alpha_j^{\,a}\log\Bigl(1-\frac1{\alpha_j}\Bigr)\qquad(\text{principal branches, no corrections}).$$
  For $k=2$, $N=4$: $\pi/(3\sqrt3)$. E.g. $k=4$, sector $2$, imaginary family: $c^{(2)}_B=\frac2{\sqrt D}\arctan\frac1{\sqrt D}$ (the $k=2$ period of the same $D$).
* **The dilogarithm period** $c^{(a)}_D$ is an explicit $\mathbf Q(\zeta_{2k},D^{1/k})$-combination of $\mathrm{Li}_2$ at the points $\alpha_j/(\alpha_j-\alpha_l)$, $(\alpha_j-1)/(\alpha_j-\alpha_l)$ and products of logarithms (derived symbolically, verified $\ge61$ digits).
* **Period space of a sector**: exactly $\mathbf Q+\mathbf Qc^{(a)}_B+\mathbf Qc^{(a)}_D$ (kernels $1,\frac1{1-t},\frac t{1-t},\log(1-t),\frac{\log(1-t)}{1-t},\frac{t\log(1-t)}{1-t}$), with elementary identities, e.g. $c^{(a)}[\log(1-t)]=\frac{kD}{Na}c^{(a)}_B-\frac{\sigma k^2}{Na^2}$ (integration by parts in $u$). **Different sectors are arithmetically independent** (no relation $c^{(b)}\in\mathbf Q+\mathbf Qc^{(a)}$ in 32 tests, 140 digits, bound $10^{12}$), so the $k-1$ statements of Theorem D per host are genuinely different.
* **Cyclotomic collapse (correcting the expectation above).** When $D=n^k$, the fold-type points $\zeta_l/(\zeta_l-\zeta_j)$ depend only on $\mu=\zeta_j/\zeta_l$ and their aggregate coefficient is $\propto\sum_l\zeta_l^a=0$: **they cancel identically** (verified to $10^{-320}$ on all collapse hosts). In particular **Catalan's constant does not appear** on the $k=4$ hosts $(1+80x)^{-1/4}$, $(1+624x)^{-1/4}$ (its coefficient is exactly $0$; the base-point Bloch–Wigner values there are not rational multiples of $G$). The only collapse host retaining an $L$-value is $(1-9x)^{-1/3}$ ($D=8$), where the base-point points collapse: the Bloch part of $c^{(1)}_D$ is $-\frac{15}{32}L(2,\chi_{-3})$ and of $c^{(2)}_D$ is $-\frac{15}{16}L(2,\chi_{-3})$ (300 digits). The correct shape of $c_D$ on collapse hosts is $\sum_z\bigl(\beta_z\,\mathrm{Re}\,\mathrm{Li}_2(z)+\gamma_z D(z)\bigr)+(\text{elementary})$ with an irreducible $\mathrm{Re}\,\mathrm{Li}_2$ part (PSLQ-excluded from $\{\pi^2,\log p\log q,\pi\log p,\theta\log p,\dots\}$ at 300 digits).
* **Mahler measures.** Every surviving Bloch–Wigner point of a collapse host is the Cassaigne–Maillot point of the triangle with vertices $1,n\zeta_l,n\zeta_j$ (sides $|n\zeta_l-n\zeta_j|,|n\zeta_l-1|,|n\zeta_j-1|$), so the Bloch–Wigner content of these periods is Mahler-measure content; for $(1-9x)^{-1/3}$ it is $-\frac5{24}\sqrt3\,\pi m(1+x+y)$ resp. $-\frac5{12}\sqrt3\,\pi m(1+x+y)$.

So Theorem D is a genuine extension to Kummer fields, but it does not produce new $L$-value statements: for $k\ge3$ the periods are mixed dilogarithms over $\mathbf Q(\zeta_{2k},D^{1/k})$, and in the cyclotomic-collapse cases the $L$-values cancel by the vanishing of the power sums of $\mu_k$ — the same mechanism ($\sum\zeta^a=0$) that makes the sectors independent.

### 9.2 Classification of rank-one hosts [proved]

**Lemma.** Let $f\in\mathbf Z[[x]]$ be an algebraic function whose only finite singularity is at one point $\delta\ne0$, with scalar local monodromy of exact order $k\ge2$ there. Then $f=R(x)(1-Nx)^{j/k}$ with $N=1/\delta\in\mathbf Z$, $\gcd(j,k)=1$, $R\in\mathbf Q(x)$ with poles only at $\delta$, and $k\,\mathrm{rad}(k)\mid N$; conversely all such series with $R(1-Nx)^{j/k}\in\mathbf Z[[x]]$ occur. *Proof.* $f^k$ is single-valued and algebraic on $\mathbf P^1\setminus\{\delta,\infty\}$, hence rational with zeros and poles only at $\delta,\infty$, so $f=c\,R_0(x)(1-x/\delta)^{j/k}$; write $\delta=p/q$ in lowest terms: the coefficient of $x^n$ in $(1-x/\delta)^{j/k}$ is $\binom{j/k}n(-q/p)^n$, whose $\ell$-adic valuation for $\ell\mid p$ is $-n\,v_\ell(p)+O(\log n)$, not compensated by a fixed rational factor, so $p=1$ and $\delta=1/N$. Then $[x^n](1-Nx)^{j/k}=\prod_{i<n}(j-ki)\,(N/k)^n/n!$ up to sign; for a prime $p\mid k$ (so $p\nmid j$) the product is a $p$-unit and $v_p$ of the coefficient is $n(v_p(N)-v_p(k))-v_p(n!)$, which at $n=p^e$ is negative unless $v_p(N)\ge v_p(k)+1$; the converse is Lemma 9.1 (the same estimate together with the arithmetic-progression argument for $p\nmid k$). $\square$ **[verified exact: the admissible $N\le400$ are exactly the multiples of $k\,\mathrm{rad}(k)$, $k\le7$.]** So the $\tau=\frac32$ theory of this note is complete: its hosts are exactly the Kummer functions $(1-Nx)^{-j/k}$, $k\,\mathrm{rad}(k)\mid N$, up to rational factors.

---

## 10. Effective irrationality measures (Theorem E)

CDT's quantitative refinement (ICM survey, Proposition `theomainreprise`): for an Apéry limit $\eta$ realised as $g=B-\eta A$ with $A,B\in\mathbf Q[[x]]$ of the given type, if $m$ independent functions $g_i=B_i-\eta A_i$ of the type have meromorphic pullbacks under $\varphi$, of which $\gamma m$ are unconditional ($A_i=0$), and $\varphi^*A_i,\varphi^*B_i$ all converge on $|z|<\rho$, then
$$m\le\frac{\mathrm{BC}(\varphi)}{\log|\varphi'(0)|-\tau(\mathbf b)-(1-\gamma)\bigl(\tfrac2\mu-\tfrac{1-\gamma}{\mu^2}\bigr)\log\tfrac1\rho},\qquad\mu=\mu_{\rm eff}(\eta),$$
whenever the denominator is positive. (As in `HOLONOMY_LINDEP.md` §6, which reproduces CDT's $24781$ from this formula, we apply it to the exhibited $m$-tuple.)

**Set-up on a Kummer host** (sector $a$; for $k=2$ drop the superscript). Three choices of pair $(A,B)$:
* $(H^a,\ H^{(a)}_D)$: $\eta=c^{(a)}_D$; inventory $\{1,\log(1-x),\log^2(1-x),g,g\circ\iota\}$, $m=5$, $\gamma=\frac35$ (the unconditional part is exactly $\{1,\log,\log^2\}$ by CDT's Theorem 2.11.11), types $(0,1,\frac32,\frac32,\frac32)$, $\tau=\frac{69}{50}$.
* $(H^a,\ H^{(a)}_B)$: $\eta=c^{(a)}_B$; inventory $\{1,\log(1-x),g,g\circ\iota\}$, $m=4$, $\gamma=\frac12$, types $(0,1,1,1)$, $\tau=\frac{15}{16}$.
* $(H^{(a)}_D,\ H^{(a)}_B)$: $\eta=c^{(a)}_B/c^{(a)}_D$, same data as the first.

$\varphi$ is the bivalent map, $\mathrm{BC}=\log8+4G/\pi$, and $\rho$ solves $\max_{|z|=\rho}|\varphi(z)|=8\rho(1+\rho^2)/(1-\rho)^4=1/N$ (the convergence radius of $H^a$, $H^{(a)}_B$, $H^{(a)}_D$; the $\iota$-images have larger radius). Solving $(1-\gamma)(\frac2\mu-\frac{1-\gamma}{\mu^2})L=\log8-\tau-\mathrm{BC}/m$ gives $\mu_*=(1-\gamma)/(1-\sqrt{1-\mathrm{gap}/L})$; for a putative larger true dimension $m'>m$ the corresponding thresholds are smaller (checked for $m'\le40$), so $m=5$ (resp. $4$) is binding.

**Theorem E [proved modulo CDT's quantitative bound, applied to the exhibited tuple].** With $\kappa_D(N),\kappa_B(N)$ as in the table, $\mu_{\rm eff}\bigl(c^{(a)}_D\bigr)\le\kappa_D(N)$, $\mu_{\rm eff}\bigl(c^{(a)}_B/c^{(a)}_D\bigr)\le\kappa_D(N)$ and $\mu_{\rm eff}\bigl(c^{(a)}_B\bigr)\le\kappa_B(N)$ for every sector:

| host | $N$ | $\rho$ | $L=\log(1/\rho)$ | $\kappa_D$ | $\kappa_B$ |
|---|---|---|---|---|---|
| $k=2$, $m=1$ | 4 | 0.02789 | 3.580 | **56.73** | 10.57 |
| $k=2$, $m=2$ | 8 | 0.01472 | 4.218 | 66.89 | 12.51 |
| $k=2$, $m=3$ | 12 | 0.01001 | 4.605 | 73.03 | 13.68 |
| $k=2$, $m=5$ | 20 | 0.00610 | 5.100 | 80.90 | 15.18 |
| $k=2$, $m=11$ | 44 | 0.00281 | 5.875 | 93.23 | 17.52 |
| $k=2$, $m=41$ | 164 | 0.00076 | 7.182 | 114.0 | 21.48 |
| $k=3$, $N=9$ | 9 | 0.01317 | 4.330 | 68.66 | 12.85 |
| $k=3$, $N=27$ | 27 | 0.00455 | 5.394 | 85.57 | 16.06 |
| $k=4$, $N=80$ | 80 | 0.00155 | 6.468 | 102.7 | 19.31 |
| $k=5$, $N=25$ | 25 | 0.00490 | 5.318 | 84.37 | 15.84 |

In particular $\mu_{\rm eff}\bigl(L(2,\chi_{-3})-\frac{\pi}{3\sqrt3}\log3\bigr)\le56.73$ (CDT's measure for the pure $L(2,\chi_{-3})$ is $24781$; no measure was known for the mixed number), $\mu_{\rm eff}\bigl(\frac{\theta_2\log\frac72-2\mathrm{Cl}_2(\pi-\theta_2)}{\sqrt7}\bigr)\le66.89$, and $\mu_{\rm eff}\bigl(\frac2{\sqrt7}\arctan\frac1{\sqrt7}\bigr)\le12.51$. For $m=1$ the elementary number is $\pi/(3\sqrt3)$, where Salikhov's $\mu(\pi/\sqrt3)\le4.6016$ is far better; for the mixed dilogarithm periods and for $\theta_m/\sqrt{D'}$ with $m\ge2$ we know of no previous explicit measure. The real family and all Kummer sectors have the same $\kappa$'s (they depend only on $N$). Script: `lattice/mixed_periods/effective_measures.py`.

---

## 11. Catalan and the angle–ratio law [proved]

On a rank-one host $H_A=(1-4mx)^{-1/2}$ (fold $\delta=\frac1{4m}$) with a kernel puncture at $s$ (kernel $\frac1{1-t/s}$, which needs $1/s\in\mathbf Z$), the substitution $u=\sqrt{1-4mx}$ gives $1-t/s=\frac{4ms-1+u^2}{4ms}$, so the dilogarithm points have angle $\theta$ with
$$\cos\theta=\frac{r-2}{r},\qquad r:=\frac{s}{\delta}\ \text{(puncture/fold ratio)},$$
and the bivalent count's entry for the mixed statement is $\log(2r)-\frac{69}{50}$ (in units $\delta=\frac14$: the bivalent map onto $\mathbf C\setminus\{s\}$ has $\varphi'(0)=8s$). Hence $\chi_{-3}$ ($\theta=\pi/3$) forces $r=4$ and entry $+0.70$; **Catalan ($\theta=\pi/2$) forces $r=2$ and entry $\log4-1.38=+0.006$** — the threshold angle is $\theta_*=\arccos(1-4e^{-69/50})\approx90.36^\circ$, and Catalan sits $0.36^\circ$ inside it. With $\mathrm{BC}=\log4+4G/\pi=2.552$ the count for "$G-\frac\pi4\log2\notin\mathbf Q+\mathbf Q\pi$" needs $m\ge405$ functions of type $[1..n][1..n/2]$ on $\mathbf P^1\setminus\{0,\frac12,\infty\}$, against the $5$ that exist. The same ratio $2$ is $|t_2|/|t_{\rm fold}|=(1/4)/(1/8)$ on the modular level-$8$ host, and the width law of `CATALAN_OBSTRUCTION.md` is the same statement in modular language. Equivalently, the host that would give Catalan the $\chi_{-3}$ geometry is $(1-2x)^{-1/2}$, which is not $2$-integral (Lemma 9.2 with $k=2$, $N=2$). Rank-one hosts being classified (§9.2), a rescue inside the holonomy method needs a rank-$\ge2$ host for the $\chi_{-4}$ class with ratio $\ge4$ outside the Beauville/four-term/AL/cyclotomic families already scanned.
