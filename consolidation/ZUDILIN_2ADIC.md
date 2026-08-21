# $\xi_2^{\rm Zud}=\zeta_2(2)$: **proved**

*Closes the last gap in the $5{:}8$ Catalan bridge.  Scripts: `lattice/zudilin_2adic/`.
Supersedes `SOURCES_S18_ZUDILIN.md` §6.3–§6.5 (see §8 below for the corrections).*

---

## 0. Verdict

**Theorem.** Let $(Q_m,P_m)$ be Zudilin's Catalan row (arXiv:math/0201024, Thm. 1).  Then
$$\boxed{\ \lim_{m\to\infty}\frac{P_m}{Q_m}=\zeta_2(2)\quad\text{in }\Q_2,\qquad
v_2\Bigl(\zeta_2(2)-\frac{P_m}{Q_m}\Bigr)=8m-1-4s_2(m)\ \text{exactly}. }$$

So $\xi_2^{\rm Zud}=\zeta_2(2)=2\,\xi_2^{\mathbf E}$, and Lemma "common $2$-adic Catalan
period" of `catalan-2-row-denominators/5-8 theorem/CATALAN_0902526_SHORT_PROOF.tex`
— its equation `\eqref{eq:Ztail}`, $v_2(G_2-P_m/Q_m)=8m-1-4s_2(m)$ — is now a **theorem**,
not a measurement.  The $5{:}8$ diagonal no longer rests on numerics on the Zudilin side.

The proof is *not* the argument sketched in the $5{:}8$ paper (that argument is invalid
as written, for the reason given in `SOURCES_S18_ZUDILIN.md` §6.4 — but the reason is
repairable, and the repair is an **exact identity**, found here, not an estimate).

**The missing identity.**  With $x_m:=\tfrac12-m$, $G_m:=\sum_{j=1}^m\frac{(-1)^{j-1}}{(2j-1)^2}$
and $(p_n,q_n)$ Beukers' Padé pair (§1),
$$\boxed{\ Q_m=q_m(x_m),\qquad P_m=\frac{(-1)^m}{8}\,p_m(x_m)+G_m\,q_m(x_m)\ }\tag{B}$$
for **all** $m\ge0$ (verified exactly for $m\le160$, `verify_bridge.gp`; proved in §4).
`SOURCES_S18_ZUDILIN.md` §6.3 recorded "$P_m\neq p_m(x_m)$, correction of bounded $v_2$"
and stopped there.  The correction is the *partial Catalan sum* $G_mQ_m$, and the sign
$(-1)^m$ is exactly the sign produced by iterating Beukers' shift equation from $x_0=\tfrac12$
down to $x_m$.  Once (B) is in hand the whole moving-point difficulty evaporates: the
target does move, but it moves by an **explicitly known elementary rational**, and $P_m$
already contains the compensating term.

Everything below is exact.  Two published inputs are used as black boxes: Beukers'
$\Theta_2(1/2)=8\zeta_2(2)$ (Prop. 5.1) and his closed form for the Padé remainder
(Prop. 6.1).  No numerical identification is used anywhere.

---

## 1. Notation

Beukers, *Irrationality of some $p$-adic $L$-values*, Acta Math. Sinica **24** (2008)
663–686 (author version `webspace.science.uu.nl/~beuke106/padicL.pdf`).  In
$K=\Q((1/x))$ put $\left[\frac nx\right]:=\dfrac{n!}{x(x+1)\cdots(x+n)}$ (Diamond) and
$$\Theta(x):=-\sum_{n\ge0}\Bigl[\frac nx\Bigr]\Bigl[\frac n{1-x}\Bigr]
=\frac1{x^2}+\frac1{x^3}+\frac1{x^4}+\cdots\in K .$$
Its Stieltjes continued fraction has convergents $p_n/q_n$ with
$$(n+1)^2w_{n+1}=(x^2-x+2n(n+1)+1)\,w_n-n^2w_{n-1},\qquad
\begin{aligned}q_0&=1,&q_1&=x^2-x+1,\\ p_0&=0,&p_1&=1,\end{aligned}\tag{1.1}$$
$\deg q_n=2n$, $\deg p_n=2n-2$, and (Beukers, Prop. 6.1, by a Zeilberger certificate)
$$r_n(x):=p_n(x)-q_n(x)\Theta(x)=(-1)^n\sum_{k\ge n}\binom kn\Bigl[\frac kx\Bigr]\Bigl[\frac k{1-x}\Bigr]
=O(x^{-2n-2}).\tag{1.2}$$
Also $q_n(x)=\sum_k\binom nk\binom{-x}k\binom{k-x}k$, and both families are invariant
under $x\mapsto1-x$ (the recursion sees only $x^2-x$), so $q_n,p_n\in\Q[x^2-x]$.

**Sign conventions (settled by exact computation, `tseries.gp`).**  The author version of
Beukers' paper is internally inconsistent by one global sign in §4–§5.  The self-consistent
normalisation — the one forced by $\Theta=\lim p_n/q_n$ and by the Diamond series above — is
$$\Theta(x)+\Theta(x+1)=+\frac2{x^2},\qquad \Theta_2(1/2)=+8\,\zeta_2(2),\qquad
\Theta(x)=-\sum_{n\ge0}(2^{n+1}-2)B_n(-1/x)^{n+1}\ (B_1=-\tfrac12).$$
The printed Cor. 4.5(3) ($-2/x^2$) and Prop. 5.1(1) ($-8\zeta_2(2)$) carry the opposite
sign; note that Beukers' own *proof* of Prop. 5.1(1) ends with "$=8\zeta_2(2)$".  Both
statements above were verified here to $2^{200}$ resp. $2^{261}$.  This is the $\mathcal T$
of the $5{:}8$ paper: $\mathcal T=\Theta$ in the corrected normalisation.

Zudilin's row: $p(n)=20n^2-8n+1$, $q(n)=3520n^6+5632n^5+2064n^4-384n^3-156n^2+16n+7$,
$$(2n{+}1)^2(2n{+}2)^2p(n)\,u_{n+1}=q(n)\,u_n+(2n{-}1)^2(2n)^2p(n{+}1)\,u_{n-1},\tag{1.3}$$
$Q_0=1,Q_1=\tfrac74$; $P_0=0,P_1=\tfrac{13}8$.

---

## 2. Lemma A (mixed Casoratian): where Zudilin's $p(n)$ comes from

Set, for $n\ge1$,
$$\Omega_n(x):=p_n(x)q_{n-1}(x{+}1)+p_{n-1}(x{+}1)q_n(x)-\frac2{x^2}\,q_n(x)q_{n-1}(x{+}1).$$

> **Lemma A.** $\displaystyle \Omega_n(x)=-\frac{(x-n)^2+n^2}{n^2x^2}$, i.e.
> $x^2\Omega_n(x)=-\Bigl(\frac{x^2}{n^2}-\frac{2x}{n}+2\Bigr)$.

*Proof.* By the shift equation, $\frac2{x^2}=\Theta(x)+\Theta(x+1)$, so
$$\Omega_n=\bigl[p_n(x)-q_n(x)\Theta(x)\bigr]q_{n-1}(x{+}1)+\bigl[p_{n-1}(x{+}1)-q_{n-1}(x{+}1)\Theta(x{+}1)\bigr]q_n(x)$$
$$=r_n(x)\,q_{n-1}(x{+}1)+r_{n-1}(x{+}1)\,q_n(x).\tag{2.1}$$
By (1.2) the first product is $O(x^{-2n-2})O(x^{2n-2})=O(x^{-4})$ and the second is
$O(x^{-2n})O(x^{2n})=O(1)$; hence $\Omega_n=O(1)$ at $\infty$.  On the other hand
$\Omega_n$ is a rational function whose only pole is $x=0$, of order $\le2$.  Therefore
$x^2\Omega_n$ is a polynomial of degree $\le2$, and three evaluations determine it.

*Constant term.* $x^2\Omega_n|_{x=0}=-2q_n(0)q_{n-1}(1)=-2$, since $q_m(0)=1$ (only $k=0$
survives in $\sum_k\binom mk\binom{-x}k\binom{k-x}k$) and $q_{m}(1)=q_m(0)=1$.

*Linear term.* Differentiating (only the $-2q_nq_{n-1}$ part survives at $x=0$):
$-2\bigl[q_n'(0)q_{n-1}(1)+q_n(0)q_{n-1}'(1)\bigr]$.  From the sum,
$q_m'(0)=\sum_{k\ge1}\binom mk\frac{(-1)^k}{k}=-H_m$, and $q_m'(1)=-q_m'(0)=H_m$.  So the
coefficient is $-2(-H_n+H_{n-1})=2/n$.

*Quadratic term* $=\lim_{x\to\infty}\Omega_n$.  By (1.2) the leading term of $r_m$ is the
$k=m$ term, $(-1)^m\frac{m!}{x^{m+1}}\cdot\frac{m!}{(-x)^{m+1}}=-\dfrac{(m!)^2}{x^{2m+2}}$,
while $q_n$ has leading coefficient $1/(n!)^2$.  Hence by (2.1)
$\lim\Omega_n=-((n-1)!)^2/(n!)^2=-1/n^2$. $\square$

At the moving point $x=x_n=\tfrac12-n$ one has $(x_n-n)^2+n^2=\tfrac14(20n^2-8n+1)$ and
$x_n^2=\tfrac14(2n-1)^2$, so
$$\boxed{\ \Omega_n(x_n)=-\frac{p(n)}{n^2(2n-1)^2}.\ }\tag{2.2}$$
*Zudilin's mysterious quadratic $p(n)=20n^2-8n+1$ is $4\bigl[(x-n)^2+n^2\bigr]$ at $x=x_n$.*

*(Verified symbolically as a polynomial identity in $x$ for $1\le n\le26$, `caso2fit.gp`.)*

---

## 3. Lemma D (skew Casoratian): where Zudilin's $q(n)$ comes from

Set $\Omega'_n(x):=p_{n+1}(x{-}1)q_{n-1}(x{+}1)-p_{n-1}(x{+}1)q_{n+1}(x{-}1)
-\Bigl(\frac2{(x-1)^2}-\frac2{x^2}\Bigr)q_{n+1}(x{-}1)q_{n-1}(x{+}1)$
and $W_n(x):=x^2(x-1)^2\,\Omega'_n(x)$.

> **Lemma D.** $W_n$ is a polynomial of degree $6$, namely
> $$n^2(n{+}1)^2W_n(x)=x^6-(4n{+}5)x^5+10(n{+}1)^2x^4-(2n{+}1)(2n{+}3)(4n{+}3)x^3$$
> $$\qquad+(2n{+}1)(2n{+}3)(4n^2{+}6n{+}1)x^2-2n(n{+}1)(4n^3{+}14n^2{+}14n{+}3)x+2n^2(n{+}1)^2(2n{+}3),$$
> and consequently
> $$\boxed{\ \Omega'_n(x_n)=\frac{q(n)}{4n^2(n+1)^2(2n-1)^2(2n+1)^2},\qquad
> W_n(x_n)=\frac{q(n)}{64\,n^2(n+1)^2}.\ }\tag{3.1}$$

*Proof.* Subtracting the two shift equations gives
$\Theta(x-1)-\Theta(x+1)=\frac2{(x-1)^2}-\frac2{x^2}$, whence
$$\Omega'_n=r_{n+1}(x{-}1)\,q_{n-1}(x{+}1)-r_{n-1}(x{+}1)\,q_{n+1}(x{-}1)
=O(x^{-6})-O(x^{-2n})O(x^{2n+2}),$$
so $\Omega'_n=O(x^2)$; its only poles are $x=0,1$, each of order $\le2$; hence
$\deg W_n\le6$.  Seven conditions determine it, and each is elementary:

1. **$x^6,x^5$** (asymptotics).  By (1.2), $r_{m}=-\frac1{(m+1)^2q_{m+1}}+O(x^{-2m-6})$
   (telescope the Casoratian of (1.1), which is $p_nq_{n-1}-q_np_{n-1}=1/n^2$), so
   $\Omega'_n=\frac{q_{n+1}(x-1)}{n^2q_n(x+1)}+O(x^{-2})$.  Writing $q_m=\tilde q_m(y)$,
   $y=x^2-x$, $\tilde q_m=\frac1{(m!)^2}(y^m+\sigma_my^{m-1}+\cdots)$, one gets
   $\frac{q_{n+1}(x-1)}{q_n(x+1)}=\frac{(x^2-3x+2)}{(n+1)^2}\Bigl(\frac{x^2-3x+2}{x^2+x}\Bigr)^{n}(1+O(x^{-2}))
   =\frac{x^2-(4n+3)x+O(1)}{(n+1)^2}$, and multiplying by $x^2(x-1)^2$ gives
   $w_6=1$, $w_5=-(4n+5)$ (after $n^2(n+1)^2$-normalisation).
2. **$W_n(0)=2q_{n+1}(-1)q_{n-1}(1)=2(2n+3)$**, using $q_m(-1)=2m+1$ (only $k=0,1$
   survive) and $q_m(1)=1$.  [Here $x^2(x-1)^2\bigl(\frac2{(x-1)^2}-\frac2{x^2}\bigr)=4x-2$.]
3. **$W_n(1)=-2q_{n+1}(0)q_{n-1}(2)=-2(2n-1)$**, using $q_m(2)=q_m(-1)=2m+1$.
4. **$W_n'(0)=-(2n+3)\bigl[4+\frac2n+\frac2{n+1}\bigr]$** and
   **$W_n'(1)=(2n-1)\bigl[-4+\frac2n+\frac2{n+1}\bigr]$**, using in addition
   $q_m'(-1)=-(2m+1)H_m$ and $q_m'(2)=-q_m'(-1)$, and $H_{n+1}-H_{n-1}=\frac1n+\frac1{n+1}$.
5. **$W_n(1/2)=\dfrac{8n^2+8n+7}{64n^2(n+1)^2}$.**  At $x=\frac12$ the elementary term
   vanishes and, by the $x\mapsto1-x$ symmetry, $\Omega'_n(\tfrac12)=
   p_{n+1}(-\tfrac12)q_{n-1}(-\tfrac12)-p_{n-1}(-\tfrac12)q_{n+1}(-\tfrac12)$, which by (1.1)
   equals $\frac{x^2-x+2n^2+2n+1}{(n+1)^2}\bigl(p_nq_{n-1}-q_np_{n-1}\bigr)
   =\frac{2n^2+2n+7/4}{n^2(n+1)^2}$ at $x=-\tfrac12$.

Seven independent conditions, so $W_n$ is as stated; (3.1) is then a polynomial identity
in $n$. $\square$

*(The closed form and (3.1) verified symbolically for $1\le n\le26$, `caso2fit.gp`;
the $q$-value closed forms of steps 2–4 in `qvals.gp`.)*

---

## 4. Proof of the bridge identity (B)

### 4.1 $Q_m=q_m(x_m)$ (Rivoal's identification — here proved)

Write $\hat q_k:=q_k(x_k)$, $\hat p_k:=p_k(x_k)$; note $x_k-1=x_{k+1}$, $x_k+1=x_{k-1}$.
Lemma A at index $n$ (variable $x_n$), Lemma A at index $n+1$ (variable $x_{n+1}=x_n-1$),
and Lemma D at $x_n$ are three **linear** equations for $(\hat p_{n-1},\hat p_n,\hat p_{n+1})$
with coefficient matrix
$$\begin{pmatrix}\hat q_n&\hat q_{n-1}&0\\ 0&\hat q_{n+1}&\hat q_n\\ -\hat q_{n+1}&0&\hat q_{n-1}\end{pmatrix},$$
whose determinant vanishes identically.  Its left null vector is
$(\hat q_{n+1},-\hat q_{n-1},\hat q_n)$, so the right-hand sides must satisfy the same
relation.  In that combination **all cubic terms $\hat q_{n+1}\hat q_n\hat q_{n-1}$ cancel**
(the elementary coefficients $\frac2{x_n^2}$, $\frac2{x_{n+1}^2}$ combine exactly with the
$\frac2{x_{n+1}^2}-\frac2{x_n^2}$ of Lemma D), leaving the *linear* three-term relation
$$\frac{p(n)}{n^2(2n-1)^2}\,\hat q_{n+1}
=\Omega'_n(x_n)\,\hat q_n+\frac{p(n+1)}{(n+1)^2(2n+1)^2}\,\hat q_{n-1}.$$
Multiplying by $4n^2(n+1)^2(2n-1)^2(2n+1)^2$ and using (3.1) turns this into **exactly**
Zudilin's recursion (1.3).  Since $\hat q_0=1=Q_0$ and $\hat q_1=q_1(-\tfrac12)=\tfrac74=Q_1$,
$$\hat q_m=Q_m\quad\text{for all }m\ge0. \qquad\square$$

*(This is Rivoal's "moving Padé specialisation"; the proof above is self-contained.
Verified exactly for $m\le160$, `verify_bridge.gp`.)*

### 4.2 Zudilin's Casoratian

From (1.3), $\mathrm{Cas}_{m+1}=-\frac{(2m-1)^2(2m)^2p(m+1)}{(2m+1)^2(2m+2)^2p(m)}\mathrm{Cas}_m$
with $\mathrm{Cas}_m:=P_mQ_{m-1}-P_{m-1}Q_m$ and $\mathrm{Cas}_1=\tfrac{13}8$.  All three
factors telescope, giving
$$\mathrm{Cas}_m=\frac{(-1)^{m-1}(20m^2-8m+1)}{8m^2(2m-1)^2}\tag{4.1}$$
— the Wronskian quoted in the $5{:}8$ paper, now proved.

### 4.3 The identity

Put $\hat P_m:=\frac{(-1)^m}{8}\hat p_m+G_m\hat q_m$ with $G_m=\sum_{j\le m}\frac{(-1)^{j-1}}{(2j-1)^2}$,
$g_m:=G_m-G_{m-1}=\frac{(-1)^{m-1}}{(2m-1)^2}$.  Then
$$\hat P_mQ_{m-1}-\hat P_{m-1}Q_m
=\frac{(-1)^m}{8}\Bigl[\hat p_m\hat q_{m-1}+\hat p_{m-1}\hat q_m-\frac2{x_m^2}\hat q_m\hat q_{m-1}\Bigr]
=\frac{(-1)^m}{8}\,\Omega_m(x_m),$$
because $g_m=-\frac{(-1)^m}{8}\cdot\frac2{x_m^2}$.  By (2.2) this equals
$\frac{(-1)^{m-1}p(m)}{8m^2(2m-1)^2}=\mathrm{Cas}_m$.  Hence $\phi_m:=\hat P_m-P_m$ satisfies
$\phi_mQ_{m-1}=\phi_{m-1}Q_m$ for all $m\ge1$, so $\phi_m/Q_m$ is constant $=\phi_0/Q_0=0$
($Q_m>0$ by Zudilin).  Therefore $P_m=\hat P_m$, which is (B). $\square$

---

## 5. The $2$-adic evaluation

Take $p=2$, $F=2$, $a_m=1-2m$, so $x_m=a_m/F$ with $2\nmid a_m$, $|x_m|_2=2$, $r=1$.

**(i) The target moves, explicitly.**  The identity $\Theta(x)+\Theta(x+1)=2/x^2$ holds in
$K$; all its terms, evaluated at any $x$ with $|x|_2>1$, tend to $0$ $2$-adically, so the
rearrangement implicit in expanding $\Theta(x+1)$ is legitimate and
$\Theta_2(x)+\Theta_2(x+1)=2/x^2$ for $|x|_2>1$.  Since $|x_m|_2=2$ for every $m$, iterating
from $x_0=\tfrac12$ with $\Theta_2(1/2)=8\zeta_2(2)$ (Beukers, Prop. 5.1(1)) gives
$$\boxed{\ \Theta_2(x_m)=(-1)^m\,8\bigl(\zeta_2(2)-G_m\bigr).\ }\tag{5.1}$$
*This is the whole content of the "moving target".*  `SOURCES_S18_ZUDILIN.md` §6.4 is right
that $\Theta_2(x_m)\not\to\Theta_2(1/2)$ — the drift is $2^{-2-v_2(m)}$ — but wrong to
conclude the passage fails: the drift is the elementary sequence $8(-1)^{m+1}G_m$, and (B)
says $P_m$ already carries $+G_mQ_m$.

**(ii) The bridge.**  Combining (B) with (5.1),
$$\frac{P_m}{Q_m}-\zeta_2(2)=\frac{(-1)^m}{8}\cdot\frac{p_m(x_m)-\Theta_2(x_m)q_m(x_m)}{q_m(x_m)}
=\frac{(-1)^m}{8}\cdot\frac{r_m(x_m)}{Q_m}.\tag{5.2}$$

**(iii) The two valuations, exactly.**  For $x=x_m$ every factor $x_m+j=\frac{2j+1-2m}2$ and
$1-x_m+j=\frac{2m+1+2j}2$ is a half-odd-integer, so
$\bigl[\frac k{x_m}\bigr]\bigl[\frac k{1-x_m}\bigr]=(k!)^2\,2^{2k+2}/\text{odd}$.  Hence in
(1.2) the $k$-th term has
$v_2=v_2\binom km+2\bigl(k+v_2(k!)\bigr)+2\ \ge\ 2\bigl(k+v_2(k!)\bigr)+2$, with equality at
$k=m$.  Since $k\mapsto k+v_2(k!)=2k-s_2(k)$ is **strictly increasing**, the $k=m$ term is
the unique minimiser and
$$v_2\bigl(r_m(x_m)\bigr)=2\bigl(m+v_2(m!)\bigr)+2=4m+2-2s_2(m).\tag{5.3}$$
The same argument applied to $q_m(x)=\sum_k\binom mk\binom{-x}k\binom{k-x}k$ at $x=x_m$
(where $\binom{m-1/2}k$ and $\binom{k+m-1/2}k$ each have $v_2=-k-v_2(k!)$) gives
$$v_2(Q_m)=-2\bigl(m+v_2(m!)\bigr)=-4m+2s_2(m),\tag{5.4}$$
the denominator law of `SOURCES_S18_ZUDILIN.md` §5 — also now proved rather than measured.

**(iv) Conclusion.**  By (5.2)–(5.4),
$$v_2\Bigl(\zeta_2(2)-\frac{P_m}{Q_m}\Bigr)=\bigl(4m+2-2s_2(m)\bigr)-3+\bigl(4m-2s_2(m)\bigr)
=8m-1-4s_2(m)\ \longrightarrow\ \infty .$$
$\blacksquare$

*(Even without (5.3) the theorem follows from Beukers' Prop. 7.1(4) alone —
$|p_n(a/F)-\Theta_p(a/F)q_n(a/F)|_p\le p^2n^2p^{-2n(r+1/(p-1))}$, **uniformly in $a$**, so
one may take $n=m$, $a=1-2m$ — giving $v_2\ge 8m-5-2s_2(m)-2\log_2m$.  (5.3) sharpens the
uniform bound to an equality.  All of $m\le60$ checked at precision $2^{700}$,
`final_check.gp`.)*

---

## 6. What this buys

* **$\xi_2^{\rm Zud}=\zeta_2(2)$ is a theorem.**  With $\xi_2^{\mathbf E}=\frac12\zeta_2(2)$
  (Theorem F, `EULER_CRITERION.md` §4.1, modulo its hypothesis (b)), the ratio
  $\xi_2^{\rm Zud}=2\xi_2^{\mathbf E}$ is proved, and with it the *identification* half of
  the common-period lemma of the $5{:}8$ paper.
* **The exact Zudilin tail** $v_2(G_2-P_m/Q_m)=8m-1-4s_2(m)$ — `\eqref{eq:Ztail}` of the $5{:}8$
  tex, previously "given the two limits are identified" — is proved outright, including
  $v_2(Q_m)=-4m+2s_2(m)$.  The $\mathbf E$-side tail `\eqref{eq:Etail}` of that paper is untouched.
* **Independent of modularity.**  `SOURCES_S18_ZUDILIN.md` §5 shows Zudilin's row has no
  modular source (order-$4$ operator, non-integral mirror map).  The proof above needs
  none: the $2$-adic period arrives through Beukers' $p$-adic Padé theory, not through an
  overconvergent Eisenstein family.  Route 2 of the task (Hadamard/quadratic transform to
  the $\mathbf E$ row) is therefore **unnecessary**, and the negative verdict of §5 stands.
* **Lean.**  `RequestProject/Limits.lean` proves $\mathrm{GZ2}$ exists from strictly
  increasing $v_2$; (5.2)+(5.3) now give its *value* and the exact rate, so the
  `wZfun` bound can be replaced by the closed form $8m-1-4s_2(m)$.

---

## 7. Remaining dependencies (all published, none numerical)

1. Beukers Prop. 5.1(1): $\Theta_2(1/2)=8\zeta_2(2)$ (his proof, via
   $\Theta_2(1/2)=-8\bigl(H_2(2,1,4)+H_2(2,3,4)\bigr)$ and the Appendix).  Independently
   re-verified here to $2^{261}$.
2. Beukers Prop. 6.1: the closed form (1.2) of the remainder (a Zeilberger certificate,
   printed in his paper).  Used for the leading coefficient in Lemma A/D and for (5.3).
3. Zudilin Thm. 1: the recursion (1.3) with the stated initial data (definition of the row).

Everything else — Lemmas A, D, §4, §5 — is proved here and machine-checked exactly.

---

## 8. Corrections to `SOURCES_S18_ZUDILIN.md`

* **§6.3** "the numerators are *not* the specialised Padé numerators … a bounded,
  never-vanishing correction": correct as far as it goes, but the correction is identified
  in (B): $P_m-\frac{(-1)^m}{8}p_m(x_m)=G_mQ_m$, the $m$-th partial sum of Catalan's series
  times $Q_m$.  (The comparison made there, $P_m-p_m(x_m)$, omits both the sign $(-1)^m$ and
  the factor $\frac18$, which is why nothing was visible.)
* **§6.4** "The moving-point passage is invalid as written": the *stated reason* (the target
  drifts by $2^{-2-v_2(m)}$) is correct, and so is the observation that
  $p_m(x_m)/Q_m$ does not converge.  The *conclusion* that the passage cannot be repaired is
  wrong: §5(i) above carries out the repair, and the "explicit elementary rational" it asked
  for is $8(-1)^{m+1}G_m$.  The "$(-1)^m$ and the exact factor $2$ (not $\frac12$)" are
  explained by (5.1) and by $\Theta_2(1/2)=8\zeta_2(2)$ (not $-8$).
* **§6.5 / §0.6** "$\xi_2^{\rm Zud}=\zeta_2(2)$ … remains **numerical**. **Open.**":
  **closed**; see §0.
* **§5** "$v_2(Q_m)=-4m+2s_2(m)$ (`zud_check.gp`, $m\le14$)": proved in (5.4).
* Unaffected: the verdict that Zudilin's row has no modular source (§5), the $s_{18}$
  analysis (§1–§4), and the normalisation discussion (§6.1).
