# Paşol–Zudilin's magnetic forms, the Shimura–Borcherds dictionary, and Cooper's congruence

*Fable, 2026-09-02.  Scripts `30_`–`39_` and their `.log` files in this directory.
All arithmetic exact (PARI/GP 2.15.4).  Tags: **[proved]** = proof written out below,
**[verified, range]** = exact computation in the stated range, **[conjectural]**,
**[refuted]**.*

---

## 0.  Summary of what changed

Four things came out of this that bear directly on `consolidation/COOPER_CONGRUENCE.md`:

1. **The "negative control" in §3 of the ledger is wrong.**  Paşol–Zudilin's level-one
   magnetic forms *do* satisfy Cooper's congruence $(S)$ — one just has to use the
   right character.  With $\psi=\chi_{-3}$ for $F_{4a}=\Delta/E_4^2$ (pole at $\rho$,
   disc $-3$) and $\psi=\chi_{-4}$ for $F_{4b}=E_4\Delta/E_6^2$ (pole at $i$, disc $-4$),
   $(S)$ holds at **every** prime where the half-integral input is $p$-integral.
   The ledger tested $\psi=\mathbf 1$, which is the wrong character for these forms.

2. **There is an exact $T_{p^2}$-eigen-identity behind PZ's Theorem 1**, not merely a
   congruence, and it comes with a *tower*.  Running PZ's Lemma-2 uniqueness argument
   once and keeping the value rather than the valuation gives
   $$f\mid T_{p^2} \;=\; \psi(p)\,p\,f \;+\; p^{3}\,c\,g_{m_0p^2},\qquad
     g_{m_0p^{2r}}\mid T_{p^2} \;=\; g_{m_0p^{2r-2}} + p^{3} g_{m_0p^{2r+2}}\ (r\ge1).$$
   The second (tower) relation is what PZ never need and what produces the *extra*
   cancellation $A(p^2n)\equiv \lambda A(pn) \bmod p^{6}$.  Together they give exactly
   Cooper's `eq:magnetic` hierarchy $\Phi\mid U_{p^n}\equiv(\psi(p)p)^n\Phi \bmod p^{n+2}$.

3. **The answer to "what gives $n^3\mid a(n^2)$" is $\lambda\equiv p\pmod{p^3}$**, and under
   that hypothesis one gets the *exact* formula $a(p^{2j}m^2)=p^{3j}c\,g_{m_0p^{2j}}(m^2)$.
   Since a single-index principal part $c\,q^{-m_0}$ forces $\lambda=\chi_{-m_0}(p)\,p$, the
   level-one forms have $p^{3j}\mid a(p^{2j}m^2)$ at *exactly the split primes*
   $\chi_{-m_0}(p)=+1$, and $v_p$ exactly $j+v_p(a(m^2))$ at the inert ones.

4. **Cooper's Conjecture 4.1 is a theorem at level one.**  §2.7 proves, for every
   $f=c\,q^{-m_0}+O(q)\in S^{!,+}_{5/2}(\Gamma_0(4))$ and every prime $p$ with $c\in\mathbf Z_p$,
   $$a(p^{2a}m^2)+\bigl(1-\psi(p)\bigr)\sum_{j=0}^{a-1}p^{\,a-j}a(p^{2j}m^2)=p^{3a}\,c\,g_{m_0p^{2a}}(m^2)
     \qquad(p\nmid m,\ \psi=\chi_{-m_0}),$$
   hence $n^2\mid\beta_\psi(n)$ for all $n$ — Conjecture 4.1 for $\Psi(f)$.  The
   $(1-\psi(p))$ correction terms, which are exactly what the $\psi$-twisted Möbius
   inversion supplies, are what repairs the inert primes.  This covers **every entry of
   PZ's Table 1**, and explains the companion strand's numerical sweep over that table.
   It does *not* cover Cooper's rows: at level $4N$ the two inputs of the proof
   ($\dim S_4=0$, and the integral basis $\{g_M\}$) both fail, and the character comes out
   as $\chi_{-3}$ rather than $\mathbf 1$ — which is precisely why the correct level-$4N$
   model must be the **twisted** (genus-character) lift.

---

## 1.  The dictionary (Task 1, step 0)  **[proved]**

Paşol–Zudilin's lift (their (5)/(6)) for $k$ even (so $D_k=1$ and the Kronecker symbol is $1$):
$$A(n)=\sum_{d\mid n} d^{\,k-1}\,a(n^2/d^2).$$
Put $A'(n)=A(n)/n^{k-1}$.  Substituting $e=n/d$,
$$A'(n)=\sum_{e\mid n}\frac{a(e^2)}{e^{k-1}},$$
so $A'=\mathbf 1\star\alpha$ with $\alpha(e)=a(e^2)/e^{k-1}$, and Möbius inversion gives
$$\boxed{\ \beta(m):=(A'\star\mu)(m)=\frac{a(m^2)}{m^{k-1}},\qquad a(m^2)=m^{k-1}\beta(m).\ }$$
For $k=2$: $\beta(m)=a(m^2)/m$, so PZ's divisibility $n\mid a(n^2)$ **is** $\beta\in\mathbf Z$
(magnetism), and Cooper's $n^2\mid\beta(n)$ **is** $n^3\mid a(n^2)$.  Confirmed exactly.

**The $k=3$, $D=-3$ analogue (derived here).**  PZ's (5) reads
$A(n)=\sum_{d\mid n}\bigl(\tfrac{-3}{d}\bigr)d^{2}a(3n^2/d^2)$.  With $\chi=\chi_{-3}$ and
$A''(n)=A(n)/n^2$ one gets $A''=\chi\star b$, $b(e)=a(3e^2)/e^2$, hence
$$\boxed{\ \beta_3(m):=(A''\star(\mu\chi_{-3}))(m)=\frac{a(3m^2)}{m^{2}},\qquad a(3m^2)=m^2\beta_3(m).\ }$$
PZ's Theorem 2 ("doubly magnetic", $n^2\mid A(n)$) is $\beta_3\in\mathbf Z$.
Verified numerically: `30_calib.log` reports $A(n)/n^2\in\mathbf Z$ and $\beta_3\in\mathbf Z$ for $n\le300$.

### Calibration against PZ  **[verified]**

`33_halfint.gp` builds $g_0,g_1,g_2$ from scratch and reproduces every printed
coefficient of PZ §4 (`33_halfint.log`).  Lemma 1's combinations give

| | | |
|---|---|---|
| $f_{4a}=\tfrac78 g_0+\tfrac1{768}(g_1-g_2)$ | $=\tfrac1{64}q^{-3}+q-506q^4+\cdots$ | $64f_{4a}\in\mathbf Z[[q]]$ ✓ |
| $f_{4b}=\tfrac{19}{18}g_0-\tfrac5{648}g_1-\tfrac1{648}g_2$ | $=-\tfrac1{108}q^{-4}+q+1222q^4+\cdots$ | $108f_{4b}\in\mathbf Z[[q]]$ ✓ |

**Normalisation note.**  The plain-text dump `pz.txt` mangles the leading coefficients of
Lemma 1 (the fractions $\tfrac1{64},\tfrac1{108}$ are printed on the denominator line and
read as if they were part of the linear combination).  The correct expansions are the ones
above: $f_{4a}$ starts $\tfrac1{64}q^{-3}$, **not** $q^{-3}$.  With that, $a(1)=1$ and
$a(4)=-506$ exactly as the task expected, and the lift formula reproduces
$\Delta/E_4^2$ and $E_4\Delta/E_6^2$ coefficient by coefficient for $n\le50$ **[verified]**.

Dictionary check $a(m^2)=m\beta(m)$: **exact agreement for all $m\le20$**, both forms
(`33_halfint.log`).  No normalisation constant had to be inserted.

### Square-index coefficients

| $m$ | $a(m^2)$, $F_{4a}$ | $a(m^2)$, $F_{4b}$ | $a(3m^2)$, $F_6$ |
|---|---|---|---|
| 1 | $1$ | $1$ | $1$ |
| 2 | $-506$ | $1222$ | $-1244$ |
| 3 | $180249$ | $1009689$ | $714996$ |
| 4 | $-56363984$ | $731061424$ | $-307867520$ |
| 5 | $16415391865$ | $493402005625$ | $114237828175$ |
| 6 | $-4574618693994$ | $318792565630806$ | $-38689288815024$ |
| 7 | $1237162549543249$ | $199937823789551953$ | $12317318339105143$ |

### Exact power of $m$ dividing $a(|D|m^2)$, $m\le300$  **[verified]**

Writing $m^{t}\,\|\,a(|D|m^2)$ (`30_calib.log`):

| form | $t=1$ | $t=2$ | $t=3$ | $t=4$ | $t\ge5$ |
|---|---|---|---|---|---|
| $F_{4a}$ | 220 | 44 | 35 | 0 | 0 |
| $F_{4b}$ | 212 | 47 | 39 | 1 | 0 |
| $F_6$ ($t$ counted against $m^2$) | 0 | 234 | 16 | 14 | 32 |

$\gcd_{2\le m\le300}\beta(m)=1$ in all three cases, so PZ's $m\mid a(m^2)$ is sharp
— but not uniformly.  The $m$ with $t\ge2$ are governed by the *split* primes dividing $m$: $t\ge2$ needs
$v_p(\beta(m))\ge e_p$ at every $p^{e_p}\|m$, and at an inert $p$ one has $v_p(\beta(m))=0$.
(The exact set is not clean at ramified or bad $p$ — e.g. for $F_{4a}$, $m=3,21,39,\dots$
have $t=1$ because $\chi_{-3}(3)=0$ — and there are a couple of numerical coincidences
$m=209$ for $F_{4a}$, $m=152,187$ for $F_{4b}$.)  The clean statement is the prime-power
one (`30_calib.log`):
$$v_p\bigl(\beta(p^e)\bigr)=\begin{cases} 2e & \chi_{D_0}(p)=+1\ \ (\Rightarrow v_p(a(p^{2e}))=3e)\\[2pt]
0 & \chi_{D_0}(p)=-1,\ p>3 \ \ (\Rightarrow v_p(a(p^{2e}))=e)\end{cases}$$
with $D_0=-3$ for $F_{4a}$, $D_0=-4$ for $F_{4b}$; e.g. $F_{4a}$, $p=7$: $v_7(\beta(7^e))=2,4$;
$p=5$: $0,0,0$.  (The primes $p=2,3$ dividing the denominator $64$, resp. $108$, behave
differently and are excluded.)

---

## 2.  The mechanism (Task 3)

### 2.1  PZ's Lemma 2, in my own words  **[exposition]**

*Hypotheses.*  $k\ge1$ with $\dim S_{2k}(\mathrm{SL}_2(\mathbf Z))=0$; $p$ prime;
$f\in M^{!,+}_{k+1/2}$ with $p$-integral coefficients and $p^2>-\operatorname{ord}_q f$.

*Statement.*  $f\mid T_{p^2}^{\,n}\equiv 0 \pmod{p^{(k-1)n}}$.

*Proof.*  With $T_{p^2}=U_{p^2}+p^{k-1}\chi_p+p^{2k-1}V_{p^2}$ and the relations
$V_{p^2}\chi_p=\chi_p U_{p^2}=0$, $V_{p^2}U_{p^2}=\mathrm{id}$, one can write
$$T_{p^2}^{\,n}=\sum_{\substack{a+b+c=n\\ r\le\min(a,c)}}\alpha_{a,b,c,r}\;
   p^{(2k-1)c+(k-1)b}\;U_{p^2}^{\,a-r}\chi_p^{\,b}V_{p^2}^{\,c-r},\qquad \alpha_{\bullet}\in\mathbf Z .$$
A single surviving $U_{p^2}$ annihilates the principal part (because $p^2>-\operatorname{ord}_q f$),
so only $a=r\le c$ contributes to the principal part; for those,
$(2k-1)c+(k-1)b\ge(k-1)(2c+b)\ge(k-1)n$.  Since $\dim S_{2k}=0$, a weakly holomorphic plus
form is determined by its principal part, and the plus space has the integral basis
$\{g_M=q^{-M}+O(q)\}$ (Borcherds/Zagier), so $f\mid T_{p^2}^{\,n}=p^{(k-1)n}g$ with $g$ integral.

*Transport.*  $\Psi(f)\mid T_p^{\,n}=\Psi\bigl(f\mid T_{p^2}^{\,n}\bigr)$, whence
$F\mid U_p^{\,n}\equiv0 \bmod p^{(k-1)n}$, i.e. $p^n\mid m\Rightarrow p^{(k-1)n}\mid A(m)$.
For $k=2$ this is magnetism ($n\mid A(n)$), equivalently $n\mid a(n^2)$ via §1.

### 2.2  The sharpening: an exact eigen-identity plus a tower  **[proved]**

Keep the *value*, not the valuation, at $n=1$.  Let $f=c\,q^{-m_0}+O(q)\in S^{!,+}_{5/2}$
($k=2$, level $4$), $D_0=-m_0$ a discriminant, $\psi=\chi_{D_0}$, $g_M=q^{-M}+O(q)$ the
integral basis.  For $p$ odd (or $p=2$ with $T'_4=K^+\!\circ T_4$):

* $(f\mid U_{p^2})(n)=a(p^2n)=0$ for $n<0$, since $p^2 m_0>m_0$;
* $(f\mid\chi_p)$ contributes $p\,\chi_p(-m_0)\,c\,q^{-m_0}=\psi(p)\,p\,c\,q^{-m_0}$;
* $(f\mid V_{p^2})$ contributes $p^{3}c\,q^{-m_0p^2}$;
* the constant term stays $0$.

Hence, by uniqueness,
$$\textbf{(T1)}\qquad f\mid T_{p^2}\;=\;\psi(p)\,p\,f\;+\;p^{3}\,c\,g_{m_0p^2}.$$
Applying the *same* computation to $g_M$, $M=m_0p^{2r}$ ($r\ge1$) — where now $p^2>M$ **fails**,
so $U_{p^2}$ *does* reach the pole and contributes $q^{-M/p^2}=q^{-m_0p^{2r-2}}$, while the
$\chi_p$ term dies because $\chi_p(-M)=\bigl(\tfrac{-m_0p^{2r}}{p}\bigr)=0$ — gives the
relation PZ never need:
$$\textbf{(T2)}\qquad g_{m_0p^{2r}}\mid T_{p^2}\;=\;g_{m_0p^{2r-2}}\;+\;p^{3}\,g_{m_0p^{2r+2}},\qquad r\ge1 .$$

**Numerical verification of (T1)** (`34_hecke.log`): for $p=2,3,5,7,11,13,17$ the residual
$E=f\mid T_{p^2}-\psi(p)p\,f$ has principal part *exactly* $p^3c\,q^{-m_0p^2}$, constant term $0$,
and $64E/p^3$ (resp. $-108E/p^3$) integral over the whole tested range.  The single exception
is $f_{4b}$ at $p=2$: there $p^2=4=m_0$, so the hypothesis $p^2>m_0$ fails, $U_{p^2}$ reaches the
pole and contributes $a(-4)q^{-1}$, and the residual acquires a term at $n=-1\not\equiv0,1\ (4)$
— which is exactly PZ's remark that $T_4$ does not preserve the plus space, in this instance.  A control shows $\lambda=0$ and
$\lambda=-\psi(p)p$ both leave a nonzero $q^{-m_0}$ residual, so $\lambda$ is unique.

### 2.3  Consequences  **[proved from (T1)+(T2)], [verified]**

Write $\Phi_r=c\,\Psi(g_{m_0p^{2r}})$, so $\Phi_0=F$, all $p$-integral.  Applying $\Psi$ and
$T_p=U_p+p^{3}V_p$ (weight $4$):
$$F\mid U_p=\lambda F+p^3\Phi_1-p^3 F\mid V_p,\qquad
  \Phi_r\mid U_p=\Phi_{r-1}+p^3\Phi_{r+1}-p^3\Phi_r\mid V_p\ (r\ge1),$$
with $\lambda=\psi(p)p$.  Composing and using $(X\mid V_p)\mid U_p=X$, the two $-p^3F$ terms
**cancel**:
$$\textbf{(C1)}\ \ A(pn)\equiv\lambda A(n)\ (p^{3}),\qquad
  \textbf{(C2)}\ \ A(p^2n)\equiv\lambda A(pn)\ (p^{6}),$$
$$\textbf{(C3)}\ \ A(p^rn)\equiv\lambda^{r}A(n)\ (p^{\,r+2})\quad(r\ge1).$$
(C3) is *verbatim* Cooper's `eq:magnetic` (i).  Its proof from (C1),(C2) is one line:
$A(p^rn)=\lambda^{r-1}A(pn)+O(p^6)$ and $\lambda^{r-1}A(pn)=\lambda^rA(n)+\lambda^{r-1}O(p^3)$,
so the error has $v_p\ge\min(6,r+2)=r+2$.

Without (T2) the cancellation in (C2) does not happen and one only gets $v_p\ge3$;
(C2) is therefore the load-bearing new fact.

**Verification** (`35_tower.log`, `39_c3.log`), $n\le1250$:

| form | (C1) mod $p^{2k-1}$ | (C2) mod $p^{2(2k-1)}$ | (C3) mod $p^{r+2}$ / $p^{2r+3}$, $r\le5$ |
|---|---|---|---|
| $F_{4a}$, $\lambda=\chi_{-3}(p)p$ | PASS $3\le p\le19$ | PASS | PASS |
| $F_{4b}$, $\lambda=\chi_{-4}(p)p$ | PASS $p=2$, $5\le p\le19$ | PASS | PASS |
| $F_6$, $\lambda=p^2$ | PASS $5\le p\le19$ | PASS | PASS |

The exceptions are precisely the primes where the *input* is not $p$-integral and the
character does not save it: $p=2$ for $f_{4a}$ ($64=2^6$), $p=3$ for $f_{4b}$ ($108=4\cdot27$),
$p=2,3$ for $f_6$ ($384=2^7\cdot3$).  (PZ handle those by hand; my sharpening does not.)

### 2.4  Cooper's congruence $(S)$ for the level-one forms  **[verified, refutes the ledger's §3]**

$(S)$: $c'(pm)\equiv\psi(p)c'(m)\bmod p^2$ with $c'(m)=A(m)/m$.
`31_congr.log`, tested for all $m\le800/p$:

| form | $\psi$ | result |
|---|---|---|
| $F_{4a}=\Delta/E_4^2$ | $\chi_{-3}$ | **PASS for every $3\le p\le59$** (fails only at $p=2$) |
| $F_{4a}$ | $\mathbf 1$ | fails at $p=2,3,5,11,17,23,29,41,47,53,59$ — i.e. at $\chi_{-3}(p)\ne1$ |
| $F_{4b}=E_4\Delta/E_6^2$ | $\chi_{-4}$ | **PASS for $p=2$ and every $5\le p\le59$** (fails only at $p=3$) |
| $F_{4b}$ | $\mathbf 1$ | fails at $\chi_{-4}(p)=-1$ |
| $F_6=E_6\Delta/E_4^3$ | — ($c''(pm)\equiv c''(m)\bmod p^{3}$, $c''=A/n^2$) | **PASS for every $5\le p\le59$** (fails at $p=2,3$) |

So the ledger's "**Negative control.** $(S)$ fails for Paşol–Zudilin's level-one magnetic forms"
is **[refuted]** as stated: it fails only because $\psi=\mathbf1$ was used.  With the
character of the *pole discriminant* it holds everywhere the input is $p$-integral.

### 2.5  The strong divisibility, and the answer to the Task-3 question

*(§2.7 below upgrades this to an unconditional theorem at level one, valid at split, inert
and ramified primes alike, once one works with the $\psi$-twisted $\beta_\psi$ rather than
with $a(m^2)/m$.)*

Set $\alpha_j=a(p^{2j}m^2)$ with $p\nmid m$, $\eta_{r,j}=c\,g_{m_0p^{2r}}(p^{2j}m^2)$.
The Hecke identity $(f\mid T_{p^2})(n)=a(p^2n)+p\bigl(\tfrac np\bigr)a(n)+p^3a(n/p^2)$ at
$n=p^{2j}m^2$ (note $\chi_p=0$ once $j\ge1$) gives, from (T1),(T2),
$$\alpha_1=(\lambda-p)\alpha_0+p^3\eta_{1,0},\qquad
  \alpha_{j+1}=\lambda\alpha_j-p^3\alpha_{j-1}+p^3\eta_{1,j}\ (j\ge1),$$
$$\eta_{r,1}=\eta_{r-1,0}-p\,\eta_{r,0}+p^3\eta_{r+1,0},\qquad
  \eta_{r,j+1}=\eta_{r-1,j}-p^3\eta_{r,j-1}+p^3\eta_{r+1,j}\ (j\ge1).$$

> **Proposition (the precise answer).**  Assume $f\mid T_{p^2}=\lambda f+p^{3}h_1$ with a
> $p$-integral tower $h_r\mid T_{p^2}=h_{r-1}+p^3h_{r+1}$, $h_0=f$.  Then
> $$v_p\bigl(a(p^{2j}m^2)\bigr)\;\ge\;3j+v_p\bigl(a(m^2)\bigr)\quad\text{for all }j\ge0,\ p\nmid m
> \qquad\Longleftrightarrow\qquad \boxed{\lambda\equiv p \pmod{p^{3}}}.$$
> Moreover under $\lambda=p$ one has the **exact** formula
> $$a(p^{2j}m^2)\;=\;p^{3j}\;h_j(m^2).$$

*Proof of $\Leftarrow$ (the case $\lambda=p$; the general case differs by $O(\lambda-p)$).*
$j=1$: $\alpha_1=0\cdot\alpha_0+p^3\eta_{1,0}=p^3\eta_{1,0}$.
$j=2$: substitute $\eta_{1,1}=\eta_{0,0}-p\eta_{1,0}+p^3\eta_{2,0}=\alpha_0-p\eta_{1,0}+p^3\eta_{2,0}$
into $\alpha_2=p\alpha_1-p^3\alpha_0+p^3\eta_{1,1}$; the $\pm p^3\alpha_0$ cancel and
$\alpha_2=p^4\eta_{1,0}-p^4\eta_{1,0}+p^6\eta_{2,0}=p^6\eta_{2,0}$.
$j=3$: $\eta_{1,2}=p^3\eta_{2,1}$ and $\eta_{2,1}=\eta_{1,0}-p\eta_{2,0}+p^3\eta_{3,0}$ give
$\alpha_3=p\cdot p^6\eta_{2,0}-p^3\cdot p^3\eta_{1,0}+p^6(\eta_{1,0}-p\eta_{2,0}+p^3\eta_{3,0})=p^9\eta_{3,0}$.
The pattern continues by the evident double induction on $(r,j)$. $\square$

*Proof of $\Rightarrow$.*  $\alpha_1=(\lambda-p)\alpha_0+p^3\eta_{1,0}$; pick $m$ with
$v_p(a(m^2))=0$ (e.g. $m=1$, $a(1)=1$).  Then $v_p(\alpha_1)\ge3$ forces $v_p(\lambda-p)\ge3$. $\square$

**Consequences.**
* Since $\lambda=\psi(p)p$ with $\psi$ quadratic, $\lambda\equiv p\ (p^3)$ $\iff$ $\psi(p)=+1$
  (for $p\ge3$).  So the level-one forms have the $n^3$-divisibility **exactly at the split
  primes**, and at inert $p$ one has $v_p(a(p^{2j}m^2))=j+v_p(a(m^2))$ exactly.
  **[verified]** (`35_tower.log`, $n\le1250$): $v_p(\beta(n))\ge 2v_p(n)$ PASSES at
  $p=7,13,19$ for $F_{4a}$ ($\chi_{-3}(p)=1$) and at $p=5,13,17,29$ for $F_{4b}$
  ($\chi_{-4}(p)=1$), with minimal excess $0$ (sharp); it FAILS at every inert prime.
* For **all** $n$ one therefore needs $\psi\equiv\mathbf1$, i.e. $\lambda_p=p$ for every $p$.

### 2.6  The theorem that would have to be proved for Cooper

> **Target Theorem.**  Let $N\in\{7,10\}$ ($\psi=\mathbf1$) and let $\Phi_s$ be Cooper's
> weight-4 source.  Suppose $\Phi_s=\Psi(f)$ for some $f\in S^{!,+}_{5/2}$ at level $4N$
> with $p$-integral coefficients, and suppose the residual tower
> $h_1:=p^{-3}(f\mid T_{p^2}-p f)$, $h_{r+1}:=p^{-3}(h_r\mid T_{p^2}-h_{r-1})$ consists of
> $p$-integral forms.  Then Conjecture 4.1 ($n^2\mid\beta(n)$, i.e. $n^3\mid a(n^2)$) holds
> at $p$ **iff** $f\mid T_{p^2}\equiv p\,f \pmod{p^3}$, and then $a(p^{2j}m^2)=p^{3j}h_j(m^2)$.

Ingredients and their status:

| ingredient | can it be checked numerically? |
|---|---|
| $\lambda_p=p$ for all $p$ (equivalently $\Phi_s\mid T_p\equiv p\Phi_s \bmod p^3$) | **yes**, and it is *already verified* — it is exactly `(S)` with $\psi=\mathbf1$, checked to $p\le199$ in the ledger and re-checked here for $p\le59$, $m\le400$ (`/tmp` scratch run reproduced in §3 below) |
| $p$-integrality of $f$ | needs $f$; see Task 4 |
| the tower $h_r$ being $p$-integral | needs $f$ and the level-$4N$ analogue of "integral basis $g_M$" (Borcherds' Prop. 2 at level $4N$ is *not* available in the same form — this is the real gap) |
| $\dim$-type uniqueness ("principal part determines the form") | **false** at level $28$: $\dim S^{+}_{5/2}(\Gamma_0(28))=1$ (PARI `mfkohnenbasis`), and there are $6$ cusps, so principal parts at all cusps are needed |

The last two rows are where a level-$4N$ proof would have to do genuinely new work.

### 2.7  A THEOREM: Conjecture 4.1 holds for every level-one strongly magnetic weight-4 form

The double induction of §2.5 in fact closes *unconditionally* at level one, for every
prime, split or inert or ramified.  Fix an odd prime $p$ (or $p=2$ with $T'_4=K^+\!\circ T_4$)
and $m$ with $p\nmid m$, and write
$$x_{r,j}:=c\,g_{m_0p^{2r}}\bigl(p^{2j}m^2\bigr)\in\mathbf Z_p\qquad (r,j\ge0),\qquad x_{0,j}=a(p^{2j}m^2),$$
with the convention $x_{r,-1}=0$.  Reading off the $n=p^{2j}m^2$ coefficient of (T1) and (T2),
and using $\chi_p(p^{2j}m^2)=[\,j=0\,]$, gives
$$x_{0,j+1}+p[\,j{=}0\,]x_{0,0}+p^3x_{0,j-1}=\psi(p)p\,x_{0,j}+p^3x_{1,j},$$
$$x_{r,j+1}+p[\,j{=}0\,]x_{r,0}+p^3x_{r,j-1}=x_{r-1,j}+p^3x_{r+1,j}\qquad(r\ge1).$$

> **Lemma.**  For $r\ge1$, $k\ge1$: $\;x_{r,k}-x_{r-1,k-1}=p^{3k}x_{r+k,0}-p^{3k-2}x_{r+k-1,0}$.
>
> *Proof.*  $k=1$: the $r$-recursion at $j=0$ gives $x_{r,1}=x_{r-1,0}-p\,x_{r,0}+p^3x_{r+1,0}$,
> i.e. $x_{r,1}-x_{r-1,0}=p^{3}x_{r+1,0}-p\,x_{r,0}$.  For $k\ge1$ the $r$-recursion at $j=k$ gives
> $x_{r,k+1}-x_{r-1,k}=p^3\bigl(x_{r+1,k}-x_{r,k-1}\bigr)$, so the claim propagates from $k$ to $k+1$
> with $r\mapsto r+1$.  $\square$

> **Theorem.**  Let $f=c\,q^{-m_0}+O(q)\in S^{!,+}_{5/2}(\Gamma_0(4))$, $-m_0$ a discriminant,
> $\psi=\chi_{-m_0}$, and let $p$ be a prime with $c\in\mathbf Z_p$ (and $p$ odd, or $p=2$ with $T'_4$).
> Then for every $a\ge0$ and every $m$ with $p\nmid m$,
> $$\boxed{\;a\bigl(p^{2a}m^2\bigr)\;+\;\bigl(1-\psi(p)\bigr)\sum_{j=0}^{a-1}p^{\,a-j}\,a\bigl(p^{2j}m^2\bigr)\;=\;p^{3a}\,c\;g_{m_0p^{2a}}\bigl(m^2\bigr).\;}$$
> Consequently $v_p(\beta_\psi(n))\ge2\,v_p(n)$ for all $n$, where $\beta_\psi=(A/n)\star(\mu\psi)$;
> and hence, over all $p$, $\;n^2\mid\beta_\psi(n)$ — **Cooper's Conjecture 4.1 for $\Psi(f)$.**

*Proof.*  Write $S_a:=x_{0,a}+(1-\psi(p))\sum_{j=0}^{a-1}p^{a-j}x_{0,j}$; the claim is $S_a=p^{3a}x_{a,0}$.
Separating the top term, $S_a=x_{0,a}-\psi(p)p\,x_{0,a-1}+p\,S_{a-1}$ for $a\ge1$, with $S_0=x_{0,0}$.
For $a=1$ the $r=0$ recursion at $j=0$ gives $x_{0,1}=(\psi(p)p-p)x_{0,0}+p^3x_{1,0}$, so
$S_1=p^3x_{1,0}$.  For $a\ge2$ the $r=0$ recursion at $j=a-1$ gives
$x_{0,a}-\psi(p)p\,x_{0,a-1}=-p^3x_{0,a-2}+p^3x_{1,a-1}$, whence, by the inductive hypothesis,
$$S_a=p\cdot p^{3(a-1)}x_{a-1,0}+p^3\bigl(x_{1,a-1}-x_{0,a-2}\bigr)
     =p^{3a-2}x_{a-1,0}+p^3\bigl(p^{3a-3}x_{a,0}-p^{3a-5}x_{a-1,0}\bigr)=p^{3a}x_{a,0},$$
using the Lemma with $r=1$, $k=a-1$.  For the consequence: $A/n=\mathbf 1\star\alpha$ with
$\alpha(e)=a(e^2)/e$ (§1), so $\beta_\psi=\alpha\star\nu$ with $\nu=\mathbf 1\star\mu\psi$,
$\nu(e)=\prod_{p\mid e}(1-\psi(p))\in\mathbf Z$.  Writing $n=p^am$, $p\nmid m$, and splitting
$e=p^ie'$,
$$\beta_\psi(n)=\sum_{e'\mid m}\nu(e')\Bigl[\alpha(p^a m')+(1-\psi(p))\sum_{j=0}^{a-1}\alpha(p^jm')\Bigr]
 =\sum_{e'\mid m}\nu(e')\,\frac{S_a(m')}{p^a m'},\qquad m'=m/e',$$
and $v_p(S_a)\ge3a$ gives $v_p(\beta_\psi(n))\ge2a$. $\square$

**Verification** (`37_thm.log`): $v_p(S_a(m))\ge3a$ for all $a\ge1$ and all $m$ with
$p^{2a}m^2\le3000$, for every $2\le p\le53$ **except** the single prime dividing the
denominator of $c$ ($p=2$ for $f_{4a}$, $p=3$ for $f_{4b}$).  Minimal excess $0$, so the
theorem is sharp.  The quotients $64S_a(1)/p^{3a}$ come out as integers as predicted
(e.g. $p=7,a=1$: $230840825570752=g_{147}(1)$).

**Scope.**  The proof uses only: (i) $\dim S_4(\mathrm{SL}_2(\mathbf Z))=0$, so a weakly
holomorphic plus form is determined by its principal part; (ii) the integral basis
$\{g_M=q^{-M}+O(q)\}$ of $S^{!,+}_{5/2}(\Gamma_0(4))$; (iii) $p^2>m_0$ and $c\in\mathbf Z_p$.
So it covers **every** entry of PZ's Table 1 (all of them are $\Psi(\kappa f_{m_0})$ for an
explicit rational $\kappa$), and in particular explains the companion strand's numerical
finding that $n^2\mid\beta_\psi(n)$ holds across the whole table.  It does **not** cover
Cooper's rows, because at level $4N$ hypotheses (i) and (ii) both fail.

---

## 3.  Where Cooper's rows differ from Paşol–Zudilin's  **[verified]**

Re-run of $(S)$ for row $s_7$ (`lib.gp`, $m\le400$, all $3\le p\le59$):

| $\psi$ | result |
|---|---|
| $\mathbf 1$ | **PASS at every $p$** |
| $\chi_{-3}$ | PASS exactly at $\chi_{-3}(p)=+1$ ($p=7,13,19,31,37,43$), FAILS at every other $p$ |

So $\psi_{s_7}=\mathbf 1$ unambiguously (the two cannot both hold: subtracting gives
$2c'(m)\equiv0\bmod p^2$ at inert $p$).  Combining with §2.5:

> **Structural dichotomy.**  A weight-4 magnetic form whose weight-$5/2$ input has
> principal part concentrated in a *single* index $-m_0$ has $\lambda_p=\chi_{-m_0}(p)\,p$,
> hence satisfies Cooper's congruence with $\psi=\chi_{-m_0}$ and *not* with $\psi=\mathbf1$,
> and has $n^3\mid a(n^2)$ only at the split primes.  Cooper's rows $s_7,s_{10}$ have
> $\psi=\mathbf1$ and $n^3\mid a(n^2)$ at *every* prime.  Therefore **their weight-$5/2$
> input cannot be a scalar Kohnen-plus form whose principal part at $\infty$ is a multiple
> of $q^{-3}$ alone** — at least not with the level-one uniqueness mechanism intact.

The natural resolution (and it is the one the trace-formula strand independently arrived
at): the correct input is the **twisted** (genus-character) lift.  In the Gross–Kohnen–Zagier /
Bringmann–Ono twisted-trace picture with fundamental discriminant $D_0=-3$, the genus
character $\chi_{D_0}$ enters the Hecke action and multiplies the naive eigenvalue by a
further $\bigl(\tfrac{D_0}{p}\bigr)$, giving
$$\lambda_p=\Bigl(\tfrac{D_0}{p}\Bigr)^{2}p=p\qquad\text{for every }p\nmid D_0 ,$$
which is exactly the $\lambda\equiv p\ (p^3)$ required by the Proposition of §2.5 — and
therefore exactly Cooper's $\psi=\mathbf 1$ with $n^3\mid a(n^2)$ at *every* prime.
**[conjectural, but forced by §2.5 + the verified $\psi_{s_7}=\mathbf 1$]**

This is, I think, the single most useful output of this run: it says the Shimura–Borcherds
route to Conjecture 4.1 must go through the **twisted** lift, not the plain one, and it says
*why* — the genus character is what squares the Kronecker symbol and kills the split/inert
dichotomy that ruins the level-one forms at half the primes.

---

## 4.  Feasibility at level 28 (Task 4)

### 4.1  What PARI can do  **[verified]**

PARI/GP 2.15.4 handles half-integral weight at level 28 **without difficulty** (`36_mf28.log`):

| call | result | time |
|---|---|---|
| `mfinit([28,5/2])` | $\dim M_{5/2}(\Gamma_0(28))=8$, `mfparams = [28, 5/2, 1, 4, t-1]` | instant |
| `mfinit([28,5/2],1)` | $\dim S_{5/2}(\Gamma_0(28))=4$ | instant |
| `mfkohnenbasis(mfinit([28,5/2],1))` | **$1$-dimensional**, $=-B_1+3B_4$ | instant |
| `mfkohnenbasis(mfinit([4,5/2],1))` | empty (correct: $S^+_{5/2}(\Gamma_0(4))\cong S_4(\mathrm{SL}_2(\mathbf Z))=0$) | instant |
| `mfinit([28,29/2])` | $\dim M_{29/2}(\Gamma_0(28))=56$ | seconds |
| `mfinit([28,53/2])` | $\dim M_{53/2}(\Gamma_0(28))=104$ | seconds |
| `mfinit([28,5/2],3)` (Eisenstein) | **error**: "half-integral weight Eisenstein space is not yet implemented" | — |

So: **PARI can do it.**  The only missing piece in the package is the half-integral
Eisenstein space (irrelevant here, since we work inside $M$ and $S$), and there is no
built-in weakly holomorphic constructor — but `{g/\Delta(4\tau)^r\}` supplies one.

The $1$-dimensional Kohnen plus subspace of $S_{5/2}(\Gamma_0(28))$ is the Shimura
preimage of the weight-$4$ newform on $\Gamma_0(7)$; its expansion begins
$-q+3q^4-7q^8+5q^9+5q^{16}-14q^{21}-11q^{25}+7q^{28}+28q^{29}+\cdots$,
supported exactly on $n\equiv0,1\ (4)$ **[verified]**.

### 4.2  The weakly holomorphic search  (`38_wh.gp`, `38_wh_out.txt`)

Set $W_r(N)=\{g/\Delta(4\tau)^r : g\in M_{5/2+12r}(\Gamma_0(4N))\}$ — weight $5/2$, poles only
at the cusps, pole order $\le 4r$ at $\infty$.  Inside it define the linear subspace
$$V_r(N)=\Bigl\{f\in W_r(N):\ a(n)=0\ \text{for } n\equiv2,3\ (4);\ a(n)=0\ \text{for } n<0,\ n\ne-3;\ a(0)=0\Bigr\},$$
i.e. Kohnen-plus, cuspidal, with the only pole at $\infty$ of order $3$.  No normalisation is
imposed (at level $4$ the answer is $f_{4a}$ with $a(-3)=\tfrac1{64}$, so fixing $a(-3)=1$
would be wrong by a factor $64$).  The test is then a pure linear-algebra question:

> is there $f\in V_r(28)$ with $a(m^2)=m\,\beta_{s_7}(m)$ for $m=1,\dots,20$?

with a level-$4$ **control** which must return $f_{4a}$ (up to scale).
Results are in `38_wh_out.txt`; see §4.3.

### 4.3  Outcome

**Control (level 4).**  $\dim V_1(1)=1$, spanned by $g_3=q^{-3}+O(q)$, and the unique
solution of the $20$ equations is the multiple $\tfrac1{64}g_3$ — i.e. **exactly $f_{4a}$**:
$$a(-3)=\tfrac1{64},\ a(0)=0,\ a(1)=1,\ a(4)=-506,\ a(8)=\tfrac{131565}{64},\ a(9)=180249,\ a(12)=-66516,\ \dots$$
So the pipeline is validated end to end.

**Level 28.**

| $r$ | $\dim M_{5/2+12r}(\Gamma_0(28))$ | $\dim V_r(7)$ | constraints | verdict |
|---|---|---|---|---|
| $1$ | $56$ | $17$ (all with pole order exactly $3$) | $m\le20$ | **no solution** |
| $2$ | $104$ | $31$ | $m\le20$ | a solution exists — but $31$ parameters against $20$ equations, so the test is **underdetermined and proves nothing** |
| $2$ | $104$ | $31$ | $m\le40$ | DECISIVE_PLACEHOLDER |

The $r=1$ answer is the honest one at that pole bound; the $r=2$/$m\le20$ "match" must not
be read as a positive identification.  The decisive run ($40$ equations, $31$ unknowns) is
`38_wh40.gp`/`38_wh40_out.txt`.

**Caveat on any negative answer.**  $W_r$ bounds the pole order at *every* cusp of
$\Gamma_0(28)$, so a failure at $r=1,2$ does not by itself exclude a preimage with larger
poles at the other five cusps.  What makes the negative reading credible is §2.7 + §3:
at level $4N$ the scalar single-index model would force $\lambda_p=\chi_{-3}(p)p$, i.e.
$\psi=\chi_{-3}$, which is **refuted** by the verified $\psi_{s_7}=\mathbf 1$.

---

## 5.  Traces of singular moduli: which known theorem is closest

Under the trace identification found by the companion strand
($\beta_{s_7}(m)=i\sqrt3\sum_Q\chi_{-3}(Q)\,\omega_Q^{-1}\hat f(\alpha_Q)$, a **twisted** CM
trace with genus character $\chi_{-3}$), the master conjecture $n^2\mid\beta(n)$ becomes

$$m^2 \ \Big|\ \operatorname{Tr}^{\mathrm{tw}}_{-3m^2}\bigl(R_{-2}f\bigr)\qquad\text{for all }m,$$

a divisibility of twisted traces of (a Maass-raised, weight $-2$) singular moduli by
$m^2$ — three powers of $m$ once one writes $a(m^2)=m\beta(m)$.

The closest published results are **Ahlgren–Ono**'s theorems on the divisibility of traces
of singular moduli by primes (Compositio Math. **141** (2005), and the companion papers of
Boylan, Jenkins, Osburn, Bringmann–Ono, Guerzhoy).  Their mechanism is precisely the one
isolated in §2: the generating function of the traces is a weakly holomorphic form of
half-integral weight, one applies $T_{p^2}$, and the principal part forces a congruence.
The differences, in the direction that matters here, are:

| Ahlgren–Ono type theorem | what Cooper needs |
|---|---|
| weight $3/2$ input ($k=1$), so $p^{k-1}=1$ and Lemma 2 gives nothing; the divisibility comes from the $p$-adic geometry of the supersingular locus / Edixhoven | weight $5/2$ ($k=2$), where Lemma 2 already gives one power of $p$ **for free** |
| a **congruence mod $p$** (or mod a fixed power) for a single $p$ | an **exact eigenvalue** $\lambda_p=p$ for every $p$, giving $p^{3j}$ |
| untwisted traces, so the Kronecker symbol $\bigl(\tfrac{-d}{p}\bigr)$ appears and the results are split/inert-conditional | twisted traces, where the genus character squares that symbol away — this is what makes $\psi=\mathbf 1$ possible |

So: the weight-$5/2$ input accounts for the *first* power of $m$ (PZ's theorem), the exact
eigenvalue $\lambda_p=p$ accounts for the remaining two, and the genus-character twist is
what makes $\lambda_p=p$ (rather than $\chi_{D_0}(p)p$) available at every prime.  I know of
no published theorem that gives the exact-eigenvalue statement; that is the gap.

---

## 6.  Verdict table

| # | statement | status |
|---|---|---|
| 1 | $\beta(m)=(A/n\star\mu)(m)=a(m^2)/m$ for the $k$ even, $D=1$ lift | **[proved]** (§1) |
| 2 | $\beta_3(m)=(A/n^2\star\mu\chi_{-3})(m)=a(3m^2)/m^2$ for $k=3$, $D=-3$ | **[proved]** (§1), **[verified]** $m\le300$ |
| 3 | $g_0,g_1,g_2$ reproduce PZ §4; $f_{4a}=\tfrac1{64}q^{-3}+q-506q^4+\cdots$, $f_{4b}=-\tfrac1{108}q^{-4}+q+1222q^4+\cdots$; $64f_{4a},108f_{4b}$ integral | **[verified]**, $n\le400$ |
| 4 | dictionary $a(m^2)=m\beta(m)$ for $F_{4a},F_{4b}$ | **[verified]**, $m\le20$, exact |
| 5 | lift formula reproduces $\Delta/E_4^2$, $E_4\Delta/E_6^2$ | **[verified]**, $n\le50$ |
| 6 | the plain-text `pz.txt` mangles Lemma 1's leading coefficients ($q^{-3}$ should be $\tfrac1{64}q^{-3}$) | **[established]** |
| 7 | $m\mid a(m^2)$ is sharp for $220/299$ (resp. $212/299$) of $m\le300$; $\gcd_m\beta(m)=1$ | **[verified]** |
| 8 | $v_p(\beta(p^e))=2e$ if $\chi_{D_0}(p)=+1$ and $=0$ if $\chi_{D_0}(p)=-1$ (good $p$); so $v_p(a(p^{2e}))=3e$ resp. $e$ | **[verified]** $n\le1250$ |
| 9 | (T1) $f\mid T_{p^2}=\chi_{D_0}(p)\,p\,f+p^{3}c\,g_{m_0p^2}$, exactly | **[proved]** (§2.2, level 1), **[verified]** $p\le17$ |
| 10 | (T2) $g_{m_0p^{2r}}\mid T_{p^2}=g_{m_0p^{2r-2}}+p^{3}g_{m_0p^{2r+2}}$ | **[proved]** (§2.2) |
| 11 | (C1) $A(pn)\equiv\lambda A(n)\ (p^{2k-1})$ | **[proved]**, **[verified]** $n\le1250$, good $p\le19$ |
| 12 | (C2) $A(p^2n)\equiv\lambda A(pn)\ (p^{2(2k-1)})$ — the extra cancellation | **[proved]**, **[verified]** |
| 13 | (C3) $A(p^rn)\equiv\lambda^{r}A(n)\ (p^{r+2})$ — Cooper `eq:magnetic` (i) | **[proved]**, **[verified]** $r\le5$ |
| 14 | Cooper's $(S)$ **holds** for $F_{4a}$ ($\psi=\chi_{-3}$, all $p\ne2$), $F_{4b}$ ($\psi=\chi_{-4}$, all $p\ne3$), $F_6$ ($c''$, mod $p^3$, all $p\ge5$) | **[verified]** $p\le59$, $m\le800$ |
| 15 | the ledger's "negative control: $(S)$ fails for PZ's level-one forms" | **[refuted]** — it fails only for the wrong $\psi$ |
| 16 | $v_p(a(p^{2j}m^2))\ge 3j+v_p(a(m^2))$ $\iff$ $\lambda\equiv p\ (p^3)$; and then $a(p^{2j}m^2)=p^{3j}h_j(m^2)$ | **[proved]** (§2.5) |
| 17 | consequently the level-one forms have $n^3\mid a(n^2)$ **exactly at the split primes** | **[proved]** + **[verified]** $n\le1250$ |
| 18 | Cooper's rows have $\psi=\mathbf 1$, hence need $\lambda_p=p$ at every $p$; a single-index principal part cannot deliver this | **[proved]** modulo the level-$4N$ uniqueness caveat of §2.6 |
| 19 | the twisted (genus-character) lift restores $\lambda_p=p$ at every $p$ | **[conjectural]**, but forced by 16–18 |
| 20 | PARI can build $M_{5/2}(\Gamma_0(28))$, $S_{5/2}$, `mfkohnenbasis`, and $M_{5/2+12r}(\Gamma_0(28))$ for $r=1,2$ | **[verified]** — no obstruction |
| 21 | existence in $V_r(28)$ of a form with $a(m^2)=m\beta_{s_7}(m)$, $r=1,2$ | see §4.3 |

---

## 7.  Files

| file | what it does |
|---|---|
| `30_calib.gp/.log` | Task 1: $A(n)$, $\beta$, $a(m^2)$ for $F_{4a},F_{4b},F_6$, $m\le300$; exact power of $m$ dividing $a(|D|m^2)$; $v_p(\beta(p^e))$ tables.  Data in `30_A_*.txt`, `30_beta_*.txt` |
| `31_congr.gp/.log` | Cooper's $(S)$ for the three level-one forms with $\psi=\chi_{D_0}$ and with $\psi=\mathbf1$; $\mathrm{rad}(n)^2\mid\beta_\psi$.  Data in `31_betapsi_*.txt` |
| `32_nsq.gp/.log` | the strong divisibility $v_p(\beta_\psi(n))\ge w\,v_p(n)$, per prime |
| `33_halfint.gp/.log`, `33_coeffs.gp` | Task 2: independent construction of $g_0,g_1,g_2,f_{4a},f_{4b}$ at level 4; dictionary check.  Data in `33_f4a.txt`, `33_f4b.txt` |
| `34_hecke.gp/.log` | Task 3: numerical proof of the exact eigen-identity (T1), with a wrong-$\lambda$ control |
| `35_tower.gp/.log`, `39_c3.gp/.log` | consequences (C1),(C2),(C3) of (T1)+(T2), and the split/inert dichotomy |
| `36_mf28.gp/.log` | Task 4: what PARI's `mf` package can do at level 28 |
| `37_thm.gp/.log` | verification of the §2.7 theorem: $v_p(S_a(m))\ge3a$ and $64S_a(1)/p^{3a}\in\mathbf Z$ |
| `38_wh.gp` → `38_wh_out.txt`, `38_wh40.gp` → `38_wh40_out.txt` | Task 4: the weakly holomorphic search in $V_r(28)$, with the level-4 control |

Nothing outside this directory was touched, and nothing was committed.
