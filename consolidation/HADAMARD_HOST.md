# The Hadamard host of (Zudilin, Nesterenko), and the arithmetic holonomy bound on it

*Fable (Opus 5), 2026-08-22.  Scripts: `lattice/hadamard_host/`
(`00_csv2gp.py`, `01_rows.gp`, `02_bridge.gp`, `lib_fit.gp`, `03_rowrecur.gp`,
`04_hadamard.gp`, `05_fitprod.gp`, `06_operator.gp`, `07_denoms.gp`,
`08_odeorder.gp`, `08b_odeorder.gp`, `09_indep.gp`, `10_folds.gp`,
`11_inventory.gp`, `12_adelic.py`, `13_baselines.gp`, `14_pairs.gp`,
`15_conformal.py`).  This is the computation set in `CDT_UNPACKED.md` §6.
Inputs and normalisations from `ADELIC_HOLONOMY.md` §2 (+ §7), `CDT_FINDER.md` §§1–3,
`CDT_NONCONGRUENCE.md` §§1–3, `CATALAN_AUDIT.md`, `CATALAN_EXPLICIT.md`,
`POSITIVITY_PROGRAM.md` §4.3.  Tags: **[proved] [verified] [measured] [open]**.
**No irrationality claim is made anywhere.***

---

## 0. Verdict

| claim | verdict |
|---|---|
| The Hadamard host of the two rows is $\mathbf P^1\setminus\{0,\ s_1t_1,\ s_1t_2,\ s_2t_1,\ s_2t_2,\ \infty\}$ with the four products of the rows' singular points, and its module has recurrence order exactly $4$ | **[verified]** §2 |
| $W:=A_Z\odot B_N-A_N\odot B_Z$ is unconditionally regular at the innermost point $s_1t_1$, and is the **unique** such $\mathbf Q$-rational element of the module (unless $1,G,G^2$ are $\mathbf Q$-dependent) | **[proved]**+**[verified]** §4 |
| $W$'s $2$-adic slope is exactly $-2$: $v_2(w_n)=-2n+O(1)$, i.e. the $24n$ cross-divisibility minus the $26n$ geometric normalisation | **[verified]** $n\le400$, §3 |
| $W$'s LCM denominator type is $[1,\dots,6n]^2$ — **the Hadamard product does not double $k$** (the naive bound $k=4$ is wrong; $k=4$ occurs only for $B_Z\odot B_N$ and for the doubly-conditional function) | **[verified]** $n\le400$, §3 |
| The **pure polylogarithm module is unavailable**: $\mathrm{Li}_j\odot\mathrm{Li}_k=\mathrm{Li}_{j+k}$ at the product point, and the four product points form a **single Galois orbit** (the leading polynomial is irreducible over $\mathbf Q$), so every $\mathbf Q$-rational pure or algebraic function branches at *all four*, including the one $W$ removes | **[proved]** §5 |
| Hence at most **one** function of the inventory (the constant) is denominator-free, and $\tau^\flat\in[10.67,12]$ for every inventory on this host | **[proved]** §6 |
| **The entry test fails, and by a wide margin.** Best configuration $(K)\ \{1,W,\theta W\}$: $\log|\varphi'(0)|=+0.4891$, $\tau^\flat=10.6667$, $\gamma_2=-0.9242$, **entry $=-11.10$**, margin $-33.79$, deficit $-16.90$ per function | **[verified]** §6 |
| Even at the hard Landau/Kodaira ceiling $\log(16\,d_{\min})=+2.340$ the entry deficit is $-9.25$: **no contour design can rescue it** | **[verified]** §6 |
| The dominant term is $\tau^\flat$, and it is inherited from the **rows**, not created by the Hadamard product: the $\mathrm{lcm}$ rate $b=6$ comes from $D_{6n}^2$ in both aligned rows. The un-sectioned Zudilin row alone scores $-4.20$; the $3$-section $-15.59$; Nesterenko $-16.89$; the Hadamard host $-16.90$ | **[verified]** §7 |
| The `POSITIVITY_PROGRAM.md` adjacent pair $(j_0,j_0+1)$ is **markedly better** — its Hadamard minor has type $[1,\dots,6n]^{\mathbf 1}$ ($\sigma_m=6$, half) and $2$-adic slope $0$ — but still fails entry by $\ge2.6$ (ceiling) / $\ge4.0$ (Koebe) | **[measured]** $n\le24$, §8 |

One sentence: *the two-host object $W$ is exactly the doubly-small function `CDT_UNPACKED.md` §5 predicted — unconditional, overconvergent by $14.44$ nats, $2$-adic slope $24$ before normalisation — and it is still $16.9$ per function short of a CDT contradiction, because the hypergeometric rows carry $\mathrm{lcm}(1,\dots,6n)^2$ where the modular rows carry $\mathrm{lcm}(1,\dots,n)^2$, and because the Hadamard host's Galois-irreducible singular set destroys the pure module that CDT's architecture depends on.*

---

## 1. Set-up: the two rows, exactly

Normalisation (both rows scaled so that $b_n/a_n\to G$):

$$a^Z_n:=Q_{3n},\quad b^Z_n:=P_{3n}\qquad\text{(Zudilin, at the aligned index }m=3n),$$
$$a^N_n:=4B_n,\quad b^N_n:=C_n\qquad\text{(Nesterenko }(4,7)).$$

$Q_m,P_m$ from Zudilin's second-order recurrence; $B_n=\sum_{j=0}^{3n}A_2(n,j)$ from the
closed form (2.11); $C_n=U_n/(4^{7n}D_{6n}^2)$ from the published exact integer row.
**[verified]** (`01_rows.gp`) All of
$$X_n=2^{e(3n)}D_{6n}^2Q_{3n},\quad Y_n=2^{e(3n)}D_{6n}^2P_{3n},\quad V_n=4^{7n+1}D_{6n}^2B_n$$
agree **digit-for-digit** with `catalan-2-row-denominators/{zudilin,nesterenko}_rows.csv`
for all $n\le98$ ($e(m)=\min(6m,4m+3+\lfloor\log_2(2m-1)\rfloor)$, $D_N=\mathrm{lcm}(1..N)$).
Rows built to $n=400$ in exact rationals; $C_n$ extended from $98$ to $400$ by the
row's own recurrence (§2).

**Measured rates** (`01_rows.gp`, `13_baselines.gp`), $n=400$ against the predicted limits:

| quantity | $n=400$ | limit |
|---|---|---|
| $\tfrac1n\log|a^Z_n|$ | $7.19605$ | $15\log\varphi=7.218177$ |
| $\tfrac1n\log|a^Z_nG-b^Z_n|$ | $-7.23872$ | $-15\log\varphi$ |
| $\tfrac1n\log|a^N_n|$ | $7.63349$ | $2\log T=7.650713$, $T=\tfrac{3303+437\sqrt{57}}{144}$ |
| $\tfrac1n\log|a^N_nG-b^N_n|$ | $-7.32625$ | $2\log T^-=-7.310915$ |

---

## 2. The $24n$ bridge, reproduced; and the host

### 2.1 The cross-divisibility **[verified]**

`02_bridge.gp`, exact integers, $h_n:=X_nU_n-V_nY_n$:

$$\min_{1\le n\le98}\bigl(v_2(h_n)-24n\bigr)=+14,\qquad
\text{least squares on }n\in[50,98]:\ v_2(h_n)=24.131\,n+35.82 .$$

$v_2(h_n)-24n$ at $n=20,40,60,80,98$ is $33,38,45,43,49$.  This reproduces
`CATALAN_AUDIT.md` §4(d) ($24.06n$) from independently rebuilt rows.
**It remains an exact-per-$n$ observation with no proved lower bound.** **[open]**

The identity tying it to the function-theoretic object is exact:
$$h_n=2^{\,e(3n)+14n}\,D_{6n}^4\;w_n,\qquad w_n:=a^Z_nb^N_n-a^N_nb^Z_n=Q_{3n}C_n-4B_nP_{3n},$$
verified for every $n\le98$.  Since $e(3n)+14n=26n+O(\log n)$,
$$v_2(w_n)=v_2(h_n)-26n+O(\log n)=\mathbf{-2n}+O(\log n).$$
Measured directly: $v_2(w_n)=-199,-399,-598,-799$ at $n=100,200,300,400$, i.e.
$v_2(w_n)=-(2n-1)$ on the nose. **[verified]**

### 2.2 The two rows' own operators **[verified]** (`03_rowrecur.gp`)

Minimal recurrences over $\mathbf Q$ (fitted mod $p$, reconstructed by CRT, verified
exactly on all indices, holding out the last $60$–$100$):

| row | $(r,d)$ | leading polynomial $P_{\rm gf}(x)$ | singular points |
|---|---|---|---|
| $(a^Z,b^Z)$ joint | $(2,22)$ | $x^2+1364x-1$ | $s_1=\varphi^{-15}=7.3313744\cdot10^{-4}$, $s_2=-\varphi^{15}=-1364.0007331$ |
| $a^N$ | $(2,18)$ | $2^{16}x^2-98077689\,x+6^6$ | $t_1=T^{-2}=4.7570467\cdot10^{-4}$, $t_2=(T^-)^{-2}=1496.5462925$ |

$b^N=C_n$ satisfies the **same** recurrence as $a^N$ (failures $=0$ on $n\le96$), which
is what licenses the extension of $C_n$ to $n=400$.

### 2.3 The Hadamard host **[verified]** (`05_fitprod.gp`, `06_operator.gp`)

All four products $A_Z\odot A_N$, $A_Z\odot B_N$, $A_N\odot B_Z$, $B_Z\odot B_N$ and
$W$ satisfy one recurrence, and the **minimal joint order is exactly $4$**:

* no order $\le3$ recurrence of degree $\le300$ exists (mod-$p$ nullity $0$);
* at $(r,d)=(4,96)$ the kernel is $1$-dimensional; reconstructed exactly over $\mathbf Q$
  ($16$ primes), verified with $0$ failures on $n=1,\dots,396$ for all four products *and*
  for $W$;
* $W$ alone has no operator of order $\le3$ and degree $\le94$.

$$P_{\rm gf}(x)=2^{32}x^4+8767272897478656\,x^3-3930485805307953\,x^2-6241544865490176\,x+6^{12},$$
**irreducible over $\mathbf Q$**, with roots

| root | value | $=$ |
|---|---|---|
| $s_1t_1$ | $3.4875690\cdot10^{-7}$ | $\varphi^{-15}T^{-2}=e^{-14.86889}$ |
| $s_1t_2$ | $+1.0971741$ | $\varphi^{-15}(T^-)^{-2}=e^{+0.09274}$ |
| $s_2t_1$ | $-0.64886152$ | $-\varphi^{15}T^{-2}=-e^{-0.43254}$ |
| $s_2t_2$ | $-2041290.2401$ | $-\varphi^{15}(T^-)^{-2}=-e^{+14.52909}$ |

**Hadamard's theorem is confirmed exactly**: the singular set is precisely the set of
products, to all printed digits. **[verified]**

*Order–degree curve* (`08_odeorder.gp`, `08b`): $(4,96),(5,63),(6,52),(7,46),(8,43),
(9,41),(10,39),(11,38),(12,37),(13,36),(15,35),(17,34),(24,33),(30,32),(36,32)$; the
data ($4\times392$ rows) cannot reach $r\ge44$.  So the **differential** order of the
module is $\le32$; a recurrence of order $4$ is the minimal-order representative.
The exact number is **[open]** and does not affect anything below (§6).

---

## 3. Denominator types and slopes, measured exactly

`07_denoms.gp`, `11_inventory.gp`, exact rational coefficients, $n\le400$.
"$k$" is the **sharp** exponent with $D_{6n}^k\,c_n\in\mathbf Z[\tfrac12]$;
$b$ is read off from the largest prime dividing the odd denominator
($P_{\max}/n=5.9975$ at $n=400$ throughout, so $b=6$).

| function | $\tfrac1n\log|c_n|$ ($n{=}400$) | limit | $\varsigma_2=v_2(c_n)/n$ | $k$ (rate $b=6$) |
|---|---|---|---|---|
| $A_Z\odot A_N$ | $14.82954$ | $14.86889$ | $-26$ | $\mathbf 0$ (denominator-free!) |
| $A_Z\odot B_N$ | $14.82932$ | $14.86889$ | $-26$ | $2$ |
| $A_N\odot B_Z$ | $14.82932$ | $14.86889$ | $-26$ | $2$ |
| $B_Z\odot B_N$ | $14.82911$ | $14.86889$ | $-26$ | $\mathbf 4$ |
| $\mathbf W$ | $\mathbf{0.39477}$ | $\mathbf{0.43254}$ | $\mathbf{-2}$ | $\mathbf 2$ |
| conditional $b(A_Z\odot B_N)-a(A_Z\odot A_N)$ | — | $-0.09274$ | $-26$ | $2$ |
| doubly conditional $(bB_Z{-}aA_Z)\odot(bB_N{-}aA_N)$ | — | $-14.52909$ | $-26$ | $4$ |

Source-row denominators, for the accounting: $\mathrm{den}(Q_{3n})$ and
$\mathrm{den}(B_n)$ are **pure powers of $2$** ($\varsigma_2=-12$ and $-14$);
$\mathrm{den}(P_{3n})=2^{12n}[1..6n]^2$, $\mathrm{den}(C_n)=2^{14n}[1..6n]^2$.
Hence $A_Z\odot A_N$ is denominator-free away from $2$, and $k$ doubles to $4$ only
in $B_Z\odot B_N$.  **The caution of `CDT_UNPACKED.md` §5(i) — "Hadamard products
square the denominator exponents" — is false for $W$ on this host**: $W$ has $k=2$,
the same as each row.

**$v_p$ profile at $n=400$** (exponent of $p$ in the denominator, against
$\lfloor\log_p6n\rfloor$, `07_denoms.gp` §4): e.g. for $W$,
$3\!:\!11/7$, $5\!:\!8/4$, $7\!:\!3/3$, $11\!:\!4/3$, $13\!:\!4/3$, $17\!:\!3/2$,
$19\!:\!4/2$, $23\!:\!2/2$, $47\!:\!3/2$ — everywhere $\le2\lfloor\log_p6n\rfloor$
and attaining $2$ at the top layer, confirming type $[1..6n]^2$ and no more.
For $B_Z\odot B_N$ the same primes read $26/7,16/4,9/3,10/3,11/3,7/2,8/2,\dots$,
confirming $k=4$.

**The conditional function's slope is that of $A_Z\odot A_N$**, as the task predicted:
both are $-25.9575$ at $n=400$ ($\to-26$). $2$-adically $b\,b^N_n-a\,a^N_n\approx
(bG_2-a)a^N_n$ with $G_2\ne a/b$, so nothing cancels — the same mechanism as
`ADELIC_HOLONOMY.md` §3.

---

## 4. Fold-regularity on this host **[proved]**+**[verified]**

By Hadamard's theorem, $\mathrm{sing}(f\odot g)\subseteq\{\text{products}\}$, so a
factor that is regular at $t_1$ deletes the whole column $\{s_1t_1,s_2t_1\}$.
Measured coefficient growth (`10_folds.gp`, $n\le400$, $G$ exact) versus
$-\log|{\rm nearest\ singularity}|$:

| function | $\tfrac1n\log|c_n|$ at $n=400$ | limit | regular at | singular at |
|---|---|---|---|---|
| $A_Z\odot A_N$ | $14.82954$ | $14.86889$ | — | all four |
| $\mathbf W$ (**unconditional**) | $0.39477$ | $0.43254$ | $s_1t_1$ | $s_1t_2,\,s_2t_1,\,s_2t_2$ |
| $A_Z\odot(B_N-GA_N)$ (**conditional**) | $-0.13020$ | $-0.09274$ | $s_1t_1,\,s_2t_1$ | $s_1t_2,\,s_2t_2$ |
| $A_N\odot(B_Z-GA_Z)$ (**conditional**) | $0.39477$ | $0.43254$ | $s_1t_1,\,s_1t_2$ | $s_2t_1,\,s_2t_2$ |
| $(B_Z-GA_Z)\odot(B_N-GA_N)$ (**doubly**) | $-14.56497$ | $-14.52909$ | $s_1t_1,s_1t_2,s_2t_1$ | $s_2t_2$ |

So on this host:
* **the hypothesis $G\in\mathbf Q$ removes a whole *column* of the product grid**
  (one factor becomes fold-regular), and using it on both factors removes three of the
  four points;
* **$W$ removes the single innermost point $s_1t_1$ unconditionally.**  Its
  overconvergence over $A_Z\odot A_N$ is $14.86889-0.43254=\mathbf{14.436}$ nats
  — exactly the prediction of `CDT_UNPACKED.md` §5.
* nothing removes $s_2t_2$.

> **Proposition (uniqueness of $W$).** Write $a^Z_n\sim\alpha\lambda_Z^n$,
> $b^Z_n\sim G\alpha\lambda_Z^n$, $a^N_n\sim\beta\lambda_N^n$,
> $b^N_n\sim G\beta\lambda_N^n$ ($\lambda_Z=\varphi^{15}$, $\lambda_N=T^2$).  For
> $c\in\mathbf Q^4$ the combination
> $c_1A_Z\odot A_N+c_2A_Z\odot B_N+c_3A_N\odot B_Z+c_4B_Z\odot B_N$ has coefficient
> growth $\lambda_Z\lambda_N$ unless $c_1+G(c_2+c_3)+G^2c_4=0$.  If $1,G,G^2$ are
> $\mathbf Q$-linearly independent this forces $c_1=c_4=0$, $c_3=-c_2$, i.e.
> **the combination is a multiple of $W$**. **[proved]**
> (If they are dependent, $G$ is algebraic of degree $\le2$ — a far stronger statement
> than anything at issue.)

So the unconditional "doubly-small" module of `CATALAN_THREE_PERIOD.md` type on this
host is exactly $\mathbf Q\cdot W$ plus its $\theta$-orbit — one generator, as on the
modular hosts.

**$\mathbf Q(x)$-independence [verified]** (`09_indep.gp`, rank of
$\{x^jf_i\}$ mod $2^{61}-1$, series to $x^{400}$): full rank at every degree $\le5$ for
$\{1,W\}$, $\{1,W,\theta W\}$, $\{1,W,\theta W,\theta^2W\}$,
$\{1,W,\theta W,\theta^2W,\theta^3W\}$, $\{1,W,\theta W,\mathrm{COND},\theta\,\mathrm{COND}\}$
and $\{1,A_Z\odot A_N,W,\theta W\}$.  (A check of, not a substitute for, the hypothesis
of CDT's Lemma 12.1.1. **[open]**)

---

## 5. Why there is no pure module here **[proved]**

On a one-row host the pure inventory is $\mathrm{Li}_j(\lambda_2x)$ — denominator-free,
which is what makes CDT's $u_j\ge1$ (and $u_2=3$) possible.  On a Hadamard host:

$$\mathrm{Li}_j(x/s)\odot\mathrm{Li}_k(x/t)=\sum_{n\ge1}\frac{(st)^{-n}}{n^{j+k}}x^n=\mathrm{Li}_{j+k}\!\left(\frac x{st}\right),$$

so the pure Hadamard orbit is just the polylogarithm orbit **at the product points**.
But $P_{\rm gf}$ is **irreducible over $\mathbf Q$** (§2.3): the four product points form a
single orbit of $\mathrm{Gal}(\mathbf Q(\sqrt5,\sqrt{57})/\mathbf Q)\cong(\mathbf Z/2)^2$.
Any $\mathbf Q$-rational combination $\sum c_{ij}\mathrm{Li}_w(x/(s_it_j))$ therefore has
Galois-conjugate coefficients, hence branches at **every** $s_it_j$ — in particular at
$s_1t_1$, the one point $W$ was built to remove.  The same argument applies to algebraic
functions: a $\mathbf Q$-rational algebraic function ramified at some $s_it_j$ is ramified
at all four, and one ramified only at $\{0,\infty\}$ is $x^{1/d}$, not a power series.

**Consequence.** Adding any pure function collapses the domain from
$\rho=1.63$ to $\rho=1.4\cdot10^{-6}$ (§6, last row of the table) — a loss of $13.97$
nats, far more than the $\tau^\flat$ it buys.  So on this host **$u_j=1$: only the
constant is denominator-free.**  This is the structural difference from CDT's host,
and it is not repairable by working harder.

---

## 6. The adelic bound: entry and margin **[verified]**

`12_adelic.py` (reusing `lattice/adelic_holonomy/adelic_bound.py`, calibrated to CDT's
own $13.9938/+0.0053$) and `15_conformal.py`.

### 6.1 The domains **[proved]**

Every relevant domain is the complement of two collinear rays,
$\Omega(a,b)=\mathbf C\setminus\bigl((-\infty,-a]\cup[b,+\infty)\bigr)$.  With
$M(w)=-\frac{w+a}{w-b}$ (which sends the two rays plus $\infty$ to $(-\infty,0]$) and
the Koebe map of the slit plane at $a/b$,
$$\varphi=M^{-1}\!\circ\!\Bigl(z\mapsto \tfrac ab\bigl(\tfrac{1+z}{1-z}\bigr)^2\Bigr):\ \mathbf D\to\Omega(a,b),\qquad
\boxed{\ |\varphi'(0)|=\frac{4ab}{a+b}\ }$$
($=4a$ for $b=\infty$: Koebe).  $\varphi$ is univalent (a composition of univalent maps),
so $\mathrm{BC}(\varphi)=\log|\varphi'(0)|$ by Grunsky — the convention of
`CDT_NONCONGRUENCE.md` §1.4.  Verified numerically in `15_conformal.py` (derivative and
injectivity on a $4\times240$ grid).

| module | slits at | $|\varphi'(0)|$ | $\log|\varphi'(0)|$ | Landau ceiling $\log(16d_{\min})$ |
|---|---|---|---|---|
| $W$, $\theta W,\dots$ | $-0.648862$, $+1.097174$ | $1.6309268$ | $+0.48915$ | $+2.34005$ |
| conditional $A_Z\odot(bB_N{-}aA_N)$ | $-2041290.24$, $+1.097174$ | $4.3886941$ | $+1.47903$ | $+2.86533$ |
| doubly conditional | $-2041290.24$ only | $8165160.96$ | $+15.91539$ | $+17.30168$ |
| everything (incl. $A_Z\odot A_N$) | $-0.648862$, $+3.48757\cdot10^{-7}$ | $1.3950\cdot10^{-6}$ | $-13.48260$ | $-12.09630$ |

*Kodaira/three-point map.* Architecture (D) of `CDT_NONCONGRUENCE.md` §3.2 needs the
module to have **one** finite nonzero branch point.  On this host that holds only for the
doubly-conditional function; $W$ has three and the conditional function two, so the
$\lambda$-map is **not admissible** for them and the Landau bound $16d_{\min}$ (the
$3$-punctured-sphere uniformisation, which ignores the extra punctures) is quoted as a
hard **ceiling**, not an achievable value.

*Rearrangement numerator* (CDT 6.0.15), computed on $|z|=r\to1$ (`15_conformal.py`):
$+0.0240$ for $W$, $+0.6265$ for the conditional, $+12.826$ for the doubly conditional,
$-14.335$ for the full module — in every case **below** $\log|\varphi'(0)|$, so CDT's
refined numerator is a genuine improvement of the *bound*.  It does not enter the
**entry** condition, which is what fails.

### 6.2 The table

$\tau^\sharp=0$ everywhere (no Eichler integrals).  Denominator layers: $b_j=6$,
$u_j=1$; $k$ layers where $k$ is the measured exponent.  Slopes as measured in §3.
"entry (no 2-adic)" is **not a valid bound** — it is what one would get if the geometric
$2$-power denominators did not exist, quoted only to separate the terms.

| inventory (architecture (K), univalent) | $m$ | $\sigma_m$ | $\tau^\flat$ | $\gamma_2$ | **entry** | entry (no 2-adic) | **margin** |
|---|---|---|---|---|---|---|---|
| $\{1,W,\theta W\}$ | $3$ | $12$ | $10.6667$ | $-0.9242$ | $\mathbf{-11.1017}$ | $-10.1775$ | $\mathbf{-33.794}$ |
| $\{1,W,\theta W,\theta^2W\}$ | $4$ | $12$ | $11.2500$ | $-1.0397$ | $-11.8006$ | $-10.7609$ | $-47.691$ |
| $\{1,W,\theta W,\theta^2W,\theta^3W\}$ | $5$ | $12$ | $11.5200$ | $-1.1090$ | $-12.1399$ | $-11.0309$ | $-61.189$ |
| $\{1,\mathrm{COND},\theta\,\mathrm{COND}\}$ | $3$ | $12$ | $10.6667$ | $-12.0146$ | $-21.2022$ | $-9.1876$ | $-65.086$ |
| $\{1,\mathrm{DBL},\theta\,\mathrm{DBL}\}$ | $3$ | $24$ | $21.3333$ | $-12.0146$ | $-17.4325$ | $-5.4179$ | $-68.213$ |
| $\{1,W,\theta W,\mathrm{COND},\theta\,\mathrm{COND}\}$ | $5$ | $12$ | $11.5200$ | $-11.7558$ | $-22.7866$ | $-11.0309$ | $-114.422$ |
| $\{1,W,\theta W,\mathrm{DBL},\theta\,\mathrm{DBL}\}$ | $5$ | $24$ | $19.2000$ | $-11.7558$ | $-30.4666$ | $-18.7109$ | $-152.822$ |
| $\{1,A_Z\!\odot\!A_N,W,\theta W\}$ | $4$ | $12$ | $11.2500$ | $-9.3575$ | $-34.0901$ | $-24.7326$ | $-122.878$ |

Architecture (D), available only for the doubly-conditional module:
$\log|\varphi'(0)|=+16.5252$, $\tau=21.3333$, $\gamma_2=-12.0146$, entry $-16.8227$
(archimedean-only $-4.8082$), margin $-67.975$.

At the **hard Landau ceiling** (upper bound on any admissible $\varphi$):

| module | $\log(16d_{\min})$ | $\tau$ | entry (adelic) | entry (archimedean only) |
|---|---|---|---|---|
| $\{1,W,\theta W\}$ | $+2.3401$ | $10.6667$ | $\mathbf{-9.2508}$ | $-8.3266$ |
| $\{1,\mathrm{COND},\theta\,\mathrm{COND}\}$ | $+2.8653$ | $10.6667$ | $-19.8159$ | $-7.8013$ |
| $\{1,\mathrm{DBL},\theta\,\mathrm{DBL}\}$ | $+17.3017$ | $21.3333$ | $-16.0462$ | $-4.0317$ |

### 6.3 Which term dominates

For the best configuration,
$$\text{entry}=\underbrace{+0.4891}_{\log|\varphi'(0)|}-\underbrace{10.6667}_{\tau^\flat}
+\underbrace{(-0.9242)}_{\gamma_2}=-11.1017 .$$

* $\tau^\flat$ is **$96\%$ of the deficit**.  It is $\sigma_m-\frac1{m^2}\sum_ju_j^2b_j
  =12-\frac{12}{m^2}$, and since $u_j=1$ is forced (§5), $\tau^\flat\ge10.667$ for every
  inventory, with equality at $m=3$.
* The $2$-adic term is a **penalty** of $0.92$ (the geometric $4^n$ denominator of $W$),
  not a gain.  Its scale-covariant twin: in the integral coordinate $x\mapsto x/4$,
  $W(4x)$ has integer numerators over $[1..6n]^2$, $\rho$ drops to $0.40773$, entry
  $=-11.5638$ and the margin is **identical** ($-33.794$) — Lemma G of
  `CDT_NONCONGRUENCE.md` §2, verified here.
* The archimedean budget is at most $+2.340$ (Landau).  **So the entry test fails by
  $\ge9.25$ for every contour and every inventory on this host.**

What would be needed: $\tau^\flat<2.340$ with $\sigma_m=12$ requires
$u/m>\sqrt{(12-2.34)/12}=0.897$, i.e. **at least $90\%$ of the inventory must be
denominator-free** — precisely the pure module that §5 shows does not exist over
$\mathbf Q$ on this host.

---

## 7. The cost of the Hadamard construction, isolated **[verified]**

`13_baselines.gp` + `12_adelic.py` run the *same* machine on each source row alone,
with $m=3=\{1,H,\theta H\}$, $H$ the row's conditional form (regular at the fold, so the
domain is the single-slit plane at $|c|=1/|\lambda_2|$):

| row | type | $\varsigma_2$ | $|c|$ | $\tau$ | $\gamma_2$ | (K) entry | (D) entry | (D) margin | **(D) deficit** |
|---|---|---|---|---|---|---|---|---|---|
| Zudilin, index $m$ | $[1..2m]^2$ | $-4$ | $11.09017$ | $3.5556$ | $-1.8484$ | $-1.6116$ | $-1.0018$ | $-8.389$ | $\mathbf{-4.195}$ |
| Zudilin $3$-section | $[1..6n]^2$ | $-12$ | $1364.00073$ | $10.6667$ | $-5.5452$ | $-7.6074$ | $-6.9976$ | $-31.189$ | $-15.594$ |
| Nesterenko $(4,7)$ | $[1..6n]^2$ | $-14$ | $1496.54629$ | $10.6667$ | $-6.4694$ | $-8.4388$ | $-7.8291$ | $-33.776$ | $-16.888$ |
| **Zudilin $\odot$ Nesterenko, $W$** | $[1..6n]^2$ | $-2$ | (two slits) | $10.6667$ | $-0.9242$ | $-11.1017$ | n/a | $-33.794$ | $\mathbf{-16.897}$ |

Read this carefully — it is the answer to "what does the Hadamard construction cost":

1. **The Hadamard product costs essentially nothing in $\tau$.**  $W$ has the *same*
   type $[1..6n]^2$ as each aligned row.  The naive "denominators multiply" would have
   given $[1..6n]^4$; it does not, because $Q_{3n}$ and $B_n$ have purely $2$-power
   denominators.  ($B_Z\odot B_N$ and the doubly-conditional function *do* pay $k=4$.)
2. **The Hadamard product buys a lot at $2$**: $\varsigma_2$ improves from $-12$/$-14$ to
   $-2$ ($\gamma_2$ from $-5.5$/$-6.5$ to $-0.9$), which is the $24n$ bridge in
   function-theoretic dress — worth $+5.5$ of entry over the Nesterenko row.
3. **But it loses exactly as much archimedean-ly**: the domain shrinks from
   $|c|\approx1400$ to a two-slit domain of conformal radius $1.63$, because $W$ removes
   only *one* of the four product singularities.  Net: $-16.897$ against $-16.888$.
   **The Hadamard host is a wash against the worse of its two rows.**
4. **The real damage was done before the Hadamard product**, by the alignment: the
   $3$-section of Zudilin's row triples the $\mathrm{lcm}$ rate ($b:2\to6$) and the
   geometric $2$-power ($4m\to12n$) while the Koebe/Landau constants do not scale, costing
   $-4.195\to-15.594$.  The un-sectioned Zudilin row is by far the best of these objects.

---

## 8. The adjacent pair $(j_0,j_0+1)$ of `POSITIVITY_PROGRAM.md` §4.3 **[measured]**

Fixed rule (no selection bias) $m=3n$, $j_0=\mathrm{round}(0.30\,m)$; the two rows are the
Beukers moments $\mathrm{mom}(3n,j_0)$ and $\mathrm{mom}(3n,j_0+1)$ (`14_pairs.gp`,
exact rationals, $n\le24$; the $O(n^3)$ big-rational solve is the limit).
$w^{\rm pair}_n:=A^{(1)}_nB^{(2)}_n-A^{(2)}_nB^{(1)}_n$.

| quantity | value |
|---|---|
| type of $A^{(i)},B^{(i)}$ | $[1..6n]^2$, $\varsigma_2\approx-13.4$ |
| **type of $w^{\rm pair}_n$** | $\mathbf{[1..6n]^1}$ for all $5\le n\le24$ — $\sigma_m=6$, **half** the $Z\times N$ value |
| **$\varsigma_2(w^{\rm pair})$** | $\mathbf 0$: $v_2(w^{\rm pair}_n)=O(1)$ (measured $+4$ at $n=24$), i.e. the pair's $2$-adic bridge is *complete* |
| $\tfrac1n\log|w^{\rm pair}_n|$ | $-0.0266,-0.0139,-0.0022,+0.0088,+0.0189$ at $n=20,\dots,24$; still rising, extrapolated limit $\approx+0.26$ (range $[0.02,0.40]$) |

So the adjacent-pair host has **$\tau^\flat=6-\frac69=5.3333$ and $\gamma_2=0$** — a
$5.33$-nat improvement in $\tau$ and a $0.92$-nat improvement at $2$ over the
$Z\times N$ host.  With $d=e^{-\text{growth}}$:

| growth | $d$ | Koebe $\log(4d)$ | entry $\le$ | Landau $\log(16d)$ | entry $\le$ |
|---|---|---|---|---|---|
| $+0.02$ | $0.9802$ | $+1.3663$ | $-3.967$ | $+2.7526$ | $\mathbf{-2.581}$ |
| $+0.26$ | $0.7711$ | $+1.1263$ | $-4.207$ | $+2.5126$ | $-2.821$ |
| $+0.40$ | $0.6703$ | $+0.9863$ | $-4.347$ | $+2.3726$ | $-2.961$ |

**The entry test still fails**, by $2.6$–$4.3$; the corresponding deficit (margin$/(m-1)$) is
$-5.25$ to $-5.63$ at the Landau ceiling and $-6.63$ to $-7.01$ with the Koebe map.  Caveats, stated plainly: the growth limit is
pre-asymptotic and the singular set of this host was **not** computed (too few terms to
fit a recurrence), so $\log(4d)$ is used as an *upper* bound for the true conformal
radius (the actual domain has at least two slits, hence a smaller $\rho$); and
$\mathbf Q(x)$-independence was not checked. **[open]**

Still, the direction is informative and matches `POSITIVITY_PROGRAM.md`: the interior
adjacent pairs are genuinely arithmetically better than the
$\{$Zudilin$\}\times\{$Nesterenko$\}$ corner, by $\approx10$ nats of entry — and the
gain comes from the denominators ($k:2\to1$), exactly as `CATALAN_EXPLICIT.md` §5
reported for the lattice covolume.

---

## 9. Comparison with the modular hosts

In the per-function units of `CDT_NONCONGRUENCE.md` §5.3 (deficit $=$ margin$/(m-1)$
in (K)/(D), margin$/13$ in (S)):

| host / architecture | deficit | status |
|---|---|---|
| Beukers row, (D) | $+0.978$ | a proof |
| CDT's own $\sqrt{s_7}$ host, (S) | $+0.0004$ | a proof |
| **Catalan, level $8/16$, (S)$+$adelic at the ceiling** | $\mathbf{-0.613}$ | the project's best |
| Catalan, level $8$, (S) with the **measured** slopes (`ADELIC_HOLONOMY.md` §7.3) | $-1.725$ | corrected |
| Zudilin's hypergeometric row alone, (D) | $-4.195$ | this note |
| adjacent pair $(j_0,j_0{+}1)$ Hadamard host, (K) | $-5.3$ (Landau) to $-7.0$ (Koebe) | this note, **[measured]** |
| Zudilin $3$-section alone, (D) | $-15.594$ | this note |
| Nesterenko $(4,7)$ alone, (D) | $-16.888$ | this note |
| **Zudilin $\odot$ Nesterenko Hadamard host, $W$, (K)** | $\mathbf{-16.897}$ | this note |

**The Hadamard host is $16.3$ per function worse than the modular Catalan host.**  The
answer to `CDT_UNPACKED.md` §5's open question — *can one bad conditional function be
carried by unconditional functions of slope $24$?* — is: **no, and not for the reason
anticipated.**  The $2$-adic slope $24$ is real and does exactly what the adelic
theorem says it should ($\gamma_2$ improves from $-6.5$ to $-0.9$), but

1. it is paid for one-for-one at the archimedean place (the domain collapses from
   $|c|\approx1400$ to $\rho=1.63$), and
2. the binding constraint was never the conditional function: it is
   $\tau^\flat=10.67$, i.e. $\mathrm{lcm}(1,\dots,6n)^2$ in the hypergeometric rows
   against $\mathrm{lcm}(1,\dots,n)^2$ in the modular rows, together with the
   Galois-irreducibility of the Hadamard singular set, which deletes the pure module
   that is the only device CDT have for pushing $\tau^\flat$ down.

`CDT_UNPACKED.md` §7's closing sentence stands, unchanged and now with the two-row
route measured rather than guessed: *Catalan needs a host for the $\chi_{-4}$ class with
$|t_2|>e^2/16$, and the two-host Hadamard object does not supply one.*

---

## 10. Honest ledger

**[proved]**
* $|\varphi'(0)|=4ab/(a+b)$ for $\mathbf C\setminus((-\infty,-a]\cup[b,\infty))$; univalence.
* Uniqueness of $W$ in the $4$-dimensional Hadamard module (§4), modulo $1,G,G^2$ being $\mathbf Q$-independent.
* $\mathrm{Li}_j\odot\mathrm{Li}_k=\mathrm{Li}_{j+k}$ at the product point, and the Galois obstruction to a $\mathbf Q$-rational pure module (§5).
* Lemma G equivalence of the rescaled and adelic accountings (re-verified numerically here).

**[verified]** (exact arithmetic, stated range)
* All rows against the published integer rows, $n\le98$; rows to $n=400$.
* $v_2(h_n)\ge24n+14$ for $n\le98$; $h_n=2^{e(3n)+14n}D_{6n}^4w_n$ exactly.
* Minimal recurrences $(2,22)$, $(2,18)$, $(4,96)$ with $0$ failures on the full range including held-out indices; $P_{\rm gf}$ irreducible; roots $=$ the four products.
* All denominator types and $2$-adic slopes of §3, $n\le400$.
* $\mathbf Q(x)$-independence to degree $5$, mod $2^{61}-1$.
* Every entry/margin number of §§6–7 (from the CDT-calibrated `adelic_bound.py`).

**[measured]** (numerical, not converged)
* Archimedean growth limits (agree with the predicted $\varphi,T$ values to $3$–$4$ digits at $n=400$; the convergence is $O(1/n)$).
* Everything in §8 (adjacent pair, $n\le24$, pre-asymptotic).
* The rearrangement numerators of §6.1.

**[open]**
* A proof of $v_2(h_n)\ge24n-O(\log n)$ (unchanged from `CATALAN_AUDIT.md`).
* The exact differential order of the Hadamard module ($\le32$, $\ge4$).
* $\mathbf Q(x)$-independence as a theorem (CDT's Lemma 12.1.1 hypothesis).
* The singular set of the adjacent-pair host.

**No irrationality claim is made, and none is possible from this note: the entry
condition of the arithmetic holonomy bound fails on every configuration examined.**
