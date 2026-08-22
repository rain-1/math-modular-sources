# CDT_FINDER — scoring modular Apéry hosts against the Calegari–Dimitrov–Tang
# arithmetic holonomy bound

*Scripts: `lattice/cdt_finder/` (`cdt_bound.py`, `conformal.py`, `hosts.py`,
`indep_check2.py`, `x15_sym2.py`, `x15_margin.py`, `final_table.py`).
Literature: F. Calegari, V. Dimitrov, Y. Tang, **The linear independence of
$1,\zeta(2),L(2,\chi_{-3})$**, arXiv:2408.15403v2 (218 pp.; LaTeX source read, not OCR)
and **Arithmetic holonomy bounds and effective Diophantine approximation**,
arXiv:2510.04156. All equation/theorem numbers below are theirs.*

**No theorem is claimed here.** Everything is either (i) a reproduction of CDT's own
numbers, (ii) an exact computation on this project's rows, or (iii) an explicitly
flagged *estimate* obtained by transporting CDT's contour and function inventory to
another host. §8 is the honest computed/estimated ledger.

---

## 1. The exact inequality

**Theorem 2.5.1 (basic form) / Theorem 6.0.2 (fine form) / Theorem 7.0.1
(Bost–Charles form).** Let $m,r\in\mathbb N_{>0}$, $\mathbf e=(e_1,\dots,e_m)\in\mathbb N^m$,
and let $\mathbf b=(b_{i,j})$ be an $m\times r$ array of non-negative reals whose
columns have the *step shape*
$$0=b_{1,j}=\dots=b_{u_j,j}<b_{u_j+1,j}=\dots=b_{m,j}=:b_j,\qquad u_j\in\{0,\dots,m\}.$$
Put $\sigma_i:=b_{i,1}+\dots+b_{i,r}$ and
$$\boxed{\ \tau^{\flat}(\mathbf b):=\frac1{m^2}\sum_{i=1}^m(2i-1)\sigma_i
 =\sigma_m-\frac1{m^2}\sum_{j=1}^r u_j^2 b_j\in[0,\sigma_m]\ }\tag{6.0.4}$$
$$\boxed{\ \tau^{\sharp}(\mathbf e):=\frac2{m^2}\min_{\xi\in[0,m]}
 \Bigl\{\xi\sum_{i=1}^m e_i+\bigl(\max_i e_i\bigr)\,I_\xi^m(\xi)\Bigr\},\qquad
 \tau(\mathbf b;\mathbf e):=\tau^{\flat}(\mathbf b)+\tau^{\sharp}(\mathbf e)\ }\tag{6.0.5–6}$$
with (Definition 6.0.1) for $0\le\max\{u,1\}\le v$, $w\le v$:
$$I_u^v(w)=\int_{\min(u,1)}^1\!\!\max\{t-w,0\}\,dt
+\int_{\max(u,1)}^v\Bigl\{\!\!\sum_{h=1}^{\lfloor (t-1)/\max(1,w)\rfloor}\!\!\tfrac1h\Bigr\}dt
+\int_{\max(u,1)}^v\max\Bigl\{\tfrac{t}{\lfloor (t+\max(0,w-1))/\max(1,w)\rfloor}-w,0\Bigr\}dt .$$

Suppose $f_1,\dots,f_m\in\mathbb Q[\![x]\!]$ are $\mathbb Q(x)$-linearly independent,
**holonomic**, of denominator type
$$f_i(x)=a_{i,0}+\sum_{n\ge1}a_{i,n}\,\frac{x^n}{n^{e_i}\,[1,\dots,b_{i,1}n]\cdots[1,\dots,b_{i,r}n]},
\qquad a_{i,n}\in\mathbb Z, \tag{6.0.9}$$
and let $\varphi:(\mathbb D,0)\to(\mathbb C,0)$ be holomorphic with
$$\textbf{ENTRY CONDITION}\qquad \log|\varphi'(0)|>\tau(\mathbf b;\mathbf e)$$
and $\varphi^*f_i$ meromorphic on $\mathbb D$ for every $i$. Then
$$\boxed{\ m\;\le\;
\frac{\displaystyle\iint_{\mathbb T^2}\log|\varphi(z)-\varphi(w)|\,\mu_{\mathrm{Haar}}(z)\,\mu_{\mathrm{Haar}}(w)}
{\log|\varphi'(0)|-\tau(\mathbf b;\mathbf e)}\ }\tag{7.0.1 = BCbound}$$
The numerator may be replaced by the rearrangement integral
$\int_0^1 2t\,(\log|\varphi(e^{2\pi it})|)^*dt$ (6.0.15), and improved by the convexity
corrections of §7 (multi-radius $\hat T(r,\varphi)$) and §8.

We write
$$\boxed{\ \textbf{margin}\;:=\;m\bigl(\log|\varphi'(0)|-\tau(\mathbf b;\mathbf e)\bigr)-\mathrm{BC}(\varphi)\ }$$
so that **margin > 0** means the bound is violated, i.e. the hypothesised rational
relation is contradicted. If the $f_i$ are *not* assumed holonomic a priori the entry
condition must be strengthened to $\log|\varphi'(0)|>\max\{\sigma_m,\tau\}$, and then
holonomicity is a *conclusion* (Corollary 2.6.1).

Univalent special case: if $\varphi$ is the Riemann map of a contractible
$\Omega\ni0$ of conformal radius $\rho$, then $\mathrm{BC}=\log\rho$ and
$m\le\log\rho/(\log\rho-\tau)$.

### The five inputs
1. **the holonomic module** — $m$ functions, $\mathbb Q(x)$-linearly independent, with a
   common ODE/level (independence is a *hypothesis of the theorem*, to be proved
   separately: CDT's Lemma 12.1.1);
2. **denominators** — the array $\mathbf b$ (the rates $b_{i,j}$ with $[1,\dots,b_{i,j}n]$)
   *and* the integration vector $\mathbf e$ (the $n^{e_i}$);
3. **the domain** — one $\Omega$ *common to all $m$ functions*, and $\varphi$ with
   $\varphi^*f_i\in\mathcal M(\mathbb D)$; only $|\varphi'(0)|$ and the boundary
   log-energy of $\varphi$ enter;
4. **the integration profile** — $\mathbf e$ again, through $\tau^{\sharp}$ only;
5. **the dimension count** $m$.

---

## 2. Calibration on CDT's own $L(2,\chi_{-3})$ proof — all inputs reproduced

Everything in this block is *recomputed from scratch* and agrees with CDT to the last
printed digit (`cdt_bound.py`, `indep_check2.py`).

| input | CDT (§11–13, App. A) | reproduced here |
|---|---|---|
| host, coordinate | $Y_0(6)$, $x=q\prod\frac{(1-q^n)^4(1-q^{6n})^8}{(1-q^{2n})^8(1-q^{3n})^4}$, cusps $\mapsto x=0,\tfrac19,1,\infty$ | — (cited) |
| conditional ODE | $x(1-x)(1-9x)y''+(1-20x+27x^2)y'+3(3x-1)y=b+\tfrac{c}{1-x}$ (Prop. 11.1.4) | solved: $H_A=1{+}3x{+}15x^2{+}93x^3$, $H_B=x{+}\tfrac{23}4x^2$, $H_C=x{+}6x^2{+}\tfrac{343}9x^3$ ✅ |
| normaliser descent | $w(x)=\frac{x}{x-1}$, $y=x+w(x)=\frac{x^2}{x-1}$; $G=\mathrm{Sym}^+H$ | ✅ (`to_y`); $B_4=\mathrm{Sym}^-\mathrm{Li}_2=-4y+\tfrac49y^2+\tfrac{31}{900}y^3$ ✅ |
| singularities | $H$ on $\mathbf P^1\!\setminus\!\{0,\tfrac19,1,\infty\}$; $G$ on $\mathbf P^1\!\setminus\!\{0,4,\infty,-\tfrac1{72}\}$, holomorphic on $\mathbb C\!\setminus\![4,\infty)$, $\mathbb Z/2$ monodromy at $4$ | ✅ ($x=\tfrac19\mapsto y=-\tfrac1{72}$) |
| $m$ | $14=7$ pure $+\,7$ conditional | ✅ |
| $\mathbf b$ | $14\times2$, $u_1=1,u_2=3$, $b_1=b_2=2$, $\sigma_m=4$ | ✅ |
| $\mathbf e$ | $(0,0,1;0^6;1,1,1,1,1)$, $\sum e_i=6$, $\max e_i=1$ | ✅ |
| $\tau^{\flat}$ | $191/49=3.897959$ | **$191/49$** ✅ |
| $I_2^{14}(2)$ | (implicit) | **$21.075$ exactly** ✅ |
| $\tau^{\sharp}$ | $27/80=0.3375$, minimiser $\xi\in[2,13/6]$ | **$0.3375$ at $\xi\approx2.16$** ✅ |
| $\tau(\mathbf b;\mathbf e)$ | $16603/3920=4.235459$ | ✅ |
| $\varphi$ | $\varphi=h\circ\psi$, $h=\lambda+\frac{\lambda}{\lambda-1}=-256\,\Delta(2\tau)/\Delta(\tau)$ the $\Gamma_0(2)$ hauptmodul; $\psi:\mathbb D\to\Omega\subset\mathbb D$ a 4-slit+lune Riemann map (Lemma A.4.4) | — (cited) |
| $|\psi'(0)|$ | $\tfrac{5448339453535586608\cdot10^9}{8658833407565631122430056127}=0.6292232680$ | ✅ |
| $|\varphi'(0)|$ | $256\cdot0.62922=161.0812$, $\log=5.081908$ | ✅ |
| $\mathrm{BC}$ integral | $11.845$ | — (cited; not recomputed) |
| **bound** | $m\le 13.9938<14$ (A.5.1) | **13.9938** ✅ |
| signed margin | — | $14\cdot0.846449-11.845=\mathbf{+0.0053}$ |

Improved numerators CDT also give: $13.730$ (Thm 6.0.2, $l=1$), $13.7206$, $13.678$,
and best $13.621$ (four radii) → best margin $\approx+0.32$.

**Independence check (their Lemma 12.1.1), verified numerically here.** All 14
functions were rebuilt exactly over $\mathbb Q$ — $G$ from a generic relation
$(a,b,c)=(1,-3,5)$, and $B_1=1$, $B_2,B_3,B_5,B_6$ from their factorial series,
$B_4=\mathrm{Sym}^-\mathrm{Li}_2$, $B_7=\int B_4\,dy/y$ — and the rank of
$\{y^jf_i\}_{i\le14,\,j\le D}$ was computed mod $p=2^{61}-1$:

| $\deg P_i\le$ | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| rank / needed | 14/14 | 28/28 | 42/42 | 56/56 | 70/70 | 84/84 |

i.e. **no $\mathbb Q(y)$-relation of degree $\le5$**; series order up to $y^{90}$.
(`indep_check2.py`.) This is a check of, not a substitute for, CDT's proof.

---

## 3. What actually varies from host to host — the reduction to three numbers

CDT's architecture, read structurally:

* the conditional function $H$ lives on $\mathbf P^1\setminus\{0,t_1,t_2,\infty\}$
  where $t_1,t_2$ are the singular $t$-values of the Picard–Fuchs operator
  (reciprocals of the characteristic roots $\lambda_1,\lambda_2$ of the Apéry
  recurrence). CDT's $x$ has $t_1=\tfrac19$, $t_2=1$;
* the **pure** module is the polylogarithm module on $\mathbf P^1\setminus\{0,s,\infty\}$
  where $s$ = the **outer** singularity $=1/\lambda_2$. Integrality of the pure
  functions $\mathrm{Li}_j(x/s)=\sum(\lambda_2x)^n/n^j$ forces $\lambda_2$ to be an
  **algebraic integer** (CDT: $\lambda_2=1$);
* the **normaliser descent** is the involution of $\mathbf P^1\setminus\{0,s,\infty\}$
  fixing $0$ and swapping $s\leftrightarrow\infty$, i.e. $w(x)=\frac{sx}{x-s}$,
  $y=x+w(x)=\frac{x^2}{x-s}$; the fixed point $x=2s$ becomes the $\mathbb Z/2$ branch
  point $y=4s$, and the inner singularity $t_1$ becomes the one *extra* point
  $y(t_1)$ that $\varphi$ must have a single preimage of. Denominators double
  ($[1..n]^k\to[1..2n]^k$, i.e. $b_j=2$);
* **uniformisation.** $\mathbf P^1\setminus\{0,s,\infty\}$ is uniformised by
  $x=s\,\lambda(z)$, $\lambda=16z-128z^2+\dots$; the symmetrised coordinate is then
  $y=-256\,s\,z^2+\dots$. Hence
  $$\boxed{\ |\varphi'(0)|\le 256\,s=256/\lambda_2\ }$$
  (CDT: $256$). This ceiling is attained only up to the loss forced by deleting the
  preimages of $y(t_1)$; CDT realise the factor $0.62922$.
* since the ambient orbifold is the **same** for every host and only the scale $s$
  changes, a contour of CDT's shape rescales by $\varphi\mapsto s\varphi$, giving
  $$\log|\varphi'(0)|=5.081908+\log s,\qquad \mathrm{BC}=11.845+\log s .$$

So a host contributes exactly three numbers: **$\lambda_2$** (through $s=1/\lambda_2$),
**$k$** (the sharp Eichler denominator exponent, $=w+1$ unless a free integration
drops it), and the **pure inventory** (which fixes $m$ and the $u_j$).

**Number fields.** When $\lambda_2\notin\mathbb Q$ (a conjugate pair, $\lambda_1\lambda_2=c$)
the construction lives over $K=\mathbb Q(\sqrt{\;})$ and there are two archimedean
places, with $s$ and $\bar s$. The LCM denominators are *rational* integers, so their
local contributions occur at every place of $K$; the correctly normalised (per-degree)
accounting therefore uses
$$\log|\varphi'(0)|\;\rightsquigarrow\;\tfrac1{[K:\mathbb Q]}\sum_{v\mid\infty}d_v\log|\varphi'(0)|_v
\;=\;\log\bigl(256\,|N(s)|^{1/[K:\mathbb Q]}\bigr),$$
with $\tau(\mathbf b;\mathbf e)$ unchanged. This is exactly the bookkeeping in Beukers'
own 1987 proof of Theorem 4, where the product of the two conjugate linear forms is
weighed against $[1,\dots,n]^{6}=[1,\dots,n]^{2k}$ — the LCM cost appearing once per
place. We therefore define
$$\lambda_2^{\mathrm{norm}}:=|N(\lambda_2)|^{1/[K:\mathbb Q]}=|c|^{1/2}\ \text{for a conjugate pair},
\qquad \text{budget}=\log(256/\lambda_2^{\mathrm{norm}}).$$
*(The unnormalised alternative — summing $\log|\varphi'(0)|_v$ over places against a
single $\tau$ — is reported alongside in §6; we believe it to be wrong, but CDT state
the plain bound only over $\mathbb Q$, so this is the one genuinely uncertain input.)*

---

## 4. Ranked table

`final_table.py`. Columns: $\lambda_2^{\mathrm{norm}}$; the field the pure module lives
over; $\tau=\tau^{\flat}+\tau^{\sharp}$; **ceil** $=\log(256/\lambda_2^{\mathrm{norm}})$
(the hard uniformisation ceiling); **entryC/entryR** $=$ ceil$-\tau$ resp.
(ceil$+\log0.62922)-\tau$ (CDT's realised contour loss transported); **margin**
$=m(\text{entryR})-\mathrm{BC}$ with $\mathrm{BC}=11.845+\log s$.

### A. CDT-identical architecture ($m=14$, $p_0=7$, CDT's pure inventory $u=(1,3,\dots)$)

| host | lvl | $w$ | $k$ | $\lambda_2^{\rm norm}$ | field | $\tau$ | ceil | entryC | entryR | **margin** | period |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **Zagier C (10,3,9)** | 6 | 1 | 2 | 1 | $\mathbb Q$ | 4.235 | 5.545 | +1.310 | +0.846 | **+0.005** | $L(2,\chi_{-3})/2$ ← **CDT** |
| Zagier A (7,2,−8) | 6 | 1 | 2 | 1 | $\mathbb Q$ | 4.235 | 5.545 | +1.310 | +0.846 | **+0.005** | $\zeta(2)/4$ |
| Cooper $s_7$ | 7 | 2 | 2 | 1 | $\mathbb Q$ | 4.235 | 5.545 | +1.310 | +0.846 | **+0.005** | $\zeta(2)/7$ (free integration) |
| Zagier D (11,3,−1) | 11 | 1 | 2 | 1 | $\mathbb Q(\sqrt5)$ unit | 4.235 | 5.545 | +1.310 | +0.846 | **+0.005** | $\zeta(2)/5$ |
| scan2 #942–4 | 16 | 1 | 2 | 2 | $\mathbb Q$ | 4.235 | 4.852 | +0.617 | +0.153 | −9.01 | ? |
| **Zagier E (12,4,32)** | 8 | 1 | 2 | 4 | $\mathbb Q$ | 4.235 | 4.159 | **−0.077** | −0.540 | −18.02 | **Catalan $G$** |
| Cooper $s_{10}$ | 10 | 2 | 2 | 4 | $\mathbb Q$ | 4.235 | 4.159 | −0.077 | −0.540 | −18.02 | $\zeta(2)/5$ |
| Apéry (17,5,1) | 5/6 | 2 | 3 | 1 | $\mathbb Q(\sqrt2)$ unit | 5.980 | 5.545 | −0.435 | −0.898 | −24.42 | $\zeta(3)/6$ |
| **$X_1(5)\,\mathrm{Sym}^2$** | 5 | 2 | 3 | 1 | $\mathbb Q(\sqrt5)$ unit | 5.980 | 5.545 | −0.435 | −0.898 | −24.42 | $8\zeta(3)-5\sqrt5L(3,\chi_5)$ |
| Zagier F (17,6,72) | 6/8 | 1 | 2 | 8 | $\mathbb Q$ | 4.235 | 3.466 | −0.770 | −1.233 | −27.03 | $\tfrac58L(2,\chi_{-3})$ |
| Cooper $s_{18}$ | 18 | 2 | 2 | 12 | $\mathbb Q$ | 4.235 | 3.060 | −1.175 | −1.638 | −32.30 | $L(2,\chi_{-3})/2$ |
| $T$ (12,4,16) | 6 | 2 | 3 | 4 | $\mathbb Q(\sqrt2)$ | 5.980 | 4.159 | −1.821 | −2.285 | −42.45 | $7\zeta(3)/32$ |
| Domb (10,4,64) | 6 | 2 | 3 | 4 | $\mathbb Q$ | 5.980 | 4.159 | −1.821 | −2.285 | −42.45 | $7\zeta(3)/24$ |
| AZ (9,3,−27) | 9 | 2 | 3 | 5.196 | $\mathbb Q(\sqrt3)$ | 5.980 | 3.897 | −2.083 | −2.546 | −45.85 | $L(3,\chi_{-3})/3$ |
| $\zeta(7)$ level 24 | 24 | 6 | 7 | 4.828 | $\mathbb Q(\sqrt2)$ | 9.695 | 3.971 | −5.724 | −6.187 | −96.89 | $\tfrac{1463}{13824}\zeta(7)$ |

**Reading.** With CDT's own function inventory the *only* positive margins are the four
hosts with $\lambda_2^{\rm norm}=1$ — the Apéry-perfect / Fricke-palindromic hosts — and
they all give the *same* margin $+0.005$ as CDT's, because the architecture is exactly
scale-covariant and CDT sit at $\lambda_2=1$. Three of the four carry $\zeta(2)$ (known
irrational); the fourth is CDT's own theorem. **The finder produces no new
$(\text{host},\text{weight})$ with positive margin in this architecture.** Catalan
(Zagier E) is the nearest miss and misses the *entry test alone* by $0.0765$ at the
hard ceiling.

### B–C. Best conceivable pure inventory (whole pure orbit denominator-free, $u_j=m/2$)

$\tau^{\flat}\to\tfrac34\sigma_m$ ($=3$ for $k=2$, $4.5$ for $k=3$).

| host | $k$ | $\tau$ | entryR | margin @ $m=14$ | margin @ $m=50$ |
|---|---|---|---|---|---|
| $\lambda_2^{\rm norm}=1$ hosts | 2 | 3.34 | +1.744 | +12.6 | +85.0 |
| Zagier E (Catalan) | 2 | 3.34 | +0.358 | −5.4 | **+17.1** (positive from $m\approx30$) |
| Apéry / $X_1(5)\,\mathrm{Sym}^2$ | 3 | 4.84 | +0.244 | −8.4 | **+10.0** (positive from $m\approx49$) |
| $T$, Domb | 3 | 4.84 | −1.14 | — | — |
| AZ(9,3,−27) → $L(3,\chi_{-3})$ | 3 | 4.84 | −1.40 | — | — |
| $\zeta(7)$ level 24 | 7 | 10.84 | −7.33 | — | — |

So the single lever that matters is **how many pure functions of the host's denominator
type are denominator-free** (CDT have exactly one, $B_1=1$, hence $u_1=1$).

---

## 5. Eisenstein extension classes and the periods they carry

For $\Phi=P(V)E_k^{\psi,\varphi}$ on $\Gamma_0(N)$ (Theorem A / Theorem B\*),
$L(\Phi,k-1)=P(k-1)L(\psi,k-1)L(\varphi,0)$; the two orientations are

* **inner** ($\varphi=\mathbf 1$, $\psi=\chi$ odd for $k=3$): period $-\tfrac12P(k-1)L(k-1,\chi)$;
* **outer** ($\psi=\mathbf 1$, $\varphi=\chi$): period $P(k-1)\zeta(k-1)L(\chi,0)$, non-zero
  only for $\chi$ odd (so $L(\varphi,0)\ne0$).

Hence at **$k=3$ ($w=1$)** the available period pairs on $\Gamma_0(N)$ are
$\{\zeta(2)\}\cup\{L(2,\chi_{-D}):D\mid N\}$; at **$k=4$ ($w=2$)** they are
$\{\zeta(3)\}\cup\{L(3,\chi_{-D})\}$ (both inducing characters odd), and the
even-character directions carry period $0$. CDT retain both orientations
$(1-8V_2)E_{3,\chi_{-3},1}$ and $(1-V_2)E_{3,1,\chi_{-3}}$ deliberately
(`book/v8/08a`, `book/v11`): full purification would leave a single period and only
prove irrationality, not linear independence.

For the priority host, $\dim M_4^{\mathrm{Eis}}(\Gamma_1(5))=4$: the pairs $(\mathbf1,\mathbf1)$
at $d=1,5$ and $(\chi_5,\mathbf1)$, $(\mathbf1,\chi_5)$ at $d=1$ (the two quartic
characters mod 5 are **odd**, hence absent in weight 4). Period vector at $s=3$:
$$\bigl(\zeta(3)\zeta(0),\ \zeta(3)\zeta(0),\ L(3,\chi_5)\zeta(0),\ \zeta(3)L(0,\chi_5)=0\bigr)
\;=\;\bigl(-\tfrac12\zeta(3),\,-\tfrac12\zeta(3),\,-\tfrac12L(3,\chi_5),\,0\bigr),$$
matching Beukers' $L(F,s)=10(1-5^{2-s})\zeta(s)\zeta(s-3)+\zeta(s)L(s-3,\chi)-5\sqrt5\,\zeta(s-3)L(s,\chi)$.
The Galois trace $F+\bar F$ is defined over $\mathbb Q$ but **kills the
$L(3,\chi_5)$ coordinate** (its coefficient is $\mp5\sqrt5$); therefore any construction
reaching $L(3,\chi_5)$ on this host is intrinsically over $\mathbb Q(\sqrt5)$.

---

## 6. Priority case: $X_1(5)$, $\mathrm{Sym}^2$ (Beukers 1987 Theorem 4)

Target: *"$1,\zeta(3),L(3,\chi_5)$ are linearly independent over $\mathbb Q(\sqrt5)$"*
(which would contain $L(3,\chi_5)\notin\mathbb Q$), by CDT-amplifying Beukers'
$8\zeta(3)-5\sqrt5L(3,\chi_5)\notin\mathbb Q(\sqrt5)$.

**Computed exactly** (`x15_sym2.py`, `x15_margin.py`):

* weight-1 layer = Zagier D, $(n{+}1)^2a_{n+1}=(11n^2{+}11n{+}3)a_n+n^2a_{n-1}$,
  $a=1,3,19,147,1251,\dots$; characteristic roots $\tfrac{11\pm5\sqrt5}{2}=\varphi^5,-\varphi^{-5}$;
* **singular $t$-values $t_1=\varphi^{-5}=0.0901699\ldots$, $t_2=-\varphi^{5}=-11.0901699\ldots$,
  $t_1t_2=-1$** — Apéry-perfect fold, and $N_{\mathbb Q(\sqrt5)/\mathbb Q}(t_2)=-1$: $t_2$ is a **unit**;
* $\mathrm{Sym}^2$ row $A_n=[t^n]F(t)^2=1,6,47,408,3745,35598,346583,3433776,\dots\in\mathbb Z$;
* **no free integration**: $(n+1)\nmid A_n$ (fails first at $n=2$);
* minimal recurrence: 5-term, cubic coefficients
  $P_0=(n{+}2)^3$, $P_1=22n^3{+}165n^2{+}419n{+}360$, $P_2=119n^3{+}1071n^2{+}3222n{+}3240$,
  $P_3=-(22n^3{+}231n^2{+}815n{+}966)$, $P_4=(n{+}4)^3$ — verified to $n=295$; characteristic
  roots $\varphi^5,\varphi^5,-\varphi^{-5},-\varphi^{-5}$, so the singular set is unchanged;
* **denominator exponent $k=3$, sharp**, verified $d_n^3B_n\in\mathbb Z$ and $d_n^2B_n\notin\mathbb Z$
  for the census companion ($b_0{=}0,b_1{=}1$) at every $n\le139$.

**Geometry** (`conformal.py`, exact via the modular $\lambda$ function):

| construction | at $v_1$ ($\sqrt5\to+2.236$) | at $v_2$ ($\sqrt5\to-2.236$) | normalised |
|---|---|---|---|
| univalent slit plane $\mathbb C\setminus(-\infty,t_2]$, $\rho=4|t_2|$ | $44.361$ ($\log=3.792$) | $0.361$ ($\log=-1.020$) | $\log 4=1.386$ |
| thrice-punctured $\mathbb C\setminus\{t_1,t_2\}$ (exact) | $1.3630$ | $1.3630$ | $\log=0.310$ |
| modular ceiling, unsymmetrised $16|s|$ | $177.44$ ($5.179$) | $1.443$ ($0.366$) | $\log16=2.773$ |
| **modular ceiling, symmetrised $256|s|$** | $2839.08$ ($7.951$) | $23.08$ ($3.139$) | $\log256=\mathbf{5.5452}$ |

(the normalised column is exactly $\log 4/16/256$ because $|N(t_2)|=1$; the involution
fixing the fold is $w(x)=\frac{t_2x}{x-t_2}$ with quotient coordinate $y=\frac{x^2}{x-t_2}$,
branch point $y=4t_2$, extra singularity $y(t_1)=\frac{t_1^2}{t_1-t_2}=-4.06\cdot10^{-4}$.)

**Verdict.**

| accounting | budget | $\tau$ ($k=3$, CDT-proportional inventory) | entry | margin ($m=14$) |
|---|---|---|---|---|
| normalised, ceiling | 5.5452 | 5.980 | **−0.435** | — |
| normalised, CDT contour loss | 5.0819 | 5.980 | **−0.898** | −24.4 |
| normalised, best inventory $u_j=m/2$ | 5.0819 | 4.837 | **+0.244** | −8.4 ($m=14$); **+10.0** at $m=50$ |
| *unnormalised* (sum over places) | 10.164 | 5.980 | +4.183 | +34.9 |

So: **the $X_1(5)\ \mathrm{Sym}^2$ host is exactly as good analytically as CDT's level-6
host** ($\lambda_2^{\rm norm}=1$ in both cases — this is the Apéry-perfect condition
$t_1t_2=\pm1$ read arithmetically as $|N(t_2)|=1$), **but it costs one more Eichler
integration ($k=3$ instead of $2$), and that extra $\approx1.74$ of $\tau$ is exactly what
the budget cannot absorb.** Under the CDT-proportional pure inventory it fails the
entry test by $0.435$ even at the hard uniformisation ceiling.

The gap is quantitative and specific: entry needs $\tau(\mathbf b;\mathbf e)<5.082$,
i.e. $\frac2{m^2}\sum_ju_j^2\ge1.24$, i.e. $\sum_j(u_j/m)^2\ge0.62$ — roughly *46 % of the
functions must be free of all three LCM layers*, against CDT's $(1/14,3/14,\cdot)$.
And a positive margin then needs $m\gtrsim49$, i.e. $p_0\gtrsim25$ per orbit.

---

## 7. Candidates, and what would have to be proved

**No $(\text{host},\text{weight})$ other than CDT's own has a positive margin under the
architecture that CDT actually used.** The three other $\lambda_2^{\rm norm}=1$ hosts
tie CDT's $+0.005$ but their periods are $\zeta(2)$ (already irrational). The
substantive candidates, all *conditional on a richer pure inventory*, are:

**(C1) Catalan $G$ on $\Gamma_0(8)$ (Zagier E), $k=2$.** Misses the entry test by
$0.0765$ at the ceiling, $0.540$ with CDT's contour. Would need: (a) a pure module on
$\mathbf P^1\setminus\{0,\tfrac14,\infty\}$ of type $[1..2n]^2$ with $u_1/m,u_2/m$
raised from CDT's $(0.07,0.21)$ to $\gtrsim(0.3,0.3)$ — i.e. $\gtrsim6$ *denominator-free*
pure functions where CDT have one; (b) $m\approx30$ functions with proved
$\mathbb Q(y)$-independence; (c) sharp $d_n^2$ arithmetic (known: $D_n^2B_n\in\mathbb Z$
sharp); (d) a contour with $|\varphi'(0)|>e^{\tau}\cdot$, i.e. losing less than CDT's
$0.629$ against the ceiling $64$.

**(C2) $X_1(5)\ \mathrm{Sym}^2$ over $\mathbb Q(\sqrt5)$ → $1,\zeta(3),L(3,\chi_5)$
independent over $\mathbb Q(\sqrt5)$.** Missing pieces, in order of difficulty:
1. *the number-field holonomy bound* — CDT state the plain bound (Thm 2.5.1/7.0.1) over
   $\mathbb Q$ only; arXiv:2510.04156 Thm 2.5 is adelic but is an approximation-measure
   statement. A clean "$m\le\mathrm{BC}/(\overline{\log|\varphi'(0)|}-\tau)$ over $K$" with
   the normalisation argued in §3 is needed. **This is the single load-bearing unproved
   input.**
2. *a pure inventory with $\gtrsim46\%$ denominator-free functions* of type $[1..2n]^3$ on
   $\mathbf P^1\setminus\{0,4t_2,\infty\}$, and $p_0\gtrsim25$;
3. *$\mathbb Q(\sqrt5)(y)$-linear independence* of the resulting $m\approx50$ functions;
4. *the conditional construction*: from $a+b\zeta(3)+cL(3,\chi_5)=0$ ($a,b,c\in\mathbb Q(\sqrt5)$)
   build $H$ and its conjugate $\bar H$ with $d_n^3$ arithmetic in $\mathbb Z[\varphi]$
   (Beukers proves $[1,\dots,n]^3c_n$ is an algebraic integer — this input is *available*);
5. *the contour*: a $\psi$ on the $\Gamma_0(2)$-uniformised $y$-disc avoiding the preimages
   of $y(t_1)=-4.06\cdot10^{-4}$, at both places.

**(C3) $\zeta(3)$ via Apéry's own host, $k=3$, $\lambda_2^{\rm norm}=1$** — numerically
identical to (C2) and equally short; it is the natural *control* for (C2), since the
answer there is known.

Everything with $\lambda_2^{\rm norm}\ge4$ and $k\ge3$ ($T$, Domb, AZ$(9,3,-27)\to L(3,\chi_{-3})$,
$\zeta(7)$ level 24) is out of reach by margins of $-26$ to $-370$ and is not a candidate.

---

## 8. Computed vs. estimated — the honest ledger

**Computed exactly (this task).** $\tau^\flat,\tau^\sharp,I_u^v(w)$ and CDT's
$191/49$, $21.075$, $27/80$, $16603/3920$, $|\varphi'(0)|=161.081$, bound $13.9938$;
CDT's $H_A,H_B,H_C$ from their ODE; the symmetrisation map and $B_4$'s published series;
the rank/independence table for CDT's 14 functions; the $X_1(5)$ $\mathrm{Sym}^2$ row,
its minimal recurrence, its singular values, the absence of a free integration, and the
sharp $k=3$; all conformal radii (thrice-punctured sphere exactly via $\lambda$;
`conformal.py` reproduces the archive's level-10 $\tau_0$ to 22 digits and $r_u=1.30651$).

**Cited, not recomputed.** CDT's Bost–Charles integral $11.845$ and its convexity
improvements; the $\Gamma_0(6)$ hauptmodul and cusp values; Beukers' $L(F,s)$ and his
$[1,\dots,n]^3c_n$ integrality; the Zagier/Cooper census rows $(\lambda_1,\lambda_2,k)$
(from `paper/tables/census.tex` and `lattice/sporadic_scan2/table.json`).

**Estimated — flagged everywhere above.**
1. *Transport of CDT's contour.* We assume a host of the same shape admits a contour with
   the same $0.62922$ loss and the same $\mathrm{BC}-\log|\varphi'(0)|=6.763$ shape term.
   Untested; the extra singularity $y(t_1)$ sits at a different hyperbolic depth for each
   host, so this could move either way by several tenths.
2. *The pure inventory* $(m,p_0,u_j,\mathbf e)$ for every host other than CDT's. Scenario A
   transports CDT's $(14,7,(1,3),6)$ literally; scenarios B–C are upper-bound fantasies.
   Nothing here identifies actual pure functions on any other host.
3. *The number-field normalisation* (§3). Argued from Beukers' own $[1..n]^{2k}$ accounting
   and from the fact that LCMs are rational integers; not proved, and the two conventions
   differ by a factor 2 in the budget — enough to flip (C2) from $-0.90$ to $+4.18$.
4. *$m=2p_0$.* Follows the coordinator's rule (one hypothesised relation ⇒ one conditional
   function ⇒ one conditional orbit) and CDT's $7+7$; not derived.

**Not done.** No contour was designed or its Bost–Charles integral computed for any host
other than CDT's; no pure module was constructed on any host other than CDT's; no
independence check was run for any candidate; the Eisenstein census was carried out
character-theoretically (§5), not by an $N\le60$ machine sweep of $M_3/M_4$ on every
genus-zero $\Gamma_0(N)$ and Atkin–Lehner quotient.
