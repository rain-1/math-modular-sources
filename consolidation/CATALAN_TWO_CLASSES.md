# Two fold-regular $\chi_{-4}$ classes on one host, and the doubly-small function they produce

*Fable (Opus 5), 2026-08-22.  Scripts and logs: `lattice/catalan_two_classes/`
(`00_setup.gp`, `02_classes.gp`, `03_level16.gp`, `04_symmetrise.gp`,
`05_margin.py`, `06_level32src.gp`, `07_ceiling.py`, `08_canon.gp`; logs
`02run.log`–`08run.log`).  Answers the "next task" of
`consolidation/ADELIC_HOLONOMY.md` §4.3/§6.3.  Conventions from
`consolidation/CDT_FINDER.md`, `consolidation/THEOREM_B_EXACT.md`,
`consolidation/EULER_CRITERION.md`.  Tags **[verified]** = exact or numerical
computation in this task; **[estimate]** = transported input.
**No irrationality claim is made anywhere.**

---


> **Retraction (2026-08-22, `CATALAN_THREE_PERIOD.md`).** §2.4 and §4 below are wrong: the outer orientation *does* contain an unconditional doubly-small class, $\Phi_0^{\rm out}=(1+3V_2-4V_4)T$. §4 tested outer classes against Theorem B*'s period formula, which does not hold at level 16; the correct outer period functional is $\xi^{\rm out}_\infty=\tfrac32\zeta(2)(P(1)-P(2))$ with fold-regularity $\iff P(0)=0$. The doubly-small orbit has dimension 2. The margins of §6 also used an incorrect slope floor in the $y$-coordinate; see `CATALAN_THREE_PERIOD.md` for the corrected table.

## 0. Verdict

| claim | verdict |
|---|---|
| **A genus-zero host carrying two independent $\chi_{-4}$ classes with fold-regular companions exists: it is $X_0(16)$ with $x=\eta(2\tau)\eta(16\tau)^2/(\eta(\tau)^2\eta(8\tau))$, $F_{16}=\eta(2\tau)^4/\eta(4\tau)^2$.** The two classes are the host's own canonical source $E=E_3^{\chi_{-4},\mathbf1}$ ($\xi_\infty=-G/2$) and $V_2(1-8V_2)E$ ($\xi_\infty=+G/8$) | **[verified]** §2 |
| The fold-regular subspace of $\langle E,V_2E,V_4E\rangle$ is **exactly** $\{c_1E+c_2V_2(1-8V_2)E\}$, cut out by the single linear condition $8c_2+c_4=0$ | **[verified]** §2.3 |
| Hence a **target-zero fold-regular class** $\Phi_0=(1+4V_2-32V_4)E$ exists, and its companion $B_0$ is an **unconditional doubly-small function**: rational coefficients, $d_n^2B_0\in\mathbf Z$ (sharp), archimedean radius $\sqrt2/4$ against the host's $1/4$, $2$-adic slope $3/2$ in $x$ (measured $1.484$ at $n=128$), $1$ in the symmetrised coordinate | **[verified]** §3 |
| **This is the object `ADELIC_HOLONOMY.md` §4.3 proved does not exist at level 8.** The level-8 obstruction was dimensional: $P(2)=P(3)=0$ is impossible for $P$ of degree $1$ in $2^{-s}$; at level 16 the degree is $2$ and both conditions can be met | **[verified]** §1.3 |
| The doubly-small orbit has dimension **exactly $1$**: level-32 sources ($V_8E$) are *not* fold-regular on this host, and the outer orientation $E_3^{\mathbf1,\chi_{-4}}$ contributes no target-zero fold-regular direction | **[verified]** §2.4, §4 |
| **The archimedean geometry does not improve.** The only normaliser involution preserving the host's singular set is $\sigma(x)=-x/(4x+1)$, i.e. $s=-1/4$ — the **fold**, not the outer singularity. Ceiling $\log(256|s|)=4.158883$, *identical to level 8*; entry and margins therefore repeat level 8's, improved only by the one extra function | **[verified]** §5 |
| Adelic margin at the ceiling: $-7.966$ (level 8, $m=14$) $\to$ $\mathbf{-7.551}$ ($m=15$ with the doubly-small function). Reaching $0$ would need $\approx14$ doubly-small functions; the Eisenstein space supplies one | **[verified]** §6 |
| **New, quantified obstruction.** In the symmetrised coordinate the pure module reaches $|y|=1$ but the conditional/doubly-small functions only $|y|=1/2$: the extra singularity sits at relative radius $1/2$, against $1/8$ on the level-8 host and $1/288$ for CDT. Transporting CDT's contour loss $0.62922$ is not defensible here | **[verified]** §5.3 |
| Level 144 is **not** a genus-zero host (genus 13; the $W_{144}$ quotient model has genus 5), has no explicit Picard–Fuchs operator or $(\lambda_1,\lambda_2)$, and its own CDT entry test fails $1.0708$ vs $1.1506$ | **[verified from the archive]** §7 |

One sentence: *the two-class host exists and is level 16; it produces exactly one
unconditional doubly-small function, of $2$-adic slope $3/2$; and it buys $+0.42$
of margin against a deficit of $-8$, because the archimedean side of the level-16
geometry is no better than level 8's and its contour is demonstrably worse.*

---

## 1. The criterion, sharpened

### 1.1 What "fold-regular" is, operationally

For a host $(\Gamma,t,F)$ with $A(t)=F$, $B_\Psi(t)=$ ($t$-expansion of $F\cdot\mathcal D^{-2}\Psi$),
write $\lambda_1$ for the reciprocal of the fold (first singularity) and $\lambda_2$
for the next.  Three inequivalent conditions must be kept apart:

1. $b_n/a_n$ **converges**;
2. $b_n/a_n\to\xi$ and $B-\xi A$ has **radius enlarged to $1/\lambda_2$**
   (*fold-regular*: trivial local monodromy at the fold);
3. $\xi=0$ (*target zero*).

(1) does **not** imply (2).  Measured on the level-16 host, the outer class
$T-V_2T$ has $b_n/a_n\to\tfrac34\cdot\tfrac{\zeta(2)}2$ to $20$ digits at $n=168$
yet $|c_n|^{1/n}=3.7735\to4$ — the log survives with a factor $(1-\lambda_1t)\log(1-\lambda_1t)$,
which only improves $t_c^{-n}/n$ to $t_c^{-n}/n^2$.  Fold-regularity is the
statement that $H_\xi=F(\Theta-\xi)$ is *invariant* under the parabolic
$\gamma_1$ generating the stabiliser of the fold cusp:
$$H_\xi\circ\gamma_1=H_\xi\iff r_{\gamma_1}(\tau)=\xi\,(j(\gamma_1,\tau)-1),$$
$r_{\gamma_1}=\alpha+\beta\tau$ the Eichler period polynomial of $\Phi$ at that
cusp.  This is **one** linear condition $\alpha(\Phi)=0$ on $\Phi$ (then $\xi$ is
determined by $\beta$).  Hence:

> **Counting rule.**  On a host whose admissible $\chi_{-4}$ Eisenstein space has
> dimension $d$, the fold-regular subspace has dimension $\ge d-1$.  Two
> fold-regular classes therefore require $d\ge3$, i.e. **level $\ge16$**.

### 1.2 Why level 8 cannot work (the level-8 obstruction, restated)

On the level-8 host the fold $t=1/8$ *is* the cusp $0$, so
Theorem B\* applies verbatim: for $\Phi=P(V)E_3^{\chi_{-4},\mathbf1}$,
fold-regularity $\iff$ endpoint annihilation $P(3)=0$, and then
$\xi_\infty=L(\Phi,2)=-\tfrac12P(2)G$.  With $P(s)=c_1+c_22^{-s}$ the two
conditions $P(3)=0$ and $P(2)=0$ have only the zero solution.  The measured
divergence rate on the level-8 host is affine in $P(3)$ — `ADELIC_HOLONOMY.md`'s
"$0.0518\,(8-c)$" is exactly $P(3)=(8-c)/8$ — which identifies the mechanism.

### 1.3 Level-16 sources on the level-8 host: killed by ramification

$t=x(1+2x)/(1+4x)^2$ **[verified exactly to $q^{40}$]**, and $1-8t=(1+4x)^{-2}$,
so $t=1/8$ is the cusp $x=\infty$ with **ramification index 2**: near the fold the
level-16 coordinate is $1/x\sim4\sqrt{1-8t}$.  A level-16 source therefore
acquires a $\sqrt{1-8t}$ branch there.  Measured on the level-8 host:

| class | $P(3)$ | $P(2)$ | $b_n/a_n$ | $|c_n|^{1/n}$, $n=108$ |
|---|---|---|---|---|
| $(1-8V_2)E$ | $0$ | $-1$ | $\to G/2$ (18 digits) | $3.797\to4$ ✅ |
| $(1-8V_2)(1-4V_2)E$ | $0$ | $0$ | $\to0$ like $n^{-1/2}$ | $7.409\to8$ ❌ |
| $V_2(1-8V_2)E$ | $0$ | $-1/4$ | $\to G/8$ like $n^{-1/2}$ | $7.315\to8$ ❌ |

The $n^{-1/2}$ decay is the signature of the square-root fold.  **The host must
itself be level 16.**

---

## 2. The level-16 host

### 2.1 Data (all rebuilt from the eta quotients, `00_setup.gp`)

$$x=\frac{\eta(2\tau)\eta(16\tau)^2}{\eta(\tau)^2\eta(8\tau)}=q+2q^2+4q^3+8q^4+14q^5+24q^6+\cdots,
\qquad F_{16}=\frac{\eta(2\tau)^4}{\eta(4\tau)^2}=1-4q^2+4q^4-\cdots$$
$$\boxed{F_{16}\,\mathcal D_qx=E:=E_3^{\chi_{-4},\mathbf 1}}\qquad\textbf{[verified to }q^{40}\textbf{]}$$
so the host's **canonical source is the bare Eisenstein series** ($P\equiv1$) — unlike
level 8, where it is $(1-8V_2)E$.
$$A(x)=1-4x^2+16x^3-44x^4+96x^5-176x^6+320x^7-860x^8+3520x^9-\cdots\in\mathbf Z[[x]].$$

**Singular set** (the six cusps of $\Gamma_0(16)$, images under $t$):
$$x=0,\quad x=-\tfrac12\ (t=0),\quad x=\infty\ (t=\tfrac18),\quad x=-\tfrac14\ (t=\infty),\quad x=\tfrac{-1\pm i}4\ (t=\tfrac14).$$
Ordering by modulus: **fold $x_1=-1/4$** ($\lambda_1=4$; measured $|a_n|^{1/n}=3.769$
at $n=168$), **next $x_2=(-1\pm i)/4$**, $|x_2|=\sqrt2/4$, $\lambda_2^{\rm norm}=2\sqrt2$
(measured $|c_n|^{1/n}=2.75\to2\sqrt2$), then $-1/2$, then $\infty$.

Note the *inversion* relative to level 8: on the level-16 host the fold is the cusp
over $t=\infty$, and the level-8 fold $t=1/8$ has moved to $x=\infty$.

### 2.2 The two fold-regular classes

$$\boxed{\Phi_1=E,\quad \xi_\infty=-\tfrac12G;\qquad
\Phi_2=V_2(1-8V_2)E=(V_2-8V_4)E,\quad \xi_\infty=+\tfrac18G.}$$

`03run.log`, $N=170$ terms, $n=168$:

| class | $b_n/a_n$ | predicted $-\tfrac12P(2)G$ | $|c_n|^{1/n}$ |
|---|---|---|---|
| $\Phi_1=E$ | $-0.45798279708860950752730175765$ | $-0.45798279708860950752730175747$ | $2.577$ |
| $\Phi_2=(V_2-8V_4)E$ | $\ \ 0.11449569927215237688183115$ | $\ \ 0.11449569927215237688182544$ | $2.741$ |
| $3\Phi_1+5\Phi_2$ | $-0.80146989490506663817274951$ | $-0.80146989490506663817277808$ | $2.767$ |

$19$–$20$ correct digits, and $|c_n|^{1/n}$ climbing to $2\sqrt2=2.8284$ (the
next singularity), against $4$ for every non-fold-regular class.  So
$\xi_\infty=L(\Phi,2)=-\tfrac12P(2)G$ continues to hold on this host even though
the fold is **not** the cusp $0$ — the classes are fold-regular *and* their
periods are the same critical values.

### 2.3 The fold-regular subspace, exactly

Writing $\Phi=(c_1+c_2V_2+c_4V_4)E$:
$$\boxed{\ \Phi\ \text{fold-regular}\iff 8c_2+c_4=0\iff \Phi\in\langle E,\ V_2(1-8V_2)E\rangle\ }$$
**[verified]**, `04run.log`: eight test classes, every $\ell=8c_2+c_4=0$ case has
$b_n/a_n$ matching $-\tfrac12P(2)G$ to $\ge17$ digits and $|c_n|^{1/n}\approx2.72$;
every $\ell\ne0$ case has the wrong limit and $|c_n|^{1/n}\approx3.70\to4$.
Structurally: *the fold-regular classes are the host's own canonical source
together with $V_2$ applied to the level-8 canonical source.*

$$\xi_\infty(c)=-\tfrac12\Bigl(c_1+\tfrac{c_2}4+\tfrac{c_4}{16}\Bigr)G
\ \stackrel{c_4=-8c_2}{=}\ -\tfrac12\Bigl(c_1-\tfrac{c_2}4\Bigr)G .$$

### 2.4 No third direction

Adding $V_8E$ (level-32 sources) gives nothing: `06run.log` shows
$(V_4-8V_8)E$ and $V_2\Phi_0$ both have $|c_n|^{1/n}=3.57,3.69\to4$ and limits
inconsistent with $-\tfrac12P(2)G$.  $X_0(32)\to X_0(16)$ ramifies over the fold.
$X_0(32)$ has genus 1 and $X_0(64)$ genus 3, so the genus-zero ladder stops at 16.

---

## 3. The doubly-small function

$$\boxed{\ \Phi_0=\Phi_1+4\Phi_2=(1+4V_2-32V_4)E,\qquad \xi_\infty=0,\ \xi_2=0\ }$$

Equivalently, in `ADELIC_HOLONOMY.md` §3's notation with $r_1=-\tfrac12$,
$r_2=+\tfrac18$: $\ r_2B_1-r_1B_2=\tfrac18(B_1+4B_2)=\tfrac18B_{\Phi_0}$.
$\xi_2=0$ is *forced* by Theorem F: $\chi_{-4}(2)=0$ so $\mathcal E_2\equiv1$,
$Q=P$ and $\xi_2=-P(2)\kappa_2$ — the same rational functional $P(2)$ that
kills the archimedean period.  This is Conjecture D doing exactly the work
§3 of `ADELIC_HOLONOMY.md` predicted it would.

$$B_{\Phi_0}(x)=x-\tfrac{64}9x^3+\tfrac{80}3x^4-\tfrac{12416}{225}x^5+\tfrac{1792}{45}x^6
+\tfrac{2103296}{11025}x^7-\tfrac{160768}{175}x^8+\tfrac{201539584}{99225}x^9-\cdots$$

Measured (`03run.log`, `04run.log`), $N=170$:

| property | measurement | reading |
|---|---|---|
| rational coefficients | exact | unconditional ✅ |
| $b_n/a_n$ | $2.29\cdot10^{-23}$ at $n=168$ | $\xi_\infty=0$ ✅ |
| $|b_n|^{1/n}$ | $2.688,\,2.732,\,2.763$ at $n=60,100,168$ | $\to2\sqrt2$: radius $\sqrt2/4>1/4$ ✅ |
| $v_2(b_n)$ | $88,\,148,\,250$ at $n=60,100,168$ | slope exactly $3/2$ ($\Delta v_2/\Delta n=1.5$ on both gaps) ✅ |
| $d_n^kB_0\in\mathbf Z$ | $k=2$ yes, $k=1$ fails at $n=3$ | $k=2$ **sharp**, same type as level 8 ✅ |

For comparison the pure module $\mathrm{Li}_j(-4x)$ on this host has $2$-adic
slope $2=v_2(\lambda_1)$ (measured $1.96$ at $n=100$).

**In the symmetrised coordinate** $y=x^2/(x-s)=4x^2/(4x+1)$, $s=-1/4$
(`04_symmetrise.gp`, exact series composition through $x(v)=\tfrac12(v^2+v\sqrt{1+v^2})$,
$y=v^2$; the odd $v$-part of $\mathrm{Sym}^+$ vanishes identically, as it must):

| function | $v_2(d_k)/k$ | $|d_k|^{1/k}$ | reading |
|---|---|---|---|
| $\mathrm{Sym}^+B_{\Phi_0}$ | $0.95,\,0.975,\,1.079$ ($k=20,40,63$) | $1.80,\,1.87,\,1.91$ | slope $\varsigma_y=1$, singular at $|y|=1/2$ |
| $\mathrm{Sym}^+\mathrm{Li}_2(-4x)$ | $1.70,\,1.80,\,1.90$ | $0.82,\,0.88,\,0.91$ | slope $\varsigma_y=2$, singular at $|y|=1$ |

The rule $\varsigma_y=2(\varsigma_x-1)$ (for $v_2(s)=-2$) reproduces
`ADELIC_HOLONOMY.md`'s level-8 measurement $\varsigma_x=\varsigma_y=2$ for the
pure orbit, and gives $\varsigma_y=1$ here.

---

## 4. The outer orientation contributes nothing

$E_3^{\mathbf1,\chi_{-4}}$ has $\xi_2=0$ always (Theorem F, $\varphi\ne\mathbf1$),
so a fold-regular outer class with $\xi_\infty=0$ would be a second doubly-small
direction.  It does not occur.  On the level-16 host (`03run.log`):
$T$, $V_2T$, $V_4T$, $V_2T-8V_4T$, $T-5V_2T+4V_4T$ all have
$|c_n|^{1/n}\to4$; the one class whose $b_n/a_n$ converges to $20$ digits,
$T-V_2T$ ($\to\tfrac38\zeta(2)$), is *not* fold-regular in the sense of §1.1
(radius stays $1/4$); and $T-5V_2T+4V_4T$, which has $P(0)=P(2)=0$ and would be
target-zero if Theorem B\* applied, instead converges to $-\tfrac34\zeta(2)$.
The $\pi^2$ orientation stays $\pi^2$.

---

## 5. Geometry: the normaliser descent, and where the level-16 host loses

### 5.1 The involution

The singular set $\{-\tfrac14,-\tfrac12,\tfrac{-1\pm i}4,\infty\}$ is preserved by
exactly one CDT-type involution fixing $0$:
$$\sigma(x)=\frac{sx}{x-s}=\frac{-x}{4x+1},\quad s=-\tfrac14;\qquad
y=x+\sigma(x)=\frac{x^2}{x-s}=\frac{4x^2}{4x+1}.$$
$\sigma$ swaps $-\tfrac14\leftrightarrow\infty$ and $\tfrac{-1+i}4\leftrightarrow\tfrac{-1-i}4$,
and fixes $0$ and $-\tfrac12$.  Images: $y(0)=0$ and $y(-\tfrac12)=-1=4s$ are the two
$\mathbf Z/2$ **branch points**; $y=\infty$ carries the fold and the cusp $0$;
$y(\tfrac{-1\pm i}4)=-\tfrac12$ is the one **extra singularity**.
(The Fricke involution $x\mapsto1/(8x)$ of the level-16 notes is *not* usable here:
it does not fix $0$, it glues the base point to the cusp $x=\infty$.)

$$\textbf{ceiling}\quad |\varphi'(0)|\le256|s|=64,\qquad \log=4.158883$$
— numerically **identical to the level-8 host**, because $|s|=1/4$ in both.

### 5.2 The inversion of roles

CDT symmetrise about the *outer* singularity $s=1/\lambda_2$ and the fold becomes
the extra point.  Here the only admissible $s$ is $1/\lambda_1$ — the **fold** —
and the *outer* pair becomes the extra point.  Consequence: the pure module
$\mathrm{Li}_j(x/s)=\mathrm{Li}_j(-4x)$ is singular at the fold, hence at $y=4s$;
the conditional and doubly-small functions are singular at $y=-1/2$.

### 5.3 The new, quantified obstruction

$$\frac{|y_{\rm extra}|}{|4s|}:\qquad \text{CDT level 6}=\tfrac1{288}=0.00347,\quad
\text{level 8}=\tfrac18=0.125,\quad \boxed{\text{level 16}=\tfrac12}$$
**[verified]**, and confirmed independently by the measured $y$-radii of §3:
the pure module reaches $|y|=1$, the conditional/doubly-small only $|y|=1/2$.
CDT's realised contour keeps $0.62922$ of the ceiling with the extra point at
relative radius $1/288$; transporting that number to relative radius $1/2$ is
not defensible, and `CDT_FINDER.md` §8 estimate 1 was already flagged as the
weakest input.

*If* $\varphi$ must avoid $y=-\tfrac12$ (the singularity there is the logarithmic
one of a cusp, not a pole), then the thrice-punctured bound applies exactly:
with $\Sigma=\{-\tfrac12,-1,\infty\}$, base $0$, the Möbius $M(y)=-2y-1$ sends
$\Sigma\mapsto\{0,1,\infty\}$ and $0\mapsto-1=\lambda(1+i)$, so
$$|\varphi'(0)|\ \le\ \tfrac12\cdot2\,\mathrm{Im}(1+i)\,|\lambda'(1+i)|
=\frac{\pi^2}{\Gamma(3/4)^4}=4.376879,\qquad \log=1.476336,$$
**[verified]** `07_ceiling.py` (using $\lambda'=i\pi\lambda(1-\lambda)\theta_3^4$,
$\lambda(1+i)=-1$), giving entry $=-2.759$ and no bound at all.  Whether CDT's
architecture tolerates a single preimage of the extra point (their Lemma A.4.4
slit map suggests it does) is **the** decisive open input, exactly as at level 8.

---

## 6. The adelic margin on the level-16 host

`05_margin.py`, reusing `lattice/adelic_holonomy/adelic_bound.py`.
$\mathrm{BC}=11.845+\log(1/4)=10.458706$; ceiling $4.158883$; "realised"
$=$ ceiling $+\log0.62922=3.695614$ **[estimate, see §5.3]**.

| accounting | $m$ | $\tau$ | $\gamma_2$ | entry | margin |
|---|---|---|---|---|---|
| archimedean only, ceiling | 14 | 4.2355 | 0 | $-0.0766$ | $-11.531$ |
| adelic (pure $\varsigma_y=2$), ceiling | 14 | 4.2355 | $+0.2546$ | $+0.1780$ | $-7.966$ |
| **adelic $+$ the doubly-small function ($\varsigma_y=1$), ceiling** | **15** | **4.2269** | $\mathbf{+0.2619}$ | $\mathbf{+0.1938}$ | $\mathbf{-7.551}$ |
| same, realised contour | 15 | 4.2269 | $+0.2619$ | $-0.2694$ | $-14.500$ |
| hypothetical $d$ doubly-small: $d=4$ | 18 | 4.2029 | $+0.2910$ | $+0.2469$ | $-6.014$ |
| hypothetical $d=14$ | 28 | 4.1464 | $+0.3855$ | $+0.3980$ | $+0.685$ |

So the doubly-small function is worth $+0.415$ of margin, and one would need
$\approx14$ of them to reach $0$ — against a supply of exactly $1$.  The
`ADELIC_HOLONOMY.md` §4.3 dimension cap is not binding here (at $\varsigma_y=1$
it permits many more than one); the binding constraint is the **dimension of the
$\chi_{-4}$ Eisenstein oldspace**, which is $3$ at level 16, minus one condition
for fold-regularity, minus one for target zero.

The inventory lever remains an order of magnitude larger: with the
"best conceivable" pure inventory ($u_1=u_2=m/2$) the same host gives
$+9.33$ at $m=14$ and $+63.5$ at $m=50$ — and nothing in this task exhibits a
single denominator-free pure function on this host.

---

## 7. Level 144, tested against the same criteria

Extracted from `packages/apery systems 3/LEVEL_144_CATALAN_SYSTEM_COMPLETE_DOSSIER_v2.docx`
and the ChatGPT archive (the `LEVEL144_..._NEXT_STEP_PACKAGE.zip` contains none
of the mathematics).

* **It is not a genus-zero host.**  $g(X_0(144))=13$; the continuation surface is
  $X_0(144)$ minus 20 bad cusp classes; the $W_{144}$-quotient model
  $\mathbf Q(t,x)$, $\deg_xP=8$, has genus 5.  Global monodromy of the degree-8
  cover is $S_8$; the "$S_3$" is a local suborbit.  There is no
  $(\lambda_1,\lambda_2)$, no Picard–Fuchs operator over $\mathbf Q(t)$, and no
  $a_n/b_n$ recurrence — all three are open tasks in the source notebook.
* **Its source already satisfies our conditions at the Mellin level.**
  $\Phi_{144}=g(2\tau)-30g(6\tau)+81g(18\tau)$, $g=E_3^{\chi_{-4},\mathbf1}$ in
  Yang's normalisation, has $P(s)=2^{-s}(1-3^{1-s})(1-3^{3-s})$, so
  $P(1)=P(3)=0$, $P(2)=-1/3$, $L(\Phi_{144},2)=-G/6$.  $P(3)=0$ is the endpoint
  condition; $P(2)\neq0$, so this is a *period-carrying* class, not a target-zero
  one.  The five-dimensional "zero-target correction space" on the $d\mid36$
  scales (basis $v_1,\dots,v_5$ in the dossier) is cut out by $P(1)=P(2)=0$ plus
  two odd-scale conditions — **none of which is the fold condition**, and the
  fold there is an interior branch point $\alpha\in Q_{18}$ (a root of an explicit
  degree-9 polynomial), $|\alpha|=0.4537$, not a cusp.  So the level-144 system
  has no "two fold-regular classes" statement to test in our sense.
* **Its own CDT entry test fails.**  $m=28$ gives $\tau=55/392$, requiring
  $|\varphi'(0)|>1.150626$; the hyperbolic comparison gives
  $|\varphi'(0)|\le1.070788$.  $\ge58$ rows would be needed merely to enter.

Level 144 is therefore *not* an alternative to level 16 for this question.

---

## 8. What is now open

1. **The contour.**  §5.3.  Whether CDT's architecture admits an extra
   singularity at relative radius $1/2$, and what $|\varphi'(0)|$ is actually
   realisable on $\mathbf P^1\setminus\{0,-\tfrac12,-1,\infty\}$ with the
   $\mathbf Z/2$ structure at $0$ and $-1$.  Until this is settled every
   "realised" row above is decoration.  **[open, decisive]**
2. **A pure inventory.**  No denominator-free pure function has been exhibited on
   any Catalan host.  This remains the only lever of the right size.  **[open]**
3. **Other coordinates on $X_0(16)$.**  The Hauptmodul is fixed only up to
   $x\mapsto x/(1+cx)$ (normalised $t=q+O(q^2)$); different $c$ move which cusp is
   the fold and change the singular configuration, hence the ceiling and possibly
   the fold-regular subspace, while keeping integrality for $c\in\mathbf Z$.
   Not scanned.  **[open]**
4. **$\mathbf Q(y)$-independence** of any candidate module on this host.  Not
   checked (as at level 8).  **[open]**
5. **A proof of the counting rule** of §1.1 (that fold-regularity is one linear
   condition, with the period read off the second coefficient of the period
   polynomial at the fold cusp).  Verified here on 20+ classes across two hosts;
   not proved.  **[open]**

**No irrationality claim.**  The one new unconditional object is
$B_{\Phi_0}\in\mathbf Q[[x]]$ of §3: archimedean radius $\sqrt2/4$ against the
host's $1/4$, $2$-adic slope $3/2$, denominator type $d_n^2$.  Nothing here bears
on the irrationality of $G$.
