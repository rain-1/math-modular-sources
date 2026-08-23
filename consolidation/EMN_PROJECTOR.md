# The EMN adelic projector: the monodromy calculus, the search, and an honest entry

*Fable, 2026-08-23.  Scripts: `lattice/emn_projector/`
(`emn_core.py`, `emn_span.py`, `emn_family.py`, `emn_cand.py`, `emn_table.py`,
`emn_slopes.py`, `emn_fold.py`, `emn_tau.py`, `emn_ledger.py`, `emn_opt.py`, `emn_bc.py`).
Source: the ChatGPT note "A Non-Eisenstein Catalan Programme from the
Eskandari–Murty–Nemoto Motive", §§2–5, 13, 17.3.  Framework: `ADELIC_HOLONOMY.md` §2.6,
`CDT_UNPACKED.md`, `CATALAN_OBSTRUCTION.md`, `CATALAN_MU4.md`, `SOL_NOTES_DIGEST.md` §2.*

**No theorem about $G$ is claimed.**  §§1–4 are proved or exactly computed; §5 is measured
to $n=1000$–$1300$; §§6–8 are a pricing in the standard units of the ledger.

---

## 0. Verdict in one line

The systematic search for §17.3's adelic projector **terminates, and returns nothing new**:
after the Niven collapse (§2) the whole shift lattice $(\pi/12)\mathbf Z$ degenerates to
$b\in\{0,\pm\pi/3\}$, and the cancellation systems (§4) have a *one-dimensional* optimal
solution space — the generator $H$ itself, with every $N\le24$ multiplication map giving back
$\pm3^{k}H+\text{const}$.  What the search **did** produce is a repair of the note's fatal
step: there **is** a fold, not at the Catalan point $z=1$ (where the log coefficient is
$\tfrac12$, $G$-free — the F3 failure recorded in `SOL_NOTES_DIGEST.md` §2.4) but at the
*other* puncture $z=2$ — the second preimage of the cusp of the $y$-line, at which
$\mathrm{Hf}$ is analytic with $\mathrm{Hf}(2)=2G$, an honest Apéry limit.  Priced honestly this gives

$$\boxed{\ \text{entry}=-0.85\ }\qquad(\text{crude }\tau:\ -0.89;\ \text{best CDT }\tau\text{ at }m=9\text{–}14:\ -0.85\text{ to }-0.89)$$

against benchmarks CDT's own $L(2,\chi_{-3})$ $+0.77$, level 8 symmetrised $-0.077$,
$\mu_4$ host A $-0.563$, level 16 $-1.46$.  So the moving EMN period is **not** the
counterfactual $+1.27$ host: it is a real CDT architecture with a real fold, and it is
worse than level 8 by $0.77$ nats.  **What binds is the conditional row's lcm weight
$k=3$** (measured true size $2.97$), not the conformal size; and $k$ is inflated by exactly
$\log 2$ because — uniquely among all hosts in this ledger — **the fold lies farther from
the expansion point than the surviving singularity does**.

---

## 1. The trigonometric calculus, closed

Write $\mathrm{ang}$ for the note's $S$ and $\mathrm{Hf}$ for its $H$ (PARI-safe names; the
scripts use `ang`, `Hfun`).  With $z=1-\cos(\mathrm{ang})$ and
$\mathscr H(\mathrm{ang})=\mathrm{Hf}(z)$,

$$\mathscr H(\mathrm{ang})=\sum_{m\ge0}\frac{(-1)^m}{(2m+1)^2}\bigl(1-\cos((2m+1)\mathrm{ang})\bigr).$$

**(1.1) The generator is a double antiderivative of a secant.**  Term-by-term,
$\sum_m(-1)^m\cos((2m+1)u)=\tfrac12\sec u$ (Abel), so

$$\boxed{\ \mathscr H''(\mathrm{ang})=\tfrac12\sec(\mathrm{ang}),\qquad \mathscr H(0)=\mathscr H'(0)=0\ }$$

and $\mathscr H'(\mathrm{ang})=\tfrac12\log\bigl|\tan(\tfrac\pi4+\tfrac{\mathrm{ang}}2)\bigr|$
(the note's formula, with the modulus that its §3 omits).

**(1.2) Clausen form.**  Since $\chi_{-4}(n)=\sin(n\pi/2)$ for all $n\ge1$,

$$\boxed{\ \mathscr H(\mathrm{ang})=G-\tfrac12\Bigl[\operatorname{Cl}_2\bigl(\tfrac\pi2+\mathrm{ang}\bigr)+\operatorname{Cl}_2\bigl(\tfrac\pi2-\mathrm{ang}\bigr)\Bigr].\ }$$

Equivalently, with $w=e^{i\,\mathrm{ang}}$, $\mathscr H=G-\tfrac12\operatorname{Im}\bigl[\operatorname{Li}_2(iw)+\operatorname{Li}_2(i/w)\bigr]$: **the EMN moving period is the
$\operatorname{Im}\operatorname{Li}_2$ module at $\mu_4$ of `CATALAN_MU4.md`, pushed to the
quotient $w\mapsto1/w$.**  That identification is the structural content of the whole note
and it is what makes every number below predictable from the $\mu_4$ ledger.

**(1.3) Special values (all verified to 30 digits, `emn_fold.py`).**
$\mathscr C(b):=\sum_{n\ge1}\chi_{-4}(n)\cos(nb)/n^2$, so $\mathscr H=G-\mathscr C$:
$$\mathscr C(0)=G,\quad \mathscr C(\pi/3)=\tfrac23G,\quad\mathscr C(\pi/2)=0,\quad
\mathscr C(2\pi/3)=-\tfrac23G,\quad\mathscr C(\pi)=-G,$$
giving $\mathrm{Hf}(0)=0$, $\mathrm{Hf}(\tfrac12)=\tfrac13G$, $\mathrm{Hf}(1)=G$,
$\mathrm{Hf}(\tfrac32)=\tfrac53G$, $\boxed{\mathrm{Hf}(2)=2G}$ — the note's §5 values, all
confirmed.  $\mathscr C(b)\in\mathbf QG$ **only** for $b\in\tfrac\pi3\mathbf Z\cup\tfrac\pi2\mathbf Z$;
e.g. $\mathscr C(\pi/6)=\tfrac{\sqrt3}2L(2,\chi_{12})\notin\mathbf QG$.

**(1.4) Two exact functional equations.**  $\mathscr H(u+\pi)=2G-\mathscr H(u)$ and
$\mathscr H(u+\tfrac\pi2)+\mathscr H(u-\tfrac\pi2)=2G$ (both immediate from (1.2)); the
first is the note's $\mathrm{Hf}(z)+\mathrm{Hf}(2-z)=2G$.

**(1.5) Master ODE.**  With $\Psi:=d^2F/d\,\mathrm{ang}^2$ expressed in $z$,
$$z(2-z)F''+(1-z)F'=\Psi(z),\qquad
n(2n-1)f_n=[z^{n-1}]\Psi+(n-1)^2f_{n-1}.$$
For $\mathrm{Hf}$, $\Psi=\tfrac1{2(1-z)}$; the recursion
$h_n=\bigl(\tfrac12+(n-1)^2h_{n-1}\bigr)/(n(2n-1))$ gives
$\tfrac12,\tfrac16,\tfrac7{90},\tfrac3{70},\dots$ and is what all the scripts run.
**$f_0$ is the only free constant**: $\Psi$ determines everything from $n\ge1$ on.

---

## 2. Rationality, and the Niven collapse of the shift lattice

Take $F=\sum_\nu c_\nu\mathscr H(a_\nu\mathrm{ang}+b_\nu)$, $a_\nu\in\mathbf Z_{>0}$,
$c_\nu\in\mathbf Q$.

**(2.1) Evenness.**  $z$-monodromy at $0$ sends $\mathrm{ang}\mapsto-\mathrm{ang}$, so $F$ is a
function of $z$ iff it is even; the $b\ne0$ terms must be paired,
$\mathfrak h_{a,b}:=\mathscr H(a\,\mathrm{ang}+b)+\mathscr H(a\,\mathrm{ang}-b)$.

**(2.2) The $\Psi$ of a pair.**  With $T_a$ the Chebyshev polynomial and $x=1-z=\cos(\mathrm{ang})$,
$$\mathfrak h_{a,b}''=\tfrac{a^2}2\bigl[\sec(a\,\mathrm{ang}+b)+\sec(a\,\mathrm{ang}-b)\bigr]
=\frac{a^2\cos b\;T_a(x)}{T_a(x)^2-\sin^2b}.$$
Hence $\Psi\in\mathbf Q(z)$ **iff $\cos b\in\mathbf Q$** (and then $\sin^2b=1-\cos^2b\in\mathbf Q$
automatically).  Different $a$ cannot conspire: the Galois group of
$\mathbf Q(\zeta_{24})^+=\mathbf Q(\sqrt2,\sqrt3)$ acts on each $a$-block separately (the
$T_a$ are defined over $\mathbf Q$ and the pole divisors of distinct $a$-blocks are distinct),
and on the four elements $b\in\{\pi/12,\pi/6,\pi/4,5\pi/12\}$ it acts by
$\sigma(\sqrt2\mapsto-\sqrt2)$: $P_{\pi/6}\mapsto P_{\pi/6}$, $P_{\pi/4},P_{\pi/12},P_{5\pi/12}\mapsto-$themselves;
$\tau(\sqrt3\mapsto-\sqrt3)$: $P_{\pi/6}\mapsto-P_{\pi/6}$, $P_{\pi/4}\mapsto P_{\pi/4}$,
$P_{\pi/12}\leftrightarrow-P_{5\pi/12}$.  The invariant subspace is $0$.

> **Theorem 1 (Niven collapse).**  *For $b$ a rational multiple of $\pi$, the pair
> $\mathfrak h_{a,b}$ has rational Taylor coefficients in $z$ iff $\cos b\in\{0,\pm\tfrac12,\pm1\}$
> (Niven's theorem), i.e. iff $b\in\tfrac\pi3\mathbf Z\cup\tfrac\pi2\mathbf Z$.  The case
> $b\equiv\tfrac\pi2$ is degenerate ($\Psi\equiv0$, by (1.4)), and $b\equiv\pi$ is $b\equiv0$
> up to sign and a constant.  Hence the entire shift lattice $(\pi/12)\mathbf Z$ — or
> $(2\pi/N)\mathbf Z$ for any $N$ — collapses to*
> $$b\in\{0\}\cup\{\pm\pi/3\}\pmod\pi .$$

So the search space is exactly two families,
$$g_{a,0}:=\mathscr H(a\,\mathrm{ang}),\qquad
g_{a,\pi/3}:=\mathfrak h_{a,\pi/3},$$
with $\Psi_{a,0}=\tfrac{a^2}{2T_a(1-z)}$ and $\Psi_{a,\pi/3}=\tfrac{a^2}2\cdot\frac{T_a(1-z)}{T_a(1-z)^2-3/4}$.
**No new shifts exist.**  (The $\pi/3$ shift is not new either: it is the cubic
distribution's own coset.)

**(2.3) Where $G$ can enter.**  $\mathfrak h_{a,b}(0)=2G-2\mathscr C(b)$, and by (1.3) this lies
in $\mathbf QG$ for every admissible $b$.  Since $\Psi\in\mathbf Q(z)$ fixes all coefficients
$n\ge1$, **the hypothesis $G\in\mathbf Q$ affects only $F(0)$** — a constant, already in every
module.  *Rationality of the coefficients is therefore never the conditional mechanism here.*
This is the precise form of the F3 objection, and §6 is how it is evaded.

---

## 3. The log-monodromy calculus

Singularities of $F$ in the $\mathrm{ang}$-line are the poles of $\Psi$, i.e.
$a\,\mathrm{ang}+b\equiv\tfrac\pi2\pmod\pi$.  Write $\mathrm{ang}=\pi\theta$, $\theta\in(0,1)$,
$z_0=1-\cos(\pi\theta)$.

> **Lemma (local monodromy).**  If $a\pi\theta\pm b=\tfrac\pi2+k\pi$ then that branch
> contributes to the $\mathrm{ang}$-residue of $\Psi$ the amount $\tfrac a2(-1)^{k+1}$, and
> $$\Delta_{z_0}F=2\pi i\,r(\theta)\,\bigl(\mathrm{ang}-\pi\theta\bigr)+O\bigl((z-z_0)^2\bigr),
> \qquad r(\theta)=\operatorname*{Res}_{\mathrm{ang}=\pi\theta}\Psi .$$
> The $z$-residue is $r(\theta)\sin(\pi\theta)$; both preimages $\pm\theta$ give the same
> $z$-pole.  **$F$ is regular at $z_0$ iff $r(\theta)=0$ — a $\mathbf Q$-linear condition on
> $(c_\nu)$**, and the conditions are automatically Galois-stable.

Two structural consequences, both used below.

* **$z=0$ and $z=2$ are never poles of $\Psi$.**  Evenness forces $r=-r$ at $w=\pm1$.  So
  $\mathrm{Hf}$ and every $F$ in the span are *analytic* at $z=2$, with
  $F(2)=\sum_\nu c_\nu\bigl(2G-2(-1)^{a_\nu}\mathscr C(b_\nu)\bigr)\in\mathbf QG$.
* **Parity.**  $T_a(-x)=(-1)^aT_a(x)$, so every basis $\Psi$ is an eigenvector of $x\mapsto-x$
  ($z\mapsto2-z$).  For an even (resp. odd) $\Psi$, $\operatorname{res}(-x_0)=-\operatorname{res}(x_0)$
  (resp. $+$).  Cancelling a pole at $x_0>0$ leaves its mirror at $-x_0$ alive with residue
  $\pm2\operatorname{res}_{\text{odd}}$: **poles can only be deleted in mirror pairs unless the
  two parities are mixed.**

Calibration of the Lemma against the note's §13: $K=\mathscr H(3\,\mathrm{ang})+3\mathscr H(\mathrm{ang})$
has, at $\theta=\tfrac12$, $r=\tfrac32(-1)^{1+1}+3\cdot\tfrac12(-1)^{0+1}=\tfrac32-\tfrac32=0$ —
the note's cancellation, reproduced.  Its remaining poles are
$\theta=\tfrac16,\tfrac56$, i.e. $z=1\mp\tfrac{\sqrt3}2$ — the note's obstruction, reproduced.

---

## 4. The linear systems, and the rigidity theorem

`emn_span.py` builds the exact rational residue matrix over the basis
$\{g_{a,0},g_{a,\pi/3}\}_{a\le A}$ and its pole set
$\theta\in\{(2k+1)/2a\}\cup\{(6k\pm1)/6a\}$, and sweeps the threshold $\theta^*$: cancel
every pole with $\theta<\theta^*$, then read off $\theta_{\min}$ of the survivors.  Since
$z_0=1-\cos(\pi\theta)$ is increasing, $\theta_{\min}$ is exactly the *outer singularity*
$t_2$ of the resulting host.

> **Theorem 2 (rigidity of the projector).**  *For $A\le24$ (all multiplication maps $N\le24$
> and both shift classes):*
> 1. $\max\theta_{\min}=\tfrac12$, i.e. $\ \boxed{t_2\le1}$ *for every nonconstant $F$ in the span;*
> 2. *the solution space with $t_2=1$ has dimension $9$ over $\mathbf Q$ ($A=24$; $5$ for
>    $A=12$), of which $6$ (resp. $2$) members are constants ($\Psi\equiv0$) and the
>    remaining $3$ are*
>    $$F_1=\mathscr H(\mathrm{ang}),\quad
>      F_3=\mathscr H(3\,\mathrm{ang})-3\mathfrak h_{1,\pi/3},\quad
>      F_9=\mathscr H(9\,\mathrm{ang})-3\mathfrak h_{3,\pi/3}+9\mathfrak h_{1,\pi/3};$$
> 3. *and these are not independent:* $F_3=-3F_1$, $F_9=9F_1$ *exactly, coefficient by
>    coefficient (`emn_cand.py`, $n\le600$).*
>
> *Hence the space of adelic projectors with the best possible geometry is
> $\mathbf Q\cdot\mathrm{Hf}\oplus(\text{constants})$: the search returns the generator.*

The mechanism behind (3) is the cubic distribution itself, in the form
$$\boxed{\ \mathscr H(3u)+3\mathscr H(u)-3\mathfrak h_{1,\pi/3}(u)=-2G\ }$$
(equivalent, via $\mathscr H(v+\pi)=2G-\mathscr H(v)$, to the note's
$\mathscr H(3u)+3\sum_{j}\mathscr H(u+\tfrac{2\pi j}3)=10G$ — **the note's $10G$ is confirmed**),
applied at levels $3$ and $9$.  The $\tfrac{10}9=1-\chi_{-4}(3)/9$ is exactly $\mathscr C(\pi/3)=\tfrac23G$
of (1.3).

**Calibration required by the brief, both directions.**
* the *full cubic trace* $\mathscr H(3u)+3\mathscr H(u)-3\mathfrak h_{1,\pi/3}$ has $\Psi\equiv0$
  and all Taylor coefficients $0$ — **it collapses to a constant** ✅;
* $K=\mathscr H(3u)+3\mathscr H(u)$ has singular set $\{1\pm\tfrac{\sqrt3}2\}$, value $4G$ at
  $z=1$, and $v_2([z^n]K)=n+1$ *exactly* for all $n\le1200$ — **$\varsigma_2=1$, the note's
  §13 claim, reproduced on the nose** ✅ (and the note's "$O(\log n)$" is in fact $O(1)$).

---

## 5. Arithmetic of the span (measured to $n\le1300$)

> **Theorem 3 (denominator type of the generator).**  $\ n^2\binom{2n}{n}h_n\in\mathbf Z$ for
> all $n\le1300$ (`quick_type.py`), with
> $$c_n:=n^2\tbinom{2n}nh_n=1,4,14,48,166,584,2092,7616,28102,104824,\dots$$
> with the closed form
> $$c_n=\frac{n\binom{2n}n}{2^n}\sum_{k=0}^{n-1}\frac{\binom{n-1}k}{2k+1}
>       =\frac{n\binom{2n}n}{2^n}\int_0^1(1+t^2)^{n-1}dt,$$
> and $h_n\sim\tfrac1{2n^2}$ exactly ($n^2h_n\to\tfrac12$ to 4 digits at $n=3000$).  Hence
> $$\tfrac1n\log\operatorname{den}(h_n)\longrightarrow 2\log2=1.38629\quad(\text{measured }1.3872\text{ at }n=2000),$$
> and the prime law is $p\mid\operatorname{den}(h_n)\iff\{n/p\}\ge\tfrac12$, whose density
> constant is $\int_0^2\chi(\{1/u\}\ge\tfrac12)\,du=1+\tfrac16+\bigl(2\log2-1-\tfrac16\bigr)=2\log2$.

This is the *arcsin$^2$ type*: $\bigl(2\arcsin\sqrt{z/2}\bigr)^2=2\sum_n(2z)^n/(n^2\binom{2n}n)$
carries the same denominator, as it must — $\binom{2n}n$ is the price of the quadratic
coordinate $z=1-\cos(\mathrm{ang})$.  In CDT's shape language the cheapest legal description is
$n^2[1..2n]$ ($\sigma=2$); the *true* size is $1.386$.  Note $\operatorname{den}(h_n)\nmid[1..2n]$
(e.g. $5^6\,\|\,\operatorname{den}(h_{1200})$ while $v_5[1..2400]=4$), so the $n^2$ layer is
not optional.

**Uniformity.**  Every one of the $24$ basis elements $g_{a,0},g_{a,\pi/3}$, $a\le12$, has
$\tfrac1n\log\operatorname{den}=1.386\pm0.006$ at $n=1000$ (`emn_slopes.py`).  *The
denominators do not depend on the projector at all.*

**$p$-adic slope census** (`slopes.json`, $n=1000$; slope $=v_p/n$):

| function | $v_2$ slope | $v_3$ | $v_5$ | $v_7$ |
|---|---|---|---|---|
| $\mathscr H(a\,\mathrm{ang})$, $a$ odd | $0$ | $0$ | $0$ | $0$ |
| $\mathscr H(a\,\mathrm{ang})$, $a=2,6,10$ | $+\tfrac12$ | $0$ | $0$ | $0$ |
| $\mathscr H(a\,\mathrm{ang})$, $a=4,12$ | $+\tfrac34$ | $0$ | $0$ | $0$ |
| $\mathscr H(8\,\mathrm{ang})$ | $+\tfrac78$ | $0$ | $0$ | $0$ |
| $\mathfrak h_{a,\pi/3}$, **all** $a\le12$ | $+1$ | $0$ | $0$ | $0$ |
| $K=\mathscr H(3u)+3\mathscr H(u)$ | $+1$ | $0$ | $0$ | $0$ |

The empirical law is $\varsigma_2\bigl(\mathscr H(a\,\mathrm{ang})\bigr)=1-2^{-v_2(a)}$, and it is
*exactly the rescaling content*: $\mathscr H(a\,\mathrm{ang})=\mathrm{Hf}(R_a(z))$ with
$R_a(z)=1-T_a(1-z)$, and $R_2(z)=2z(2-z)$ carries the $2$.  **There is no $3$-adic slope
anywhere in the span** — in particular the note's §7 $3$-adic realisation
$\Phi_3(3/2)=L_3(2,\chi_{12})$, whatever its status, buys no $\gamma_3$.  And the $2$-adic
slopes cannot be harvested: `ADELIC_HOLONOMY.md` §2.6 requires $(\ast_p)$ for **all** $m$
functions, the pure module has slope $0$, and rescaling $z\mapsto z/2$ converts $\log R_2$
into exactly the same $\log2$ of conformal size.  $\sum_p\gamma_p=0$.

---

## 6. The fold that does exist: $z=2$

This is the one genuinely new object in this note, and it repairs the F3 failure of
`SOL_NOTES_DIGEST.md` §2.4 ("no conditional function at all").

**The geometry.**  The local system of $\mathrm{Hf}$ is $\{1,\ \mathrm{ang}(z)=2\arcsin\sqrt{z/2},\ \mathrm{Hf}\}$
with $\gamma_1:\mathrm{Hf}\mapsto\mathrm{Hf}-\pi i(\mathrm{ang}-\tfrac\pi2)$, $\gamma_0,\gamma_2:\mathrm{ang}\mapsto-\mathrm{ang},\,2\pi-\mathrm{ang}$
(both fix $\mathrm{Hf}$).  Since $\gamma_1\gamma_0\gamma_1^{-1}\gamma_0^{-1}(\mathrm{Hf})=\mathrm{Hf}-2\pi i\,\mathrm{ang}$,
the honest host is
$$\Sigma=\{0,1,2,\infty\},\qquad X=\mathbf P^1\setminus\Sigma .$$
$\mathrm{Hf}$ is **analytic** at $z=0$ and at $z=2$; those two points are punctures only because
the pure companion $\mathrm{ang}(z)$ branches there ($\gamma_0,\gamma_2$ have order $2$).

**Exact uniformisation of $X$** (this answers the brief's "compute the uniformisation").
$y=2z-z^2=\sin^2(\mathrm{ang})$ maps $X\to Y=\mathbf P^1\setminus\{0,1,\infty\}$; the fibres are
$y=0\mapsto\{0,2\}$, $y=1\mapsto\{1\}$, $y=\infty\mapsto\{\infty\}$, so on the *punctured*
surfaces this is an **unramified** degree-$2$ cover ($\chi=-2=2\cdot(-1)$), and the cusp $z=0$
has width $1$.  With the $\lambda$-template $y=16q+O(q^2)$ and $z=\tfrac y2+\tfrac{y^2}8+\cdots$,
$$\boxed{\ z(q)=8q+O(q^2)\ \Longrightarrow\ |\varphi'(0)|\le8\ \text{ on }X\ }$$
— exactly half of $16$, and (contrast `CATALAN_MU4.md` §4) **not** $16^{1/2}=4$: the
involution $z\mapsto2-z$ *moves* the cusp instead of fixing it, so no square root is taken.
This is the sense in which the $z$-line beats $\mu_4$ host A by $\log2$.

**The fold.**  $\mathrm{Hf}(2)=2G$ (1.3) and $\mathrm{Hf}$ is regular there.  Put
$$A(z)=-\log\bigl(1-\tfrac z2\bigr),\qquad
\boxed{\ B(z)=\int_0^z\frac{\mathrm{Hf}(t)}{2-t}\,dt=\int_0^{\mathrm{ang}}\mathscr H(s)\tan\tfrac s2\,ds\ }$$
(the second form is the substitution $\sin s/(1+\cos s)=\tan\frac s2$, and shows the
continuation past $z=1$ is along the real axis, $\mathrm{Hf}$ being continuous there).  Then
$\Delta_{z=2}B=-4\pi iG$, $\Delta_{z=2}A=-2\pi i$, so

> **Theorem 4 (the EMN fold).**  $B-2G\,A$ *is regular at $z=2$; under $G=p/q$ the series*
> $qB-2pA\in\mathbf Q[[z]]$ *is a conditional function in CDT's sense, living on
> $\mathbf P^1\setminus\{0,1,\infty\}$ rather than on $X$.*

**Verified numerically** (`emn_fold.py`, 30 digits): $B+2G\log(2-z)$ tends to
$-0.3606366\ldots$ as $z\to2^-$ with increments $O(2-z)$ — a smooth finite limit, i.e.
$B-2GA\to-0.360637-2G\log2=-1.630343$.

**The caveat, identical to `CATALAN_MU4.md` §3.**  Only the principal branch is killed:
$\Delta_2\Delta_1B\in\mathbf Q^\times\pi^3$, so full descent needs the *three-period*
hypothesis that $1,\pi^2,G$ are $\mathbf Q$-dependent; with $G\in\mathbf Q$ alone one is in
CDT's `stacky overconvergence` regime and must impose $\#\varphi^{-1}(2)=1$.

**Every projector has this fold**, since $F(2)\in\mathbf QG$ for all $F$ in the span (§3).  So
the fold is a property of the *motive*, not of the projector — which is why the search of §4
cannot improve it.

---

## 7. The honest CDT entry

Work in $Z=z/2$ (fold at $Z=1$, so the kernel $1/(1-T)$ is integral); the surviving
singularity is $t_2=\tfrac12$ and the width-law ceiling after the fold is deleted is
$\log(16|t_2|)=\log8=2.0794$.  Equivalently in $z$: ceiling $\log16=2.7726$ and the
conditional row carries an extra $2^{-n}$.  Both accountings agree throughout.

**Denominators of the conditional row** (`emn_meas.py`, `quick_Btype.py`, $n\le600$):
$$\tfrac1n\log\operatorname{den}(B_n)\to2.97\quad(2.905,\,2.994,\,2.948,\,2.967,\,2.971,\,2.992,\,2.976
\text{ at }n=100,200,400,600,800,1200,1600),$$
with **no** cancellation against $\operatorname{lcm}_{k<n}\operatorname{den}(\widetilde h_k)$
(measured $2.912,2.968,2.944,2.957$): the partial sums $\sum_{k<n}h_k2^k$ accumulate the whole
lcm.  The cheapest valid CDT shape is
$$\boxed{\ B\ \text{has type }\ n\,[1..n]\,[1..2n],\qquad k=\sigma_B=3\ }$$
(verified integral for **all** $n\le220$, `quick_shape.py`; $n[1..n]^2$ fails at $216$ of
$220$ indices, $n[1..2n]$ at $216$, $[1..n][1..2n]$ at $159$, $n[1..n][1..\tfrac32n]$ at $213$
— so $k=3$ is the exact shape weight, and it matches the measured true size $2.95$ to
within the fluctuation).
In $z$: $\sigma_B=3.64$.  Compare CDT's own conditional row $n[1..n]^2$, $k=2$.

**Entry ledger** (`emn_ledger.py`; pure module = the CDT five functions transported to
$\mathbf P^1\setminus\{0,\tfrac12,\infty\}$, i.e. $1,\log(1-2Z),\operatorname{Li}_2(2Z),\log^2(1-2Z),J(2Z)$):

| description | $c$ | $m$ | $\sigma_m$ | $\tau^\flat$ | $\tau^\sharp$ | $\tau$ | ceiling | **entry** |
|---|---|---|---|---|---|---|---|---|
| measured types | 4 | 9 | 3.000 | 2.5802 | 0.8745 | 3.4547 | $+2.0794$ | $-1.375$ |
| measured types | 9 | 14 | 3.000 | 2.8265 | 0.6950 | 3.5216 | $+2.0794$ | $-1.442$ |
| **relaxed pure** | 4 | 9 | 3.000 | 2.9259 | 0 | 2.9259 | $+2.0794$ | $\mathbf{-0.847}$ |
| **relaxed pure** | 7 | 12 | 3.000 | 2.9583 | 0 | 2.9583 | $+2.0794$ | $\mathbf{-0.879}$ |
| relaxed pure | 9 | 14 | 3.000 | 2.9694 | 0 | 2.9694 | $+2.0794$ | $-0.890$ |
| crude Liouville ($\tau=$ true size $2.97$) | — | any | — | — | — | 2.970 | $+2.0794$ | $-0.891$ |

**Best entry $=-0.85$.**  (Minimising $\tau$ over *all* inventories, including long polylog
towers, gives $-0.42$ at $m=3$ — meaningless, since a proof needs $m>\mathrm{BC}/\text{entry}\approx6$;
entry is monotone decreasing in $m$ on this host, as on every Catalan host.)

**Candidate table** (`emn_table.py`, $n=400$; entry $=\log(16|t_2|)-\sigma_B^{(z)}$, crude):

| projector | $t_2$ | #poles | $\sigma_F$ | $\sigma_{B_F}$ | $\varsigma_2(F)$ | ceiling | **entry** |
|---|---|---|---|---|---|---|---|
| $\mathrm{Hf}=\mathscr H(\mathrm{ang})$ | $1$ | 1 | 1.379 | 3.641 | $0$ | $+2.7726$ | $\mathbf{-0.868}$ |
| $\mathscr H(3u)-3\mathfrak h_{1,\pi/3}=-3\mathrm{Hf}$ | $1$ | 1 | 1.376 | 3.638 | $0$ | $+2.7726$ | $-0.865$ |
| $\mathscr H(2u)$ | $0.2929$ | 2 | 1.388 | 2.942 | $+\tfrac12$ | $+1.5446$ | $-1.397$ |
| $\mathscr H(2u)\pm2\mathscr H(u)$ | $0.2929$ | 3 | 1.37 | 3.642 | $0$ | $+1.5446$ | $-2.097$ |
| $K=\mathscr H(3u)+3\mathscr H(u)$ (note §13) | $0.1340$ | 2 | 1.386 | 3.646 | $+1$ | $+0.7625$ | $-2.884$ |
| $\mathfrak h_{1,\pi/3}$ | $0.1340$ | 2 | 1.388 | 3.649 | $+1$ | $+0.7625$ | $-2.886$ |
| $\mathscr H(4u)+2\mathscr H(2u)+4\mathscr H(u)$ | $0.0761$ | 7 | 1.388 | 3.634 | $0$ | $+0.1972$ | $-3.437$ |
| $\mathscr H(5u)$ | $0.0489$ | 5 | 1.382 | 3.632 | $0$ | $-0.2445$ | $-3.877$ |
| full cubic trace | — | 0 | — | — | — | — | **constant** |

**Bost–Charles numerators** (`emn_bc.py`, $\varphi=\lambda/2$ on $|q|<r$, $N=2048$; sanity
$\mathrm{BC}(\rho z)=\log\rho$ to 13 digits):

| $r$ | $\log|\varphi'(0)|$ | $\mathrm{BC}$ | $\mathrm{BC}-\log|\varphi'(0)|$ |
|---|---|---|---|
| 0.20 | $+0.4700$ | $+0.4700$ | $0$ (univalent) |
| 0.40 | $+1.1632$ | $+1.9080$ | $+0.745$ |
| 0.60 | $+1.5686$ | $+3.3247$ | $+1.756$ |
| 0.75 | $+1.7918$ | $+4.6556$ | $+2.864$ |

Since entry $<0$ the margin $m\cdot\text{entry}-\mathrm{BC}$ is bookkeeping only
($\approx-13$ at $m=9$, $r=0.75$).

**Does the width law bind?**  Yes, unchanged by the trigonometric (non-modular) structure —
with one new term.  The fold cusp has width $w=1$ (the descent $z\mapsto2z-z^2$ is unramified
at $0$ and $2$), $t_2=1$ in the coordinate in which the primary row is integral, and $k=3.64$:
$$\text{entry}\ \le\ \tfrac1w\log16+\log|t_2|-k\ =\ 2.7726+0-3.641\ =\ -0.868 .$$
The novelty is *where* the $k$ comes from.  On every previous host the fold was the **nearest**
singularity ($1/9$ for CDT, $\tfrac18$ for Zagier $\mathbf E$, $\pm1$ for $\mu_4$), so its
position was a *credit*.  Here it is the **farthest**, and the law needs the extra term
$$k=k_0+\log\bigl|t_{\text{fold}}/t_2\bigr|,\qquad k_0=2.95,\ \log|t_{\text{fold}}/t_2|=\log2 .$$
Counterfactuals (`emn_ledger.py`): with $k=2$ (a CDT-quality conditional row) entry
$=-0.49$; with the fold moved onto $t_2$ (so $\log|t_{\rm fold}/t_2|=0$) entry $=-0.70$.
(The note's own $a_n=\bigl(2^{n+1}(n+1)\bigr)^{-1}\sum_k\binom nk/(2k+1)$ was checked against the
ODE recursion for $n<60$: identical.)
**Neither repair alone, nor the two together at the crude $\tau$ ($\log8-2=+0.08$), survives
the rearranged $\tau$.**

---

## 8. Where this sits, and what binds

$$\text{CDT }L(2,\chi_{-3})\ +0.77\ \gg\ \text{level 8 sym }-0.077\ >\ \mu_4\text{ host A }-0.563\ >\ \boxed{\text{EMN }z\text{-line }-0.85}\ >\ \text{level 16 }-1.46 .$$

The comparison with $\mu_4$ is exact and instructive, because by (1.2) the two are the same
module on the two sides of the degree-$2$ quotient $w\mapsto1/w$:

| | ceiling | conditional $k$ | entry |
|---|---|---|---|
| $\mu_4$ host A ($w$-line, $\operatorname{Ti}_2$ lcm-free) | $\log4=1.386$ | $2$ | $-0.61$ crude, $-0.563$ |
| EMN $z$-line (quotient) | $\log16=2.773$ | $2.95+\log2=3.64$ | $-0.87$ crude, $-0.85$ |

Descending to the quotient buys $\log4=1.386$ of conformal size and costs
$0.95$ (the $\binom{2n}n$ of the quadratic coordinate leaking through the integration)
$+\ \log2=0.693$ (the fold doubles its distance).  **Net $-0.26$: the quotient is a net
loss, for the same reason `CATALAN_MU4.md` §5 gives — the involution that gains the size
also moves the fold.**

**What binds, in order.**
1. **The conditional row's weight $k=3$** (needs $<2.08$).  It is $2.95$ because the fold
   companion is a *partial sum* $\sum_{k<n}h_k2^k$ and the $h_k$ have coprime
   $\binom{2k}k$-type denominators, so the lcm is taken with no cancellation.  This is
   intrinsic: $B'=\mathrm{Hf}/(2-z)$ forces it.
2. **The fold's position**, $\log2$, novel to this host.
3. Not binding, and this is the news: the **conformal size is fine** ($t_2=1$, $|\varphi'(0)|=16$
   after the fold — the first Catalan host in the ledger with $|t_2|=O(1)$ *and* a genuine
   fold), the **denominators of the primary row are fine** ($1.386$, better than $[1..n]^2$),
   and $\tau^\sharp$ is avoidable.
4. $\sum_p\gamma_p=0$: no $p$-adic rescue (§5), exactly as `CATALAN_OBSTRUCTION.md` (2.3)
   predicts, and here for the sharper reason that the slopes present are pure rescaling.

**What would evade.**  A conditional function on this host of lcm weight $<2.08$.  Equivalently:
an EMN Apéry *pair* $(A_n,B_n)$ with $B_n/A_n\to2G$ realised by a recurrence rather than by
$\int\mathrm{Hf}\,dt/(2-t)$ — i.e. `CATALAN_OBSTRUCTION.md` (E3) with the extra requirement that
the fold companion be a $q$-expansion-type object, not an antiderivative.  The note's §17.1
("compute the integral structure of the EMN motive itself") is exactly the right instinct, and
this document makes it quantitative: **one needs to lose $0.9$ nats on $B$, not on $\mathrm{Hf}$.**

**Ledger amendment.**  `SOL_NOTES_DIGEST.md` §2.4's entry for the moving EMN period should
read: *fold* **yes** (at $z=2$, $\mathrm{Hf}(2)=2G$, via $\int\mathrm{Hf}\,dt/(2-t)$), *ceiling*
$\log 8$ before / $\log16$ after the fold, $\tau=2.95$–$3.0$ for the conditional row,
**entry $-0.85$** — not $+1.27$ counterfactual and not F3.  And
`CATALAN_OBSTRUCTION.md`'s standing question ("a host where $G$ is an Apéry limit at a fold
*and* $|t_2|=O(1)$; the two have never yet co-occurred") is **answered in the affirmative
here** — and the answer is still not enough, because the co-occurrence costs $k$.

---

## 9. Not done

* $\mathbf Q(z)$-independence of the transported pure module and of $\{B,\int B,\dots\}$ is
  assumed, not checked (it is checked for the same functions in `CATALAN_MU4.md` §2).
* No slit/gobble contour was designed; only concentric $|q|<r$ (§7 BC table).
* Theorem 2 is a finite verification for $A\le24$, not a proof for all $N$; the parity
  argument of §3 explains why mirror pairs cannot be split, but not why mixing parities
  never wins.
* $c_n=1,4,14,48,166,584,2092,7616,28102,104824,394404,1494240$ has the closed form above;
  a hypergeometric/combinatorial identification (and a direct proof of its integrality, which
  here is only verified to $n=1300$) is open.  $c_n-4c_{n-1}=0,-2,-8,-26,-80,-244,-752,\dots$
* The note's §7 claim $\Phi_3(3/2)=L_3(2,\chi_{12})$ was not re-verified here; §5 shows it is
  in any case slope-free.
