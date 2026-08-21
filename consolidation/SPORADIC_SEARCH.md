# A mechanical search for sporadic Apéry-like sequences from eta-quotient data on genus-zero $X_0(N)$

*Claude (Fable), 2026-08-21. Scripts: `lattice/sporadic_search/{search.py,run_scan2.py,analyze.py,limits.py,make_report.py}`.
Raw output: `lattice/sporadic_search/{hits.json,table.json,table.md,scan_all.log,lindep.out}`.*

This executes item 3 of `ZETA3_TWO_LATTICE.md` §7.4 ("systematically enumerating
$(N,\text{cusp},w)$ and testing integrality is a well-defined search for new sporadic
sequences"), with the §5.1 design rule of §7.4(1) as the ranking criterion.

**Headline.** The construction was run exhaustively over all $14$ genus-zero levels;
**520 005** pairs $(t,F)$ were tested. It reproduces **all six** third-order
Almkvist–Zudilin rows and **five of the six** Zagier second-order rows, together with
their Apéry limits ($\tfrac16\zeta(3)$, $\tfrac7{24}\zeta(3)$, $\tfrac7{32}\zeta(3)$,
$\tfrac12 L(2,\chi_{-3})$, $\tfrac13 L(3,\chi_{-3})$, $\tfrac14\zeta(2)$, $\tfrac12 G$),
which is the acceptance test.  Beyond them it produces **no new sporadic row in Zagier's
normalisation**: every additional hit is either a two-term (hypergeometric) degeneration,
a lacunary/quadratic degeneration, a classical convolution sequence, or a Möbius
reparametrisation of a known row.  The single object that is not obviously accounted for
is one **weight-one row over the Domb parameter** with $\Lambda=16$ and budget $+0.773$
(§5).  **The search does not find any row with $\Lambda$ larger than Apéry's own
$\Lambda=12+8\sqrt2+\dots=33.97$; §7.4(2)'s hope of buying a larger budget at higher
level is not realised inside this family.**

---

## 1. The construction, exactly as enumerated

For each genus-zero $N\in\{1,\dots,10,12,13,16,18,25\}$ and each divisor list $D=\{d:d\mid N\}$:

**Parameter.** $t=\prod_{d\mid N}\eta(d\tau)^{r_d}$ with
$$\sum_d r_d=0,\qquad \tfrac1{24}\sum_d d\,r_d=1\quad(\Rightarrow t=q+O(q^2)),$$
subject to the Ligozat conditions $\sum_d d\,r_d\equiv\sum_d (N/d)r_d\equiv0\pmod{24}$
and integrality of every cusp order
$$\operatorname{ord}_c(f)=\frac{N}{24\gcd(c^2,N)}\sum_{d\mid N}\frac{\gcd(c,d)^2}{d}r_d ,\qquad c\mid N .$$
Box $|r_d|\le24$ ($\le16$ when $N$ has $\ge5$ divisors).

*Note.* An earlier pass additionally demanded that $t$ be a **Hauptmodul** (single simple
pole).  That is too strict and it **loses Domb and Apéry**: the Domb parameter
$t=(\eta_2\eta_6/\eta_1\eta_3)^6$ has *two* simple poles on $\Gamma_0(6)$ (it is a
Hauptmodul for the Fricke quotient, not for $\Gamma_0(6)$).  The final scan therefore
imposes only $\operatorname{ord}_\infty t=1$.  This is worth recording: **the sporadic
rows do not all come from Hauptmoduls of $\Gamma_0(N)$ itself.**

**Form.** $F=\prod_d\eta(d\tau)^{s_d}$ of weight $w\in\{1,2\}$:
$$\sum_d s_d=2w,\qquad \sum_d d\,s_d=0\ (\Rightarrow F=1+O(q)),\qquad \operatorname{ord}_c F\ge0\ \forall c .$$
Nebentypus is allowed (the extra Ligozat "$\prod d^{r_d}$ is a square" condition is *not*
imposed — imposing it kills level 5 and hence $(11,5,125)$).

**Sequence.** $a_n$ defined by $F(q)=\sum_n a_n t(q)^n$, obtained by triangular peeling
against the powers of $t$.  Since $t=q+O(q^2)$ with integer coefficients, $a_n\in\mathbf Z$
**automatically**; integrality is therefore *not* a discriminating test in this family
(all 520 005 candidates were integral).  The discriminating test is the recurrence.

**Test.** Fit $p_2(n)a_{n+1}+p_1(n)a_n+p_0(n)a_{n-1}=0$ with $\deg p_i\le3$ (12 unknowns)
by exact linear algebra; pre-filter by rank over $\mathbf F_p$, $p=2^{61}-1$, on 20
equations at scan precision 28; then **re-derive the sequence to $q^{66}$ and verify the
recurrence on every $n\le64$**.  Survivors are then matched against
$(n+1)^2a_{n+1}=(an^2+an+b)a_n-cn^2a_{n-1}$ (order 2) and
$(n+1)^3a_{n+1}=(2n+1)(an^2+an+b)a_n-cn^3a_{n-1}$ (order 3).
Both signs of $t$ are tried ($a_n\mapsto(-1)^na_n$).

## 2. Scan size

| $N$ | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 12 | 13 | 16 | 18 | 25 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| $t$ candidates | 1 | 1 | 3 | 1 | 37 | 1 | 21 | 5 | 120 | 1647 | 1 | 127 | 4164 | 9 |
| $F$, $w=1$ | 0 | 0 | 0 | 0 | 3 | 0 | 3 | 1 | 0 | 15 | 0 | 7 | 16 | 0 |
| $F$, $w=2$ | 0 | 0 | 2 | 1 | 6 | 0 | 6 | 2 | 10 | 60 | 0 | 19 | 78 | 2 |

$N=1$ admits no eta quotient with $\operatorname{ord}_\infty t=1$.
Total pairs scanned **520 005**; 99 were rank-deficient at scan precision; after
de-duplication and exact verification to $n=64$, **50 distinct sequences** survive.
Levels $10,13,18$ produced **no** surviving sequence at all.

## 3. Sanity check: the known lists are recovered

**Almkvist–Zudilin / third order — 6 of 6.**

| $(a,b,c)$ | $N$ | $t$ | $F$ | $a_n$ | Apéry limit (recomputed, `lindep`) |
|---|---|---|---|---|---|
| $(17,5,1)$ | 6 | $(\eta_1\eta_6/\eta_2\eta_3)^{12}$ | $\eta_2^7\eta_3^7/\eta_1^5\eta_6^5$ | 1,5,73,1445,33001 | $\tfrac16\zeta(3)$ |
| $(12,4,16)$ | 8 | $\eta_1^{-8}\eta_2^{16}\eta_4^{-16}\eta_8^{8}$ | $\eta_1^4\eta_2^{-6}\eta_4^{10}\eta_8^{-4}$ | 1,4,40,544,8536 | $\tfrac7{32}\zeta(3)$ |
| $(10,4,64)$ (Domb) | 6 | $(\eta_2\eta_6/\eta_1\eta_3)^{6}$ | $(\eta_1\eta_3)^4(\eta_2\eta_6)^{-2}$ | 1,4,28,256,2716 | $\tfrac7{24}\zeta(3)$ |
| $(9,3,-27)$ | 9 | $\eta_1^6\eta_3^{-12}\eta_9^6$ | $\eta_1^{-3}\eta_3^{10}\eta_9^{-3}$ | 1,3,27,309,4059 | $\tfrac13 L(3,\chi_{-3})$ |
| $(11,5,125)$ | 5 | $(\eta_5/\eta_1)^6$ | $\eta_1^5\eta_5^{-1}$ | 1,5,35,275,2275 | $1.129718\ldots$ (not in basis) |
| $(7,3,81)$ | 6 | $\eta_1^{-4}\eta_2^{-4}\eta_3^4\eta_6^4$ | $\eta_1^3\eta_2^3\eta_3^{-1}\eta_6^{-1}$ | 1,3,9,3,−279 | $0.327349\ldots$ (not in basis) |

**Zagier / second order — 5 of 6.**

| $(a,b,c)$ | $N$ | $a_n$ | limit |
|---|---|---|---|
| $(17,6,72)$ | 6, 12 | 1,6,42,312,2394 | $0.488314\ldots$ |
| $(10,3,9)$ | 6, 12 | 1,3,15,93,639 | $\tfrac12 L(2,\chi_{-3})$ |
| $(12,4,32)$ | 8 | 1,4,20,112,676 | $\tfrac12 G$ (Catalan) |
| $(7,2,-8)$ | 6 | 1,2,10,56,346 | $\tfrac14\zeta(2)$ |
| $(9,3,27)$ | 9 | 1,3,9,21,9 | $0.634021\ldots$ |
| $(11,3,-1)$ | — | **not produced** | — |

**Why $(11,3,-1)$ is missed, and why that is a property of the construction, not a bug.**
Its parameter is level $5$ and its $F$ must have weight $1$.  On $\Gamma_0(5)$ a
weight-one eta quotient would need $s_1+s_5=2$, $s_1+5s_5=0$, i.e. $s_5=-\tfrac12$:
no integral solution.  On $\Gamma_0(10)$ the exhaustive search ($|s_d|\le24$) finds no
weight-one eta quotient holomorphic at every cusp with $\operatorname{ord}_\infty=0$.
The relevant weight-one form (the $\chi_5$ Eisenstein series) is simply **not an eta
quotient**.  So the honest statement of scope is:

> *The eta-quotient construction on genus-zero $\Gamma_0(N)$ recovers the whole
> Almkvist–Zudilin third-order list and exactly those Zagier rows whose weight-one form
> is an eta quotient; $(11,3,-1)$ falls outside because its form is not.*

The recomputation of the Apéry limits is an independent check on the whole chain: the
project's ledger values $\tfrac7{24}\zeta(3)$ (Domb) and $\tfrac7{32}\zeta(3)$ (T-row)
came out of `lindep` on the companion recurrence with no input from the ledger.

## 4. Everything else the scan found, and why it is not new

The full table is `lattice/sporadic_search/table.md` (50 rows, sorted by budget).
The non-known rows fall into four groups.

1. **Two-term degenerations ($p_0\equiv0$).** $(32,8,0)$ at $N=4$, which is
   $a_n=\binom{2n}n^3$; $(16,4,0)$ at $N=8$, which is $a_n=\binom{2n}n^2$; and two
   $\deg p_2=2$ rows at $N=8,9$ with $p_1=4(4n+1)^2$, $3(3n+1)^2$.  These have
   $\Lambda=64,\,27$ and hence formally huge budget ($+2.16$, $+1.30$), but a two-term
   recurrence has **no second solution**, so there is no companion, no extension class and
   no resource to spend.  They are the classical hypergeometric families, not sporadic.
2. **Lacunary degenerations ($p_1\equiv0$).** $(0,0,\pm16)$, $(0,0,64)$: $a_{2k+1}=0$;
   these are pullbacks along $t\mapsto t^2$ of the rows above.  Apéry limit $0$.
3. **A classical convolution.** $(16,8,256)$ at $N=4$, $a_n=1,8,88,1088,14296,\dots$
   Verified identity
   $$a_n=\sum_k\binom{2k}k^2\binom{2n-2k}{n-k}^2 ,$$
   the self-convolution of $\binom{2n}n^2$.  Characteristic polynomial $x^2-32x+256$ has a
   **double** root $16$, so $a_n\sim16^n\cdot\text{poly}$ and the budget is negative.
   Well known; not sporadic.
4. **Non-Zagier-normal rows (29 of them).** These have a genuine three-term recurrence of
   degree $\le3$ but with $p_2=(n+1)^2(n+2)$ or similar, i.e. they are **Möbius
   reparametrisations $t\mapsto t/(1+\alpha t)$** of the known rows (they arise because
   the scan admits every $t$ with $\operatorname{ord}_\infty t=1$, not only the extremal
   ones).  Every one of them has a **rational** Apéry limit — the diagnostic that the
   extension class has been trivialised by the coordinate change — except the four
   discussed in §5.

## 5. The one row that is not accounted for

$N=12$ (equivalently the level-6 Domb parameter):
$$t=-\Bigl(\frac{\eta_2\eta_6}{\eta_1\eta_3}\Bigr)^{6}\ \ (\text{the Domb parameter}),\qquad
F=\frac{\eta_1^2\eta_3^2}{\eta_2\eta_6}\quad(\text{weight }1).$$

Note $F^2=\eta_1^4\eta_3^4/\eta_2^2\eta_6^2$ **is exactly the Domb weight-two form** of the
table in §3, over the same $t$.  So this row is literally the square root of the Domb row's
modular data — the $\operatorname{Sym}^1$ to Domb's $\operatorname{Sym}^2$.
(The same sequence reappears at $N=12$ from
$t=(\eta_1\eta_3\eta_4\eta_{12}/\eta_2\eta_6)^{6}$, $F=\eta_2^5\eta_6^5/\eta_1^2\eta_3^2\eta_4^2\eta_{12}^2$,
without the sign flip.)
$$a_n=1,\;2,\;12,\;104,\;1078,\;12348,\;150528,\;1914432,\;25108668,\;337111736,\dots$$
$$\boxed{(n+1)^2a_{n+1}=(20n^2+10n+2)\,a_n-16(2n-1)^2\,a_{n-1}}$$
**verified directly from the $q$-expansion for every $n\le208$.**
Characteristic polynomial $x^2-20x+64$: roots $\Lambda=16$, $\lambda_2=4$ — *identical to
Domb's*.  Slope $\sigma_2=v_2(64)=6$, again as Domb.  Budget
$$\log\Lambda-k=\log16-2=0.7726>0.53 .$$
The companion $b_0=0,b_1=1$ satisfies $d_n^2b_n\in\mathbf Z$ ($d_n=\operatorname{lcm}(1..n)$,
checked to $n=40$), so it is a *bona fide* Apéry construction, and
$$\lim_{n\to\infty}\frac{b_n}{a_n}=0.737292996185596240176426197802293697960919143276529185222868885113\ldots$$
`lindep` at 60 digits against $[1,\zeta(2),\zeta(3),G,L(2,\chi_{-3}),L(3,\chi_{-3}),\pi^3,\pi^4]$,
and at 55 digits against a 15-element enlarged basis (adding $\log2$, $\pi\log2$,
$\log^22$, $\zeta(3)/\pi^2$, $\pi^2\log2$, $\sqrt2$, $\sqrt3$), **finds nothing**.
No simple binomial-sum closed form was found either (a brute-force scan over products of
$\binom nk,\binom{n+k}k,\binom{2k}k,\binom{2n-2k}{n-k},\binom{2k}n,\binom{3k}k$ with
exponents $\le2$ and total degree $\le4$ returned no match).

**Honest assessment.**  This row is *not* a counterexample to Zagier's classification and
should not be advertised as "a seventh sporadic sequence": Zagier searched the **symmetric**
normalisation $p_1=-(an^2+an+b)$, $p_0=cn^2$, and this recurrence is not of that shape
($p_1=-(20n^2+10n+2)$ is not symmetric under $n\mapsto-n-1$, and $p_0=16(2n-1)^2\ne cn^2$).
It is therefore *outside the class he enumerated*, exactly as it should be.  What it is, in
the §7 language, is the **$\operatorname{Sym}^1$ (weight-one) companion of the Domb
$\operatorname{Sym}^2$ local system on the same parameter** — same $t$, same characteristic
polynomial, same $2$-adic slope, different weight.  Its interest is that its extension class
appears not to be a $\mathbf Q$-multiple of anything in the standard period basis, and that
it carries budget $0.773$ at $k=2$ — the largest in the whole scan among rows with a
genuine companion.  Two further rows are in the same situation (unidentified, non-rational
limit): $N=16$, $w=1$, $\Lambda=16$ (double root $-16$), limit $0.865086\ldots$; and $N=9$, $w=2$,
$\Lambda=27$ (double root $-27$), limit $-0.424611\ldots$.
Both have a repeated characteristic root, so their $a_n$ grow like $\Lambda^n$ times a
polynomial and the two-lattice machinery has nothing to separate.

## 6. Ranking by the §5.1 budget, and the negative conclusion

Top of `table.md` restricted to rows with a genuine second solution:

| rank | row | order $k$ | $\Lambda$ | budget $\log\Lambda-k$ | limit |
|---|---|---|---|---|---|
| 1 | §5 row, $N=12$, $w=1$ | 2 | 16 | **+0.7726** | unidentified |
| 2 | $N=16$, $w=1$, double root | 2 | 16 | +0.7726 | unidentified |
| 3 | **$(17,5,1)$ Apéry**, $N=6$ | 3 | $33.9706$ | **+0.5255** | $\tfrac16\zeta(3)$ |
| 4 | $N=9$, $w=2$, double root | 3 | 27 | +0.2958 | unidentified |
| 5 | $(17,6,72)$, $(10,3,9)$ | 2 | 9 | +0.1972 | $0.4883\ldots$, $\tfrac12L(2,\chi_{-3})$ |
| 6 | $(12,4,16)$ (the T row) | 3 | $12+8\sqrt2$ | +0.1490 | $\tfrac7{32}\zeta(3)$ |
| 7 | $(7,2,-8)$, $(12,4,32)$ | 2 | 8 | +0.0794 | $\tfrac14\zeta(2)$, $\tfrac12 G$ |
| 8 | $(10,4,64)$ Domb | 3 | 16 | −0.2274 | $\tfrac7{24}\zeta(3)$ |

The $\zeta(3)$ row of `ZETA3_TWO_LATTICE.md` is $(12,4,16)$ at $+0.1490$; the Apéry row
itself sits at $+0.5255$; the threshold quoted in the brief ($0.53$, Apéry's $\zeta(3)$ row)
is essentially the Apéry row's own value.  **Nothing in the scan clears it except the two
$\Lambda=16$, $k=2$ rows of §5, and those have unidentified periods.**  In particular:

> **Negative result (the useful one).** Over all $14$ genus-zero levels and all $520\,005$
> eta-quotient pairs, the largest characteristic root attached to a *known-period*
> Apéry-like row is $\Lambda=12+8\sqrt2+\dots=33.97$, i.e. Apéry's own $(17,5,1)$.
> Raising the level does **not** raise $\log\Lambda$: levels $10$, $13$, $18$ contribute
> nothing at all, and levels $12$, $16$ only reproduce rows already present at $4,6,8,9$.
> §7.4(2)'s suggestion — "look for a cusp of large width at high level" — is, within the
> eta-quotient family, **false**: the cusp widths grow but the $q$-expansions that survive
> the recurrence test do not.  The ceiling $\gamma\log\Lambda_{\rm dec}$ of §7.5 cannot be
> raised this way.

## 7. What was not done

- **Fricke quotients $\Gamma_0(N)+$** were not enumerated separately.  Partially moot: the
  Domb and Apéry parameters, which *are* Fricke Hauptmoduls, are already in the scan
  because the only condition imposed on $t$ is $\operatorname{ord}_\infty t=1$.
- **Weights $w\ge3$** (fourth-order recurrences, Cooper's cubic systems) were not scanned.
- **Non-eta-quotient forms.**  This is the real gap, and $(11,3,-1)$ is the proof that it
  matters: an entire Zagier row is invisible to the eta-quotient parametrisation.  A scan
  over bases of $M_w(\Gamma_0(N),\chi)$ rather than eta quotients would be strictly larger
  and is the natural next step.
- The period of the §5 row is unidentified; it is the one loose end worth pulling.

## 8. Reproduction

```
python3 lattice/sporadic_search/run_scan2.py            # full scan -> hits.json (~4 min)
python3 lattice/sporadic_search/analyze.py              # roots, budget -> table.json
python3 lattice/sporadic_search/limits.py               # Apery limits -> _lindep.gp
gp -q lattice/sporadic_search/_lindep.gp > lindep.out   # identification
python3 lattice/sporadic_search/make_report.py          # -> table.md
```

---

## ADDENDUM (parent agent): the unidentified limit is a cusp-form $L$-value

The row reported above as "unaccounted for" —
$$(n+1)^2a_{n+1}=(20n^2+10n+2)a_n-16(2n-1)^2a_{n-1},\qquad a_n=1,2,12,104,1078,12348,\dots$$
arising over the Domb parameter with the **weight-1** form $F=\eta_1^2\eta_3^2/(\eta_2\eta_6)$
(so $F^2=F_{\text{Domb}}$) — has now been identified.

$$\boxed{\ \lim_{n\to\infty}\frac{b_n}{a_n}=L(f,2),\qquad
f=\eta(2\tau)^3\eta(6\tau)^3\in S_3\bigl(\Gamma_0(12),\chi_{-3}\bigr)\ }$$

the weight-3 CM newform of level 12, $f=q-3q^3+2q^7+9q^9-22q^{13}+\cdots$.

*Verified*: agreement to **60 digits** ( `lindep` relation exactly $[-1,0,1]$, i.e. the limit is
$L(f,2)$ **on the nose**, with no rational factor); the eigenform was produced independently
by `mfinit([12,3,-3],0)` and confirmed equal to `mffrometaquo([2,3;6,3])`.
$L(f,2)=0.737292996185596240176426197802293697960919143276529185222869\ldots$

### Why this matters: it has the largest budget in the census, and the reason is $k=2$

Characteristic roots of $x^2-20x+64$ are $16$ and $4$, so $\Lambda=16$, $\lambda_2=4$,
$c=\Lambda\lambda_2=64$, $\sigma_2=v_2(64)=6$. Denominators: $d_n^2b_n\in\mathbf Z$
[verified $n\le200$] — **$k=2$, not 3**.

$$\text{score}=\log\tfrac14-2=-3.3863,\qquad
\text{budget}=\log16-2=+0.7726,\qquad \text{headroom}=6\log2=4.1589 .$$

Compare `certificates/documents/main(6).tex`, which manufactures a cusp-form apparatus on
the *same* Domb curve but with the **weight-4** level-6 newform $f_6=(\eta_1\eta_2\eta_3\eta_6)^2$,
limit $L(f_6,3)/2$, $d_n^3B^*_n\in\mathbf Z$, rate $4^{-n}$. Its abstract does the bookkeeping
honestly: "$4^n$ against $d_n^3\sim e^{3n}$ loses, so no irrationality follows" — budget
$\log16-3=-0.2274$, **negative**.

Our row is the **weight-3** analogue on the same curve with the same $4^{-n}$ rate, and it
costs only $d_n^2$. That single unit of denominator flips the sign:

| | source weight | limit | $k$ | budget |
|---|---|---|---|---|
| `main(6).tex` row | 4 ($f_6$, level 6) | $L(f_6,3)/2$ | 3 | $-0.2274$ |
| **this row** | 3 ($f$, level 12, $\chi_{-3}$) | $L(f,2)$ | **2** | $\mathbf{+0.7726}$ |

So "the cusp form costs nothing" (main(6).tex §§195, 413) is *more* true at weight 3, and
this is the largest budget of any row in the whole census — larger than Apéry's $+0.5255$.

### The honest caveat, stated plainly

**Budget is a ceiling, not an achievement.** This row has $\lambda_2=4>1$: its linear form
**grows**, and its own score is $-3.386$. Alone it proves nothing whatsoever. It is an
*engine* in the §5 sense, and a good one — $\sigma_{\rm eng}=6$, $k=2$, $\log\rho_2=\log4$,
strictly better than Catalan's E row ($\sigma=5$) — but an engine is useless without a
decayer **sharing its period**. The period here is $L(f,2)$ for a weight-3 CM cusp form,
and no Apéry-like row with that period and $|\lambda_2|<1$ is known to exist.

By the §5.1 design rule with this engine ($k=2$, $\log\rho_2^{\rm eng}=\log4$, $\sigma_{\rm eng}=6$):

| $\sigma_{\rm dec}$ | $r=\sigma_{\rm dec}/6$ | need $\Lambda_{\rm dec}>$ |
|---|---|---|
| 3 | $1/2$ | $14.8$ |
| 4 | $2/3$ | $18.6$ |
| 6 | $1$ | $29.6$ |

**New search target, sharper than anything previously on the list: an Apéry-like row with
period $L(f,2)$, $f=\eta_2^3\eta_6^3$, decaying linear form, and $\Lambda>14.8$.**
This is a well-posed question that did not exist before this session. It is also entirely
possible no such row exists — the census's other clean negative (raising the level does not
raise $\log\Lambda$ inside the eta-quotient family) is evidence against.

**Correction (Fable, review pass):** the row's trailing coefficient is $16(2n-1)^2$, not $64n^2$, so its Casoratian is $\propto 4^n\binom{2n}{n}^2$ and the 2-adic slope is **σ₂ = 2**, not 6 (measured 2 by the census agent; confirmed by hand). The budget $\log16-2=+0.77$ is unaffected, but the engine parameters and the partner thresholds quoted above must use σ_eng = 2 (Λ_dec > 14.8, 29.6, 160 for σ_dec = 1, 2, 3).
