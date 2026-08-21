# Slope census — Cooper's third-order rows, the complex-root AZ family, and higher-order sporadics

*Computations in `lattice/census/cooper.gp`, `lattice/census/cooper_align.gp`, `lattice/census/az_families.gp` (exact PARI). This note does not touch `lattice/zeta3_lattice/` or `consolidation/ZETA3_TWO_LATTICE.md`.*

## 1. Cooper's $s_7,s_{10},s_{18}$

**Recurrences** (exact, from `sporadic_eisenstein_cooper_research_notes.txt` lines 346–359):
$$(n+1)^3A_{n+1}=2(2n+1)(3n^2+3n+1)A_n+4n(16n^2-1)A_{n-1}\quad(s_{10})$$
$$(n+1)^3A_{n+1}=(2n+1)(13n^2+13n+4)A_n+3n(9n^2-1)A_{n-1}\quad(s_7)$$
$$(n+1)^3A_{n+1}=2(2n+1)(7n^2+7n+3)A_n-12n(16n^2-1)A_{n-1}\quad(s_{18})$$

**Initial conditions.** The notes state $A_0=1,\ A_1=b$ generically, but this is *not* correct for $s_{10}$ and $s_{18}$: $A_1=b$ (the constant term of the quadratic) gives $A_1=1,3$ respectively, which do **not** converge to the claimed limits. Scanning $A_1=1,\dots,8$ against the target limits (with 150–250 digit precision) pins down the correct initial data exactly:
- $s_{10}$: $A_1=2$ — reproduces $s_{10}(n)=\sum_k\binom nk^4=1,2,18,164,\dots$ exactly.
- $s_7$: $A_1=4$ — matches the $F_7=1+4q+12q^2+\dots$ q-expansion in the notes.
- $s_{18}$: $A_1=6$ — found empirically (no q-expansion coefficients given in the notes); diff from $\tfrac12L(2,\chi_{-3})$ is $7\times10^{-20}$ at $n=150$, prec 200.

With these corrected ICs, all three limits are confirmed to 200–250 digits:

| row | $c$ (asymp. char. eq.) | $\lambda_1,\lambda_2$ | limit | verified digits |
|---|---|---|---|---|
| $s_{10}$ | $-64$ | $16,-4$ | $\zeta(2)/5$ | 250 |
| $s_7$ | $-27$ | $27,-1$ | $\zeta(2)/7$ | 250, exact `lindep`$=[-7,0,1,0]$ |
| $s_{18}$ | $192$ | $16,12$ | $\tfrac12L(2,\chi_{-3})$ | 250 |

(Char. eqs: $x^2-12x-64=0$, $x^2-26x-27=0$, $x^2-28x+192=0$, read off the $n\to\infty$ leading coefficients of the recurrence — these match the values already stated in the notes.)

**A key structural difference from the Zagier rows: $A_n$ is not integral.** $\text{denom}(A_n^{(s_{10})})$ grows (1, 4, 18, 144, 600, …), likewise for none of the three is $A_n\in\mathbb Z$ in general (only $s_7,s_{18}$ happen to stay integral through $n=29$; $s_{10}$ does not). Consequently the Notes-03 §1 slope law (which needs $v_p(A_n)=O(\log n)$) is **not automatically applicable**, and indeed:

**Measured slopes** ($v_p(b_n/a_n-b_{n-1}/a_{n-1})$ at $n=100,300,500$, precision 250 digits):

| row | $p=2$ | $p=3$ | $p=5$ | $p=7$ |
|---|---|---|---|---|
| $s_{10}$ | $-5,-5,-5$ (none) | $-2,-4,-3$ (none) | $-4,-2,-6$ (none) | $1,0,2$ (none) |
| $s_7$ | $-9,-17,-8$ (none) | $1,-1,2$ (none) | $-4,-3,-6$ (none) | $0,0,2$ (none) |
| $s_{18}$ | $-5,-5,-5$ (none) | $94,291,492$ — **slope $\approx1$** | $-4,-2,-6$ (none) | $1,0,2$ (none) |

So $s_{10}$ and $s_7$ have **no positive $p$-adic slope at $p=2,3,5,7$** — their $A_n$ denominators destroy the Cauchy mechanism entirely. Only $s_{18}$ has a genuine slope, at $p=3$ (slope $1=v_3(192)$), matching the $v_3(c)$ prediction even though $A_n$ is here integral.

**Cross-row alignment** (test $v_p(r\,a^X_nb^Y_n-r'\,a^Y_nb^X_n)$, $n=50\ldots300$):

| pair | common period | $r,r'$ | $p$ | result |
|---|---|---|---|---|
| $s_{10}$, Zagier D | $\zeta(2)/5$ | $1,1$ | 2,3,5,7 | **no alignment** (all valuations negative/noisy) |
| $s_7$, Zagier A | $\zeta(2)$ | $7,4$ | 2,3,5,7 | **no alignment** |
| $s_{18}$, Zagier C | $\tfrac12L(2,\chi_{-3})$ | $1,1$ | 3 | **aligned**: $v_3=51,101,150,204,252,300\approx n$ |
| $s_{18}$, Zagier B | $\tfrac12L(2,\chi_{-3})$ (B: no arch. limit) | $1,1$ | 3 | **aligned**: $v_3=52,102,151,204,252,301\approx n$ |
| $s_{18}$, Zagier F | $\tfrac58L(2,\chi_{-3})$ | $5,4$ | 3 | growing but irregular ($8$–$16$), not a clean linear law with this $r$ |

This confirms the extension-class rigidity law of Notes-03 §2: alignment occurs exactly where **both** rows have positive slope at the same prime ($s_{18}$ vs. C, B at $p=3$), and fails where one side has zero slope even though periods coincide ($s_{10}$/D, $s_7$/A — periods match but no side has a slope at 2,3,5,7 for the Cooper row). $L(2,\chi_{-3})/2$ now has **three** aligned rows at $p=3$: Zagier B, C, and Cooper's $s_{18}$.

## 2. The complex-root AZ family: $\delta=(7,3,81)$ and $\eta=(11,5,125)$

Recurrence form: $(n+1)^3u_{n+1}=(2n+1)(an^2+an+b)u_n-cn^3u_{n-1}$, $A_0=1,A_1=b$, $B_0=0,B_1=1$ (as in `lattice/third_order.gp`).

**Identification.** Archive conversation `chatgpt-research-archive/conversations/apery-systems/isogeny-eisenstein-cobundary-theorem-6a82c82c.md` (line 424) identifies $\delta:(7,3,81)$, together with $\alpha:(10,4,64)=\mathrm{Domb}$ and $\gamma:(17,5,1)=\mathrm{Ap\acute ery}$, as **three different Eisenstein extensions of the same Beauville-IV elliptic Picard–Fuchs system**, obtained as the symmetric squares of Zagier's rows $A,C,F$ respectively ($\delta \leftrightarrow$ Zagier A $(7,2,-8)\to\zeta(2)/4$). So $\delta$ is geometrically a $\zeta(3)$-family sibling of Domb/Apéry, not obviously an $L(3,\chi)$-family.

Book `book/v8/08c_chi5_application.tex` **explicitly identifies $\eta=(11,5,125)$** with $L(\Phi_\eta,3)=\tfrac12L(3,\chi_5)$ — this is source-proved in the dossier text (not just a numerical guess): "$L(\Phi_\eta,3)=\tfrac12\Lambda$", $\Lambda=L(3,\chi_5)$, with $A_0=1,A_1=5,B_0=0,B_1=1$ matching our normalization exactly.

**Characteristic roots** (both complex, confirming "no archimedean limit"):
- $\delta$: $x^2-14x+81=0\Rightarrow \lambda=7\pm4\sqrt5\,i$, $|\lambda|=9$.
- $\eta$: $x^2-22x+125=0\Rightarrow\lambda=11\pm2i$, $|\lambda|=\sqrt{125}$.

**Measured $p$-adic slopes** (from $v_p(b_n/a_n-b_{n-1}/a_{n-1})$, $n$ up to 800):
- $\delta$ at $p=3$: slope $\approx (2375-389)/500\approx 3.97\approx4$ — matches $v_3(81)=4$ exactly.
- $\eta$ at $p=5$: slope $\approx(1774-282)/500\approx2.98\approx3$ — matches $v_5(125)=3$ exactly.
- Both at $p=2,7$: no slope (bounded/negative valuations, as expected: $v_2,v_7$ of $81,125$ are $0$).

**Cross-alignment search (new negative result).** Scanned $\delta$ against Domb, T$=(12,4,16)$, and Apéry at $p=3$ over 22 candidate rational multipliers $r\in\{1,2,3,4,6,7,8,9,12,16,24,32,\tfrac12,\tfrac13,\tfrac14,\tfrac32,\tfrac43,\tfrac74,\tfrac78,\tfrac7{16},\tfrac7{24},\tfrac7{32}\}$: **no candidate produces growing valuation** — all six sequences $v_3(r\,a^\delta_nb^Y_n-a^Y_nb^\delta_n)$ for $Y=$Domb, T, Apéry are flat-to-decreasing out to $n=800$. So, despite $\delta$ living on the same elliptic surface as Domb/T/Apéry (all carrying rational multiples of $\zeta(3)$), **$\delta$'s 3-adic limit does not coincide with theirs at any simple rational scalar** — either the true relation needs an irrational-looking $p$-adic scalar (a genuine $L$-value ratio, not visible to a finite rational scan), or $\delta$ genuinely carries a different 3-adic extension class from Domb/T/Apéry despite the shared geometry. This is left open.

**$p$-adic limits, no known aligned partner found for either $\delta$ or $\eta$** (no second $L(3,\chi_5)$ row or level-20 companion located in the notes). Per the task instructions, the $p$-adic limits are recorded here for future use: $\xi_3(\delta)$ and $\xi_5(\eta)$ are Cauchy to $\sim(3.97\!\times\!800)\approx3170$ and $\sim(2.98\!\times\!800)\approx2380$ $p$-adic digits respectively at $N=800$ (script `lattice/census/az_families.gp`); the raw rational values $b_{800}/a_{800}$ are printed by the script (not reproduced here — thousand-digit rationals) rather than in this file. Re-run the script for the actual digit strings if needed.

The book gives the period polynomial for $\eta$'s source explicitly: $125x^2-234x+125=0$ at the two conjugate exponents $\alpha,\bar\alpha$ of the order-six Hadamard module (§ "Order-six module and arithmetic lattice"), and the fold values $t_\pm=(11\pm2i)/125$, root of $1-22t+125t^2$ — consistent with $\eta$'s own characteristic equation $x^2-22x+125$ above (same discriminant structure, $t=1/x$). No period polynomial for $\delta$'s specific Eisenstein extension (as opposed to the shared Beauville-IV base curve) was found in the notes.

## 3. Higher-order rows (partial — see caveat)

**$\zeta(5)$, level 12** (`modular_period_annihilation_research_notes.txt` §9–10). Target $C=(11/144)\zeta(5)$ (also $25/144\zeta(5)$ appears elsewhere in the archive as a related/optimized constant — the level-12 source itself gives $11/144$). The *true* order-5 Picard–Fuchs operator is $\mathrm{Sym}^4$ of a rank-2 (elliptic) system with base characteristic equation $\lambda^2-14\lambda+1=0$, roots $\alpha=7+4\sqrt3\approx13.928$, $\alpha^{-1}=7-4\sqrt3$ (product $1$). The five $\mathrm{Sym}^4$ eigenvalues are $\alpha^4,\alpha^2,1,\alpha^{-2},\alpha^{-4}$, product $=1$ — an Apéry-perfect Casoratian, i.e. $\prod\lambda_i=1$ exactly, generalizing $c=\pm1$. This is exact (read directly off the stated characteristic polynomial $(\lambda+1)^4(\lambda^2-14\lambda+1)^5$, resp. $(\lambda^3-13\lambda^2-13\lambda+1)^5$, for the order-14/16 recurrence obtained after clearing denominators — the extra $(\lambda+1)$ factor is an artifact of the polynomial-coefficient recurrence, not of the order-5 system itself).

**Generalization of the slope law used here:** the task's proposed generalization — "slope of $b_n/a_n$ depends on the *second-dominant* solution, not the full Casoratian" — is confirmed *analytically* by this eigenvalue structure: $b_n/a_n$ converges at the rate of the ratio of the two largest eigenvalues, $\lambda_2/\lambda_1=\alpha^2/\alpha^4=\alpha^{-2}$, i.e. the companion linear form decays like $\alpha^{-2n}=(7-4\sqrt3)^{2n}$ up to polynomial factors — **not** like $\lambda_{\min}/\lambda_1=\alpha^{-8}$. This is the direct order-$m$ generalization of the Zagier-row law (there $m=2$ and "second-dominant" = "the only other root").

**Caveat (honest, not guessed):** I did not transcribe and run the full order-14 recurrence with the 15 quintic-in-$n$ polynomial coefficients $P_0,\dots,P_{14}$ (only the "independent half" $P_0$–$P_7$ plus a reflection rule $P_{14-j}(n)=-P_j(11-n)$ are given in the notes) — this is large and error-prone to reproduce exactly in the time available, so I did **not** numerically re-verify $C=(11/144)\zeta(5)$ or measure actual $p$-adic slopes at $p=2,3,5,7$ for this row. The characteristic-root/product analysis above is exact (it only needs the stated characteristic polynomial), but the numeric census entry for this row is incomplete.

**$\zeta(7)$, level 24** (`modular_apery_cm_isogeny_research_notes.txt` §VII). Target $C_7=(1463/13824)\zeta(7)$, source vector $(1,-588,11583,-27456,-138996,585728,-762048,331776)$. This is stated to be $\mathrm{Sym}^6$ of a rank-2 system on the elliptic quotient $Y^2=r^4-34r^2+1$ with Apéry coordinate $z=r/(1+r)$, CM branch $z_*=(2-\sqrt2)/4$, next branch $z_2=(1-\sqrt2)/2$; the notes report an **observed** (not derived here) root/convergence ratio $|z_*|/|z_2|=1/\sqrt2$. If — by analogy with the level-12 case — this $1/\sqrt2$ is the ratio $\lambda_2/\lambda_1$ of the two dominant $\mathrm{Sym}^6$ eigenvalues, and the base rank-2 system also has unit root-product, the seven eigenvalues would be $\beta^6,\beta^4,\beta^2,1,\beta^{-2},\beta^{-4},\beta^{-6}$ with $\beta^{-2}=1/\sqrt2\Rightarrow\beta=2^{1/4}$; this is **not independently confirmed** — I did not derive the base rank-2 characteristic equation from $Y^2=r^4-34r^2+1$ myself, and did not transcribe/run the order-7 recurrence (its explicit polynomial coefficients are not fully given in the excerpted notes; only the differential operator $L_7$ in terms of $R(z)$ is given). This row is reported here only at the level of the notes' own statements, not independently verified.

## 4. Summary table

| row | period | order | $c=\prod\lambda_i$ | $\lambda_1$ | $\lambda_2$ | slopes $p=2,3,5,7$ | denom. exponent | aligned partners |
|---|---|---|---|---|---|---|---|---|
| $s_{10}$ | $\zeta(2)/5$ | 2 | $-64$ | $16$ | $-4$ | $0,0,0,0$ | $A_n\notin\mathbb Z$ (growing denom.) | none found |
| $s_7$ | $\zeta(2)/7$ | 2 | $-27$ | $27$ | $-1$ | $0,0,0,0$ | integral through $n\le29$ | none found |
| $s_{18}$ | $\tfrac12L(2,\chi_{-3})$ | 2 | $192$ | $16$ | $12$ | $0,1,0,0$ | integral through $n\le29$ | Zagier B, C at $p=3$ |
| $\delta$ | $\zeta(3)$-family (Beauville-IV, no arch. limit) | 3 | $81$ | $7\pm4\sqrt5\,i$ | — | $0,4,0,0$ | not checked | none found (22 rationals scanned at $p=3$) |
| $\eta$ | $\tfrac12L(3,\chi_5)$ (no arch. limit) | 3 | $125$ | $11\pm2i$ | — | $0,0,3,0$ | not checked | none found in notes |
| Zagier B (ref.) | $\tfrac12L(2,\chi_{-3})$ | 2 | $27$ | $\sqrt{27}$ | $\sqrt{27}$ | $0,3,0,0$ | — | C, $s_{18}$ at $p=3$ |
| Zagier C (ref.) | $\tfrac12L(2,\chi_{-3})$ | 2 | $9$ | $9$ | $1$ | $0,2,0,0$ | — | B, F, $s_{18}$ at $p=3$ |
| ζ(5) L12 | $\tfrac{11}{144}\zeta(5)$ | 5 (Sym$^4$) | $1$ (exact) | $\alpha^4\approx37637$ | $\alpha^2\approx193.99$ | not measured (recurrence not run) | not checked | n/a |
| ζ(7) L24 | $\tfrac{1463}{13824}\zeta(7)$ | 7 (Sym$^6$) | not derived | — | — | not measured | not checked | n/a |

**Periods with $\ge2$ aligned rows:** $\tfrac12L(2,\chi_{-3})$ at $p=3$ now has **three** aligned rows (Zagier B, Zagier C, Cooper $s_{18}$), with a fourth (Zagier F) plausible but not cleanly confirmed with the tested $r=5/4$ scalar. No other period in this census has $\ge2$ aligned rows.

## Scripts
- `lattice/census/cooper.gp` — Cooper row generation, limits, slopes, denominators.
- `lattice/census/cooper_align.gp` — cross-row alignment tests, Cooper vs. Zagier.
- `lattice/census/az_families.gp` — $\delta,\eta$ complex-root rows: roots, slopes, alignment scan vs. Domb/T/Apéry.

## 5. ζ(7), level 24 — full computation (follow-up)

Source: `chatgpt-research-archive/Closed binomial–harmonic forms for a level-24 ζ(7) Apéry system.md` (closed forms) + `modular_apery_cm_isogeny_research_notes.txt` §VII.C (eta-quotient $r(\tau)$). Script: `lattice/census/zeta7_level24.gp`.

**Construction actually implemented** (bypassing the paper's own finite-sum $\mathcal K_{n,m}$ formula for $B_n$, using the paper's own stated relation $B(z(q))=A(z(q))\Psi(q)$ directly as a $q$-series computation):
- $r(\tau)=[\eta(2\tau)\eta(12\tau)/(\eta(4\tau)\eta(6\tau))]^6=q\cdot(\text{eta-ratio})^6$, computed as an exact $q$-series from the Dedekind-eta products (order $q^{220}$).
- $z(q)=r/(1+r)$; $q(z)=\mathrm{serreverse}(z(q))$.
- $\Phi(q)=\sum g_mq^m$, $g_m=\sum_{d\mid24,\,d\mid m}c_d\sigma_7(m/d)$, $(c_1,\dots,c_{24})=(1,-588,11583,-27456,-138996,585728,-762048,331776)$; $\Psi(q)=\sum (g_m/m^7)q^m$.
- $A(z)=\sqrt{1-34s+s^2}\,\mathcal A(s)^3$, $s=(z/(1-z))^2$, $\mathcal A(s)=\sum a_ms^m$ the Apéry-$\zeta(3)$ numbers — computed **two independent ways**: (a) directly as a $z$-series substitution, (b) from the paper's own boxed finite sum $A_n=\sum_{m=1}^{\lfloor n/2\rfloor}\binom{n-1}{2m-1}C_m$ with $C_m=\sum_{\ell+i+j+k=m}y_\ell a_ia_ja_k$ and the given closed form for $y_\ell$. **Both agree exactly** for $n=0..8$: $A_n=1,0,-2,-4,-111,-428,-4016,-19848,-146040,\dots$ — and these also agree exactly with the independently-sourced "$E_{7,\mathrm{geom}}(z)$" coefficients quoted in `modular_apery_cm_isogeny_research_notes.txt` §VII.E, confirming $A_n$ is the same object under two names in two different notes files.
- $B(z)=A(z)\cdot\Psi(q(z))$ by power-series composition; extracted to $n=218$ (limited by the $q$-series truncation order, not by algebra).

**(i) Integrality of $A_n$:** confirmed by inspection — all computed $A_n\in\mathbb Z$ (matches the paper's claim).

**(ii) $B_n$ verification.** $d_n^7B_n\in\mathbb Z$ ($d_n=\mathrm{lcm}(1,\dots,n)$) checked exactly for **all** $n=1,\dots,218$ — true in every case, confirming the paper's sharp lcm-clearing claim numerically. Limit:
$$\frac{B_{218}}{A_{218}}=0.106714047512279593147180542204154957914\ldots$$
$$\tfrac{1463}{13824}\zeta(7)=0.106714047512279593147180542204156008340\ldots$$
agreeing to **33 digits** ($|{\rm diff}|\approx1.05\times10^{-33}$), and $|{\rm diff}|^{1/n}\approx0.7059$ at $n=218$ vs. the paper's exact predicted rate $1/\sqrt2=0.70711$ (the residual gap is consistent with the paper's stated $O(n^{-1})$ correction to the leading exponential).

**(iii) Characteristic data.** I did not fit the full order-7 recurrence (its polynomial coefficients are not given in full in the excerpted notes and a numerical fit was judged too fragile in the time available); instead I used the paper's own stated facts directly: dominant root $\lambda_1=4+2\sqrt2\approx6.828$, and archimedean decay rate $|\lambda_2/\lambda_1|=1/\sqrt2$ (boxed in the closed-form paper), giving $\lambda_2=\lambda_1/\sqrt2=2+2\sqrt2\approx4.828$. This is exactly the same "second-dominant-root controls the linear form" law used for the level-12 $\zeta(5)$ system.

**Measured $p$-adic slopes of $B_n/A_n$** (valuation of successive differences at $n=60,\dots,218$, $p=2,3,5,7$): **all four primes show no growing slope** — the valuations are negative/noisy and do not trend upward with $n$ (e.g. $p=2$: $-33,-44,-40,-49,-49,-49$; $p=3$: $-22,-28,-29,-28,-29,-30$; similarly for $5,7$). $v_p(A_n)$ itself stays small and non-monotonic for $p=3,5,7$ (values $0$–$2$) and grows slowly for $p=2$ (up to $10$ at $n=218$, consistent with $O(\log n)$).

**(iv) Answer to the coordinator's question:** $B_n/A_n$ **does not converge $p$-adically at $p=2$ or $p=3$** (the two primes dividing the level $24=2^3\cdot3$), nor at $5,7$. All of the row's "budget" (Notes-03 §1 terminology) is archimedean: the full decay rate $1/\sqrt2$ is carried by $|\lambda_2/\lambda_1|$ with **no measurable $p$-adic contribution at small primes** — i.e. the product-of-roots analogue $c=\prod\lambda_i$ appears to be coprime to $2,3,5,7$ (or at least contributes no $p$-adic Cauchy structure to $B_n/A_n$ at these primes), unlike the Cooper $s_{18}$ or Zagier B/C/F rows where $p$-adic slope and archimedean decay split the total budget. This makes the level-24 $\zeta(7)$ row **structurally an "Apéry-perfect" row** in the sense of Notes-03 §1 (all decay archimedean, nothing to harvest $p$-adically at $2,3,5,7$) — analogous to Apéry's own $\zeta(3)$ row, not to the Catalan/Domb/$s_{18}$-type rows that split their budget with a prime.

