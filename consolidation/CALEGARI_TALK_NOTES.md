# Notes on Calegari's talk "The arithmetic of some Dirichlet L-values"

*Fable, 2026-08-21. Source: `chatgpt-research-archive/calegari-talk.txt` (transcript reformatted by River). What the talk says, and where it touches our programme.*

## Content of the talk (compressed)

1. **The list of known irrationality results** for $L(m,\chi)$, $\chi$ quadratic: Lambert ($\pi$), Legendre ($\pi^2$), Lindemann (all right-parity values, via $\pi$), Gelfond–Schneider ($m=1$, all characters, via logarithms), Apéry ($\zeta(3)$), then nothing for 45 years (Rivoal–Zudilin give infinitely many but none specific). Main theorem (CDT): $L(2,\chi_{-3})\notin\mathbb Q$, indeed $1,\pi^2,L(2,\chi_{-3})$ linearly independent over $\mathbb Q$. Corollary: $\sum_{n\equiv\pm1\ (3)}\pm n^{-2}$ irrational. Figure-eight volume $=\tfrac{3\sqrt3}{2}L(2,\chi_{-3})$: first contact with Thurston's last open problem (volumes not all rationally related).
2. **Explicit Catalan remark:** "The next smallest discriminant would be $\mathbb Q(\sqrt{-1})$, and you would have Catalan's constant here, but alas our argument does not work in this case."
3. **Apéry via Beukers/Siegel $G$-functions:** recurrence $\Rightarrow$ order-3 ODE with singularities $0,\infty,\alpha,\beta$ ($\alpha=(\sqrt2-1)^4$ small); codimension-one holomorphy at $\alpha$ $\Rightarrow$ some combination $A-\zeta(3)B$ overconverges to $\beta$. Three ingredients: (i) small denominators, (ii) why $\zeta(3)$, (iii) "why is 33 bigger than 20". (i),(ii) understood conjecturally via Bombieri–Dwork ("integer coefficients $\Rightarrow$ from geometry"; Legendre family example); (iii) is the real mystery; generalisation attempts have always failed on (iii): the second singularity $\beta$ must exceed $e^{k}$ and usually doesn't.
4. **Exploit more than the disc.** Pólya: integer coefficients + conformal radius $>1$ $\Rightarrow$ rational ($4\beta$ for the slit plane). André: integer coefficients + $|\varphi'(0)|>1$ for a non-univalent uniformising map $\Rightarrow$ algebraic; Kodaira's modular map gives $16\beta$ on $\mathbb P^1\setminus\{0,\beta,\infty\}$, optimal ($\sum\binom{2n}n^2x^n$, singular at $1/16$, is transcendental). Proof sketch: Siegel's lemma + Cauchy estimate of the leading coefficient.
5. **Denominators.** André-type theorems extend to finitely many primes ("conformal derivative times the product of $p$-adic radii of convergence $>1$") but fail for $\log(1-x)$ (all $p$-adic radii $=1$, $16>1$, yet transcendental). Fix: replace powers by derivatives (holonomy, not algebraicity), denominators of LCM type $d_n^M$. Then $16>e$ gives finite dimension for $d_n$-denominators on $\mathbb P^1\setminus\{0,1,\infty\}$.
6. **Quantitative holonomy bounds, three versions.** (I) $\dim<\infty$ if $\log|\varphi'(0)|>M$. (II) replace max by average: many variables, concentration of measure, Vandermonde determinant kills the small-measure bad set (2012–13; used for unbounded denominators). (III) $M\to M^\flat$: exploit a subset of functions with better denominators — "more of an arithmetic nature". Specialising (II) to univalent $\varphi$, $M=0$ recovers Pólya with constant exactly $1$ (Jensen), so it is sharp there.
7. **Qualitative is useless here:** the functions already satisfy ODEs; the contradiction must be "ODE of order $\le r$" against a monodromy computation. Hence the need for near-optimal constants.

## Where it touches our programme

| talk | ours |
|---|---|
| Catalan: "our argument does not work" | `ADELIC_HOLONOMY.md`, paper §6.5: margin $-8.0$ (level 8), $-7.55$ (level 16, two classes), structural; entry condition met for the first time via the $2$-adic slope. No number is given in the talk; ours is the quantitative form. |
| ingredient (ii), Bombieri–Dwork, "periods of degenerate fibres" | Theorems A/B*/F: the Apéry limit is the regulator of the Eisenstein class at the place where the cusp disc has a free boundary (`THEORY_NOTES_05/06`). Integral MUM variation = "comes from geometry". |
| ingredient (iii), the coincidence | the score $\log(1/\lambda_2)-k$ and River's near-miss principle; scalar barrier. |
| finitely-many-prime $p$-adic radii version of André | Theorem (adelic holonomy bound) $\tau_{\rm ad}=\tau-\sum_p\gamma_p$: LCM denominators *and* $p$-adic overconvergence at the bad primes; $\gamma_p$ = log of a $p$-adic radius. Consistent with the talk. |
| version (III) $M^\flat$ | our $\tau^\flat/\tau^\sharp$ split in `CDT_FINDER.md`; $\sum_p\gamma_p$ is a fourth, $p$-adic, refinement of the same arithmetic kind. |
| Kodaira $16\beta$; $\lambda(\tau)$ geometry | level-16 Catalan host: singularity at relative radius $1/2$ (the $1/16$ of $\lambda$); no CDT contour loss transports (`ADELIC_HOLONOMY.md`). |
| dimension bound vs monodromy | our margin in log units $\approx$ log of the surplus dimension; "fourteen functions needed where Eisenstein supplies one". |
| Thurston's problem | Catalan $=\tfrac14\mathrm{vol}$(Whitehead link complement); same flavour. One sentence for §7. |

Not in the talk: $L(3,\chi_5)$, Beukers' mixed statements, simultaneous/vector-valued approximation, anything $p$-adic beyond the radii remark.

## Two things to take from it

1. In the paper, describe $\sum_p\gamma_p$ as "a refinement in the direction of CDT's $\tau^\flat$: rewarding functions that are better at a $p$-adic place instead of at the archimedean denominator count". That is the honest placement.
2. Calegari's framing of (iii) as the sole obstruction supports spending effort on *finding* hosts with larger second singularity relative to Hodge depth (the MUM survey / non-congruence scan), rather than on strengthening the bound further for fixed hosts: for Catalan the deficit is $7.5$, and bound-side refinements have so far produced $+0.45$.
