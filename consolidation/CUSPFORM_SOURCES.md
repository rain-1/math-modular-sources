# Cusp-form sources on the Fricke hosts: the Fricke-sign criterion, five new Apéry-like rows, no irrationality

*Fable, 2026-09-02. Closure of the "cusp-form sources" probe (bucket 5, item 3). Full report and 22 scripts: `lattice/cuspform_sources/REPORT.md`.*

## 1. The question

The companion formula $B_n=\sum_{m\le n}c(m)m^{-r}e_{n,m}$ gives $d_n^{\,r}B_n\in\mathbf Z$ for **any** source with integer coefficients; Eisenstein sources were only ever needed for fold-regularity. So: do weight-four cusp forms on the Fricke hosts ($F=D\log u$ of weight 2, $r=3$) give Apéry-like linear forms in $1$ and $L(f,3)$, a critical value of a rank-two motive whose period is not known to be irrational?

## 2. Theorem A (the criterion) [exact; calibrated on Apéry's host]

Let $u|W_N=1/(Cu)$, $x=u/(1+Bu+Cu^2)$, $F=D\log u$, so $x$ is $W_N$-invariant with a critical point (the fold) at $\tau_c=i/\sqrt N$, $F|_2W_N=-F$, and the canonical source $\Phi_0=F\cdot Dx=xF^2\sqrt{P(x)}$ satisfies $\Phi_0|_4W_N=-\Phi_0$. For $\Phi\in M_4(\Gamma_0(N))$ with $a_0(\Phi)=0$ and $\Theta=D^{-3}\Phi$: near $\tau_c$ the singular part of $F$ is its odd part in $s=\tau-\tau_c$, so $B-\xi F=F(\Theta-\xi)$ is analytic at $t_1$ iff $F(\Theta-\xi)$ is $W_N$-invariant, i.e.
$$(\star)\qquad \Theta|_{-2}W_N+\Theta=\xi\,(N\tau^2+1).$$
With $\Theta|_{-2}W_N=\varepsilon(\Theta+4\pi^3iR)$ for $\Phi|_4W_N=\varepsilon\Phi$, $R(\tau)=i\Lambda(1)\tau^2+2\Lambda(2)\tau-i\Lambda(3)$, $\Lambda(j)=\int_0^\infty y^{j-1}\Phi(iy)\,dy$:
* $\varepsilon=+1$ makes $(\star)$ read $2\Theta+P=$ quadratic, impossible for $\Phi\ne0$;
* $\varepsilon=-1$ makes $(\star)$ read $\Lambda(2)=0$ and $\Lambda(1)=-N\Lambda(3)$, both automatic from $\Lambda(s)=\varepsilon N^{2-s}\Lambda(4-s)$, and then $\xi=4\pi^3\Lambda(3)=L(\Phi,3)$.

**So $\Phi$ is fold-regular on every Fricke host of level $N$ iff $\Phi|_4W_N=-\Phi$, with Apéry limit $L(\Phi,3)$.** Calibration: Apéry's source $E_4-28E_4(2\tau)+63E_4(3\tau)-36E_4(6\tau)$ has $W_6=-1$ and $L(\Phi,3)=\zeta(3)/6$ (60 digits). For $\varepsilon=+1$ the limit still exists but the convergence is $\xi+c/n+O(n^{-2})$ (measured ratio $0.99000$ against $((n-1)/n)^2=0.990025$ at $n=200$ on all twelve hosts), the limit is not $L(\Phi,3)$, and no Eisenstein admixture repairs it because $M_4=\mathrm{Eis}\oplus S_4$ is $W_N$-stable.

**Corollary (host-independence).** $F'(\tau_c)=N\tau_cF(\tau_c)$, hence $\xi=\Theta(q_c)+(2\pi/\sqrt N)\,D\Theta(q_c)$ with $q_c=e^{-2\pi/\sqrt N}$: the host fixes rate and denominators only. For $\varepsilon=-1$ this is the classical Fricke-symmetric series $L(\Phi,3)=\sum a(m)m^{-3}e^{-2\pi m/\sqrt N}+\tfrac{2\pi}{\sqrt N}\sum a(m)m^{-2}e^{-2\pi m/\sqrt N}$ (96 digits against `lfun`).

## 3. The census [verified]

$\dim S_4(\Gamma_0(N))^{W_N=-1}=0$ for $N=5,6,7,8,9$: the newforms $5.4.a.a,\ldots,9.4.a.a$ all have Fricke sign $+1$. **No cusp form is fold-regular on either Apéry-perfect host** (Apéry's $N=6$, Cooper's $N=7$). Five fold-regular (host, form) pairs exist, all oldform combinations $f-d^2f(d\tau)$:

| host | $N,B,C$ | source | limit $\xi$ | rate | $K_-$ |
|---|---|---|---|---|---|
| Cooper $s_{10}$ | $10,6,25$ | $f_5-4f_5(2\tau)$ | $L(5.4.a.a,3)/2$ | $1/4$ | not identified |
| Domb | $12,10,9$ | $f_6-4f_6(2\tau)$ | $L(6.4.a.a,3)/2=0.36455672878339462975728560797\ldots$ | $1/4$ | $L(6.4.a.a,2)/\sqrt\pi$ |
| level 12, $C=1$ | $12,34,1$ | $f_6-4f_6(2\tau)$ | $L(6.4.a.a,3)/2$ | $8/9$ | not identified |
| Cooper $s_{18}$ | $18,14,1$ | $f_6-9f_6(3\tau)$ | $2L(6.4.a.a,3)/3$ | $3/4$ | $\sqrt3\,L(6.4.a.a,2)/\sqrt\pi$ |
| Cooper $s_{18}$ | $18,14,1$ | $f_9-4f_9(2\tau)$ | $L(9.4.a.a,3)/2$ | $3/4$ | not identified |

Domb host row: $A_n=1,4,28,256,2716,\ldots$ (Domb numbers), $B_n=0,1,\tfrac{37}4,\tfrac{818}9,\tfrac{141587}{144},\ldots$, $d_n^{\,3}B_n\in\mathbf Z$ with $k=3$ sharp (no free integration anywhere, even on the free-integration hosts), $A_n\xi-B_n\sim K_-\lambda_2^nn^{-3/2}$, 144 digits at $n=240$. Where $K_-$ is identified it is the **central** value $L(f,2)$ of the newform, nonzero because the newform's own sign is $+1$. The far singularity is a square-root branch at the far fold with coefficient $-2\sqrt\pi K_-$. $B_n$ satisfies no three-term recurrence (minimal order 3 or 4, characteristic roots $\lambda_1,\lambda_2,\lambda_2$), so the Casoratian identity $K_+K_-(\lambda_1-\lambda_2)=\kappa$ of `ASYMPTOTIC_CONSTANTS.md` fails for these rows. Closed form in three of five cases: $\Phi=xF^2\sqrt{1-\lambda_1x}=\Phi_0/\sqrt{1-\lambda_2x}$, equivalently the host's recurrence leaves the exact residual $\binom{2n+2}{n+1}(\lambda_2/4)^{n+1}$.

**Mixed periods.** The full fold-regular source space is one-dimensional on Apéry's host ($\zeta(3)/6$ only), three-dimensional at level 12 ($\zeta(3)$, $L(6.4.a.a,3)$) and five-dimensional at level 18 ($\zeta(3)$, $L(3,\chi_{-3})$, $L(6.4.a.a,3)$, $L(9.4.a.a,3)$).

## 4. Diophantine verdict

None. Only Apéry's host has $\log|\lambda_2|+3<0$ (namely $-0.525$), and it carries no cusp source; the best host that does (Domb, $s_{10}$) misses the entry condition by $4.386$ nats. The objects are new and explicit (Apéry-like approximations with lcm-cubed denominators to critical values of non-CM weight-four newforms), and the arithmetic side is untested: hypothesis (a) of the companion Lucas theorem is false for a Hecke eigenform, so the Lucas law is expected to fail.

## 5. Not done

The weight-one hosts (cusps, not folds) need a different argument, and the only level-7 weight-one four-point host has complex conjugate dominant singularities, so the CM form $\eta(\tau)^3\eta(7\tau)^3$ has nowhere to act; three of five $K_-$ unidentified; meromorphic $W_N$-anti-invariant sources outside Theorem A; no $p$-adic analysis.
