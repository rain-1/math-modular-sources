# The identification of $\gamma$ (row $s_7$): a twisted CM trace of a Maass-raised
# weight $-2$ form

## The form

On $\Gamma_0(7)$ put, with $u=(\eta(7\tau)/\eta(\tau))^4$, $g(u)=1+13u+49u^2$, $x=u/g(u)$,
$F=D\log u\in M_2(\Gamma_0(7))$,
$$\boxed{\;f:=\frac1{x\,F}=\frac{g(u)}{u\,F}\;}
= q^{-1}+5+19q+52q^2+137q^3+316q^4+695q^5+1440q^6+2887q^7+5564q^8+\cdots\in q^{-1}\mathbf Z[[q]].$$

**Lemma.** $f$ is a *weakly holomorphic modular form of weight $-2$ on $\Gamma_0(7)$*:
it is holomorphic on $\mathbf H$ and $f|_{-2}W_7=-f$.  Indeed $x$ has a pole of order
$\nu=3$ (in $\tau$) at the two CM points $\tau_0$ of discriminant $-3$, while $F$ vanishes
there to order $2$ (the divisor of $F$ on $X_0(7)$ has degree $2\cdot8/12=4/3$ and is
supported on the two elliptic points of order $3$); so $f=(1/x)(1/F)$ vanishes to order
$3-2=1$ at $\tau_0$ and has no other zero or pole on $\mathbf H$.  $x|W_7=x$ and
$F|_2W_7=-F$ give the Fricke sign.  Up to scalars $f$ is the *unique* weight-$(-2)$
weakly holomorphic form on $\Gamma_0(7)$ with pole order $\le1$ at both cusps. $\square$

Let $R_{-2}=2i\,\partial_\tau-2y^{-1}$ be the Maass raising operator ($y=\operatorname{Im}\tau$)
and put
$$\widehat f:=-\frac{R_{-2}f}{4\pi}=Df+\frac{f}{2\pi y},\qquad D=q\frac{d}{dq},$$
a **real-analytic $\Gamma_0(7)$-invariant function of weight $0$**
(checked numerically to 50 digits against $\binom{1\ 0}{7\ 1}$ and $\binom{8\ 1}{7\ 1}$;
`t3.gp`).  Its CM values are algebraic (Shimura–Maass); e.g.
$\widehat f\bigl(\tfrac{-5+\sqrt{-3}}{14}\bigr)=-\sqrt{-3}$.

## The identity

For $d=-3m^2$ let $\mathcal Q_m$ be the set of $\Gamma_0(7)$-classes of Heegner forms
$[A,B,C]$, $7\mid A$, $B\equiv5m\pmod{14}$, $B^2-4AC=-3m^2$ (the GKZ convention
$\beta=\delta_0\sqrt d$ with $\delta_0=5$, $\delta_0^2\equiv-3\bmod 28$), let
$\chi_{-3}$ be the genus character and $\omega_Q$ the order of the stabiliser of
$\alpha_Q$ in $\overline{\Gamma_0(7)}$.  Define the twisted trace
$$\operatorname{Tr}_m(\widehat f):=\sum_{Q\in\mathcal Q_m}\chi_{-3}(Q)\,\frac{\widehat f(\alpha_Q)}{\omega_Q}.$$

$$\boxed{\;\beta_{s_7}(m)\;=\;\sqrt{-3}\;\operatorname{Tr}_m(\widehat f)\;}
\qquad\text{equivalently}\qquad
\gamma_{s_7}(m)=\frac{\sqrt{-3}}{m^2}\operatorname{Tr}_m(\widehat f),\quad
a(m^2)=m\,\sqrt{-3}\operatorname{Tr}_m(\widehat f).$$

**[verified]** for every $m\le45$ with $7\nmid m$, to $60$ decimals (worst absolute error
$4\cdot10^{-61}$ against values of size up to $2.3\cdot10^{15}$, i.e. $\ge70$ significant
digits); the number of $\Gamma_0(7)$-Heegner classes found equals the number of
$\mathrm{SL}_2(\mathbf Z)$-classes of discriminant $-3m^2$ in every case
(`17_ident7wide.gp/.log`, and `13_ident7.gp/.log`).  For $7\mid m$ the two $\beta$-classes
merge, the trace becomes real, and $\beta(m)=0$ by Theorem 6.3 of round 1.

$\widehat f$ was computed two independent ways which agree to 50 digits (`maass.gp`,
`maass2.gp`): a closed form
$\widehat f=-f^2(F\,Dx+x\,DF)+f/(2\pi y)$ with $Dx=Fu(1-Cu^2)/g(u)^2$,
$DF=\sum_d\lambda_dd\,(E_2(d\tau)^2-E_4(d\tau))/12$ and $E_2^*$, $E_4$ evaluated by
reduction to the $\mathrm{SL}_2(\mathbf Z)$ fundamental domain; and the $q$-expansions of
$f$, $Df$ evaluated after reduction in $\Gamma_0(7)+7$ (used at the zeros of $F$, i.e. at
the polar CM points).

## What it buys

1. **$\gamma$ is identified.**  $K=\sum\gamma(m)q^m$ is $\sum_m m^{-2}\sqrt{-3}\operatorname{Tr}_m(\widehat f)q^m$.
2. **The master conjecture becomes a divisibility of traces of singular moduli:**
   $$n^2\mid\beta(n)\iff n^2\ \Bigm|\ \sqrt{-3}\operatorname{Tr}_n(\widehat f),$$
   i.e. exactly the Ahlgren–Ono/Edixhoven phenomenon, for the *first Shimura–Maass
   derivative of a weight $-2$ form* rather than for a modular function.  Route (B) of the
   task is now a precisely stated problem, and the $p$-adic input is available: on the
   ordinary locus the Shimura–Maass derivative $R_{-2}$ agrees, at CM points, with the
   Katz $\theta$-operator under the Katz $p$-adic period comparison, and the disc-$D_0p^{2k}$
   CM points lie in the $p$-adic residue discs of the disc-$D_0$ points.
3. **It explains the failure of the naive route** (`NOTES_TRACE.md`): the raising operator
   contributes the term $f/(2\pi y)$, which at the top Heegner point is
   $q^{-1}/(2\pi m\operatorname{Im}\tau_0)=q^{-1}\kappa/m$ — exactly the $1/m$ tail that no
   modular *function* can produce.
4. **It matches the Shimura–Borcherds picture** (`COOPER_CONGRUENCE.md` §4): the theta lift
   of a weakly holomorphic form of weight $2-2k$ has weight $k+1/2$; for $k=2$ the input has
   weight $-2$ and the output weight $5/2$.  So $\sum_d\operatorname{Tr}_d(\widehat f)q^d$
   should be (a normalisation of) the weight-$5/2$ Kohnen-plus input $f_{5/2}$ at level $28$
   whose Shimura–Borcherds lift is $\Phi$, and its square-index coefficients are
   $a(m^2)=m\beta(m)$.
