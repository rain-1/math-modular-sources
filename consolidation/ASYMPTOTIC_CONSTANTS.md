# Asymptotic constants of modular Apéry sequences: elementary unless the companion vanishes at the fold

*Fable, 2026-09-02. Theory here; census numbers in `scratchpad/asymK/REPORT.md` (agent, pending) and then in §4. Generalises the closed form of the level-12 $\zeta(5)$ constant $K$ (`ZETA5_K_CLOSED_FORM.md`).*

## 1. Setting

A modular Apéry row $(x,F)$: $x(q)=q+O(q^2)$ an integral Hauptmodul-type coordinate on a genus-zero group $\Gamma\supseteq\Gamma_0(N)$, $F$ a modular form of weight $k$ (integral $q$-expansion, $F=1+O(q)$), $A(x)=F(q(x))=\sum a_nx^n$. Suppose the dominant singularity of $A$ is a **simple fold** of $x$ at the Fricke point $\tau_c=i/\sqrt N$: $x'(q_c)=0\ne x''(q_c)$, $q_c=e^{-2\pi/\sqrt N}$, $x_+=x(q_c)$, and $F$ is analytic at $q_c$ with $F'(q_c)\ne0$. Then (fold lemma, `paper/sections/02_sources.tex` Lemma fold) $A=g+h\sqrt{1-x/x_+}$ near $x_+$ with $g,h$ analytic, and singularity analysis gives
$$a_n\sim-\frac{h(x_+)}{2\sqrt\pi}\,x_+^{-n}n^{-3/2},\qquad h(x_+)=F'(q_c)\sqrt{\frac{2x_+}{-x''(q_c)}} .$$
With $D=q\,d/dq$ (so $F'(q_c)=DF(\tau_c)/q_c$ and $x''(q_c)=D^2x(\tau_c)/q_c^2$ at the fold):
$$\boxed{\ K:=\lim_n|a_n|\,x_+^{\,n}n^{3/2}=\frac{|DF(\tau_c)|}{2\sqrt\pi}\sqrt{\frac{2x_+}{-D^2x(\tau_c)}},\qquad K^2=\frac{x_+\,DF(\tau_c)^2}{2\pi\,(-D^2x(\tau_c))}\ }\tag{1.1}$$
(For a companion-type sequence $U=RF/S$ with $R,S$ rational functions of $x$, replace $F$ by $G=R(x)F/S(x)$; at the fold $x'=0$ so $DG(\tau_c)=R(x_+)DF(\tau_c)/S(x_+)$.)

## 2. The derivative at the Fricke point

**Lemma 2.1 [proved].** Let $f$ be a meromorphic modular form of weight $k$ on $\Gamma_0(N)$ with $f|_kW_N=\varepsilon f$, $\varepsilon=\pm1$, analytic at $\tau_c=i/\sqrt N$. Then
$$i^k\varepsilon=+1\ \Longrightarrow\ Df(\tau_c)=\frac{k\sqrt N}{4\pi}\,f(\tau_c);\qquad i^k\varepsilon=-1\ \Longrightarrow\ f(\tau_c)=0 .$$
*Proof.* $f|_kW_N=\varepsilon f$ reads $f(-1/(N\tau))=\varepsilon(\sqrt N\tau)^kf(\tau)$. Differentiate in $\tau$ and put $\tau=\tau_c$ (a fixed point, $\sqrt N\tau_c=i$, $1/(N\tau_c^2)=-1$): $-f'(\tau_c)=\varepsilon\bigl[k\sqrt N\,i^{k-1}f(\tau_c)+i^kf'(\tau_c)\bigr]$. If $i^k\varepsilon=1$ this is $2f'(\tau_c)=-k\sqrt N\,i^{k-1}\varepsilon f(\tau_c)=ik\sqrt N f(\tau_c)$, and $Df=f'/(2\pi i)$. If $i^k\varepsilon=-1$ it is $k\sqrt N\,i^{k-1}\varepsilon f(\tau_c)=0$. $\square$

(For weight $2$ and $\varepsilon=-1$, or weight $4$ and $\varepsilon=+1$: the value is free and the derivative is elementary times it. For weight $4$ and $\varepsilon=-1$: the form vanishes at the Fricke point — this is the level-12 $\zeta(5)$ companion $E$, and also the mechanism of the Domb apparatus in `paper/main.tex` Thm oldspace.)

## 3. The dichotomy

**Theorem 3.1 [proved modulo the fold hypotheses of §1].** (a) *Elementary case.* If $i^k\varepsilon_F=+1$ then
$$K=\frac{k\sqrt N}{8\pi^{3/2}}\;\bigl|F(\tau_c)\bigr|\sqrt{\frac{2x_+}{-D^2x(\tau_c)}} .$$
If moreover $x=r(h)$ is a rational function of a weight-$0$ Hauptmodul $h$ of $\Gamma_0(N)$ (so that $D^2x(\tau_c)=r''(h_c)\,(Dh(\tau_c))^2$ at the fold $r'(h_c)=0$, $h_c=h(\tau_c)$) and $k=2$, then $F/Dh$ is a modular function and
$$K=\frac{\sqrt N}{4\pi^{3/2}}\Bigl|\frac{F}{Dh}(\tau_c)\Bigr|\sqrt{\frac{2x_+}{-r''(h_c)}}\ \in\ \overline{\mathbf Q}\cdot\pi^{-3/2},$$
with all three factors $(F/Dh)(\tau_c)$, $x_+$, $r''(h_c)$ algebraic (values of modular functions and of rational functions at a CM point).
(b) *CM case.* If $i^k\varepsilon_F=-1$ then $F(\tau_c)=0$, $DF(\tau_c)=\vartheta_kF(\tau_c)$ is a genuine weight-$(k+2)$ CM value (no $E_2$ anomaly), $D^2x(\tau_c)$ a weight-$4$ CM value, and by Chowla–Selberg
$$K\in\overline{\mathbf Q}\cdot\frac{\Omega_{\tau_c}^{\,k}}{\sqrt\pi},$$
where $\Omega_{\tau_c}^2$ is the weight-$2$ period of the CM point ($\Gamma(1/3)^6/\pi^4$ for $\mathbf Q(\sqrt{-3})$, $\Gamma(1/4)^4/\pi^3$ for $\mathbf Q(i)$, …). For the level-12 $\zeta(5)$ companion ($k=4$, $\tau_c=i/\sqrt{12}$): $K\in\overline{\mathbf Q}\,\Gamma(1/3)^{12}\pi^{-17/2}$, as found.

*Proof.* (1.1) with Lemma 2.1; in (a) the anomaly of $DF$ is absorbed into the elementary factor $k\sqrt N/(4\pi)$, and $Dh$ is a genuine weight-$2$ meromorphic form because $\sum_dr_d=0$ for a weight-$0$ eta quotient; in (b) $F(\tau_c)=0$ removes the anomaly term $\frac k{12}E_2F$ from $DF=\vartheta_kF+\frac k{12}E_2F$. $\square$

**Remark 3.2.** The classical asymptotic of Apéry's numbers, $a_n\sim(1+\sqrt2)^{4n+2}/(2^{9/4}\pi^{3/2}n^{3/2})$, is case (a) at $N=6$; the $\pi^{-3/2}$ is the $1/(4\pi)$ of Lemma 2.1 times the $1/(2\sqrt\pi)$ of singularity analysis, and $(1+\sqrt2)^2/2^{9/4}$ is $\frac{\sqrt6}{4}\,|F/Dh|(\tau_c)\sqrt{2x_+/(-r''(h_c))}$ — to be confirmed by the census. The same fold lemma gives the Apéry limit (Theorem B$^*$: $\xi=\Theta(q_c)+F\Theta'/F'$) and the asymptotic constant: both are values at the Fricke point, the limit of the Eichler integral and the constant of the form's derivative.

**Remark 3.3 (what decides the case).** For the second-order sporadic rows the fold is a cusp (log singularity) and the constant is elementary for a different reason (cusp width and the constant term of $F|S$). For the third-order rows $F$ has weight $2$; a fold at the Fricke point with $F(\tau_c)\ne0$ forces $\varepsilon_F=-1$ (Fricke-odd), case (a). CM periods can only enter through companions whose Fricke sign is "wrong" for their weight: weight $4$ anti-Fricke ($\zeta(5)$ constructions, Domb's $f^*$), weight $2$ Fricke-even, weight $6$ Fricke-even, etc.

## 4. Census

*(to be filled from the agent report: rows α, γ, ε, ζ, Cooper's s₇, s₁₀, s₁₈; the level-12 and level-16 ζ(5) systems; the level-12 ζ(7) parent.)*
