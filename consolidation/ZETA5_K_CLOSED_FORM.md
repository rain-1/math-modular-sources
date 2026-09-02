# The level-12 $\zeta(5)$ asymptotic constant $K$: closed form

*Fable, 2026-09-02. Closes the "open piece" of `CLAUDE_FINDINGS.md` (2026-08-21): the constant $K$ in $c_n\sim-K\,(7+4\sqrt3)^n n^{-3/2}$ for the level-12 CM-isogeny $\zeta(5)$ system (`packages/phase 1 zeta math package/modular_apery_cm_isogeny_research_notes.txt`, §II.G–L). Scripts: `lattice/zeta5_K/`.*

**Setting.** $h_{12}=\eta_1^3\eta_4\eta_6^2/(\eta_2^2\eta_3\eta_{12}^3)=q^{-1}-3+\dots$, $x=h_{12}/((h_{12}+3)(h_{12}+4))=q-4q^2+\dots$ (Fricke-invariant, branch values $x_\pm=7\mp4\sqrt3$ at $h_{12}=\pm2\sqrt3$, the Fricke point $\tau_*=i/\sqrt{12}$, $q_c=e^{-\pi/\sqrt3}$), $E=(-9E_4(3\tau)+16E_4(4\tau))/7$, $Q(x)=143x^3+189x^2+21x+7$, $U=7E/Q(x)=\sum c_nx^n$ with $c_n\in\mathbf Z$ ($1,-3,-18,-194,\dots$, regenerated from the definition and checked against the 14-term recurrence of §II.J).

**Value [verified].** From $c_n$ to $n=4200$ (exact recurrence) and Richardson extrapolation in $1/n$ (16 points, stable to $47$ digits):
$$K=0.68537184849053440595143004488054241747925748918627\ldots$$

**Fold formula [proved, verified $10^{-51}$].** $x$ has a simple fold at $q_c$ ($x'(q_c)=0$, verified), $U\circ x=G(q)=7E/Q(x)$ is analytic at $q_c$, and singularity analysis (as in the fold lemma of `paper/main.tex` §2) gives $c_n\sim-\frac{h(x_+)}{2\sqrt\pi}x_+^{-n}n^{-3/2}$ with $h(x_+)=G'(q_c)\sqrt{2x_+/(-x''(q_c))}$, $G'(q_c)=7E'(q_c)/Q(x_+)$ (since $x'(q_c)=0$). In terms of $D=q\,d/dq$ and $D^2x(\tau_*)=q_c^2x''(q_c)$:
$$\boxed{\ K^2=\frac{49\,x_+\,(DE)(\tau_*)^2}{2\pi\,Q(x_+)^2\,\bigl(-D^2x(\tau_*)\bigr)}\ }$$
Numerically $K^2=0.46973457070333204886965161962712295965\ldots$ on both sides.

**Two structural facts.** (i) $E(\tau_*)=0$ (to $10^{-77}$): the anti-Fricke weight-4 companion vanishes at the Fricke point, so $DE(\tau_*)=\vartheta E(\tau_*)$ is a genuine weight-6 CM value with no $E_2$ anomaly. (ii) $D^2x(\tau_*)=x''(h_c)\,(Dh_{12})(\tau_*)^2$ with $h_c=2\sqrt3$, $x''(h_c)=-2h_c/((h_c+3)(h_c+4))^2$, and $Dh_{12}=h_{12}\cdot\frac1{24}\sum_dr_d\,d\,E_2(d\tau)$ is a genuine weight-2 meromorphic form ($\sum r_d=0$). Hence $K^2\in\overline{\mathbf Q}\cdot\Omega^{8}/\pi$ with $\Omega^2:=\Gamma(1/3)^6/\pi^4$ the weight-2 Chowla–Selberg unit of $\mathbf Q(\sqrt{-3})$.

**CM values [identified, algdep at 70–130 digits, then verified].** With $W^2=\Gamma(1/3)^6/\pi^4$:
$$\frac{DE(\tau_*)}{W^6}=u=\frac{243\,(23-15\sqrt3)}{14336}\quad(2^{21}\!\cdot\!49\,u^2-7\!\cdot\!2^{11}\!\cdot\!3^5\!\cdot\!23\,u-3^{10}\!\cdot\!73=0),\qquad \Bigl(\frac{Dh_{12}(\tau_*)}{W^2}\Bigr)^3=a^3=-\frac{81\,(45+26\sqrt3)}{2048}\quad(2^{22}a^6+2^{12}3^65\,a^3-3^9=0).$$
Also $Q(x_+)=96(2205-1273\sqrt3)$, $x_+/(-x''(h_c))=12+7\sqrt3$.

**Closed form [verified to $51$ digits].**
$$\boxed{\ K=\frac{\Gamma(1/3)^{12}}{\pi^{17/2}}\cdot\frac{3^{8/3}\,(15\sqrt3-23)\,\sqrt{(12+7\sqrt3)/2}}{2^{37/3}\,(2205-1273\sqrt3)\,(45+26\sqrt3)^{1/3}}\ }$$
i.e. $K=t\,\Gamma(1/3)^{12}\pi^{-17/2}$ with $t=0.0843620623288518148561713194046687969100356126578958\ldots$, an algebraic number of degree $>12$ (it involves $\sqrt3$, a square root of $(12+7\sqrt3)/2$ and the cube root of $45+26\sqrt3$). Consequently the full linear-form asymptotic of §II.L is now closed:
$$\frac{d_n}{c_n}-\frac{11}{144}\zeta(5)\ \sim\ \frac{13}{96K}(-1)^n(7-4\sqrt3)^n\,n^{1/2}(\log n)^4,\qquad \frac{13}{96K}=0.19758130854792016\ldots$$

**Method note.** The identification needed the right transcendental unit: for $\mathbf Q(\sqrt{-3})$ the weight-$k$ CM values are algebraic multiples of $(\Gamma(1/3)^6/\pi^4)^{k/2}$ (not $\Gamma(1/3)^6/\pi^2$); with the wrong unit `algdep` returns noise. The conductor-$4$ point $\tau_*=i/\sqrt{12}$ (discriminant $-48$) contributes the cube root.
