# Working note: the CM-trace obstruction (row s7)

Scripts `02_trace7.gp` … `07_obstr7.gp`; helper `heeg.gp` (Heegner forms on
$\Gamma_0(N)$, genus character, stabiliser orders) and `e2.gp` ($E_2^*$ by reduction
to the $\mathrm{SL}_2(\mathbf Z)$ fundamental domain).

Setting: $N=7$, $D_0=-3$, $\tau_0=(5+\sqrt{-3})/14$, $R=e^{\pi\sqrt3/7}=2.175682\ldots$,
$\kappa=1/(2\pi\operatorname{Im}\tau_0)=14/(2\pi\sqrt3)=1.286009\ldots$.
For $d=D_0m^2$ take the Heegner forms $[A,B,C]$, $7\mid A$, $B\equiv 5m\pmod{14}$
(the GKZ convention $\beta=\delta_0\sqrt d$ with $\delta_0^2\equiv-3\ (28)$, $\delta_0=5$),
the genus character $\chi_{-3}$, and the Fricke-anti-invariant Hauptmodul combination
$H=1/u-49u$.  Then
$$\operatorname{Tr}^{\mathrm{tw}}_{m}(H):=\sum_{Q}\chi_{-3}(Q)\,\frac{H(\alpha_Q)}{\omega_Q}
= i\sqrt3\;\tau_H(m),\qquad \tau_H(m)\in\mathbf Z .$$
$\tau_H(m)$, $m=1,\dots,20$ (the entries with $7\mid m$ are anomalous, the two
$\beta$-classes merging there):
$1,-6,5,18,-57,90,\ast,-458,1227,-1138,-2609,12570,-22061,\ast,104595,-283366,274837,596598,-2922599,5100094$.

**Finding.** $\beta(m)$ and $\tau_H(m)$ have the same sign pattern and
$\beta(m)/\tau_H(m)\to3$, but not exactly: numerically (`07_obstr7.log`, $m\le34$)
$$\frac{m\,\beta(m)}{(m-\kappa)\,\tau_H(m)}\ \longrightarrow\ 3
\qquad(3.000000000012\ \text{at}\ m=33,\ \ 2.999999967\ \text{at}\ m=27),$$
i.e. $\;a(m^2)=m\beta(m)=3(m-\kappa)\tau_H(m)+O(mR^{\theta m})$, $\theta<1$.
Since $\kappa\notin\mathbf Q$ (it is $7\sqrt3/\pi\cdot\frac13$), $\beta$ is **not** a rational
multiple of $\tau_H$; and an exact three-parameter fit
$\beta=\lambda_1\tau_H+\lambda_2\tau_J+\lambda_3\tau_1$ against the full basis
$\{1/u-49u,\;1/u+49u,\;1\}$ of functions with poles of order $\le1$ at the two cusps,
solved at $m=22,23,24$, misses $m=25,\dots,34$ by about $1\%$ — a relative error of
size $\kappa/m$, not of size $R^{-m/2}$.

**Obstruction (proof).** Every twisted CM trace of a modular function $G$ on $X_0(N)$
with poles only at the cusps satisfies
$\operatorname{Tr}^{\mathrm{tw}}_m(G)=\sum_{j\ge1}2\operatorname{Re}(c_jq_0^{-jm})+O(R^{\theta m})$
with **constant** $c_j$ (the principal-part coefficients of $G$), because the Heegner point
of largest height is $\alpha=-m\bar\tau_0$ (height $m\operatorname{Im}\tau_0$) together with
its Fricke partner, all other classes having height $\le\tfrac12 m\operatorname{Im}\tau_0$.
By `cooper_sources` §2.4 the polar data of $\Phi$ give
$\beta(m)=c'(m)+O(R^{\theta m})=-\nu^2\cdot2\operatorname{Re}\!\bigl(q_0^{-m}/g'(u_0)\bigr)\bigl(1-\kappa/m\bigr)+O(R^{\theta m})$.
The $1/m$ tail $\;\kappa\cdot 2\operatorname{Re}(\cdots)R^m/m\;$ is of size $R^m/m$, which
dominates $R^{\theta m}$ and is not an exponential polynomial with constant coefficients.
Hence $\beta$ is not a finite $\mathbf Q$-linear combination of twisted CM traces of
modular functions at the discriminants $D_0m^2$.  $\square$

The $1/m$ is exactly the Eichler-integral signature ($\Xi=D^{-1}\Phi$, and $\Phi$ has a
**double** pole).  A CM-trace identification of $\gamma$, if one exists, must therefore use a
second-order object (a higher Green's function $G_2$, or the CM values of the
non-holomorphic weight-2 completion), not the CM values of a modular function.

**A control that failed.** The natural real-analytic weight-0 invariant $E_2^*/F$
($F=D\log u\in M_2(\Gamma_0(7))$) does not repair this: numerically
$\sum_Q\chi(Q)\,\omega_Q^{-1}(E_2^*H/F)(\alpha_Q)=-3\operatorname{Tr}^{\mathrm{tw}}_m(H)+(\text{small})$,
so it produces no new sequence (`06_e2trace7.log`).
