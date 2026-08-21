# `lattice/sources_s18_zud/` — sources of Cooper $s_{18}$ and Zudilin's Catalan row

Companion scripts to `consolidation/SOURCES_S18_ZUDILIN.md`.
Run from the repository root with stdin closed, e.g.

    gp -q lattice/sources_s18_zud/s18_struct.gp < /dev/null

(all `read()`s use absolute paths; gp enters its interactive loop after a script,
hence `< /dev/null`).

| script | what it does |
|---|---|
| `s18_setup.gp` | shared: $x_{18},t_{18},F_{18},\Phi=F\thq t$ as $q$-series (`NQ` = precision) |
| `s18_order.gp` | $F=\sum A_nt^n$; shows $[t^n]F\thq^{-k}\Phi=B_n$ iff $k=3$ ⇒ $r=3$, $w=2$ |
| `s18_order2test.gp` | no order-2 annihilator (deg ≤10); prints the order-3 operator |
| `s18_struct.gp` | $Z=-F\sqrt{1-16t}$, $(x-3)/(x+3)=\sqrt{1-12t}$, $c(m)/m\in\Z$, $\sqrt{F}\notin\Z[[q]]$ |
| `s18_mf.gp`, `s18_mf2.gp` | `mfinit`/`mftobasis`: $W=-F\thq\log x\in M_4(\Gamma_0(18))$ and its exact Eisenstein + CM-cuspidal split |
| `s18_pole.gp` | numerical: $x_{18}((3+i)/6)=-3$, $q_0=-e^{-\pi/3}$ |
| `s18_padic.gp` | $\xi_3=\tfrac12\zeta_3(2)$ to $3^{219}$; refutes $\tfrac12\zeta_3(3)$ |
| `ode_fit.gp` | generic minimal-ODE fitter `odefit(y,order,degree)` |
| `zud_row.gp` | Zudilin's $(Q_m,P_m)$ and Beukers' $q_n(x)$ |
| `zud_check.gp` | $Q_n=q_n(-n+\tfrac12)$ exactly; $q_n(-n)=$ Apéry $\zeta(2)$; $v_2(Q_m)=-4m+2s_2(m)$ |
| `zud_ode.gp`, `zud_ode2.gp` | minimal ODE has order 4; symbol $(1-11T-T^2)(1+10T+T^2)$ |
| `zud_mirror.gp` | mirror map $t(q)$ — non-integral, no Hauptmodul |
| `zud_padic.gp` | $\xi_2=\zeta_2(2)$ to $2^{398}$ |
| `zud_beukers.gp` | $P_m$ vs $p_m(-m+\tfrac12)$; Beukers' fixed-$x$ convergents at $x=1/2$ |
| `zud_moving.gp` | the moving Padé quotient does **not** converge 2-adically ($v_2$ of differences $\equiv4$) |
