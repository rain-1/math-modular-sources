# The $\Gamma_1(5)$ theorem does not follow: closure of bucket 4, Track D

*Fable, 2026-09-02. Full report with scripts: `lattice/gamma15/REPORT.md` (800 lines, three sub-reports in `task2/`, `task3/`). This note records the statement, the verdict, the two structural reasons, and the by-products.*

## 1. The target

On Zagier's row **D** (Apéry's $\zeta(2)$ row on $\Gamma_1(5)$; $A_n=\sum_k\binom nk^2\binom{n+k}k$, $\lambda_1\lambda_2=-1$, $k=2$) the fold-regular weight-three Eisenstein space is two-dimensional over $K=\mathbf Q(\sqrt5)$ (`lattice/hostscan/REPORT.md` §10). Besides $\zeta(2)/5$ it carries
$$\xi=\varphi^5\,\mathrm{Im}\,L(2,\psi_4)-\mathrm{Re}\,L(2,\psi_4)=0.65563418884065676633\ldots,$$
$\psi_4$ the odd quartic character mod 5, source $\Phi_{\rm new}=(1+i\varphi^5)E_3^{\psi_4,\mathbf 1}+(1-i\varphi^5)E_3^{\bar\psi_4,\mathbf 1}\in\mathbf Z[\varphi][[q]]$. Since the host, $k$, and (as Task 3 verified) the whole denominator array are CDT's own, the question was whether CDT's construction, run over $K$ with the averaged number-field bound (`NUMBER_FIELD_HOLONOMY.md`, which is CDT's own Remark BCboundK), proves

> **(T)** $1,\ \zeta(2),\ \xi$ are linearly independent over $\mathbf Q(\sqrt5)$.

## 2. Verdict: no, by 1.44 nats at best, and (T) is not even the hypothesis it refutes

| | hypothesis refuted | $\bar L$ | $\overline{\mathrm{BC}}$ | entry $\bar L-\tau$ | bound | margin $14(\bar L-\tau)-\overline{\mathrm{BC}}$ |
|---|---|---|---|---|---|---|
| A | (T′), see §3.2 | 4.90502 | 10.81633 | $+0.670$ | $m\le16.15$ | $\mathbf{-1.44}$ |
| B | (T) itself | 4.59778 | 9.47631 | $+0.362$ | $m\le26.15$ | $\mathbf{-4.40}$ |

with $\tau=16603/3920$ unchanged, $m=14$ optimal (every other $m$ is worse; `INVENTORY_BOUND.md` R1 says 14 is the complete supply at $\max e_i\le1$), entry passing at both branches. The verdict is robust: CDT's convexity gain transported (factor 0.973) leaves $m\le15.7$; dropping every unproved exclusion changes the margin by $<0.02$ (§6(b′) of the report); the conformal radii of the whole slit-and-lune family saturate at $0.823$ / $0.407$, far from the $0.485$ the second place would need.

## 3. The two reasons, both at the second real place

**3.1 Geometry.** The pure module's puncture is $s=t_2=-\varphi^5$, a unit. At $v_1$ the removed fold $t_1=\varphi^{-5}$ is hyperbolically *deeper* than CDT's $-1/72$ (normalised $Y_1=-\varphi^{-15}/(5\sqrt5)=-6.6\cdot10^{-5}$; deepest bad $h$-preimage at $|z|=0.536$ against CDT's $0.402$), and $v_1$ beats CDT by $0.58$ nats. At $v_2$ the roles of $t_1,t_2$ swap: $Y_2=-\varphi^{15}/(5\sqrt5)=-122$, the fold sits at $30$ orbifold radii, the deepest bad preimage at $|z|=0.214$, and $v_2$ loses $3.51$ nats. Halved, $-1.46$. Since $Y_1Y_2=1/125$ is forced by $t_1t_2=-1$, no contour or inventory can repair this; only a different host could, and by the host sweep there is none.

**3.2 Hypothesis.** A single $K$-relation (i) $a+b\,\zeta(2)/5+c\,\xi=0$ removes the fold singularity of the conditional function $H$ at $v_1$ only. At $v_2$ the conjugate $\sigma H$ must be regular at the *far* cusp, where the periods are not the conjugates of the near-cusp ones: the far-cusp period of Zagier D's companion is $-\tfrac{11}{5}\zeta(2)$ (not $\zeta(2)/5$), while $\Phi_{\rm new}$'s is exactly $\sigma(\xi)=-\mathrm{Re}\,L-\varphi^{-5}\mathrm{Im}\,L$ (monodromy around $t_2$, 174 digits; Abel summation cross-check). So the second place needs
$$\text{(ii)}\qquad \sigma(a)-11\,\sigma(b)\,\zeta(2)/5+\sigma(c)\,\sigma(\xi)=0,$$
which (i) does not imply, and which is never proportional to (i) when $c\ne0$. The construction therefore refutes **(T′)**: no $(a,b,c)\in K^3\setminus0$ satisfies (i) and (ii) simultaneously. That is branch A. Branch B keeps (T) honest by making the $v_2$ contour avoid the principal preimage as well, and pays three more nats.

*Why $-11$.* It is the orientation, not Galois. The outer direction is $D^{-2}\Phi_D=\sum_d w(d)\mathrm{Li}_2(q^d)$ and $\mathrm{Li}_2$ takes the value $\zeta(2)$ at $q\to+1$ but $-\zeta(2)/2$ on odd $d$ at $q\to-1$, turning the character sum $\sum_j j\,w(j)=-1$ (factor $1/5$) into $22$ (factor $-11/5$). The inner direction is a Lambert series whose cusp weight is $-\tfrac12$ uniformly, so on the $\xi$-line equivariance is exact. Note $11=\mathrm{Tr}_{K/\mathbf Q}\varphi^5=-(t_1+t_2)$; the coincidence is not claimed to be structural.

## 4. By-products worth keeping

1. **Two corrections to CDT's Appendix A** (their arXiv:2408.15403v2, `papers/cdt/cdt2/L2chi.tex`). (i) The published slit parameters $\theta_i$ do not satisfy their Lemma A.4.4: $-1/72$ has six, not four, non-principal $h$-preimages inside $|z|<77/100$ (the pair at $|z|=0.753925$ was missed; argument principle gives 7 solutions in the disc), and with the printed $\theta_i$ five of the six lie inside their $\Omega$. (ii) The Bost–Charles integral of the published $\varphi$ is $12.12$, not $11.845$. **The theorem is unaffected**: a successive-slitting design with their own $R=77/100$, $c=15/2$, shrink $995/1000$ and six slits gives $|\psi'(0)|=0.62934$ (theirs $0.62922$), $\mathrm{BC}=11.8449$ (theirs $11.845$), bound $13.9906$ (theirs $13.9938$); only the transcription of the slit placement is wrong. Optimising $R$ in the same family gives margin $+0.021$ against their $+0.005$.
2. **A Jensen reduction of the Bost–Charles integral.** For $\varphi=h\circ\Phi$ with $\Phi$ univalent, $\mathrm{BC}(\varphi)=\log|\varphi'(0)|+\int_{\mathbf T}\sum_{w\in\mathbf D,\ \varphi(w)=\varphi(z)}\log\frac1{|w|}\,d\mu(z)$: the excess over $\log|\varphi'(0)|$ is a pure multivalency term, zero iff $\varphi$ is univalent (recovering CDT's univalent case). Proved; the interior solutions are the $\Gamma_0(2)$-orbit points of $\Phi(z)$ lying in $\Omega$.
3. **A weight-three identity.** $\mathrm{Re}\,L(3,\psi_4)=\varphi^5\,\mathrm{Im}\,L(3,\psi_4)$ (212 digits in the report; re-checked independently at 60 digits with PARI's `lfun`). Since $s=3$ is critical for an odd character, both sides are $\pi^3$ times algebraic numbers and the identity is in principle a Gauss-sum/Bernoulli computation; it is the $s=3$ shadow of the fold-regularity ratio $-\varphi^5$ at $s=2$, which is the non-critical, non-trivial one.
4. **The monodromy criterion for bad preimages.** Only those $h$-preimages of the fold at which the continued conditional function acquires its logarithm need be excluded; for a deck element $\gamma$ this is the cocycle condition $\rho(\gamma)H-H\in\langle u_{t_1}\rangle$. For the decisive pair (bottom rows $(2,\pm1)$, the parabolic generator at the cusp $y=\infty$, whose loop encircles $t_2$) it fails because the fixed lines $\langle u_{t_1}\rangle\ne\langle u_{t_2}\rangle$ (distinct cusps; numerically $|u_{t_1}\wedge u_{t_2}|=0.318$ on $\Gamma_1(5)$, $0.737$ on CDT's host). Had it held, the margin would be $+0.51$ here and $+0.65$ on CDT's own host. Shallower preimages are moot (§6(b′)).

## 5. What is left

* A host whose **outer** direction has uniform cusp weights (Lambert rather than dilogarithm orientation) would remove the hypothesis tax and put branch A back on (T). That is a concrete search criterion, orthogonal to the geometry, but the geometry ($Y_2=-\varphi^{15}/5\sqrt5$) already loses 1.46 nats, and the host sweep says level 5 and level 6 are the only Apéry-perfect $k=2$ hosts.
* Assumed, not transported: unbounded-degree $K(y)$-independence of the fourteen (CDT's Lemma 12.1.1; degree $\le5$ certified at a split and an inert prime); the dictionary "Apéry limit = far-cusp period = log coefficient" at the far cusp (verified to 212 digits at the near cusp).
* The Appendix-A corrections should be communicated to the authors.

## 6. Bucket 4 as a whole

| track | question | outcome |
|---|---|---|
| A | can a larger lcm-free inventory rescue Catalan / $X_1(5)\,\mathrm{Sym}^2$? | no: CDT's 14 are the complete supply, best margins $-8.04$ / $-6.99$ (`INVENTORY_BOUND.md` §4) |
| B | is there another Apéry-perfect host, $N\le120$? | no; but a new period $\xi$ on $\Gamma_1(5)$ (`lattice/hostscan/REPORT.md` §10) |
| C | do split-denominator two-variable companions exist? | no: tail obstruction and Landau ceiling (`TWO_VARIABLE_HOLONOMY.md` §6) |
| D | does CDT's construction over $\mathbf Q(\sqrt5)$ prove (T)? | no: $-1.44$ nats on the doubled hypothesis, $-4.40$ on (T); second real place responsible twice |

Every lever the finder identified has now been pulled and measured. The arithmetic-holonomy architecture at weight $\le3$ proves exactly what CDT proved, and the reasons it proves nothing else are quantitative and recorded.
