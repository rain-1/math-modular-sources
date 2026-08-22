# lattice/catalan_three_period

Scripts for `consolidation/CATALAN_THREE_PERIOD.md`.  All GP scripts use absolute
paths, `default(parisizemax,...)`, and end with `quit`.

| script | what it does | log |
|---|---|---|
| `10_periods.gp` | Eichler period polynomials $(\alpha,\beta)$ of the six weight-3 $\chi_{-4}$ classes at all five non-$\infty$ cusps of $\Gamma_0(16)$, from the numerical transformation law at the symmetric point of each horocycle. Two-point consistency of the affine fit is $\sim10^{-133}$ at every cusp and class. | `10run.log` |
| `11_a0.gp` | three-parameter fit $(\alpha,\beta,a_0)$; recovers $a_0=0$ (inner) and $a_0=-1/4$ (outer) to 130 digits, confirming the constant-term correction term | (stdout) |
| `12_exact.gp` | recognises $(\alpha,\beta,a_0)$ at the fold cusp $1/2$ exactly via `bestappr`; prints the two linear functionals and $\xi=-\beta/16$ | `12run.log` |
| `13_outer_dbl.gp` | numerical test of the outer prediction on the actual power series: $b_n/a_n$, $\|c_n\|^{1/n}$, $v_2(b_n)$, $d_n^k$ denominators, $N=220$ | `13run.log` |
| `14_symm.gp` | normaliser descent $y=4x^2/(4x+1)$; 2-adic slope and radius in $y$ of $A$, $\mathrm{Li}_j(-4x)$, and the two doubly-small companions | `14run.log` |
| `15_rank.gp` | $\mathbf Q(y)$-independence (rank mod $2^{61}-1$) of the three small generators and of their CDT 7-member $\theta$-orbits | `15run.log` |
| `17_conditions.gp` | reduces the exact $(\alpha,\beta)$ at each of the five non-$\infty$ cusps to the linear condition on $(c_1,c_2,c_4)$: $(64,8,1)$, $(0,8,1)$, $(0,0,1)$, $(0,0,1)$, vacuous — identifying the cusp $1/2$ as the fold | `17run.log` |
| `16_margin.py` | first-pass adelic margins with the corrected *inventory* but the previously **assumed** slope profile; reuses `lattice/adelic_holonomy/adelic_bound.py` | `16run.log` |
| `18_slopes.gp` | 2-adic slope in $y$ of every object in the inventory, including the host $A$ and the raw companions — this is where the conditional slope $-2$ is found | `18run.log` |
| `19_margin3.py` | **definitive** margins, both slope profiles (assumed / measured), with a scale-covariance check | `19run.log` |

Reads `lattice/catalan_two_classes/00_setup.gp` for the hosts and companions.
