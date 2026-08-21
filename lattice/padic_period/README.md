# lattice/padic_period — scripts for consolidation/PADIC_PERIOD.md

Run with `timeout 110 gp -q <script>` from the repo root (they `read("lattice/rigidity/lib.gp")`).

- `involution.gp`     gamma = [8,-3;12,-4] = W_4 on X_0(12): t_C(gamma tau)=ubar,
                      t_C(2 gamma tau)=vbar, t_F invariant, F_F and Phi_F have
                      W_4-eigenvalue +1 (37 digits)
- `arch_period.gp`    the archimedean W-period identity: P_gamma(tau) =
                      (5/4) xi_inf (6tau-3) - i pi^2/(6 sqrt3) (6tau-1), 7 points, 37 digits
- `padic_point.gp`    3-adic pointwise test of the W-functional equation at three
                      rational points, precision 3^700; and recovery of xi_3 from it
- `padic_point_hi.gp` same at precision 3^1600 (N=1800, K=900); ~0.4 s
