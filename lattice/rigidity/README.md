# lattice/rigidity — scripts for consolidation/RIGIDITY_PROOF.md

Run with `timeout 560 gp -q <script>` (PARI/GP). `lib.gp` and `setup.gp` are read
by the others; `cover.gp` sets `PS` (q-precision) before reading `setup.gp`.

- `setup.gp`      eta quotients t_X, F_X and sources Phi_X = F_X D t_X for X = B, C, F
- `cover_checks.gp`  source identities and the C/F cover identities (all EXACT to O(q^220))
- `psi0_test.gp`  3-adic limit of the (1-4V_2)S companion on the C-curve  (= 0, slope 2)
- `galois_test.gp`  branch-involution invariance of the transfer identity, 3-adically
- `bcase_test.gp` Theta_B decomposition; C-curve rows = Zagier rows; failure of naive B expansion
- `btest.gp`      t_B is not in Q(t_C, t_C(2tau)); cubic for t_B over Q(t_F)
- `padic_limits.gp`  xi_3 for rows B, C, F to N = 400
- `radius.gp`     kappa_3 = 0 evidence (max v_3(a_n) for n <= 400)
