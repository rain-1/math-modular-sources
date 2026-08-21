# lattice/adelic_pade — scripts for consolidation/ADELIC_PADE.md

Every script is run by prepending the shared row builders, e.g.

    cat lattice/catalan_positivity/rows_common.gp lattice/catalan_explicit/moments.gp \
        lattice/adelic_pade/g2.gp lattice/adelic_pade/lawcheck.gp > run.gp
    gp -q run.gp

* `g2.gp`      — G_2 (2-adic Catalan period, repo normalisation = 8*zeta_2(2)) from the
                 Zudilin row: G2rat(M) returns [P_M/Q_M, 8M-1-4 s_2(M)] = value and precision.
                 Also `v2()`.  Verifies v_2(G_2 - P_m/Q_m) = 8m-1-4 s_2(m) for m<=20.
* `tj.gp`      — tau_j = v_2(A_{n,j} G_2 - B_{n,j}) for the whole moment family, and the
                 exact ultrametric identity for the mixed minor h_{j0,j1}.  n<=4, ~90 s.
* `lawcheck.gp`— verifies the three closed forms (v2A, v2B, tau) for all (n,j), n<=5. ~110 s.
* `sieve.gp`   — 2-adic margin v_2(G_2 - p/q) of the short vectors of the two-row lattice
                 K_n at k=22.4, over a ball of 2600 lattice points.  n<=16, ~100 s.
* `eps.gp`     — cancellation excess eps(c) and the exact identity margin = qZ - eps. n<=20.
