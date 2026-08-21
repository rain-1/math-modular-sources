Scripts for consolidation/CATALAN_DIRECTIONAL.md (evaluation of Sol's directional proposal).
gp does not pick up read("rows_common.gp") reliably here; run each as

  cat rows_common.gp directional.gp > /tmp/run.gp; echo '\q' >> /tmp/run.gp; gp -q /tmp/run.gp

  rows_common.gp  exact Zudilin (index 3n) and Nesterenko (4,7) integer rows in the
                  Lean/paper normalisation (from lattice/catalan_audit/{rows.gp,nest2.gp})
  directional.gp  per-n directional mass log(Phi_n)/n, rigorous vs optimistic, net of
                  lattice-forced divisibility; plus the beta_Z/beta_N/F/H exponent table
  absorb.gp       the absorption check (v_l(Y)=0, v_l(U)<=v_l(S)) and the window-modulus variant
  detail.gp       per-prime valuation table at n=20
  ceiling.gp      absolute ceiling: odd gcd(V,U) over odd gcd(X,Y), all primes incl. l>6n
