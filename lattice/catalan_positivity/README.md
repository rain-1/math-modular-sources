Scripts for consolidation/CATALAN_POSITIVITY.md.
Run each as:  cat rows_common.gp SCRIPT.gp > run.gp; gp -q run.gp   (each <15 s)
  rows_common.gp  exact Zudilin(3n)/Nesterenko(4,7) integer rows (copied from lattice/catalan_directional)
  diag.gp         2-adic valuations of the rows; index of K_n vs M_n
  integrals.gp    Mechanism A: verifies K_N(n)/K_Z(3n) = (xy/(1-xy))^n
  cone_honest.gp  honest lattice {c : q,p both integral}: index, first minimum, positive-cone minimum
  cone_sweep.gp   main table, n <= 44
  control.gp      same, with G replaced by the rational G* = bestappr(G,10^320)
  character.gp    Mechanism C: sign vs congruence class of (c_Z,c_N)
