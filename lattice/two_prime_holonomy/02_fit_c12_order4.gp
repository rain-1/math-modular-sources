default(parisizemax, 8000000000);
\\ 02_fit_c12_order4.gp -- direct joint fit for the conductor-12 product-form
\\ family c12rowB (longer sequences, b <= 340).  The minimal joint operator turns
\\ out to have ORDER 3 (not 4): Q and P both lie in Sym^2 of the common
\\ 2-dimensional half-row solution space, not in the full tensor product.
\\ Predicted from the halves:
\\   alpha=2 : chi(x) = (x - 1048576)(x + 1024)(x - 1)
\\   alpha=1 : balance polynomial with roots {(17+12 sqrt2)^2, 1, (17-12 sqrt2)^2}
\\             = {577+408 sqrt2, 1, 577-408 sqrt2}, scaling u_(n+1)/u_n ~ rho n^-4.
read("/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/02_fit_recur.gp");
qsweep(nm, rmax, dmax) =
{ my(fn = concat(concat(outdir, nm), ".seq"), NA, NF, nl);
  read(fn); NA = #QV - 1; NF = NA - 10;
  print("--- sweep ", nm, " joint  NALL=", NA);
  nl = vector(NF - rmax, k, k);
  sweep([QV, PV], rmax, dmax, nl, nextprime(2^61));
};
qsweep("c12_a2_long", 3, 75);
qsweep("c12_a1_long", 3, 75);
run("c12_a2_long", 3, 75, 0);
run("c12_a1_long", 3, 75, 0);
quit;
