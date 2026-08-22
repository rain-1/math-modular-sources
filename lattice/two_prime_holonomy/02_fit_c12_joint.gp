default(parisizemax, 8000000000);
\\ 02_fit_c12_joint.gp -- minimal joint (Q,P) operator for the conductor-12
\\ product-form family c12rowB, on the long sequences (b <= 340).
\\ alpha=2 was already settled at (r,d) = (3,60) from the b<=200 data:
\\   chi(x) = (x - 1048576)(x + 1024)^2.
\\ Here we search harder for alpha=1.
read("/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/02_fit_recur.gp");
dsweep(nm, r, dlo, dhi) =
{ my(fn = concat(concat(outdir, nm), ".seq"), NA, nl);
  read(fn); NA = #QV - 1;
  nl = vector(NA - 10 - r, k, k);
  print("--- ", nm, " r=", r, "  rows=", 2*#nl);
  for(d = dlo, dhi,
    if((r+1)*(d+1) > 2*#nl - 20, print("   out of data at d=", d); break);
    my(nu = #matker(buildmat([QV,PV], r, d, nl, nextprime(2^61))));
    if(nu > 0, print("   HIT d=", d, " nullity=", nu); break));
  print("   (sweep finished)");
};
dsweep("c12_a2_long", 3, 58, 62);
dsweep("c12_a1_long", 3, 60, 150);
dsweep("c12_a1_long", 4, 60, 120);
quit;
