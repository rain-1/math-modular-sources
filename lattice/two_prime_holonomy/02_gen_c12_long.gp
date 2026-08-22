default(parisizemax, 8000000000);
\\ 02_gen_c12_long.gp -- longer conductor-12 product-form sequences, for the
\\ direct order-4 joint recurrence fit.
read("/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/lib12.gp");
NMAX = 340;
outdir = "/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/out/";
dump(nm, f) =
{ my(qv = vector(NMAX+1), pv = vector(NMAX+1), z, fn);
  gettime();
  for(b = 0, NMAX, z = f(b); qv[b+1] = z[1]; pv[b+1] = z[2]);
  fn = concat(concat(outdir, nm), ".seq");
  write(fn, "QV = ", qv); write(fn, "PV = ", pv);
  print(nm, " done  time=", gettime(), "ms"); };
dump("c12_a2_long", (b) -> c12rowB(2*b, 2*b, b));
dump("c12_a1_long", (b) -> c12rowB(b, b, b));
quit;
