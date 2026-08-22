default(parisizemax, 8000000000);
\\ 02_gen_halves.gp -- the two conductor-12 HALF rows c12half(a,b,1), c12half(a,b,5).
\\ c12rowB(a1,a5,b) = [160*s1*s5, s5*P1 + s1*P5], so the (Q,P) pair of the product
\\ form lies in the 4-dimensional shift-stable span of {s1 s5, s1 P5, P1 s5, P1 P5}.
read("/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/lib12.gp");
NMAX = 200;
outdir = "/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/out/";
dump(nm, f) =
{ my(qv = vector(NMAX+1), pv = vector(NMAX+1), z, fn);
  gettime();
  for(b = 0, NMAX, z = f(b); qv[b+1] = z[1]; pv[b+1] = z[2]);
  fn = concat(concat(outdir, nm), ".seq");
  write(fn, "QV = ", qv); write(fn, "PV = ", pv);
  print(nm, " done  time=", gettime(), "ms");
};
dump("h12_a2_r1", (b) -> c12half(2*b, b, 1));
dump("h12_a2_r5", (b) -> c12half(2*b, b, 5));
dump("h12_a1_r1", (b) -> c12half(b, b, 1));
dump("h12_a1_r5", (b) -> c12half(b, b, 5));
quit;
