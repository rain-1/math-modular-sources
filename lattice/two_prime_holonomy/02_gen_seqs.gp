default(parisizemax, 8000000000);
\\ 02_gen_seqs.gp -- generate exact (Q_b,P_b) sequences for the conductor-6 and
\\ conductor-12 well-poised families, alpha = 2 and alpha = 1.
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/lib.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/lib12.gp");

NMAX = 200;
outdir = "/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/out/";

dump(nm, fq, fp) =
{ my(qv = vector(NMAX+1), pv = vector(NMAX+1), z, fn);
  gettime();
  for(b = 0, NMAX,
    z = fq(b);
    qv[b+1] = z[1]; pv[b+1] = z[2]);
  fn = concat(concat(outdir, nm), ".seq");
  write(fn, "QV = ", qv);
  write(fn, "PV = ", pv);
  print(nm, "  done, NMAX=", NMAX, "  time=", gettime(), "ms  ",
        "size(Q_last)=", #digits(numerator(qv[NMAX+1])), " digits");
};

dump("c6_a2",  (b) -> chi6row(2*b, b), 0);
dump("c6_a1",  (b) -> chi6row(b, b), 0);
dump("c12_a2", (b) -> c12rowB(2*b, 2*b, b), 0);
dump("c12_a1", (b) -> c12rowB(b, b, b), 0);
quit;
