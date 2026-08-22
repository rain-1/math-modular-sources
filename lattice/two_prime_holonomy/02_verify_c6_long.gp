default(parisizemax, 8000000000);
\\ 02_verify_c6_long.gp -- regenerate the conductor-6 alpha=2 and alpha=1 rows out
\\ to b = 320, refit the minimal joint operator from b <= 192 only, and verify it
\\ exactly over Q on all indices up to 320 (i.e. ~128 fresh ones).
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/lib.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/02_fit_recur.gp");
NL = 320;
chk(nm, f, r, d) =
{ my(qv = vector(NL+1), pv = vector(NL+1), z, seqsS, seqsL, nl, ek, cs, den, gg);
  gettime();
  for(b = 0, NL, z = f(b); qv[b+1] = z[1]; pv[b+1] = z[2]);
  print("--- ", nm, "  generated b<=", NL, " in ", gettime(), "ms");
  seqsS = [vector(201, k, qv[k]), vector(201, k, pv[k])];
  nl = vector(192 - r, k, k);
  ek = exactker(seqsS, r, d, nl, 40);
  cs = tocoefs(ek[1], r, d);
  den = lcm(vector(#ek[1], t, denominator(ek[1][t])));
  cs = vector(r+1, i, cs[i]*den);
  gg = 0; for(i = 1, r+1, gg = gcd(gg, content(cs[i])));
  cs = vector(r+1, i, cs[i]/gg);
  seqsL = [qv, pv];
  print("   fit used n <= ", 192-r, ";  total failures on n=1..", NL-r,
        " = ", verify(seqsL, cs, r, 1, NL - r));
  print("   failures on FRESH n = 200..", NL-r, " : ", verify(seqsL, cs, r, 200, NL - r));
};
chk("c6_a2", (b) -> chi6row(2*b, b), 2, 15);
chk("c6_a1", (b) -> chi6row(b, b),   2, 13);
quit;
