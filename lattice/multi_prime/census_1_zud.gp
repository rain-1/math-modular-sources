default(parisizemax, 8000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/lib.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/census_util.gp");
fn = "/home/ubuntu/code/math-modular-sources/lattice/multi_prime/out/census_1_zud.log";
MM = 300;
zz = zudrow(MM);
qv = zz[1]; pv = zz[2];

logit(fn, Str("### Zudilin Catalan row zudrow(", MM, ") ###"));

/* (a) is den(Q_m) a pure power of 2, with v_2 = 4m - 2 s_2(m)? */
s2(m) = { my(t = 0, x = m); while(x > 0, t += x % 2; x = x \ 2); t; };
{
  my(badpure = 0, badfor = 0, firstbad = -1);
  for(m = 0, MM,
    my(d = denominator(qv[m+1]), v = valuation(denominator(qv[m+1]), 2));
    if(d / 2^v != 1, badpure++; if(firstbad < 0, firstbad = m));
    if(v != 4*m - 2*s2(m), badfor++));
  logit(fn, Str("  # m<=", MM, " with den(Q_m) NOT a pure power of 2 : ", badpure,
                "  (first bad m=", firstbad, ")"));
  logit(fn, Str("  # m<=", MM, " with v_2(den Q_m) != 4m - 2 s_2(m)  : ", badfor));
}
/* same for P */
{
  my(badpure = 0);
  for(m = 0, MM, my(d = denominator(pv[m+1]), v = valuation(d, 2));
    if(d / 2^v != 1, badpure++));
  logit(fn, Str("  # m<=", MM, " with den(P_m) NOT a pure power of 2 : ", badpure));
}
{ logit(fn, Str("  v_2(den Q_m) at m=100,200,299: ", [valuation(denominator(qv[101]),2), valuation(denominator(qv[201]),2), valuation(denominator(qv[300]),2)])); }
logit(fn, "");

ns = [50, 100, 150, 200, 250, 299];
pl = primes(15);
logit(fn, Str("primes tested: ", pl));
report(fn, "zudrow", qv, pv, pl, ns);

/* fine increments of the 2-adic slope, last 20 n */
{
  my(s = "");
  for(m = MM-20, MM-1, s = Str(s, " ", valuation(pv[m+1]/qv[m+1] - pv[m]/qv[m], 2)));
  logit(fn, Str("  fine v_2(incr), m=", MM-20, "..", MM-1, ":", s));
}
{
  my(s = "");
  for(m = MM-20, MM-1, s = Str(s, " ", valuation(pv[m+1]/qv[m+1] - pv[m]/qv[m], 2)
                                  - valuation(pv[m]/qv[m] - pv[m-1]/qv[m-1], 2)));
  logit(fn, Str("  first differences of v_2(incr):", s));
}
logit(fn, "DONE");
quit;
