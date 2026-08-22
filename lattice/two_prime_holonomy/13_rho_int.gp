default(parisizemax, 8000000000);
default(realprecision,60);
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/lib.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/lib12.gp");
\\ rho_int = (1/b) log | Q_b * den(Q_b) |  -- the growth rate of the integral normalisation.
\\ Both log Lambda and nu drift (log b); their sum should not.
ri(q, bb) = my(dn = denominator(q)); log(abs(q*1.0)*dn)/bb;
{
print("=== rho_int stability:  (1/b) log |Q_b * den(Q_b)| ===");
print("family/alpha       b=20      b=40      b=60      b=80");
foreach([[1,4],[1,2],[1,1],[2,1]], fr,
  my(v = vector(4));
  for(k=1,4, my(bb = 20*k, aa); if(bb % fr[2], v[k]=0.0, aa = fr[1]*bb/fr[2];
     v[k] = ri(chi6row(aa,bb)[1], bb)));
  printf("c6  alpha=%d/%d  %9.4f %9.4f %9.4f %9.4f\n", fr[1],fr[2],v[1],v[2],v[3],v[4]));
foreach([[1,4],[1,2],[1,1],[2,1]], fr,
  my(v = vector(4));
  for(k=1,4, my(bb = 20*k, aa); if(bb % fr[2], v[k]=0.0, aa = fr[1]*bb/fr[2];
     v[k] = ri(c12rowB(aa,aa,bb)[1], bb)));
  printf("c12 alpha=%d/%d  %9.4f %9.4f %9.4f %9.4f\n", fr[1],fr[2],v[1],v[2],v[3],v[4]));
}
quit;
