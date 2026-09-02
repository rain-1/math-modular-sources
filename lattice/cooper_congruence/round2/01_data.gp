\\ 01_data.gp -- c, c', beta, gamma for the three rows, n <= M.  Exact.
read("lib.gp");
M = 600;
{
for(k=1,3,
  my(cp, b, g, bad);
  cp = CPvec(k,M);
  b  = Bvec(k,cp);
  g  = vector(M);
  bad = 0;
  for(n=1,M, if(b[n]%(n^2)!=0, bad=n; break, g[n]=b[n]/n^2));
  print("row ", NAM[k], "  n^2 | beta(n) for n<=", M, ": ", if(bad==0,"YES",concat("FAIL at ",bad)));
  print("  gamma(1..14) = ", vector(14,i,g[i]));
  write(concat(concat("gamma_",NAM[k]),".txt"), g);
  write(concat(concat("beta_",NAM[k]),".txt"), b);
  write(concat(concat("cp_",NAM[k]),".txt"), cp);
);
}
quit;
