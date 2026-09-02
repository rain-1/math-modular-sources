\\ 20_extend.gp -- extend c', beta, gamma to n <= M (exact) and write data files.
read("lib.gp");
M = 12000;
gettime();
{
for(k=1,3,
  my(cp, b, g, bad);
  cp = CPvec(k,M);
  b  = Bvec(k,cp);
  g  = vector(M);
  bad = 0;
  for(n=1,M, if(b[n]%(n^2)!=0, bad=n; break, g[n]=b[n]/n^2));
  print("row ", NAM[k], "  n^2 | beta(n) for n<=", M, ": ", if(bad==0,"YES",concat("FAIL at ",bad)), "   [", gettime(), " ms]");
  write(concat(concat("20c_gamma_",NAM[k]),".txt"), g);
  write(concat(concat("20c_beta_",NAM[k]),".txt"), b);
  write(concat(concat("20c_cp_",NAM[k]),".txt"), cp);
);
}
quit;
