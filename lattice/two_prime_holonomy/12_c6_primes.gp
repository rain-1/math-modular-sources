default(parisizemax, 8000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/lib.gp");
{
print("=== conductor-6 family (MULTI_PRIME lib): prime support of den(Q_b), den(P_b) ===");
foreach([[1,4,20],[1,4,40],[1,4,60],[2,1,20],[2,1,40]], v,
  my(bb = v[3], aa = v[1]*bb/v[2], r = chi6row(aa,bb),
     dq = denominator(r[1]), dp = denominator(r[2]), f1, f2, n2, n3);
  f1 = if(dq==1, 1, my(f=factor(dq)); f[#f[,1],1]);
  f2 = if(dp==1, 1, my(f=factor(dp)); f[#f[,1],1]);
  n2 = valuation(dq,2); n3 = valuation(dq,3);
  printf("alpha=%d/%d b=%d: kap2=%.3f kap3=%.3f nu=%.3f  max prime den(Q)=%d den(P)=%d  (6b+5=%d)\n",
    v[1],v[2],bb, n2*1.0/bb, n3*1.0/bb,
    log((dq/2^n2/3^n3)*1.0)/bb, f1, f2, 6*bb+5));
}
quit;
