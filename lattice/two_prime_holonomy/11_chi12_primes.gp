default(parisizemax, 8000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/lib12.gp");
{
print("=== product form: largest prime in den(Q_b), den(P_b) vs 12b+11 ===");
foreach([[1,2,20],[1,1,20],[1,2,40],[1,1,40],[1,2,60],[1,1,60]], v,
  my(bb = v[3], aa = v[1]*bb/v[2], r = c12rowB(aa,aa,bb),
     dq = denominator(r[1]), dp = denominator(r[2]), f1, f2);
  f1 = if(dq==1, 1, my(f=factor(dq)); f[#f[,1],1]);
  f2 = if(dp==1, 1, my(f=factor(dp)); f[#f[,1],1]);
  printf("alpha=%d/%d b=%d: max prime den(Q) = %d, den(P) = %d,  12b+11 = %d  -> %s\n",
    v[1],v[2],bb,f1,f2,12*bb+11, if(f1<=12*bb+11 && f2<=12*bb+11,"LCM-type OK","EXCEEDS")));
}
quit;
