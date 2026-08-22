default(parisizemax, 8000000000);
default(realprecision, 60);
read("/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/lib12.gp");
\\ Is den(Q_b), den(P_b) of LCM (partial-fraction) type, i.e. rate bounded in b?
\\ Provable partial-fraction bound for a conductor-12 half row: 12^* * D_{12b+11}^2.
{
print("=== conductor-12 product form: denominator type test ===");
print("alpha  b     nu(Q)     nu(P)     eta      does den(Q)*  |  D_{12b+11}^2 ?   (prime-to-6 parts)");
foreach([[1,2],[1,1]], fr,
  foreach([10,20,30,40,50,60], bb,
    my(nu0=fr[1],de0=fr[2]); if(bb % de0, next);
    my(aa = nu0*bb/de0, r = c12rowB(aa,aa,bb), dq = denominator(r[1]),
       dp = denominator(r[2]), dd, q6, p6, dv);
    q6 = dq/2^valuation(dq,2)/3^valuation(dq,3);
    p6 = dp/2^valuation(dp,2)/3^valuation(dp,3);
    dd = lcm(vector(12*bb+11,i,i));   \\ D_{12b+11}
    dv = if(q6 == 1, "yes", if((dd^2) % q6 == 0, "yes", if((dd^4) % q6 == 0, "D^4", "no")));
    printf("%d/%d  %3d  %8.3f  %8.3f  %8.3f   %s\n", nu0,de0,bb,
      log(q6*1.0)/bb, log(p6*1.0)/bb, log(lcm(dq,dp)*1.0)/bb, dv)));
}
quit;
