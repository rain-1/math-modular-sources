default(parisizemax, 8000000000);
default(realprecision, 80);
read("/home/ubuntu/code/math-modular-sources/lattice/euler_criterion/lp.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/lib12.gp");
raw(aa, bb) =
{ my(r = c12rowB(aa,aa,bb), q = r[1], p = r[2], dn = denominator(q),
     v2 = valuation(dn,2), v3 = valuation(dn,3), dp = denominator(p));
  [ log(abs(q*1.0)), log(abs(q*Catalan - p)), v2, v3,
    log((dn/2^v2/3^v3)*1.0), log(lcm(dn,dp)*1.0) ];
};
{
print("=== conductor-12 Catalan row (product form): wide-window rates, b: 30 -> 60 ===");
print("alpha    logLam    loglam    kap2     kap3     nu       eta      largest prime in den(Q) at b=60");
foreach([[1,4],[1,2],[3,4],[1,1],[5,4],[3,2],[7,4],[2,1]], fr,
  my(nu0 = fr[1], de0 = fr[2], b1 = 30, b2 = 60, u1, u2, d, lp0);
  if(b1 % de0 || b2 % de0, next);
  u1 = raw(nu0*b1/de0, b1); u2 = raw(nu0*b2/de0, b2); d = b2-b1;
  lp0 = my(dn = denominator(c12rowB(nu0*b2/de0,nu0*b2/de0,b2)[1]));
  lp0 = if(dn==1, 1, my(f=factor(dn)); f[#f[,1],1]);
  printf("%d/%d    %8.4f  %8.4f  %7.3f  %7.3f  %7.3f  %7.3f   %d   (12b+11 = %d)\n", nu0,de0,
     (u2[1]-u1[1])/d, (u2[2]-u1[2])/d, (u2[3]-u1[3])/d, (u2[4]-u1[4])/d,
     (u2[5]-u1[5])/d, (u2[6]-u1[6])/d, lp0, 12*b2+11));
print("");
print("=== p-adic limits vs Theorem F ===");
print("  prediction: xi_2 = zeta_2(2)   (E_2(2)=1, chi_-4 omega^-1 = 1 at p=2)");
print("              xi_3 = (9/10) L_3(2,chi_12)   (E_3(2)=10/9, chi_-4 omega^-1 = chi_12 at p=3)");
my(PR2 = 1200, PR3 = 500, z2 = Lp(2,triv,2,PR2), l3 = Lp(3,chi12,2,PR3));
print("alpha  b    v2(xi-zeta_2(2))  Cauchy2    v3(xi-(9/10)L_3(2,chi12))  Cauchy3");
foreach([[1,2],[1,1],[3,2],[2,1]], fr,
  foreach([30,45,60], bb,
    my(nu0=fr[1], de0=fr[2]); if(bb % de0, next);
    my(aa = nu0*bb/de0, r = c12rowB(aa,aa,bb), rm, x, xm);
    x = r[2]/r[1];
    rm = c12rowB(nu0*(bb-de0)/de0, nu0*(bb-de0)/de0, bb-de0); xm = rm[2]/rm[1];
    printf("%d/%d  %3d   %8d  %8d    %8d  %8d\n", nu0,de0,bb,
      valuation((x+O(2^PR2)) - z2, 2), valuation(x-xm,2)/de0*de0,
      valuation((x+O(3^PR3)) - (9/10)*l3, 3), valuation(x-xm,3))));
print("");
print("=== controls at alpha=1, b=60 ===");
my(r = c12rowB(60,60,60), x = r[2]/r[1]);
foreach([1/2,1,2,4,5/4,4/5,8/5,9/10,10/9,3/2,2/3,5/8], sc,
  printf("  c=%6s : v2(xi - c*zeta_2(2)) = %5d   v3(xi - c*L_3(2,chi12)) = %5d\n",
    sc, valuation((x+O(2^PR2))-sc*z2,2), valuation((x+O(3^PR3))-sc*l3,3)));
}
quit;
