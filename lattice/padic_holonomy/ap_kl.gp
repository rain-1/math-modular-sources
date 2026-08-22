/* ap_kl.gp -- TASK 1d: cross-check eta against Kubota-Leopoldt.
   zeta*(1-2j) = -(1 - 5^(2j-1)) B_{2j}/(2j)  (the 5-stabilised Riemann zeta).
   Taking 2j = -2 mod 4*5^r  (so 1-2j = 3 mod 4*5^r, branch 3 mod 4) the values
   converge 5-adically to zeta_5(3).  Claim to test: 2*eta = zeta_5(3).
   Exact rational Bernoulli numbers -> exact 5-adic comparison.
*/
default(parisizemax, 12000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/eta500.gp");
print("eta (from N=500 run) mod 5^12 = ", eta500 % 5^12);
te = 2*eta500;
{
for(r = 1, 7,
  my(md, jj, bz, zs, df, agree, t0);
  md = 4*5^r;
  jj = md - 2;              /* 2j = -2 mod 4*5^r, smallest positive */
  t0 = gettime();
  bz = bernfrac(jj);
  zs = -(1 - 5^(jj-1))*bz/jj;
  /* 5-adic comparison */
  df = te - zs;
  agree = if(df == 0, oo, valuation(numerator(df),5) - valuation(denominator(df),5));
  print("r=", r, "  2j=", jj, "  v_5(denom B_2j)=", valuation(denominator(bz),5),
        "  v_5(2*eta - zeta*(1-2j)) = ", agree, "   (bernfrac ms=", gettime(), ")");
  print("      zeta*(1-2j) mod 5^", r+2, " = ", lift(Mod(numerator(zs),5^(r+2))/Mod(denominator(zs),5^(r+2))),
        "     2*eta mod 5^", r+2, " = ", te % 5^(r+2));
);
}
quit;
