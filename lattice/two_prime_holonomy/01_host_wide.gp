default(parisizemax, 8000000000);
\\ PART 3: the host recurrence (Zagier row C, (a,b,c)=(10,3,9))
\\ (n+1)^2 u_{n+1} = (10n^2+10n+3) u_n - 9 n^2 u_{n-1}
NMAX = 1700;
av = vector(NMAX+2); bv = vector(NMAX+2);
av[1] = 1; av[2] = 3; bv[1] = 0; bv[2] = 1;
{
for(nn = 1, NMAX,
  av[nn+2] = ((10*nn^2+10*nn+3)*av[nn+1] - 9*nn^2*av[nn])/(nn+1)^2;
  bv[nn+2] = ((10*nn^2+10*nn+3)*bv[nn+1] - 9*nn^2*bv[nn])/(nn+1)^2;
);
}
print("a_0..a_5 = ", vector(6, ii, av[ii]));
print("b_0..b_5 = ", vector(6, ii, bv[ii]));
print("");

print("=== 3(i) increments  d_n = b_n/a_n - b_{n-1}/a_{n-1} ===");
dv = vector(NMAX+1);
{ for(nn = 1, NMAX, dv[nn+1] = bv[nn+2-1]/av[nn+2-1] - bv[nn+1-1]/av[nn+1-1]); }
{
foreach([2,3], pp,
  print("  p=",pp,":  v_",pp,"(d_n)  and v_",pp,"(d_n)/n");
  foreach([10,25,50,100,150,200,250,300,350,400,500,600,700,800], nn,
    my(t = valuation(dv[nn+1], pp));
    print("     n=", nn, "  v=", t, "  v/n=", strprintf("%.6f", t*1.0/nn));
  );
);
}
print("");
print("  exact: v_3(d_n) - 2n  for n=100..110: ", vector(11, ii, valuation(dv[100+ii], 3) - 2*(99+ii)));
print("  exact: v_2(d_n)       for n=100..110: ", vector(11, ii, valuation(dv[100+ii], 2)));
print("");

\\ ---- 3(ii)/(iii): slope of v_p of coefficients of B - r A ----
\\ slope measured on n in [N0,N1] by least squares and by (v(N1)-v(N0))/(N1-N0)
N0 = 200; N1 = 800;
{
slopes(rr, pp) =
  my(cnt=0, sx=0, sy=0, sxx=0, sxy=0, t, lo, hi, v0=0, v1=0);
  lo = 10^9; hi = -10^9;
  for(nn = N0, N1,
    t = bv[nn+1] - rr*av[nn+1];
    if(t == 0, next);
    t = valuation(t, pp);
    cnt++; sx += nn; sy += t; sxx += nn^2; sxy += nn*t;
    if(t/nn < lo, lo = t/nn); if(t/nn > hi, hi = t/nn);
    if(nn == N0, v0 = t); if(nn == N1, v1 = t);
  );
  if(cnt < 2, return([0,0,0,0]));
  [(cnt*sxy-sx*sy)/(cnt*sxx-sx^2)*1.0, (v1-v0)/(N1-N0)*1.0, lo*1.0, hi*1.0];
}

{
foreach([3,2], pp,
  print("=== 3(", if(pp==3,"ii","iii"), "):  p=", pp, "  slopes of v_", pp, "(b_n - r a_n), n in [",N0,",",N1,"] ===");
  my(best = -10^9, bestr = 0, res, tab = List());
  forstep(qq = 1, 6, 1,
    for(numr = -6, 6,
      if(gcd(abs(numr), qq) != 1 && numr != 0, next);
      if(numr == 0 && qq != 1, next);
      my(rr = numr/qq);
      res = slopes(rr, pp);
      listput(tab, [res[1], rr, res[2], res[3], res[4]]);
      if(res[1] > best, best = res[1]; bestr = rr);
    );
  );
  tab = vecsort(Vec(tab), 1, 4);
  print("   top 8 rationals r by LS slope:");
  for(ii = 1, 8,
    print("     r=", tab[ii][2], "   LS=", strprintf("%+.6f", tab[ii][1]),
          "  wide=", strprintf("%+.6f", tab[ii][3]),
          "  min v/n=", strprintf("%+.6f", tab[ii][4]),
          "  max v/n=", strprintf("%+.6f", tab[ii][5]));
  );
  print("   BEST rational r = ", bestr, "  LS slope = ", strprintf("%+.6f", best));
  \\ now r = xi_p computed p-adically as lim b_n/a_n
  my(xi = bv[NMAX+1]/av[NMAX+1], res2 = slopes(xi, pp));
  print("   r = xi_",pp," := b_",NMAX,"/a_",NMAX," (p-adic limit surrogate):  LS=",
        strprintf("%+.6f", res2[1]), "  wide=", strprintf("%+.6f", res2[2]),
        "  min v/n=", strprintf("%+.6f", res2[3]), "  max v/n=", strprintf("%+.6f", res2[4]));
  print("   v_",pp,"(b_n - xi a_n) at n=100,...,800: ",
        vector(8, ii, valuation(bv[100*ii+1] - xi*av[100*ii+1], pp)));
  print("");
);
}

\\ 3-adic digits of xi_3 = lim b_n/a_n, and its stability
{
my(xi8 = bv[NMAX+1]/av[NMAX+1], xi4 = bv[401]/av[401], xi2 = bv[201]/av[201]);
print("xi_3 = lim b_n/a_n, 3-adically:");
print("   b_800/a_800 + O(3^40) = ", xi8 + O(3^40));
print("   b_400/a_400 + O(3^40) = ", xi4 + O(3^40));
print("   v_3(b_800/a_800 - b_400/a_400) = ", valuation(xi8-xi4, 3));
print("   v_3(b_400/a_400 - b_200/a_200) = ", valuation(xi4-xi2, 3));
print("   v_2(b_800/a_800 - b_400/a_400) = ", valuation(xi8-xi4, 2));
print("v_2(a_n) n=100,200,400,800: ", [valuation(av[101],2),valuation(av[201],2),valuation(av[401],2),valuation(av[801],2)]);
print("v_3(a_n) n=100,200,400,800: ", [valuation(av[101],3),valuation(av[201],3),valuation(av[401],3),valuation(av[801],3)]);
print("v_2(b_n) n=100,200,400,800: ", [valuation(bv[101],2),valuation(bv[201],2),valuation(bv[401],2),valuation(bv[801],2)]);
print("v_3(b_n) n=100,200,400,800: ", [valuation(bv[101],3),valuation(bv[201],3),valuation(bv[401],3),valuation(bv[801],3)]);
}
