default(parisizemax, 8000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/01_lib.gp")
\\ PART 4: scan rational (a,b,c) of small height for p-adic overconvergence of
\\ H = a H_A + b H_B + c H_C in the x-coordinate.
PX = 430; N0 = 100; N1 = 400; HT = 6;
HA = solveode(0,0,PX); HB = solveode(1,0,PX); HC = solveode(0,1,PX);
if(HA[4]!=93 || HB[3]!=23/4 || HC[4]!=343/9, error("checks failed"));

\\ integral representatives: dn[n] = common denominator at index n
dn = vector(N1+1, ii, lcm([denominator(HA[ii]), denominator(HB[ii]), denominator(HC[ii])]));
Ai = vector(N1+1, ii, HA[ii]*dn[ii]);
Bi = vector(N1+1, ii, HB[ii]*dn[ii]);
Ci = vector(N1+1, ii, HC[ii]*dn[ii]);
vd = [vector(N1+1, ii, valuation(dn[ii],2)), vector(N1+1, ii, valuation(dn[ii],3))];

{
scan(aa, bb, cc, pk) =
  my(pp = if(pk==1,2,3), cnt=0, sx=0,sy=0,sxx=0,sxy=0, t, lo=10^9, hi=-10^9, v0=0,v1=0, z);
  for(nn = N0, N1,
    z = aa*Ai[nn+1] + bb*Bi[nn+1] + cc*Ci[nn+1];
    if(z == 0, next);
    t = valuation(z, pp) - vd[pk][nn+1];
    cnt++; sx += nn; sy += t; sxx += nn^2; sxy += nn*t;
    if(t/nn < lo, lo = t/nn); if(t/nn > hi, hi = t/nn);
    if(nn == N0, v0 = t); if(nn == N1, v1 = t);
  );
  if(cnt < 2, return([-10^9, 0, 0, 0, cnt]));
  [(cnt*sxy-sx*sy)/(cnt*sxx-sx^2)*1.0, (v1-v0)/(N1-N0)*1.0, lo*1.0, hi*1.0, cnt];
}

{
foreach([1,2], pk,
  my(pp = if(pk==1,2,3), tab = List(), res);
  for(aa = -HT, HT, for(bb = -HT, HT, for(cc = -HT, HT,
    if(aa==0 && bb==0 && cc==0, next);
    if(gcd([aa,bb,cc]) != 1, next);
    \\ up to sign: keep first nonzero positive
    if(aa < 0 || (aa==0 && (bb < 0 || (bb==0 && cc<0))), next);
    res = scan(aa, bb, cc, pk);
    listput(tab, [res[1], [aa,bb,cc], res[2], res[3], res[4], res[5]]);
  )));
  tab = vecsort(Vec(tab), 1, 4);
  print("=== PART 4: p=", pp, "  H = a H_A + b H_B + c H_C,  x-coefficients, n in [",N0,",",N1,"] ===");
  print("    scanned ", #tab, " primitive triples with entries in [-",HT,",",HT,"] (up to sign)");
  print("    top 10 by least-squares slope of v_", pp, "(H_n):");
  for(ii = 1, 10,
    print("      (a,b,c)=", tab[ii][2], "  LS=", strprintf("%+.6f", tab[ii][1]),
          "  wide=", strprintf("%+.6f", tab[ii][3]),
          "  min v/n=", strprintf("%+.6f", tab[ii][4]),
          "  max v/n=", strprintf("%+.6f", tab[ii][5]),
          "  nonzero=", tab[ii][6]);
  );
  print("    MAX least-squares slope over all triples: ", strprintf("%+.6f", tab[1][1]), " at ", tab[1][2]);
  \\ also rank by min v/n (a rigorous lower bound on the slope over the window)
  my(tb2 = vecsort(Vec(tab), 4, 4));
  print("    top 5 by min_{n in window} v_", pp, "(H_n)/n:");
  for(ii = 1, 5,
    print("      (a,b,c)=", tb2[ii][2], "  min v/n=", strprintf("%+.6f", tb2[ii][4]),
          "  LS=", strprintf("%+.6f", tb2[ii][1]));
  );
  print("    MAX wide-window slope over all triples: ",
        strprintf("%+.6f", vecsort(Vec(tab), 3, 4)[1][3]), " at ", vecsort(Vec(tab), 3, 4)[1][2]);
  print("");
);
}
