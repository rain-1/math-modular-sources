\\ 75_dcd.gp -- THE CLEAN STATEMENT:  d | c(d)  for every admissible d.
\\ It contains Conjecture 4.1 as the case d = m^2, since c(m^2) = beta_{s7}(m) and
\\ m^2 | beta(m) is exactly m^2 | c(m^2).
L = read("/home/ubuntu/code/math-modular-sources/lattice/cooper_congruence/round2/73_cd_s7_2500.txt");
{ my(tot=0, bad=List(), sq=0, mx=0, G=List());
  for(i=1,#L,
    my(d=L[i][1], v=L[i][2]);
    if(type(v)!="t_INT", next);
    tot++;
    if(v % d != 0, listput(bad,d), listput(G,[d, v/d]));
    if(issquare(d), sq++));
  print("admissible d tested: ", tot, "   (of which squares: ", sq, ")   max d = ", L[#L][1]);
  print("d | c(d) failures: ", #bad, if(#bad>0, concat("  ", Str(Vec(bad))), ""));
  print("");
  print("the quotient  G(d) := c(d)/d,  first 40:");
  print(vector(min(40,#G), i, Vec(G)[i]));
  \\ sharpness: is c(d)/d ever a p-adic unit for each p?
  my(s="");
  forprime(p=2,29,
    my(mn=99);
    for(i=1,#G, my(e=Vec(G)[i]); if(e[2]!=0 && e[1]%p==0, my(w=valuation(e[2],p)); if(w<mn,mn=w)));
    s = concat(s, concat(concat(" ",p),concat(":",if(mn==99,"-",mn)))));
  print("");
  print("min_{p | d} v_p(c(d)/d)  (0 means d | c(d) is sharp at p): ", s);
  \\ and the square case reproduces gamma
  print("");
  print("check: c(m^2)/m^2 = gamma_{s7}(m) for m <= 20:");
  my(g=read("/home/ubuntu/code/math-modular-sources/lattice/cooper_congruence/round2/gamma_s7.txt"), ok=1);
  for(m=1,20, if(m%7==0, next);
    my(f=0); for(i=1,#G, if(Vec(G)[i][1]==m^2, f=Vec(G)[i][2]));
    if(f!=g[m], ok=0; print("   mismatch at m=",m,": ",f," vs ",g[m])));
  print("   ", if(ok,"all agree","MISMATCH"));
}
quit;
