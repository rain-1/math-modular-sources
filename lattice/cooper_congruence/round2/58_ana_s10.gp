\\ 58_ana_s10.gp -- divisibility analysis of c(d) for s10.  Raw data 57_cd_s10.txt holds
\\ A(d) = T(d)/i.  Test both normalisations c = -A (lambda=i) and c = -2A (lambda=2i).
read("50_lib.gp");
V = read("57_cd_s10.txt");
bet = read("20_beta_s10.txt");
print("s10: ", #V, " admissible d");
{ print("  non-integer / NOREP / DEGEN entries: ",
    select(x->type(x[4])!="t_INT", V)); }
W = select(x->type(x[4])=="t_INT", V);
print("  usable: ", #W);
print();
print("== calibration on squares: c(m^2) vs beta(m) ==");
{ for(i=1,#W, my(d=W[i][1], A=W[i][4]);
   if(!issquare(d), next);
   my(m=sqrtint(d));
   print("   d=",d," m=",m,"   -A=",-A,"   -2A=",-2*A,"   beta(m)=",bet[m],
         "   -A=beta? ",-A==bet[m],"   -2A=beta? ",-2*A==bet[m])); }
print();
print("== divisibility d | c(d) ==");
{ for(lam=1,2,
   my(fails=List());
   for(i=1,#W, my(d=W[i][1], c=-lam*W[i][4]);
     if(c%d != 0, listput(fails,[d,c])));
   print("  lambda = ",if(lam==1,"i","2i"),"  (c = ",if(lam==1,"-A","-2A"),")   tests=",#W,
         "   failures=",#fails, if(#fails>0, concat("  first: ",Str(Vec(fails)[1])), "")));
}
print();
print("== with lambda = 2i:  G(d) = c(d)/d, and sharpness min_{p|d} v_p(G(d)) ==");
{ my(G=List(), DS=List());
  for(i=1,#W, my(d=W[i][1], c=-2*W[i][4]); if(c%d==0, listput(G,c/d); listput(DS,d)));
  G=Vec(G); DS=Vec(DS);
  print("  G(d) for the first 24 admissible d:");
  print("   d = ", vector(min(24,#DS),i,DS[i]));
  print("   G = ", vector(min(24,#G),i,G[i]));
  forprime(p=2,31,
    my(mn=oo, at=0, n=0);
    for(i=1,#DS, if(DS[i]%p!=0 || G[i]==0, next); n++;
      my(v=valuation(G[i],p)); if(v<mn, mn=v; at=DS[i]));
    print("   p=",p,"   min_{p|d, G!=0} v_p(G(d)) = ",mn,"  (at d=",at,")   [",n," tests]"));
}
print();
print("== p-part form: p^(2a) | c(p^(2a) d) ==");
{ my(M=Map(), tot=0, bad=List());
  for(i=1,#W, mapput(M, W[i][1], -2*W[i][4]));
  for(i=1,#W, my(d=W[i][1]);
    forprime(p=2,23,
      my(e=p^2*d);
      if(!mapisdefined(M,e), next);
      tot++;
      if(mapget(M,e) % (p^2*if(d%p==0,p^0,1)) != 0, listput(bad,[p,d]))));
  print("  tests=",tot,"   failures=",if(#bad==0,"NONE",Vec(bad)));
}
quit;
