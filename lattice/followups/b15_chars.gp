/* b15_chars.gp -- xi_infty(207), 1139 digits, against ALL Dirichlet characters of
   conductor dividing 8*13*17 (and a few others), m=2..6, real and imaginary parts. */
default(realprecision, 1150);
xi = eval(Str(read("b14_xi_arch_1200.txt")));
default(realprecision, 1100); xi = xi*1.0;
print("xi to 1100 digits loaded.");
NM=List(); VL=List(); ad(n,v)={if(abs(v)>1e-900, listput(NM,n); listput(VL,v));}
{ my(FS=[3,4,5,7,8,11,12,13,16,17,20,24,26,32,34,39,40,51,52,68,104,136,221,272,442,884,1768]);
  foreach(FS, f,
    my(G=znstar(f,1), CH=chargalois(G));
    for(i=1,#CH,
      my(chi=CH[i], L);
      if(zncharconductor(G,chi)!=f, next);
      L = lfuncreate([G,chi]);
      for(m=2,6,
        my(v = lfun(L, m));
        ad(Str("ReL(",f,"#",i,",",m,")"), real(v));
        if(abs(imag(v))>1e-900, ad(Str("ImL(",f,"#",i,",",m,")"), imag(v)))))); }
V=Vec(VL); NA=Vec(NM);
print("targets: ", #V);
print("\n-- singles: xi = q*T --");
{ for(j=1,#V, my(q=bestappr(xi/V[j],10^25));
   if(q!=0 && abs(xi/V[j]-q)<1e-900, print("  HIT xi = ",q," * ",NA[j]))); }
print("\n-- pairs --");
{ for(j=1,#V, for(l=j+1,#V, my(r=lindep([xi,V[j],V[l]],900));
   if(type(r)=="t_COL" && #r==3 && r[1]!=0 && vecmax(abs(r))<10^20,
     print("  PAIR ",NA[j],",",NA[l]," : ",r~)))); }
print("\n-- with Pi powers and zeta(3) --");
{ my(ex=[["1",1],["z3",zeta(3)],["P2",Pi^2],["P3",Pi^3],["P4",Pi^4],["P",Pi],["s17",sqrt(17)],["z5",zeta(5)]]);
  for(j=1,#V, for(e=1,#ex, my(r=lindep([xi,V[j],ex[e][2]],900));
   if(type(r)=="t_COL" && #r==3 && r[1]!=0 && vecmax(abs(r))<10^20,
     print("  ",NA[j]," , ",ex[e][1]," : ",r~)))); }
quit
