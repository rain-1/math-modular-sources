default(parisize, 10000000000);
default(realprecision, 330);
xi = eval(Str(read("b14_xi_arch_1200.txt")))*1.0;
NM=List(); VL=List(); ad(n,v)={if(abs(v)>1e-250, listput(NM,n); listput(VL,v));}
{ my(FS=[3,4,5,7,8,11,12,13,16,17,20,24,26,32,34,39,40,51,52,68,104,136,221,272]);
  foreach(FS, f,
    my(G=znstar(f,1), CH=chargalois(G));
    for(i=1,#CH,
      my(chi=CH[i]);
      if(zncharconductor(G,chi)!=f, next);
      my(L = lfuncreate([G,chi]));
      for(m=2,6,
        my(v = lfun(L, m));
        ad(Str("ReL(f",f,"#",i,",",m,")"), real(v));
        if(abs(imag(v))>1e-250, ad(Str("ImL(f",f,"#",i,",",m,")"), imag(v)))))); }
V=Vec(VL); NA=Vec(NM);
print("targets: ", #V);
print("-- singles --");
{ for(j=1,#V, my(q=bestappr(xi/V[j],10^12));
   if(q!=0 && abs(xi/V[j]-q)<1e-250, print("  HIT xi = ",q," * ",NA[j]))); }
print("-- pairs --");
{ for(j=1,#V, for(l=j+1,#V, my(r=lindep([xi,V[j],V[l]],280));
   if(type(r)=="t_COL" && #r==3 && r[1]!=0 && vecmax(abs(r))<10^15,
     print("  PAIR ",NA[j],",",NA[l]," : ",r~)))); }
print("-- with constants --");
{ my(ex=[["1",1],["z3",zeta(3)],["P2",Pi^2],["P3",Pi^3],["P4",Pi^4],["P",Pi],["s17",sqrt(17)],["z5",zeta(5)],["s13",sqrt(13)]]);
  for(j=1,#V, for(e=1,#ex, my(r=lindep([xi,V[j],ex[e][2]],280));
   if(type(r)=="t_COL" && #r==3 && r[1]!=0 && vecmax(abs(r))<10^15,
     print("  ",NA[j]," , ",ex[e][1]," : ",r~)))); }
quit
