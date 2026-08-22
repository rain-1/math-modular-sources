/* b16_mf.gp -- newform L-values, weights 2,3,4,6, levels from {2,13,17}. */
default(parisize, 10000000000);
default(realprecision, 1100);
xi = eval(Str(read("b14_xi_arch_1200.txt")))*1.0;
LV=List(); LN=List();
{ my(NS=[2,4,8,13,16,17,26,32,34,52,64,68,104,128,136,208,221,256,272,416,442,544]);
  foreach(NS, N, foreach([2,3,4,6], k,
    my(M = mfinit([N,k],0));
    if(mfdim(M)==0, next);
    my(EB = mfeigenbasis(M));
    for(i=1,#EB,
      my(LL = lfunmf(M, EB[i]));
      if(type(LL)!="t_VEC", LL=[LL]);
      for(t=1,#LL,
        my(Ld);
        Ld = LL[t];
        for(s=1,k-1,
          my(v = lfun(Ld, s));
          if(abs(v)>1e-900,
            listput(LN, Str("L(",N,".",k,".",i,".",t,",",s,")")); listput(LV, v))))))); }
V=Vec(LV); NA=Vec(LN);
print("newform L-values: ", #V);
{ for(j=1,#V, my(q=bestappr(xi/V[j],10^25));
   if(q!=0 && abs(xi/V[j]-q)<1e-900, print("  HIT xi = ",q," * ",NA[j]))); }
print("-- pairs with Pi powers --");
{ my(ex=[["1",1],["z3",zeta(3)],["P2",Pi^2],["P3",Pi^3],["P4",Pi^4],["P",Pi],["s17",sqrt(17)],["s13",sqrt(13)]]);
  for(j=1,#V, for(e=1,#ex, my(r=lindep([xi,V[j],ex[e][2]],900));
   if(type(r)=="t_COL" && #r==3 && r[1]!=0 && vecmax(abs(r))<10^20,
     print("  ",NA[j]," , ",ex[e][1]," : ",r~)))); }
print("-- pairs among newform values --");
{ for(j=1,#V, for(l=j+1,#V, my(r=lindep([xi,V[j],V[l]],900));
   if(type(r)=="t_COL" && #r==3 && r[1]!=0 && vecmax(abs(r))<10^20,
     print("  ",NA[j],",",NA[l]," : ",r~)))); }
quit
