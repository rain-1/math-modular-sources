default(parisize, 4000000000);
default(realprecision, 105);
xi = eval(Str(read("b14_xi_arch_1200.txt")))*1.0;
LV=List(); LN=List();
getvals(Ld, N,k,i,t) = {
  for(s=1,k-1,
    my(v = iferr(lfun(Ld,s), E, 0));
    if(type(v)=="t_REAL" || type(v)=="t_COMPLEX",
      if(abs(v)>1e-80, listput(LN, Str("L(",N,".",k,".",i,".",t,",",s,")")); listput(LV, real(v)))));
}
{ my(NS=[2,4,8,13,16,17,26,32,34,52,68,104,136,221,272]);
  foreach(NS, N, foreach([2,4], k,
    my(M = iferr(mfinit([N,k],0), E, 0));
    if(M==0, next);
    if(mfdim(M)==0, next);
    my(EB = mfeigenbasis(M));
    for(i=1,#EB,
      my(LL = lfunmf(M, EB[i]));
      my(ok = iferr(lfun(LL,2);1, E, 0));
      if(ok, getvals(LL,N,k,i,1),
             for(t=1,#LL, getvals(LL[t],N,k,i,t)))))); }
V=Vec(LV); NA=Vec(LN);
print("newform L-values: ", #V);
print("-- singles --");
{ for(j=1,#V, my(q=bestappr(xi/V[j],10^8));
   if(q!=0 && abs(xi/V[j]-q)<1e-80, print("  HIT xi = ",q," * ",NA[j]))); }
print("-- with constants --");
{ my(ex=[["1",1],["z3",zeta(3)],["P2",Pi^2],["P3",Pi^3],["P4",Pi^4],["P",Pi],["s17",sqrt(17)],["s13",sqrt(13)]]);
  for(j=1,#V, for(e=1,#ex, my(r=lindep([xi,V[j],ex[e][2]],280));
   if(type(r)=="t_COL" && #r==3 && r[1]!=0 && vecmax(abs(r))<10^15,
     print("  ",NA[j]," , ",ex[e][1]," : ",r~)))); }
print("-- pairs --");
{ for(j=1,#V, for(l=j+1,#V, my(r=lindep([xi,V[j],V[l]],280));
   if(type(r)=="t_COL" && #r==3 && r[1]!=0 && vecmax(abs(r))<10^15,
     print("  ",NA[j],",",NA[l]," : ",r~)))); }
quit
