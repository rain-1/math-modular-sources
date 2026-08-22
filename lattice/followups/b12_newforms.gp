/* b12_newforms.gp -- test xi_infty(207) against L(f,s) for newforms of weight 2 and 4
   at levels built from the bad primes 2, 13, 17.                                   */
default(parisize, 8000000000);
default(realprecision, 340);
xi0 = eval(Str(read("b2_xi_arch.txt"))); xi = xi0*1.0;
LV = List(); LN = List();
{ my(NS=[2,4,8,13,16,17,26,32,34,52,64,68,104,128,136,208,221,256,272,416,442,544,832,884]);
  foreach(NS, N,
    foreach([2,4], k,
      my(M);
      M = mfinit([N,k],0);
      if(mfdim(M)==0, next);
      my(EB = mfeigenbasis(M));
      for(i=1,#EB,
        my(f=EB[i], L);
        L = lfunmf(M,f);
        if(type(L)=="t_VEC" && #L>0 && type(L[1])=="t_VEC", L=L, L=[L]);
        for(t=1,#L,
          for(s=1,k-1,
            my(v = lfun(L[t], s));
            if(abs(v)>1e-200,
              listput(LN, Str("L(f_",N,"_",k,"_",i,"_",t,",",s,")"));
              listput(LV, v)))))));}
V=Vec(LV); NA=Vec(LN);
print("newform L-values collected: ", #V);
print("\n-- singles --");
{ for(j=1,#V, my(q=bestappr(xi/V[j],10^12));
   if(q!=0 && abs(xi/V[j]-q)<1e-280, print("  xi = ",q," * ",NA[j]))); }
print("\n-- pairs among newform values --");
{ for(j=1,#V, for(l=j+1,#V, my(r=lindep([xi,V[j],V[l]],290));
   if(type(r)=="t_COL" && #r==3 && r[1]!=0 && vecmax(abs(r))<10^8,
     print("  ",NA[j]," , ",NA[l]," : ",r~)))); }
print("\n-- newform value with zeta(3), Pi^2, Pi^4, 1 --");
{ my(ex=[["1",1],["zeta(3)",zeta(3)],["Pi^2",Pi^2],["Pi^4",Pi^4],["Pi^3",Pi^3],["s17",sqrt(17)],["Pi",Pi]]);
  for(j=1,#V, for(e=1,#ex, my(r=lindep([xi,V[j],ex[e][2]],290));
   if(type(r)=="t_COL" && #r==3 && r[1]!=0 && vecmax(abs(r))<10^8,
     print("  ",NA[j]," , ",ex[e][1]," : ",r~)))); }
quit
