/* b6_arch2.gp -- broader archimedean search: algdep, logs of singular values,
   sqrt(17)/sqrt(13) twists, 3-term combos over a curated weight-<=4 basis.  */
default(realprecision, 400);
xi0 = eval(Str(read("b2_xi_arch.txt")));
default(realprecision, 340);
xi = xi0*1.0;
s17 = sqrt(17); s13=sqrt(13); s2=sqrt(2);
zp = (349+85*s17)/131072; zm = (349-85*s17)/131072; zc = 1/53248.;
print("singular z: ", zm, " ", zc, " ", zp);
print("1/z: ", 1/zm, " ", 1/zc, " ", 1/zp);
print("\n--- algdep tests ---");
for(d=1,8, my(P=algdep(xi,d)); if(poldegree(P)<=d && vecmax(abs(Vec(P)))<10^25, print("  d=",d,": ",P)));
print("\n--- curated basis, 3-term lindep ---");
NM=List(); VL=List();
ad(n,v)={listput(NM,n);listput(VL,v);}
ad("zeta(3)",zeta(3)); ad("Pi^2",Pi^2); ad("Pi^4",Pi^4); ad("Pi^3",Pi^3);
ad("zeta(5)",zeta(5)); ad("Pi^2*zeta(3)",Pi^2*zeta(3));
ad("log2",log(2)); ad("log2^2",log(2)^2); ad("log2^3",log(2)^3); ad("log2^4",log(2)^4);
ad("log13",log(13)); ad("log17",log(17));
ad("s17",s17); ad("s17*zeta(3)",s17*zeta(3)); ad("s17*Pi^2",s17*Pi^2);
ad("s17*Pi^4",s17*Pi^4); ad("s17*Pi^3",s17*Pi^3);
ad("s13",s13); ad("s13*zeta(3)",s13*zeta(3));
ad("L(17,2)",lfun(17,2)); ad("L(17,3)",lfun(17,3)); ad("L(17,4)",lfun(17,4));
ad("L(13,2)",lfun(13,2)); ad("L(13,3)",lfun(13,3)); ad("L(13,4)",lfun(13,4));
ad("L(-4,2)",lfun(-4,2)); ad("L(-4,3)",lfun(-4,3)); ad("L(-4,4)",lfun(-4,4));
ad("L(8,2)",lfun(8,2)); ad("L(8,3)",lfun(8,3)); ad("L(8,4)",lfun(8,4));
ad("L(-8,2)",lfun(-8,2)); ad("L(-8,3)",lfun(-8,3)); ad("L(-8,4)",lfun(-8,4));
ad("L(-68,3)",lfun(-68,3)); ad("L(-68,2)",lfun(-68,2)); ad("L(-68,4)",lfun(-68,4));
ad("L(-52,3)",lfun(-52,3)); ad("L(-52,2)",lfun(-52,2)); ad("L(-52,4)",lfun(-52,4));
ad("L(221,3)",lfun(221,3)); ad("L(221,2)",lfun(221,2)); ad("L(221,4)",lfun(221,4));
ad("L(104,3)",lfun(104,3)); ad("L(-104,3)",lfun(-104,3));
ad("L(136,3)",lfun(136,3)); ad("L(-136,3)",lfun(-136,3));
ad("Li4(1/2)",polylog(4,1/2)); ad("Catalan",Catalan);
ad("log(zm_abs)",log(abs(zm))); ad("log(zp)",log(zp)); ad("log(zc)",log(zc));
V=Vec(VL); NA=Vec(NM);
print("basis size ", #V);
{ for(j=1,#V, for(l=j+1,#V,
    my(r=lindep([xi,V[j],V[l]],300));
    if(type(r)=="t_COL" && #r==3 && r[1]!=0 && vecmax(abs(r))<10^10,
       print("  PAIR ",NA[j]," , ",NA[l]," : ",r~)))); }
print("\n--- triples over a small core ---");
{ my(core=[1,2,3,5,7,13,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46]);
  for(a=1,#core, for(b=a+1,#core, for(c=b+1,#core,
    my(j=core[a],l=core[b],k=core[c]);
    my(r=lindep([xi,V[j],V[l],V[k]],280));
    if(type(r)=="t_COL" && #r==4 && r[1]!=0 && vecmax(abs(r))<10^6,
       print("  TRIPLE ",NA[j],",",NA[l],",",NA[k]," : ",r~))))); }
quit
