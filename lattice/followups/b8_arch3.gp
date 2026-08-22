/* b8_arch3.gp -- archimedean sweep including the constant 1 and 3-term combos. */
default(realprecision, 400);
xi0 = eval(Str(read("b2_xi_arch.txt")));
default(realprecision, 340);
xi = xi0*1.0;
NM=List(); VL=List(); ad(n,v)={listput(NM,n);listput(VL,v);}
ad("1",1);
ad("zeta(3)",zeta(3)); ad("zeta(5)",zeta(5)); ad("zeta(7)",zeta(7));
ad("Pi^2",Pi^2); ad("Pi^3",Pi^3); ad("Pi^4",Pi^4); ad("Pi^5",Pi^5); ad("Pi^6",Pi^6);
ad("Pi^2*zeta(3)",Pi^2*zeta(3)); ad("zeta(3)^2",zeta(3)^2);
ad("log2",log(2)); ad("log2^2",log(2)^2); ad("log2^3",log(2)^3); ad("log2^4",log(2)^4);
ad("log13",log(13)); ad("log17",log(17)); ad("Li4(1/2)",polylog(4,1/2));
ad("zeta(3)log2",zeta(3)*log(2)); ad("Pi^2log2^2",Pi^2*log(2)^2);
ad("s17",sqrt(17)); ad("s13",sqrt(13)); ad("s2",sqrt(2));
ad("s17*zeta(3)",sqrt(17)*zeta(3)); ad("s17*Pi^2",sqrt(17)*Pi^2);
ad("s17*Pi^3",sqrt(17)*Pi^3); ad("s17*Pi^4",sqrt(17)*Pi^4);
ad("Catalan",Catalan); ad("G*Pi",Catalan*Pi);
{ my(DS=[-3,-4,5,8,-7,-8,12,13,17,-11,-19,-20,-24,21,24,29,33,-39,-40,-51,-52,
         -56,-68,-88,-104,-119,-136,-152,-187,-221,-232,-264,-296,-312,
         44,56,60,65,68,69,76,85,88,89,92,101,104,105,113,120,136,140,156,
         168,172,184,204,209,221,232,248,264,272,273,280,312,408,424,440,884,1768]);
  DS=Set(DS);
  for(i=1,#DS, my(D=DS[i]); if(isfundamental(D),
    for(m=2,5, ad(Str("L(",D,",",m,")"), lfun(D,m))))); }
V=Vec(VL); NA=Vec(NM);
print("basis size ", #V, "   (index 1 = constant 1)");
print("\n--- pairs including 1 ---");
{ for(l=2,#V,
    my(r=lindep([xi,1,V[l]],300));
    if(type(r)=="t_COL" && #r==3 && r[1]!=0 && vecmax(abs(r))<10^12,
       print("  ", NA[l], " : ", r~))); }
print("\n--- all pairs ---");
{ for(j=2,#V, for(l=j+1,#V,
    my(r=lindep([xi,V[j],V[l]],300));
    if(type(r)=="t_COL" && #r==3 && r[1]!=0 && vecmax(abs(r))<10^9,
       print("  ",NA[j]," , ",NA[l]," : ",r~)))); }
print("\n--- triples with 1 ---");
{ for(j=2,#V, for(l=j+1,#V,
    my(r=lindep([xi,1,V[j],V[l]],290));
    if(type(r)=="t_COL" && #r==4 && r[1]!=0 && vecmax(abs(r))<10^7,
       print("  1,",NA[j],",",NA[l]," : ",r~)))); }
quit
