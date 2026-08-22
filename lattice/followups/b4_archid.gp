/* b4_archid.gp -- identification sweep for xi_infty(AESZ 207). */
default(realprecision, 400);
xi0 = eval(Str(read("b2_xi_arch.txt")));
default(realprecision, 345);
xi = xi0 * 1.0;
print("xi  = ", xi);
print("1/xi= ", 1/xi);
print("digits available ~348");

NAMES = List(); VALS = List();
add(nm, v) = {listput(NAMES,nm); listput(VALS,v);}

/* zeta / pi powers */
add("zeta(3)", zeta(3)); add("zeta(5)", zeta(5)); add("zeta(7)", zeta(7));
add("Pi^2", Pi^2); add("Pi^3", Pi^3); add("Pi^4", Pi^4); add("Pi^5", Pi^5); add("Pi^6", Pi^6);
add("zeta(3)^2", zeta(3)^2); add("Pi^2*zeta(3)", Pi^2*zeta(3));
add("log(2)", log(2)); add("log(13)", log(13)); add("log(17)", log(17));
add("zeta(3)*log2", zeta(3)*log(2)); add("Pi^2*log2^2", Pi^2*log(2)^2);
add("log2^4", log(2)^4); add("Li4(1/2)", polylog(4,1/2)); add("Li5(1/2)", polylog(5,1/2));
add("Li4(-1)", polylog(4,-1)); add("Li3(1/2)", polylog(3,1/2));
add("Catalan", Catalan);

/* Dirichlet L-values, quadratic characters */
{ my(DS=[-3,-4,5,8,-7,-8,12,13,17,-11,-15,-19,-20,-23,-24,21,24,28,29,33,
         -35,-39,-40,-43,-51,-52,-55,-56,-68,-84,-88,-104,-119,-136,-152,
         -187,-221,-232,-260,-264,-296,-312,-408,-424,-440,44,56,60,65,68,
         69,76,85,88,89,92,93,101,104,105,109,113,120,136,140,156,168,172,
         184,204,209,221,232,248,264,272,273,280,312,408,424,440,884,1768,
         -884,-1768,-136,-1768]);
  DS = Set(DS);
  for(i=1,#DS, my(D=DS[i]);
    if(isfundamental(D),
      for(m=2,5, add(Str("L(",D,",",m,")"), lfun(D,m))))); }

VV = Vec(VALS); NN = Vec(NAMES);
print("catalogue size = ", #VV);

print("\n--- single: xi = q * C, q rational of small height ---");
{ for(j=1,#VV,
    if(VV[j]==0, next);
    my(r = xi/VV[j], q = bestappr(r, 10^12));
    if(q!=0 && abs(r-q) < abs(r)*1e-300, print("  HIT  xi = ", q, " * ", NN[j]))); }

print("\n--- pairs: xi = q1*C1 + q2*C2 ---");
{ for(j=1,#VV, for(l=j+1,#VV,
    my(v = lindep([xi, VV[j], VV[l]], 300));
    if(type(v)=="t_COL" && #v==3 && v[1]!=0 && vecmax(abs(v))<10^9,
       print("  PAIR ", NN[j], " , ", NN[l], "  : ", v~)))); }
quit
