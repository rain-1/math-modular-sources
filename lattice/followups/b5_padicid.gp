/* b5_padicid.gp -- 2-adic identification sweep for xi_2(AESZ 207). */
default(parisize, 12000000000);
read("../mum_survey/apery.gp"); read("../mum_survey/ops.gp"); read("../mum_survey/lpgen.gp");
OP=0; for(i=1,#OPS, if(OPS[i][1]=="207", OP=OPS[i]));
NN=400;                      /* sigma_2=12 => ~4800 2-adic digits */
pr=aperyPair(OP[4],NN); A=pr[1]; B=pr[2];
PR = 1200;
xi2 = B[NN+1]/A[NN+1] + O(2^PR);
cau = valuation(B[NN+1]/A[NN+1] - B[NN]/A[NN], 2);
print("N=",NN,"  Cauchy 2-adic digits = ", cau, "   (using PR=",PR,")");
print("v_2(xi_2) = ", valuation(xi2,2));
print("xi_2 * 2^4 mod 2^40 = ", lift((xi2*16)+O(2^40)));

/* character catalogue: quadratic, conductor built from {2,13,17} and small others */
{ DS = [1,-3,-4,5,8,-7,-8,12,13,17,-11,-15,-19,-20,-23,-24,21,24,28,29,33,
        -35,-39,-40,-43,-51,-52,-55,-56,-68,-84,-88,-104,-119,-136,-152,
        -187,-221,-232,-260,-264,-296,-312,-408,-424,-440,44,56,60,65,68,
        69,76,85,88,89,92,93,101,104,105,109,113,120,136,140,156,168,172,
        184,204,209,221,232,248,264,272,273,280,312,408,424,440,884,1768,
        -884,-1768];
  DS = Set(DS); }
CH = List(); CN = List();
listput(CN,"1"); listput(CH, [1,[1]]);
{ for(i=1,#DS, my(D=DS[i]);
    if(D!=1 && isfundamental(D),
      listput(CN, Str("chi",D));
      listput(CH, [abs(D), vector(abs(D), a, kronecker(D,a))]))); }
CH=Vec(CH); CN=Vec(CN);
print("quadratic characters: ", #CH);

TN = List(); TV = List();
{ for(u=1,#CH, for(m=2,5,
    my(v = LpG(2, CH[u], 0, m, PR));
    if(v!=0 && valuation(v,2) < PR-40,
       listput(TN, Str("L_2(",m,",",CN[u],")")); listput(TV, v)))); }
TN=Vec(TN); TV=Vec(TV);
print("nonzero 2-adic targets: ", #TV);

print("\n--- single: xi_2 = q * T ---");
{ for(j=1,#TV,
    my(rel = lindep([xi2, TV[j]]));
    if(#rel==2 && rel[1]!=0 && abs(rel[1])<10^8 && abs(rel[2])<10^8,
      my(q=-rel[2]/rel[1]);
      if(valuation(xi2 - q*TV[j], 2) > PR-25,
        print("  HIT  xi_2 = ", q, " * ", TN[j])))); }

print("\n--- pairs: xi_2 = q1*T1 + q2*T2 ---");
{ for(j=1,#TV, for(l=j+1,#TV,
    my(rel = lindep([xi2, TV[j], TV[l]]));
    if(#rel==3 && rel[1]!=0 && vecmax(abs(rel))<10^6,
      my(q1=-rel[2]/rel[1], q2=-rel[3]/rel[1]);
      if(valuation(xi2 - q1*TV[j] - q2*TV[l], 2) > PR-40,
        print("  PAIR xi_2 = ", q1, "*", TN[j], " + ", q2, "*", TN[l]))))); }
quit
