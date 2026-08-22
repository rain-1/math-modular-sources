default(parisize, 6000000000);
read("../mum_survey/apery.gp"); read("../mum_survey/ops.gp"); read("../mum_survey/lpgen.gp");
OP=0; for(i=1,#OPS, if(OPS[i][1]=="207", OP=OPS[i]));
NN=500; pr=aperyPair(OP[4],NN); A=pr[1]; B=pr[2];
PR=5500;
xi2 = B[NN+1]/A[NN+1] + O(2^PR);
print("Cauchy = ", valuation(B[NN+1]/A[NN+1]-B[NN]/A[NN],2), " ; using PR=",PR);
r17 = sqrt(17+O(2^(PR+5)));
print("\n-- xi_2 in Q(sqrt17)? --");
{ my(r=lindep([xi2,1,r17]));
  print("  lindep(xi2,1,sqrt17) = ", r~);
  if(#r==3 && r[1]!=0, my(q1=-r[2]/r[1], q2=-r[3]/r[1]);
    print("  residual v_2 = ", valuation(xi2-q1-q2*r17,2))); }
print("\n-- algdep over Q_2, degrees 1..6 (residual valuation reported) --");
{ for(d=1,6, my(P=algdep(xi2,d));
    print("  d=",d,"  height=", vecmax(abs(Vec(P))), "  v_2(P(xi2)) = ", valuation(subst(P,x,xi2),2))); }
print("\n-- xi_2 vs 2-adic zeta values with 2-power rational coefficients --");
{ my(T=[["z2(2)",LpG(2,triv,0,2,PR)],["z2(3)",LpG(2,triv,0,3,PR)],["z2(4)",LpG(2,triv,0,4,PR)],
        ["z2(5)",LpG(2,triv,0,5,PR)],["z2(6)",LpG(2,triv,0,6,PR)],["z2(7)",LpG(2,triv,0,7,PR)],
        ["L2(2,chi8)",LpG(2,[8,[1,0,0,0,-1,0,0,0]],0,2,PR)],
        ["L2(3,chi8)",LpG(2,[8,[1,0,0,0,-1,0,0,0]],0,3,PR)],
        ["L2(2,chi17)",LpG(2,[17,vector(17,a,kronecker(17,a))],0,2,PR)],
        ["L2(3,chi17)",LpG(2,[17,vector(17,a,kronecker(17,a))],0,3,PR)],
        ["L2(4,chi17)",LpG(2,[17,vector(17,a,kronecker(17,a))],0,4,PR)],
        ["L2(2,chi13)",LpG(2,[13,vector(13,a,kronecker(13,a))],0,2,PR)],
        ["L2(3,chi13)",LpG(2,[13,vector(13,a,kronecker(13,a))],0,3,PR)],
        ["L2(4,chi13)",LpG(2,[13,vector(13,a,kronecker(13,a))],0,4,PR)]]);
  for(j=1,#T, if(T[j][2]==0, print("  ",T[j][1]," = 0"); next);
    print("  ",T[j][1],"  v_2=",valuation(T[j][2],2),
          "   lindep(xi2,T) = ", lindep([xi2,T[j][2]])~));
  print("\n  pairs with sqrt17 coefficients:");
  for(j=1,#T, if(T[j][2]==0,next);
    my(r=lindep([xi2, T[j][2], r17*T[j][2]]));
    if(#r==3 && r[1]!=0 && vecmax(abs(r))<10^12, print("   *** ",T[j][1]," : ",r~)));
}
quit
