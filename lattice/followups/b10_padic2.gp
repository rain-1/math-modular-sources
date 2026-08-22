/* enriched 2-adic basis: zeta_2 odd, logs of units, Gamma_2 values, products. */
default(parisize, 12000000000);
read("../mum_survey/apery.gp"); read("../mum_survey/ops.gp"); read("../mum_survey/lpgen.gp");
OP=0; for(i=1,#OPS, if(OPS[i][1]=="207", OP=OPS[i]));
NN=400; pr=aperyPair(OP[4],NN); A=pr[1]; B=pr[2];
PR=900;
xi2 = B[NN+1]/A[NN+1] + O(2^PR);
print("v_2(xi2)=",valuation(xi2,2));
print("\n-- is xi_2 rational? --");
print("  lindep([xi2,1]) = ", lindep([xi2,1]));
print("  algdep deg<=4: ", algdep(xi2,4));
NM=List(); VL=List(); ad(n,v)={if(v!=0, listput(NM,n);listput(VL,v));}
ad("zeta_2(3)", LpG(2,triv,0,3,PR));
ad("zeta_2(5)", LpG(2,triv,0,5,PR));
ad("zeta_2(7)", LpG(2,triv,0,7,PR));
ad("L_2(2,chi-4)", LpG(2,chim4,0,2,PR));
ad("L_2(4,chi-4)", LpG(2,chim4,0,4,PR));
ad("L_2(2,chi8)", LpG(2,[8,[1,0,0,0,-1,0,0,0]],0,2,PR));
{ my(u=[3,5,7,11,13,17,29,349,-1,-3,-5]);
  for(i=1,#u, my(x=u[i]);
    if(x%2!=0, ad(Str("log2(",x,")"), log(x+O(2^PR))))); }
/* Gamma_2 via Morita on Z_2: Gamma_2(x) for x in Z_2 -- use the standard product */
gam2(x, PR) = { my(s=1+O(2^PR), y=x+O(2^PR)); 0; }
V=Vec(VL); NA=Vec(NM);
print("\nbasis: ", NA);
print("\n-- singles --");
{ for(j=1,#V, my(rel=lindep([xi2,V[j]]));
   if(#rel==2 && rel[1]!=0 && vecmax(abs(rel))<10^8,
     my(q=-rel[2]/rel[1]); if(valuation(xi2-q*V[j],2)>PR-25, print("  xi2 = ",q," * ",NA[j])))); }
print("\n-- pairs --");
{ for(j=1,#V, for(l=j+1,#V, my(rel=lindep([xi2,V[j],V[l]]));
   if(#rel==3 && rel[1]!=0 && vecmax(abs(rel))<10^6,
     my(q1=-rel[2]/rel[1],q2=-rel[3]/rel[1]);
     if(valuation(xi2-q1*V[j]-q2*V[l],2)>PR-40, print("  xi2 = ",q1,"*",NA[j]," + ",q2,"*",NA[l]))))); }
print("\n-- products log*log, log*zeta --");
NM2=List(); VL2=List();
{ for(j=1,#V, for(l=j,#V, listput(NM2,Str(NA[j],"*",NA[l])); listput(VL2,V[j]*V[l]))); }
V2=Vec(VL2); NA2=Vec(NM2);
{ for(j=1,#V2, my(rel=lindep([xi2,V2[j]]));
   if(#rel==2 && rel[1]!=0 && vecmax(abs(rel))<10^8,
     my(q=-rel[2]/rel[1]); if(valuation(xi2-q*V2[j],2)>PR-25, print("  xi2 = ",q," * ",NA2[j])))); }
quit
