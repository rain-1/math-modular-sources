/* b17_gamma.gp -- the Gamma-class ansatz for xi(AESZ 207) at the dominant conifold
   z_- = (349-85 sqrt17)/2^17 (exponents 0,1,1,2), over Q and over Q(sqrt17),
   at both places.  17 = 1 mod 8 so sqrt17 lies in Q_2 as well.                */
default(parisize, 6000000000);
default(realprecision, 1100);
xi = eval(Str(read("b14_xi_arch_1200.txt")))*1.0;
s17 = sqrt(17.);
zm = (349-85*s17)/131072;
Lam = log(abs(1/zm));
Z2 = zeta(2); Z3 = zeta(3); P2 = Pi^2;
tryv(nm, v, hb) = { my(r = lindep(v, 1000));
  if(type(r)=="t_COL" && r[1]!=0 && vecmax(abs(r)) < hb,
     print("  *** ", nm, " : ", r~), print("  ", nm, " : none (hb=",hb,")")); }
print("=== archimedean, 1100 digits ===");
B7 = [1, Lam, Lam^2, Lam^3, P2, P2*Lam, Z3];
tryv("Q-span of {1,L,L^2,L^3,pi^2,pi^2 L,zeta(3)}", concat([xi], B7), 10^30);
B14 = concat(B7, s17*B7);
tryv("Q(sqrt17)-span of the same", concat([xi], B14), 10^25);
tryv("Q-span + zeta(5),pi^4", concat([xi], concat(B7,[zeta(5),Pi^4])), 10^25);
B8 = concat(B7,[Pi^4]);
tryv("Q(sqrt17)-span of {..,pi^4}", concat([xi], concat(B8, s17*B8)), 10^20);

print("\n=== 2-adic ===");
read("../mum_survey/apery.gp"); read("../mum_survey/ops.gp"); read("../mum_survey/lpgen.gp");
OP=0; for(i=1,#OPS, if(OPS[i][1]=="207", OP=OPS[i]));
NN=400; pr=aperyPair(OP[4],NN); A=pr[1]; B=pr[2];
PR=1400;
xi2 = B[NN+1]/A[NN+1] + O(2^PR);
r17 = sqrt(17+O(2^(PR+5)));
print("  sqrt17 in Q_2:  ", r17+O(2^20));
al = 349+85*r17; be = 349-85*r17;
print("  v_2(alpha)=",valuation(al,2),"  v_2(beta)=",valuation(be,2),"  v_2(alpha*beta)=",valuation(al*be,2));
zm2 = be/131072; print("  v_2(z_-)=",valuation(zm2,2));
ual = al/2^valuation(al,2); ube = be/2^valuation(be,2);
La = log(ual); Lb = log(ube);
z2p = LpG(2,triv,0,2,PR);  z3p = LpG(2,triv,0,3,PR);
print("  v_2(L_2(2,1))=",valuation(z2p,2),"  v_2(L_2(3,1))=",valuation(z3p,2));
print("  v_2(log2(alpha_unit))=",valuation(La,2),"  v_2(log2(beta_unit))=",valuation(Lb,2));
t2(nm,v,hb) = { my(r=lindep(v));
  if(#r==#v && r[1]!=0 && vecmax(abs(r))<hb,
     my(ok=1); print("  *** ", nm, " : ", r~), print("  ", nm, " : none")); }
G7a = [1, La, La^2, La^3, z2p, z2p*La, z3p];
t2("Q-span {1,La,La^2,La^3,z2p,z2p*La,z3p}", concat([xi2], G7a), 10^30);
G7b = [1, Lb, Lb^2, Lb^3, z2p, z2p*Lb, z3p];
t2("Q-span with Lb", concat([xi2], G7b), 10^30);
t2("Q(sqrt17)-span with La", concat([xi2], concat(G7a, r17*G7a)), 10^22);
t2("Q(sqrt17)-span with Lb", concat([xi2], concat(G7b, r17*G7b)), 10^22);
t2("{1,La,Lb,z3p,z2p}", concat([xi2], [1,La,Lb,z3p,z2p]), 10^40);
t2("{1,La,Lb,La^2,Lb^2,La^3,Lb^3,z3p,z2p*La,z2p*Lb,z2p}",
   concat([xi2], [1,La,Lb,La^2,Lb^2,La^3,Lb^3,z3p,z2p*La,z2p*Lb,z2p]), 10^25);
quit
