default(parisizemax, 6000000000);
LOG = "/home/ubuntu/code/math-modular-sources/lattice/multislope/row2_kl.log";
W(s) = write(LOG, s);
read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row2_data.txt");
rB = R2X1/R2A; rC = R2X2/R2A; rD = R2X3/R2A;
KTR = 400; PRE = 800;
om(a) = if(a%4==1, 1, -1);
tw(a) = a/om(a);
/* Washington Thm 5.11, chi = 1, p = 2, F = 4 : zeta_2(s) = L_2(s,1) */
z2(s) = { (1/(4*(s-1))) * sum(i=1,2, my(a=[1,3][i]); (tw(a)+O(2^PRE))^(1-s) * sum(j=0, KTR, binomial(1-s,j)*((4/a)+O(2^PRE))^j*bernfrac(j))) };
/* validation */
{ my(ok=1);
  for(m=1,6, my(nn=2*m, lhs, rhs);
    lhs = (1/(4*((1-nn)-1)))*sum(i=1,2, my(a=[1,3][i]); tw(a)^(1-(1-nn))*sum(j=0,nn+2, binomial(1-(1-nn),j)*(4/a)^j*bernfrac(j)));
    rhs = -(1-2^(nn-1))*bernfrac(nn)/nn;
    if(lhs != rhs, ok = 0));
  W(Str("=== KL zeta_2 validation (exact interpolation L_2(1-n,1) = -(1-2^{n-1})B_n/n): ", if(ok,"PASS","FAIL"), " ===")); }
ZV = vector(6, i, z2(i+1));   /* zeta_2(2)..zeta_2(7) */
{ for(i=1,6, W(Str("  v_2(zeta_2(", i+1, ")) = ", valuation(ZV[i],2)))); }
NM = ["B","C","D"]; RV = [rB,rC,rD];
W("");
W("--- is xi_2(X) a rational multiple of zeta_2(m)?  (lindep at 700 2-adic digits) ---");
{ for(s=1,3, my(x = RV[s] + O(2^700));
    for(i=1,6,
      my(L = lindep([x, ZV[i]+O(2^700)]), h = max(#Str(abs(L[1])), #Str(abs(L[2]))));
      W(Str("  ", NM[s], " vs zeta_2(", i+1, "):  coeff sizes ", h, " digits", if(h<25, Str("  -> ", L), "")))));
}
W("");
W("--- xi_2(X) in Q + Q*zeta_2(5)?  and in Q*zeta_2(3)+Q*zeta_2(5)? ---");
{ for(s=1,3, my(x = RV[s] + O(2^700), L1, L2);
    L1 = lindep([x, 1+O(2^700), ZV[4]+O(2^700)]);
    L2 = lindep([x, ZV[2]+O(2^700), ZV[4]+O(2^700)]);
    W(Str("  ", NM[s], ": [xi,1,z2(5)] coeff sizes ", vector(3,i,#Str(abs(L1[i]))), "   [xi,z2(3),z2(5)] ", vector(3,i,#Str(abs(L2[i])))))); }
W("");
W("--- is xi_2 rational? algebraic of degree <= 4? ---");
{ for(s=1,3, my(x = RV[s]+O(2^900), L);
    W(Str("  ", NM[s], ": bestappr-style lindep([xi,1]) coeff sizes ", vector(2,i,#Str(abs(lindep([x,1+O(2^900)])[i]))),
          "   algdep deg 4 height ", #Str(vecmax(abs(Vec(algdep(x,4))))))) ); }
W("DONE"); quit;
