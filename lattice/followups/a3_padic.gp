default(parisize, 4000000000);
read("../mum_survey/apery.gp"); read("../mum_survey/ops.gp"); read("../mum_survey/lpgen.gp");
O184=0; for(i=1,#OPS, if(OPS[i][1]=="184", O184=OPS[i]));
rowR3(a,b,c,N)={my(A=vector(N+1),B=vector(N+1)); A[1]=1;A[2]=b;B[1]=0;B[2]=1;
 for(n=1,N-1, A[n+2]=((2*n+1)*(a*n^2+a*n+b)*A[n+1]-c*n^3*A[n])/(n+1)^3;
              B[n+2]=((2*n+1)*(a*n^2+a*n+b)*B[n+1]-c*n^3*B[n])/(n+1)^3); [A,B];}
N=1050;
pr=aperyPair(O184[4],N); A4=pr[1]; B4=pr[2];
pe=rowR3(11,5,125,N);   Ae=pe[1]; Be=pe[2];
PR=3050;
chi5=[5,[1,-1,-1,1,0]];
T1 = LpG(5, triv, -2, 3, PR);   \\ L_5(3, omega^-2)
T2 = LpG(5, chi5, -2, 3, PR);   \\ L_5(3, chi5*omega^-2) = L_5(3, 1)
T3 = LpG(5, triv,  0, 3, PR);   \\ L_5(3, 1) via trivial char, no twist
T4 = LpG(5, chi5,  0, 3, PR);   \\ L_5(3, chi5)
print("v_5(T2-T3) = ", valuation(T2-T3,5), "   (both should be L_5(3,1))");
print("v_5(T1-T4) = ", valuation(T1-T4,5), "   (both should be L_5(3,chi5))");
print("T1 mod 5^12 = ", lift(T1+O(5^12)));
print("T2 mod 5^12 = ", lift(T2+O(5^12)));
{ my(Nn=1000, x4=B4[Nn+1]/A4[Nn+1]+O(5^PR), xe=Be[Nn+1]/Ae[Nn+1]+O(5^PR));
  print("Cauchy(184) = ", valuation(x4-(B4[Nn]/A4[Nn]+O(5^PR)),5));
  for(j=1,4, my(T=[T1,T2,T3,T4][j]);
    print("  target ",j,": v_5(xi184 - (1/4)T) = ", valuation(x4-T/4,5),
          "    v_5(xi_eta - (1/2)T) = ", valuation(xe-T/2,5)));
  print("  lindep(xi184, T2) = ", lindep([x4,T2]));
  print("  lindep(xi_eta,T2) = ", lindep([xe,T2]));
}
quit
