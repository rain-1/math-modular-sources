/* a2_align.gp -- exact Hadamard identity, 5-adic limits, cross-determinant. */
default(parisize, 4000000000);
read("../mum_survey/apery.gp"); read("../mum_survey/ops.gp"); read("../mum_survey/lpgen.gp");
default(realprecision, 60);
O184=0; for(i=1,#OPS, if(OPS[i][1]=="184", O184=OPS[i]));
rowR3(a,b,c,N)={my(A=vector(N+1),B=vector(N+1)); A[1]=1;A[2]=b;B[1]=0;B[2]=1;
 for(n=1,N-1, A[n+2]=((2*n+1)*(a*n^2+a*n+b)*A[n+1]-c*n^3*A[n])/(n+1)^3;
              B[n+2]=((2*n+1)*(a*n^2+a*n+b)*B[n+1]-c*n^3*B[n])/(n+1)^3); [A,B];}
N=1200;
pr=aperyPair(O184[4],N); A4=pr[1]; B4=pr[2];
pe=rowR3(11,5,125,N);   Ae=pe[1]; Be=pe[2];

print("=== exact Hadamard identities, n=0..",N," ===");
print("  #{n: A^184_n != binom(2n,n) A^eta_n} = ", sum(n=0,N, A4[n+1]!=binomial(2*n,n)*Ae[n+1]));
print("  #{n: B^184_n != (1/2) binom(2n,n) B^eta_n} = ", sum(n=0,N, B4[n+1]!=binomial(2*n,n)*Be[n+1]/2));
print("  => B^184_n/A^184_n = (1/2) B^eta_n/A^eta_n  identically.");
print("  cross-determinant  B^eta_n A^184_m - 2 A^eta_n B^184_m  at n=m: max |.| over n<=",N," = ",
      vecmax(vector(N+1,i, abs(Be[i]*A4[i]-2*Ae[i]*B4[i]))));

print("\n=== 5-adic limits ===");
PR = 3000;
z5 = LpG(5, triv, -2, 3, PR);   \\ zeta_5(3) = L_5(3, omega^-2)
L5chi = LpG(5, [5,[1,-1,-1,1,0]], -2, 3, PR); \\ L_5(3, chi5 omega^-2) -- chi5 quadratic mod 5
print("  v_5( L_5(3,chi5 om^-2) - zeta_5(3) ) = ", valuation(L5chi-z5,5), "   (chi5*om^-2 = trivial at p=5)");
{ my(Nn=1000);
  my(x4 = B4[Nn+1]/A4[Nn+1] + O(5^PR), xe = Be[Nn+1]/Ae[Nn+1] + O(5^PR));
  my(c4 = B4[Nn]/A4[Nn] + O(5^PR),  ce = Be[Nn]/Ae[Nn] + O(5^PR));
  print("  N=",Nn);
  print("  Cauchy precision  184: v_5(xi_N - xi_{N-1}) = ", valuation(x4-c4,5));
  print("  Cauchy precision  eta: v_5(xi_N - xi_{N-1}) = ", valuation(xe-ce,5));
  print("  v_5( xi_184 - (1/4) zeta_5(3) ) = ", valuation(x4 - z5/4, 5));
  print("  v_5( xi_eta - (1/2) zeta_5(3) ) = ", valuation(xe - z5/2, 5));
  print("  v_5( xi_eta - 2 xi_184 )        = ", valuation(xe-2*x4, 5), "  (exact identity => oo)");
}

print("\n=== cross-determinant on a grid (n,m), D = B^eta_n A^184_m - 2 A^eta_n B^184_m ===");
{ my(ns=[100,200,300,400,600,800,1000]);
  for(i=1,#ns, for(j=1,#ns,
    my(n=ns[i], m=ns[j], D = Be[n+1]*A4[m+1] - 2*Ae[n+1]*B4[m+1]);
    printf("  n=%4d m=%4d  v_5(D)=%s\n", n, m, if(D==0,"+oo (D=0)",Str(valuation(D,5))))));
}
quit
