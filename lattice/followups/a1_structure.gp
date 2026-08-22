/* a1_structure.gp -- AESZ 184 basic structure: Hadamard relation to eta=(11,5,125),
   characteristic roots, integrality, denominator exponent, kappa_5, slopes.        */
default(parisize, 2000000000);
read("../mum_survey/apery.gp");
read("../mum_survey/ops.gp");
default(realprecision, 120);

/* locate AESZ 184 */
O184 = 0;
for(i=1,#OPS, if(OPS[i][1]=="184", O184=OPS[i]));
print("AESZ184 nn=", O184[2], " degz=", O184[3]);
print("  P_i(X) = ", O184[4]);
print("  disc/sing = ", O184[7]);

N = 400;
pr = aperyPair(O184[4], N);  A4 = pr[1]; B4 = pr[2];

/* eta row (11,5,125), order-3 Zagier-type: (n+1)^3 u_{n+1}=(2n+1)(a n^2+a n+b)u_n - c n^3 u_{n-1} */
rowR3(a,b,c,N) =
{ my(A=vector(N+1), B=vector(N+1));
  A[1]=1; A[2]=b; B[1]=0; B[2]=1;
  for(n=1,N-1,
    A[n+2] = ((2*n+1)*(a*n^2+a*n+b)*A[n+1] - c*n^3*A[n])/(n+1)^3;
    B[n+2] = ((2*n+1)*(a*n^2+a*n+b)*B[n+1] - c*n^3*B[n])/(n+1)^3);
  [A,B];
}
pe = rowR3(11,5,125,N); Ae = pe[1]; Be = pe[2];

print("\n--- A_n of AESZ184 vs binom(2n,n)*A_n(eta) ---");
bad = 0;
for(n=0,N, if(A4[n+1] != binomial(2*n,n)*Ae[n+1], bad++; if(bad<4,print("  MISMATCH n=",n))));
print("  mismatches among n=0..",N,": ", bad);

print("\n--- B_n of AESZ184 vs binom(2n,n)*B_n(eta) ---");
bad2 = 0; for(n=0,N, if(B4[n+1] != binomial(2*n,n)*Be[n+1], bad2++));
print("  mismatches: ", bad2);
print("  first 8 B(184): ", vector(8,i,B4[i]));
print("  first 8 binom*B(eta): ", vector(8,i,binomial(2*(i-1),i-1)*Be[i]));

print("\n--- integrality ---");
print("  A_n(184) all integers? ", sum(n=0,N, denominator(A4[n+1])!=1)==0);
print("  A_n(eta) all integers? ", sum(n=0,N, denominator(Ae[n+1])!=1)==0);
print("  kappa_5(184) = ", -valuation(A4[N+1],5)*1.0/N, "   v_5(A_N)=",valuation(A4[N+1],5));
print("  kappa_2(184) = ", -valuation(A4[N+1],2)*1.0/N);
print("  max v_5(A_n(184)), n<=400: ", vecmax(vector(N,n,valuation(A4[n+1],5))));
print("  max v_5(A_n(eta)),  n<=400: ", vecmax(vector(N,n,valuation(Ae[n+1],5))));

print("\n--- denominator exponent k (minimal with d_n^k B_n in Z) ---");
print("  k(184) up to n=200: ", denomExp(B4, 200));
print("  k(eta)  up to n=200: ", denomExp(Be, 200));
/* sharpness: is k-1 enough? */
{ for(kk=1,4,
  my(ok=1, dn=1);
  for(n=1,200, dn=lcm(dn,n); if(denominator(B4[n+1]*dn^kk)!=1, ok=0; break));
  print("   d_n^",kk," B_n(184) integral: ", ok));
}
{ for(kk=1,4,
  my(ok=1, dn=1);
  for(n=1,200, dn=lcm(dn,n); if(denominator(Be[n+1]*dn^kk)!=1, ok=0; break));
  print("   d_n^",kk," B_n(eta) integral: ", ok));
}

print("\n--- characteristic roots ---");
print("  184: x^2 - 88x + 2000 roots = ", polroots(x^2-88*x+2000));
print("       |rho| = ", abs(polroots(x^2-88*x+2000)[1]), "  = 20*sqrt5 = ", 20*sqrt(5));
print("  eta: x^2 - 22x + 125 roots = ", polroots(x^2-22*x+125));
print("       |rho| = ", abs(polroots(x^2-22*x+125)[1]), " = sqrt125 = ",sqrt(125));
print("  c(184) = prod roots = 2000 = 2^4*5^3 ; v_5=",valuation(2000,5)," v_2=",valuation(2000,2));

print("\n--- measured lambda1 = |A_n|^{1/n}, and the oscillation ---");
for(j=1,8, my(n=50*j); print("  n=",n,"  |A_n|^(1/n) = ", exp(log(abs(1.0*A4[n+1]))/n)));
print("\n--- B_n/A_n (archimedean): does it converge? ---");
for(j=1,8, my(n=50*j); print("  n=",n,"  B_n/A_n = ", 1.0*B4[n+1]/A4[n+1]));
print("  eta:");
for(j=1,8, my(n=50*j); print("  n=",n,"  B_n/A_n = ", 1.0*Be[n+1]/Ae[n+1]));

print("\n--- Casoratian valuations (w = lim v_p(Cas_n)/n) ---");
{ for(j=1,8, my(n=50*j, C4=A4[n+1]*B4[n]-A4[n]*B4[n+1], Ce=Ae[n+1]*Be[n]-Ae[n]*Be[n+1]);
    print("  n=",n," v_5(Cas184)=",valuation(C4,5)," rate=",valuation(C4,5)*1.0/n,
          "   v_2(Cas184)=",valuation(C4,2),
          "   v_5(Cas_eta)=",valuation(Ce,5)," rate=",valuation(Ce,5)*1.0/n)); }
quit
