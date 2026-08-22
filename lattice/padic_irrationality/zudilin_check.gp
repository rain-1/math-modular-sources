\\ zudilin_check.gp -- the kappa-corrected criterion on a NON-integral row.
\\ Zudilin's Catalan row (arXiv math/0201024 Thm 1), recurrence as in
\\ consolidation/ZUDILIN_2ADIC.md (1.3).  Measures kappa_2 = -v_2(Q_m)/m,
\\ the quality slope sigma_2 = v_2(xi_2 - P_m/Q_m)/m, the archimedean rate
\\ log|Q_m|/m, and the prime-to-2 clearing exponent.
default(parisizemax,"6G");
default(realprecision,40);
MZ = 420;

pz(n) = 20*n^2 - 8*n + 1;
qz(n) = 3520*n^6 + 5632*n^5 + 2064*n^4 - 384*n^3 - 156*n^2 + 16*n + 7;

{ zudrow(MZ) =
  my(Q=vector(MZ+1), P=vector(MZ+1));
  Q[1]=1; Q[2]=7/4; P[1]=0; P[2]=13/8;
  for(n=1,MZ-1,
    Q[n+2] = (qz(n)*Q[n+1] + (2*n-1)^2*(2*n)^2*pz(n+1)*Q[n]) / ((2*n+1)^2*(2*n+2)^2*pz(n));
    P[n+2] = (qz(n)*P[n+1] + (2*n-1)^2*(2*n)^2*pz(n+1)*P[n]) / ((2*n+1)^2*(2*n+2)^2*pz(n));
  );
  [Q,P];
}
{
my(R=zudrow(MZ), Q=R[1], P=R[2], xi=P[MZ+1]/Q[MZ+1]);
print("Q_0..Q_5 = ", vector(6,j,Q[j]));
print("P_0..P_5 = ", vector(6,j,P[j]));
print("\n m | v_2(Q_m)  (theory -4m+2s_2(m)) | v_2(xi-P_m/Q_m) (theory 8m-1-4s_2(m)) | log|Q_m|/m (-> 5 log phi = ",5*log((1+sqrt(5))/2),")");
for(j=1,6, my(m=[64,100,128,200,300,400][j]);
  printf(" %3d | %6d  (%6d) | %6d  (%6d) | %.6f\n",
    m, valuation(Q[m+1],2), -4*m+2*hammingweight(m),
    valuation(xi-P[m+1]/Q[m+1],2), 8*m-1-4*hammingweight(m),
    log(abs(Q[m+1]*1.))/m));
print("\nprime-to-2 denominator of P_m: needs lcm(1..2m-1)^2 ?");
for(j=1,4, my(m=[20,40,60,80][j], D=lcm(vector(2*m-1,i,i)), x=P[m+1]*2^(-valuation(P[m+1],2)));
  printf(" m=%3d : D_{2m-1}^2 P_m integral (up to 2-part)? %d ;  log(D_{2m-1}^2)/m = %.5f\n",
     m, denominator(x*D^2)==1, 2*log(D*1.)/m));
print("\n=> kappa_2 = 4, sigma_2 = 8 = v_2(c) + 2 kappa_2 with v_2(c)=0,");
print("   k = 4 nats/index, log lambda_1 = 5 log phi = ", 5*log((1+sqrt(5))/2));
my(phi=(1+sqrt(5))/2, sg=8, kap=4, kk=4, l1=5*log(phi));
printf("   S_2 = %.6f   theta_2 = %.6f\n", sg*log(2)-kk-kap*log(2)-l1, sg*log(2)/(kk+kap*log(2)+l1));
}
quit;
