/* 08_census_check.gp -- alignment-prime census cross-check (SLOPE_CENSUS.md) */
N = 600;
{ zpair(aa,bb,dd,NM) = my(u=vector(NM+2), w=vector(NM+2)); u[1]=1; u[2]=bb; w[1]=0; w[2]=1;
  for(n=1,NM, my(Pn=aa*(n^2+n)+bb, Qn=dd*n^2);
      u[n+2]=(Pn*u[n+1]-Qn*u[n])/(n+1)^2; w[n+2]=(Pn*w[n+1]-Qn*w[n])/(n+1)^2);
  [vector(NM+1,k,u[k]), vector(NM+1,k,w[k])]; }
/* Cooper s_18 : (n+1)^3 A_{n+1} = 2(2n+1)(7n^2+7n+3) A_n - 12 n(16n^2-1) A_{n-1} */
{ s18(NM) = my(u=vector(NM+2), w=vector(NM+2)); u[1]=1; u[2]=6; w[1]=0; w[2]=1;
  for(n=1,NM, my(Pn=2*(2*n+1)*(7*n^2+7*n+3), Qn=12*n*(16*n^2-1));
      u[n+2]=(Pn*u[n+1]-Qn*u[n])/(n+1)^3; w[n+2]=(Pn*w[n+1]-Qn*w[n])/(n+1)^3);
  [vector(NM+1,k,u[k]), vector(NM+1,k,w[k])]; }
RB=zpair(9,3,27,N); RC=zpair(10,3,9,N); RF=zpair(17,6,72,N); RE=zpair(12,4,32,N); RD=zpair(11,3,-1,N); RA=zpair(7,2,-8,N);
S18=s18(N);
xB=RB[2][N+1]/RB[1][N+1]; xC=RC[2][N+1]/RC[1][N+1]; xF=RF[2][N+1]/RF[1][N+1];
xE=RE[2][N+1]/RE[1][N+1]; xD=RD[2][N+1]/RD[1][N+1]; xA=RA[2][N+1]/RA[1][N+1];
x18=S18[2][N+1]/S18[1][N+1];
print("N = ", N);
print("=== alignment primes sigma_p = v_p(d) for Zagier's six ===");
nm=["A","B","C","D","E","F"]; dd=[-8,27,9,-1,32,72];
for(j=1,6, print("  ",nm[j],"  d=",dd[j],"   sigma_2=",valuation(dd[j],2),"  sigma_3=",valuation(dd[j],3),"  sigma_5=",valuation(dd[j],5)));
print();
print("=== p=3 cluster:  B, C, F, and Cooper s_18 ===");
print("v_3(xi_B - xi_C)        = ", valuation(xB-xC,3));
print("v_3(4 xi_F - 5 xi_C)    = ", valuation(4*xF-5*xC,3));
print("v_3(xi_s18 - xi_C)      = ", valuation(x18-xC,3));
print("v_3(4 xi_F - 5 xi_s18)  = ", valuation(4*xF-5*x18,3));
print("s_18 Cauchy: v_3(xi_s18(N)-xi_s18(N-1)) = ", valuation(x18 - S18[2][N]/S18[1][N],3), "  [N = ",N,"]");
print();
print("=== p=2 cluster:  A, E, F ===");
print("v_2(xi_A)               = ", valuation(xA,2), "   (predicted 0, at rate 3n)");
print("v_2(xi_E)               = ", valuation(xE,2));
print("v_2(xi_F)               = ", valuation(xF,2));
print("v_2(xi_E - xi_F)        = ", valuation(xE-xF,2), "   (different L-values: zeta_2(2) vs L_2(2,chi_12))");
print("v_2(xi_A - xi_E)        = ", valuation(xA-xE,2));
print();
print("=== D: no alignment prime at all ===");
print("v_2(xi_D(N)-xi_D(N-1))  = ", valuation(xD - RD[2][N]/RD[1][N],2));
print("v_3(xi_D(N)-xi_D(N-1))  = ", valuation(xD - RD[2][N]/RD[1][N],3));
print("v_5(xi_D(N)-xi_D(N-1))  = ", valuation(xD - RD[2][N]/RD[1][N],5));
quit;
