/* ============================================================
   05_padic.gp -- p-adic Apery limits of the cusp-move orbit {A, C, F}
   ============================================================ */
N = 800;
{ pair(aa,bb,dd,NM) =
  my(u=vector(NM+2), w=vector(NM+2)); u[1]=1; u[2]=bb; w[1]=0; w[2]=1;
  for(n=1,NM, my(Pn=aa*(n^2+n)+bb, Qn=dd*n^2);
      u[n+2]=(Pn*u[n+1]-Qn*u[n])/(n+1)^2;
      w[n+2]=(Pn*w[n+1]-Qn*w[n])/(n+1)^2);
  [vector(NM+1,k,u[k]), vector(NM+1,k,w[k])];
}
RA = pair(7,2,-8,N); RC = pair(10,3,9,N); RF = pair(17,6,72,N); RFp = pair(-17,-6,72,N);
xA = RA[2][N+1]/RA[1][N+1];
xC = RC[2][N+1]/RC[1][N+1];
xF = RF[2][N+1]/RF[1][N+1];
xFp= RFp[2][N+1]/RFp[1][N+1];
xAm= RA[2][N]/RA[1][N]; xCm= RC[2][N]/RC[1][N]; xFm= RF[2][N]/RF[1][N];

print("N = ", N);
print();
print("=== A. slopes sigma_p = v_p(d) and Cauchy precision actually attained ===");
print("row A  d=-8 : sigma_2=", valuation(-8,2), "  sigma_3=", valuation(-8,3));
print("row C  d= 9 : sigma_2=", valuation(9,2),  "  sigma_3=", valuation(9,3));
print("row F  d=72 : sigma_2=", valuation(72,2), "  sigma_3=", valuation(72,3));
print("v_3(xi_C(N)-xi_C(N-1)) = ", valuation(xC-xCm,3), "   [2N = ", 2*N, "]");
print("v_3(xi_F(N)-xi_F(N-1)) = ", valuation(xF-xFm,3), "   [2N = ", 2*N, "]");
print("v_2(xi_A(N)-xi_A(N-1)) = ", valuation(xA-xAm,2), "   [3N = ", 3*N, "]");
print("v_3(xi_A(N)-xi_A(N-1)) = ", valuation(xA-xAm,3), "   [no 3-adic slope]");
print();
print("=== B. the predicted cusp-move ratio  4 xi_3^F = 5 xi_3^C ===");
print("v_3(xi_C) = ", valuation(xC,3), ",  v_3(xi_F) = ", valuation(xF,3));
print("v_3(4 xi_F - 5 xi_C)      = ", valuation(4*xF-5*xC,3));
print("v_3(4 xi_F' + 5 xi_C)     = ", valuation(4*xFp+5*xC,3));
print("v_3(xi_F' + xi_F)         = ", valuation(xFp+xF,3), "   (exact: F' row = (-1)^n twist)");
print("v_3(xi_F - xi_C)          = ", valuation(xF-xC,3), "   (control: should be small)");
print("v_3(3 xi_F - 4 xi_C)      = ", valuation(3*xF-4*xC,3), "   (control)");
print();
print("=== C. 3-adic digits of xi_3^C, xi_3^F ===");
{ p3(z,pp,k) = my(v=valuation(z,pp), y=z/pp^v);
  Str(pp,"^",v," * ", lift(Mod(numerator(y),pp^k)*Mod(denominator(y),pp^k)^-1)); }
print("xi_3^C mod 3^40 : ", p3(xC,3,40));
print("xi_3^F mod 3^40 : ", p3(xF,3,40));
print("(5/4)xi_3^C     : ", p3(5*xC/4,3,40));
print("  [CONJ_D_PROOF records xi_3^C = 3^-1 * 8386265965554334030 mod 3^39]");
print();
print("=== D. row A at p = 2 : the limit is ZERO ===");
for(j=1,8, my(n=100*j); print("  n=", n, "  v_2(b^A_n/a^A_n) = ", valuation(RA[2][n+1]/RA[1][n+1],2), "   (3n = ", 3*n, ")"));
print();
print("=== E. row F at p = 2 (second alignment prime of F; A's partner would be here) ===");
for(j=1,4, my(n=200*j); print("  n=", n, "  v_2(b^F_n/a^F_n) = ", valuation(RF[2][n+1]/RF[1][n+1],2)));
print("v_2(xi_F(N)-xi_F(N-1)) = ", valuation(xF-xFm,2), "   [3N = ", 3*N, "]");
print("v_2(xi_A - xi_F)  = ", valuation(xA-xF,2), " ; v_2(xi_A) = ", valuation(xA,2), " ; v_2(xi_F) = ", valuation(xF,2));
quit;
