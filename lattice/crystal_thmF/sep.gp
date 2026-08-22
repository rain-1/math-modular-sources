/* sep.gp -- clean p-adic separation census.  mu1 = v_p(a_N)/N, muC = v_p(Cas_N)/N,
   sigma_p = muC - 2 mu1  (exact rates up to O(log N / N)).  Also the slope-prime
   tower eigenvalue (should be 1) for Zudilin's row at p=2. */
default(parisizemax, 6000000000);
read("lattice/euler_criterion/rows.gp");
N = 1500;
{ ROWS = [
 ["A     (7,2,-8)  R2", rowR2(7,2,-8,N), -8],
 ["B     (9,3,27)  R2", rowR2(9,3,27,N), 27],
 ["C     (10,3,9)  R2", rowR2(10,3,9,N), 9],
 ["D z2  (11,3,-1) R2", rowR2(11,3,-1,N), -1],
 ["E     (12,4,32) R2", rowR2(12,4,32,N), 32],
 ["F     (17,6,72) R2", rowR2(17,6,72,N), 72],
 ["gamma z3        R3", rowR3(17,5,1,N), 1],
 ["alpha Domb      R3", rowR3(10,4,64,N), 64],
 ["delta (7,3,81)  R3", rowR3(7,3,81,N), 81],
 ["eps   (12,4,16) R3", rowR3(12,4,16,N), 16],
 ["zeta  (9,3,-27) R3", rowR3(9,3,-27,N), -27],
 ["eta   (11,5,125)R3", rowR3(11,5,125,N), 125],
 ["s18 Cooper      R3", rowS18(N), 192],
 ["cusp L(f,2) N=12  ", rowCusp(N), 64]
]; }
PR = [2,3,5,7,11,13];
vp(x,p) = if(x==0, 0, valuation(x,p));
print("row                 |  p | mu1=v(a_N)/N | sigma_p (meas) | v_p(c) | agree?");
{ for(r=1,#ROWS,
  my(nm=ROWS[r][1], A=ROWS[r][2][1], B=ROWS[r][2][2], c=ROWS[r][3], Cas=vector(N));
  for(n=1,N-1, Cas[n]=A[n+2]*B[n+1]-A[n+1]*B[n+2]);
  for(j=1,#PR, my(p=PR[j]);
    my(m1=vp(A[N+1],p)/N, mC=vp(Cas[N-1],p)/(N-1), sg=mC-2*m1, pc=valuation(c,p));
    print(nm," | ",p," | ",1.0*m1," | ",1.0*sg," | ",pc," | ",
          if(abs(1.0*sg-pc)<0.05,"YES","** NO **")));); }

print(); print("=== Zudilin Catalan row at the slope prime p=2: tower eigenvalues ===");
MZ = 2200; Z = rowZud(MZ);
{ for(a=1,3,
   my(s=0,L=[],LB=[]);
   while(a*2^(s+1)<=MZ,
     my(n0=a*2^s,n1=a*2^(s+1));
     L=concat(L,[valuation(Z[1][n1+1]/Z[1][n0+1]-1,2)]);
     LB=concat(LB,[valuation(Z[2][n1+1]/Z[2][n0+1]-1,2)]); s++);
   print(" a=",a,"  v(rhoA-1)=",L); print("        v(rhoB-1)=",LB)); }
print();
print("=== Zudilin: v_2(P_m/Q_m - xi) vs predicted 8m-1-4 s_2(m) ===");
{ my(xi = Z[2][2000+1]/Z[1][2000+1]);
  for(i=1,10, my(m=[5,10,17,32,33,64,100,127,128,200][i]);
    print("  m=",m,"  v_2(P/Q - xi_ref)=",valuation(Z[2][m+1]/Z[1][m+1]-xi,2),
          "   8m-1-4s_2(m)=",8*m-1-4*hammingweight(m))); }
