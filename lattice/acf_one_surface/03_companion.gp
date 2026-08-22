/* ============================================================
   03_companion.gp -- how the companion B = F * D^{-2}Phi moves
   under the cusp change, and the transfer identity.
   ============================================================ */
read("lib.gp");
E1 = etaprod(1); E2 = etaprod(2); E3 = etaprod(3); E6 = etaprod(6);
tC = 'q * E1^4 * E6^8 / (E2^8 * E3^4);
FC = E2^6 * E3 / (E1^3 * E6^2);
tA  = tC/(1-tC);      FA  = (1-tC)*FC;
tFp = tC/(1-9*tC);    FFp = (1-9*tC)*FC;

Es = Ess(); Et = Ett();
PhiC  = Es - 8*Vop(2,Es);
PhiFp = Es +   Vop(2,Es);
PhiA  = Et -   Vop(2,Et);
Psi0  = Es - 4*Vop(2,Es);

ThC = Dinv2(PhiC); ThFp = Dinv2(PhiFp); ThA = Dinv2(PhiA); ThPsi = Dinv2(Psi0);

BC  = FC  * ThC;      /* = sum b_n^C  t_C^n   */
BFp = FFp * ThFp;     /* = sum b_n^F' t_F'^n  */
BA  = FA  * ThA;
XI  = FC  * ThPsi;    /* the Psi_0-companion on the C-curve */

print("=== A. transfer identity for the companion (q-series, exact) ===");
print("B^{F'} + (5/4)(1-9t_C) B^C - (9/4)(1-9t_C) XI : ", BFp + 5*(1-9*tC)*BC/4 - 9*(1-9*tC)*XI/4);
print("A^{F'} - (1-9t_C) A^C                          : ", FFp - (1-9*tC)*FC);
print();
print("=== B. t-expansions: rows from the modular construction ===");
NN = 60;
cf(G,t) = my(H = texp(G,t)); vector(NN+1, k, polcoef(H, k-1, 'q));
aC = cf(FC ,tC ); bC = cf(BC ,tC );
aFp= cf(FFp,tFp); bFp= cf(BFp,tFp);
aA = cf(FA ,tA ); bA = cf(BA ,tA );
psn= cf(XI ,tC );
gc(v,n) = v[n+1];
print("a^C  : ", vector(6,k,gc(aC ,k-1)));
print("b^C  : ", vector(6,k,gc(bC ,k-1)));
print("a^F' : ", vector(6,k,gc(aFp,k-1)));
print("b^F' : ", vector(6,k,gc(bFp,k-1)));
print("a^A  : ", vector(6,k,gc(aA ,k-1)));
print("b^A  : ", vector(6,k,gc(bA ,k-1)));
print("psi_n: ", vector(6,k,gc(psn,k-1)));

/* recurrence rows, second solution normalised u_0=0,u_1=1 */
{ pair(aa,bb,dd,N) =
  my(u=vector(N+2), w=vector(N+2)); u[1]=1; u[2]=bb; w[1]=0; w[2]=1;
  for(n=1,N, my(Pn=aa*(n^2+n)+bb, Qn=dd*n^2);
      u[n+2]=(Pn*u[n+1]-Qn*u[n])/(n+1)^2;
      w[n+2]=(Pn*w[n+1]-Qn*w[n])/(n+1)^2);
  [vector(N+1,k,u[k]), vector(N+1,k,w[k])];
}
RC = pair(10,3,9,NN); RA = pair(7,2,-8,NN); RFp = pair(-17,-6,72,NN); RF = pair(17,6,72,NN);
er(v,w) = my(m=min(#v,#w)); vecmax(vector(m,k,abs(v[k]-w[k])));
print();
print("max |a^C - recurrence|  = ", er(vector(NN+1,k,gc(aC ,k-1)), RC[1]));
print("max |b^C - recurrence|  = ", er(vector(NN+1,k,gc(bC ,k-1)), RC[2]));
print("max |a^A - recurrence|  = ", er(vector(NN+1,k,gc(aA ,k-1)), RA[1]));
print("max |b^A - recurrence|  = ", er(vector(NN+1,k,gc(bA ,k-1)), RA[2]));
print("max |a^F'- recurrence|  = ", er(vector(NN+1,k,gc(aFp,k-1)), RFp[1]));
print("max |b^F'- recurrence|  = ", er(vector(NN+1,k,gc(bFp,k-1)), RFp[2]));
print("max |a^F' - (-1)^n a^F| = ", er(RFp[1], vector(NN+1,k,(-1)^(k-1)*RF[1][k])));
print("max |b^F' + (-1)^n b^F| = ", er(RFp[2], vector(NN+1,k,-(-1)^(k-1)*RF[2][k])));
print();
print("=== C. 3-adic size of the Psi_0-companion on the C-curve (P(1/4)=0 class) ===");
for(j=1,6, my(n=10*j); print("  n=",n,"  v3(psi_n/a^C_n) = ", valuation(gc(psn,n),3)-valuation(RC[1][n+1],3), "   (2n = ",2*n,")"));
print();
print("=== D. Theta_A is NOT in the span of Theta_C, Theta_{Psi_0} ===");
Ecal = Dinv2(Es); EcalT = Dinv2(Et);
mat4 = matrix(8,4,i,j, my(L=[Ecal, Vop(2,Ecal), EcalT, Vop(2,EcalT)]); polcoef(L[j],i,'q));
print("rank of {Ecal, V2 Ecal, EcalT, V2 EcalT} on q^1..q^8 : ", matrank(mat4));
mat3 = matrix(8,3,i,j, my(L=[ThC, ThPsi, ThA]); polcoef(L[j],i,'q));
print("rank of {Th_C, Th_Psi0, Th_A}                        : ", matrank(mat3));
quit;
