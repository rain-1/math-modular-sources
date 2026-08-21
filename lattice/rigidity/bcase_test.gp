\\ B-case: expand the source Phi_B on C's curve (t_C, F_C) and measure the 3-adic limit.
PS=70; read("cover.gp");
chi(n)=kronecker(-3,n);
S = sum(m=1,PREC-1, sumdiv(m,d, chi(m/d)*d^2)*q^m) + O(q^PREC);
PhiBs = S - 6*Vd(S,2) - 8*Vd(S,4);
print("Phi_B(eta) == (1-6V2-8V4)S : ", PhiB-PhiBs == O(q^PREC));
Psi0 = S - 4*Vd(S,2);
ThB = Dqinv(PhiBs,2); ThC = Dqinv(PhiC,2); ThP = Dqinv(2*Vd(Psi0,2),2);
print("Theta_B = Theta_C + D^-2(2 V2 Psi0) : ", ThB-ThC-ThP == O(q^PREC));
rev = serreverse(tC);
Aser = subst(FC,q,rev);
Bser = subst(FC*ThC,q,rev);
Bbser= subst(FC*ThB,q,rev);          \\ companion of Phi_B ON C's CURVE
Pser = subst(FC*ThP,q,rev);          \\ companion of 2V2(1-4V2)S on C's curve
\\ compare with the genuine Zagier rows
N=65; rows(a,b,c)={my(A=vector(N+1),B=vector(N+1));A[1]=1;A[2]=b;B[1]=0;B[2]=1;
 for(n=1,N-1,A[n+2]=((a*n^2+a*n+b)*A[n+1]-c*n^2*A[n])/(n+1)^2;B[n+2]=((a*n^2+a*n+b)*B[n+1]-c*n^2*B[n])/(n+1)^2);[A,B]};
RC=rows(10,3,9); RB=rows(9,3,27);
print("A^C from eta == Zagier row C a_n (n<=90): ", vector(61,i,polcoeff(Aser,i-1))==vector(61,i,RC[1][i]));
print("B^C from eta == Zagier row C b_n (n<=90): ", vector(61,i,polcoeff(Bser,i-1))==vector(61,i,RC[2][i]));
xiC = RC[2][N+1]/RC[1][N+1]; xiB = RB[2][N+1]/RB[1][N+1];
print("n   v3(psi_n/a_n)   v3(bB_n/a_n - xi_C)   [psi = 2V2(1-4V2)S companion on C-curve]");
for(k=1,6, n=10*k; a=polcoeff(Aser,n); d=polcoeff(Bbser,n)/a-xiC; print(n,"   ",valuation(polcoeff(Pser,n)/a,3),"   ",if(d==0,"ZERO",valuation(d,3))));
