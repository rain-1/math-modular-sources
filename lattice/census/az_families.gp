\\ AZ complex-root third-order rows delta=(7,3,81) [Beauville-IV, symmetric square of Zagier A]
\\ and eta=(11,5,125) [carries L(3,chi_5)/2, book v8 08c_chi5_application.tex].
\\ Recurrence form (from lattice/third_order.gp): (n+1)^3 u_{n+1}=(2n+1)(a n^2+an+b)u_n - c n^3 u_{n-1}.
default(realprecision,400);
N=800;
rows(a,b,c)={my(A=vector(N+1),B=vector(N+1));A[1]=1;A[2]=b;B[1]=0;B[2]=1;
 for(n=1,N-1,A[n+2]=((2*n+1)*(a*n^2+a*n+b)*A[n+1]-c*n^3*A[n])/(n+1)^3;
            B[n+2]=((2*n+1)*(a*n^2+a*n+b)*B[n+1]-c*n^3*B[n])/(n+1)^3);[A,B]};

Rdelta=rows(7,3,81); \\ delta
Reta=rows(11,5,125); \\ eta
RDomb=rows(10,4,64);
RT=rows(12,4,16);
RApery=rows(17,5,1);

print("=== characteristic roots (complex) ===");
print("delta(7,3,81): ",polroots(x^2-14*x+81)~);
print("eta(11,5,125): ",polroots(x^2-22*x+125)~);

print("\n=== 3-adic limit xi_3 of delta (from b_n/a_n, p-adic precision from slope~4n) ===");
{my(A=Rdelta[1],B=Rdelta[2]);
 my(x1=B[N+1]/A[N+1]);
 print("delta: 3-adic digits of b_n/a_n at n=",N,", precision ~",4*N," (padicprec cap; printing valuation of successive diffs)");
 for(k=1,6,my(n=100*k);print("  n=",n," v_3(b_n/a_n - b_{n-1}/a_{n-1})=",valuation(B[n+1]/A[n+1]-B[n]/A[n],3)));}

print("\n=== 5-adic limit xi_5 of eta ===");
{my(A=Reta[1],B=Reta[2]);
 for(k=1,6,my(n=100*k);print("  n=",n," v_5(b_n/a_n - b_{n-1}/a_{n-1})=",valuation(B[n+1]/A[n+1]-B[n]/A[n],5)));}

print("\n=== other primes for delta,eta (p=2,7 sanity) ===");
{my(A=Rdelta[1],B=Reta[1]);
 for(k=1,4,my(n=200*k);
  print("  delta n=",n," v_2=",valuation(Rdelta[2][n+1]/Rdelta[1][n+1]-Rdelta[2][n]/Rdelta[1][n],2),
        " v_7=",valuation(Rdelta[2][n+1]/Rdelta[1][n+1]-Rdelta[2][n]/Rdelta[1][n],7));
  print("  eta   n=",n," v_2=",valuation(Reta[2][n+1]/Reta[1][n+1]-Reta[2][n]/Reta[1][n],2),
        " v_7=",valuation(Reta[2][n+1]/Reta[1][n+1]-Reta[2][n]/Reta[1][n],7)));}

print("\n=== cross-determinant test: delta vs Domb, delta vs T at p=3, scanning small rational r ===");
{my(cands=[1,2,3,4,6,7,8,9,12,16,24,32,1/2,1/3,1/4,3/2,4/3,7/4,7/8,7/16,7/24,7/32]);
 print("delta vs Domb: v_3(r*a_d(n)*b_D(n) - a_D(n)*b_d(n)) at n=200,400,600,800 for each r");
 for(ci=1,#cands,my(r=cands[ci]);
   my(v=vector(4,k,my(n=200*k);valuation(r*Rdelta[1][n+1]*RDomb[2][n+1]-RDomb[1][n+1]*Rdelta[2][n+1],3)));
   print("  r=",r," : ",v));
 print("delta vs T: v_3(r*a_d(n)*b_T(n) - a_T(n)*b_d(n)) at n=200,400,600,800 for each r");
 for(ci=1,#cands,my(r=cands[ci]);
   my(v=vector(4,k,my(n=200*k);valuation(r*Rdelta[1][n+1]*RT[2][n+1]-RT[1][n+1]*Rdelta[2][n+1],3)));
   print("  r=",r," : ",v));
 print("delta vs Apery: v_3(r*a_d(n)*b_Ap(n) - a_Ap(n)*b_d(n)) at n=200,400,600,800 for each r");
 for(ci=1,#cands,my(r=cands[ci]);
   my(v=vector(4,k,my(n=200*k);valuation(r*Rdelta[1][n+1]*RApery[2][n+1]-RApery[1][n+1]*Rdelta[2][n+1],3)));
   print("  r=",r," : ",v));}
\q
