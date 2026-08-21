par=[[9,3,27],[10,3,9],[17,6,72]]; N=300;
rows(a,b,c)={my(A=vector(N+1),B=vector(N+1));A[1]=1;A[2]=b;B[1]=0;B[2]=1;
 for(n=1,N-1,A[n+2]=((a*n^2+a*n+b)*A[n+1]-c*n^2*A[n])/(n+1)^2;B[n+2]=((a*n^2+a*n+b)*B[n+1]-c*n^2*B[n])/(n+1)^2);[A,B]};
RB=rows(9,3,27);RC=rows(10,3,9);RF=rows(17,6,72);
{for(k=1,6,n=50*k;
 print("n=",n,"  v3(4 aC bF - 5 aF bC)=",valuation(4*RC[1][n+1]*RF[2][n+1]-5*RF[1][n+1]*RC[2][n+1],3),
 "   v3(4 aB bF - 5 aF bB)=",valuation(4*RB[1][n+1]*RF[2][n+1]-5*RF[1][n+1]*RB[2][n+1],3),
 "   v3(aB bC - aC bB)=",valuation(RB[1][n+1]*RC[2][n+1]-RC[1][n+1]*RB[2][n+1],3)));}
\q
