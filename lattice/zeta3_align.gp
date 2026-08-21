N=400;
rows(a,b,c)={my(A=vector(N+1),B=vector(N+1));A[1]=1;A[2]=b;B[1]=0;B[2]=1;
 for(n=1,N-1,A[n+2]=((2*n+1)*(a*n^2+a*n+b)*A[n+1]-c*n^3*A[n])/(n+1)^3;B[n+2]=((2*n+1)*(a*n^2+a*n+b)*B[n+1]-c*n^3*B[n])/(n+1)^3);[A,B]};
D=rows(10,4,64); T=rows(12,4,16); P=rows(17,5,1);
{for(k=1,8,n=50*k;
 print("n=",n,"  v2(3 aT bD - 4 aD bT)=",valuation(3*T[1][n+1]*D[2][n+1]-4*D[1][n+1]*T[2][n+1],2),
 "   v2(4 aP bD - 7 aD bP)=",valuation(4*P[1][n+1]*D[2][n+1]-7*D[1][n+1]*P[2][n+1],2),
 "   v2(a^T_n) ",valuation(T[1][n+1],2)," v2(a^D_n) ",valuation(D[1][n+1],2)));}
\q
