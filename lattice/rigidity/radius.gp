N=400;
rows(a,b,c)={my(A=vector(N+1),B=vector(N+1));A[1]=1;A[2]=b;B[1]=0;B[2]=1;
 for(n=1,N-1,A[n+2]=((a*n^2+a*n+b)*A[n+1]-c*n^2*A[n])/(n+1)^2;B[n+2]=((a*n^2+a*n+b)*B[n+1]-c*n^2*B[n])/(n+1)^2);[A,B]};
RB=rows(9,3,27); RC=rows(10,3,9); RF=rows(17,6,72); RR=[RB,RC,RF]; NM=["B","C","F"];
for(j=1,3, R=RR[j]; c0=sum(n=1,N, valuation(R[1][n+1],3)==1); mx=vecmax(vector(N+1,i,valuation(R[1][i],3))); print(NM[j],":  #{1<=n<=",N,": v3(a_n)=1} = ",c0,"    max v3(a_n) = ",mx));
