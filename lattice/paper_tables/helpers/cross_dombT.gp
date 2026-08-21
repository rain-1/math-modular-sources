default(realprecision,80);
rows3(a,b,c,N)={my(A=vector(N+1),B=vector(N+1));A[1]=1;A[2]=b;B[1]=0;B[2]=1;
 for(n=1,N-1,A[n+2]=((2*n+1)*(a*n^2+a*n+b)*A[n+1]-c*n^3*A[n])/(n+1)^3;
            B[n+2]=((2*n+1)*(a*n^2+a*n+b)*B[n+1]-c*n^3*B[n])/(n+1)^3);[A,B]};
N=400;
RD=rows3(10,4,64,N);
RT=rows3(12,4,16,N);
ns=[50,100,150,200,250,300,350,400];
vals=vector(#ns,k,my(n=ns[k]); valuation(3*RT[1][n+1]*RD[2][n+1]-4*RD[1][n+1]*RT[2][n+1],2));
print(vals);
\q
