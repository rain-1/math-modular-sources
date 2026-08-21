default(realprecision,250);
rowsZ(a,b,c,N)={my(A=vector(N+1),B=vector(N+1));A[1]=1;A[2]=b;B[1]=0;B[2]=1;
 for(n=1,N-1,A[n+2]=((a*n^2+a*n+b)*A[n+1]-c*n^2*A[n])/(n+1)^2;
            B[n+2]=((a*n^2+a*n+b)*B[n+1]-c*n^2*B[n])/(n+1)^2);[A,B]};
N=300;
RB=rowsZ(9,3,27,N); RC=rowsZ(10,3,9,N); RF=rowsZ(17,6,72,N);
crossdet(X,Y,r,rp,p,ns)={vector(#ns,k,my(n=ns[k]);valuation(r*X[1][n+1]*Y[2][n+1]-rp*Y[1][n+1]*X[2][n+1],p))};
ns=[50,100,150,200,250,300];
print("B,C r=1,1 p=3: ",crossdet(RB,RC,1,1,3,ns));
print("B,F (4 xiF = 5 xiB) p=3: ",crossdet(RB,RF,4,5,3,ns));
print("C,F (4 xiF = 5 xiC) p=3: ",crossdet(RC,RF,4,5,3,ns));
\q
