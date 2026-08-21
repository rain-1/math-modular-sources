\\ L(f,2) row: (n+1)^2 a_{n+1} = (20n^2+10n+2) a_n - 16(2n-1)^2 a_{n-1}; companion b_0=0,b_1=1.
N=600; A=vector(N+1);B=vector(N+1);A[1]=1;A[2]=2;B[1]=0;B[2]=1;
for(n=1,N-1,A[n+2]=((20*n^2+10*n+2)*A[n+1]-16*(2*n-1)^2*A[n])/(n+1)^2;B[n+2]=((20*n^2+10*n+2)*B[n+1]-16*(2*n-1)^2*B[n])/(n+1)^2);
print("a_n: ",vector(8,i,A[i]));
print("v2(b_n/a_n) at n=50..600 step 50: ",vector(12,k,valuation(B[50*k+1]/A[50*k+1],2)));
print("v2(a_n) max: ",vecmax(vector(N,n,valuation(A[n+1],2))),"  v3 slope check v3(b/a diff)/n at N: ",valuation(B[N+1]/A[N+1]-B[N]/A[N],3)*1./N);
\q
