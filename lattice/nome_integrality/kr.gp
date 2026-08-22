N=60;
A=vector(N+2); A[1]=1; A[2]=5;
{for(n=1,N, A[n+2] = ((2*n+1)*(17*n^2+17*n+5)*A[n+1] - n^3*A[n])/(n+1)^3);}
B=vector(N+2); B[1]=0;
{for(n=1,N+1, my(m=n-1,bm1,bm2,am,am1,am2); bm1=B[n]; bm2=if(n>=2,B[n-1],0); am=A[n+1]; am1=A[n]; am2=if(n>=2,A[n-1],0); B[n+1] = ( (2*n-1)*(17*m^2+17*m+5)*bm1 - m^3*bm2 - 3*n^2*am + (102*m^2+102*m+27)*am1 - 3*m^2*am2 )/n^3 );}
H(k) = sum(i=1,k,1/i);
KRB(n) = 2*sum(j=0,n, binomial(n,j)^2*binomial(n+j,j)^2*(H(n+j)-H(n-j)));
bad=0; {for(n=0,N, if(KRB(n)!=B[n+1], bad++; if(bad<4,print("MISMATCH n=",n," ODE=",B[n+1]," KR=",KRB(n)))));}
print("KR B_{2,2} vs ODE log-solution: mismatches over n<=",N," = ",bad);
print("B_0..B_5 = ", vector(6,i,B[i]));
\\ also check GKZ identity  C(k,j)C(k+j,j) = (k+j)!/(j!^2 (k-j)!)
bad=0; {for(k=0,40, for(j=0,k, if(binomial(k,j)*binomial(k+j,j) != (k+j)!/(j!^2*(k-j)!), bad++)));}
print("GKZ multinomial identity mismatches: ",bad);
quit;
