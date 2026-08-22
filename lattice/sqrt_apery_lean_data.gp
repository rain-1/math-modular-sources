default(parisizemax,2000000000); default(realprecision,80);
N=60; A=vector(N+1);B=vector(N+1);A[1]=1;A[2]=10;B[1]=0;B[2]=1;
for(n=1,N-1,A[n+2]=((136*n^2+68*n+10)*A[n+1]-4*(2*n-1)^2*A[n])/(n+1)^2;B[n+2]=((136*n^2+68*n+10)*B[n+1]-4*(2*n-1)^2*B[n])/(n+1)^2);
print("a_0..a_10: ",vector(11,i,A[i]));
print("b_1..b_8: ",vector(8,i,B[i+1]));
dn=1;print("d_n^2 b_n, n=1..8: ",vector(8,n,dn=lcm(dn,n);dn^2*B[n+1]));
print("Casoratian a_n b_{n-1}-a_{n-1} b_n, n=1..6: ",vector(6,n,A[n+1]*B[n]-A[n]*B[n+1]),"  Catalan^2: ",vector(6,n,(binomial(2*n-2,n-1)/n)^2));
print("ratios a_{n+1}/a_n: ",vector(8,n,A[n+2]*1./A[n+1]));
print("min over n>=3 of a_{n+1}/a_n (n<=59): ",vecmin(vector(57,n,A[n+4]*1./A[n+3])));
print("xi=",B[N+1]*1./A[N+1]);
\\ Apery numbers and check a_n = 4^n [t^n] sqrt(sum A t^n) via Catalan convolution formula
Ap=vector(12);Ap[1]=1;Ap[2]=5;for(n=1,10,Ap[n+2]=((2*n+1)*(17*n^2+17*n+5)*Ap[n+1]-n^3*Ap[n])/(n+1)^3);
print("Apery A_n: ",Ap);
\\ F sqrt(1-34t+t^2) integrality check
s=Ser(Ap,x)+O(x^12); print("F*sqrt(1-34t+t^2): ",Vec(s*sqrt(1-34*x+x^2+O(x^12))));
\q
