\\ Integer model: a_n integer;  b_n = c_n/(n!)^3 with
\\   c_{n+1} = (2n+1)(a n^2+a n+b) c_n - C n^6 c_{n-1},  c_0=0, c_1=1.
NMAX=800;
rowsI(a,b,C)={my(A=vector(NMAX+1),Cc=vector(NMAX+1));A[1]=1;A[2]=b;Cc[1]=0;Cc[2]=1;
 for(n=1,NMAX-1,
   A[n+2]=((2*n+1)*(a*n^2+a*n+b)*A[n+1]-C*n^3*A[n])/(n+1)^3;
   Cc[n+2]=(2*n+1)*(a*n^2+a*n+b)*Cc[n+1]-C*n^6*Cc[n];
 );[A,Cc]};
v2fact(n)=n-hammingweight(n);
