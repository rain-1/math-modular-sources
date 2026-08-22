/* End-to-end check of the Apery criterion for the row
     (n+1)^2 a_{n+1} = (234n^2-39n-78) a_n + (27n^2-63n+30) a_{n-1},
     a_0=1, a_1=-78 ;  b_0=0, b_1=1.                                        */
default(parisizemax, 8000000000);
default(realprecision, 400);
{
N=600;
a=vector(N+2); b=vector(N+2); a[1]=1;a[2]=-78; b[1]=0;b[2]=1;
for(n=1,N, my(P=234*n^2-39*n-78, Q=-27*n^2+63*n-30);
  a[n+2]=(P*a[n+1]-Q*a[n])/(n+1)^2; b[n+2]=(P*b[n+1]-Q*b[n])/(n+1)^2);
print("a_n in Z for n<=",N+1,": ", vecmax(vector(N+2,i,denominator(a[i])))==1);
D=vector(N+2); D[1]=1; for(i=2,N+2, D[i]=lcm(D[i-1],i-1));
print("d_n^2 b_n in Z: ", vecmax(vector(N+2,i,denominator(D[i]^2*b[i])))==1);
print("d_n^1 b_n in Z: ", vecmax(vector(N+2,i,denominator(D[i]*b[i])))==1);
W=vector(N+1,i,a[i]*b[i+1]-a[i+1]*b[i]);
print("W_n = 0 for some n<=",N,"? ", vecmin(vector(N,i,if(W[i]==0,0,1)))==0);
S=sum(m=1,N, W[m]/(a[m]*a[m+1]));   xi=S*1.0;
print("xi = ", xi);
print("");
print("  n     (1/n)log|a_n|    (1/n)log|r_n|    (1/n)log(d_n^2|r_n|)");
for(j=1,5, my(n=100*j);
  my(T=sum(m=n,N-1, W[m]/(a[m]*a[m+1])), r=a[n+1]*T);
  print("  ",n,"   ", log(abs(a[n+1]*1.0))/n, "   ", log(abs(r*1.0))/n, "   ",
        (2*log(D[n+1]*1.0)+log(abs(r*1.0)))/n));
print("");
print("log lam1 = ", log((234+sqrt(234^2+4*27))/2), "   log lam2 = ", log(abs((234-sqrt(234^2+4*27))/2)));
print("score = ", -log(abs((234-sqrt(234^2+4*27))/2))-2);
print("mu <= ", 1 + (2+log((234+sqrt(234^2+4*27))/2))/(-log(abs((234-sqrt(234^2+4*27))/2))-2));
}
