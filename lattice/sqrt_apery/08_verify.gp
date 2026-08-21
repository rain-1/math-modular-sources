default(parisizemax,6000000000);
{print("[L1] binom(1/2,k)*4^k == (-1)^(k-1)*2*Cat(k-1), k=1..80: ",
  vector(80,k, binomial(1/2,k)*4^k) == vector(80,k,(-1)^(k-1)*2*binomial(2*k-2,k-1)/k));}
N=700; a=vector(N+2); b=vector(N+2); a[1]=1;a[2]=10;b[1]=0;b[2]=1;
{for(n=1,N, a[n+2]=((136*n^2+68*n+10)*a[n+1]-4*(2*n-1)^2*a[n])/(n+1)^2;
            b[n+2]=((136*n^2+68*n+10)*b[n+1]-4*(2*n-1)^2*b[n])/(n+1)^2);}
{print("[C] Casoratian == Catalan_n^2 for n<=",N,": ",
  vector(N,n,a[n]*b[n+1]-a[n+1]*b[n])==vector(N,n,(binomial(2*n-2,n-1)/n)^2));}
{print("[P] a_n in Z, a_n>=10a_{n-1}, n<=",N,": ",
  vector(N+2,i,denominator(a[i]))==vector(N+2,i,1), " ",
  vector(N,n,a[n+2]>=10*a[n+1])==vector(N,n,1));}
dn=1; ok2=1; {for(n=1,N, dn=lcm(dn,n); if(denominator(dn^2*b[n+1])!=1, ok2=0));}
print("[D] d_n^2 b_n in Z to ",N,": ",ok2);
\\ exact positive-series linear form r_n = a_n * sum_{m>=n} C_m^2/(a_m a_{m+1})  -- no cancellation
default(realprecision,80);
r(n) = a[n+1]*sum(m=n,N-1,(binomial(2*m,m)/(m+1))^2/(a[m+1]*a[m+2])*1.0);
dn=1; L=List();
{for(n=1,600, dn=lcm(dn,n); if(n%100==0, listput(L,[n, log(r(n))/n, log(dn^2*r(n))/n, log(dn^2*a[n+1])/n])));}
print("n | log r_n /n | log(d_n^2 r_n)/n | log(d_n^2 a_n)/n");
for(i=1,#L, print("  ",L[i][1],"  ",L[i][2],"  ",L[i][3],"  ",L[i][4]));
print("targets: log lambda2=",log(4*(17-12*sqrt(2))),"  2+log lambda2=",2+log(4*(17-12*sqrt(2))),"  2+log lambda1=",2+log(4*(17+12*sqrt(2))));
\q
