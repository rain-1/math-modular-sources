default(parisizemax,6000000000);
default(realprecision,1100);
N=420;
a=vector(N+2); b=vector(N+2); a[1]=1;a[2]=10;b[1]=0;b[2]=1;
{for(n=1,N, a[n+2]=((136*n^2+68*n+10)*a[n+1]-4*(2*n-1)^2*a[n])/(n+1)^2;
            b[n+2]=((136*n^2+68*n+10)*b[n+1]-4*(2*n-1)^2*b[n])/(n+1)^2);}
\\ exact Casoratian check
print("Casoratian a_n b_{n+1}-a_{n+1} b_n == Catalan_n^2 for n<=100? ",
  vector(100,n,a[n]*b[n+1]-a[n+1]*b[n])==vector(100,n,(binomial(2*(n-1),n-1)/n)^2));
xi = sum(m=0,N, (binomial(2*m,m)/(m+1))^2/(a[m+1]*a[m+2])*1.0);
print("xi = ",xi);
write("lattice/sqrt_apery/xi.txt", xi);
\\ sanity vs b/a
print("xi - b_N/a_N = ",xi - b[N+1]*1.0/a[N+1]);
\q
