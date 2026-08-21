default(parisizemax,8000000000);
N=1500; a=vector(N+2); b=vector(N+2); a[1]=1;a[2]=10;b[1]=0;b[2]=1;
{for(n=1,N, a[n+2]=((136*n^2+68*n+10)*a[n+1]-4*(2*n-1)^2*a[n])/(n+1)^2;
            b[n+2]=((136*n^2+68*n+10)*b[n+1]-4*(2*n-1)^2*b[n])/(n+1)^2);}
print("[P] a_n in Z to ",N,": ",vector(N+2,i,denominator(a[i]))==vector(N+2,i,1));
dn=1; ok=1; {for(n=1,N, dn=lcm(dn,n); if(denominator(dn^2*b[n+1])!=1, ok=0; print("FAIL n=",n)));}
print("[D] d_n^2 b_n in Z to ",N,": ",ok);
{print("[C] Casoratian ok to ",N,": ",vector(N,n,a[n]*b[n+1]-a[n+1]*b[n])==vector(N,n,(binomial(2*n-2,n-1)/n)^2));}
\q
