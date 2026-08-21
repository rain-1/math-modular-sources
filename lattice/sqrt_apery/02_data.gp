default(parisizemax,4000000000);
N=200;
a=vector(N+1); b=vector(N+1); a[1]=1;a[2]=10;b[1]=0;b[2]=1;
{for(n=1,N-1, a[n+2]=((136*n^2+68*n+10)*a[n+1]-4*(2*n-1)^2*a[n])/(n+1)^2;
             b[n+2]=((136*n^2+68*n+10)*b[n+1]-4*(2*n-1)^2*b[n])/(n+1)^2);}
print("a_n n=0..8: ",vector(9,i,a[i]));
print("b_n n=0..8: ",vector(9,i,b[i]));
print("den(b_n) n=1..24: ",vector(24,i,denominator(b[i+1])));
dn(n)=lcm(vector(max(n,1),i,i));
print("factor den(b_n) n=1..18:");
for(n=1,18, print("  n=",n,"  den=",factor(denominator(b[n+1]))));
print("d_n^2 b_n integral to ",N,"? ", vector(N,n,denominator(dn(n)^2*b[n+1]))==vector(N,n,1));
print("d_n b_n integral? ", vector(N,n,denominator(dn(n)*b[n+1]))==vector(N,n,1));
print("den(b_n) | d_n^2 quotient n=1..20: ",vector(20,n,dn(n)^2/denominator(b[n+1])));
\q
