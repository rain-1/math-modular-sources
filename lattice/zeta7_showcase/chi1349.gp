default(parisize,"40G");
cn=read("cn.txt");
NC=900;
p=2^61-1;
cp=vector(NC+1,i,Mod(cn[i],p));
r=13; D=49;
sz=(r+1)*(D+1);
nrows=min(NC-r+1,sz+25);
M=matrix(nrows,sz,a,b,my(n=r+a-1, j=(b-1)\(D+1), k=(b-1)%(D+1)); Mod(n,p)^k*cp[n-j+1]);
K=matker(M);
print("dim = ",#K);
lead=vector(r+1,jj,K[(jj-1)*(D+1)+D+1,1]);
print("leading vector p_j (mod p), normalised by p_0:");
print(vector(r+1,i,lift(lead[i]/lead[1])));
{for(a=0,r, if((r-a)%2==1,next); my(b=(r-a)/2, pol=(t+1)^a*(t^2-14*t+1)^b, v=vector(r+1,i,Mod(polcoeff(pol,r-i+1),p)));
  my(ok=1); for(i=1,r+1, if(lead[i]*v[1]!=v[i]*lead[1], ok=0));
  if(ok, print("CLEAN chi = (L+1)^",a,"(L^2-14L+1)^",b)));}
\\ also print P_0's structure: is n^7 | P_0 ?
P0=sum(k=0,D, lift(K[k+1,1]/lead[1])*n^k);
print("P_0 (normalised) divisible by n^7 mod p? ", Mod(1,p)*P0 == Mod(1,p)*(n^7*(P0\n^7)));
print("valuation of P_0 at n=0 (mod p): ",valuation(Mod(1,p)*P0,n));
quit;
