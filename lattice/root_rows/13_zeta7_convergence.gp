read("lattice/root_rows/lib.gp");
default(realprecision,120);
N=600;
ap=vector(N+1,i,my(m=i-1);sum(k=0,m,binomial(m,k)^2*binomial(m+k,k)^2));
Aser=sum(m=0,N,ap[m+1]*'s^m)+O('s^(N+1));
At=Aser^3*sqrt(1-34*'s+'s^2+O('s^(N+1)));
Av=vector(N+1,i,polcoeff(At,i-1));
rr=rootrow(Av,6,N); a=rr[2];
r=[1,-612,93798,-49572,6561];  \\ placeholder
mr=minrec(a,4,4); rq=mr[3];
qs=vector(5,j,subst(rq[j],'n,'n+(j-1)));
b=vector(N+1); b[1]=0; b[2]=1;
{for(m=2,N, my(t=0); for(j=1,4, if(m-j>=0, t+=subst(qs[j+1],'n,m-j)*b[m-j+1])); b[m+1]=-t/m^2);}
print("a_n integral to ",N,": ",sum(i=1,N+1,denominator(a[i])!=1)==0);
{my(dn=1,k2=1,k1=1); for(n=1,N,dn=lcm(dn,n);
  if(denominator(dn^2*b[n+1])!=1,k2=0); if(denominator(dn*b[n+1])!=1,k1=0));
 print("d_n^2 b_n in Z to n=",N,": ",k2,"    d_n b_n in Z: ",k1);}
{my(prev=0); for(j=1,10, my(n=60*j);
  my(d=b[n+1]/a[n+1]-b[n]/a[n]);
  print("n=",n,"  b_n/a_n = ",b[n+1]/a[n+1]*1.0,"   diff=",d*1.0,
        if(prev!=0, Str("   log-ratio exponent = ", log(abs(d/prev))/log(n/(n-60.))), ""));
  prev=d);}
print("xi ~ ",b[N+1]/a[N+1]*1.0);
\q
