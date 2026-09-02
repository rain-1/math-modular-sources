default(parisize,"24G");
BASE="/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/9a849c0a-95f8-4d19-b342-98033d0d9c03/scratchpad/zeta7/";
DN=read(concat(BASE,"dn900.txt"));
ND=900;
p=2^61-1;
DR=vector(ND+1,i,lift(Mod(DN[i],p)));
kd(r,D)=my(n0=r+1,nr=ND-r,sz=(r+1)*(D+1),PW,M);if(nr<sz+20,return(-1));PW=matrix(nr,D+1,a,b,lift(Mod(n0+a-1,p)^(b-1)));M=matrix(nr,sz,a,b,my(j=(b-1)\(D+1),k=(b-1)%(D+1),nn=n0+a-1);PW[a,k+1]*DR[nn-j+1]%p);#matker(M*Mod(1,p));
pred(r,D)=max(0,(r-13)*(D-6)-42);
gettime();
{foreach([[17,15],[18,15],[19,13],[19,14],[20,12],[20,13],[21,11],[21,12],[22,10],[22,11],[23,10],[24,10],[27,9],[28,9],[34,8],[35,8],[55,7],[56,7],[57,7],[18,14],[26,12],[30,11]],rd,my(r=rd[1],D=rd[2],z=kd(r,D));print("(",r,",",D,") dim=",z,"   predicted (r-13)(D-6)-42 -> ",pred(r,D),"   match=",z==pred(r,D),"   t=",gettime()));}
quit;
