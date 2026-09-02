default(parisize,"20G");
BASE="/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/9a849c0a-95f8-4d19-b342-98033d0d9c03/scratchpad/zeta7/";
DN=read(concat(BASE,"dn.txt"));
ND=400;
print("#dn = ",#DN);
p=2^61-1;
badden=0;
{for(i=1,ND+1, if(Mod(denominator(DN[i]),p)==0, badden=badden+1));}
print("denominators divisible by p: ",badden);
DR=vector(ND+1,i,lift(Mod(DN[i],p)));
kd(r,D)=my(n0=r+1,nr=ND-n0+1,sz=(r+1)*(D+1),PW,M);if(nr<sz+15,return(-1));PW=matrix(nr,D+1,a,b,lift(Mod(n0+a-1,p)^(b-1)));M=matrix(nr,sz,a,b,my(j=(b-1)\(D+1),k=(b-1)%(D+1),nn=n0+a-1);PW[a,k+1]*DR[nn-j+1]%p);#matker(M*Mod(1,p));
gettime();
{for(r=8,44, my(v=List(),any=0); for(D=0,28, my(z=kd(r,D)); listput(v,z); if(z>0,any=1)); print("r'=",r," D'=0..28 dims: ",Vec(v),"   t=",gettime()));}
quit;
