default(parisize,"20G");
cn=read("cn.txt");
NC=900;
p=2^61-1;
cp=vector(NC+1,i,Mod(cn[i],p));
trial(r,D)=
{
  my(sz=(r+1)*(D+1), nrows, M);
  nrows=min(NC-r+1, sz+30);
  if(nrows<sz+10, return(-1));
  M=matrix(nrows,sz,a,b,my(n=r+a-1, j=(b-1)\(D+1), k=(b-1)%(D+1)); Mod(n,p)^k*cp[n-j+1]);
  return(#matker(M));
}
gettime();
{for(D=7,9, my(v=List()); for(r=13,70, listput(v,trial(r,D))); print("D=",D," r=13..70 dims: ",Vec(v)); print("  time ",gettime()));}
quit;
