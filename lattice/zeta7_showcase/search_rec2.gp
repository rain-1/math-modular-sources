default(parisize,"16G");
cn=read("cn.txt");
NC=900;
p=2^61-1;
cp=vector(NC+1,i,Mod(cn[i],p));
gettime();
trial(r,D)=
{
  my(sz=(r+1)*(D+1), nrows, M);
  nrows=min(NC-r+1, sz+40);
  if(nrows<sz+10, return(-1));
  M=matrix(nrows,sz,a,b,my(n=r+a-1, j=(b-1)\(D+1), k=(b-1)%(D+1)); Mod(n,p)^k*cp[n-j+1]);
  return(#matker(M));
}
{for(r=6,30, my(row=vector(14,ii,my(D=ii+2, sz=(r+1)*(D+1)); if(sz>430,-9,trial(r,D)))); print("r=",r,"  dims D=3..16: ",row));}
print("time ",gettime());
quit;
