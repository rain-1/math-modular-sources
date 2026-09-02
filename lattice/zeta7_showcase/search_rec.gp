default(parisize,"16G");
cn=read("cn.txt");
NC=900;
p=2^61-1;
cp=vector(NC+1,i,Mod(cn[i],p));
gettime();
trial(r,D)=
{
  my(sz=(r+1)*(D+1), nrows, M, K);
  nrows=min(NC-r+1, sz+40);
  if(nrows<sz+5, return([-1]));
  M=matrix(nrows,sz,a,b,my(n=r+a-1, j=(b-1)\(D+1), k=(b-1)%(D+1)); Mod(n,p)^k*cp[n-j+1]);
  K=matker(M);
  return(K);
}
found=0;
{for(r=6,32, for(D=3,12, if(found,next); my(K=trial(r,D)); if(type(K)=="t_VEC", next); if(#K>0, print("FOUND r=",r," D=",D," kernel dim=",#K," time=",gettime()); found=[r,D];)));}
print("result: ",found);
quit;
