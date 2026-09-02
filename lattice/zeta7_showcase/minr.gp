default(parisize,"40G");
cn=read("cn.txt");
NC=900;
p=2^61-1;
cp=vector(NC+1,i,Mod(cn[i],p));
trial(r,D)=
{
  my(sz=(r+1)*(D+1), nrows, M);
  nrows=min(NC-r+1, sz+25);
  if(nrows<sz+10, return(-1));
  M=matrix(nrows,sz,a,b,my(n=r+a-1, j=(b-1)\(D+1), k=(b-1)%(D+1)); Mod(n,p)^k*cp[n-j+1]);
  return(#matker(M));
}
gettime();
{foreach([[16,16],[16,17],[16,18],[15,20],[15,21],[15,22],[14,27],[14,28],[14,29],[13,48],[13,49],[13,50],[12,60],[12,80]],rd,
  my(d=trial(rd[1],rd[2]), pred=max(0,(rd[1]-12)*(rd[2]-6)-42));
  print("(r,D)=(",rd[1],",",rd[2],")  dim=",d,"  predicted=",pred,"  match=",d==pred,"  t=",gettime()));}
quit;
