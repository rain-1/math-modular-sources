default(parisize,"20G");
cn=read("cn.txt");
NC=900;
p=2^61-1;
cp=vector(NC+1,i,Mod(cn[i],p));
ker(r,D)=
{
  my(sz=(r+1)*(D+1), nrows, M);
  nrows=min(NC-r+1, sz+30);
  M=matrix(nrows,sz,a,b,my(n=r+a-1, j=(b-1)\(D+1), k=(b-1)%(D+1)); Mod(n,p)^k*cp[n-j+1]);
  return(matker(M));
}
lead(K,r,D)=vector(r+1,jj,K[(jj-1)*(D+1)+D+1]);
test(r,D)=
{
  my(K=ker(r,D), nb=#K, L, target, res=List(), MM, kk);
  if(nb==0, return([]));
  L=matrix(r+1,nb,i,j,lead(K[,j],r,D)[i]);
  for(a=0,r, if((r-a)%2==1, next); my(b=(r-a)/2, pol=(t+1)^a*(t^2-14*t+1)^b);
    target=vector(r+1,i,Mod(polcoeff(pol,r-i+1),p))~;
    MM=matconcat([L,-target]); kk=matker(MM);
    if(#kk>0, my(ok=0); for(s=1,#kk, if(kk[nb+1,s]!=0, ok=1)); if(ok, listput(res,[a,b]))));
  return(Vec(res));
}
{foreach([[17,15],[18,14],[19,13],[20,12],[21,11],[22,11],[23,10],[24,10],[27,9],[34,8],[55,7]],rd,
  print("r=",rd[1]," D=",rd[2]," dim=",#ker(rd[1],rd[2]),"  clean char polys (a,b): ",test(rd[1],rd[2])));}
quit;
