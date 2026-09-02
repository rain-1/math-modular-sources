default(parisize,"40G");
cn=read("cn.txt");
NC=900;
p=2^61-1;
cp=vector(NC+1,i,Mod(cn[i],p));
ker(r,D)=
{
  my(sz=(r+1)*(D+1), nrows=min(NC-r+1,(r+1)*(D+1)+25), M);
  M=matrix(nrows,sz,a,b,my(nn=r+a-1, j=(b-1)\(D+1), k=(b-1)%(D+1)); Mod(nn,p)^k*cp[nn-j+1]);
  matker(M);
}
test(r,D)=
{
  my(K=ker(r,D), nb=#K, L, res=List(), tg, MM, kk);
  if(nb==0, return([]));
  L=matrix(r+1,nb,i,j,K[(i-1)*(D+1)+D+1,j]);
  for(a=0,r, if((r-a)%2==1,next); my(b=(r-a)/2, pol=(t+1)^a*(t^2-14*t+1)^b);
    tg=vector(r+1,i,Mod(polcoeff(pol,r-i+1),p))~;
    MM=matconcat([L,-tg]); kk=matker(MM);
    if(#kk>0, my(o=0); for(s=1,#kk, if(kk[nb+1,s]!=0, o=1)); if(o, listput(res,[a,b]))));
  Vec(res);
}
{foreach([[13,49],[14,28],[15,21],[16,17]],rd, print("(r,D)=",rd,"  dim=",#ker(rd[1],rd[2]),"  clean (a,b): ",test(rd[1],rd[2])));}
quit;
