default(parisize,"40G");
cn=read("cn.txt");
NC=900;
r=13; D=49;
sz=(r+1)*(D+1);
nrows=min(NC-r+1,sz+25);
getPJ(p)=
{
  my(cp=vector(NC+1,i,Mod(cn[i],p)), M, K, nrm);
  M=matrix(nrows,sz,a,b,my(nn=r+a-1, j=(b-1)\(D+1), k=(b-1)%(D+1)); Mod(nn,p)^k*cp[nn-j+1]);
  K=matker(M);
  if(#K!=1, return(0));
  K=K[,1]; nrm=K[D+1];
  vector(r+1,jj,sum(k=0,D, (K[(jj-1)*(D+1)+k+1]/nrm)*n^k));
}
{foreach([2^61-1, precprime(2^61-100), precprime(2^60)],p,
  my(PJ=getPJ(p), lead, ok=0);
  lead=vector(r+1,i,polcoeff(PJ[i],D));
  for(a=0,r, if((r-a)%2==1,next); my(b=(r-a)/2, pol=(t+1)^a*(t^2-14*t+1)^b, v=vector(r+1,i,Mod(polcoeff(pol,r-i+1),p)), o=1);
    for(i=1,r+1, if(lead[i]*v[1]!=v[i]*lead[1], o=0));
    if(o, print("p=",p,"  CLEAN chi = (L+1)^",a,"(L^2-14L+1)^",b); ok=1));
  if(!ok, print("p=",p,"  no clean chi"));
  print("   valuation of P_0 at n=0: ",valuation(PJ[1],n)));}
p=2^61-1; PJ=getPJ(p);
found=List();
{for(A2=-160,240, my(A=A2/2); for(e=0,1, my(ep=if(e==0,1,-1), ok=1);
  for(j=0,r, if(PJ[r-j+1]-ep*subst(PJ[j+1],n,A-n)!=0, ok=0; break));
  if(ok, listput(found,[A,ep]))));}
print("symmetries with A in (1/2)Z, |A|<=120: ",Vec(found));
{my(ok=0); for(A3=-360,360, my(A=A3/3, o=1); for(j=0,r, if(PJ[r-j+1]+subst(PJ[j+1],n,A-n)!=0, o=0;break)); if(o, print("eps=-1 symmetry at A=",A); ok=1)); if(!ok, print("no eps=-1 symmetry with A in (1/3)Z, |A|<=120"));}
quit;
