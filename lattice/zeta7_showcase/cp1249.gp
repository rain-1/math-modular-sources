default(parisize,"48G");
C=read("cn3000.txt");
NN=3000;
kerp(r,D,p)=my(cp=vector(NN,i,Mod(C[i+1]+C[i],p)));my(sz=(r+1)*(D+1));my(nrows=min(NN-r,sz+25));my(M=matrix(nrows,sz,a,b,my(nn=r+a);my(j=(b-1)\(D+1));my(k=(b-1)%(D+1));Mod(nn,p)^k*cp[nn-j]));matker(M);
p=2^61-1;
{foreach([[11,49],[11,60],[12,48],[12,49],[12,50],[13,28],[54,7],[53,7]],rd,print("(r,D)=",rd,"  dim=",#kerp(rd[1],rd[2],p),"  predicted (r-11)(D-6)-42 = ",max(0,(rd[1]-11)*(rd[2]-6)-42)));}
r=12;D=49;
{foreach([2^61-1,precprime(2^61-100),precprime(2^60)],pp,
 my(K=kerp(r,D,pp));
 print("p=",pp," dim=",#K);
 my(v=vector(r+1,jj,K[(jj-1)*(D+1)+D+1,1]));
 print("  leading p_j/p_0: ",vector(r+1,i,my(z=lift(v[i]/v[1]));if(z>pp/2,z-pp,z)));
 my(ok=0);
 for(a=0,r,if((r-a)%2==1,next);my(b=(r-a)/2);my(pol=(t+1)^a*(t^2-14*t+1)^b);my(w=vector(r+1,i,Mod(polcoeff(pol,r-i+1),pp)));my(o=1);for(i=1,r+1,if(v[i]*w[1]!=w[i]*v[1],o=0));if(o,print("  CLEAN chi = (L+1)^",a,"(L^2-14L+1)^",b);ok=1));
 if(!ok,print("  no clean chi"));
 my(P0=sum(k=0,D,(K[k+1,1]/K[D+1,1])*n^k));
 print("  valuation of P_0 at n=0: ",valuation(P0,n)));}
quit;
