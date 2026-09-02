default(parisize,"48G");
default(realprecision,600);
C=read("cn3000.txt");
NN=3000;
s3=sqrt(3);xp=7-4*s3;
K7=3^(27/4)*(739-356*s3)/2^(77/6)*gamma(1/3)^12/Pi^(19/2);
Kn(nn)=C[nn+1]*xp^nn*nn^(3/2);
rich(ns)=my(m=#ns);my(V=matrix(m,m,i,j,1.0/ns[i]^(j-1)));my(b=vector(m,i,Kn(ns[i])));(matsolve(V,b~))[1];
print("K7 closed form = ",precision(K7,60));
{foreach([10,20,30,40,50],m,my(ns=vector(m,i,3000-45*(i-1)));my(v=rich(ns));print("  m=",m,"  rel.diff to closed form = ",precision(v/K7-1,6)));}
print();
print("=== symmetry scan for the primed (12,49) recurrence");
p=2^61-1;
cp=vector(NN,i,Mod(C[i+1]+C[i],p));
r=12;D=49;
sz=(r+1)*(D+1);
nrows=min(NN-r,sz+25);
M=matrix(nrows,sz,a,b,my(nn=r+a);my(j=(b-1)\(D+1));my(k=(b-1)%(D+1));Mod(nn,p)^k*cp[nn-j]);
K=matker(M)[,1];
PJ=vector(r+1,jj,sum(k=0,D,(K[(jj-1)*(D+1)+k+1]/K[D+1])*n^k));
found=List();
{for(A2=-240,240,my(A=A2/2);for(e=0,1,my(ep=if(e==0,1,-1));my(ok=1);for(j=0,r,if(PJ[r-j+1]-ep*subst(PJ[j+1],n,A-n)!=0,ok=0;break));if(ok,listput(found,[A,ep]))));}
print("symmetries P_{12-j}(n)=eps P_j(A-n), A in (1/2)Z, |A|<=120: ",Vec(found));
quit;
