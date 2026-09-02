default(parisize,"48G");
DN=read("dn3000.txt");
NN=3000;
p=2^61-1;
dv=vector(NN+1,i,Mod(DN[i],p));
trial(r,D)=my(sz=(r+1)*(D+1));my(nrows=min(NN-r+1,sz+25));if(nrows<sz+10,return(-1));my(M=matrix(nrows,sz,a,b,my(nn=r+a-1);my(j=(b-1)\(D+1));my(k=(b-1)%(D+1));Mod(nn,p)^k*dv[nn-j+1]));#matker(M);
print("denominators cleared mod p? ",vecmin(vector(NN+1,i,denominator(DN[i])%p!=0)));
{foreach([[13,49],[13,60],[14,48],[14,49],[14,50],[15,28],[15,29],[18,15],[17,15]],rd,print("d_n: (r,D)=",rd,"  dim=",trial(rd[1],rd[2]),"  predicted (r-13)(D-6)-42 = ",max(0,(rd[1]-13)*(rd[2]-6)-42)));}
\\ primed d' row
dp=vector(NN,i,Mod(DN[i+1]+DN[i],p));
trialp(r,D)=my(sz=(r+1)*(D+1));my(nrows=min(NN-r,sz+25));if(nrows<sz+10,return(-1));my(M=matrix(nrows,sz,a,b,my(nn=r+a);my(j=(b-1)\(D+1));my(k=(b-1)%(D+1));Mod(nn,p)^k*dp[nn-j]));#matker(M);
{foreach([[12,49],[13,49],[13,48],[17,15],[18,15]],rd,print("d'_n: (r,D)=",rd,"  dim=",trialp(rd[1],rd[2]),"  predicted (r-12)(D-6)-42 = ",max(0,(rd[1]-12)*(rd[2]-6)-42)));}
quit;
