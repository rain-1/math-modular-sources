default(parisize,"48G");
C=read("cn3000.txt");
NN=3000;
p=2^61-1;
cp=vector(NN,i,Mod(C[i+1]+C[i],p));
trial(r,D)=my(sz=(r+1)*(D+1));my(nrows=min(NN-r,sz+25));if(nrows<sz+10,return(-1));my(M=matrix(nrows,sz,a,b,my(nn=r+a);my(j=(b-1)\(D+1));my(k=(b-1)%(D+1));Mod(nn,p)^k*cp[nn-j]));#matker(M);
gettime();
print("dimension grid for c'_n = c_n + c_{n-1}   [predicted analogue: (r-13)(D-6)-42 ?]");
{foreach([[13,49],[14,49],[14,28],[15,28],[15,21],[16,21],[16,17],[17,17],[17,15],[18,15],[19,13],[20,13],[20,12],[21,12],[23,10],[24,10],[27,9],[28,9],[55,7],[56,7]],rd,
 my(d=trial(rd[1],rd[2]));
 print("  (r,D)=(",rd[1],",",rd[2],")  dim=",d,"   (r-12)(D-6)-42=",(rd[1]-12)*(rd[2]-6)-42,"   (r-13)(D-6)-42=",(rd[1]-13)*(rd[2]-6)-42,"   t=",gettime()));}
quit;
