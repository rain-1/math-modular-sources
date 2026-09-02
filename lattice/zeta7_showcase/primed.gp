default(parisize,"48G");
default(realprecision,3600);
C=read("cn3000.txt");
DN=read("dn3000.txt");
NN=3000;
xi=(209/1728)*zeta(7);
s3=sqrt(3);xp=7-4*s3;
K7=3^(27/4)*(739-356*s3)/2^(77/6)*gamma(1/3)^12/Pi^(19/2);
K7p=(1+xp)*K7;
al=14161/5040;
cp=vector(NN,i,C[i+1]+C[i]);
dp=vector(NN,i,DN[i+1]+DN[i]);
print("K7  = ",precision(K7,40));
print("K7' = (8-4sqrt3) K7 = ",precision(K7p,40));
Kn(nn)=cp[nn]*xp^nn*nn^(3/2);
rich(ns)=my(m=#ns);my(V=matrix(m,m,i,j,1.0/ns[i]^(j-1)));my(b=vector(m,i,Kn(ns[i])));(matsolve(V,b~))[1];
{foreach([10,16,22,28,34],m,my(ns=vector(m,i,3000-60*(i-1)));my(v=rich(ns));print("  Richardson m=",m," K7' = ",precision(v,50),"   rel.diff = ",precision(v/K7p-1,6)));}
print();
print("=== primed digit table");
lead(nn)=(7*al/K7p)*(-1)^(nn+1)*xp^nn*nn^(1/2)*log(nn)^6;
print("n | digits zeta(7) | |d'/c'-xi|^(1/n) | err/leading");
{foreach([5,10,20,50,100,200,400,900,1500,3000],nn,my(rr=dp[nn]/cp[nn]);my(er=rr-xi);print(nn," | ",precision(-log(abs((1728/209)*rr-zeta(7)))/log(10),8)," | ",precision(abs(er)^(1/nn),8)," | ",precision(er/lead(nn),8)));}
print();
print("=== unprimed row for comparison (same n)");
leadu(nn)=(al/K7)*(-1)^(nn+1)*xp^nn*nn^(3/2)*log(nn)^7;
{foreach([100,400,900,1500,3000],nn,my(rr=DN[nn+1]/C[nn+1]);my(er=rr-xi);print(nn," | ",precision(-log(abs((1728/209)*rr-zeta(7)))/log(10),8)," | ",precision(er/leadu(nn),8)));}
quit;
