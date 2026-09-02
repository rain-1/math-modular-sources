default(parisize,"20G");
default(realprecision,420);
s3=sqrt(3);
fv=read("foldvals.txt");
K7fold=fv[1];
kapfold=fv[6];
kap=3^(27/4)*(739-356*s3)/2^(77/6);
print("kappa (new closed form 3^(27/4)(739-356 sqrt3)/2^(77/6)) = ",precision(kap,60));
print("kappa (fold, 400 d) = ",precision(kapfold,60));
print("closed/fold - 1 = ",precision(kap/kapfold-1,10));
K7=kap*gamma(1/3)^12/Pi^(19/2);
print("K7 closed = ",precision(K7,60));
print("K7 fold   = ",precision(K7fold,60));
print("K7 closed/fold - 1 = ",precision(K7/K7fold-1,10));
print("norm(739-356 sqrt3) = ",739^2-3*356^2,"  = ",factor(739^2-3*356^2));
k6=Mod(3^40*w*(739-356*w)^6,w^2-3)/2^77;
print("kappa^6 exact = ",lift(k6));
print("  numeric = ",precision(subst(lift(k6),w,s3),40),"   kappa^6 = ",precision(kap^6,40));
print("minpoly(kappa^6) = ",minpoly(k6));
mp=subst(minpoly(k6),x,y^6);
print("minpoly(kappa) deg 12 = ",mp);
print("  content: ",content(mp),"  primitive: ",mp/content(mp));
print("  check kappa root: ",precision(subst(mp/content(mp),y,kap),10));
\\ Richardson
cn=read("cn.txt");
xp=7-4*s3;
Kn(n)=cn[n+1]*xp^n*n^(3/2);
rich(ns)={my(m=#ns,V=matrix(m,m,i,j,1.0/ns[i]^(j-1)),b=vector(m,i,Kn(ns[i])));(matsolve(V,b~))[1];}
print("=== Richardson from exact c_n (n<=900), nodes n=900-25(i-1)");
{for(m=4,26,my(ns=vector(m,i,900-25*(i-1)),val=rich(ns)); print("m=",m,"  K7 ~ ",precision(val,55),"   rel.diff to closed = ",precision(val/K7-1,6)));}
quit;
