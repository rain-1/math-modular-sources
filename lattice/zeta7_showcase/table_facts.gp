default(parisize,"40G");
default(realprecision,620);
cn=read("cn.txt");
dn=read("dn900.txt");
xi=(209/1728)*zeta(7);
s3=sqrt(3); xp=7-4*s3;
K7=3^(27/4)*(739-356*s3)/2^(77/6)*gamma(1/3)^12/Pi^(19/2);
al=2023/720;
print("=== ITEM 6: Dirichlet polynomial");
cv=[1,-572,11583,-36608,46332,-20736]; dv=[1,2,3,4,6,12];
P(s)=sum(i=1,6,cv[i]*dv[i]^(-s));
{for(j=0,8, print("  P(",j,") = ",P(j)));}
print("  P(7) = ",P(7),"   -209/864 = ",-209/864,"   equal? ",P(7)==-209/864);
print("  xi = P(7)*zeta(0) = -P(7)/2 = ",209/1728," * zeta(7)");
print();
print("=== ITEM 2: exact display values");
print("  c_0..c_8 = ",vector(9,i,cn[i]));
print("  d_1..d_4 = ",vector(4,i,dn[i+1]));
print("  (1728/209) d_n/c_n for n=1..4 : ",vector(4,i,(1728/209)*dn[i+1]/cn[i+1]));
print();
print("=== digit table");
print("n | digits of zeta(7) | |d/c-xi|^(1/n) | err/leading | log10|err|");
lead(n)=(al/K7)*(-1)^(n+1)*xp^n*n^(3/2)*log(n)^7;
{foreach([5,10,20,50,100,200,400,900],n,
  my(r=dn[n+1]/cn[n+1], er=r-xi);
  print(n," | ",precision(-log(abs((1728/209)*r-zeta(7)))/log(10),8),
        " | ",precision(abs(er)^(1/n),8),
        " | ",precision(er/lead(n),8),
        " | ",precision(log(abs(er))/log(10),8)));}
print();
print("K7 = ",precision(K7,50));
print("alpha/K7 = ",precision(al/K7,30));
print("7-4sqrt3 = ",precision(xp,20),"  log(1/(7-4sqrt3)) = ",precision(-log(xp),20));
quit;
