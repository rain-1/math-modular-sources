default(parisize,"40G");
default(realprecision,1200);
cn=read("cn.txt");
dn=read("dn900.txt");
xi=(209/1728)*zeta(7);
s3=sqrt(3); xp=7-4*s3;
K7=3^(27/4)*(739-356*s3)/2^(77/6)*gamma(1/3)^12/Pi^(19/2);
al=2023/720;
NN=900;
Lg=log(1+x+O(x^(NN+2)));
S7=Lg^7/(1+x+O(x^(NN+2)));
lead(n)=(al/K7)*(-1)^(n+1)*xp^n*n^(3/2)*log(n)^7;
lead2(n)=if(n<7,1,(al/K7)*polcoeff(S7,n)*xp^n*n^(3/2));
print("n | digits zeta(7) | |d/c-xi|^(1/n) | err/leading(log^7) | err/leading(exact S7)");
{foreach([5,10,20,50,100,200,400,900],n,
  my(r=dn[n+1]/cn[n+1], er=r-xi);
  print(n," | ",precision(-log(abs((1728/209)*r-zeta(7)))/log(10),8),
        " | ",precision(abs(er)^(1/n),8),
        " | ",precision(er/lead(n),8),
        " | ",precision(er/lead2(n),8)));}
print();
print("check e_n = d_n - xi c_n against alpha*[x^n](log^7(1+x)/(1+x)):");
{foreach([100,200,400,700,900],n, my(en=dn[n+1]-xi*cn[n+1]); print("  n=",n,"  e_n/(alpha*S7_n) = ",precision(en/(al*polcoeff(S7,n)),15)));}
quit;
