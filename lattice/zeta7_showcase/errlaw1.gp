default(parisize,"20G");
default(realprecision,700);
cn=read("cn400.txt");
dn=read("dn.txt");
xi=(209/1728)*zeta(7);
e=vector(401,i,dn[i]-xi*cn[i]);
print("n | (-1)^n e_n | e_n*n/(log n)^6 | e_n/(log n)^7 | e_n*n/(log n)^7");
{foreach([20,40,60,80,100,150,200,250,300,350,400],n,
 my(s=(-1)^n, en=e[n+1]);
 print(n,"  ",precision(s*en,12),"   ",precision(s*en*n/log(n)^6,12),"   ",precision(s*en/log(n)^7,12),"   ",precision(s*en*n/log(n)^7,12)));}
print();
print("sign of e_n for n=380..400: ",vector(21,i,sign(e[380+i])));
quit;
