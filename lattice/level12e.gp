getA(M)={default(seriesprecision,M+5);my(et(k)=eta(q^k+O(q^(M+5))));
 my(u=3*q*et(2)^2*et(12)^4/(et(4)^4*et(6)^2), ps1=et(2)^6*et(3)/(et(1)^3*et(6)^2));
 my(F12=ps1*(3-4*u-3*u^2)/(3*(1+u)^2), tt=serreverse(u/3), Ft=subst(F12,q,tt));
 vector(M,n,polcoeff(Ft,n-1))};
A1=getA(60); A2=getA(150);
print("agree on first 60? ", vector(60,i,A1[i])==vector(60,i,A2[i]));
print("first mismatch: ", for(i=1,60,if(A1[i]!=A2[i],print1(i," ");break)));
print(vector(10,i,A2[50+i]));
\q
