read("fitrec.gp");
M=150; default(seriesprecision,M+5);
et(k)=eta(q^k+O(q^(M+5)));
u = 3*q*et(2)^2*et(12)^4/(et(4)^4*et(6)^2);
ps1 = et(2)^6*et(3)/(et(1)^3*et(6)^2);
F12 = ps1*(3-4*u-3*u^2)/(3*(1+u)^2);
tt = serreverse(u/3); Ft = subst(F12, q, tt);
A = vector(M, n, polcoeff(Ft, n-1));
{for(ord=2,5,for(d=1,5,r=fitrec(A,ord,d);if(r!=0,print("order ",ord," deg ",d,": ",r);break(2))));}
print("done");

