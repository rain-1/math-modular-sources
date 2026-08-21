M=200; default(seriesprecision,M+5);
et(k)=eta(q^k+O(q^(M+5)));
u = 3*q*et(2)^2*et(12)^4/(et(4)^4*et(6)^2);
ps1 = et(2)^6*et(3)/(et(1)^3*et(6)^2);
F = ps1*(3-4*u-3*u^2)/(3*(1+u)^2);
t = u/3;
Dq(f)=q*deriv(f,q);
th(f)=t*Dq(f)/Dq(t);          \\ theta_t = t d/dt
F0=F; F1=th(F0); F2=th(F1); F3=th(F2);
\\ find polynomials P_j(t) deg<=D with sum P_j(t) theta^j F = 0, order 2 then 3
{for(ord=2,3, for(D=2,12,
  my(cols=[]); for(j=0,ord, for(e=0,D, cols=concat(cols,[t^e*[F0,F1,F2,F3][j+1]])));
  my(Mx=matrix(M-2,#cols,i,k,polcoeff(cols[k],i-1))); my(K=matker(Mx));
  if(#K>0, print("order ",ord," deg ",D," kerdim ",#K); my(v=K[,1]); v=v/content(v);
    for(j=0,ord,print("  P",j," = ",factor(Polrev(vector(D+1,e,v[j*(D+1)+e]),x)))); break(2))));}
\q
