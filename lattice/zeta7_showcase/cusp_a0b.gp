default(parisize,"20G");
default(realprecision,220);
M=9000;
et(k)=prod(m=1,M\k+1,1-'q^(k*m)+O('q^(M+2)));
h=et(1)^3*et(4)*et(6)^2/(et(2)^2*et(3)*et(12)^3)/'q;
x=h/((h+3)*(h+4));
E8(k)=1+480*sum(m=1,M\k+1,sigma(m,7)*'q^(k*m))+O('q^(M+2));
Phi=(E8(1)-572*E8(2)+11583*E8(3)-36608*E8(4)+46332*E8(6)-20736*E8(12))/480;
U=Phi/('q*deriv(x,'q));
xs=truncate(x); Ps=truncate(Phi); Us=truncate(U);
evq(s,tau)=subst(s,'q,exp(2*Pi*I*tau));
rich(ys,w,fun)=
{
  my(m=#ys,V=matrix(m,m,i,j,exp(2*Pi*I*(I*ys[i])/w)^(j-1)),rhs=vector(m,i,fun(I*ys[i]))~);
  matsolve(V,rhs)[1];
}
print("=== cusp 1/2 (c=2, w=3): predicted a_0 = -119/81");
f1(tp)=my(tau=tp/(2*tp+1)); evq(Ps,tau)/(2*tp+1)^8;
f2(tp)=my(tau=tp/(2*tp+1)); evq(Us,tau)/(2*tp+1)^6*(1+evq(xs,tau));
{my(ys=vector(11,k,2.0+0.25*(k-1)));
 print("  a_0        = ",precision(rich(ys,3,f1),30));
 print("  -119/81    = ",precision(-119/81.,30));
 print("  lim U|g(1+x) = ",precision(rich(ys,3,f2),30));
 print("  -119/27      = ",precision(-119/27.,30));}
print();
print("=== cusp 1/3 (c=3, w=4): predicted a_0 = 6273/256");
f3(tp)=my(tau=tp/(3*tp+1)); evq(Ps,tau)/(3*tp+1)^8;
{my(ys=vector(9,k,2.4+0.3*(k-1)));
 print("  a_0     = ",precision(rich(ys,4,f3),25));
 print("  6273/256= ",precision(6273/256.,25));}
quit;
