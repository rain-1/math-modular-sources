default(parisize,"16G");
default(realprecision,150);
M=4000;
et(k)=prod(m=1,M\k+1,1-'q^(k*m)+O('q^(M+2)));
h=et(1)^3*et(4)*et(6)^2/(et(2)^2*et(3)*et(12)^3)/'q;
x=h/((h+3)*(h+4));
E8(k)=1+480*sum(m=1,M\k+1,sigma(m,7)*'q^(k*m))+O('q^(M+2));
Phi=(E8(1)-572*E8(2)+11583*E8(3)-36608*E8(4)+46332*E8(6)-20736*E8(12))/480;
U=Phi/('q*deriv(x,'q));
xs=truncate(x); Ps=truncate(Phi); Us=truncate(U);
evq(s,tau)=subst(s,'q,exp(2*Pi*I*tau));
a=1;b=0;c=2;d=1;w=3;
print("cusp 1/2, gamma=[1,0;2,1], width 3");
print("predicted a_0(Phi) = -119/81 = ",precision(-119/81.,20));
print("predicted w*a_0 = -119/27 = ",precision(-119/27.,20));
{foreach([1.6,1.9,2.2,2.5,2.8],Y,
  my(tp=I*Y, tau=tp/(2*tp+1), jf=(2*tp+1), qp=exp(2*Pi*I*tp/w));
  my(P8=evq(Ps,tau)/jf^8, U6=evq(Us,tau)/jf^6, xx=evq(xs,tau));
  print("  Y'=",Y,"  q'=",precision(abs(qp),6),
        "  Phi|_8g = ",precision(real(P8),18),
        "   U|_6g*(1+x) = ",precision(real(U6*(1+xx)),18)));}
print();
print("Richardson in q' on Phi|_8gamma:");
{my(ys=[2.0,2.3,2.6,2.9,3.2,3.5],m=6,V,rhs);
 V=matrix(m,m,i,j,exp(2*Pi*I*(I*ys[i])/w)^(j-1));
 rhs=vector(m,i,my(tp=I*ys[i],tau=tp/(2*tp+1)); evq(Ps,tau)/(2*tp+1)^8)~;
 print("  a_0 = ",precision(matsolve(V,rhs)[1],25));}
{my(ys=[2.0,2.3,2.6,2.9,3.2,3.5],m=6,V,rhs);
 V=matrix(m,m,i,j,exp(2*Pi*I*(I*ys[i])/w)^(j-1));
 rhs=vector(m,i,my(tp=I*ys[i],tau=tp/(2*tp+1)); evq(Us,tau)/(2*tp+1)^6*(1+evq(xs,tau)))~;
 print("  lim U|_6g*(1+x) = ",precision(matsolve(V,rhs)[1],25));}
print("alpha = a_0^2 w^8/7! = ",precision((119/81)^2*3^8/5040.,25),"   2023/720 = ",precision(2023/720.,25));
quit;
