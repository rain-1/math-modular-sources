default(parisize,"16G");
default(realprecision,120);
M=3000;
et(k)=prod(m=1,M\k+1,1-'q^(k*m)+O('q^(M+2)));
h=et(1)^3*et(4)*et(6)^2/(et(2)^2*et(3)*et(12)^3)/'q;
x=h/((h+3)*(h+4));
E8(k)=1+480*sum(m=1,M\k+1,sigma(m,7)*'q^(k*m))+O('q^(M+2));
Phi=(E8(1)-572*E8(2)+11583*E8(3)-36608*E8(4)+46332*E8(6)-20736*E8(12))/480;
Dx='q*deriv(x,'q);
hs=truncate(h); xs=truncate(x); Ps=truncate(Phi); Ds=truncate(Dx);
evq(s,tau)=subst(s,'q,exp(2*Pi*I*tau));
\\ cusp a/c, gamma=[a,b;c,d]
{foreach([[1,0,2,1],[1,0,3,1],[1,0,4,1],[1,0,6,1],[0,-1,1,0]],G,
  my(a=G[1],b=G[2],c=G[3],d=G[4],w);
  w=12/gcd(c^2,12);
  print("cusp ",a,"/",c,"  width w=",w);
  foreach([1.5,2.0,2.5],Y,
    my(tp=I*Y, tau=(a*tp+b)/(c*tp+d), hh=evq(hs,tau), xx=evq(xs,tau));
    print("   Y'=",Y," Im(tau)=",precision(imag(tau),8)," |q|=",precision(abs(exp(2*Pi*I*tau)),8),
          "  h=",precision(hh,14),"  x=",precision(xx,14))));}
quit;
