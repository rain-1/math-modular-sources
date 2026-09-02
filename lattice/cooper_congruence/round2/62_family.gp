\\ 62_family.gp -- is the master congruence a GENERAL property of weight-4 meromorphic
\\ modular forms with a double pole on a Heegner divisor?  Family: for a discriminant D<0
\\ with Hilbert class polynomial H_D, put  Phi_D = E_4 / H_D(j)^2  (weight 4, holomorphic
\\ at the cusp, double poles exactly at the disc-D CM points).  Test magnetism and
\\ n^2 | beta(n) with psi = chi_D, after the natural integral normalisation.
M = 260;
E4s = 1 + 240*sum(n=1,M,sigma(n,3)*'q^n) + O('q^(M+1));
E6s = 1 - 504*sum(n=1,M,sigma(n,5)*'q^n) + O('q^(M+1));
Dl  = 'q*prod(n=1,M,(1-'q^n+O('q^(M+1)))^24);
js  = E4s^3/Dl;
{ testD(D) = my(H, Ph, c, g, cp, b, sc, bad, s);
  H = polclass(D);
  Ph = E4s/subst(H,variable(H),js)^2;
  c = vector(M,m,polcoeff(Ph,m));
  g = content(vector(M,m,c[m]));       \\ clear the common denominator/content
  sc = 1/g;
  c = vector(M,m,sc*c[m]);
  \\ magnetism
  bad = 0; for(m=1,M, if(c[m]%m!=0, bad=m; break));
  print("D=",D,"  h=",poldegree(H),"  scale=",sc,"   c(1..4)=",vector(4,i,c[i]));
  print("    m | c(m) for m<=",M,"? ", if(bad==0,"YES",concat("NO at ",bad)));
  if(bad!=0, return(0));
  cp = vector(M,m,c[m]/m);
  b = vector(M);
  for(n=1,M, my(t=0); fordiv(n,d, t += moebius(d)*kronecker(D,d)*cp[n/d]); b[n]=t);
  bad = List(); for(n=1,M, if(b[n]%(n^2)!=0, listput(bad,n)));
  print("    n^2 | beta(n) with psi=chi_",D,"? ", if(#bad==0,"YES",concat(concat("NO, ",#bad)," failures: ")));
  if(#bad>0, print("       ",vector(min(12,#bad),i,Vec(bad)[i])));
  s=""; forprime(p=2,23, my(mn=999); for(n=1,M, if(n%p==0, my(v=valuation(b[n],p)-2*valuation(n,p)); if(v<mn,mn=v))); s=concat(s,concat(concat(" ",p),concat(":",mn))));
  print("    min v_p(beta(n))-2v_p(n): ",s);
  print("    beta(1..6)=",vector(6,i,b[i]));
}
{ for(i=1,8, my(DD=[-7,-8,-11,-19,-43,-67,-163,-15][i]); testD(DD)); }
quit;
