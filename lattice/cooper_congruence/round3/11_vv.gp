\\ 11_vv.gp -- the vector-valued transformation law.
\\ Components beta in Z/14, exponents d/28 with d = beta^2 (28); a(-beta,d) = -a(beta,d)
\\ (verified in 09_beta.log), so the beta=0,7 components vanish identically.
\\ Test  G_beta(-1/tau) = C tau^{5/2} sum_gamma sin(pi beta gamma/7) G_gamma(tau)
\\ with unknowns C and the polar coefficient p at d=-3 (component beta=5).
default(realprecision, 60);
DMAX = 2500;
DAT = read("../round2/73_cd_s7_2500.txt");
CD = vector(DMAX); BU = vector(DMAX);
{ for(i=1,#DAT, my(e=DAT[i], d=e[1], v=e[2], bs=List(), bt);
    if(d>DMAX, next);
    for(b=0,13, if((b^2+3*d)%28==0, listput(bs,b)));
    bs = Vec(bs);
    bt = if(issquare(d), (5*sqrtint(d))%14, bs[1]);
    if(!setsearch(Set(bs), bt), bt = bs[1]);
    BU[d] = (3*bt)%14;
    CD[d] = if(type(v)=="t_STR", 0, v)); }
\\ a(beta,d) for beta=1..6
{ aa(beta,d) = if(d<1 || d>DMAX || CD[d]==0, return(0));
  if(BU[d]==beta, return(CD[d]));
  if((14-BU[d])%14==beta, return(-CD[d]));
  0; }
print("check: d=1 -> beta=",BU[1]," a=",CD[1],"   d=4 -> beta=",BU[4]," a=",CD[4]);
\\ G_beta(tau) = sum_{d>0} a(beta,d) e(d tau/28)   (known part A)
{ Aser(beta,tau) = my(s=0., w=exp(2*Pi*I*tau/28)); for(d=1,DMAX, if(CD[d]!=0, s += aa(beta,d)*w^d)); s; }
\\ polar basis vector: d=-3 sits in components beta = 5 and 9 (5^2=25=-3 mod 28)
{ Bser(beta,tau) = my(w=exp(2*Pi*I*tau/28));
  if(beta==5, return(w^(-3))); if(beta==9, return(-w^(-3))); 0; }
{ S(b,g) = sin(Pi*b*g/7); }
{ test(tau) = my(t2=-1/tau, M, Y, V, res);
  M = matrix(6,3); Y = vector(6);
  for(b=1,6,
    my(sa = sum(g=1,6, S(b,g)*Aser(g,tau)), sb = sum(g=1,6, S(b,g)*Bser(g,tau)), pw = exp(5/2*log(tau)));
    M[b,1] = Bser(b,t2);       \\ coefficient of p
    M[b,2] = -pw*sa;           \\ coefficient of u = C
    M[b,3] = -pw*sb;           \\ coefficient of v = C p
    Y[b] = -Aser(b,t2));
  my(RR=[5,1,2], MM=matrix(3,3,i,j,M[RR[i],j]), YY=vector(3,i,Y[RR[i]]));
  V = matsolve(MM, YY~);
  res = M*V - Y~;
  print("tau=",tau);
  print("   p = ", V[1]);
  print("   C = ", V[2]);
  print("   v/u = ", V[3]/V[2], "   (should equal p)");
  print("   residual (all 6 eqs, rel): ", vector(6,i,abs(res[i])/(1e-99+abs(Y[i]))));
}
test(0.3+0.6*I);
test(-0.2+0.55*I);
test(0.41+0.72*I);
quit;
