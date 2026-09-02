\\ 12_vv2.gp -- systematic vector-valued S-transformation test.
\\ Components beta=1..6 (antisymmetric; beta=0,7 vanish), exponent EF*d/28,
\\ coefficients d^s * c(d), weight k, unknown polar coefficients at d in NEG.
\\ Law:  G_b(-1/tau) = C (tau/i)^k sum_{g=1}^{6} sin(pi b g/7) G_g(tau),  |C| should be 1/sqrt(14).
default(realprecision, 60);
DMAX = 2500;
DAT = read("../round2/73_cd_s7_2500.txt");
CD = vector(DMAX); BT = vector(DMAX);
{ for(i=1,#DAT, my(e=DAT[i], d=e[1], v=e[2], bs=List(), bt);
    if(d>DMAX, next);
    for(b=0,13, if((b^2+3*d)%28==0, listput(bs,b)));
    bs = Vec(bs);
    bt = if(issquare(d), (5*sqrtint(d))%14, bs[1]);
    if(!setsearch(Set(bs), bt), bt = bs[1]);
    BT[d] = bt; CD[d] = if(type(v)=="t_STR", 0, v)); }
NEG = [-3,-12,-19,-20,-24,-27];
TAUS = [0.3+0.6*I, -0.2+0.55*I, 0.41+0.72*I, 0.13+0.83*I];
{ lab(EF,d) = if(EF==1, (3*BT[d])%14, BT[d]); }
{ labneg(EF,d) = my(bt=0); for(b=0,13, if((b^2+3*d)%28==0, bt=b; break));
   if(EF==1, (3*bt)%14, bt); }
{ ampl(beta, lb) = if(lb==beta, 1, if((14-lb)%14==beta, -1, 0)); }
{ Aser(beta,tau,EF,s) = my(sm=0., w=exp(2*Pi*I*EF*tau/28), e);
   for(d=1,DMAX, if(CD[d]!=0, e=ampl(beta,lab(EF,d)); if(e!=0, sm += e*d^s*CD[d]*w^d))); sm; }
{ Bser(beta,tau,EF,dneg) = my(w=exp(2*Pi*I*EF*tau/28), e=ampl(beta,labneg(EF,dneg)));
   if(e==0, 0., e*w^dneg); }
{ SS(b,g) = sin(Pi*b*g/7); }
{ run(EF,k,s) =
  my(nn=#NEG, nu=2*nn+1, nr=6*#TAUS, M=matrix(nr,nu), Y=vector(nr), V, res, C, ok);
  for(ti=1,#TAUS, my(tau=TAUS[ti], t2=-1/tau, pw=exp(k*log(tau/I)));
    for(b=1,6, my(row=6*(ti-1)+b, sa);
      sa = sum(g=1,6, SS(b,g)*Aser(g,tau,EF,s));
      for(i=1,nn, M[row,i] = Bser(b,t2,EF,NEG[i]));
      M[row,nn+1] = -pw*sa;
      for(i=1,nn, M[row,nn+1+i] = -pw*sum(g=1,6, SS(b,g)*Bser(g,tau,EF,NEG[i])));
      Y[row] = -Aser(b,t2,EF,s)));
  V = matsolve(M[1..nu,], Y[1..nu]~);
  res = M*V - Y~;
  C = V[nn+1];
  ok = vector(nn, i, abs(V[nn+1+i] - C*V[i]));
  print("EF=",EF," k=",k," s=",s,"  maxres=", vecmax(vector(nr,i,abs(res[i]))),
        "  |C|=", abs(C), "  (1/sqrt14=", 1/sqrt(14), ")  max|w-Cp|=", vecmax(ok));
}
{ foreach([1,3], EF, foreach([1/2,3/2,5/2,7/2], k, foreach([-1/2,0,1/2,1], s, run(EF,k,s)))); }
quit;
