\\ 02_trace7.gp -- is beta_{s7}(m) a twisted trace of H = 1/u - 49u over
\\ Heegner points of discriminant -3m^2 on X_0(7)?
read("lib.gp"); read("heeg.gp");
default(realprecision, 80);
N = 7; C7 = 49;
bet = read("beta_s7.txt");
uval(t) = (eta(7*t,1)/eta(t,1))^4;
Hval(t) = my(w=uval(t)); 1/w - C7*w;
{
for(m=1,18,
  my(d=-3*m^2, beta, RF, T, S, nf, rep, al, om);
  beta = (5*m)%14;
  RF = redforms(d);
  T = 0.; nf = 0;
  for(i=1,#RF,
    rep = heegrep(RF[i],N,beta,40);
    if(rep==0, print("  m=",m," NO REP for ",RF[i]); next);
    nf++;
    om = omeg(rep);
    al = (-rep[2] + I*m*sqrt(3))/(2*rep[1]);
    T += Hval(al)/om);
  print("m=",m,"  #classes=",#RF,"  #found=",nf,"  T=",T);
  print("     beta=",bet[m],"   beta/T = ", bet[m]/T);
);
}
quit;
