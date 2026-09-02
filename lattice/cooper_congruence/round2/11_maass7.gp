\\ 11_maass7.gp -- THE TEST: twisted CM trace of the weight-0 invariant
\\   fhat = D f + f/(2 pi y),   f = 1/(x F)  weakly holomorphic of weight -2 on Gamma_0(7),
\\ over the Heegner points of discriminant -3m^2, against beta_{s7}(m).
read("lib.gp"); read("heeg.gp"); read("wt2.gp");
default(realprecision, 60);
N=7; D0=-3;
{ tr(m) = my(d=-3*m^2, bt, RF, rep, al, ch, om, t=0., nf=0);
  bt = (5*m)%14;
  RF = redforms(d);
  for(i=1,#RF,
    rep = heegrep(RF[i],N,bt,60);
    if(rep==0, next);
    ch = genchar(rep,D0); om = omeg(rep); nf++;
    al = (-rep[2] + I*m*sqrt(3))/(2*rep[1]);
    t += ch*fhat7(al)[2]/om);
  [t,nf];
}
bet = read("beta_s7.txt");
{
for(m=1,18,
  my(T=tr(m));
  print("m=",m,"  #",T[2],"  T=",T[1],"  beta=",bet[m],"  beta/T=",bet[m]/T[1]);
);
}
quit;
