\\ 05_wtrace7.gp -- A-weighted twisted traces on X_0(7); test beta = l0*S0 + l1*S1/m
read("lib.gp"); read("heeg.gp");
default(realprecision, 90);
N = 7; D0 = -3;
{ tr(m) = my(d=-3*m^2, bt, RF, rep, al, ch, om, w, h, s0=0., s1=0., s2=0.);
  bt = (5*m)%14;
  RF = redforms(d);
  for(i=1,#RF,
    rep = heegrep(RF[i],N,bt,60);
    if(rep==0, next);
    ch = genchar(rep,D0); om = omeg(rep);
    al = (-rep[2] + I*m*sqrt(3))/(2*rep[1]);
    w = (eta(7*al,1)/eta(al,1))^4;
    h = 1/w - 49*w;
    s0 += ch*h/om;
    s1 += ch*rep[1]*h/om;
    s2 += ch*rep[1]^2*h/om);
  [s0,s1,s2];
}
bet = read("beta_s7.txt");
S = sqrt(3)*I;
{
for(m=1,20,
  my(T=tr(m));
  print("m=",m,"  s0=",T[1]/S,"  s1=",T[2]/S,"  s2=",T[3]/S,"  beta=",bet[m]);
);
}
quit;
