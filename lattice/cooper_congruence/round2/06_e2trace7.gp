\\ 06_e2trace7.gp -- twisted CM traces on X_0(7) of the real-analytic weight-0
\\ invariants H and (E_2^*/F)*H, H = 1/u - 49u.  Fit beta against them.
read("lib.gp"); read("heeg.gp"); read("e2.gp");
default(realprecision, 60);
N=7; D0=-3;
{ tr(m) = my(d=-3*m^2, bt, RF, rep, al, ch, om, w, h, Fv, t0=0., t1=0.);
  bt = (5*m)%14;
  RF = redforms(d);
  for(i=1,#RF,
    rep = heegrep(RF[i],N,bt,60);
    if(rep==0, next);
    ch = genchar(rep,D0); om = omeg(rep);
    al = (-rep[2] + I*m*sqrt(3))/(2*rep[1]);
    w = (eta(7*al,1)/eta(al,1))^4;
    h = 1/w - 49*w;
    Fv = (7*E2star(7*al) - E2star(al))/6;
    t0 += ch*h/om;
    t1 += ch*h*E2star(al)/Fv/om);
  [t0,t1];
}
bet = read("beta_s7.txt");
S = sqrt(3)*I;
{
for(m=1,14,
  my(T=tr(m));
  print("m=",m,"  t0=",T[1]/S,"   t1=",T[2]/S,"   beta=",bet[m]);
);
}
quit;
