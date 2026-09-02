\\ 14_ident10.gp -- identification for row s10 (N=10, D_0=-4, disc -4m^2), several variants.
read("lib.gp"); read("heeg.gp"); read("maass2.gp");
default(realprecision, 60);
N=10; initfser(2,700);
{ tr(m, b0, D0) = my(d=-4*m^2, bt, RF, rep, al, ch, om, t=0., nf=0);
  bt = (b0*m)%(2*N);
  RF = redforms(d);
  for(i=1,#RF,
    rep = heegrep(RF[i],N,bt,60);
    if(rep==0, next);
    ch = if(D0==0, 1, genchar(rep,D0)); om = omeg(rep); nf++;
    al = (-rep[2] + I*m*2)/(2*rep[1]);
    t += ch*fhatR(2,al)/om);
  [t,nf,#RF];
}
bet = read("beta_s10.txt");
{
for(m=1,16,
  my(T=tr(m,6,-4), U=tr(m,6,0), V=tr(m,14,-4));
  print("m=",m," cls ",T[3],"  2i*Tr(chi,b=6m)=", 2*I*T[1]);
  print("        2i*Tr(nochi)=", 2*I*U[1], "   2i*Tr(chi,b=-6m)=", 2*I*V[1], "   beta=", bet[m]);
);
}
quit;
