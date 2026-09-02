default(parisizemax,12000000000);
read("lib.gp");
NQ = 205;
LIM = 200;
{
for(hh=1, #HOSTS,
  my(H); my(N); my(C); my(B); my(dv); my(rv); my(tag);
  H = HOSTS[hh]; N = H[1]; C = H[2]; B = H[3]; dv = H[4]; rv = H[5]; tag = H[6];
  if(N!=5 && N!=6, next);
  my(us); my(Fs); my(F2);
  us = useries(dv,rv,NQ); Fs = Fseries(dv,rv,NQ); F2 = Fs^2;
  print("### host ", tag, "  N=",N," C=",C," B=",B);
  for(Bp = -60, 60,
    my(den); my(rho); my(Phi); my(bi); my(bm);
    den = (1 + Bp*us + C*us^2)^2;
    rho = us*(1 - C*us^2)/den;
    Phi = F2*rho;
    bi = intfail(Phi, LIM); bm = magfail(Phi, LIM);
    if(bi==0, print("  Bp=",Bp,"  integral  magfail=",bm, "   c(1..8)=", vector(8,j,polcoeff(Phi,j))));
    if(bi!=0 && bi>3, print("  Bp=",Bp,"  nonintegral first at m=",bi));
  );
);
}
quit;
