default(parisizemax,12000000000);
read("lib.gp");
default(realprecision, 60);
NQ = 405;
print("### the twelve Fricke hosts: is the canonical source Phi_0 = F*Dx holomorphic?");
print("R = limsup |c(m)|^{1/m};  R=1 <=> Phi_0 holomorphic on H;  R>1 <=> double pole at tau_0, Im tau_0 = log R/(2pi)");
print("");
print("tag | C | B | roots of 1+Bu+Cu^2 | R (m=380..400) | Im tau_0 | free int (magnetic)?");
{
for(hi=1, #HOSTS,
  my(HH, N, C, B, dv, rv, us, Fs, F2, Phi, R, rts, MAG, T, L1, L2, Xs, disc);
  HH = HOSTS[hi]; N=HH[1]; C=HH[2]; B=HH[3]; dv=HH[4]; rv=HH[5];
  us = useries(dv,rv,NQ); Fs = Fseries(dv,rv,NQ); F2 = Fs^2;
  Phi = F2*us*(1-C*us^2)/(1+B*us+C*us^2)^2;
  R = 0;
  for(m=380, 400, my(t); t = abs(polcoeff(Phi,m)*1.0)^(1/m); if(t>R, R=t));
  Xs = xibasis(C, us, F2, [B^2, 2*C*B, C^2], 400);
  T = xitomat(Xs, 400);
  L1 = latmag(matrix(1,100,j,m,T[j,m])); L2 = latmag(T);
  MAG = if(abs(L1[1,1])==abs(L2[1,1]), "MAGNETIC", "not magnetic");
  disc = B^2-4*C;
  rts = if(disc>=0, Str("real: ", (-B+sqrt(disc*1.0))/(2*C), ", ", (-B-sqrt(disc*1.0))/(2*C)), Str("complex, disc=", disc));
  print(HH[6], " | ", C, " | ", B, " | ", rts, " | R=", R, " | Im=", log(R)/(2*Pi), " | ", MAG);
);
}
quit;
