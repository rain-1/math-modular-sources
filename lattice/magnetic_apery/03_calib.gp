default(parisizemax,12000000000);
read("lib.gp");
NQ = 405; M3 = 400;
{
for(hh=1,#HOSTS,
  my(HH, N, C, B, dv, rv, us, Fs, F2, qv, Xs, T, L1, L2, L3);
  HH = HOSTS[hh]; N=HH[1]; C=HH[2]; B=HH[3]; dv=HH[4]; rv=HH[5];
  us = useries(dv,rv,NQ); Fs = Fseries(dv,rv,NQ); F2 = Fs^2;
  qv = [B^2, 2*C*B, C^2];
  Xs = xibasis(C, us, F2, qv, M3);
  T = xitomat(Xs, M3);
  L1 = latmag(matrix(1,100,j,m,T[j,m]));
  L2 = latmag(matrix(1,200,j,m,T[j,m]));
  L3 = latmag(T);
  print(HH[6], "   own-Q generators at M=100,200,400: ", L1[1,1], " ", L2[1,1], " ", L3[1,1]);
);
}
quit;
