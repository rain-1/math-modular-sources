\\ 04_gkz7.gp -- GKZ twisted traces on X_0(7): D0=-3, d=m^2, disc -3m^2,
\\ beta = 5m mod 14, genus character chi_{-3}.  Candidate functions 1, u, 1/u.
read("lib.gp"); read("heeg.gp");
default(realprecision, 90);
N = 7; D0 = -3;
{ gkz7(m) = my(d=-3*m^2, bt, RF, rep, al, om, ch, t1, tu, ti, w, nf);
  bt = (5*m)%14;
  RF = redforms(d);
  t1=0.; tu=0.; ti=0.; nf=0;
  for(i=1,#RF,
    rep = heegrep(RF[i],N,bt,60);
    if(rep==0, print("  NOREP m=",m," ",RF[i]); next);
    ch = genchar(rep,D0);
    nf++;
    om = omeg(rep);
    al = (-rep[2] + I*m*sqrt(3))/(2*rep[1]);
    w = (eta(7*al,1)/eta(al,1))^4;
    t1 += ch*1./om; tu += ch*w/om; ti += ch*(1/w)/om);
  [t1,tu,ti,nf];
}
bet = read("beta_s7.txt");
{
for(m=1,16,
  my(A=gkz7(m), th);
  th = A[3] - 49*A[2];
  print("m=",m,"  T1=",A[1],"  TH=",th,"  TH/(I*sqrt(3))=", th/(I*sqrt(3)), "  beta=",bet[m]);
);
}
quit;
