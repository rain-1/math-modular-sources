\\ 17_ident7wide.gp -- the s7 identification, m <= 45, 40 digits, with class-count checks.
read("lib.gp"); read("heeg.gp"); read("maass2.gp");
default(realprecision, 60);
N=7; initfser(1,900);
{ tr(m) = my(d=-3*m^2, bt, RF, rep, al, ch, om, t=0., nf=0);
  bt = (5*m)%14;
  RF = redforms(d);
  for(i=1,#RF,
    rep = heegrep(RF[i],N,bt,80);
    if(rep==0, print("  NOREP m=",m," ",RF[i]); next);
    ch = genchar(rep,-3); om = omeg(rep); nf++;
    al = (-rep[2] + I*m*sqrt(3))/(2*rep[1]);
    t += ch*fhatR(1,al)/om);
  [t,nf,#RF];
}
bet = read("beta_s7.txt");
print("f = 1/(xF) q-expansion, coefficients of q^-1 .. q^8:");
print(vector(10,i,polcoeff(FSER[1][1],i-2)));
print("");
bad = 0; worst = 0.;
{
for(m=1,45,
  my(T, P, e);
  if(m%7==0, print("m=",m,"  [7|m: beta=",bet[m],", trace is real, skipped]"); next);
  T = tr(m);
  if(T[2]!=T[3], print("  CLASS COUNT MISMATCH at m=",m,": ",T[2],"/",T[3]); bad++);
  P = I*sqrt(3)*T[1];
  e = abs(P-bet[m]);
  if(e > worst, worst = e);
  print("m=",m,"  cls=",T[3],"  i sqrt3 Tr = ",P,"   beta = ",bet[m],"   |err| = ",e);
);
}
print("");
print("class-count mismatches: ", bad, "    worst |i sqrt3 Tr - beta| over m<=45, 7 nmid m: ", worst);
quit;
