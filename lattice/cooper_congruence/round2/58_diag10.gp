\\ 58_diag10.gp -- per-form diagnostic for row s10: content, genus character, Im(alpha), term.
read("50_lib.gp"); read("heeg.gp"); read("maass2.gp");
default(realprecision, 60);
N=10; initfser(2,700);
bet = read("20_beta_s10.txt");
cont(Q) = gcd(gcd(Q[1],Q[2]),Q[3]);
{
for(m=1,16,
  if(m%5==0, next);
  my(d=-4*m^2, bt=(6*m)%20, RF=redforms(d), rep, ch, om, al, v, T=0., Tz=0., nz=0);
  print("=== m=",m,"  d=",d,"  beta=",bt,"  #forms=",#RF,"  beta_row=",bet[m]);
  for(i=1,#RF,
    rep = heegrep(RF[i],N,bt,60);
    if(rep==0, print("   ",RF[i]," NOREP"); next);
    ch = genchar(rep,-4); om = omeg(rep);
    al = (-rep[2] + I*m*2)/(2*rep[1]);
    v = fhatR(2,al);
    if(ch==0, nz++; Tz += v/om, T += ch*v/om);
    print("   ",RF[i]," cont=",cont(RF[i])," -> ",rep," cont=",cont(rep),"  chi=",ch,"  om=",om,"  Im(al)=",imag(al),"  |fhat|=",abs(v)));
  print("   T(chi) = ",T,"   i*T = ",I*T,"   2i*T = ",2*I*T,"   beta = ",bet[m]);
  print("   #chi=0 forms: ",nz,"   their sum Tz = ",Tz,"   i*(T+Tz)=",I*(T+Tz),"  i*(T-Tz)=",I*(T-Tz));
);
}
quit;
