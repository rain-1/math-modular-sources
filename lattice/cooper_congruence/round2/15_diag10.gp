read("lib.gp"); read("heeg.gp"); read("maass2.gp");
default(realprecision, 30);
N=10; initfser(2,700);
{
for(m=2,6,
  my(d=-4*m^2, bt=(6*m)%20, RF=redforms(d), rep, ch, om, al, v);
  print("=== m=",m,"  d=",d,"  beta=",bt);
  for(i=1,#RF,
    rep = heegrep(RF[i],N,bt,60);
    if(rep==0, print("   ",RF[i]," NOREP"); next);
    ch = genchar(rep,-4); om = omeg(rep);
    al = (-rep[2] + I*m*2)/(2*rep[1]);
    v = fhatR(2,al);
    print("   ",RF[i]," -> ",rep,"  chi=",ch,"  om=",om,"  fhat=",v, "  term=",ch*v/om));
);
}
quit;
