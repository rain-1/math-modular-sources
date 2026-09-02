\\ 79_zero21.gp -- diagnostic: why is c_{s10}(21) = 0 in my pipeline?  Print the terms.
read("lib.gp"); read("heeg.gp"); read("maass2.gp");
default(realprecision, 40);
N=10; initfser(2,900);
{ show(d, bt) = my(D=-4*d, RF=redforms(D), rep, ch, om, al, v, t=0.);
  print("d=",d,"  D=",D,"  beta=",bt,"  #SL2 classes=",#RF);
  for(i=1,#RF,
    rep = heegrep(RF[i],N,bt,90);
    if(rep==0, print("   ",RF[i]," NOREP"); next);
    ch = genchar(rep,-4); om = omeg(rep);
    al = (-rep[2] + I*sqrt(-D))/(2*rep[1]);
    v = fhatR(2,al);
    t += ch*v/om;
    print("   ",RF[i]," -> A=",rep[1]," B=",rep[2]," chi=",ch," om=",om,"  alpha=",al,"  fhat=",v));
  print("   SUM = ", t, "     2i*SUM = ", 2*I*t);
}
show(21,6); show(21,14); show(29,2); show(1,6); show(41,6);
quit;
