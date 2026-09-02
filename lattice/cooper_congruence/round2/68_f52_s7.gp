\\ 68_f52_s7.gp -- construct the weight-5/2 Shimura-Borcherds input of Cooper's row s7
\\ from the theta lift:  a(d) = sqrt(d) * sqrt(-3) * Tr_{-3d}(fhat),  fhat = Df + f/(2 pi y),
\\ f = 1/(xF) on Gamma_0(7).  Prints a(d) for every admissible d, for each beta-class.
read("lib.gp"); read("heeg.gp"); read("maass2.gp");
default(realprecision, 50);
N=7; initfser(1,900);
{ trD(D, bt) = my(RF=redforms(D), rep, ch, om, al, t=0., nf=0);
  for(i=1,#RF,
    rep = heegrep(RF[i],N,bt,80);
    if(rep==0, return([0,-1]));
    ch = genchar(rep,-3); om = omeg(rep); nf++;
    al = (-rep[2] + I*sqrt(-D))/(2*rep[1]);
    t += ch*fhatR(1,al)/om);
  [t,nf];
}
{
for(d=1,40,
  my(D=-3*d, bs=List(), r);
  if((D%4)!=0 && ((D-1)%4)!=0, next);
  for(b=0,13, if((b^2-D)%28==0, listput(bs,b)));
  if(#bs==0, print("d=",d,"  (no beta mod 14)"); next);
  print("d=",d,"  D=",D,"  beta-classes ",Vec(bs));
  for(i=1,#bs,
    r = trD(D, Vec(bs)[i]);
    if(r[2]<0, print("     beta=",Vec(bs)[i],"   NOREP"); next);
    print("     beta=",Vec(bs)[i],"  #cls=",r[2],"   a(d)= sqrt(d) sqrt(-3) Tr = ", sqrt(d)*I*sqrt(3)*r[1]));
);
}
quit;
