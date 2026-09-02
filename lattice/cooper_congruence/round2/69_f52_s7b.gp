\\ 69_f52_s7b.gp -- the integral trace sequence c(d) = sqrt(-3) Tr_{-3d}(fhat) for row s7,
\\ over ALL admissible d (not just squares).  c(m^2) = beta_{s7}(m).
read("lib.gp"); read("heeg.gp"); read("maass2.gp");
default(realprecision, 50);
N=7; initfser(1,900);
{ trD(D, bt) = my(RF=redforms(D), rep, ch, om, al, t=0., nf=0);
  for(i=1,#RF,
    rep = heegrep(RF[i],N,bt,90);
    if(rep==0, return([0,-1]));
    ch = genchar(rep,-3); om = omeg(rep); nf++;
    al = (-rep[2] + I*sqrt(-D))/(2*rep[1]);
    t += ch*fhatR(1,al)/om);
  [t,nf];
}
bet = read("beta_s7.txt");
print("d    c(d) = sqrt(-3) Tr_{-3d}(fhat)   [both beta-classes]     beta(sqrt d) if square");
{
for(d=1,90,
  my(D=-3*d, bs=List(), r, v, sq);
  if((D%4)!=0 && ((D-1)%4)!=0, next);
  for(b=0,13, if((b^2-D)%28==0, listput(bs,b)));
  if(#bs==0, next);
  sq = if(issquare(d), bet[sqrtint(d)], "-");
  v = List();
  for(i=1,#bs,
    r = trD(D, Vec(bs)[i]);
    if(r[2]<0, listput(v,"NOREP"), listput(v, I*sqrt(3)*r[1])));
  print("d=",d,"  betaclasses=",Vec(bs),"   c=",Vec(v),"    beta=",sq);
);
}
quit;
