\\ 70_cd_s7.gp -- the integer sequence c(d) = sqrt(-3) Tr_{-3d}(fhat) for s7, d <= 400,
\\ written to 71_cd_s7_big.txt as a list of [d, c(d)] with the beta-class congruent to
\\ 5*sqrt(d) mod 14 on squares (and the first listed class otherwise).
read("lib.gp"); read("heeg.gp"); read("maass2.gp");
default(realprecision, 60);
N=7; initfser(1,1600);
{ trD(D, bt) = my(RF=redforms(D), rep, ch, om, al, t=0.);
  for(i=1,#RF,
    rep = heegrep(RF[i],N,bt,120);
    if(rep==0, return("NOREP"));
    ch = genchar(rep,-3); om = omeg(rep);
    al = (-rep[2] + I*sqrt(-D))/(2*rep[1]);
    t += ch*fhatR(1,al)/om);
  t;
}
L = List();
{
for(d=1,1600,
  my(D=-3*d, bs=List(), r, v, bt);
  if((D%4)!=0 && ((D-1)%4)!=0, next);
  for(b=0,13, if((b^2-D)%28==0, listput(bs,b)));
  if(#bs==0, next);
  bt = if(issquare(d), (5*sqrtint(d))%14, Vec(bs)[1]);
  if(!setsearch(Set(Vec(bs)), bt), bt = Vec(bs)[1]);
  r = trD(D, bt);
  if(type(r)=="t_STR", listput(L,[d,"NOREP"]); next);
  v = I*sqrt(3)*r;
  if(abs(imag(v)) > 1e-25*(1+abs(v)), listput(L,[d,"DEGEN"]), listput(L,[d,round(real(v))]));
);
}
write("71_cd_s7_big.txt", Vec(L));
print("wrote ", #L, " entries; first 40:");
print(vector(min(40,#L), i, Vec(L)[i]));
quit;
