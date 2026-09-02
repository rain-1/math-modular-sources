\\ 09_beta.gp -- Is the vector-valued s7 input symmetric or antisymmetric under beta -> -beta?
\\ Computes the twisted CM trace c(beta,d) for BOTH beta classes at each admissible d.
\\ RUN FROM ../round2 (its libs use relative reads):  gp -q ../round3/09_beta.gp
read("lib.gp"); read("heeg.gp"); read("maass2.gp");
default(realprecision, 50);
N=7; initfser(1,600);
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
{ for(d=1,120,
  my(D=-3*d, bs=List(), r1, r2, b1, b2);
  if((D%4)!=0 && ((D-1)%4)!=0, next);
  for(b=0,13, if((b^2-D)%28==0, listput(bs,b)));
  bs = Vec(bs);
  if(#bs==0, next);
  b1 = bs[1]; b2 = if(#bs>1, bs[2], bs[1]);
  r1 = trD(D,b1); r2 = trD(D,b2);
  if(type(r1)=="t_STR" || type(r2)=="t_STR", print(d," NOREP"); next);
  print("d=",d," betas=",bs," c(b1)=", I*sqrt(3)*r1, "   c(b2)=", I*sqrt(3)*r2);
  listput(L,[d,bs,I*sqrt(3)*r1,I*sqrt(3)*r2]);
); }
write("../round3/09_beta.txt", Vec(L));
quit;
