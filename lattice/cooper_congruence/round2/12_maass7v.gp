\\ 12_maass7v.gp -- variants of the Maass-raised CM trace on X_0(7)
read("lib.gp"); read("heeg.gp"); read("wt2.gp");
default(realprecision, 50);
N=7; D0=-3;
{ tr(m) = my(d=-3*m^2, bt, RF, rep, al, ch, om, cont, v, t=vector(4,i,0.));
  bt = (5*m)%14;
  RF = redforms(d);
  for(i=1,#RF,
    cont = gcd(gcd(RF[i][1],RF[i][2]),RF[i][3]);
    rep = heegrep(RF[i],N,bt,60);
    if(rep==0, next);
    ch = genchar(rep,D0); om = omeg(rep);
    al = (-rep[2] + I*m*sqrt(3))/(2*rep[1]);
    v = fhat7(al)[2]/om;
    t[1] += ch*v;
    t[3] += v;
    if(cont==1, t[2] += ch*v; t[4] += v));
  t;
}
bet = read("beta_s7.txt");
cp  = read("cp_s7.txt");
S = I*sqrt(3);
{
for(m=1,16,
  my(T=tr(m));
  print("m=",m,"  allchi=",T[1]*S,"  primchi=",T[2]*S,"  all=",T[3]*S,"  prim=",T[4]*S);
  print("        beta=",bet[m],"   c'=",cp[m]);
);
}
quit;
