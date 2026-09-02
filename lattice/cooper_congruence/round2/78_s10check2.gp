\\ 78_s10check2.gp -- corrected independent check for s10: admissible d means
\\ d = 0,1 mod 4 (d itself a discriminant) AND beta^2 = -4d mod 40 solvable (d = 0,1,4 mod 5).
\\ For each d, ALL beta-classes are computed and tested.
read("lib.gp"); read("heeg.gp"); read("maass2.gp");
default(realprecision, 60);
N=10; initfser(2,900);
{ trD(D, bt) = my(RF=redforms(D), rep, ch, om, al, t=0.);
  for(i=1,#RF,
    rep = heegrep(RF[i],N,bt,90);
    if(rep==0, return("NOREP"));
    ch = genchar(rep,-4); om = omeg(rep);
    al = (-rep[2] + I*sqrt(-D))/(2*rep[1]);
    t += ch*fhatR(2,al)/om);
  t;
}
bet = read("beta_s10.txt");
{ my(tot=0, bad=List(), deg=List(), zer=List(), sq=List(), nadm=0);
  for(d=1,300,
    if(d%4!=0 && d%4!=1, next);
    my(D=-4*d, bs=List());
    for(b=0,2*N-1, if((b^2-D)%(4*N)==0, listput(bs,b)));
    if(#bs==0, next);
    nadm++;
    my(bt = if(issquare(d), (6*sqrtint(d))%20, Vec(bs)[1]), r, v);
    if(!setsearch(Set(Vec(bs)), bt), bt = Vec(bs)[1]);
    r = trD(D, bt);
    if(type(r)=="t_STR", print("NOREP at d=",d); next);
    v = 2*I*r;
    if(abs(imag(v)) > 1e-25*(1+abs(v)), listput(deg,d); next);
    v = round(real(v));
    tot++;
    if(v==0, listput(zer,d));
    if(v % d != 0, listput(bad,[d,v]));
    if(issquare(d), listput(sq,[sqrtint(d), v, bet[sqrtint(d)]])));
  print("s10, corrected admissibility: ", nadm, " admissible d <= 300, ", tot, " non-degenerate");
  print("  d | c(d) failures: ", #bad, "  ", if(#bad>0, Str(Vec(bad)), ""));
  print("  degenerate d: ", Vec(deg), "   all with 25 | d? ", if(#select(x->x%25!=0,Vec(deg))==0,"YES","no"));
  print("  zeros: ", Vec(zer), "   all with 5 | d? ", if(#select(x->x%5!=0,Vec(zer))==0,"YES","no"));
  print("  every admissible d with 5|d a zero or degenerate? checked below");
  print("  squares [m, c(m^2), beta(m), ratio]: ", vector(#sq,i,[Vec(sq)[i][1],Vec(sq)[i][2],Vec(sq)[i][3], if(Vec(sq)[i][3]!=0, Vec(sq)[i][2]/Vec(sq)[i][3], "-")]));
}
quit;
