\\ 80_s10min.gp -- independent s10 check with MINIMAL-A representatives.
read("lib.gp"); read("heeg.gp"); read("heegmin.gp"); read("maass2.gp");
default(realprecision, 60);
N=10; initfser(2,900);
{ trD(D, bt) = my(RF=redforms(D), rep, ch, om, al, t=0., mnIm=9.);
  for(i=1,#RF,
    rep = heegminA(RF[i],N,bt,60);
    if(rep==0, return(["NOREP",0]));
    ch = genchar(rep,-4); om = omeg(rep);
    al = (-rep[2] + I*sqrt(-D))/(2*rep[1]);
    if(imag(al)<mnIm, mnIm=imag(al));
    t += ch*fhatR(2,al)/om);
  [t,mnIm];
}
bet = read("beta_s10.txt");
{ my(tot=0, bad=List(), deg=List(), zer=List(), sq=List(), worst=9.);
  for(d=1,700,
    if(d%4!=0 && d%4!=1, next);
    my(D=-4*d, bs=List());
    for(b=0,2*N-1, if((b^2-D)%(4*N)==0, listput(bs,b)));
    if(#bs==0, next);
    my(bt = if(issquare(d), (6*sqrtint(d))%20, Vec(bs)[1]), r, v);
    if(!setsearch(Set(Vec(bs)), bt), bt = Vec(bs)[1]);
    r = trD(D, bt);
    if(type(r[1])=="t_STR", print("NOREP at d=",d); next);
    if(r[2]<worst, worst=r[2]);
    v = 2*I*r[1];
    if(abs(imag(v)) > 1e-25*(1+abs(v)), listput(deg,d); next);
    v = round(real(v));
    tot++;
    if(v==0, listput(zer,d));
    if(v % d != 0, listput(bad,[d,v]));
    if(issquare(d), listput(sq,[sqrtint(d), v, bet[sqrtint(d)]])));
  print("s10 (minimal-A): ", tot, " non-degenerate admissible d <= 700;  min Im(alpha) over all = ", worst);
  print("  d | c(d) failures: ", #bad, "  ", if(#bad>0, Str(Vec(bad)), ""));
  print("  degenerate d: ", Vec(deg), "   all 25|d? ", if(#select(x->x%25!=0,Vec(deg))==0,"YES","no"));
  print("  zeros: ", Vec(zer));
  print("   zero set = {5|d} union {d = 5 mod 8}? ",
    if(#select(x->x%5!=0 && x%8!=5, Vec(zer))==0, "every zero is 5|d or 5 mod 8: YES", "NO"));
  my(pred=List()); for(d=1,700, if((d%4==0||d%4==1) && (d%5==0||d%5==1||d%5==4) && (d%5==0||d%8==5), listput(pred,d)));
  print("   predicted zeros (admissible with 5|d or d=5 mod 8): ", #pred, "   actual zeros: ", #zer);
  print("   predicted-but-not-zero: ", select(x->!setsearch(Set(Vec(zer)),x) && !setsearch(Set(Vec(deg)),x), Vec(pred)));
  print("  squares [m, c(m^2), beta(m), ratio]: ", vector(#sq,i,[Vec(sq)[i][1],Vec(sq)[i][2],Vec(sq)[i][3],if(Vec(sq)[i][3]!=0,Vec(sq)[i][2]/Vec(sq)[i][3],"-")]));
}
quit;
