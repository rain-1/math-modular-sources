\\ 77_s10check.gp -- INDEPENDENT check (my own heeg.gp/maass2.gp pipeline, which matched the
\\ Gamma_0(10) class counts) of the claim  d | c(d)  for row s10 with c(d) = 2i T(d).
read("lib.gp"); read("heeg.gp"); read("maass2.gp");
default(realprecision, 60);
N=10; initfser(2,900);
{ trD(D, bt) = my(RF=redforms(D), rep, ch, om, al, t=0., nf=0);
  for(i=1,#RF,
    rep = heegrep(RF[i],N,bt,90);
    if(rep==0, return(["NOREP",0]));
    ch = genchar(rep,-4); om = omeg(rep); nf++;
    al = (-rep[2] + I*sqrt(-D))/(2*rep[1]);
    t += ch*fhatR(2,al)/om);
  [t,nf];
}
bet = read("beta_s10.txt");
{ my(tot=0, bad=List(), deg=List(), zer=List(), sqok=List());
  for(d=1,300,
    my(D=-4*d, bs=List(), r, v);
    if((D%4)!=0 && ((D-1)%4)!=0, next);
    for(b=0,2*N-1, if((b^2-D)%(4*N)==0, listput(bs,b)));
    if(#bs==0, next);
    my(bt = if(issquare(d), (6*sqrtint(d))%20, Vec(bs)[1]));
    if(!setsearch(Set(Vec(bs)), bt), bt = Vec(bs)[1]);
    r = trD(D, bt);
    if(type(r[1])=="t_STR", next);
    v = 2*I*r[1];
    if(abs(imag(v)) > 1e-25*(1+abs(v)), listput(deg,d); next);
    v = round(real(v));
    tot++;
    if(v==0, listput(zer,d));
    if(v % d != 0, listput(bad,[d,v]));
    if(issquare(d), listput(sqok,[sqrtint(d), v, bet[sqrtint(d)]]));
  );
  print("s10: tested ", tot, " admissible d <= 300");
  print("  d | c(d) failures: ", #bad, if(#bad>0, Str(Vec(bad)), ""));
  print("  degenerate d (trace not purely imaginary): ", Vec(deg));
  print("  zeros of c: ", Vec(zer));
  print("  all zeros divisible by 5? ", if(#select(x->x%5!=0,Vec(zer))==0,"YES","no"));
  print("  squares [m, c(m^2), beta(m)]: ", Vec(sqok));
}
quit;
