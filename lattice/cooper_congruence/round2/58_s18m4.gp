\\ 58_s18m4.gp -- s18 with the FUNDAMENTAL discriminant: family disc = -4d on Gamma_0(18),
\\ chi_{-4}, beta^2 = -4d mod 72.  The Cooper squares sit at d = 9m^2 (disc -36m^2).
\\ Admissible: d = 0,1 mod 4 (so D_1 = d is a discriminant) and -d a QR mod 18.
read("56_cdlib.gp");
default(realprecision, 60);
initfser(3,2200);
cont(Q) = gcd(gcd(Q[1],Q[2]),Q[3]);
{ TR4(d, bt, PB) = my(D=-4*d, RF=redforms(D), t=0., rep, ch, om, al, sq=sqrt(-D), nf=0, nad=0);
  for(i=1,#RF,
    rep = heegmin2(RF[i],18,bt,PB);
    if(rep==0, rep = heegmin2(RF[i],18,bt,3*PB));
    if(rep==0, next);          \\ no Heegner representative: skip (recorded separately)
    nf++;
    om = omeg(rep); ch = genchar(rep,-4);
    al = (-rep[2] + I*sq)/(2*rep[1]);
    t += ch*fhx(3,al)/om);
  [t,nf,#RF];
}
L=List();
{
for(d=1,200,
  if(d%4!=0 && d%4!=1, next);
  my(D=-4*d, bs=betas(D,18), best=0., bb=-1, nf=0, nt=0);
  if(#bs==0, next);
  for(i=1,#bs,
    my(R=TR4(d,bs[i],40));
    if(abs(R[1])>abs(best), best=R[1]; bb=bs[i]; nf=R[2]; nt=R[3]));
  my(v=best, tag);
  tag = if(abs(v)<1e-25, 0,
          if(abs(imag(v))>1e-22*(1+abs(v)), "CPLX",
             if(abs(real(v)-round(real(v)))<1e-22, round(real(v)),
                if(abs(real(v)/sqrt(3)-round(real(v)/sqrt(3)))<1e-22, Str(concat("sqrt3*",Str(round(real(v)/sqrt(3))))), "OTHER"))));
  listput(L,[d,bb,nf,nt,tag]);
);
}
V=Vec(L);
print("s18, family -4d:  ", #V, " admissible d <= 200");
print(V);
write("58_s18m4.txt", V);
print();
print("== the Cooper squares d = 9m^2 ==");
b18 = read("20_beta_s18.txt");
{ for(i=1,#V, my(d=V[i][1]);
   if(d%9!=0 || !issquare(d/9), next);
   my(m=sqrtint(d/9));
   print("   d=",d,"  m=",m,"  value=",V[i][5],"   beta(m)=",b18[m])); }
quit;
