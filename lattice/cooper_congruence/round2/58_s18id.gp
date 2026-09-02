\\ 58_s18id.gp -- identify the s18 traces T(d) for non-square d:  in which field do they live?
read("56_cdlib.gp");
default(realprecision, 90);
initfser(3,2200);
DS=[41,42,43,44,45,46,47,48,50,51,52,53,54,55,5,7,11,13,3,6,10,12];
{
for(j=1,#DS,
  my(d=DS[j], D=-36*d, bs=betas(D,18), best=0, bb=0);
  for(i=1,#bs,
    my(R=Traw(3,d,bs[i],40));
    if(type(R)=="t_STR", next);
    if(abs(R[1])>abs(best), best=R[1]; bb=bs[i]));
  if(abs(best)<1e-30, print("d=",d,"  T = 0 for every beta"); next);
  my(T=real(best), r1=T, r2=T/sqrt(3), r3=T/sqrt(d), r4=T/sqrt(3*d));
  print("d=",d," (beta=",bb,")   T=",T);
  print("     T        -> ",if(abs(r1-round(r1))<1e-25, concat("INTEGER ",Str(round(r1))), "-"));
  print("     T/sqrt3  -> ",if(abs(r2-round(r2))<1e-25, concat("INTEGER ",Str(round(r2))), "-"));
  print("     T/sqrt(d)-> ",if(abs(r3-round(r3))<1e-25, concat("INTEGER ",Str(round(r3))), "-"));
  print("     T/sqrt(3d)->",if(abs(r4-round(r4))<1e-25, concat("INTEGER ",Str(round(r4))), "-"));
  print("     lindep[1,sqrt3,T] = ", lindep([1,sqrt(3),T],30));
);
}
quit;
