read("common.gp");
default(realprecision, 40);
NT = 700;
{
for(i=1,#ROWS,
  my(R=ROWS[i], nm=R[1], r=R[2], a=R[3], b=R[4], c=R[5]);
  if(nm!="delta" && nm!="eta" && nm!="B", next);
  my(BB=build(r,a,b,c,NT), tq=truncate(BB[1]), dt=deriv(tq,'t));
  print("=== ",nm);
  my(best=List());
  for(ix=0,50, for(iy=5,40,
    my(tau=ix/100. + I*iy/100., q0=exp(2*Pi*I*tau), v=abs(q0*subst(dt,'t,q0)));
    listput(best,[v,tau])));
  my(v=Vec(best)); v=vecsort(v,1);
  print("   16 smallest |q t'(q)|:");
  for(k=1,16, print("     ",v[k][1],"  tau=",v[k][2]));
);
}
quit;
