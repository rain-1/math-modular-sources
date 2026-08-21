read("common.gp");
default(realprecision, 30);
{
for(i=1,#ROWS,
  my(R=ROWS[i], nm=R[1], r=R[2], a=R[3], b=R[4], c=R[5]);
  if(r!=3, next);
  my(BB=build(r,a,b,c,140), tq=truncate(BB[1]), dt=deriv(tq,'t));
  my(rts=polroots(1-2*a*x+c*x^2));
  print("=== ",nm,"  t_c ",rts~);
  my(sols=List());
  for(ir=2,36, my(rr=ir/50.);
    for(ia=0,143, my(th=ia*Pi/72., q0=rr*exp(I*th), ok=1);
      for(it=1,60, my(dv=subst(deriv(dt,'t),'t,q0)); if(dv==0,ok=0;break);
          q0 = q0 - subst(dt,'t,q0)/dv;
          if(abs(q0)>0.7 || abs(q0)<0.01, ok=0; break));
      if(ok && abs(subst(dt,'t,q0))<1e-15,
        my(new=1); for(j=1,#sols, if(abs(sols[j]-q0)<1e-12, new=0));
        if(new, listput(sols,q0)))));
  for(j=1,#sols, my(q0=sols[j]);
    print("   q_c=",q0,"  |q|=",abs(q0),"  t=",subst(tq,'t,q0),"  tau=",log(q0)/(2*Pi*I)));
);
}
quit;
