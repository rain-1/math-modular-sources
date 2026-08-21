read("common.gp");
default(realprecision, 60);
NT = 300;
{
for(i=1,#ROWS,
  my(R=ROWS[i], nm=R[1], r=R[2], a=R[3], b=R[4], c=R[5]);
  my(BB=build(r,a,b,c,NT), tq=truncate(BB[1]), Fq=truncate(BB[2]));
  my(dt=deriv(tq,'t));
  my(Pc = if(r==2, 1-a*x+c*x^2, 1-2*a*x+c*x^2), rts=polroots(Pc));
  print("=== ",nm,"  r=",r,"  roots of P: ",rts~);
  \\ scan real q
  my(prev=0, found=[]);
  for(k=1,990, my(q0=k/1000., v=subst(dt,'t,q0));
     if(k>1 && sign(v)!=sign(prev), found=concat(found,[q0]));
     prev=v);
  print("   real sign changes of t'(q) at q ~ ", found);
  if(#found, my(q0=found[1]);
     for(it=1,60, q0 = q0 - subst(dt,'t,q0)/subst(deriv(dt,'t),'t,q0));
     print("   q_c = ", q0, "   t(q_c) = ", subst(tq,'t,q0), "   F(q_c)=",subst(Fq,'t,q0));
     print("   tau_c = ", log(q0)/(2*Pi)," i    (1/tau^2 -> N = ",-1/(log(q0)/(2*Pi))^2,")"));
  );
}
quit;
