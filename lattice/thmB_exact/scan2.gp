read("common.gp");
default(realprecision, 40);
NT = 300;
{
for(i=1,6,
  my(R=ROWS[i], nm=R[1], r=R[2], a=R[3], b=R[4], c=R[5]);
  my(BB=build(r,a,b,c,NT), tq=truncate(BB[1]));
  my(dt=deriv(tq,'t));
  print("=== ",nm, "  t_c roots: ", polroots(1-a*x+c*x^2)~);
  print("  t(q) at q=-0.9..0.9: ", vector(19,k, my(q0=(k-10)/10.); [q0, subst(tq,'t,q0)]));
);
}
quit;
