read("/home/ubuntu/code/math-modular-sources/lattice/sources_s18_zud/zud_row.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/sources_s18_zud/ode_fit.gp");
NT=70; r=zudrow(NT); u=r[1];
U = vector(NT+1, i, 16^(i-1)*u[i]);
y = sum(n=0,NT-2, U[n+1]*x^n);
z = odefit(y,4,4); K=z[1]; lab=z[2];
v = K[,1]; v = v/content(v);
/* assemble polynomial coefficients p_i(t) of theta^i */
{ my(P=vector(5));
  for(b=1,#lab, P[lab[b][1]+1] += v[b]*x^(lab[b][2]));
  for(i=0,4, print("p_",i,"(t) = ", P[i+1]));
  print("leading symbol p_4(t) = ", factor(P[5]));
  print("p_0 factor: ", factor(P[1]));
}
quit;
