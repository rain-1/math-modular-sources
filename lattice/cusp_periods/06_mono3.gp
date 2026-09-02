/* 06_mono3.gp -- (P4) for the third-order (Almkvist-Zudilin) rows: the two finite
   singular points of P_3 = 1 - 2a x + c x^2 are FOLDS (order-2 orbifold points of
   the Fricke host), not cusps, so the monodromy there is not unipotent.  We
   compute (M-1)y anyway and compare with the q-side cusp periods.            */
read("lib.gp");
read("mono.gp");
default(parisizemax, 20000000000);
default(realprecision, 100);
NC = 520;
A3(a,b,c) = my(v=vector(NC+1)); v[1]=1; v[2]=b; for(n=1,NC-1, v[n+2]=((2*n+1)*(a*n^2+a*n+b)*v[n+1]-c*n^3*v[n])/(n+1)^3); v;
B3r(a,b,c) = my(v=vector(NC+1)); v[1]=0; v[2]=1; for(n=1,NC-1, v[n+2]=((2*n+1)*(a*n^2+a*n+b)*v[n+1]-c*n^3*v[n])/(n+1)^3); v;
PC3(a,b,c) = [[-b, c], [1, -(6*a+2*b), 7*c], [0, 3, -9*a, 6*c], [0, 0, 1, -2*a, c]];
RC0(x) = vector(MT+1, i, 0);
RC1(x) = vector(MT+1, i, if(i==1, 1, 0));
dl3(nm, a, b, c, xs, appr0, rho, Nc) = {
  my(PC = PC3(a,b,c), Av=A3(a,b,c), Bv=B3r(a,b,c));
  my(SG = concat([0], Vec(polroots(1 - 2*a*'x + c*'x^2))));
  my(xb = xs + rho);
  my(appr = concat(appr0, [xb]));
  my(circ = vector(Nc+1, i, xs + rho*exp(2*Pi*I*(i-1)/Nc)));
  my(path = concat(appr, circ[2..Nc+1]));
  my(x0 = appr[1]);
  my(sA = runp(PC, SG, RC0, evs(Av,x0,3), appr));
  my(sB = runp(PC, SG, RC1, evs(Bv,x0,3), appr));
  my(gA = runp(PC, SG, RC0, evs(Av,x0,3), path));
  my(gB = runp(PC, SG, RC1, evs(Bv,x0,3), path));
  print("  ", nm, "  x* = ", xs);
  print("     Delta A  = ", gA[1]-sA[1]);
  print("     Delta B  = ", gB[1]-sB[1]);
  print("     ratio    = ", (gB[1]-sB[1])/(gA[1]-sA[1]));
  print("     from y'  = ", (gB[2]-sB[2])/(gA[2]-sA[2]));
  print("     from y'' = ", (gB[3]-sB[3])/(gA[3]-sA[3]));
  (gB[1]-sB[1])/(gA[1]-sA[1]);
}
print("=== Apery gamma (17,5,1): P_3 = 1-34x+x^2, roots 17 -+ 12 sqrt2 = 0.0294373, 33.9706");
g1 = dl3("near", 17,5,1, 17-12*sqrt(2), [0.005, 0.005+0.008*I, 0.037+0.008*I], 0.008, 48);
print("     q-side Pi(cusp 0) for gamma = ", cperiod(srcbyname("gamma"),0,1)[1], "  = zeta(3)/6 = ", zeta(3)/6);
g2 = dl3("far ", 17,5,1, 17+12*sqrt(2), [0.005, 0.005+3*I, 45+3*I], 11.0, 64);
print("     q-side cusp values of gamma: 1/1 ", cperiod(srcbyname("gamma"),1,1)[1], "   1/2 ", cperiod(srcbyname("gamma"),1,2)[1]);
print("                                  1/3 ", cperiod(srcbyname("gamma"),1,3)[1], "   1/6 ", cperiod(srcbyname("gamma"),1,6)[1]);
print();
print("=== Domb alpha (10,4,64): P_3 = (1-4x)(1-16x), roots 1/16, 1/4");
a1 = dl3("near", 10,4,64, 1/16., [0.008, 0.008+0.012*I, 0.075+0.012*I], 0.012, 48);
a2 = dl3("far ", 10,4,64, 1/4., [0.008, 0.008+0.09*I, 0.33+0.09*I], 0.08, 48);
print("     q-side alpha: 1/1 ", cperiod(srcbyname("alpha"),1,1)[1], "   1/3 ", cperiod(srcbyname("alpha"),1,3)[1]);
print("                   1/4 ", cperiod(srcbyname("alpha"),1,4)[1], "   1/6 ", cperiod(srcbyname("alpha"),1,6)[1], "  1/12 ", cperiod(srcbyname("alpha"),1,12)[1]);
print();
print("=== AZ epsilon (12,4,16): P_3 = 1-24x+16x^2, roots (3 -+ 2 sqrt2)/4 = 0.0428932, 1.4571");
e1 = dl3("near", 12,4,16, (3-2*sqrt(2))/4, [0.008, 0.008+0.012*I, 0.055+0.012*I], 0.012, 48);
e2 = dl3("far ", 12,4,16, (3+2*sqrt(2))/4, [0.008, 0.008+0.6*I, 2.0+0.6*I], 0.5, 48);
print("     q-side eps: 1/1 ", cperiod(srcbyname("eps"),1,1)[1], "   1/2 ", cperiod(srcbyname("eps"),1,2)[1]);
print("                 1/4 ", cperiod(srcbyname("eps"),1,4)[1], "   1/8 ", cperiod(srcbyname("eps"),1,8)[1]);
quit;
