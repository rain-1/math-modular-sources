/* 09_cooper.gp -- bonus: Cooper's three meromorphic-source rows s7, s10, s18,
   and Zagier B / row zeta, by analytic continuation only (no q-side formula
   exists for a meromorphic source).  Also fixes the far fold of row zeta.    */
read("lib.gp");
read("mono.gp");
default(parisizemax, 20000000000);
default(realprecision, 100);
NC = 520;
A3d(a,b,c,d) = my(v=vector(NC+1)); v[1]=1; v[2]=b; for(n=1,NC-1, v[n+2]=((2*n+1)*(a*n^2+a*n+b)*v[n+1]-n*(c*n^2+d)*v[n])/(n+1)^3); v;
B3d(a,b,c,d) = my(v=vector(NC+1)); v[1]=0; v[2]=1; for(n=1,NC-1, v[n+2]=((2*n+1)*(a*n^2+a*n+b)*v[n+1]-n*(c*n^2+d)*v[n])/(n+1)^3); v;
PC3d(a,b,c,d) = [[-b, c+d], [1, -(6*a+2*b), 7*c+d], [0, 3, -9*a, 6*c], [0, 0, 1, -2*a, c]];
RC0(x) = vector(MT+1, i, 0);
RC1(x) = vector(MT+1, i, if(i==1, 1, 0));
dl(nm, a, b, c, d, xs, appr0, rho, Nc) = my(PC=PC3d(a,b,c,d)); my(Av=A3d(a,b,c,d)); my(Bv=B3d(a,b,c,d)); my(SG=concat([0], Vec(polroots(1 - 2*a*'x + c*'x^2)))); my(appr=concat(appr0,[xs+rho])); my(circ=vector(Nc+1,i, xs+rho*exp(2*Pi*I*(i-1)/Nc))); my(path=concat(appr, circ[2..Nc+1])); my(x0=appr[1]); my(sA=runp(PC,SG,RC0,evs(Av,x0,3),appr)); my(sB=runp(PC,SG,RC1,evs(Bv,x0,3),appr)); my(gA=runp(PC,SG,RC0,evs(Av,x0,3),path)); my(gB=runp(PC,SG,RC1,evs(Bv,x0,3),path)); print("  ",nm," x* = ",xs); print("     A_0..A_4 = ", vector(5,i,Av[i])); print("     period  = ", (gB[1]-sB[1])/(gA[1]-sA[1])); print("     from y' = ", (gB[2]-sB[2])/(gA[2]-sA[2])); (gB[1]-sB[1])/(gA[1]-sA[1]);
G3 = znstar(3,1);
L2m3 = lfun(lfuncreate([G3,[1]]),3-1);
print("=== Cooper s7 (13,4,-27,3):  P_3 = 1-26x-27x^2, near 1/27, far -1");
c71 = dl("near", 13,4,-27,3, 1/27., [0.004, 0.004+0.008*I, 0.045+0.008*I], 0.008, 48);
print("     zeta(2)/7 = ", zeta(2)/7);
c72 = dl("far ", 13,4,-27,3, -1., [0.004, 0.004+0.3*I, -0.5+0.3*I], 0.5, 48);
print();
print("=== Cooper s10 (6,2,-64,4):  P_3 = 1-12x-64x^2, near 1/16, far -1/4");
cA1 = dl("near", 6,2,-64,4, 1/16., [0.006, 0.006+0.012*I, 0.0745+0.012*I], 0.012, 48);
print("     zeta(2)/5 = ", zeta(2)/5);
cA2 = dl("far ", 6,2,-64,4, -0.25, [0.006, 0.006+0.12*I, -0.13+0.12*I], 0.12, 48);
print();
print("=== Cooper s18 (14,6,192,-12): P_3 = 1-28x+192x^2, near 1/16, far 1/12");
cB1 = dl("near", 14,6,192,-12, 1/16., [0.005, 0.005+0.01*I, 0.0725+0.01*I], 0.01, 48);
print("     L(2,chi-3)/2 = ", L2m3/2);
cB2 = dl("far ", 14,6,192,-12, 1/12., [0.005, 0.005+0.012*I, 0.0953+0.012*I], 0.012, 48);
print();
print("=== row zeta (9,3,-27,0): P_3 = 1-18x-27x^2, near 0.0515669, far -0.718234");
z1 = dl("near", 9,3,-27,0, (-9+6*sqrt(3))/27, [0.005, 0.005+0.008*I, 0.0596+0.008*I], 0.008, 48);
z2 = dl("far ", 9,3,-27,0, (-9-6*sqrt(3))/27, [0.005, 0.005+0.3*I, -0.4+0.3*I], 0.3, 48);
print("     q-side zeta row: cusp 0 ", cperiod(srcbyname("zeta"),0,1)[1], "  cusp 1/3 ", cperiod(srcbyname("zeta"),1,3)[1], "  cusp 1/9 ", cperiod(srcbyname("zeta"),1,9)[1]);
quit;
